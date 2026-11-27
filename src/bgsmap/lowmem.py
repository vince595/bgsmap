"""Analytic GSS enrichment with a bounded queue of parallel spot blocks."""
from __future__ import annotations

from collections import deque
from concurrent.futures import ThreadPoolExecutor
import logging
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from . import __version__
from .gss_stream import analysis_source
from .resources import plan_blocks
from .io import load_burden, load_coords, read_gene_map
from .pipeline import _add_covariates, _align_spot_ids, _write_results
from .statistics import (
    EnrichmentResult,
    _analytic_tables,
    _burden_weights,
    _gene_contributions,
    _permutation_variance,
    _profile_variance,
    _result_metadata,
    _variance_factors,
)

LOG = logging.getLogger(__name__)


def _parallel_blocks(blocks, compute, threads, name, n_spots, metrics=None):
    """Read ahead by at most one block and reduce worker results in input order."""
    started = perf_counter()
    read_seconds = 0.0
    done = 0
    worker_seconds = 0.0
    n_blocks = 0
    progress = max(1, n_spots // 10)
    next_log = progress
    pending = deque()
    iterator = iter(blocks)

    def timed_compute(start, block):
        then = perf_counter()
        result = compute(start, block)
        return result, perf_counter() - then

    def report(task):
        nonlocal done, next_log, worker_seconds, n_blocks
        start, size, future = task
        result, seconds = future.result()
        worker_seconds += seconds
        n_blocks += 1
        if start != done:
            raise RuntimeError("Unexpected GSS block order")
        done += size
        if done >= next_log or done == n_spots:
            LOG.info("%s %d/%d spots (%.1f s)", name, done, n_spots, perf_counter() - started)
            next_log = done + progress
        return start, result

    with ThreadPoolExecutor(max_workers=threads) as executor:
        while True:
            then = perf_counter()
            try:
                start, block = next(iterator)
            except StopIteration:
                break
            read_seconds += perf_counter() - then
            if len(pending) == threads:
                yield report(pending.popleft())
            pending.append((start, len(block), executor.submit(timed_compute, start, block)))
        while pending:
            yield report(pending.popleft())
    if done != n_spots:
        raise RuntimeError(f"Read only {done}/{n_spots} GSS spots during {name}")
    if metrics is not None:
        metrics.update(wall_seconds=perf_counter() - started, read_seconds=read_seconds,
                       worker_seconds=worker_seconds, n_blocks=n_blocks, threads=threads)
    LOG.info("%s finished in %.2f s, including %.2f s pulling GSS blocks",
             name, perf_counter() - started, read_seconds)


def _gene_moments(blocks, n_genes, n_spots, threads, metrics=None):
    # Shifting by the first spot avoids cancellation for nearly invariant
    # nonzero genes. Each worker returns only vectors of length n_genes.
    anchor = np.empty(n_genes, dtype=np.float64)

    def with_anchor():
        for start, block in blocks:
            if start == 0:
                anchor[:] = block[0]
                anchor[~np.isfinite(anchor)] = 0.0
            yield start, block

    def moments(start, block):
        invalid = ~np.isfinite(block)
        centered = np.subtract(block, anchor, dtype=np.float64)
        centered[invalid] = 0
        return (invalid.any(axis=0), centered.sum(axis=0),
                np.einsum("ij,ij->j", centered, centered))

    sums = np.zeros(n_genes, dtype=np.float64)
    squares = np.zeros_like(sums)
    invalid = np.zeros(n_genes, dtype=bool)
    for _, (bad, block_sum, block_squares) in _parallel_blocks(
        with_anchor(), moments, threads, "Pass 1/2 — gene QC", n_spots, metrics
    ):
        invalid |= bad
        sums += block_sum
        squares += block_squares
    mean = sums / n_spots
    sd = np.sqrt(np.maximum(squares / n_spots - mean * mean, 0.0))
    good = ~invalid & (sd > 1e-9)
    return good, (anchor[good] + mean[good]).astype(np.float32), sd[good].astype(np.float32)


def _spot_statistics(blocks, means, stds, w, strata, codes, n_domains, threads, metrics=None):
    n_spots = len(codes)
    factors = _variance_factors(w, strata)
    weights = w.astype(np.float32)

    def score(start, block):
        block -= means
        block /= stds
        obs = (np.einsum("ij,j->i", block, weights, dtype=np.float64) / len(w)).astype(np.float32)
        var = _profile_variance(block, factors)
        labels = codes[start:start + len(block)]
        order = np.argsort(labels, kind="stable")
        ordered_labels = labels[order]
        groups = np.r_[0, np.flatnonzero(ordered_labels[1:] != ordered_labels[:-1]) + 1]
        if len(groups) == 1:
            partial = block.sum(axis=0, dtype=np.float64)[None, :]
        else:
            grouped = np.split(block[order], groups[1:], axis=0)
            partial = np.stack([group.sum(axis=0, dtype=np.float64) for group in grouped])
        return obs, var, ordered_labels[groups], partial

    obs = np.empty(n_spots, dtype=np.float32)
    var = np.empty(n_spots, dtype=np.float64)
    sums = np.zeros((n_domains, len(w)), dtype=np.float64)
    for start, (block_obs, block_var, present, partial) in _parallel_blocks(
        blocks, score, threads, "Pass 2/2 — spot statistics", n_spots, metrics
    ):
        obs[start:start + len(block_obs)] = block_obs
        var[start:start + len(block_var)] = block_var
        sums[present] += partial
    counts = np.bincount(codes, minlength=n_domains)
    profile = (sums / counts[:, None]).astype(np.float32)
    overall = (sums.sum(axis=0) / n_spots).astype(np.float32)
    # Same domain-vs-rest contrast as bgsmap 1.0.0. A single domain has
    # no comparison set and therefore a zero competitive profile.
    comp = np.vstack([
        row - ((n_spots * overall - count * row) / (n_spots - count))
        if count < n_spots else np.zeros_like(row)
        for row, count in zip(profile, counts)
    ]).astype(np.float32)
    return obs, var, profile, comp, counts


def run_pipeline_fast(args):
    started = perf_counter()
    with analysis_source(
        args.gss, gene_col=args.gss_gene_col, cache=getattr(args, "gss_cache", None),
        temp_dir=getattr(args, "temp_dir", None),
        projection_size=getattr(args, "projection_size", None),
        memory_limit=getattr(args, "memory_limit", None),
    ) as source:
        import_seconds = perf_counter() - started
        return _run_source(args, source, import_seconds)


def _run_source(args, source, import_seconds):
    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    burden, info = load_burden(args.burden, args.gene_col, args.p_col, args.mask_col,
                               args.mask, args.trait_col, args.trait, args.duplicate_policy)
    gene_map = read_gene_map(args.gene_map)
    mapped = pd.Index([gene_map.get(gene) if gene_map is not None else gene
                       for gene in source.genes], name="gene")
    if mapped.dropna().duplicated().any():
        raise ValueError("Gene map yielded duplicated target genes")
    burden = burden.set_index("gene", drop=False)
    rows = np.flatnonzero(mapped.isin(burden.index))
    if len(rows) < args.min_genes:
        raise ValueError(f"Only {len(rows)} matched genes (<{args.min_genes})")
    coords = load_coords(coords=args.coords, h5ad=args.h5ad, spatial_key=args.spatial_key,
                         annotation=args.annotation)
    all_spots = pd.Index(source.spots, name="spot")
    indices, coords = _align_spot_ids(
        all_spots, coords, h5ad=args.h5ad, spatial_key=args.spatial_key, annotation=args.annotation,
    )
    spots = all_spots[indices]
    n_spots = len(spots)
    LOG.info("Matched %d genes and %d spots before QC", len(rows), n_spots)

    labels = coords.loc[spots, "domain"].to_numpy(dtype=str)
    domains = pd.Index(pd.unique(labels), name="domain")
    settings = plan_blocks(
        source, len(rows), n_spots, len(domains),
        threads=getattr(args, "threads", None), block_size=getattr(args, "block_size", None),
        projection_size=getattr(args, "projection_size", None),
        memory_limit=getattr(args, "memory_limit", None),
    )
    settings.update(storage=source.storage, gss_cache=getattr(args, "gss_cache", None))
    LOG.info("Analysis resource plan %s", settings)
    performance = {"source_setup_seconds": import_seconds,
                   "alignment_seconds": perf_counter() - started,
                   "gene_qc": {}, "spot_statistics": {}}

    def blocks(genes):
        return source.blocks(indices, genes, settings["spot_block"], settings["projection_size"])

    good, means, stds = _gene_moments(
        blocks(rows), len(rows), n_spots, settings["threads"], performance["gene_qc"],
    )
    genes = rows[good]
    if len(genes) < args.min_genes:
        raise ValueError(f"Only {len(genes)} finite/nonconstant genes after QC (<{args.min_genes})")
    LOG.info("QC retained %d/%d genes", len(genes), len(rows))
    burden = _add_covariates(burden.loc[mapped[genes]].reset_index(drop=True), args)
    if not args.match_cols:
        LOG.warning("No gene-level power covariates matched. Unstratified gene-label permutation is exploratory.")
    w_raw, w, strata, matching = _burden_weights(
        burden, args.score_mode, args.binary_threshold, args.p_cap, args.match_cols, args.match_bins,
    )
    codes = pd.Categorical(labels, categories=domains).codes
    obs, var, profile, comp, counts = _spot_statistics(
        blocks(genes), means, stds, w, strata, codes, len(domains), settings["threads"],
        performance["spot_statistics"],
    )
    source.check_source()
    w32 = w.astype(np.float32)
    out_sp, out_reg = _analytic_tables(
        coords, spots, domains, counts, obs, var, (profile @ w32) / len(w),
        _permutation_variance(profile, w, strata), (comp @ w32) / len(w),
        _permutation_variance(comp, w, strata),
    )
    top_domain = str(out_reg.iloc[0].domain)
    out_gene = _gene_contributions(burden, w_raw, w, top_domain, domains, comp)
    meta = _result_metadata(
        n_spots, len(domains), w, 0, args.seed, args.score_mode, args.p_cap,
        args.binary_threshold, matching, top_domain,
    )
    details = {
        "tool": "bgsmap", "version": __version__, "input_gss": str(Path(args.gss).resolve()),
        "input_burden": str(Path(args.burden).resolve()),
        "input_spatial": str(Path(args.h5ad or args.coords).resolve()),
        "input_burden_column_mapping": info, "genes_before_qc": len(rows),
        "genes_after_qc": len(genes), "excluded_invalid_or_invariant_genes": len(rows) - len(genes),
        "n_input_gss_spots": n_spots, "trait_name": args.trait_name, **meta,
        "low_memory": settings, "performance": performance,
    }
    _write_results(EnrichmentResult(out_sp, out_reg, out_gene, meta), burden, details, out, args,
                   skip_plots=getattr(args, "skip_plots", False))
    return details
