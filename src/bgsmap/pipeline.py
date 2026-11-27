"""End-to-end data alignment, QC, enrichment, persistence, and visualization."""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .io import (
    canonical_gene, load_burden, load_coords, load_coords_anndata,
    load_gss, read_gene_map, read_table,
)
from .statistics import run_enrichment
from .plotting import write_plots

LOG = logging.getLogger(__name__)


def _align_spot_ids(spots, coords, h5ad=None, spatial_key="spatial", annotation=None):
    """Return index positions in GSS order and coordinates in spot-ID space.

    Do not infer alignment from row order: mixed-order H5AD and Feather inputs
    must always be matched using exact spot identifier strings.
    """
    keep_spots=np.flatnonzero(spots.isin(coords.index))
    if len(keep_spots)<5 and h5ad:
        # The index observed through a low-level H5AD decoder can occasionally
        # differ from AnnData's obs_names, even though the same source file was
        # used to produce gsMap Feather GSS. Recover the canonical obs_names
        # instead of editing IDs or assuming row order.
        LOG.warning(
            "Fast H5AD reader matched only %d/%d GSS spot IDs; "
            "retrying with AnnData obs_names.", len(keep_spots), len(spots)
        )
        try:
            alternate = load_coords_anndata(
                h5ad, spatial_key=spatial_key, annotation=annotation
            )
            better = np.flatnonzero(spots.isin(alternate.index))
            if len(better)>len(keep_spots):
                coords, keep_spots = alternate, better
                LOG.info("AnnData reader recovered %d/%d GSS spot IDs", len(better), len(spots))
        except (ImportError, ValueError, KeyError) as exc:
            LOG.warning("AnnData H5AD fallback unavailable: %s", exc)
    if len(keep_spots)<5:
        raise ValueError(
            f"Only {len(keep_spots)} GSS spots match spatial coordinates "
            f"(GSS={len(spots)}, spatial={len(coords)}). "
            f"First GSS spot IDs: {list(spots[:5])}; "
            f"first spatial IDs: {list(coords.index[:5])}. "
            "If the IDs are identical in AnnData, ensure anndata is installed; "
            "otherwise check the --gss / --h5ad input combination."
        )
    LOG.info("Spot alignment: %d/%d GSS spots matched %d spatial IDs",
             len(keep_spots), len(spots), len(coords))
    return keep_spots, coords


def run_pipeline(args):
    out=Path(args.out).expanduser().resolve()
    out.mkdir(parents=True,exist_ok=True)
    LOG.info("Output directory: %s", out)
    LOG.info("Reading burden summary: %s", args.burden)
    burden, info=load_burden(args.burden,args.gene_col,args.p_col,args.mask_col,args.mask,
                             args.trait_col,args.trait,args.duplicate_policy)
    LOG.info("Loaded %d burden genes. Reading GSS matrix: %s", len(burden), args.gss)
    gene_ix,spots,gss=load_gss(args.gss,args.gss_gene_col)
    LOG.info("Loaded GSS for %d genes and %d spots", len(gene_ix), len(spots))
    gene_map=read_gene_map(args.gene_map)
    if gene_map is not None:
        keep=np.array([g in gene_map for g in gene_ix]); gss=gss[keep]
        gene_ix=pd.Index([gene_map[g] for g in gene_ix[keep]],name="gene")
    if gene_ix.duplicated().any():
        raise ValueError("Gene map yielded duplicated target genes")
    LOG.info("Reading spatial coordinates: %s", args.h5ad or args.coords)
    coords=load_coords(coords=args.coords,h5ad=args.h5ad,spatial_key=args.spatial_key,
                       annotation=args.annotation)
    keep_spots,coords=_align_spot_ids(
        spots,coords,h5ad=args.h5ad,spatial_key=args.spatial_key,annotation=args.annotation
    )
    spots=spots[keep_spots]; gss=gss[:,keep_spots]
    burden=burden.set_index("gene",drop=False)
    keep_genes=np.flatnonzero(gene_ix.isin(burden.index))
    gss=gss[keep_genes];gene_ix=gene_ix[keep_genes]
    if len(gene_ix)<args.min_genes:
        raise ValueError(f"Only {len(gene_ix)} intersecting genes < --min-genes {args.min_genes}; check Ensembl/symbol and ortholog mapping")
    LOG.info("Checking GSS quality for %d intersecting genes", len(gene_ix))
    # Input GSS missing/constant markers do not encode localized information.
    finite=np.isfinite(gss).all(axis=1)
    nonconstant=np.std(np.where(np.isfinite(gss),gss,0),axis=1)>1e-9
    good=finite & nonconstant
    if np.count_nonzero(good) < args.min_genes:
        raise ValueError("Insufficient nonconstant genes after QC")
    skipped=int(np.sum(~good))
    gss=gss[good];gene_ix=gene_ix[good]
    burden=burden.loc[list(gene_ix)].reset_index(drop=True)
    burden = _add_covariates(burden, args)
    # Matching null is only plausible when power-related gene covariates are supplied.
    if not args.match_cols:
        LOG.warning("No gene-level power covariates matched. Unstratified gene-label permutation is exploratory.")
    LOG.info("Matched %d genes and %d spots; %d genes excluded for invalid/invariant GSS",len(gene_ix),len(spots),skipped)
    r=run_enrichment(gss, gene_ix,spots,burden,coords,
                     permutations=args.permutations,seed=args.seed,chunk_size=args.chunk_size,
                     p_cap=args.p_cap,score_mode=args.score_mode,binary_threshold=args.binary_threshold,
                     matching_cols=args.match_cols,match_bins=args.match_bins)
    details={"tool":"bgsmap","version":__version__,"input_gss":str(Path(args.gss).resolve()),
             "input_burden":str(Path(args.burden).resolve()),
             "input_spatial":str(Path(args.h5ad or args.coords).resolve()),
             "input_burden_column_mapping":info,"genes_before_qc":len(keep_genes),
             "genes_after_qc":len(gene_ix),"excluded_invalid_or_invariant_genes":skipped,
             "n_input_gss_spots":len(keep_spots),"trait_name":args.trait_name,
             **r.meta}
    _write_results(r, burden, details, out, args, skip_plots=getattr(args, 'skip_plots', False))
    return details


def _add_covariates(burden, args):
    if args.covariates:
        LOG.info("Reading gene matching covariates: %s", args.covariates)
        cv=read_table(args.covariates)
        cv_gene=args.cov_gene_col
        if cv_gene not in cv:
            raise ValueError(f"Missing covariate gene column '{cv_gene}'")
        cv["gene"]=canonical_gene(cv[cv_gene]);
        if cv.gene.duplicated().any():
            raise ValueError("Duplicate covariate gene identifiers")
        match_cols=args.match_cols or []
        for col in match_cols:
            if col not in cv:
                raise ValueError(f"Covariate file missing '{col}'")
        if set(match_cols) & set(burden.columns):
            raise ValueError("Matching column appears in both burden and covariate inputs; rename ambiguous column")
        burden=burden.merge(cv[["gene"]+match_cols],on="gene",how="left",validate="one_to_one")
    return burden


def _write_results(r, burden, details, out, args, skip_plots=False):
    LOG.info("Writing result tables and metadata to %s", out)
    for src,target,compression in ((r.spot,"spot_results.tsv.gz","gzip"),
                                    (r.region,"region_results.tsv",None),
                                    (r.gene,"gene_contributions.tsv.gz","gzip")):
        src.to_csv(out/target,sep="\t",index=False,compression=compression)
    if r.spot_permutation is not None:
        r.spot_permutation.to_csv(out/"spot_permutation_diagnostics.tsv.gz",sep="\t",index=False,compression="gzip")
        r.region_permutation.to_csv(out/"region_permutation_diagnostics.tsv",sep="\t",index=False)
    burden.to_csv(out/"aligned_burden.tsv.gz",sep="\t",index=False,compression="gzip")
    (out/"run_metadata.json").write_text(json.dumps(details,indent=2,ensure_ascii=False),encoding="utf-8")
    if skip_plots:
        LOG.info("Result tables complete")
        return
    LOG.info("Rendering HTML reports and vector PDFs")
    write_plots(r.spot,r.region,r.gene,out/"plots",trait=args.trait_name,
                invert_y=not args.no_invert_y)
    LOG.info("Complete. HTML: %s | PDF: %s",out/"plots"/"report.html",out/"plots"/"spatial_enrichment.pdf")


def replot(out, trait_name="Trait", invert_y=True):
    out=Path(out)
    LOG.info("Reading result tables from %s", out)
    spot=pd.read_csv(out/"spot_results.tsv.gz",sep="\t")
    region=pd.read_csv(out/"region_results.tsv",sep="\t")
    gene=pd.read_csv(out/"gene_contributions.tsv.gz",sep="\t")
    if not {'logp_analytic','p_analytic'}.issubset(spot.columns) or not {'logp_competitive_analytic','p_competitive_analytic'}.issubset(region.columns):
        raise ValueError("These result tables lack analytic P values. Re-run 'bgsmap run' to generate analytic results before plotting.")
    LOG.info("Rendering HTML reports and vector PDFs in %s", out / "plots")
    write_plots(spot,region,gene,out/"plots",trait=trait_name,invert_y=invert_y)
