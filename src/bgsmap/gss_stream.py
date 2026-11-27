"""Bounded Feather V2 projections and a reusable spot-major GSS cache."""
from __future__ import annotations

import errno
import json
import logging
from pathlib import Path
import shutil
import tempfile
from time import perf_counter

import numpy as np
import pandas as pd

from .io import _choose, canonical_gene

LOG = logging.getLogger(__name__)
CACHE_DTYPE = np.dtype("<f4")
CACHE_FORMAT = 1
GENE_COLUMNS = ("HUMAN_GENE_SYM", "gene", "gene_symbol", "GENE", "index")


def _source_identity(path):
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


class FeatherGSS:
    def __init__(self, path, gene_col=None):
        import pyarrow as pa

        self.path = Path(path).expanduser().resolve()
        if self.path.suffix.lower() not in (".feather", ".arrow"):
            raise ValueError("The fast engine requires Feather V2 / Arrow IPC file GSS")
        self.identity = _source_identity(self.path)
        with pa.memory_map(str(self.path), "r") as source:
            names = pa.ipc.open_file(source).schema.names
            self.gene_col = _choose(names, gene_col, GENE_COLUMNS, "gss-gene")
            self.gene_field = names.index(self.gene_col)
            self.spots = [name for i, name in enumerate(names) if i != self.gene_field]
            if len(set(names)) != len(names):
                raise ValueError("Duplicate Feather field names")
            reader = self._project(source, [self.gene_field])
            columns = (reader.get_batch(i).column(0).to_pylist()
                       for i in range(reader.num_record_batches))
            self.genes = list(canonical_gene(pd.Series(
                [gene for column in columns for gene in column], dtype="string")))
        if pd.isna(self.genes).any() or len(set(self.genes)) != len(self.genes):
            raise ValueError("Missing or duplicate GSS gene IDs after normalization")

    @staticmethod
    def _project(source, fields):
        from pyarrow import ipc

        opts = ipc.IpcReadOptions(included_fields=fields, use_threads=False)
        return ipc.open_file(source, options=opts)

    def blocks(self, spots, rows, block_size, projection_size):
        import pyarrow as pa

        # These indexes follow source order after barcode alignment. Arrow returns
        # included fields in schema order. Each mapping is released per projection.
        for start in range(0, len(spots), projection_size):
            positions = spots[start:start + projection_size]
            fields = (positions + (positions >= self.gene_field)).tolist()
            with pa.memory_map(str(self.path), "r") as source:
                reader = self._project(source, fields)
                batches = [reader.get_batch(i) for i in range(reader.num_record_batches)]
                try:
                    for offset in range(0, len(fields), block_size):
                        stop = min(offset + block_size, len(fields))
                        result = np.empty((stop - offset, len(rows)), dtype=np.float32)
                        for i in range(offset, stop):
                            if len(batches) == 1:
                                column = batches[0].column(i)
                            else:
                                column = pa.concat_arrays([batch.column(i) for batch in batches])
                            values = column.to_numpy(zero_copy_only=False)
                            result[i - offset] = values if len(rows) == len(values) else values[rows]
                        del column, values
                        yield start + offset, result
                finally:
                    del batches, reader

    def check_source(self):
        if self.identity != _source_identity(self.path):
            raise ValueError(f"GSS changed during reading: {self.path}")


class CachedGSS:
    def __init__(self, path, gss, gene_col=None):
        self.path = Path(path).expanduser().resolve()
        self.source = Path(gss).expanduser().resolve()
        try:
            manifest = json.loads((self.path / "manifest.json").read_text(encoding="utf-8"))
            if manifest["format"] != CACHE_FORMAT:
                raise ValueError("unsupported cache format")
            self.identity = manifest["source"]
            if gene_col and manifest["gene_column"] != gene_col:
                raise ValueError("the GSS gene column differs")
            self.genes = manifest["genes"]
            self.spots = manifest["spots"]
            self.matrix = self.path / "gss.f32"
            expected = CACHE_DTYPE.itemsize * len(self.genes) * len(self.spots)
            if not self.genes or not self.spots or self.matrix.stat().st_size != expected:
                raise ValueError("invalid cache dimensions or matrix size")
            self.check_source()
        except (ValueError, OSError, KeyError, TypeError) as exc:
            raise ValueError(f"Cannot use GSS cache {self.path}: {exc}. Rebuild in a new directory.") from exc

    def check_source(self):
        if self.identity != _source_identity(self.source):
            raise ValueError(f"Cache source does not match current GSS: {self.source}")

    def blocks(self, spots, rows, block_size, projection_size):
        width = len(self.genes)
        stride = width * CACHE_DTYPE.itemsize
        with self.matrix.open("rb") as source:
            for start in range(0, len(spots), block_size):
                positions = spots[start:start + block_size]
                breaks = np.flatnonzero(np.diff(positions) != 1) + 1
                runs = np.split(positions, breaks)
                # Common case uses one sequential read. Larger gaps are skipped
                # without allocating a dense matrix across non-aligned spots.
                result = np.empty((len(positions), len(rows)), dtype=np.float32) if len(runs) > 1 else None
                offset = 0
                for run in runs:
                    target = int(run[0]) * stride
                    if source.tell() != target:
                        source.seek(target)
                    raw = np.fromfile(source, dtype=CACHE_DTYPE, count=len(run) * width)
                    if raw.size != len(run) * width:
                        raise ValueError(f"Incomplete GSS matrix in {self.path}")
                    data = raw.reshape(len(run), width)
                    if len(rows) != width:
                        data = np.ascontiguousarray(data[:, rows], dtype=np.float32)
                    if result is None:
                        result = data
                    else:
                        result[offset:offset + len(run)] = data
                    offset += len(run)
                yield start, result


def build_cache(gss, out, gene_col=None, projection_size=4096):
    """Publish the complete cache atomically after reading one Feather pass."""
    target = Path(out).expanduser().resolve()
    if target.exists():
        CachedGSS(target, gss, gene_col)
        LOG.info("Reusing GSS cache %s", target)
        return target
    started = perf_counter()
    source = FeatherGSS(gss, gene_col)
    rows = np.arange(len(source.genes))
    spots = np.arange(len(source.spots))
    if not len(rows) or not len(spots):
        raise ValueError("GSS cache requires genes and spot columns")
    target.parent.mkdir(parents=True, exist_ok=True)
    required = CACHE_DTYPE.itemsize * len(rows) * len(spots)
    if shutil.disk_usage(target.parent).free < required:
        raise OSError(errno.ENOSPC, f"GSS cache requires {required:,} bytes", str(target.parent))
    LOG.info("Building GSS cache %s (%d genes × %d spots, %.2f GiB)",
             target, len(rows), len(spots), required / 2**30)
    temp = Path(tempfile.mkdtemp(prefix=f".{target.name}.", dir=target.parent))
    try:
        progress = max(projection_size, len(spots) // 10)
        next_log = progress
        with (temp / "gss.f32").open("wb") as dest:
            for start, block in source.blocks(spots, rows, 512, projection_size):
                block.astype(CACHE_DTYPE, copy=False).tofile(dest)
                seen = start + len(block)
                if seen >= next_log or seen == len(spots):
                    LOG.info("Cache build %d/%d spots (%.1f s)", seen, len(spots),
                             perf_counter() - started)
                    next_log = seen + progress
        source.check_source()
        manifest = {"format": CACHE_FORMAT, "source": source.identity,
                    "gene_column": source.gene_col, "genes": source.genes, "spots": source.spots}
        (temp / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False), encoding="utf-8")
        try:
            temp.rename(target)
        except OSError as exc:
            if exc.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                raise
            CachedGSS(target, gss, gene_col)
    finally:
        if temp.exists():
            shutil.rmtree(temp)
    LOG.info("GSS cache ready in %.1f s: %s", perf_counter() - started, target)
    return target
