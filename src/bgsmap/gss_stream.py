"""Bounded GSS readers and reusable disk caches for spatial datasets."""
from __future__ import annotations

import errno
from contextlib import contextmanager
import json
import logging
import mmap
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
    storage = "columnar"
    def __init__(self, path, gene_col=None):
        import pyarrow as pa

        self.path = Path(path).expanduser().resolve()
        if self.path.suffix.lower() not in (".feather", ".arrow"):
            raise ValueError("The fast engine requires Feather V2 / Arrow IPC file GSS")
        self.identity = _source_identity(self.path)
        with self.path.open("rb") as source:
            if source.read(4) == b"FEA1":
                raise ValueError("Fast computation requires Feather V2. Use --engine legacy "
                                 "for Feather V1, or convert the source to Feather V2.")
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
        _validate_genes(self.genes)

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


def _validate_genes(genes):
    if pd.isna(genes).any() or len(set(genes)) != len(genes):
        raise ValueError("Missing or duplicate GSS gene IDs after normalization")


class ParquetGSS:
    storage = "columnar"

    def __init__(self, path, gene_col=None):
        import pyarrow.parquet as pq

        self.path = Path(path).expanduser().resolve()
        self.identity = _source_identity(self.path)
        with pq.ParquetFile(self.path) as source:
            metadata = source.schema_arrow.pandas_metadata or {}
            index_fields = {field for field in metadata.get("index_columns", [])
                            if isinstance(field, str)}
            names = [field for field in source.schema_arrow.names if field not in index_fields]
            self.gene_col = _choose(names, gene_col, GENE_COLUMNS, "gss-gene")
            if len(set(names)) != len(names):
                raise ValueError("Duplicate Parquet field names")
            self.spots = [name for name in names if name != self.gene_col]
            self.genes = list(canonical_gene(source.read(
                columns=[self.gene_col], use_threads=False).to_pandas()[self.gene_col]))
        _validate_genes(self.genes)

    def blocks(self, spots, rows, block_size, projection_size):
        import pyarrow.parquet as pq

        with pq.ParquetFile(self.path) as source:
            for start in range(0, len(spots), projection_size):
                names = [self.spots[i] for i in spots[start:start + projection_size]]
                table = source.read(columns=names, use_threads=False)
                for offset in range(0, len(names), block_size):
                    stop = min(offset + block_size, len(names))
                    result = np.empty((stop - offset, len(rows)), dtype=np.float32)
                    for i in range(offset, stop):
                        values = table.column(i).combine_chunks().to_numpy(zero_copy_only=False)
                        result[i - offset] = values[rows]
                    del values
                    yield start + offset, result
                del table

    def check_source(self):
        if self.identity != _source_identity(self.path):
            raise ValueError(f"GSS changed during reading: {self.path}")


class TextGSS:
    storage = "text"

    def __init__(self, path, gene_col=None):
        self.path = Path(path).expanduser().resolve()
        self.identity = _source_identity(self.path)
        name = self.path.name.lower()
        self.separator = "," if name.endswith((".csv", ".csv.gz", ".csv.xz", ".csv.bz2")) else "\t"
        names = list(pd.read_csv(self.path, sep=self.separator, nrows=0).columns)
        self.gene_col = _choose(names, gene_col, GENE_COLUMNS, "gss-gene")
        self.spots = [name for name in names if name != self.gene_col]
        self.genes = list(canonical_gene(pd.read_csv(
            self.path, sep=self.separator, usecols=[self.gene_col])[self.gene_col]))
        _validate_genes(self.genes)

    def write_cache(self, dest, memory_limit):
        # One decompression/numeric parse into sequential gene-major disk data.
        # The next two analysis passes read only the selected spot/gene tiles.
        from .resources import available_memory, memory_bytes

        budget = memory_bytes(memory_limit) or min(4 * 2**30, available_memory() // 2)
        chunk_rows = max(1, min(1024, budget // max(1, len(self.spots) * 64)))
        offset = 0
        with pd.read_csv(self.path, sep=self.separator, chunksize=chunk_rows) as reader:
            for chunk in reader:
                genes = list(canonical_gene(chunk[self.gene_col]))
                if genes != self.genes[offset:offset + len(chunk)]:
                    raise ValueError("GSS gene rows changed during cache creation")
                values = chunk[self.spots].to_numpy(dtype=CACHE_DTYPE)
                values.tofile(dest)
                offset += len(chunk)
                LOG.info("Cache build %d/%d genes", offset, len(self.genes))
        if offset != len(self.genes):
            raise ValueError("Incomplete text GSS during cache creation")

    def check_source(self):
        if self.identity != _source_identity(self.path):
            raise ValueError(f"GSS changed during reading: {self.path}")


def open_gss(path, gene_col=None):
    name = str(path).lower()
    if name.endswith((".feather", ".arrow")):
        return FeatherGSS(path, gene_col)
    if name.endswith(".parquet"):
        return ParquetGSS(path, gene_col)
    return TextGSS(path, gene_col)


class CachedGSS:
    storage = "cache"
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
            self.layout = manifest.get("layout", "spot-major")
            if self.layout not in ("spot-major", "gene-major"):
                raise ValueError("unsupported matrix layout")
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
        if self.layout == "gene-major":
            yield from self._gene_major_blocks(spots, rows, block_size)
            return
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

    def _gene_major_blocks(self, spots, rows, block_size):
        shape = (len(self.genes), len(self.spots))
        with self.matrix.open("rb") as source:
            for start in range(0, len(spots), block_size):
                # Close each mapping before yielding. Only the copied tile
                # remains resident while workers process it.
                with mmap.mmap(source.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
                    matrix = np.ndarray(shape, dtype=CACHE_DTYPE, buffer=mapped)
                    tile = matrix[np.ix_(rows, spots[start:start + block_size])]
                    result = np.ascontiguousarray(tile.T, dtype=np.float32)
                    del matrix, tile
                yield start, result


def build_cache(gss, out, gene_col=None, projection_size=None, memory_limit=None):
    """Publish a reusable cache atomically after a bounded source conversion."""
    target = Path(out).expanduser().resolve()
    if target.exists():
        CachedGSS(target, gss, gene_col)
        LOG.info("Reusing GSS cache %s", target)
        return target
    started = perf_counter()
    source = open_gss(gss, gene_col)
    rows = np.arange(len(source.genes))
    spots = np.arange(len(source.spots))
    if not len(rows) or not len(spots):
        raise ValueError("GSS cache requires genes and spot columns")
    target.parent.mkdir(parents=True, exist_ok=True)
    required = CACHE_DTYPE.itemsize * len(rows) * len(spots)
    if shutil.disk_usage(target.parent).free < required:
        raise OSError(errno.ENOSPC, f"GSS cache requires {required:,} bytes", str(target.parent))
    if source.storage == "text":
        block_size = 512
        projection_size = 512
        layout = "gene-major"
    else:
        from .resources import plan_blocks

        settings = plan_blocks(source, len(rows), len(spots), 0, threads=1,
                               projection_size=projection_size, memory_limit=memory_limit)
        block_size = settings["spot_block"]
        projection_size = settings["projection_size"]
        layout = "spot-major"
    LOG.info("Building GSS cache %s (%d genes × %d spots, %.2f GiB)",
             target, len(rows), len(spots), required / 2**30)
    temp = Path(tempfile.mkdtemp(prefix=f".{target.name}.", dir=target.parent))
    try:
        progress = max(projection_size, len(spots) // 10)
        next_log = progress
        with (temp / "gss.f32").open("wb") as dest:
            if source.storage == "text":
                source.write_cache(dest, memory_limit)
            else:
                for start, block in source.blocks(spots, rows, block_size, projection_size):
                    block.astype(CACHE_DTYPE, copy=False).tofile(dest)
                    seen = start + len(block)
                    if seen >= next_log or seen == len(spots):
                        LOG.info("Cache build %d/%d spots (%.1f s)", seen, len(spots),
                                 perf_counter() - started)
                        next_log = seen + progress
        source.check_source()
        manifest = {"format": CACHE_FORMAT, "source": source.identity,
                    "gene_column": source.gene_col, "genes": source.genes, "spots": source.spots,
                    "layout": layout}
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


@contextmanager
def analysis_source(gss, *, gene_col=None, cache=None, temp_dir=None,
                    projection_size=None, memory_limit=None):
    """Provide a projected source or a temporary/persistent numeric cache."""
    if cache:
        build_cache(gss, cache, gene_col, projection_size, memory_limit)
        yield CachedGSS(cache, gss, gene_col)
        return
    if str(gss).lower().endswith((".feather", ".arrow", ".parquet")):
        yield open_gss(gss, gene_col)
        return
    with tempfile.TemporaryDirectory(prefix="bgsmap-gss-", dir=temp_dir) as temp:
        target = Path(temp) / "cache"
        LOG.info("Importing text GSS into temporary disk cache %s", target)
        build_cache(gss, target, gene_col, projection_size, memory_limit)
        yield CachedGSS(target, gss, gene_col)
