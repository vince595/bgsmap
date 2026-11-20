"""Validated, unambiguous input adapters for gene-level burden and gsMap GSS."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    name = path.name.lower()
    if name.endswith((".feather", ".arrow")):
        try:
            return pd.read_feather(path)
        except Exception as first_error:
            # Some cluster builds can read Feather V2 through Arrow IPC while
            # pandas.read_feather fails; keep this fallback explicit.
            try:
                import pyarrow as pa
                with pa.memory_map(str(path), "r") as source:
                    return pa.ipc.open_file(source).read_all().to_pandas()
            except Exception as second_error:
                raise RuntimeError(
                    "Failed reading GSS Feather via pandas AND Arrow IPC. "
                    "Install a compatible pyarrow/libstdc++ (e.g. pip install 'bgsmap[feather]') "
                    "or convert GSS to TSV/CSV in a working environment. "
                    f"pandas: {first_error}; ipc: {second_error}"
                ) from second_error
    if name.endswith(".parquet"):
        try:
            return pd.read_parquet(path)
        except ImportError as e:
            raise ImportError("Parquet input requires pyarrow") from e
    sep = "," if any(name.endswith(x) for x in (".csv", ".csv.gz", ".csv.xz", ".csv.bz2")) else "\t"
    return pd.read_csv(path, sep=sep, compression="infer", low_memory=False)


def feather_metadata(path: str | Path) -> dict:
    """Inspect Arrow IPC/Feather schema without loading a giant GSS matrix.

    This matters for MOSTA E16.5: its GSS has >120,000 spot columns.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        import pyarrow as pa
        import pyarrow.ipc as ipc
        with pa.memory_map(str(path), "r") as source:
            reader = ipc.open_file(source)
            columns = list(reader.schema.names)
            # Arrow file-level record batches expose metadata without converting
            # the entire gene x spot table to pandas.
            return {"n_columns": len(columns), "columns": columns,
                    "n_record_batches": reader.num_record_batches}
    except ImportError as e:
        raise ImportError("Feather metadata inspection requires pyarrow") from e


def _validated_coordinates(tab: pd.DataFrame) -> pd.DataFrame:
    tab = tab.copy()
    tab["spot"] = tab["spot"].astype(str)
    tab["domain"] = tab["domain"].fillna("unassigned").astype(str)
    if tab.spot.duplicated().any():
        raise ValueError("Duplicate spatial spot identifiers")
    if not np.isfinite(tab[["x", "y"]].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite spatial coordinates")
    return tab.set_index("spot", drop=False)


def load_coords_anndata(h5ad, spatial_key="spatial", annotation=None):
    """Read coordinates using AnnData's own obs_names decoder.

    Used as a safety fallback if the fast h5py reader and Feather spot IDs
    disagree. AnnData is already an upstream gsMap dependency and understands
    all supported H5AD index encodings. The expression matrix is *not* read.
    """
    try:
        import anndata as ad
    except ImportError as e:
        raise ImportError(
            "AnnData-based H5AD spot ID recovery needs anndata. "
            "Install anndata or pass a verified --coords TSV instead."
        ) from e
    obj = ad.read_h5ad(h5ad, backed="r")
    try:
        spots = np.asarray(obj.obs_names.astype(str))
        if spatial_key not in obj.obsm:
            raise ValueError(f"H5AD missing obsm['{spatial_key}']")
        xy = np.asarray(obj.obsm[spatial_key], dtype=float)
        if xy.ndim != 2 or xy.shape[1] < 2 or xy.shape[0] != len(spots):
            raise ValueError("AnnData spatial coordinates must be N x 2+ with N obs_names")
        if annotation:
            if annotation not in obj.obs.columns:
                raise ValueError(f"H5AD annotation '{annotation}' missing in obs")
            domains = obj.obs[annotation].astype("string").fillna("unassigned").astype(str).to_numpy()
        else:
            domains = np.repeat("all", len(spots))
        tab = pd.DataFrame({"spot": spots, "x": xy[:, 0], "y": xy[:, 1], "domain": domains})
    finally:
        obj.file.close()
    return _validated_coordinates(tab)


def canonical_gene(s: pd.Series) -> pd.Series:
    # Strip Ensembl version only for Ensembl identifiers; do not munge gene symbols.
    z = s.astype("string").str.strip()
    z = z.str.replace(r"^(ENS[A-Z]*G\d+)\.\d+$", r"\1", regex=True)
    return z.str.upper()


def _choose(columns, requested, candidates, label):
    if requested:
        if requested not in columns:
            raise ValueError(f"Missing {label} column '{requested}'. Available: {list(columns)[:30]}")
        return requested
    lookup = {str(x).lower(): x for x in columns}
    for x in candidates:
        if x.lower() in lookup:
            return lookup[x.lower()]
    raise ValueError(f"Cannot detect {label} column. Pass --{label}-col. Available: {list(columns)[:30]}")


def load_burden(
    path, gene_col=None, p_col=None, mask_col=None, mask=None,
    trait_col=None, trait=None, duplicate_policy="error",
):
    df = read_table(path)
    g = _choose(df.columns, gene_col, ["gene", "gene_symbol", "genename", "symbol", "GENE", "HUMAN_GENE_SYM", "gene_id"], "gene")
    p = _choose(df.columns, p_col, ["p", "pvalue", "p_value", "pval", "Pvalue", "P_VALUE", "p.value", "P_BURDEN", "pval_burden"], "p")
    if mask is not None:
        m = _choose(df.columns, mask_col, ["mask", "mask_name", "annotation", "gene_mask"], "mask")
        df = df.loc[df[m].astype(str) == str(mask)].copy()
    if trait is not None:
        t = _choose(df.columns, trait_col, ["trait", "phenotype", "pheno", "phenotype_id"], "trait")
        df = df.loc[df[t].astype(str) == str(trait)].copy()
    if df.empty:
        raise ValueError("No burden rows after selecting trait/mask")
    input_rows = len(df)
    df = df.copy()
    df["gene"] = canonical_gene(df[g])
    df["p"] = pd.to_numeric(df[p], errors="coerce")
    df = df.loc[df.gene.notna() & df.p.notna()].copy()
    if (df.p < 0).any() or (df.p > 1).any():
        raise ValueError("Burden p values must be in [0,1]")
    # Preserve literal input P (including underflowed zeros) in the output.
    # Numerical clipping is performed only inside the derived -log10 score.
    df["p_original"] = df["p"]
    duplicates = df.duplicated("gene", keep=False)
    if duplicates.any():
        if duplicate_policy == "error":
            examples = df.loc[duplicates, "gene"].unique()[:8]
            raise ValueError(f"Duplicate burden genes {list(examples)}. Filter one mask/trait or use --duplicate-policy bonferroni")
        if duplicate_policy != "bonferroni":
            raise ValueError("Unknown duplicate policy")
        # Conservative correction for selecting a minimum among multiple gene tests.
        counts = df.groupby("gene")["p"].size()
        df = df.sort_values("p").drop_duplicates("gene").copy()
        df["p"] = np.minimum(1.0, df["p"] * df["gene"].map(counts))
    return df.reset_index(drop=True), {"gene_column": str(g), "p_column": str(p), "input_rows": input_rows, "retained_rows": len(df)}


def load_gss(path, gene_col=None):
    df = read_table(path)
    g = _choose(df.columns, gene_col, ["HUMAN_GENE_SYM", "gene", "gene_symbol", "GENE", "index"], "gene")
    genes = canonical_gene(df[g]); ids = [str(x) for x in df.columns if x != g]
    if genes.duplicated().any():
        raise ValueError("Duplicate GSS gene IDs after normalization. Resolve ortholog/alias collisions.")
    if pd.isna(genes).any():
        raise ValueError("Missing GSS gene IDs")
    if len(ids) < 5:
        raise ValueError("GSS must contain at least five spot columns")
    try:
        mat = df[ids].to_numpy(dtype=np.float32)
    except (TypeError, ValueError) as e:
        raise ValueError("GSS spot columns must be numeric") from e
    return pd.Index(genes.astype(str), name="gene"), pd.Index(ids, name="spot"), mat


def load_coords(coords=None, h5ad=None, spatial_key="spatial", annotation=None):
    if bool(coords) == bool(h5ad):
        raise ValueError("Supply exactly one of --coords or --h5ad")
    if coords:
        tab = read_table(coords)
        needed = ["spot", "x", "y"]
        if not all(k in tab for k in needed):
            raise ValueError("Coordinates table must have columns: spot, x, y [and optionally domain]")
        cols = needed + ([annotation] if annotation else (["domain"] if "domain" in tab else []))
        tab = tab[cols].copy()
        if annotation:
            tab = tab.rename(columns={annotation: "domain"})
        if "domain" not in tab:
            tab["domain"] = "all"
    else:
        try:
            import h5py
        except ImportError as e:
            raise ImportError("H5AD requires pip install 'bgsmap[h5ad]'") from e
        with h5py.File(h5ad, "r") as f:
            if "obs" not in f or "obsm" not in f or spatial_key not in f["obsm"]:
                raise ValueError(f"H5AD missing /obs, /obsm/{spatial_key}")
            obs = f["obs"]
            index_key = obs.attrs.get("_index", "_index")
            if isinstance(index_key, bytes):
                index_key = index_key.decode()
            if index_key not in obs:
                raise ValueError(f"H5AD missing obs index {index_key}")
            spots = _decode_h5(obs[index_key][...])
            xy = np.asarray(f["obsm"][spatial_key][...], dtype=float)
            if xy.ndim != 2 or xy.shape[1] < 2:
                raise ValueError("Spatial coordinates must be an Nx2+ matrix")
            if annotation:
                if annotation not in obs:
                    raise ValueError(f"H5AD annotation '{annotation}' missing in obs")
                node = obs[annotation]
                if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
                    codes = np.asarray(node["codes"][...], dtype=int)
                    categories = _decode_h5(node["categories"][...])
                    domains = np.array([categories[i] if 0 <= i < len(categories) else "unassigned" for i in codes])
                elif isinstance(node, h5py.Dataset):
                    domains = _decode_h5(node[...]).astype(str)
                else:
                    raise ValueError("Unsupported H5AD obs annotation representation")
            else:
                domains = np.repeat("all", len(spots))
            tab = pd.DataFrame({"spot":spots, "x":xy[:,0], "y":xy[:,1], "domain":domains})
    return _validated_coordinates(tab)


def _decode_h5(x):
    # Explicitly handle bytes and NumPy's bytes_ scalar (both common H5AD
    # encodings) without accidentally stringifying them as "b'123_456'".
    return np.array([
        v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v)
        for v in x
    ], dtype=str)


def read_gene_map(path):
    if not path:
        return None
    m = read_table(path)
    if not {"st_gene", "burden_gene"}.issubset(m):
        raise ValueError("Gene map must contain 'st_gene' and 'burden_gene' (one-to-one orthologs)")
    m = m[["st_gene", "burden_gene"]].dropna().copy()
    m["st_gene"] = canonical_gene(m["st_gene"])
    m["burden_gene"] = canonical_gene(m["burden_gene"])
    if m.duplicated("st_gene").any() or m.duplicated("burden_gene").any():
        raise ValueError("Gene map must be one-to-one, no duplicated source or target genes")
    return dict(zip(m.st_gene, m.burden_gene))
