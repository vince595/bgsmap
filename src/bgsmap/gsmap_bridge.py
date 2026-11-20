"""Optional integration with the *unmodified* upstream gsMap GNN/GSS steps.

Only runs the original expression-side steps. Does not invoke SNP annotations,
LD score construction or spatial S-LDSC on gene-level burden statistics.
"""
from __future__ import annotations

import logging
from pathlib import Path

LOG=logging.getLogger(__name__)


def precompute_gss(input_h5ad, workdir, sample_name, data_layer="X", annotation=None,
                   latent_key=None, homolog_file=None, epochs=300, skip_existing=True):
    try:
        from gsMap.config import FindLatentRepresentationsConfig, LatentToGeneConfig
        from gsMap.find_latent_representation import run_find_latent_representation
        from gsMap.latent_to_gene import run_latent_to_gene
    except ImportError as e:
        raise ImportError("Original gsMap package is required for GNN/GSS generation: pip install gsMap. "
                          "Alternatively supply an existing --gss file to 'bgsmap run'.") from e
    base=Path(workdir).expanduser().resolve()/sample_name
    h5ad_out=base/"find_latent_representations"/f"{sample_name}_add_latent.h5ad"
    gss_out=base/"latent_to_gene"/f"{sample_name}_gene_marker_score.feather"
    if not Path(input_h5ad).is_file():
        raise FileNotFoundError(input_h5ad)
    if latent_key is None and not (skip_existing and h5ad_out.exists()):
        LOG.info("Training original gsMap latent representation on %s",input_h5ad)
        conf=FindLatentRepresentationsConfig(
            input_hdf5_path=str(input_h5ad),workdir=str(workdir),sample_name=sample_name,
            annotation=annotation,data_layer=data_layer,epochs=epochs)
        run_find_latent_representation(conf)
    else:
        LOG.info("Using existing/provided latent representation")
    if not (skip_existing and gss_out.exists()):
        LOG.info("Computing gene specificity scores with upstream gsMap")
        kwargs={"workdir":str(workdir),"sample_name":sample_name,
                "annotation":annotation,"homolog_file":homolog_file}
        if latent_key is not None:
            kwargs["latent_representation"]=latent_key
            kwargs["input_hdf5_path"]=str(input_h5ad)
        conf=LatentToGeneConfig(**kwargs)
        run_latent_to_gene(conf)
    else:
        LOG.info("Using existing gene specificity scores: %s", gss_out)
    if not gss_out.is_file() or not h5ad_out.is_file():
        raise RuntimeError("Upstream gsMap GSS/latent outputs missing; inspect gsMap logs")
    return str(gss_out),str(h5ad_out)
