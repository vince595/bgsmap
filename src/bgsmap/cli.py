"""Command-line entry point for reproducible burden × gsMap GSS analysis."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from time import perf_counter

from . import __version__
from .pipeline import run_pipeline, replot
from .demo import generate_demo

LOG = logging.getLogger(__name__)


def build_parser():
    parser=argparse.ArgumentParser(prog="bgsmap",description="Spatial gene-level burden enrichment with gsMap GSS (NOT SNP-level LDSC)")
    parser.add_argument("--version",action="version",version=f"bgsmap {__version__}")
    subs=parser.add_subparsers(dest="command",required=True)
    r=subs.add_parser("run",help="End-to-end QC → enrichment → HTML/PDF plots",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    r.add_argument("--config",help="JSON config keys matching the CLI option names (hyphens become underscores)")
    r.add_argument("--gss",help="gsMap latent-to-gene marker-score feather, or TSV/CSV")
    r.add_argument("--burden",help="Gene-level burden summary (tsv/csv[.gz/.xz], parquet)")
    r.add_argument("--h5ad",help="Original gsMap spatial H5AD; read obs and obsm without scanpy")
    r.add_argument("--coords",help="Alternative coordinates table (spot,x,y,domain)")
    r.add_argument("--spatial-key",default="spatial",help="Key under H5AD /obsm/")
    r.add_argument("--annotation",help="H5AD obs annotation key or coords TSV annotation column")
    r.add_argument("--gss-gene-col",help="GSS gene ID column (auto-detect by default)")
    r.add_argument("--gene-col",help="Burden gene ID column (auto-detect)")
    r.add_argument("--p-col",help="Burden p column (auto-detect)")
    r.add_argument("--mask-col",help="Burden mask column")
    r.add_argument("--mask",help="Select exactly one burden mask")
    r.add_argument("--trait-col",help="Burden trait/phenotype column")
    r.add_argument("--trait",help="Select one phenotype in multi-trait burden table")
    r.add_argument("--trait-name",default="BurdenTrait",help="Display name of output trait")
    r.add_argument("--gene-map",help="Two-column one-to-one ortholog file st_gene,burden_gene")
    r.add_argument("--duplicate-policy",choices=["error","bonferroni"],default="error")
    r.add_argument("--covariates",help="Separate gene-level matching covariate table")
    r.add_argument("--cov-gene-col",default="gene")
    r.add_argument("--match-cols",nargs="+",help="Power-related gene covariates; gene length, MAC, variants tested etc")
    r.add_argument("--match-bins",type=int,default=4,help="Quantile bins for each continuous matching column")
    r.add_argument("--score-mode",choices=["logp","binary"],default="logp")
    r.add_argument("--binary-threshold",type=float,default=0.05)
    r.add_argument("--p-cap",type=float,default=30,help="Cap -log10(burden p) to limit extreme outliers")
    r.add_argument("--min-genes",type=int,default=100)
    r.add_argument("--permutations",type=int,default=0,help="Optional permutation diagnostics, 0=analytic only (normal-tail P values are always primary)")
    r.add_argument("--seed",type=int,default=42)
    r.add_argument("--chunk-size",type=int,default=64)
    r.add_argument("--out",help="Output directory")
    r.add_argument("--no-invert-y",action="store_true",help="Do not invert histology-oriented y-axis")
    q=subs.add_parser("quick",help="Generate upstream gsMap GSS then run burden enrichment (requires gsMap)")
    # The 'quick' parser accepts all burden-run options (plus upstream inputs).
    for action in r._actions:
        if action.dest == "help":
            continue
        from copy import copy
        cloned=copy(action)
        q._add_action(cloned)
    q.add_argument("--workdir",required=True,help="Where upstream gsMap stores GNN/GSS files")
    q.add_argument("--sample-name",required=True,help="Upstream gsMap sample name")
    q.add_argument("--data-layer",default="X",help="Original gsMap counts layer; set 'counts' for raw counts")
    q.add_argument("--latent-key",help="Existing H5AD obsm latent key; skip GNN training")
    q.add_argument("--homolog-file",help="Upstream gsMap homolog file, species_gene HUMAN_GENE_SYM")
    q.add_argument("--epochs",type=int,default=300)
    q.add_argument("--force-gss",action="store_true",help="Recompute upstream gsMap output")
    p=subs.add_parser("plot",help="Regenerate all HTML/PDF plots from an existing results folder")
    p.add_argument("--out",required=True)
    p.add_argument("--trait-name",default="BurdenTrait")
    p.add_argument("--no-invert-y",action="store_true")
    i=subs.add_parser("inspect",help="Inspect H5AD obs/obsm and burden/GSS input columns")
    i.add_argument("--h5ad")
    i.add_argument("--coords", help="Optional coordinates TSV instead of --h5ad")
    i.add_argument("--spatial-key", default="spatial")
    i.add_argument("--annotation", help="H5AD obs annotation / coords domain column")
    i.add_argument("--burden")
    i.add_argument("--gss")
    i.add_argument("--check-spots", action="store_true",
                   help="Compare GSS Feather schema columns with spatial spot IDs without reading the matrix")
    d=subs.add_parser("demo",help="Create synthetic toy GSS/coords/burden inputs")
    d.add_argument("--out",required=True)
    d.add_argument("--seed",type=int,default=19)
    return parser


def _config_apply(parser, args):
    if args.command not in ("run","quick") or not args.config:
        return args
    conf=json.loads(Path(args.config).read_text(encoding="utf-8"))
    if not isinstance(conf,dict):
        raise ValueError("Config JSON must be an object")
    allowed=set(vars(args))
    unknown=set(conf)-allowed
    if unknown:
        raise ValueError(f"Unknown config keys: {sorted(unknown)}")
    # Do not override flags actually supplied on command line. Compare to parser defaults.
    defaults={a.dest:a.default for a in parser._subparsers._group_actions[0].choices[args.command]._actions}
    for key,val in conf.items():
        if key not in ("command","config") and getattr(args,key,None)==defaults.get(key):
            setattr(args,key,val)
    return args


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.getLogger("bgsmap").setLevel(logging.INFO)
    start = perf_counter()
    LOG.info("Starting bgsmap %s (%s)", args.command, __version__)
    if args.command=="inspect":
        from .io import read_table, feather_metadata, load_coords, load_coords_anndata
        info={}
        if args.h5ad:
            LOG.info("Inspecting spatial H5AD: %s", args.h5ad)
            try:
                import h5py
            except ImportError as e:
                raise ImportError("For H5AD inspection: pip install h5py") from e
            with h5py.File(args.h5ad,"r") as f:
                obs=f["obs"]
                index_name=obs.attrs.get("_index","_index")
                if isinstance(index_name,bytes):
                    index_name=index_name.decode("utf-8")
                info["h5ad"]={"obs_columns":list(obs.keys())[:40],
                    "obsm_keys":list(f["obsm"].keys()),
                    "n_spots":int(obs[index_name].shape[0])}
        if args.coords:
            LOG.info("Inspecting spatial coordinates: %s", args.coords)
            table=read_table(args.coords)
            info["coords"]={"rows":len(table),"columns":list(map(str,table.columns[:40]))}
        if args.burden:
            LOG.info("Inspecting burden summary: %s", args.burden)
            table=read_table(args.burden)
            info["burden"]={"rows":len(table),"columns":list(map(str,table.columns[:40]))}
        if args.gss:
            LOG.info("Inspecting GSS input: %s", args.gss)
            if str(args.gss).lower().endswith((".feather",".arrow")):
                schema=feather_metadata(args.gss)
                info["gss"]={"n_columns":schema["n_columns"],
                    "first_columns":schema["columns"][:8],
                    "last_columns":schema["columns"][-8:],
                    "record_batches":schema["n_record_batches"],
                    "mode":"schema-only; numeric GSS matrix not loaded"}
            else:
                table=read_table(args.gss)
                info["gss"]={"rows":len(table),"columns":list(map(str,table.columns[:40]))}
        if args.check_spots:
            if not args.gss or bool(args.h5ad)==bool(args.coords):
                parser.error("--check-spots requires --gss and exactly one of --h5ad / --coords")
            LOG.info("Comparing GSS and spatial spot IDs")
            if str(args.gss).lower().endswith((".feather",".arrow")):
                columns=feather_metadata(args.gss)["columns"]
            else:
                columns=list(read_table(args.gss).columns)
            for gene_name in ("HUMAN_GENE_SYM", "gene", "gene_symbol", "GENE", "index"):
                if gene_name in columns:
                    columns.remove(gene_name)
                    break
            spatial=load_coords(h5ad=args.h5ad,coords=args.coords,
                                spatial_key=args.spatial_key,annotation=args.annotation)
            matched=sum(x in spatial.index for x in columns)
            reader="h5py" if args.h5ad else "coords TSV"
            if matched<5 and args.h5ad:
                try:
                    alt=load_coords_anndata(args.h5ad,spatial_key=args.spatial_key,
                                            annotation=args.annotation)
                    alt_matched=sum(x in alt.index for x in columns)
                    if alt_matched>matched:
                        spatial,matched,reader=alt,alt_matched,"anndata fallback"
                except (ImportError,ValueError,KeyError) as exc:
                    info["fallback_warning"]=str(exc)
            info["spot_alignment"]={"n_gss_spots":len(columns),
                  "n_spatial_spots":len(spatial),"matched":matched,
                  "reader":reader,"gss_examples":columns[:5],
                  "spatial_examples":list(spatial.index[:5])}
        if not info:
            parser.error("inspect requires --h5ad, --coords, --burden, or --gss")
        print(json.dumps(info,indent=2,ensure_ascii=False))
    elif args.command=="demo":
        LOG.info("Generating synthetic inputs in %s", args.out)
        dest=generate_demo(args.out,args.seed);print(f"Synthetic demo created: {dest}")
    elif args.command=="plot":
        replot(args.out,args.trait_name,invert_y=not args.no_invert_y)
        print(f"Plots written to: {Path(args.out)/'plots'}")
    else:
        args=_config_apply(parser,args)
        if args.command=="quick":
            from .gsmap_bridge import precompute_gss
            if not args.h5ad:
                parser.error("quick requires --h5ad input for original gsMap")
            args.gss,args.h5ad=precompute_gss(args.h5ad,args.workdir,args.sample_name,
                data_layer=args.data_layer,annotation=args.annotation,latent_key=args.latent_key,
                homolog_file=args.homolog_file,epochs=args.epochs,skip_existing=not args.force_gss)
        if not args.gss or not args.burden or not args.out:
            parser.error("run requires --gss, --burden and --out (or config equivalents)")
        if bool(args.h5ad)==bool(args.coords):
            parser.error("run requires exactly one of --h5ad or --coords")
        if args.min_genes<20:
            parser.error("--min-genes must be >=20")
        if args.p_cap<=0 or args.binary_threshold<=0 or args.binary_threshold>=1:
            parser.error("Invalid p-cap/binary-threshold")
        details=run_pipeline(args)
        print(json.dumps({"results":str(Path(args.out).resolve()),"n_genes":details["n_genes"],
                          "n_spots":details["n_spots"],"report":str(Path(args.out).resolve()/"plots"/"report.html")},indent=2))
    LOG.info("Finished bgsmap %s in %.1f s", args.command, perf_counter() - start)


if __name__=="__main__":
    main()
