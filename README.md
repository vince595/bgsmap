# bgsmap

[简体中文](https://github.com/vince595/bgsmap/blob/main/README_CN.md)

bgsmap integrates gene-level rare-variant burden summaries with gene specificity scores (GSS) from [gsMap](https://github.com/JianYang-Lab/gsMap). It computes enrichment per spatial spot and domain and exports TSV, interactive HTML, and vector PDF.

Primary P values use normal-tail approximations and the in-memory engine saves optional empirical permutation diagnostics in separate files. See the [statistical model](https://github.com/vince595/bgsmap/blob/main/docs/statistics.md) or the [P-value migration guide (Chinese)](https://github.com/vince595/bgsmap/blob/main/docs/PVALUES_v0.2_CN.md).

Version 1.3.0 enables automatic bounded-memory analytic computation through
`bgsmap run` for Feather V2, Parquet, TSV and CSV. Spot-block workers, adaptive
tiles, reusable disk caches and phase timings are part of the software.
See the [large-dataset performance guide](docs/PERFORMANCE_CN.md).

## Install

Python 3.10 or later. Once bgsmap is published to PyPI, install it with either command:

```bash
python -m pip install bgsmap
python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple bgsmap
```

The Tsinghua mirror makes new releases available after it syncs with PyPI. Until then, install from the release wheel or an unpacked source tree:

```bash
python -m pip install ./bgsmap-1.3.0-py3-none-any.whl
# Or run this from the unpacked source tree:
python -m pip install .
```

For Feather and H5AD support, install `'bgsmap[all]'` from PyPI or `'.[all]'` from source. The `feather` extra also enables Parquet reading. The optional AnnData spot-ID fallback requires a separate `anndata` installation. The `quick` command additionally requires `gsMap` in a compatible Python environment.

```bash
bgsmap --version
python -m bgsmap --help
```

## Workflow

```text
Spatial H5AD ── gsMap GNN latent ── gsMap latent-to-gene GSS (gene × spot)
                                                    │
Gene-level burden summary (gene, p, mask) ── gene ID alignment and QC
                                                    │
                                GSS standardization + burden score
                                                    │
                           analytic conditional null (optional permutation diagnostics)
                                                    │
                     spot / region / competitive continuous P, logP, BH, Holm
                                                    │
                         interactive HTML + vector PDF + TSV
```

Crucially, **no GWAS SNP-to-gene annotation, PLINK LD weights, LD-score generation, or S-LDSC are performed on burden summaries**.

## 1. Demonstration without real data

```bash
bgsmap demo --out examples/demo_data
bgsmap run \
  --gss examples/demo_data/GSS.tsv.gz \
  --burden examples/demo_data/burden.tsv.gz \
  --coords examples/demo_data/spatial_coords.tsv \
  --trait-name SyntheticTrait \
  --permutations 0 \
  --out results/demo
```

The demo generates 160 genes, 96 spots, and a planted Domain 3 signal.

Commands print timestamped progress to stderr. The `run` and `inspect` commands write JSON to stdout. To save both streams separately, append `> summary.json 2> run.log` to a `run` command.

## 2. From an existing gsMap GSS matrix (recommended)

Upstream gsMap's actual output format is a Feather file, with a `HUMAN_GENE_SYM` gene column and one column per spatial barcode:

```text
{workdir}/{sample}/latent_to_gene/{sample}_gene_marker_score.feather
```

Use the **matching** H5AD (same spot barcodes) and a real burden summary. Select the actual annotation and column names by inspection:

```bash
bgsmap inspect --h5ad /path/to/Visium.h5ad --burden /path/to/Finngen_gene_burden.tsv.gz

bgsmap run \
  --gss /path/to/Visium_gene_marker_score.feather \
  --h5ad /path/to/Visium_add_latent.h5ad \
  --annotation domain \
  --burden /path/to/Finngen_gene_burden.tsv.gz \
  --gene-col gene --p-col p \
  --trait-name G6_MIGRAINE \
  --permutations 0 --seed 42 \
  --out results/G6_MIGRAINE
```

The H5AD reader loads only `/obs` and `/obsm/spatial` (not the expression matrix), and supports AnnData categorical obs encoding. If your domain annotation is stored elsewhere, export `spot,x,y,domain` as TSV and use `--coords` instead of `--h5ad`. A multi-domain annotation is needed for the *competitive* region comparison.

To select one mask or trait from a multi-mask/multi-trait summary use `--mask-col Mask --mask Mask1.0.01` or `--trait-col trait --trait G6_MIGRAINE`. Ambiguous duplicate gene tests are **rejected by default**; the optional `--duplicate-policy bonferroni` conservatively corrects the selected minimum p-value for the number of tests per gene. Never mix masks without documenting that strategy.

### Matching burden test power (strongly recommended)

A gene's burden p-value distribution can depend on gene length, number of qualifying variants, carrier counts, coverage, variant mask and other factors. Permuting all gene labels without correction is exploratory; for interpretable inference, match genes on available testing-power covariates:

```bash
bgsmap run \
  --gss /path/to/GSS.feather --h5ad /path/to/Visium.h5ad --annotation domain \
  --burden /path/to/burden.tsv.gz --gene-col Gene --p-col P \
  --covariates /path/to/gene_test_covariates.tsv \
  --match-cols gene_length n_tested_variants \
  --match-bins 4 --permutations 0 \
  --out results/matched
```

`gene_test_covariates.tsv` needs `gene,gene_length,n_tested_variants`. These columns may instead be present in the burden summary itself, in which case omit `--covariates`. Each joint stratum must contain at least 12 matched genes (small strata error rather than silent fallback). For actual burden results, a **matched null is preferred to an unrestricted gene-label shuffle**, but it is not a proof of calibration; gene-gene correlation and differences in burden test power can remain.

### Mouse tissue versus human burden

Original gsMap `latent-to-gene` can convert mouse genes to human symbols with a species–human homolog file. If your existing GSS is *not* already human-annotated, pass a strictly one-to-one table:

```tsv
st_gene	burden_gene
MouseSymbol1	HumanSymbol1
MouseSymbol2	HumanSymbol2
```

and add `--gene-map /path/to/one_to_one.tsv`. This pipeline never silently transforms mouse gene symbols into human genes. The `run_metadata.json` file records the retained gene/spot counts.

## 3. From raw spatial H5AD (reuse original gsMap)

```bash
pip install gsMap

bgsmap quick \
  --h5ad /path/to/V1_Adult_Mouse_Brain_Coronal_Section.h5ad \
  --workdir work/gsMap_precompute --sample-name VisiumCoronal \
  --data-layer counts --annotation domain \
  --homolog-file /path/to/mouse_human_homolog.tsv \
  --burden /path/to/FinnGen_MIGRAINE.tsv.gz \
  --gene-col gene --p-col p \
  --permutations 0 --out results/coronal
```

`quick` uses upstream `gsMap.find_latent_representation` and `gsMap.latent_to_gene` and then calls this package's gene-level inference. **It has not been integration-tested against your server's gsMap version.** It does not call gsMap's GWAS analysis. To skip expensive GNN training when H5AD already has an embedding, pass `--latent-key latent_GVAE` (replace with your actual `obsm` key). Upstream gsMap may copy the H5AD into its workdir; keep sufficient disk space. Use `--force-gss` to recompute upstream steps.

## 4. Configuration / figures

```bash
# Edit the input and column names first
bgsmap run --config examples/visium_finngen.json
# Regenerate plots from tables that include analytic P values
bgsmap plot --out results/coronal --trait-name G6_MIGRAINE
```

Output directory:

```text
results/coronal/
├── aligned_burden.tsv.gz
├── spot_results.tsv.gz                # spot, x/y, domain, T, p_analytic, logp_analytic, q_bh_analytic, p_holm_analytic
├── region_results.tsv                # region / competitive analytic p, logp, BH, Holm
├── region_permutation_diagnostics.tsv # optional, only if permutations>0
├── spot_permutation_diagnostics.tsv.gz # optional, only if permutations>0
├── gene_contributions.tsv.gz         # descriptive gene contributions
├── run_metadata.json                 # seeds, method and gene/spots overlap
└── plots/
    ├── report.html                   # offline report, links to interactive figures
    ├── spatial_enrichment.html       # offline Plotly, interactive (self-contained)
    ├── spatial_enrichment.pdf        # vector point symbols, two spatial panels
    ├── region_enrichment.html
    ├── region_enrichment.pdf
    └── gene_contributions.pdf
```

PDFs use Matplotlib vector artists and embedded TrueType text (`pdf.fonttype=42`), editable with Adobe Illustrator (subject to font availability). HTML is standalone and can be opened locally without a web server. Copy the **entire** `plots` folder to preserve report iframe references.

## Statistical model

See [the statistical model](https://github.com/vince595/bgsmap/blob/main/docs/statistics.md) for the exact conditional variance formula, Gaussian-tail approximation, competitive test, within-run multiplicity and validity limitations. Evaluate tail calibration on the intended data before drawing research conclusions.

## Development

```bash
python -m pip install -e '.[all,dev]'
python -m pytest -q
python -m build
```

CI covers Python 3.10 and 3.12. Tests exercise the synthetic analysis and the HTML/PDF outputs.

## License and citation

bgsmap is released under the [MIT license](https://github.com/vince595/bgsmap/blob/main/LICENSE). When using GSS computed by gsMap, cite the upstream paper. [Song et al., *Nature* 641, 932–941 (2025).](https://www.nature.com/articles/s41586-025-08757-x)
