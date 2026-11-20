import json
from pathlib import Path

import pandas as pd

from bgsmap.cli import main


def test_end_to_end_cli_and_pdf(tmp_path, caplog, capsys):
    inp=tmp_path/"data"
    main(["demo","--out",str(inp)])
    capsys.readouterr()
    result=tmp_path/"results"
    with caplog.at_level("INFO", logger="bgsmap"):
        main(["run","--gss",str(inp/"GSS.tsv.gz"),"--burden",str(inp/"burden.tsv.gz"),
              "--coords",str(inp/"spatial_coords.tsv"),"--trait-name","TEST",
              "--out",str(result),"--permutations","39","--chunk-size","16"])
    summary = json.loads(capsys.readouterr().out)
    assert summary["n_genes"] == 160
    assert "Reading burden summary" in caplog.text
    assert "Permutation diagnostics: 16/39 completed" in caplog.text
    assert "Permutation diagnostics: 39/39 completed" in caplog.text
    assert "Finished bgsmap run" in caplog.text
    plots=result/"plots"
    for n in ("report.html","spatial_enrichment.html","region_enrichment.html",
              "spatial_enrichment.pdf","region_enrichment.pdf","gene_contributions.pdf"):
        assert (plots/n).stat().st_size>1000
    assert (plots/"spatial_enrichment.pdf").open("rb").read(4)==b"%PDF"
    assert "plotly" in (plots/"spatial_enrichment.html").read_text(encoding="utf-8")
    region=pd.read_csv(result/"region_results.tsv",sep="\t")
    assert str(region.loc[0,"domain"])=="3"
    assert {'p_analytic', 'logp_competitive_analytic', 'p_holm_competitive_analytic'}.issubset(region.columns)
    assert 'p_empirical' not in region.columns
    assert (result/'region_permutation_diagnostics.tsv').exists()
    assert (result/'spot_permutation_diagnostics.tsv.gz').exists()
    assert 'p_empirical' not in (plots/'region_enrichment.html').read_text()
    assert 'p_empirical' not in (plots/'spatial_enrichment.html').read_text()
    assert 'permutation p' not in (plots/'region_enrichment.html').read_text()
    assert 'p_empirical' not in (plots/'report.html').read_text()
    meta=json.loads((result/"run_metadata.json").read_text())
    assert meta["n_genes"]==160 and meta["n_spots"]==96
    main(["plot","--out",str(result),"--trait-name","New title"])
    assert "New title" in (plots/"report.html").read_text(encoding="utf-8")


def test_config_file(tmp_path):
    inp=tmp_path/"data"
    main(["demo","--out",str(inp)])
    out=tmp_path/"cfig"
    cfg={"gss":str(inp/"GSS.tsv.gz"),"burden":str(inp/"burden.tsv.gz"),
         "coords":str(inp/"spatial_coords.tsv"),"permutations":19,"out":str(out)}
    p=tmp_path/"options.json"
    p.write_text(json.dumps(cfg))
    main(["run","--config",str(p)])
    assert (out/"plots"/"report.html").exists()


def test_default_is_analytic_no_permutations(tmp_path):
    inp=tmp_path/'demo';main(['demo','--out',str(inp)])
    out=tmp_path/'no_perm'
    main(['run','--gss',str(inp/'GSS.tsv.gz'),'--burden',str(inp/'burden.tsv.gz'),
          '--coords',str(inp/'spatial_coords.tsv'),'--out',str(out)])
    assert not (out/'region_permutation_diagnostics.tsv').exists()
    region=pd.read_csv(out/'region_results.tsv',sep='\t')
    assert 'p_competitive_analytic' in region.columns
    assert 'p_competitive_empirical' not in region.columns
    assert 'logp_competitive_analytic' in region.columns
