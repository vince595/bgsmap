"""Regression coverage for the bounded Feather and cached analytic engines."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from bgsmap.cli import main
from bgsmap.lowmem import _gene_moments


@pytest.fixture
def inputs(tmp_path):
    h5py = pytest.importorskip("h5py")
    feather = pytest.importorskip("pyarrow.feather")
    rng = np.random.default_rng(31)
    genes = [f"G{i}" for i in range(84)]
    spots = [f"spot-{i}" for i in range(121)]
    matrix = rng.gamma(2, 3, size=(len(genes), len(spots))).astype(np.float16)
    matrix[:10, ::3] += 3
    matrix[2, :] = np.float16(0.125)
    matrix[3, 9] = np.nan
    matrix[5, 19] = np.inf
    gss = pd.DataFrame(matrix, columns=spots)
    gss.insert(17, "HUMAN_GENE_SYM", genes)
    feather.write_feather(gss, tmp_path / "multi.feather", compression="lz4", chunksize=13)
    burden = pd.DataFrame({"gene_symbol": genes[:-5],
                           "Pvalue_Burden": rng.uniform(0.001, 0.45, len(genes) - 5)})
    burden = pd.concat([burden, burden.iloc[[8]].assign(Pvalue_Burden=1e-5)], ignore_index=True)
    burden.to_csv(tmp_path / "burden.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": genes, "power": np.arange(len(genes)) % 2}).to_csv(
        tmp_path / "covariates.tsv", sep="\t", index=False,
    )
    keep = [i for i in rng.permutation(len(spots)) if i % 7 != 0]
    with h5py.File(tmp_path / "spatial.h5ad", "w") as h5:
        obs = h5.create_group("obs")
        obs["_index"] = np.array([spots[i].encode() for i in keep])
        annotation = obs.create_group("substructure")
        annotation["categories"] = [b"A", b"B", b"C"]
        annotation["codes"] = np.asarray(keep) % 3
        h5.create_group("obsm")["spatial"] = rng.random((len(keep), 2))
    return [
        "--gss", str(tmp_path / "multi.feather"), "--h5ad", str(tmp_path / "spatial.h5ad"),
        "--annotation", "substructure", "--burden", str(tmp_path / "burden.tsv"),
        "--gene-col", "gene_symbol", "--p-col", "Pvalue_Burden", "--min-genes", "20",
        "--duplicate-policy", "bonferroni", "--seed", "42", "--permutations", "0", "--skip-plots",
    ]


@pytest.mark.parametrize("mode", ["logp", "binary"])
def test_analytic_outputs_match_legacy(inputs, tmp_path, mode):
    common = [*inputs, "--score-mode", mode]
    if mode == "binary":
        common += ["--covariates", str(tmp_path / "covariates.tsv"), "--match-cols", "power",
                   "--binary-threshold", "0.2"]
    cache = str(tmp_path / "cache")
    main(["cache-gss", "--gss", common[1], "--out", cache, "--projection-size", "31"])
    main(["cache-gss", "--gss", common[1], "--out", cache])
    main(["run", *common, "--engine", "legacy", "--out", str(tmp_path / "legacy")])
    for name, flags in [
        ("direct", ["--engine", "fast", "--threads", "1", "--block-size", "17",
                    "--projection-size", "23"]),
        ("cached", ["--n-threads", "4", "--block-size", "11", "--gss-cache", cache]),
    ]:
        out = tmp_path / name
        main(["run", *common, *flags, "--out", str(out)], lowmem=(name == "cached"))
        for file, key in [
            ("spot_results.tsv.gz", "spot"),
            ("region_results.tsv", "domain"),
            ("gene_contributions.tsv.gz", "gene"),
            ("aligned_burden.tsv.gz", "gene"),
        ]:
            expected = pd.read_csv(tmp_path / "legacy" / file, sep="\t").set_index(key).sort_index()
            actual = pd.read_csv(out / file, sep="\t").set_index(key).sort_index()
            pd.testing.assert_frame_equal(actual, expected, rtol=3e-5, atol=2e-6)
        details = json.loads((out / "run_metadata.json").read_text())
        assert details["excluded_invalid_or_invariant_genes"] == 3
        assert details["n_genes"] == 76
        assert details["n_domains"] == 3
        assert details["random_seed"] == 42


def test_stale_or_truncated_cache_rejected(inputs, tmp_path):
    gss, cache = inputs[1], tmp_path / "cache"
    main(["cache-gss", "--gss", gss, "--out", str(cache)])
    kwargs = ["run", *inputs, "--gss-cache", str(cache), "--out", str(tmp_path / "out")]
    with (cache / "gss.f32").open("r+b") as file:
        size = file.seek(0, 2)
        file.truncate(size - 4)
    with pytest.raises(ValueError, match="matrix size"):
        main(kwargs, lowmem=True)
    main(["cache-gss", "--gss", gss, "--out", str(cache.with_name("valid"))])
    stat = os.stat(gss)
    os.utime(gss, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
    kwargs[kwargs.index("--gss-cache") + 1] = str(cache.with_name("valid"))
    with pytest.raises(ValueError, match="does not match"):
        main(kwargs, lowmem=True)


def test_qc_resolves_near_constant_genes():
    base = np.float32(1000)
    one_ulp = np.nextafter(base, np.float32(np.inf))
    matrix = np.array([[base, base, 0], [base, one_ulp, np.nan], [base, base, 0]])
    keep, means, std = _gene_moments(iter([(0, matrix[:2]), (2, matrix[2:])]), 3, 3, 2)
    np.testing.assert_array_equal(keep, [False, True, False])
    np.testing.assert_allclose(means, [base])
    np.testing.assert_allclose(std[0], np.std(matrix[:, 1].astype(np.float64)), rtol=1e-6)


def test_fast_rejects_permutation_diagnostics(inputs, tmp_path, capsys):
    args = list(inputs)
    args[args.index("--permutations") + 1] = "19"
    with pytest.raises(SystemExit, match="2"):
        main(["run", *args, "--engine", "fast", "--out", str(tmp_path / "out")])
    assert "requires --permutations 0" in capsys.readouterr().err


def test_single_domain_has_zero_competitive_signal(inputs, tmp_path):
    import h5py

    h5ad = inputs[inputs.index("--h5ad") + 1]
    with h5py.File(h5ad, "r+") as source:
        source["obs/substructure/codes"][:] = 0
    out = tmp_path / "single"
    main(["run", *inputs, "--engine", "fast", "--threads", "2",
          "--block-size", "16", "--out", str(out)])
    region = pd.read_csv(out / "region_results.tsv", sep="\t")
    assert len(region) == 1 and region["n_spots"].item() == 103
    assert region["competitive_score"].item() == 0
    assert region["null_sd_competitive_analytic"].item() == 0
    assert region["p_holm_competitive_analytic"].item() == 1


@pytest.mark.parametrize("suffix", [".parquet", ".tsv.gz", ".csv.xz"])
def test_auto_engine_formats_and_reusable_cache(inputs, tmp_path, suffix):
    from bgsmap.io import read_table

    frame = read_table(inputs[1])
    converted = tmp_path / ("input" + suffix)
    if suffix == ".parquet":
        # Several row groups and an interior gene field exercise Parquet projections.
        numeric = frame.select_dtypes(include="number").columns
        frame[numeric] = frame[numeric].astype(np.float32)
        frame.index = [f"old-row-{i}" for i in range(len(frame))]
        frame.to_parquet(converted, row_group_size=9, index=True)
    else:
        frame.to_csv(converted, sep="," if ".csv" in suffix else "\t", index=False)
    common = list(inputs)
    common[1] = str(converted)
    cache = tmp_path / "reusable"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    main(["run", *common, "--engine", "legacy", "--out", str(tmp_path / "reference")])
    for name, extra in [
        ("auto", ["--temp-dir", str(scratch)]),
        ("cached", ["--gss-cache", str(cache)]),
        ("reused", ["--gss-cache", str(cache)]),
    ]:
        out = tmp_path / name
        main(["run", *common, "--threads", "2", "--block-size", "13",
              "--projection-size", "37", "--memory-limit", "256M", "--out", str(out), *extra])
        for filename, key in [("spot_results.tsv.gz", "spot"), ("region_results.tsv", "domain"),
                              ("gene_contributions.tsv.gz", "gene"), ("aligned_burden.tsv.gz", "gene")]:
            expected = pd.read_csv(tmp_path / "reference" / filename, sep="\t").set_index(key).sort_index()
            actual = pd.read_csv(out / filename, sep="\t").set_index(key).sort_index()
            pd.testing.assert_frame_equal(actual, expected, rtol=3e-5, atol=2e-6)
        details = json.loads((out / "run_metadata.json").read_text())
        assert details["execution"] == {"requested_engine": "auto", "engine": "fast"}
        assert details["low_memory"]["memory_budget_bytes"] == 256 * 2**20
        assert details["low_memory"]["estimated_working_set_bytes"] <= 256 * 2**20
        timing = details["performance"]
        assert timing["total_seconds"] > timing["gene_qc"]["wall_seconds"] > 0
        assert timing["spot_statistics"]["n_blocks"] >= 8
        assert timing["write_tables_seconds"] > 0
    assert not list(scratch.iterdir())


def test_budget_resizes_tiles_and_rejects_oversized_explicit_request():
    from types import SimpleNamespace
    from bgsmap.resources import plan_blocks

    source = SimpleNamespace(storage="columnar", genes=range(17008))
    plan = plan_blocks(source, 17008, 563390, 81, threads=4, memory_limit="1G")
    assert plan["estimated_working_set_bytes"] <= 2**30
    assert plan["projection_size"] < 4096
    with pytest.raises(ValueError, match="exceeds --memory-limit"):
        plan_blocks(source, 17008, 563390, 81, threads=4, block_size=4096,
                    projection_size=4096, memory_limit="1G")


def test_cli_flag_overrides_config_even_when_equal_to_default(inputs, tmp_path):
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"engine": "legacy", "threads": 1, "permutations": 19,
                                  "out": str(tmp_path / "result")}))
    main(["run", *inputs, "--config", str(config), "--engine", "auto"])
    details = json.loads((tmp_path / "result" / "run_metadata.json").read_text())
    assert details["execution"]["engine"] == "fast"
    assert details["permutations"] == 0


def test_temporary_text_cache_is_cleaned_after_failure(inputs, tmp_path):
    from bgsmap.io import read_table

    text = tmp_path / "input.tsv.gz"
    read_table(inputs[1]).to_csv(text, sep="\t", index=False)
    args = list(inputs)
    args[1] = str(text)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    with pytest.raises(ValueError, match="exceeds --memory-limit"):
        main(["run", *args, "--temp-dir", str(scratch), "--memory-limit", "1M",
              "--out", str(tmp_path / "result")])
    assert not list(scratch.iterdir())


def test_python_pipeline_uses_auto_engine(inputs, tmp_path):
    from bgsmap.cli import build_parser
    from bgsmap.pipeline import run_pipeline

    args = build_parser().parse_args(["run", *inputs, "--out", str(tmp_path / "result")])
    details = run_pipeline(args)
    assert details["execution"]["engine"] == "fast"
    assert args.engine == "auto"  # Caller-owned options retain the requested value.


def test_auto_preserves_feather_v1_compatibility(inputs, tmp_path):
    import pyarrow.feather as feather
    from bgsmap.io import read_table

    frame = read_table(inputs[1])
    numeric = frame.select_dtypes(include="number").columns
    frame[numeric] = frame[numeric].astype(np.float32)
    v1 = tmp_path / "old.feather"
    feather.write_feather(frame, v1, version=1)
    args = list(inputs)
    args[1] = str(v1)
    out = tmp_path / "old-result"
    main(["run", *args, "--out", str(out)])
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["execution"]["engine"] == "legacy"
    with pytest.raises(ValueError, match="requires Feather V2"):
        main(["run", *args, "--engine", "fast", "--out", str(tmp_path / "fast")])


def test_quick_dispatches_through_auto_engine(inputs, tmp_path, monkeypatch):
    def precompute(*args, **kwargs):
        return inputs[1], inputs[inputs.index("--h5ad") + 1]

    monkeypatch.setattr("bgsmap.gsmap_bridge.precompute_gss", precompute)
    out = tmp_path / "quick"
    main(["quick", *inputs, "--workdir", str(tmp_path / "work"), "--sample-name", "generic",
          "--out", str(out)])
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["execution"] == {"requested_engine": "auto", "engine": "fast"}
