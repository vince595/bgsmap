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
    main(["run", *common, "--out", str(tmp_path / "legacy")])
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
