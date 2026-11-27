"""Compare bgsmap 1.2.0 with the original low-memory wheel on synthetic GSS."""
import argparse
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
from time import perf_counter

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as feather

# Measure each process separately. The old code's --skip-plots is an env var.
RUN = """
import json, resource, runpy, sys
module = sys.argv.pop(1)
try:
    runpy.run_module(module, run_name="__main__", alter_sys=True)
finally:
    print("PEAK_RSS", json.dumps(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
          file=sys.stderr)
"""


def measure(python, module, args, root, name, env):
    started = perf_counter()
    with (root / f"{name}.log").open("w") as log:
        subprocess.run([python, "-c", RUN, module, *args], stdout=log, stderr=log,
                       env=env, check=True)
    elapsed = perf_counter() - started
    log = (root / f"{name}.log").read_text()
    rss = int(re.search(r"PEAK_RSS (\d+)", log)[1])
    unit = 2**20 if sys.platform == "darwin" else 2**10
    phases = [list(map(float, values)) for values in re.findall(
        r"finished in ([\d.]+) s, including ([\d.]+) s pulling GSS blocks", log,
    )]
    result = {"seconds": round(elapsed, 3), "peak_rss_mib": round(rss / unit, 1),
              "passes_wall_and_read_seconds": phases}
    print(name, result, flush=True)
    return result


def make_inputs(root, genes, spots, domains):
    rng = np.random.default_rng(42)
    gene_names = [f"G{i}" for i in range(genes)]
    spot_names = [f"S{i:06}" for i in range(spots)]
    # Match gsMap's spot columns and a one-batch LZ4-compressed Feather V2.
    # Generate in spot-major order so that Arrow can use zero-copy columns
    # while writing the fixture. File creation is outside the measurements.
    markers = np.empty((spots, genes), dtype=np.float16)
    for start in range(0, spots, 256):
        tile = markers[start:start+256]
        tile[:] = rng.gamma(1.5, 2.0, size=tile.shape)
    markers[::domains, :20] += 2
    table = pa.Table.from_arrays([pa.array(gene_names), *map(pa.array, markers)],
                                 names=["HUMAN_GENE_SYM", *spot_names])
    feather.write_feather(table, root / "gss.feather", compression="lz4", chunksize=genes)
    pd.DataFrame({"spot": spot_names, "x": np.arange(spots) % 1024,
                  "y": np.arange(spots) // 1024,
                  "domain": [f"D{i % domains:02}" for i in range(spots)]}
                 ).to_csv(root / "coords.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": gene_names,
                  "Pvalue_Burden": 10.0 ** -rng.uniform(0, 5, size=genes)}
                 ).to_csv(root / "burden.tsv", sep="\t", index=False)


def compare(root, reference, other):
    max_diff = {}
    for file, key in [
        ("spot_results.tsv.gz", "spot"),
        ("region_results.tsv", "domain"),
        ("gene_contributions.tsv.gz", "gene"),
        ("aligned_burden.tsv.gz", "gene"),
    ]:
        a = pd.read_csv(root / reference / file, sep="\t").set_index(key).sort_index()
        b = pd.read_csv(root / other / file, sep="\t").set_index(key).sort_index()
        pd.testing.assert_frame_equal(a, b, rtol=5e-5, atol=3e-6)
        numeric = a.select_dtypes(include="number").columns
        max_diff[file] = float((a[numeric] - b[numeric]).abs().to_numpy().max())
    return max_diff


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="New benchmark directory")
    parser.add_argument("--genes", type=int, default=2048)
    parser.add_argument("--spots", type=int, default=40000)
    parser.add_argument("--domains", type=int, default=81)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--block-size", type=int, default=512)
    parser.add_argument("--projection-size", type=int, default=4096)
    parser.add_argument("--baseline-python", default=sys.executable,
                        help="Python with bgsmap 1.0.0 and bgsmap-lowmem 1.1.0")
    parser.add_argument("--baseline-path", help="Optional baseline-only PYTHONPATH")
    parser.add_argument("--check-legacy", action="store_true",
                        help="Also run original bgsmap 1.0.0 with plots, outside the speed comparison")
    args = parser.parse_args()
    root = args.out.resolve()
    root.mkdir(parents=True)
    make_inputs(root, args.genes, args.spots, args.domains)
    env = os.environ.copy()
    env.update(dict(OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1",
                    VECLIB_MAXIMUM_THREADS="1", MPLBACKEND="Agg",
                    BGS_LOWMEM_SPOT_BLOCK="256", BGS_LOWMEM_SKIP_PLOTS="1"))
    old_env = env.copy()
    if args.baseline_path:
        old_env["PYTHONPATH"] = args.baseline_path
    gss = str(root / "gss.feather")
    cache = str(root / "cache")
    common = ["run", "--gss", gss, "--coords", str(root / "coords.tsv"),
              "--burden", str(root / "burden.tsv"), "--gene-col", "gene_symbol",
              "--p-col", "Pvalue_Burden", "--min-genes", "20", "--permutations", "0",
              "--seed", "42", "--score-mode", "logp", "--duplicate-policy", "bonferroni"]
    meta = {"platform": platform.platform(), "python": platform.python_version(),
            "arrow": pa.__version__, "numpy": np.__version__, "spots": args.spots,
            "genes": args.genes, "domains": args.domains, "block_size": args.block_size,
            "projection_size": args.projection_size}
    runs = {}
    differences = {}
    runs["lowmem_1_1_0"] = measure(args.baseline_python, "bgsmap_lowmem",
                                   [*common, "--out", str(root / "lowmem_1_1_0")],
                                   root, "lowmem_1_1_0", old_env)
    fast = [*common, "--engine", "fast", "--skip-plots",
            "--projection-size", str(args.projection_size), "--block-size", str(args.block_size)]
    for use_cache in (False, True):
        if use_cache:
            runs["build_cache"] = measure(sys.executable, "bgsmap", [
                "cache-gss", "--gss", gss, "--out", cache,
                "--projection-size", str(args.projection_size),
            ], root, "build_cache", env)
        for workers in dict.fromkeys((1, args.threads)):
            label = f"{'cached' if use_cache else 'direct'}_{workers}"
            extra = ["--gss-cache", cache] if use_cache else []
            runs[label] = measure(sys.executable, "bgsmap",
                                  [*fast, *extra, "--threads", str(workers),
                                   "--out", str(root / label)], root, label, env)
            differences[label] = compare(root, "lowmem_1_1_0", label)
    if args.check_legacy:
        measure(args.baseline_python, "bgsmap", [*common, "--out", str(root / "legacy")],
                root, "legacy", old_env)
        differences["legacy"] = compare(root, "lowmem_1_1_0", "legacy")
    (root / "metrics.json").write_text(
        json.dumps({"fixture": meta, "runs": runs, "max_absolute_difference": differences}, indent=2))


if __name__ == "__main__":
    main()
