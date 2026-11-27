# Changelog

## 1.2.0

- Integrate the fast, bounded-memory Feather V2 analytic engine. Run it with
  `bgsmap-lowmem run` or `bgsmap run --engine fast`.
- Separate wide Feather projections from spot compute blocks, with a bounded
  pool of spot workers and separate read/compute phase logging.
- Add `cache-gss` to reuse spot-major GSS data across burden traits.
- Share burden weighting, analytic P, BH/Holm and result-table code with the
  in-memory engine. The default `bgsmap run` continues to use that engine.
- Add HESTA Slurm instructions, numeric regression tests and a benchmark
  against the `bgsmap-lowmem 1.1.0` extension wheel.

## 1.0.0 (2026-11-20)

- First stable release of the `bgsmap` package and CLI, formerly `gsmap-burden`.
- Source archives and wheels use the same version as the CLI and run metadata.

## 0.2.0

- Introduce the `bgsmap` package and CLI name and timestamped progress logging to stderr.
- Introduce primary **continuous analytic Gaussian-tail P values**, using the exact finite-population variance of the conditional within-stratum gene-label permutation distribution (new statistic-level inference; see docs/statistics.md).
- Write `p_analytic`, `logp_analytic`, `z_analytic`, analytic BH/Holm and corresponding competitive columns into the primary spot/domain result tables.
- All spatial / competitive region HTML, PDF and report tables use analytic P or `-log10(P)` (computed without float64 underflow). **No plots use `p_empirical`.** Top descriptive gene domain uses analytic competitive P.
- Optional permutations (`--permutations >=19`) go into standalone `*_permutation_diagnostics.tsv*` files, no longer limit or alter primary P. Default `--permutations 0`.
- Explicit error if attempting v0.2.0 `plot` on older empirical-only result tables; rerun the analysis first.
- Retain per-run (per phenotype) multiple testing scope; no cross-phenotype correction is added.
- Include exact-variance enumeration tests, extreme-tail log P tests, and HTML/PDF regression tests.
- Gaussian tail approximation requires empirical calibration in rare/sparse burden settings; no guarantee of exact small-P calibration is claimed.

## 0.1.1 (2026-11-12)

- Recover `obs_names` via AnnData in backed mode when fast H5AD decoding produces zero GSS spot matches.
- Preserve exact spot-ID alignment; never align GSS to coordinates by row order.
- Add informative unmatched-ID diagnostics, alignment match logging, and `inspect --check-spots`.
- Inspect very-wide Feather GSS using Arrow schema metadata without materializing numerical values.
- No changes to the burden statistic, permutation algorithm, or existing CLI options.
