# Statistical methods in v0.2.0

This implementation tests association of gene-level burden intensity with a gsMap GSS-derived *spatial* profile. It is not gsMap SNP-level spatial LDSC.

Let `G` be the number of matched genes. For each gene, standardize its GSS across spots: `X[g,s] = (GSS[g,s] - mean_spots(GSS[g]))/sd_spots(GSS[g])`. Set burden strength `b[g] = min(-log10(max(burden_p[g],1e-300)), p_cap)` (binary alternative is available). Remove within-stratum means of `b` and normalize the remaining weights to standard deviation one, obtaining `w[g]`.

For each spot, define `T_s = (1/G) Σ_g X[g,s]*w[g]`. Domain mean inference replaces `X[g,s]` by the mean across spots of a domain. Competitive domain inference replaces it by mean inside the domain minus mean outside it. All are one-sided positive enrichment tests.

The null holds the full spatial GSS and the multiset of burden weights *fixed*, permuting burden weights over genes within each matching stratum `h` (conditionally exchangeable genes). The finite-population null variance of a profile `a` is exactly

`Var_perm(T_a) = G^-2 * Σ_h [(Σ_{i∈h}(w_i - mean_h(w))²)/(n_h-1)] * Σ_{i∈h}(a_i - mean_h(a))²`.

The primary **analytic Gaussian-tail approximation** uses `Z_a = T_a / sqrt(Var_perm(T_a))` and `p_analytic = Phi(-Z_a)`. It is computed in log survival space with SciPy `log_ndtr(-Z)`, and both `p_analytic` and `logp_analytic = -log10(p_analytic)` are written. All spot maps render `logp_analytic`. All competitive domain maps render `logp_competitive_analytic`. The approximation is not quantized by a finite number of permutations; P values below the float64 minimum underflow to zero in the `p_...` column but **not** in `logp_...`.

Within each analysis/run, `q_bh_analytic` is BH-adjusted and `p_holm_analytic` is Holm-adjusted; the correction is done in log space for extreme P values. Corresponding competitive columns refer to only the *competitive domain tests within that same run*. No correction pooling all 108 phenotypes is performed. BH under spatial dependency requires caution; Holm provides FWER control *if each marginal analytic P value is valid*, not a cure for invalid normal tails.

`--permutations 0` (default) computes only analytic P. When `--permutations >=19`, synchronized gene-label shuffles produce diagnostics (empirical P, empirical BH, synchronized maxT permutation FWER) written **only** to `spot_permutation_diagnostics.tsv.gz` and `region_permutation_diagnostics.tsv`. They are never used for ranking, plotting, or replacing `p_analytic` in main tables.

**Statistical validity:** Exact conditional second moments do **not** imply exact Gaussian tails. Normal approximation may fail for very sparse/highly influential burden signals, nonexchangeable strata, extreme GSS, or gene-gene dependence. `effective_gene_weight_count = (sum_g w_g²)² / sum_g w_g⁴` is a heuristic that flags concentrated burdens (<30) but is not a proof of asymptotic validity. Before interpreting extremely small P values, assess Type-I and FWER calibration with independently replicated null simulations, matched gene detection-opportunity covariates, competitive-null simulations and sensible alternative models. Exact permutation diagnostics remain valuable at moderate P values, despite their limited tail resolution. Large `logp` should not be sold as a calibrated P without these validations.

**One-domain edge case:** competitive mean inside minus mean outside is undefined when the whole section is one domain; the code returns a zero statistic and P=1 as a non-significant placeholder. Do not interpret it as a competitive test.
