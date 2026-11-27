"""Spatial gene-label enrichment: continuous normal-approximation P values.

Conditional on a *fixed* GSS matrix and burden weights, the null randomizes
burden gene labels within exchangeable strata. The finite-population variance
is exact. Gaussian tail probabilities based on that variance are approximations,
not exact permutation P values. Simulations/QQ diagnostics are essential for
strongly sparse, correlated or heterogeneous gene signals.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
import math

import numpy as np
import pandas as pd
from scipy.special import log_ndtr

LOG = logging.getLogger(__name__)


@dataclass
class EnrichmentResult:
    spot: pd.DataFrame
    region: pd.DataFrame
    gene: pd.DataFrame
    meta: dict
    spot_permutation: pd.DataFrame | None = None
    region_permutation: pd.DataFrame | None = None


def bh(p):
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('BH requires a 1D array of valid finite P values')
    if len(p) == 0:
        return p.copy()
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p), dtype=float)
    out[order] = np.minimum(q, 1.0)
    return out


def _log10_adjust(logp, method):
    """Return -log10(adjusted P), without underflow even at P < 1e-308.

    logp is -log10(raw P). Holm controls FWER under arbitrary dependence;
    BH is a separate FDR summary and requires its usual dependence conditions.
    """
    v = np.asarray(logp, dtype=np.float64)
    if v.ndim != 1 or not np.isfinite(v).all() or (v < -1e-10).any():
        raise ValueError('Invalid -log10 P values')
    n = len(v)
    if n == 0:
        return v.copy()
    order = np.argsort(-v, kind='stable')
    seq = v[order]
    if method == 'bh':
        candidates = seq - np.log10(n / np.arange(1, n+1))
        adjusted = np.maximum.accumulate(candidates[::-1])[::-1]
    elif method == 'holm':
        candidates = seq - np.log10(np.arange(n, 0, -1))
        adjusted = np.minimum.accumulate(candidates)
    else:
        raise ValueError('Unknown adjustment: ' + method)
    ans = np.empty(n, dtype=np.float64)
    ans[order] = np.maximum(adjusted, 0.0)
    return ans


def _p_columns(table, stem, observed, variance):
    """Add continuous one-sided Gaussian-tail estimates and adjustments."""
    obs = np.asarray(observed, dtype=np.float64)
    var = np.asarray(variance, dtype=np.float64)
    if not np.isfinite(var).all() or (var < 0).any():
        raise ValueError('Non-finite / negative exact permutation variance')
    zero = var <= np.finfo(float).tiny
    z = np.zeros_like(obs)
    np.divide(obs, np.sqrt(var), out=z, where=~zero)
    # A zero-variance statistic has no positive null enrichment support.
    # Return p=1 when statistic is approximately zero, or fail loudly if
    # degenerate math is inconsistent with a nonzero observation.
    if np.any(zero & (np.abs(obs) > 1e-6)):
        raise ValueError('Observed statistic is nonzero but null variance is zero')
    logp = np.maximum(-log_ndtr(-z) / math.log(10), 0)
    logp[zero] = 0
    logq_bh = _log10_adjust(logp, 'bh')
    logp_holm = _log10_adjust(logp, 'holm')
    table[f'z_{stem}'] = z
    table[f'null_sd_{stem}'] = np.sqrt(var)
    table[f'p_{stem}'] = np.power(10., -logp)
    table[f'logp_{stem}'] = logp
    table[f'q_bh_{stem}'] = np.power(10., -logq_bh)
    table[f'logq_bh_{stem}'] = logq_bh
    table[f'p_holm_{stem}'] = np.power(10., -logp_holm)
    table[f'logp_holm_{stem}'] = logp_holm
    return table


def make_strata(table, columns=None, bins=4, min_size=12):
    """Build covariate-matched strata; reject nonexchangeable tiny groups."""
    if not columns:
        return np.zeros(len(table), dtype=np.int32), {'matched':False, 'columns':[]}
    if not set(columns).issubset(table.columns):
        raise ValueError(f'Missing matching columns: {set(columns)-set(table.columns)}')
    labels = []
    for col in columns:
        value = table[col]
        if value.isna().any():
            raise ValueError(f'Matching column {col} contains missing values')
        numeric = pd.to_numeric(value, errors='coerce')
        if numeric.notna().all():
            if numeric.nunique() > bins:
                coded = pd.qcut(numeric.rank(method='first'), q=bins, labels=False)
                labels.append(coded.astype(str).to_numpy())
            else:
                labels.append(value.astype(str).to_numpy())
        else:
            labels.append(value.astype(str).to_numpy())
    combined = np.array(['|'.join(fields) for fields in zip(*labels)], dtype=str)
    names, codes = np.unique(combined, return_inverse=True)
    counts = np.bincount(codes)
    if counts.min() < min_size:
        raise ValueError(f'Matching bins too small: smallest has {counts.min()} genes (<{min_size}); reduce --match-bins or matching columns')
    return codes.astype(np.int32), {'matched':True,'columns':list(columns),
                                     'strata_sizes':counts.tolist(), 'n_strata':len(names)}


def _permutation_variance(profiles, w, strata_indices, block_rows=256):
    """Exact Var_perm[sum(a_i*w_perm_i)/G], batched to cap RAM.

    For stratum h, Var(sum_i a_i W_perm_i) =
       sum_i(w_i-mean(w_h))**2 / (m_h-1)
       * sum_i(a_i-mean(a_h))**2.
    No independent-gene or independent-spot approximation is made for
    these second moments. Normality of the tail remains an assumption.
    """
    factors = _variance_factors(w, strata_indices)
    out = np.empty(len(profiles), dtype=np.float64)
    for start in range(0, len(profiles), block_rows):
        chunk = profiles[start:start+block_rows]
        out[start:start+len(chunk)] = _profile_variance(chunk, factors)
    return out


def _variance_factors(w, strata_indices):
    factors = []
    for ix in strata_indices:
        if len(ix) > 1:
            wi = np.asarray(w[ix], dtype=np.float64)
            wi -= wi.mean()
            factor = np.dot(wi, wi) / (len(ix) - 1)
            if factor:
                factors.append((ix, factor))
    return factors


def _profile_variance(profiles, factors):
    n_rows, n_genes = profiles.shape
    out = np.zeros(n_rows, dtype=np.float64)
    for ix, factor in factors:
        a = profiles if len(ix) == n_genes else profiles[:, ix]
        sums = a.sum(axis=1, dtype=np.float64)
        ss = np.einsum('ij,ij->i', a, a, dtype=np.float64)
        out += factor * np.maximum(ss - sums*sums/len(ix), 0.0)
    return out / float(n_genes)**2


def _burden_weights(burden, score_mode, binary_threshold, p_cap, matching_cols, match_bins):
    if score_mode not in ('logp', 'binary'):
        raise ValueError('Invalid score_mode')
    if score_mode == 'logp':
        w_raw = np.minimum(-np.log10(np.maximum(burden.p.to_numpy(dtype=float), 1e-300)), p_cap)
    else:
        w_raw = (burden.p.to_numpy(dtype=float) <= binary_threshold).astype(float)
    if not np.isfinite(w_raw).all() or np.std(w_raw) <= 1e-10:
        raise ValueError('Burden signal has no variation. Select a different score mode or p threshold')
    strata, matching = make_strata(burden, columns=matching_cols, bins=match_bins)
    strata_indices = [np.flatnonzero(strata == label) for label in np.unique(strata)]
    w = w_raw.astype(np.float64).copy()
    for ix in strata_indices:
        w[ix] -= w[ix].mean()
    if np.std(w) < 1e-10:
        raise ValueError('No within-stratum burden weight variance; relax matching strata')
    w /= np.std(w)
    return w_raw, w, strata_indices, matching


def run_enrichment(
    gss: np.ndarray, genes, spots, burden: pd.DataFrame, coord: pd.DataFrame,
    permutations=0, seed=42, chunk_size=64, p_cap=30.0,
    score_mode='logp', binary_threshold=0.05,
    matching_cols=None, match_bins=4,
) -> EnrichmentResult:
    """Main continuous P: normal tail with exact stratified randomization variance.

    Optional synchronized permutations provide *separate* empirical diagnostics
    and maxT, never drive plots or replace the primary continuous P values.
    """
    if permutations < 0 or (0 < permutations < 19) or chunk_size < 1:
        raise ValueError('Use 0 permutations (analytic only), or >=19 for diagnostics')
    if not np.isfinite(gss).all():
        raise ValueError('GSS has non-finite values')
    n_genes, n_spots = gss.shape
    if n_genes < 20 or n_spots < 5:
        raise ValueError('Too few aligned genes (<20) or spots (<5)')
    if len(genes) != n_genes or len(spots) != n_spots or len(burden) != n_genes:
        raise ValueError('GSS/genes/spots/burden mismatch')
    w_raw, w, strata_indices, matching = _burden_weights(
        burden, score_mode, binary_threshold, p_cap, matching_cols, match_bins
    )
    LOG.info("Standardizing GSS and burden weights (%d gene strata, score mode %s)",
             len(strata_indices), score_mode)
    x = gss.astype(np.float32, copy=True)
    x -= x.mean(axis=1, keepdims=True)
    std = x.std(axis=1)
    if (std < 1e-9).any():
        raise ValueError('Zero-variance GSS genes must be filtered before enrichment')
    x /= std[:,None]
    xs = np.ascontiguousarray(x.T, dtype=np.float32)
    w32 = w.astype(np.float32)
    obs = (xs @ w32) / n_genes
    labels = coord.loc[list(spots), 'domain'].to_numpy(dtype=str)
    domains = pd.Index(pd.unique(labels), name='domain')
    dom_ix = [np.flatnonzero(labels == d) for d in domains]
    dom_profile = np.vstack([xs[ix].mean(axis=0) for ix in dom_ix]).astype(np.float32)
    dom_obs = (dom_profile @ w32) / n_genes
    overall = xs.mean(axis=0)
    comp_profile = np.vstack([
        dom_profile[j] - ((n_spots*overall - len(ix)*dom_profile[j])/(n_spots-len(ix)))
        if len(ix) < n_spots else np.zeros(n_genes, dtype=np.float32)
        for j, ix in enumerate(dom_ix)
    ]).astype(np.float32)
    comp_obs = (comp_profile @ w32) / n_genes

    LOG.info("Computing analytic enrichment for %d spots and %d domains", n_spots, len(domains))
    spot_var = _permutation_variance(xs, w, strata_indices)
    dom_var = _permutation_variance(dom_profile, w, strata_indices)
    comp_var = _permutation_variance(comp_profile, w, strata_indices)

    out_sp, out_reg = _analytic_tables(
        coord, spots, domains, [len(ix) for ix in dom_ix],
        obs, spot_var, dom_obs, dom_var, comp_obs, comp_var,
    )

    spot_perm = None
    region_perm = None
    if permutations:
        LOG.info("Running %d diagnostic permutations (batch size %d, seed %d)",
                 permutations, chunk_size, seed)
        progress_interval = max(chunk_size, (permutations + 9) // 10)
        next_progress = chunk_size
        count_sp = np.zeros(n_spots, dtype=np.int64)
        count_reg = np.zeros(len(domains), dtype=np.int64)
        count_comp = np.zeros(len(domains), dtype=np.int64)
        max_sp = np.empty(permutations, dtype=np.float32)
        max_reg = np.empty(permutations, dtype=np.float32)
        max_comp = np.empty(permutations, dtype=np.float32)
        rng = np.random.default_rng(seed)
        for start in range(0, permutations, chunk_size):
            end = min(start + chunk_size, permutations)
            wp = np.broadcast_to(w32, (end-start,n_genes)).copy()
            for ix in strata_indices:
                if len(ix) > 1:
                    for k in range(end-start):
                        wp[k,ix] = w32[ix][rng.permutation(len(ix))]
            wp = wp.T.copy()
            null = (xs @ wp) / n_genes
            null_reg = (dom_profile @ wp) / n_genes
            null_comp = (comp_profile @ wp) / n_genes
            count_sp += (null >= obs[:,None]-1e-9).sum(axis=1)
            count_reg += (null_reg >= dom_obs[:,None]-1e-9).sum(axis=1)
            count_comp += (null_comp >= comp_obs[:,None]-1e-9).sum(axis=1)
            max_sp[start:end] = null.max(axis=0)
            max_reg[start:end] = null_reg.max(axis=0)
            max_comp[start:end] = null_comp.max(axis=0)
            if end >= next_progress or end == permutations:
                LOG.info("Permutation diagnostics: %d/%d completed", end, permutations)
                next_progress = end + progress_interval
        p_count = lambda count: (1+count)/(1+permutations)
        spot_p = p_count(count_sp)
        reg_p = p_count(count_reg)
        comp_p = p_count(count_comp)
        spot_perm = pd.DataFrame({'spot':out_sp.spot, 'domain':out_sp.domain,
                                  'p_empirical':spot_p,
                                  'q_bh_empirical':bh(spot_p),
                                  'p_fwer_permutation':p_count((max_sp[None,:]>=obs[:,None]-1e-9).sum(axis=1))})
        region_perm = pd.DataFrame({'domain':domains.astype(str),
                                    'p_empirical':reg_p,
                                    'q_bh_empirical':bh(reg_p),
                                    'p_fwer_permutation':p_count((max_reg[None,:]>=dom_obs[:,None]-1e-9).sum(axis=1)),
                                    'p_competitive_empirical':comp_p,
                                    'q_competitive_empirical':bh(comp_p),
                                    'p_competitive_fwer_permutation':p_count((max_comp[None,:]>=comp_obs[:,None]-1e-9).sum(axis=1))})
    top_domain = str(out_reg.iloc[0].domain)
    out_gene = _gene_contributions(burden, w_raw, w, top_domain, domains, comp_profile)
    meta = _result_metadata(
        len(spots), len(domains), w, permutations, seed, score_mode, p_cap,
        binary_threshold, matching, top_domain,
    )
    return EnrichmentResult(out_sp, out_reg, out_gene, meta, spot_perm, region_perm)


def _analytic_tables(coords, spots, domains, counts, obs, spot_var, dom_obs, dom_var,
                     comp_obs, comp_var):
    out_sp = coords.loc[list(spots), ['spot','x','y','domain']].reset_index(drop=True).copy()
    out_sp['score'] = obs
    _p_columns(out_sp, 'analytic', obs, spot_var)
    # Legacy 'logp' aliases are deliberately omitted: scripts must explicitly
    # select the inferential target they render.
    out_reg = pd.DataFrame({'domain':domains.astype(str),
                            'n_spots':counts,
                            'mean_score':dom_obs,
                            'competitive_score':comp_obs})
    _p_columns(out_reg, 'analytic', dom_obs, dom_var)
    _p_columns(out_reg, 'competitive_analytic', comp_obs, comp_var)
    out_reg = out_reg.sort_values('p_competitive_analytic', kind='stable').reset_index(drop=True)

    return out_sp, out_reg


def _gene_contributions(burden, w_raw, w, top_domain, domains, comp_profile):
    n_genes = len(w)
    top_j = domains.get_loc(top_domain)
    contrib = w * comp_profile[top_j] / n_genes
    out_gene = burden[['gene','p']].copy()
    out_gene['burden_weight'] = w_raw
    out_gene['standardized_weight'] = w
    out_gene['top_region'] = top_domain
    out_gene['region_vs_rest_gss_z'] = comp_profile[top_j]
    out_gene['contribution'] = contrib
    out_gene = out_gene.sort_values('contribution', ascending=False).reset_index(drop=True)
    return out_gene


def _result_metadata(n_spots, n_domains, w, permutations, seed, score_mode, p_cap,
                     binary_threshold, matching, top_domain):
    n_genes = len(w)
    effective_n = (np.dot(w, w)**2 / np.dot(w*w, w*w))
    meta = {'n_genes':n_genes,'n_spots':n_spots,'n_domains':n_domains,
            'permutations':permutations,'random_seed':seed,
            'score_mode':score_mode,'p_cap':p_cap,'binary_threshold':binary_threshold,
            'matching':matching,'one_sided':True,
            'primary_p_column':'p_analytic',
            'primary_competitive_p_column':'p_competitive_analytic',
            'primary_spatial_plot_column':'logp_analytic',
            'primary_region_plot_column':'logp_competitive_analytic',
            'primary_method':'normal approximation using exact conditional within-stratum permutation variance',
            'pvalue_type':'approximate, not exact permutation',
            'effective_gene_weight_count':round(float(effective_n),2),
            'normal_approximation_warning':bool(effective_n < 30),
            'permutation_diagnostics_saved':bool(permutations),
            'top_region_for_descriptive_gene_contributions':top_domain,
            'warning':'Normal-tail P values avoid finite Monte Carlo resolution, but can be miscalibrated for sparse/heavy-tailed burden weights, gene dependence, or nonexchangeable strata; test null calibration before publication.'}
    return meta
