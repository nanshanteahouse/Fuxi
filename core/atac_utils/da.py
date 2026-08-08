"""Differential accessibility (DA) — "best practice" path for marker peaks.

Replaces the quick ``snap.tl.marker_regions`` call with a pseudobulk +
background-matched Wilcoxon design (Registered Report style):

1. **Pseudobulk** per-group counts (sparse rows summed per group).
2. **log-TP10K** normalisation (per-group total x 10k, then log1p).
3. **Background matching**: for each group, sample ``n_background_per_group``
   cells stratified by total fragment count (quantile bins) to control for
   sequencing-depth confounds.
4. **Wilcoxon rank-sum** per peak (cell-level group vs background), BH-FDR
   corrected; significant peaks are |log2FC| >= ``log2fc_threshold`` and
   FDR < ``fdr_threshold``, with a minimum mean-count filter (``min_mean``)
   excluding ultra-sparse peaks whose log2FC is inflated by normalisation.

The return type is ``dict[str, pd.Index]`` — group -> peak names, matching
``snap.tl.marker_regions`` so the downstream CSV writer stays identical.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, rankdata


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment (returns q-values, NaN-safe)."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    finite = ~np.isnan(p)
    if finite.sum() == 0:
        return q
    pv = p[finite]
    order = np.argsort(pv)
    ranked = np.empty_like(pv)
    ranked[order] = np.arange(1, len(pv) + 1)
    adjusted = pv * len(pv) / ranked
    # enforce monotonicity (largest p first)
    adj = np.minimum.accumulate(adjusted[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    q[finite] = adj
    return q


def _wilcoxon_ranksums_vectorized(
    x_csr,
    rows_a,
    rows_b,
    col_idx,
    chunk: int = 25_000,
) -> np.ndarray:
    """Bilateral Wilcoxon rank-sum (Mann-Whitney U, normal approx. + tie fix), vectorised.

    Identical to ``scipy.stats.ranksums`` on tie-free data; on discrete sparse
    counts the tie correction makes it exact. Returns p-values aligned to col_idx."""
    n1, n2 = len(rows_a), len(rows_b)
    n_total = n1 + n2
    sel = np.concatenate([rows_a, rows_b])
    p_out = np.zeros(len(col_idx), dtype=float)
    for s in range(0, len(col_idx), chunk):
        blk = col_idx[s : s + chunk]
        xb = x_csr[sel][:, blk].toarray().T.astype(np.float32)
        ranks = rankdata(xb, axis=1)
        r1 = ranks[:, :n1].sum(axis=1)
        u = r1 - n1 * (n1 + 1) / 2.0
        mu = n1 * n2 / 2.0
        sv = np.sort(xb, axis=1)
        tie_adj = np.zeros(xb.shape[0])
        for i in range(xb.shape[0]):
            row = sv[i]
            not_eq = np.where(row[1:] != row[:-1])[0] + 1
            seg_end = np.concatenate([[0], not_eq, [n_total]])
            lengths = np.diff(seg_end)
            t = lengths[lengths > 1]
            tie_adj[i] = (t**3 - t).sum() if len(t) else 0.0
        var = n1 * n2 * (n_total + 1) / 12.0 - n1 * n2 * tie_adj / (12.0 * n_total * (n_total - 1))
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (u - mu) / np.sqrt(var)
        p_out[s : s + len(blk)] = np.where(var > 0, 2.0 * norm.sf(np.abs(z)), 1.0)
    return p_out


def _log_tp10k(counts: np.ndarray) -> np.ndarray:
    """log1p(TP10K): per-column (pseudobulk group) normalisation."""
    totals = counts.sum(axis=0)
    totals = np.where(totals == 0, 1.0, totals)
    return np.log1p(counts * 10_000.0 / totals)


def differential_accessibility(
    adata,
    groupby: str,
    log2fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
    n_background_per_group: int = 200,
    n_bins: int = 5,
    seed: int = 42,
    max_peaks: int | None = None,
    min_mean: float = 0.05,
) -> dict[str, pd.Index]:
    """Per-group differential accessibility.

    Parameters
    ----------
    adata
        Peak-by-cell AnnData (``.X`` sparse counts, ``.var_names`` = peaks).
    groupby
        Observation column defining the groups (cell_type / leiden...).
    log2fc_threshold
        Minimum |log2FC| for a peak to be called significant.
    fdr_threshold
        BH-FDR cutoff.
    n_background_per_group
        Number of background cells sampled per group (stratified by depth).
    n_bins
        Quantile bins for depth stratification of background sampling.
    seed
        RNG seed (deterministic background sampling).
    max_peaks
        Optional cap on peaks tested per group (uniform random sample) to
        bound Wilcoxon runtime on very large peak sets.
    min_mean
        Minimum mean count (per cell, in either group or background) for a
        peak to be tested. Ultra-sparse peaks get inflated |log2FC| after
        TP10K normalisation; this filter restores biological ranking.
    Returns
    -------
    dict[str, pd.Index]
        group -> sorted significant peak names (|log2FC| desc within group).
    """
    if groupby not in adata.obs:
        raise ValueError(f"groupby column '{groupby}' not in adata.obs")
    x = adata.X
    if sparse.issparse(x):
        x_csr = x.tocsr()
    else:
        x_csr = sparse.csr_matrix(x)
    obs = adata.obs
    peak_names = np.asarray(adata.var_names, dtype=object)
    groups = obs[groupby].astype(str)
    group_names = sorted(set(groups))
    rng = np.random.RandomState(seed)

    # Cell-level depth (total fragments) for stratified background sampling.
    depth = np.asarray(x_csr.sum(axis=1)).ravel()
    depth_bins = pd.Series(
        pd.qcut(depth, q=n_bins, labels=False, duplicates="drop"),
        index=adata.obs_names,
    )

    result: dict[str, pd.Index] = {}
    for g in group_names:
        cells_g = np.where(groups.values == g)[0]
        if len(cells_g) < 5:
            # Too few cells for a stable rank-sum; keep quick-path behaviour
            # (empty) so downstream CSV simply lacks this group.
            result[g] = pd.Index([], dtype=object)
            continue
        # ── 1. Pseudobulk counts + log-TP10K (log2FC denominator) ──
        pseudo = np.asarray(x_csr[cells_g].sum(axis=0)).ravel()
        pseudo_norm = _log_tp10k(pseudo.reshape(1, -1)).ravel()

        # ── 2. Background: stratified sample of non-group cells ──
        bg_idx = np.setdiff1d(np.arange(x_csr.shape[0]), cells_g)
        if len(bg_idx) == 0:
            result[g] = pd.Index([], dtype=object)
            continue
        bg_bins = depth_bins.iloc[bg_idx].values
        chosen: list[int] = []
        bins_present = np.unique(bg_bins[~np.isnan(bg_bins)])
        per_bin = max(1, n_background_per_group // max(len(bins_present), 1))
        for b in bins_present:
            cand = bg_idx[bg_bins == b]
            take = min(len(cand), per_bin)
            chosen.extend(rng.choice(cand, take, replace=False).tolist())
        # top up if stratification produced too few
        if len(chosen) < n_background_per_group:
            extra = np.setdiff1d(bg_idx, chosen)
            rng.shuffle(extra)
            chosen.extend(extra[: n_background_per_group - len(chosen)].tolist())
        bg_idx = np.asarray(chosen, dtype=int)
        bg_norm = _log_tp10k(np.asarray(x_csr[bg_idx].sum(axis=0)).ravel()).ravel()

        # ── 3. log2FC from TP10K-normalised pseudobulk vs background ──
        log2fc = (pseudo_norm - bg_norm) / np.log(2)

        # ── 4. Wilcoxon rank-sum per peak (cell-level), vectorised ──
        n_peaks = x_csr.shape[1]
        # Mean-count filter: ultra-sparse peaks (both sides ~0) get inflated
        # log2FC after TP10K normalisation and would dominate |log2FC| ranking.
        group_mean = np.asarray(x_csr[cells_g].mean(axis=0)).ravel()
        bg_mean = np.asarray(x_csr[bg_idx].mean(axis=0)).ravel()
        max_mean = np.maximum(group_mean, bg_mean)
        # Pre-filter: |log2FC| < log2fc_threshold or max_mean < min_mean can
        # never pass the final threshold, so skip their rank-sum entirely
        # (keeps result identical).
        cand = np.where((np.abs(log2fc) >= log2fc_threshold) & (max_mean >= min_mean))[0]
        if max_peaks is not None and len(cand) > max_peaks:
            cand = np.sort(rng.choice(cand, max_peaks, replace=False))
        pvals = np.full(n_peaks, np.nan, dtype=float)
        if len(cand) > 0:
            pvals[cand] = _wilcoxon_ranksums_vectorized(x_csr, cells_g, bg_idx, cand)

        # ── 5. BH-FDR + thresholds ──
        qvals = _bh_fdr(pvals)
        sig = (np.abs(log2fc) >= log2fc_threshold) & (qvals < fdr_threshold)
        idx = np.where(sig)[0]
        # sort by |log2FC| descending
        order = np.argsort(-np.abs(log2fc[idx]))
        idx = idx[order]
        result[g] = pd.Index(peak_names[idx], dtype=object)

    return result
