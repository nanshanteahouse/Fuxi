#!/usr/bin/env python3
"""
Cluster parameter selection methods for the Fuxi pipeline.

Provides objective, quantitative selection of the best (n_neighbors, resolution)
from a grid search summary, supporting four methods: 'pareto_elbow',
'silhouette', 'multi_metric' (default), and None (manual).
from a grid search summary, replacing naive silhouette-score-max selection.

Exports:
    select_best_params(results_summary, method, best_resolution=None, best_n_neighbors=0, multi_metric_weights=None)
        -> (best_n, best_r, method_label, reason_str)
    _compute_stability(adata, resolution, ...) -> float
    _compute_cluster_coherence(adata, cluster_key, per_cell_scores, ...) -> float
"""

import numpy as np
from scipy.spatial import ConvexHull
from typing import Literal

import logging

logger = logging.getLogger(__name__)

DEFAULT_MULTI_METRIC_WEIGHTS = {
    'silhouette': 0.2,
    'stability': 0.2,
    'cluster_coherence': 0.3,
    'splitting_gain': 0.2,
    'kb_annotatable_rate': 0.1,
}


def select_best_params(results_summary, method="pareto_elbow", best_resolution=None, best_n_neighbors=0, multi_metric_weights=None, log=None):
    """Select the best (n_neighbors, resolution) from a grid search summary.

    Parameters
    ----------
    results_summary : list of dict
        Each dict must have keys: 'n_neighbors', 'resolution',
        'n_clusters', 'silhouette_score'.
    method : str or None
        "pareto_elbow"  — Pareto frontier + normalized elbow detection
        "silhouette"    — Pick max silhouette score
        "multi_metric"  — Composite scoring: silhouette + stability + cluster coherence
        None            — Manual via best_resolution + best_n_neighbors
                           (falls back to max silhouette within matching
                           resolution if n_neighbors=0, then globally if
                           no match found)
    best_resolution : float or None
        Only used when method is None.  If present in the grid, picks
        the best silhouette among entries matching that resolution.
    best_n_neighbors : int
        Only used when method is None.  If > 0, requires an exact match
        on both resolution and n_neighbors.  Default 0 = auto-pick best
        silhouette at the given resolution.
    multi_metric_weights : dict[str, float] | None
        Only used when method is "multi_metric".  Custom metric weights.
        Defaults to :data:`DEFAULT_MULTI_METRIC_WEIGHTS` (silhouette=0.2,
        stability=0.2, cluster_coherence=0.3, splitting_gain=0.2,
        kb_annotatable_rate=0.1) when None.

    Returns
    -------
    best_n : int
        Best n_neighbors value.
    best_r : float
        Best resolution value.
    method_label : str
        Human-readable method name for logging.
    reason : str
        One-line diagnostic for the logger.
    """
    # -- Filter invalid / missing silhouette scores --
    valid = [r for r in results_summary
             if r.get('silhouette_score') is not None
             and not (isinstance(r['silhouette_score'], float)
                      and np.isnan(r['silhouette_score']))]

    if not valid:
        raise ValueError("No valid silhouette scores in results_summary")

    # ── Method dispatch ──
    if method is None:
        return _select_manual(valid, best_resolution, best_n_neighbors)
    elif method == "pareto_elbow":
        return _select_pareto_elbow(valid)
    elif method == "silhouette":
        return _select_max_silhouette(valid)
    elif method == "multi_metric":
        return _select_multi_metric(valid, weights=multi_metric_weights, log=log)
    else:
        raise ValueError(
            f"Unknown cluster_selection_method: {method!r}. "
            f"Valid options: 'pareto_elbow', 'silhouette', 'multi_metric', None"
        )


# ═══════════════════════════════════════════════════════════════════════
#  Internal selection strategies
# ═══════════════════════════════════════════════════════════════════════

def _select_max_silhouette(valid):
    """Pick the combination with the highest silhouette score."""
    best = max(valid, key=lambda r: r['silhouette_score'])
    return (
        best['n_neighbors'],
        best['resolution'],
        "silhouette",
        f"silhouette={best['silhouette_score']:.4f} k={best['n_clusters']}",
    )


def _select_manual(valid, best_resolution, best_n_neighbors=0):
    """Manual selection via best_resolution + optional best_n_neighbors.

    - If best_n_neighbors > 0: require exact match on both.
    - If best_n_neighbors == 0: filter by resolution only, pick best silhouette.
    - Falls back to max silhouette if the requested combination is not in grid.
    """
    if best_resolution is not None:
        matching = [r for r in valid if r['resolution'] == best_resolution]
        if matching:
            if best_n_neighbors and best_n_neighbors > 0:
                # Exact combination requested
                exact = [r for r in matching if r['n_neighbors'] == best_n_neighbors]
                if exact:
                    best = exact[0]
                    return (
                        best['n_neighbors'],
                        best['resolution'],
                        "manual",
                        f"n_neighbors={best_n_neighbors}, resolution={best['resolution']:.1f} "
                        f"(configured) silhouette={best['silhouette_score']:.4f} "
                        f"k={best['n_clusters']}",
                    )
                # n_neighbors not found at this resolution → fall through to
                # pick best silhouette at the given resolution
            best = max(matching, key=lambda r: r['silhouette_score'])
            return (
                best['n_neighbors'],
                best['resolution'],
                "manual",
                f"resolution={best['resolution']:.1f} (configured) "
                f"silhouette={best['silhouette_score']:.4f} k={best['n_clusters']}",
            )
        # best_resolution set but not in grid → fall through to auto
    # Fallback: max silhouette
    best = max(valid, key=lambda r: r['silhouette_score'])
    return (
        best['n_neighbors'],
        best['resolution'],
        "silhouette",
        f"best_resolution={best_resolution} not in grid, "
        f"fallback silhouette={best['silhouette_score']:.4f} k={best['n_clusters']}",
    )


def _select_pareto_elbow(valid):
    """Pareto frontier + normalized elbow detection.

    Algorithm:
      1. Compute the Pareto frontier in (n_clusters, silhouette_score) space.
         A point i is dominated if there exists j with:
           n_clusters[j] <= n_clusters[i] AND silhouette[j] >= silhouette[i]
           AND at least one strict inequality.
      2. Normalize both axes to [0, 1].
      3. Pick the Pareto point closest to the ideal point (k_min=0, s_max=1)
         in normalized space.

    Returns (best_n, best_r, method_label, reason_str).
    """
    # Build (k, ss) array
    pts = np.array([(r['n_clusters'], r['silhouette_score']) for r in valid])

    # -- Pareto frontier (O(n log n) sort+scan) --
    n = len(pts)
    # Sort by n_clusters asc, silhouette desc: for equal clusters the best silhouette comes first.
    order = np.lexsort((-pts[:, 1], pts[:, 0]))
    pts_sorted = pts[order]

    is_pareto_sorted = np.zeros(n, dtype=bool)
    best_s = -np.inf
    best_k = -1
    for i in range(n):
        k, s = pts_sorted[i]
        if s > best_s or (s == best_s and k == best_k):
            # Not dominated: either better silhouette than any earlier point,
            # or identical to an earlier Pareto point (no domination between equals)
            is_pareto_sorted[i] = True
            best_s = s
            best_k = k

    # Map Pareto flags back to original (unsorted) order
    is_pareto = np.zeros(n, dtype=bool)
    is_pareto[order] = is_pareto_sorted

    pareto_idx = np.where(is_pareto)[0]
    pareto_k = pts[pareto_idx, 0]
    pareto_s = pts[pareto_idx, 1]

    # Sort by k ascending
    sort_order = np.argsort(pareto_k)
    pareto_k = pareto_k[sort_order]
    pareto_s = pareto_s[sort_order]
    pareto_idx = pareto_idx[sort_order]

    # -- Handle single Pareto point --
    if len(pareto_k) == 1:
        best = valid[pareto_idx[0]]
        return (
            best['n_neighbors'],
            best['resolution'],
            "pareto_elbow",
            f"single_pareto_point silhouette={best['silhouette_score']:.4f} "
            f"k={best['n_clusters']}",
        )

    # -- Normalize --
    eps = 1e-10
    k_norm = (pareto_k - pareto_k.min()) / (pareto_k.max() - pareto_k.min() + eps)
    s_norm = (pareto_s - pareto_s.min()) / (pareto_s.max() - pareto_s.min() + eps)

    # Distance to ideal point (k_norm=0, s_norm=1)
    dist = np.sqrt((1.0 - s_norm)**2 + k_norm**2)
    elbow_idx = np.argmin(dist)

    best = valid[pareto_idx[elbow_idx]]

    # Supplementary: compute ΔSS/Δk for the transition into this point
    # (from the previous Pareto point, if any)
    delta_note = ""
    if elbow_idx > 0:
        dk = pareto_k[elbow_idx] - pareto_k[elbow_idx - 1]
        ds = pareto_s[elbow_idx] - pareto_s[elbow_idx - 1]
        ratio = ds / dk if dk > 0 else 0
        delta_note = f" ΔSS/Δk={ratio:.6f}"

    return (
        best['n_neighbors'],
        best['resolution'],
        "pareto_elbow",
        f"dist_to_ideal={dist[elbow_idx]:.4f} "
        f"silhouette={best['silhouette_score']:.4f} k={best['n_clusters']}"
        f"{delta_note}",
    )


# ═══════════════════════════════════════════════════════════════════════
#  Multi-metric scoring helpers
# ═══════════════════════════════════════════════════════════════════════


def _compute_stability(adata, resolution, leiden_flavor: Literal['leidenalg', 'igraph'] = 'igraph', n_seeds=5, base_seed=42):
    """Compute cross-seed clustering stability via pairwise ARI.

    Re-runs Leiden clustering with *n_seeds* different random seeds at the
    given *resolution*, then computes pairwise adjusted Rand index between
    all label sets.  Returns mean ARI as a stability score.

    Parameters
    ----------
    adata : AnnData
    resolution : float
        Resolution parameter for Leiden clustering.
    leiden_flavor : str
        Flavour of leiden algorithm (default 'igraph').
    n_seeds : int
        Number of different seeds to run.
    base_seed : int
        Starting seed value.

    Returns
    -------
    float
        Mean pairwise ARI.  1.0 when *n_seeds* <= 1 or all clusterings
        are identical.
    """
    if n_seeds <= 1:
        return 1.0

    from sklearn.metrics import adjusted_rand_score
    import scanpy as sc

    seeds = range(base_seed, base_seed + n_seeds)
    temp_keys = [f'_temp_stab_{i}' for i in range(n_seeds)]

    label_sets = []
    for i, seed in enumerate(seeds):
        key = temp_keys[i]
        try:
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=key,
                random_state=seed,
                flavor=leiden_flavor,
                n_iterations=2,
                directed=False,
            )
            label_sets.append(adata.obs[key].values)
        except Exception:
            pass
        finally:
            if key in adata.obs.columns:
                del adata.obs[key]

    # Final cleanup: remove any remaining temp columns
    for col in list(adata.obs.columns):
        if col.startswith('_temp_stab_'):
            del adata.obs[col]

    n_runs = len(label_sets)
    if n_runs <= 1:
        return 1.0

    # Compute pairwise ARI
    aris = []
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            ari = adjusted_rand_score(label_sets[i], label_sets[j])
            aris.append(ari)

    return float(np.mean(aris))


def _compute_cluster_coherence(adata, cluster_key, per_cell_scores, dominance_threshold=1.5, min_expression=0.05):
    """Per-cluster marker coherence metric that peaks at intermediate resolution.

    For each cluster label in adata.obs[cluster_key]:
      - Compute mean of per_cell_scores per cell type
      - Find best_score (max) and second_best_score
      - A cluster is "coherent" if:
        (a) best_score > min_expression (not noise)
        (b) best_score / max(second_best, 1e-10) > dominance_threshold
      - This ensures the cluster has ONE clearly dominant cell type

    coherence = n_coherent / n_total_clusters

    Parameters
    ----------
    adata : AnnData
    cluster_key : str
        Key in ``adata.obs`` for cluster labels.
    per_cell_scores : dict[str, np.ndarray]
        Cell-type → per-cell score array (shape: n_cells,).  Pre-computed
        by the caller (e.g. via ``sc.tl.score_genes()``).
    dominance_threshold : float
        Minimum ratio of best score to second-best for coherence.
    min_expression : float
        Minimum mean expression for a cluster to be considered (not noise).

    Returns
    -------
    float
        Fraction of clusters that are coherent (0.0-1.0).
        1.0 if *per_cell_scores* is empty or has no valid entries.
    """
    if not per_cell_scores:
        return 1.0

    # Check for any valid (non-None) entries
    has_valid = any(v is not None for v in per_cell_scores.values())
    if not has_valid:
        return 1.0

    cluster_labels = adata.obs[cluster_key].values
    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)

    if n_clusters <= 1:
        return 1.0

    cell_types = list(per_cell_scores.keys())
    n_cells = len(cluster_labels)
    n_coherent = 0

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        mean_scores = []
        for ct in cell_types:
            scores = per_cell_scores[ct]
            if scores is not None and len(scores) == n_cells:
                mean_scores.append((ct, float(np.mean(scores[mask]))))

        if not mean_scores:
            continue

        # Sort by mean score descending
        mean_scores.sort(key=lambda x: x[1], reverse=True)
        top_score = mean_scores[0][1]

        if top_score <= min_expression:
            continue

        if len(mean_scores) == 1:
            n_coherent += 1
        else:
            second_score = mean_scores[1][1]
            if top_score / max(second_score, 1e-10) > dominance_threshold:
                n_coherent += 1

    return n_coherent / n_clusters


def _compute_splitting_gain(valid_by_resolution: list[dict]) -> dict[float, float]:
    """Compute splitting gain per resolution.

    Splitting gain measures how many new clusters are created per unit
    resolution increase. For each resolution r_i (except the lowest):
        splitting_gain(r_i) = max(0, (n_clusters(r_i) - n_clusters(r_{i-1})) / (r_i - r_{i-1}))

    Parameters
    ----------
    valid_by_resolution : list of dict
        List of results_summary entries sorted by resolution ascending,
        all from the SAME n_neighbors group.

    Returns
    -------
    dict[float, float]
        {resolution_value: splitting_gain_value}
        Empty dict if fewer than 2 resolutions.
    """
    if len(valid_by_resolution) < 2:
        return {}

    # Ensure sorted by resolution
    sorted_entries = sorted(valid_by_resolution, key=lambda e: e.get('resolution', 0.0))

    gains = {}
    for i in range(1, len(sorted_entries)):
        prev = sorted_entries[i - 1]
        curr = sorted_entries[i]
        r_prev = prev['resolution']
        r_curr = curr['resolution']
        k_prev = prev['n_clusters']
        k_curr = curr['n_clusters']

        delta_k = k_curr - k_prev
        delta_r = r_curr - r_prev

        if delta_r > 0:
            gain = max(0.0, delta_k / delta_r)
        else:
            gain = 0.0

        gains[r_curr] = gain

    return gains

def _select_multi_metric(valid, weights=None, log=None):
    """Select best clustering via composite multi-metric scoring.

    Reads precomputed ``silhouette_score``, ``stability_score``,
    ``cluster_coherence``, and ``splitting_gain`` from each entry dict, normalises each metric to
    [0, 1] across all entries, then computes a weighted composite.
    Returns argmax.

    Parameters
    ----------
    valid : list[dict]
        Pre-filtered entries (must have valid ``silhouette_score``).
    weights : dict[str, float] | None
        Metric weights.  Defaults to :data:`DEFAULT_MULTI_METRIC_WEIGHTS`.
        Degrades from 4-metric to 3-metric when no entry has
        ``splitting_gain``, and to silhouette+stability when no entry
        has ``cluster_coherence``.

    Returns
    -------
    tuple[int, float, str, str]
        (n_neighbors, resolution, 'multi_metric', reason_str)
    """
    n = len(valid)

    # ── Gather raw scores ──
    sil_scores = np.array([r['silhouette_score'] for r in valid])
    stab_scores = np.array([r.get('stability_score', 0.0) for r in valid])

    has_coherence = any('cluster_coherence' in r for r in valid)
    coh_scores: np.ndarray = np.zeros(n)
    if has_coherence:
        coh_scores = np.array([r.get('cluster_coherence', 0.0) for r in valid])

    has_splitting_gain = any('splitting_gain' in r for r in valid)
    split_scores: np.ndarray = np.zeros(n)
    if has_splitting_gain:
        split_scores = np.array([r.get('splitting_gain', 0.0) for r in valid])
    has_kb_rate = any('kb_annotatable_rate' in r for r in valid)
    kb_scores: np.ndarray = np.zeros(n)
    if has_kb_rate:
        kb_scores = np.array([r.get('kb_annotatable_rate', 0.0) for r in valid])

    # ── Determine initial weights ──
    if weights is None:
        if has_coherence and has_splitting_gain and has_kb_rate:
            active_weights = dict(DEFAULT_MULTI_METRIC_WEIGHTS)
        elif has_coherence and has_splitting_gain:
            # Degrade to 4-metric (no kb_annotatable_rate)
            active_weights = {'silhouette': 0.2, 'stability': 0.2, 'cluster_coherence': 0.35, 'splitting_gain': 0.25}
        elif has_coherence:
            # Degrade to 3-metric (no splitting_gain)
            active_weights = {'silhouette': 0.25, 'stability': 0.25, 'cluster_coherence': 0.5}
        else:
            active_weights = {'silhouette': 0.5, 'stability': 0.5}
    else:
        active_weights = dict(weights)

    # Remove cluster_coherence from weights if entries lack it
    if not has_coherence:
        active_weights.pop('cluster_coherence', None)

    # Remove splitting_gain from weights if entries lack it
    if not has_splitting_gain:
        active_weights.pop('splitting_gain', None)

    # Remove kb_annotatable_rate from weights if entries lack it
    if not has_kb_rate:
        active_weights.pop('kb_annotatable_rate', None)

    # -- Coherence mismatch auto-degrade: if all entries have cluster_coherence < 0.1 --
    if has_coherence and float(np.max(coh_scores)) < 0.1:
        logger.warning(
            "Max cluster_coherence=%.4f < 0.1 across all entries — marker_dict may be mismatched. "
            "Degrading to silhouette+stability only.",
            float(np.max(coh_scores)),
        )
        has_coherence = False
        active_weights = {'silhouette': 0.5, 'stability': 0.5}

    # ── Low-variance guard (on raw scores) ──
    metrics_raw = {
        'silhouette': sil_scores,
        'stability': stab_scores,
    }
    if has_coherence:
        metrics_raw['cluster_coherence'] = coh_scores
    if has_splitting_gain:
        metrics_raw['splitting_gain'] = split_scores
    if has_kb_rate:
        metrics_raw['kb_annotatable_rate'] = kb_scores

    for metric_name in list(active_weights.keys()):
        scores = metrics_raw.get(metric_name)
        if scores is None:
            continue
        score_range = float(np.max(scores) - np.min(scores))
        if score_range < 0.01:
            logger.warning(
                "%s variance < 0.01 (range=%.4f) — disabling metric",
                metric_name, score_range,
            )
            del active_weights[metric_name]

    if not active_weights:
        logger.warning("All metrics dropped — falling back to silhouette only")
        active_weights = {'silhouette': 1.0}

    # ── Renormalise weights to sum=1.0 ──
    total_w = sum(active_weights.values())
    if total_w > 0:
        for k in active_weights:
            active_weights[k] /= total_w

    # ── Normalise each metric to [0, 1] ──
    def _normalize(values):
        vmin, vmax = float(np.min(values)), float(np.max(values))
        if vmax - vmin < 1e-10:
            return np.ones_like(values)
        return (values - vmin) / (vmax - vmin + 1e-10)

    norm = {}
    if 'silhouette' in active_weights:
        norm['silhouette'] = _normalize(sil_scores)
    if 'stability' in active_weights:
        norm['stability'] = _normalize(stab_scores)
    if 'cluster_coherence' in active_weights and has_coherence:
        norm['cluster_coherence'] = _normalize(coh_scores)
    if 'splitting_gain' in active_weights and has_splitting_gain:
        norm['splitting_gain'] = _normalize(split_scores)
    if 'kb_annotatable_rate' in active_weights and has_kb_rate:
        norm['kb_annotatable_rate'] = _normalize(kb_scores)

    # ── Composite score ──
    composite = np.zeros(n)
    for metric_name, w in active_weights.items():
        composite += w * norm[metric_name]

    # ── 3-tier resolution recommendation logging ──
    try:
        entries_by_resolution = sorted(
            [(valid[i], composite[i]) for i in range(n)],
            key=lambda x: x[0]['resolution']
        )
        comp_max = max(composite)
        stab_max = max(stab_scores) if 'stability' in active_weights else None

        # Coarse: lowest resolution whose composite > 0.7 of max composite
        coarse_entry = None
        for entry, comp in entries_by_resolution:
            if comp > 0.7 * comp_max:
                coarse_entry = entry
                break

        # Balanced: the best (current behavior)
        balanced_entry = valid[int(np.argmax(composite))]

        # Fine: highest resolution whose stability > 0.85 of max stability
        fine_entry = None
        if stab_max is not None and stab_max > 0:
            for entry, _ in reversed(entries_by_resolution):
                stab = entry.get('stability_score', 0.0)
                if stab > 0.85 * stab_max:
                    fine_entry = entry
                    break

        (log or logger).info(
            "[multi_metric 3-tier] coarse: r=%.2f (k=%d) / balanced: r=%.2f (k=%d) / fine: r=%.2f (k=%d)",
            coarse_entry['resolution'] if coarse_entry else float('nan'),
            coarse_entry['n_clusters'] if coarse_entry else 0,
            balanced_entry['resolution'],
            balanced_entry['n_clusters'],
            fine_entry['resolution'] if fine_entry else float('nan'),
            fine_entry['n_clusters'] if fine_entry else 0,
        )
    except Exception as exc:
        logger.debug("3-tier computation skipped: %s", exc)

    best_idx = int(np.argmax(composite))
    best = valid[best_idx]

    # ── Build reason string ──
    sil_val = best['silhouette_score']
    stab_val = best.get('stability_score', 0.0)
    coh_val = best.get('cluster_coherence', 0.0) if has_coherence else 0.0
    split_val = best.get('splitting_gain', 0.0) if has_splitting_gain else 0.0
    kb_val = best.get('kb_annotatable_rate', 0.0) if has_kb_rate else 0.0
    k = best['n_clusters']

    sil_norm_val = norm.get('silhouette', np.zeros(n))[best_idx]
    stab_norm_val = norm.get('stability', np.zeros(n))[best_idx]

    reason = (
        f"composite={composite[best_idx]:.4f} "
        f"sil={sil_val:.4f}(n={sil_norm_val:.2f}) "
        f"stab={stab_val:.3f}(n={stab_norm_val:.2f}) "
        f"coherence={coh_val:.3f} "
        f"split_gain={split_val:.3f} kb_rate={kb_val:.3f} k={k}"
    )

    return (
        best['n_neighbors'],
        best['resolution'],
        'multi_metric',
        reason,
    )



def _detect_granularity(results_summary: list[dict], cv_threshold: float = 0.05, min_clusters: int = 10) -> str:
    """Determine whether the data is tissue-level or subtype-level.

    Two-path architecture for cluster granularity detection. This function
    analyses grid-search results to decide between two downstream strategies:

    - "tissue": multiple distinct cell types → use full DE pipeline
    - "subtype": FACS-enriched / similar cells → use gated DE (Wave 4)

    The decision is based on the coefficient of variation (CV) of silhouette
    scores within the median n_neighbors group. Low CV means silhouette
    scores are flat across resolutions — typical of subtype data where all
    partitions are similarly mediocre. Combined with a low maximum cluster
    count, this signals subtype-level resolution.

    Algorithm
    ---------
    1. Group entries by n_neighbors, pick the median-size group.
    2. From that group, sort by resolution, collect silhouette scores
       and n_clusters.
    3. Compute CV = std(silhouette_scores) / mean(silhouette_scores).
    4. If CV < cv_threshold AND max_n_clusters < min_clusters: "subtype".
    5. Otherwise: "tissue" (conservative default).

    Edge cases (all return "tissue"):
    - Empty results_summary
    - Single entry
    - Mean silhouette = 0 (division guard)

    Parameters
    ----------
    results_summary : list of dict
        Grid-search summary entries with keys: n_neighbors, resolution,
        n_clusters, silhouette_score.
    cv_threshold : float
        CV threshold below which silhouette flatness signals subtype data.
    min_clusters : int
        Maximum n_clusters below which (combined with low CV) signals
        subtype data.

    Returns
    -------
    str
        "tissue" or "subtype"
    """
    # Edge case: empty or single entry → conservative default
    if not results_summary or len(results_summary) <= 1:
        return "tissue"

    # 1. Group entries by n_neighbors
    groups: dict[int, list[dict]] = {}
    for entry in results_summary:
        if 'n_neighbors' not in entry:
            continue
        nn = entry['n_neighbors']
        groups.setdefault(nn, []).append(entry)

    if not groups:
        return "tissue"

    # Pick median-size group (by number of entries)
    group_sizes = sorted(groups.items(), key=lambda kv: len(kv[1]))
    median_idx = len(group_sizes) // 2
    _median_nn, median_group = group_sizes[median_idx]

    # 2. Collect silhouette scores and n_clusters from the median group
    #    Sort by resolution for deterministic ordering
    median_group_sorted = sorted(median_group, key=lambda e: e.get('resolution', 0.0))

    silhouette_values = []
    n_clusters_values = []
    for entry in median_group_sorted:
        if 'silhouette_score' not in entry or entry['silhouette_score'] is None:
            continue
        if 'n_clusters' not in entry or entry['n_clusters'] is None:
            continue
        silhouette_values.append(entry['silhouette_score'])
        n_clusters_values.append(entry['n_clusters'])

    if not silhouette_values:
        return "tissue"

    # 3. Compute CV of silhouette scores
    sil_arr = np.array(silhouette_values, dtype=float)
    mean_sil = float(np.mean(sil_arr))
    if mean_sil == 0.0:
        logger.debug(
            "_detect_granularity: mean silhouette is 0 → conservative 'tissue' "
            f"(n_neighbors={_median_nn})"
        )
        return "tissue"

    std_sil = float(np.std(sil_arr))
    cv = std_sil / mean_sil

    # 4. & 5. Decision
    max_n_clusters = max(n_clusters_values) if n_clusters_values else 0

    if cv < cv_threshold and max_n_clusters < min_clusters:
        logger.debug(
            f"_detect_granularity: CV={cv:.5f} < {cv_threshold} AND "
            f"max_n_clusters={max_n_clusters} < {min_clusters} → 'subtype' "
            f"(n_neighbors={_median_nn}, n_entries={len(median_group)})"
        )
        return "subtype"
    else:
        logger.debug(
            f"_detect_granularity: CV={cv:.5f} (threshold={cv_threshold}), "
            f"max_n_clusters={max_n_clusters} (threshold={min_clusters}) → 'tissue' "
            f"(n_neighbors={_median_nn}, n_entries={len(median_group)})"
        )
        return "tissue"

def _select_de_gated(valid, adata, de_gate_threshold=25):
    """Select best resolution using DE-gated criterion for subtype-level data.

    For subtype-level data where silhouette scores are flat across resolutions,
    this method selects the highest resolution that maintains a minimum number of
    differentially expressed genes between every cluster pair.

    Algorithm follows Shekhar 2016 merge.clusters.DE pattern (inverted):
    select highest resolution where pairwise DE >= threshold.

    Parameters
    ----------
    valid : list[dict]
        Grid-search entries with keys: n_clusters, resolution, cluster_key.
        Must be unique by resolution.
    adata : AnnData
        Annotated data matrix.
    de_gate_threshold : int
        Minimum number of DE genes required between every cluster pair.
        Default 25.

    Returns
    -------
    tuple[int, float, str, str]
        (n_clusters, resolution, cluster_key, reason_str)
    """
    import scanpy as sc

    # Edge case: single entry -> return as-is
    if len(valid) <= 1:
        entry = valid[0]
        return (
            entry["n_clusters"],
            entry["resolution"],
            entry["cluster_key"],
            "de_gated(single_entry, min_pairwise_de=N/A)",
        )

    # Sort by resolution ascending for deterministic iteration
    sorted_entries = sorted(valid, key=lambda e: e.get("resolution", 0.0))

    # Collect (entry, min_pairwise_de) for each resolution
    candidates = []  # list of (entry, min_pairwise_de)

    for entry in sorted_entries:
        cluster_key = entry["cluster_key"]

        # Check if rank_genes_groups already computed for this groupby key
        existing_rg = adata.uns.get("rank_genes_groups")
        recompute = True
        if existing_rg is not None:
            existing_params = existing_rg.get("params", {})
            if existing_params.get("groupby") == cluster_key:
                recompute = False

        if recompute:
            try:
                sc.tl.rank_genes_groups(
                    adata,
                    groupby=cluster_key,
                    method="wilcoxon",
                    n_genes=50,
                    use_raw=True,
                )
            except Exception as e:
                logger.warning(
                    "rank_genes_groups failed for cluster_key=%s (resolution=%.2f): %s",
                    cluster_key, entry["resolution"], e,
                )
                continue

        # Extract min pairwise DE count: for each group, count genes with
        # padj < 0.05 AND log2FC > 1.0, then take the minimum across groups.
        try:
            rg = adata.uns["rank_genes_groups"]
            pvals_adj = rg["pvals_adj"]
            logfoldchanges = rg["logfoldchanges"]
            group_names = pvals_adj.dtype.names

            if group_names is None or len(group_names) == 0:
                logger.warning(
                    "No cluster groups in rank_genes_groups output for cluster_key=%s",
                    cluster_key,
                )
                continue

            de_counts = []
            for group in group_names:
                padj = pvals_adj[group]
                lfc = logfoldchanges[group]
                n_de = int(np.sum((padj < 0.05) & (lfc > 1.0)))
                de_counts.append(n_de)

            min_pairwise_de = int(min(de_counts))
        except Exception as e:
            logger.warning(
                "Failed to extract DE counts for cluster_key=%s: %s",
                cluster_key, e,
            )
            continue

        candidates.append((entry, min_pairwise_de))

    if not candidates:
        # All entries failed -> fallback to first entry
        logger.warning(
            "DE-gated selection: all rank_genes_groups calls failed -- "
            "fallback to first entry"
        )
        entry = sorted_entries[0]
        return (
            entry["n_clusters"],
            entry["resolution"],
            entry["cluster_key"],
            "de_gated(all_failed, min_pairwise_de=N/A)",
        )

    # Select entry with highest n_clusters where min_pairwise_de >= threshold
    candidates_by_n = sorted(
        candidates, key=lambda c: c[0]["n_clusters"], reverse=True,
    )
    best_entry = None
    best_min_de = -1
    for entry, min_de in candidates_by_n:
        if min_de >= de_gate_threshold:
            best_entry = entry
            best_min_de = min_de
            break

    # If no entry meets threshold, fallback to lowest resolution (most conservative)
    if best_entry is None:
        best_entry, best_min_de = candidates[0]
        logger.info(
            "DE-gated selection: best_resolution=%.2f, min_pairwise_de=%d "
            "(fallback: no entry met threshold=%d)",
            best_entry["resolution"], best_min_de, de_gate_threshold,
        )
    else:
        logger.info(
            "DE-gated selection: best_resolution=%.2f, min_pairwise_de=%d",
            best_entry["resolution"], best_min_de,
        )

    return (
        best_entry["n_clusters"],
        best_entry["resolution"],
        best_entry["cluster_key"],
        f"de_gated(min_pairwise_de={best_min_de})",
    )


def select_best_umap_params(adata, best_n, min_dist_grid, spread_grid, method, CFG, use_rep, log):
    """Sweep min_dist × spread on the best (n_neighbors) neighbor graph,
    or use manual fallback.

    The KNN graph is rebuilt once (or reused if already present).
    Selection `method` follows the same logic as cluster_selection_method:

        "convex_hull"  — auto-sweep, pick largest convex-hull area (default)
        None            — manual: use CFG.clustering.umap_min_dist / CFG.clustering.umap_spread directly

    Parameters
    ----------
    adata : AnnData
        Must already have PCA representation (use_rep).
    best_n : int
        Best n_neighbors value (from select_best_params).
    min_dist_grid : list of float or None
        Values to sweep in "convex_hull" mode.
    spread_grid : list of float or None
        Values to sweep in "convex_hull" mode.
    method : str or None
        "convex_hull" | None
    CFG : Config
    use_rep : str
        Key in adata.obsm for PCA (e.g. 'X_pca_harmony' or 'X_pca').
    log : logging.Logger

    Returns
    -------
    best_min_dist : float
    best_spread : float
    method_label : str
        Human-readable method name for logging.
    results : list of dict
        Each dict: {min_dist, spread, convex_hull_area}
        Empty list if sweep was skipped.
    """
    import scanpy as sc
    import numpy as np
    import pandas as pd

    use_paga = getattr(CFG.clustering, 'umap_paga_init', False)

    # ── Manual mode ──
    if method is None:
        md = getattr(CFG.clustering, 'umap_min_dist', 0.3)
        sp = getattr(CFG.clustering, 'umap_spread', 1.0)
        log.info("UMAP params (manual): min_dist=%.2f, spread=%.1f", md, sp)
        return md, sp, "manual", []

    if method != "convex_hull":
        raise ValueError(
            f"Unknown umap_selection_method: {method!r}. "
            f"Valid options: 'convex_hull', None"
        )

    # ── Auto-sweep: convex_hull ──
    do_sweep = True
    if min_dist_grid is None or spread_grid is None:
        do_sweep = False
    elif len(min_dist_grid) <= 1 and len(spread_grid) <= 1:
        do_sweep = False

    if not do_sweep:
        md = getattr(CFG.clustering, 'umap_min_dist', 0.3)
        sp = getattr(CFG.clustering, 'umap_spread', 1.0)
        log.info("UMAP params (convex_hull, empty grid → fallback): min_dist=%.2f, spread=%.1f",
                 md, sp)
        return md, sp, "convex_hull", []

    # ── Ensure neighbor graph exists for best_n ──
    log.info("Building KNN graph (n_neighbors=%d) for UMAP parameter sweep...", best_n)
    try:
        sc.pp.neighbors(
            adata, n_neighbors=best_n,
            n_pcs=CFG.pca.n_pcs_use, use_rep=use_rep,
            random_state=CFG.execution.random_seed,
        )
    except Exception as e:
        log.error("KNN graph build failed for UMAP sweep: %s", e)
        return (
            getattr(CFG.clustering, 'umap_min_dist', 0.3),
            getattr(CFG.clustering, 'umap_spread', 1.0),
            "convex_hull",
            [],
        )

    # ── Sweep ──
    results = []
    best_area = -1.0
    best_md = min_dist_grid[0]
    best_sp = spread_grid[0]

    for md in min_dist_grid:
        for sp in spread_grid:
            try:
                sc.tl.umap(adata, min_dist=md, spread=sp,
                           init_pos='paga' if use_paga else 'spectral',
                           random_state=CFG.execution.random_seed)
                coords = adata.obsm['X_umap']
                hull = ConvexHull(coords)
                area = float(hull.volume)  # 2D → area
                results.append({
                    'min_dist': md,
                    'spread': sp,
                    'convex_hull_area': area,
                })
                log.info("  min_dist=%.2f, spread=%.1f → convex_hull_area=%.2f",
                         md, sp, area)
                if area > best_area:
                    best_area = area
                    best_md = md
                    best_sp = sp
            except Exception as e:
                log.warning("  UMAP failed (min_dist=%.2f, spread=%.1f): %s", md, sp, e)
                results.append({
                    'min_dist': md,
                    'spread': sp,
                    'convex_hull_area': None,
                })

    log.info("Best UMAP params (convex_hull): min_dist=%.2f, spread=%.1f (area=%.2f)",
             best_md, best_sp, best_area)
    return best_md, best_sp, "convex_hull", results
