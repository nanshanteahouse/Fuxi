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
    _compute_marker_coverage(adata, cluster_key, per_cell_scores, ...) -> float
"""

import numpy as np
from scipy.spatial import ConvexHull
from typing import Literal

import logging

logger = logging.getLogger(__name__)

DEFAULT_MULTI_METRIC_WEIGHTS = {
    'silhouette': 0.3,
    'stability': 0.3,
    'marker_coverage': 0.4,
}


def select_best_params(results_summary, method="pareto_elbow", best_resolution=None, best_n_neighbors=0, multi_metric_weights=None):
    """Select the best (n_neighbors, resolution) from a grid search summary.

    Parameters
    ----------
    results_summary : list of dict
        Each dict must have keys: 'n_neighbors', 'resolution',
        'n_clusters', 'silhouette_score'.
    method : str or None
        "pareto_elbow"  — Pareto frontier + normalized elbow detection
        "silhouette"    — Pick max silhouette score
        "multi_metric"  — Composite scoring: silhouette + stability + marker coverage
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
        Defaults to :data:`DEFAULT_MULTI_METRIC_WEIGHTS` (silhouette=0.3,
        stability=0.3, marker_coverage=0.4) when None.

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
        return _select_multi_metric(valid, weights=multi_metric_weights)
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


def _compute_marker_coverage(adata, cluster_key, per_cell_scores, ratio_threshold=1.5):
    """Compute fraction of clusters with a clear marker-gene match.

    For each cluster, computes the mean per-cell score for each cell type.
    A cluster is "matched" if the top cell type's mean score is at least
    *ratio_threshold* × the second-best score.

    Parameters
    ----------
    adata : AnnData
    cluster_key : str
        Key in ``adata.obs`` for cluster labels.
    per_cell_scores : dict[str, np.ndarray]
        Cell-type → per-cell score array (shape: n_cells,).  Pre-computed
        by the caller (e.g. via ``sc.tl.score_genes()``).
    ratio_threshold : float
        Minimum ratio of top score to second-best for a match.

    Returns
    -------
    float
        Fraction of clusters with a clear marker match (0.0–1.0).
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
    n_matched = 0

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

        if top_score <= 0.0:
            continue

        if len(mean_scores) == 1:
            n_matched += 1
        else:
            second_score = mean_scores[1][1]
            if top_score >= ratio_threshold * second_score:
                n_matched += 1

    return n_matched / n_clusters


def _select_multi_metric(valid, weights=None):
    """Select best clustering via composite multi-metric scoring.

    Reads precomputed ``silhouette_score``, ``stability_score``, and
    ``marker_coverage`` from each entry dict, normalises each metric to
    [0, 1] across all entries, then computes a weighted composite.
    Returns argmax.

    Parameters
    ----------
    valid : list[dict]
        Pre-filtered entries (must have valid ``silhouette_score``).
    weights : dict[str, float] | None
        Metric weights.  Defaults to :data:`DEFAULT_MULTI_METRIC_WEIGHTS`.
        Degrades to silhouette+stability (0.5/0.5) when no entry has
        ``marker_coverage``.

    Returns
    -------
    tuple[int, float, str, str]
        (n_neighbors, resolution, 'multi_metric', reason_str)
    """
    n = len(valid)

    # ── Gather raw scores ──
    sil_scores = np.array([r['silhouette_score'] for r in valid])
    stab_scores = np.array([r.get('stability_score', 0.0) for r in valid])

    has_marker_coverage = any('marker_coverage' in r for r in valid)
    mc_scores: np.ndarray = np.zeros(n)
    if has_marker_coverage:
        mc_scores = np.array([r.get('marker_coverage', 0.0) for r in valid])

    # ── Determine initial weights ──
    if weights is None:
        if has_marker_coverage:
            active_weights = dict(DEFAULT_MULTI_METRIC_WEIGHTS)
        else:
            active_weights = {'silhouette': 0.5, 'stability': 0.5}
    else:
        active_weights = dict(weights)

    # Remove marker_coverage from weights if entries lack it
    if not has_marker_coverage:
        active_weights.pop('marker_coverage', None)

    # ── Low-variance guard (on raw scores) ──
    metrics_raw = {
        'silhouette': sil_scores,
        'stability': stab_scores,
    }
    if has_marker_coverage:
        metrics_raw['marker_coverage'] = mc_scores

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
    if 'marker_coverage' in active_weights and has_marker_coverage:
        norm['marker_coverage'] = _normalize(mc_scores)

    # ── Composite score ──
    composite = np.zeros(n)
    for metric_name, w in active_weights.items():
        composite += w * norm[metric_name]

    best_idx = int(np.argmax(composite))
    best = valid[best_idx]

    # ── Build reason string ──
    sil_val = best['silhouette_score']
    stab_val = best.get('stability_score', 0.0)
    mc_val = best.get('marker_coverage', 0.0) if has_marker_coverage else 0.0
    k = best['n_clusters']

    sil_norm_val = norm.get('silhouette', np.zeros(n))[best_idx]
    stab_norm_val = norm.get('stability', np.zeros(n))[best_idx]

    reason = (
        f"composite={composite[best_idx]:.4f} "
        f"sil={sil_val:.4f}(n={sil_norm_val:.2f}) "
        f"stab={stab_val:.3f}(n={stab_norm_val:.2f}) "
        f"marker_cov={mc_val:.3f} k={k}"
    )

    return (
        best['n_neighbors'],
        best['resolution'],
        'multi_metric',
        reason,
    )


def select_best_umap_params(adata, best_n, min_dist_grid, spread_grid, method, CFG, use_rep, log):
    """Sweep min_dist × spread on the best (n_neighbors) neighbor graph,
    or use manual fallback.

    The KNN graph is rebuilt once (or reused if already present).
    Selection `method` follows the same logic as cluster_selection_method:

        "convex_hull"  — auto-sweep, pick largest convex-hull area (default)
        None            — manual: use CFG.umap_min_dist / CFG.umap_spread directly

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

    use_paga = getattr(CFG, 'umap_paga_init', False)

    # ── Manual mode ──
    if method is None:
        md = getattr(CFG, 'umap_min_dist', 0.3)
        sp = getattr(CFG, 'umap_spread', 1.0)
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
        md = getattr(CFG, 'umap_min_dist', 0.3)
        sp = getattr(CFG, 'umap_spread', 1.0)
        log.info("UMAP params (convex_hull, empty grid → fallback): min_dist=%.2f, spread=%.1f",
                 md, sp)
        return md, sp, "convex_hull", []

    # ── Ensure neighbor graph exists for best_n ──
    log.info("Building KNN graph (n_neighbors=%d) for UMAP parameter sweep...", best_n)
    try:
        sc.pp.neighbors(
            adata, n_neighbors=best_n,
            n_pcs=CFG.n_pcs_use, use_rep=use_rep,
            random_state=CFG.random_seed,
        )
    except Exception as e:
        log.error("KNN graph build failed for UMAP sweep: %s", e)
        return (
            getattr(CFG, 'umap_min_dist', 0.3),
            getattr(CFG, 'umap_spread', 1.0),
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
                           random_state=CFG.random_seed)
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
