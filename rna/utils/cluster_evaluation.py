#!/usr/bin/env python3
"""
Cluster parameter selection methods for the Fuxi pipeline.

Provides objective, quantitative selection of the best (n_neighbors, resolution)
from a grid search summary, replacing naive silhouette-score-max selection.

Exports:
    select_best_params(results_summary, method, best_resolution=None)
        -> (best_n, best_r, method_label, reason_str)
"""

import numpy as np
from scipy.spatial import ConvexHull


def select_best_params(results_summary, method="pareto_elbow", best_resolution=None, best_n_neighbors=0):
    """Select the best (n_neighbors, resolution) from a grid search summary.

    Parameters
    ----------
    results_summary : list of dict
        Each dict must have keys: 'n_neighbors', 'resolution',
        'n_clusters', 'silhouette_score'.
    method : str or None
        "pareto_elbow"  — Pareto frontier + normalized elbow detection
        "silhouette"    — Pick max silhouette score
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
    else:
        raise ValueError(
            f"Unknown cluster_selection_method: {method!r}. "
            f"Valid options: 'pareto_elbow', 'silhouette', None"
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

    # -- Pareto frontier --
    n = len(pts)
    is_pareto = np.ones(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # j dominates i
            if (pts[j, 0] <= pts[i, 0] and pts[j, 1] >= pts[i, 1]
                    and (pts[j, 0] < pts[i, 0] or pts[j, 1] > pts[i, 1])):
                is_pareto[i] = False
                break

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
