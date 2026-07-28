"""Selection strategies for grid-search parameter evaluation.

Submodule of ``core/cluster/evaluation`` containing the top-level dispatch
functions and all selection-helper strategies.
"""

import logging

import numpy as np
from scipy import stats
from scipy.spatial import ConvexHull

from .enrichment import DEFAULT_MULTI_METRIC_WEIGHTS

logger = logging.getLogger(__name__)


def select_best_params(
    results_summary,
    method="pareto_elbow",
    best_resolution=None,
    best_n_neighbors=0,
    multi_metric_weights=None,
    log=None,
):
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
    valid = [
        r
        for r in results_summary
        if r.get("silhouette_score") is not None
        and not (isinstance(r["silhouette_score"], float) and np.isnan(r["silhouette_score"]))
    ]

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
    best = max(valid, key=lambda r: r["silhouette_score"])
    return (
        best["n_neighbors"],
        best["resolution"],
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
        matching = [r for r in valid if r["resolution"] == best_resolution]
        if matching:
            if best_n_neighbors and best_n_neighbors > 0:
                # Exact combination requested
                exact = [r for r in matching if r["n_neighbors"] == best_n_neighbors]
                if exact:
                    best = exact[0]
                    return (
                        best["n_neighbors"],
                        best["resolution"],
                        "manual",
                        f"n_neighbors={best_n_neighbors}, resolution={best['resolution']:.1f} "
                        f"(configured) silhouette={best['silhouette_score']:.4f} "
                        f"k={best['n_clusters']}",
                    )
                # n_neighbors not found at this resolution → fall through to
                # pick best silhouette at the given resolution
            best = max(matching, key=lambda r: r["silhouette_score"])
            return (
                best["n_neighbors"],
                best["resolution"],
                "manual",
                f"resolution={best['resolution']:.1f} (configured) "
                f"silhouette={best['silhouette_score']:.4f} k={best['n_clusters']}",
            )
        # best_resolution set but not in grid → fall through to auto
    # Fallback: max silhouette
    best = max(valid, key=lambda r: r["silhouette_score"])
    return (
        best["n_neighbors"],
        best["resolution"],
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
    pts = np.array([(r["n_clusters"], r["silhouette_score"]) for r in valid])

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
            best["n_neighbors"],
            best["resolution"],
            "pareto_elbow",
            f"single_pareto_point silhouette={best['silhouette_score']:.4f} "
            f"k={best['n_clusters']}",
        )

    # -- Normalize --
    eps = 1e-10
    k_norm = (pareto_k - pareto_k.min()) / (pareto_k.max() - pareto_k.min() + eps)
    s_norm = (pareto_s - pareto_s.min()) / (pareto_s.max() - pareto_s.min() + eps)

    # Distance to ideal point (k_norm=0, s_norm=1)
    dist = np.sqrt((1.0 - s_norm) ** 2 + k_norm**2)
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
        best["n_neighbors"],
        best["resolution"],
        "pareto_elbow",
        f"dist_to_ideal={dist[elbow_idx]:.4f} "
        f"silhouette={best['silhouette_score']:.4f} k={best['n_clusters']}"
        f"{delta_note}",
    )


# ═══════════════════════════════════════════════════════════════════════
#  Multi-metric scoring helpers
# ═══════════════════════════════════════════════════════════════════════


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

    # -- Single-result fast path: skip multi-metric scoring --
    if n == 1:
        r = valid[0]
        return (
            r["n_neighbors"],
            r["resolution"],
            "multi_metric",
            f"single_result k={r['n_clusters']} composite={r.get('composite_score', 0):.4f}",
        )

    # ── Gather raw scores ──
    sil_scores = np.array([r["silhouette_score"] for r in valid])
    stab_scores = np.array([r.get("stability_score") or 0.0 for r in valid])

    has_coherence = any(
        "cluster_coherence" in r and r["cluster_coherence"] is not None for r in valid
    )
    coh_scores: np.ndarray = np.zeros(n)
    if has_coherence:
        coh_scores = np.array([r.get("cluster_coherence", 0.0) for r in valid])

    has_splitting_gain = any("splitting_gain" in r for r in valid)
    split_scores: np.ndarray = np.zeros(n)
    if has_splitting_gain:
        split_scores = np.array([r.get("splitting_gain", 0.0) for r in valid])
    has_kb_rate = any("kb_annotatable_rate" in r for r in valid)
    kb_scores: np.ndarray = np.zeros(n)
    if has_kb_rate:
        kb_scores = np.array([r.get("kb_annotatable_rate", 0.0) for r in valid])

    # ── Determine initial weights ──
    if weights is None:
        if has_coherence and has_splitting_gain and has_kb_rate:
            active_weights = dict(DEFAULT_MULTI_METRIC_WEIGHTS)
        elif has_coherence and has_splitting_gain:
            # Degrade to 4-metric (no kb_annotatable_rate)
            active_weights = {
                "silhouette": 0.15,
                "stability": 0.20,
                "cluster_coherence": 0.50,
                "splitting_gain": 0.15,
            }
        elif has_coherence:
            # Degrade to 3-metric (no splitting_gain)
            active_weights = {"silhouette": 0.20, "stability": 0.25, "cluster_coherence": 0.55}
        else:
            active_weights = {"silhouette": 0.5, "stability": 0.5}
    else:
        active_weights = dict(weights)

    # Remove cluster_coherence from weights if entries lack it
    if not has_coherence:
        active_weights.pop("cluster_coherence", None)

    # Remove splitting_gain from weights if entries lack it
    if not has_splitting_gain:
        active_weights.pop("splitting_gain", None)

    # Remove kb_annotatable_rate from weights if entries lack it
    if not has_kb_rate:
        active_weights.pop("kb_annotatable_rate", None)

    # -- Coherence mismatch auto-degrade: if all entries have cluster_coherence < 0.1 --
    if (
        has_coherence
        and any(c is not None for c in coh_scores)
        and float(np.max([c for c in coh_scores if c is not None])) < 0.1
    ):
        logger.warning(
            "Max cluster_coherence=%.4f < 0.1 across all entries — marker_dict may be mismatched. "
            "Degrading to silhouette+stability only.",
            float(np.max([c for c in coh_scores if c is not None])),
        )
        has_coherence = False
        active_weights = {"silhouette": 0.5, "stability": 0.5}

    # ── Low-variance guard (on raw scores) ──
    metrics_raw = {
        "silhouette": sil_scores,
        "stability": stab_scores,
    }
    if has_coherence:
        metrics_raw["cluster_coherence"] = coh_scores
    if has_splitting_gain:
        metrics_raw["splitting_gain"] = split_scores
    if has_kb_rate:
        metrics_raw["kb_annotatable_rate"] = kb_scores

    for metric_name in list(active_weights.keys()):
        scores = metrics_raw.get(metric_name)
        if scores is None:
            continue
        score_range = float(np.max(scores) - np.min(scores))
        if score_range < 0.01:
            logger.warning(
                "%s variance < 0.01 (range=%.4f) — disabling metric",
                metric_name,
                score_range,
            )
            del active_weights[metric_name]

    if not active_weights:
        logger.warning("All metrics dropped — falling back to silhouette only")
        active_weights = {"silhouette": 1.0}

    # ── Renormalise weights to sum=1.0 ──
    total_w = sum(active_weights.values())
    if total_w > 0:
        for k in active_weights:
            active_weights[k] /= total_w

    # ── Build weight string for reason ──
    _weight_short = {
        "silhouette": "sil",
        "stability": "stab",
        "cluster_coherence": "coh",
        "splitting_gain": "split",
        "kb_annotatable_rate": "kb",
    }
    _wp: list[str] = []
    for _mk in [
        "silhouette",
        "stability",
        "cluster_coherence",
        "splitting_gain",
        "kb_annotatable_rate",
    ]:
        if _mk in active_weights:
            _wp.append(f"{_weight_short[_mk]}:{active_weights[_mk]:.2f}")
    weights_str = "weights=" + ",".join(_wp)

    # ── Normalise each metric to [0, 1] ──
    def _normalize(values):
        n = len(values)
        if n <= 1:
            return np.array([0.5])
        ranks = stats.rankdata(values)
        return (ranks - 1) / (n - 1)

    norm = {}
    if "silhouette" in active_weights:
        norm["silhouette"] = _normalize(sil_scores)
    if "stability" in active_weights:
        norm["stability"] = _normalize(stab_scores)
    if "cluster_coherence" in active_weights and has_coherence:
        norm["cluster_coherence"] = _normalize(coh_scores)
    if "splitting_gain" in active_weights and has_splitting_gain:
        norm["splitting_gain"] = _normalize(split_scores)
    if "kb_annotatable_rate" in active_weights and has_kb_rate:
        norm["kb_annotatable_rate"] = _normalize(kb_scores)

    # ── Composite score ──
    composite = np.zeros(n)
    for metric_name, w in active_weights.items():
        composite += w * norm[metric_name]

    # ── 3-tier resolution recommendation ──
    tier_str = ""
    try:
        entries_by_resolution = sorted(
            [(valid[i], composite[i]) for i in range(n)], key=lambda x: x[0]["resolution"]
        )
        comp_max = max(composite)
        stab_max = max(stab_scores) if "stability" in active_weights else None

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
                stab = entry.get("stability_score", 0.0)
                if stab > 0.85 * stab_max:
                    fine_entry = entry
                    break

        coarse_r = f"{coarse_entry['resolution']:.2f}" if coarse_entry else "nan"
        coarse_k = coarse_entry["n_clusters"] if coarse_entry else 0
        fine_r = f"{fine_entry['resolution']:.2f}" if fine_entry else "nan"
        fine_k = fine_entry["n_clusters"] if fine_entry else 0

        tier_str = (
            f"coarse: r={coarse_r} (k={coarse_k})"
            f" / balanced: r={balanced_entry['resolution']:.2f} (k={balanced_entry['n_clusters']})"
            f" / fine: r={fine_r} (k={fine_k})"
        )
        (log or logger).info(
            "[multi_metric 3-tier] %s",
            tier_str,
        )
    except Exception as exc:
        logger.debug("3-tier computation skipped: %s", exc)

    best_idx = int(np.argmax(composite))
    best = valid[best_idx]

    # ── Build reason string ──
    sil_val = best["silhouette_score"]
    stab_val = best.get("stability_score", 0.0)
    coh_val = best.get("cluster_coherence", 0.0) if has_coherence else 0.0
    split_val = best.get("splitting_gain", 0.0) if has_splitting_gain else 0.0
    kb_val = best.get("kb_annotatable_rate", 0.0) if has_kb_rate else 0.0
    k = best["n_clusters"]

    sil_norm_val = norm.get("silhouette", np.zeros(n))[best_idx]
    stab_norm_val = norm.get("stability", np.zeros(n))[best_idx]

    reason = (
        f"composite={composite[best_idx]:.4f} "
        f"sil={sil_val:.4f}(n={sil_norm_val:.2f}) "
        f"stab={stab_val:.3f}(n={stab_norm_val:.2f}) "
        f"coherence={coh_val:.3f} "
        f"split_gain={split_val:.3f} kb_rate={kb_val:.3f} k={k}"
        f" {weights_str}"
        f" {tier_str}"
    ).rstrip()

    return (
        best["n_neighbors"],
        best["resolution"],
        "multi_metric",
        reason,
    )


def select_best_umap_params(
    adata,
    best_n,
    min_dist_grid,
    spread_grid,
    method,
    CFG,  # noqa: N803
    use_rep,
    log,
    device="cpu",
    metric="trustworthiness",
):  # noqa: N803
    """Sweep min_dist × spread on the best (n_neighbors) neighbor graph,
    or use manual fallback.

    The KNN graph is rebuilt once (or reused if already present).

    Parameters
    ----------
    adata : AnnData
        Must already have PCA representation (use_rep).
    best_n : int
        Best n_neighbors value (from select_best_params).
    min_dist_grid : list of float or None
        Values to sweep in auto-sweep mode.
    spread_grid : list of float or None
        Values to sweep in auto-sweep mode.
    method : str or None
        "convex_hull" | None
    CFG : Config
    use_rep : str
        Key in adata.obsm for PCA (e.g. 'X_pca_harmony' or 'X_pca').
    log : logging.Logger
    device : str
        "cpu" or GPU device name.
    metric : str
        "trustworthiness" | "convex_hull" | "fixed"

        Controls how UMAP parameters are selected:
        - "trustworthiness": compute sklearn trustworthiness for each
          (min_dist, spread) combo, pick highest score (default).
        - "convex_hull": pick combo with largest convex hull area (legacy).
        - "fixed": skip sweep entirely, use CFG.clustering.umap_min_dist
          and umap_spread as-is (saves ~3700s on large datasets).

    Returns
    -------
    best_min_dist : float
    best_spread : float
    method_label : str
        Human-readable method name for logging.
    results : list of dict
        Each dict: {min_dist, spread, <metric>_score} plus "coords" when run.
        Empty list if sweep was skipped.
    """
    import scanpy as sc

    use_paga = getattr(CFG.clustering, "umap_paga_init", False)

    # ── Manual mode ──
    if method is None:
        md = getattr(CFG.clustering, "umap_min_dist", 0.3)
        sp = getattr(CFG.clustering, "umap_spread", 1.0)
        log.info("UMAP params (manual): min_dist=%.2f, spread=%.1f", md, sp)
        return md, sp, "manual", []

    # method is not None → auto-sweep mode. The scoring metric is controlled by
    # the `metric` parameter (umap_selection_metric in config).

    # ── Auto-sweep (empty grid → fallback) ──
    if (
        min_dist_grid is None
        or spread_grid is None
        or (len(min_dist_grid) <= 1 and len(spread_grid) <= 1)
    ):
        md = getattr(CFG.clustering, "umap_min_dist", 0.3)
        sp = getattr(CFG.clustering, "umap_spread", 1.0)
        log.info(
            "UMAP params (%s, empty grid → fallback): min_dist=%.2f, spread=%.1f", metric, md, sp
        )
        return md, sp, metric, []

    # ── Reuse neighbor graph if it already matches best_n + use_rep ──
    # grid_search_clustering already built KNN for the selected n_neighbors;
    # rebuilding here is a redundant O(n_obs × log n_obs) cost on 1M-cell datasets.
    _need_rebuild = True
    try:
        _nb_params = adata.uns.get("neighbors", {}).get("params", {}) or {}
        if (
            _nb_params.get("n_neighbors") == best_n
            and _nb_params.get("use_rep") == use_rep
            and "connectivities" in adata.obsp
            and "distances" in adata.obsp
        ):
            log.info(
                "Reusing existing KNN graph (n_neighbors=%d, use_rep=%s) for UMAP sweep",
                best_n,
                use_rep,
            )
            _need_rebuild = False
    except Exception:
        pass

    if _need_rebuild:
        log.info("Building KNN graph (n_neighbors=%d) for UMAP parameter sweep...", best_n)
        try:
            if device != "cpu":
                from core.utils._gpu import gpu_neighbors

                gpu_neighbors(
                    adata,
                    log=log,
                    device=device,
                    n_neighbors=best_n,
                    n_pcs=min(CFG.pca.n_pcs_use, adata.obsm[use_rep].shape[1]),
                    use_rep=use_rep,
                    random_state=CFG.execution.random_seed,
                )
            else:
                sc.pp.neighbors(
                    adata,
                    n_neighbors=best_n,
                    n_pcs=min(CFG.pca.n_pcs_use, adata.obsm[use_rep].shape[1]),
                    use_rep=use_rep,
                    random_state=CFG.execution.random_seed,
                )
        except Exception as e:
            log.error("KNN graph build failed for UMAP sweep: %s", e)
            return (
                getattr(CFG.clustering, "umap_min_dist", 0.3),
                getattr(CFG.clustering, "umap_spread", 1.0),
                metric,
                [],
            )
    # ── Sweep ──
    results = []
    if metric == "trustworthiness":
        best_score = -1.0
    else:
        best_area = -1.0
    best_md = min_dist_grid[0]
    best_sp = spread_grid[0]

    # Cache of the previous UMAP embedding — used to warm-start the next
    # combo instead of paying for spectral initialization every iteration.
    _prev_embedding = None
    for md in min_dist_grid:
        for sp in spread_grid:
            # ── Stage 1: UMAP embedding ──
            try:
                # UMAP init: first combo uses spectral/paga; subsequent combos
                # warm-start from the previous embedding (legal scanpy API,
                # much faster than re-running spectral initialization each time).
                _init = (
                    _prev_embedding
                    if _prev_embedding is not None
                    else ("paga" if use_paga else "spectral")
                )
                if device != "cpu":
                    from core.utils._gpu import gpu_umap

                    gpu_umap(
                        adata,
                        log=log,
                        device=device,
                        min_dist=md,
                        spread=sp,
                        init_pos=_init,
                        random_state=CFG.execution.random_seed,
                    )
                else:
                    sc.tl.umap(
                        adata,
                        min_dist=md,
                        spread=sp,
                        init_pos=_init,
                        random_state=CFG.execution.random_seed,
                    )
                coords = adata.obsm["X_umap"]
                _prev_embedding = np.asarray(coords).copy()
            except Exception as e:
                # UMAP itself broke — no embedding produced, nothing to score.
                log.warning(
                    "  UMAP step failed (min_dist=%.2f, spread=%.1f): %s",
                    md,
                    sp,
                    e,
                )
                entry = {"min_dist": md, "spread": sp}
                if metric == "trustworthiness":
                    entry["trustworthiness"] = None
                else:
                    entry["convex_hull_area"] = None
                results.append(entry)
                continue

            # ── Stage 2/3: scoring (UMAP already succeeded; coords are valid).
            # Separate try so a scoring-side failure (bad import, ConvexHull
            # degeneracy, …) is not misattributed to the UMAP step above.
            try:
                if metric == "trustworthiness":
                    # sklearn ≥1.6 renamed `t_sne.py` → `_t_sne.py` (private);
                    # `trustworthiness` is now exported directly from sklearn.manifold.
                    from sklearn.manifold import trustworthiness

                    score = float(
                        trustworthiness(
                            adata.obsm[use_rep],
                            coords,
                            n_neighbors=15,
                        )
                    )
                    results.append(
                        {
                            "min_dist": md,
                            "spread": sp,
                            "trustworthiness": score,
                            "coords": _prev_embedding,
                        }
                    )
                    log.info("  min_dist=%.2f, spread=%.1f → trustworthiness=%.4f", md, sp, score)
                    if score > best_score:
                        best_score = score
                        best_md = md
                        best_sp = sp
                else:
                    hull = ConvexHull(coords)
                    area = float(hull.volume)  # 2D → area
                    results.append(
                        {
                            "min_dist": md,
                            "spread": sp,
                            "convex_hull_area": area,
                            # Save coords so the comparison figure can reuse them
                            # without re-running UMAP (huge win on 1M-cell datasets:
                            # ~3 × 47min cold spectral = ~2.4h saved).
                            "coords": _prev_embedding,
                        }
                    )
                    log.info("  min_dist=%.2f, spread=%.1f → convex_hull_area=%.2f", md, sp, area)
                    if area > best_area:
                        best_area = area
                        best_md = md
                        best_sp = sp
            except Exception as e:
                # UMAP ran fine; only the metric computation broke. Report it as
                # a scoring failure so the log doesn't blame the embedding.
                label = (
                    "Trustworthiness scoring"
                    if metric == "trustworthiness"
                    else "ConvexHull scoring"
                )
                log.warning(
                    "  %s failed (min_dist=%.2f, spread=%.1f): %s",
                    label,
                    md,
                    sp,
                    e,
                )
                entry = {"min_dist": md, "spread": sp}
                if metric == "trustworthiness":
                    entry["trustworthiness"] = None
                else:
                    entry["convex_hull_area"] = None
                results.append(entry)

    if metric == "trustworthiness":
        log.info(
            "Best UMAP params (trustworthiness): min_dist=%.2f, spread=%.1f (score=%.4f)",
            best_md,
            best_sp,
            best_score,
        )
    else:
        log.info(
            "Best UMAP params (convex_hull): min_dist=%.2f, spread=%.1f (area=%.2f)",
            best_md,
            best_sp,
            best_area,
        )
    return best_md, best_sp, metric, results
