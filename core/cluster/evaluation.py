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

import logging
from typing import Literal

import numpy as np
from scipy import stats
from scipy.spatial import ConvexHull

logger = logging.getLogger(__name__)

DEFAULT_MULTI_METRIC_WEIGHTS = {
    "silhouette": 0.15,
    "stability": 0.20,
    "cluster_coherence": 0.45,
    "splitting_gain": 0.15,
    "kb_annotatable_rate": 0.05,
}


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


def _compute_stability(
    adata,
    resolution,
    leiden_flavor: Literal["leidenalg", "igraph"] = "igraph",
    n_seeds=5,
    base_seed=42,
    device: str = "cpu",
    cfg=None,
    n_iterations=-1,
):
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
    # Resolve n_iterations from config when available
    if cfg is not None:
        n_iterations = getattr(cfg.clustering, "stability_leiden_n_iterations", n_iterations)

    if n_seeds <= 1:
        return float("nan")

    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score

    seeds = range(base_seed, base_seed + n_seeds)
    temp_keys = [f"_temp_stab_{i}" for i in range(n_seeds)]

    label_sets, consecutive_aris = [], []
    for i, seed in enumerate(seeds):
        key = temp_keys[i]
        try:
            if device != "cpu":
                # Lazy import to avoid circular dependency at module load.
                from core.utils._gpu import gpu_leiden

                gpu_leiden(
                    adata,
                    resolution=resolution,
                    key_added=key,
                    random_state=seed,
                    # cuGraph doesn't accept flavor/directed/n_iterations;
                    # gpu_leiden strips them on the GPU path automatically.
                    flavor=leiden_flavor,
                    n_iterations=n_iterations,
                    directed=False,
                )
            else:
                sc.tl.leiden(
                    adata,
                    resolution=resolution,
                    key_added=key,
                    random_state=seed,
                    flavor=leiden_flavor,
                    n_iterations=n_iterations,
                    directed=False,
                )
            label_sets.append(adata.obs[key].values)
        except Exception as e:
            logger.warning(f"Stability seed {seed} failed: {e}")
        finally:
            if key in adata.obs.columns:
                del adata.obs[key]
        # Early termination: check if last 2 consecutive ARIs are 1.0
        if len(label_sets) >= 2:
            ari = adjusted_rand_score(label_sets[-2], label_sets[-1])
            consecutive_aris.append(ari)
            if len(consecutive_aris) >= 2 and all(a == 1.0 for a in consecutive_aris[-2:]):
                logger.info(
                    f"Early termination at seed {seed + 1}/{n_seeds}: "
                    "last 2 consecutive ARIs are 1.0"
                )
                break

    # Final cleanup: remove any remaining temp columns
    for col in list(adata.obs.columns):
        if col.startswith("_temp_stab_"):
            del adata.obs[col]

    n_runs = len(label_sets)
    if n_runs <= 1:
        return float("nan")

    # Compute pairwise ARI
    aris = []
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            ari = adjusted_rand_score(label_sets[i], label_sets[j])
            aris.append(ari)

    return float(np.mean(aris))


def _compute_cluster_coherence(
    adata,
    cluster_key,
    per_cell_scores,
    dominance_threshold=2.5,
    min_expression=0.05,
    min_cluster_size=10,
    log=None,
):
    """Per-cluster marker coherence metric that peaks at intermediate resolution.

    For each cluster label in adata.obs[cluster_key]:
      - Compute mean of per_cell_scores per cell type
      - Find best_score (max) and second_best_score
      - A cluster is "coherent" if:
        (a) best_score > min_expression (not noise)
        (b) best_score / max(second_best, 1e-10) > dominance_threshold
      - This ensures the cluster has ONE clearly dominant cell type

    coherence = n_coherent / n_valid_clusters
    where n_valid_clusters excludes clusters smaller than min_cluster_size.

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
    min_expression : float or str
        Minimum mean expression for a cluster to be considered (not noise).
        If a string starting with 'p' (e.g. 'p25'), interpreted as a percentile
        of per-cell scores for the best-matching cell type.
    min_cluster_size : int
        Minimum number of cells a cluster must have to be evaluated.
        Clusters below this threshold are skipped.
    log : Logger or None

    Returns
    -------
    float
        Fraction of clusters that are coherent (0.0-1.0).
        NaN if *per_cell_scores* is empty or has no valid entries.
    """
    if not per_cell_scores:
        return float("nan")

    # Check for any valid (non-None) entries
    has_valid = any(v is not None for v in per_cell_scores.values())
    if not has_valid:
        return float("nan")

    cluster_labels = adata.obs[cluster_key].values
    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)

    if n_clusters <= 1:
        return 1.0

    cell_types = list(per_cell_scores.keys())
    n_cells = len(cluster_labels)
    n_coherent = 0
    n_valid_clusters = 0
    _log = log or logger

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        cluster_size = int(np.sum(mask))

        # Skip clusters smaller than min_cluster_size
        if cluster_size < min_cluster_size:
            _log.debug(
                "Skipping cluster %s in coherence: size %d < min_cluster_size %d",
                cluster,
                cluster_size,
                min_cluster_size,
            )
            continue

        n_valid_clusters += 1
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
        best_ct = mean_scores[0][0]

        # Resolve threshold (float or percentile string)
        threshold = min_expression
        if isinstance(min_expression, str) and min_expression.startswith("p"):
            pct = int(min_expression[1:])
            if best_ct in per_cell_scores and per_cell_scores[best_ct] is not None:
                threshold = float(np.percentile(per_cell_scores[best_ct], pct))

        if top_score <= threshold:
            continue

        if len(mean_scores) < 2:
            return float("nan")
        else:
            second_score = mean_scores[1][1]
            if top_score / max(second_score, 1e-10) > dominance_threshold:
                n_coherent += 1

    if n_valid_clusters == 0:
        return float("nan")

    return n_coherent / n_valid_clusters


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
    sorted_entries = sorted(valid_by_resolution, key=lambda e: e.get("resolution", 0.0))

    gains = {}
    for i in range(1, len(sorted_entries)):
        prev = sorted_entries[i - 1]
        curr = sorted_entries[i]
        r_prev = prev["resolution"]
        r_curr = curr["resolution"]
        k_prev = prev["n_clusters"]
        k_curr = curr["n_clusters"]

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
    sil_scores = np.array([r["silhouette_score"] for r in valid])
    stab_scores = np.array([r.get("stability_score", 0.0) for r in valid])

    has_coherence = any("cluster_coherence" in r for r in valid)
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
    if has_coherence and float(np.max(coh_scores)) < 0.1:
        logger.warning(
            "Max cluster_coherence=%.4f < 0.1 across all entries — marker_dict may be mismatched. "
            "Degrading to silhouette+stability only.",
            float(np.max(coh_scores)),
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


def _detect_granularity(
    results_summary: list[dict], cv_threshold: float = 0.05, min_clusters: int = 10
) -> str:
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
        if "n_neighbors" not in entry:
            continue
        nn = entry["n_neighbors"]
        groups.setdefault(nn, []).append(entry)

    if not groups:
        return "tissue"

    # Pick median-size group (by number of entries)
    group_sizes = sorted(groups.items(), key=lambda kv: len(kv[1]))
    median_idx = len(group_sizes) // 2
    _median_nn, median_group = group_sizes[median_idx]

    # 2. Collect silhouette scores and n_clusters from the median group
    #    Sort by resolution for deterministic ordering
    median_group_sorted = sorted(median_group, key=lambda e: e.get("resolution", 0.0))

    silhouette_values = []
    n_clusters_values = []
    for entry in median_group_sorted:
        if "silhouette_score" not in entry or entry["silhouette_score"] is None:
            continue
        if "n_clusters" not in entry or entry["n_clusters"] is None:
            continue
        silhouette_values.append(entry["silhouette_score"])
        n_clusters_values.append(entry["n_clusters"])

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


def _compute_pairwise_de_markers(
    adata,
    cluster_key,
    n_genes=50,
    padj_threshold=0.05,
    lfc_threshold=1.0,
    max_clusters_for_pairwise=30,
    log=None,
):
    """Compute per-cluster pairwise DE marker counts (Shekhar 2016 semantics).

    For each cluster C_i, compares against every other cluster C_j individually.
    A gene is a pairwise marker for C_i if it is upregulated (padj < threshold
    AND log2FC > threshold) in ALL pairwise comparisons C_i vs C_j.

    Falls back to one-vs-rest if n_clusters > max_clusters_for_pairwise
    (unless max_clusters_for_pairwise=0, meaning no cap).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts in .raw.
    cluster_key : str
        Key in adata.obs storing cluster assignments.
    n_genes : int
        Number of top DE genes to request from rank_genes_groups.
    padj_threshold : float
        Adjusted p-value threshold for significance.
    lfc_threshold : float
        Log2 fold-change threshold for upregulation.
    max_clusters_for_pairwise : int
        Maximum number of clusters for pairwise mode. 0 = no cap.
    log : logging.Logger or None
        Logger for warnings.

    Returns
    -------
    dict[str, int]
        Cluster name → pairwise marker count.
    """
    import scanpy as sc

    unique_vals = adata.obs[cluster_key].unique()
    # Handle categorical with possibly None/unused categories
    clusters = sorted([str(v) for v in unique_vals if v is not None])
    n_clusters = len(clusters)

    # Fallback: too many clusters
    if max_clusters_for_pairwise > 0 and n_clusters > max_clusters_for_pairwise:
        if log:
            log.warning(
                "n_clusters=%d exceeds de_pairwise_max_clusters=%d — "
                "falling back to one-vs-rest DE",
                n_clusters,
                max_clusters_for_pairwise,
            )
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            n_genes=n_genes,
            use_raw=True,
            key_added="_pairwise_fallback",
        )
        result = adata.uns["_pairwise_fallback"]
        marker_counts = {}
        for cluster_name in clusters:
            pvals = result["pvals_adj"][cluster_name]
            lfcs = result["logfoldchanges"][cluster_name]
            is_sig = (pvals < padj_threshold) & (lfcs > lfc_threshold)
            marker_counts[cluster_name] = int(is_sig.sum())
        del adata.uns["_pairwise_fallback"]
        return marker_counts

    # Pairwise mode
    marker_counts = {}
    for ci in clusters:
        marker_set = None
        for cj in clusters:
            if cj == ci:
                continue
            sc.tl.rank_genes_groups(
                adata,
                groupby=cluster_key,
                groups=[ci],
                reference=cj,
                method="wilcoxon",
                n_genes=n_genes,
                use_raw=True,
                key_added="_pairwise_tmp",
            )
            tmp = adata.uns["_pairwise_tmp"]
            pvals = tmp["pvals_adj"][ci]
            lfcs = tmp["logfoldchanges"][ci]
            is_sig = (pvals < padj_threshold) & (lfcs > lfc_threshold)
            # Use gene indices as proxy for identity (names are mock-safe)
            sig_mask = np.asarray(is_sig, dtype=bool)
            sig_indices = set(np.where(sig_mask)[0])
            marker_set = sig_indices if marker_set is None else (marker_set & sig_indices)
            del adata.uns["_pairwise_tmp"]

        marker_counts[ci] = len(marker_set) if marker_set is not None else 0

    return marker_counts


def _select_de_gated(valid, adata, de_gate_threshold=25, pairwise_max_clusters=30):
    """Select best resolution using DE-gated criterion for subtype-level data.

    For subtype-level data where silhouette scores are flat across resolutions,
    this method selects the highest resolution that maintains a minimum number of
    differentially expressed genes for every cluster.

    Uses pairwise DE comparisons (Shekhar 2016 semantics): each cluster is
    compared against every other cluster individually. A gene is a pairwise
    marker for cluster C_i if it is upregulated in ALL pairwise comparisons
    C_i vs C_j. Falls back to one-vs-rest when n_clusters > pairwise_max_clusters
    (0 = always pairwise).

    Parameters
    ----------
    valid : list[dict]
        Grid-search entries with keys: n_clusters, resolution, cluster_key.
        Must be unique by resolution.
    adata : AnnData
        Annotated data matrix.
    de_gate_threshold : int
        Minimum number of pairwise DE genes required for every cluster.
        Default 25.
    pairwise_max_clusters : int
        Maximum number of clusters for pairwise comparisons. Above this
        number, falls back to one-vs-rest DE to avoid O(n²) cost.
        Default 30.

    Returns
    -------
    tuple[int, float, str, str]
        (n_clusters, resolution, cluster_key, reason_str)
    """

    # Edge case: single entry -> return as-is
    if len(valid) <= 1:
        entry = valid[0]
        return (
            entry["n_clusters"],
            entry["resolution"],
            entry["cluster_key"],
            "de_gated(single_entry, min_de=N/A)",
        )

    # Sort by resolution ascending for deterministic iteration
    sorted_entries = sorted(valid, key=lambda e: e.get("resolution", 0.0))

    # Collect (entry, min_de) for each resolution
    candidates = []  # list of (entry, min_de)

    for entry in sorted_entries:
        cluster_key = entry["cluster_key"]

        try:
            marker_counts = _compute_pairwise_de_markers(
                adata,
                cluster_key,
                max_clusters_for_pairwise=pairwise_max_clusters,
                log=logger,
            )
            min_de = min(marker_counts.values()) if marker_counts else 0
        except Exception as e:
            logger.warning(
                "pairwise DE failed for cluster_key=%s (resolution=%.2f): %s",
                cluster_key,
                entry["resolution"],
                e,
            )
            continue

        candidates.append((entry, min_de))

    # Clean up temp rank_genes_groups keys (belt-and-suspenders)
    for k in list(adata.uns):
        if k.startswith("_de_gated_") or k.startswith("_pairwise_"):
            del adata.uns[k]

    if not candidates:
        # All entries failed -> fallback to first entry
        logger.warning(
            "DE-gated selection: all pairwise DE calls failed -- fallback to first entry"
        )
        entry = sorted_entries[0]
        return (
            entry["n_clusters"],
            entry["resolution"],
            entry["cluster_key"],
            "de_gated(all_failed, min_de=N/A)",
        )

    # Select entry with highest n_clusters where min_de >= threshold
    candidates_by_n = sorted(
        candidates,
        key=lambda c: c[0]["n_clusters"],
        reverse=True,
    )
    best_entry = None
    best_min_de = -1
    for entry, min_de in candidates_by_n:
        if min_de >= de_gate_threshold:
            best_entry = entry
            best_min_de = min_de
            break

    # If no entry meets threshold, fallback to entry with highest DE count
    if best_entry is None:
        best_entry, best_min_de = max(candidates, key=lambda c: c[1])
        logger.info(
            "DE-gated selection: best_resolution=%.2f, min_de=%d "
            "(fallback: no entry met threshold=%d)",
            best_entry["resolution"],
            best_min_de,
            de_gate_threshold,
        )
    else:
        logger.info(
            "DE-gated selection: best_resolution=%.2f, min_de=%d",
            best_entry["resolution"],
            best_min_de,
        )

    return (
        best_entry["n_clusters"],
        best_entry["resolution"],
        best_entry["cluster_key"],
        f"de_gated(min_de={best_min_de})",
    )


def enrich_grid_results(
    adata,
    results_summary,
    cfg,
    compute_per_cell_scores_fn=None,
    log=None,
    use_rep=None,
):
    """Enrich *results_summary* entries with multi-metric scores.

    For each (n_neighbors, resolution) entry, computes:
      - stability_score  (via cross-seed ARI)
      - cluster_coherence (via marker gene per-cell scores)
      - kb_annotatable_rate (via tissue knowledge base, if *cfg.tissue_kb* is set)
      - splitting_gain (per n_neighbors group)

    The list is modified **in place** and returned. The caller is
    responsible for persistence (e.g. via ``param_grid_summary.csv``).

    Parameters
    ----------
    adata : AnnData
    results_summary : list[dict]
        Each dict has keys: *n_neighbors*, *resolution*, *cluster_key*.
        Enriched with *stability_score*, *cluster_coherence*, etc. IN PLACE.
    cfg : Config
    compute_per_cell_scores_fn : callable or None
        ``(adata, cfg) -> dict[str, np.ndarray]``
        Pre-computed per-cell scores for each cell type.
        When ``None`` and ``cfg.marker.marker_dict`` has entries, the
        default marker-based logic (``sc.tl.score_genes``) is used.
    log : Logger or None
    use_rep : str or None
        Key in ``adata.obsm`` for the embedding (e.g. ``'X_pca'``).
        Required when the KNN graph needs rebuilding.

    Returns
    -------
    list[dict]
        The same *results_summary* (modified in place).
    """
    import logging as _logging

    import numpy as np
    import scanpy as sc

    if log is None:
        log = _logging.getLogger(__name__)

    marker_dict = getattr(cfg.marker, "marker_dict", None) or {}
    has_markers = bool(marker_dict)
    n_stab_seeds = getattr(cfg.clustering, "stability_n_seeds", 12)
    dominance_threshold = getattr(cfg.clustering, "multi_metric_coverage_ratio_threshold", 2.5)
    leiden_flavor = getattr(cfg.clustering, "leiden_flavor", "igraph")
    device = getattr(cfg.execution, "device", "cpu")

    # Group results by n_neighbors
    by_n = {}
    for r in results_summary:
        n_ = r.get("n_neighbors")
        by_n.setdefault(n_, []).append(r)

    for n_val, group in by_n.items():
        # --- Check if KNN rebuild is needed (optimisation) ---
        _need_rebuild_knn = True
        try:
            _nb_params = adata.uns.get("neighbors", {}).get("params", {}) or {}
            if (
                _nb_params.get("n_neighbors") == n_val
                and use_rep is not None
                and _nb_params.get("use_rep") == use_rep
                and "connectivities" in adata.obsp
            ):
                _need_rebuild_knn = False
        except Exception:
            pass

        if _need_rebuild_knn:
            if use_rep is None:
                log.warning(
                    "KNN rebuild needed but use_rep is None -- skipping n_neighbors=%d",
                    n_val,
                )
                continue
            try:
                sc.pp.neighbors(
                    adata,
                    n_neighbors=n_val,
                    n_pcs=cfg.pca.n_pcs_use,
                    use_rep=use_rep,
                    random_state=cfg.execution.random_seed,
                )
            except Exception as e:
                log.warning(
                    "KNN rebuild failed for n_neighbors=%d: %s -- skipping group",
                    n_val,
                    e,
                )
                continue

        # --- Pre-compute per-cell scores ---
        per_cell_scores = {}
        if compute_per_cell_scores_fn is not None:
            try:
                scores = compute_per_cell_scores_fn(adata, cfg)
                if scores is not None:
                    per_cell_scores = scores
            except Exception as e:
                log.warning(
                    "per_cell_scores computation failed: %s -- falling back to no markers",
                    e,
                )
        elif has_markers and adata.raw is not None:
            # Default marker_dict-based scoring (RNA and Spatial common case)
            from anndata import utils as anndata_utils

            adata.raw._var.index = anndata_utils.make_index_unique(adata.raw._var.index, join="-")
            try:
                for ct, genes in marker_dict.items():
                    valid_genes = [g for g in genes if g in adata.raw.var_names]
                    if valid_genes:
                        sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f"_score_{ct}")
                        per_cell_scores[ct] = adata.obs[f"_score_{ct}"].values.copy()
                # Clean up temporary score columns
                for col in list(adata.obs.columns):
                    if col.startswith("_score_") and col in adata.obs.columns:
                        adata.obs.drop(columns=[col], inplace=True)
            except Exception as e:
                log.warning(
                    "Marker score pre-computation failed: %s -- falling back to no markers",
                    e,
                )
                per_cell_scores = {}
        elif has_markers and adata.raw is None:
            log.warning(
                "adata.raw is None -- cannot compute marker coverage. Degrading to silhouette+stability only."
            )
            has_markers = False

        # --- Compute stability + marker coverage for each entry ---
        for entry in group:
            try:
                resolution = entry["resolution"]
                ck = entry["cluster_key"]

                entry["stability_score"] = _compute_stability(
                    adata,
                    resolution=resolution,
                    leiden_flavor=leiden_flavor,
                    n_seeds=n_stab_seeds,
                    device=device,
                    cfg=cfg,
                )

                if per_cell_scores:
                    entry["cluster_coherence"] = _compute_cluster_coherence(
                        adata,
                        ck,
                        per_cell_scores,
                        dominance_threshold=dominance_threshold,
                    )

                # --- KB annotatable rate ---
                if getattr(cfg, "tissue_kb", None) and per_cell_scores:
                    labels = adata.obs[ck].values
                    unique_clusters = np.unique(labels)
                    n_total = len(unique_clusters)
                    n_annotatable = 0
                    for cl in unique_clusters:
                        mask_ = labels == cl
                        best_score = 0.0
                        for ct_ in per_cell_scores:
                            scores = per_cell_scores[ct_]
                            if scores is not None and len(scores) == len(labels):
                                mean_val = float(np.mean(scores[mask_]))
                                if mean_val > best_score:
                                    best_score = mean_val
                        if best_score > 0.5:
                            n_annotatable += 1
                    rate = n_annotatable / n_total if n_total > 0 else 0.0
                    entry["kb_annotatable_rate"] = rate
                    log.info("KB annotatable rate: %.3f", rate)

            except Exception as e:
                log.warning(
                    "Enrichment failed for n_neighbors=%d, resolution=%.1f: %s",
                    entry.get("n_neighbors"),
                    entry.get("resolution"),
                    e,
                )
                entry["stability_score"] = None
                entry["cluster_coherence"] = None
                entry["kb_annotatable_rate"] = None

        # --- Compute splitting_gain for this n_neighbors group ---
        if len(group) >= 2:
            group_sorted = sorted(group, key=lambda e: e.get("resolution", 0.0))
            gains = _compute_splitting_gain(group_sorted)
            for entry in group:
                entry["splitting_gain"] = gains.get(entry["resolution"], 0.0)

    return results_summary


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
                    n_pcs=CFG.pca.n_pcs_use,
                    use_rep=use_rep,
                    random_state=CFG.execution.random_seed,
                )
            else:
                sc.pp.neighbors(
                    adata,
                    n_neighbors=best_n,
                    n_pcs=CFG.pca.n_pcs_use,
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
