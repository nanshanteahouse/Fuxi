#!/usr/bin/env python3
"""Funnel grid search — progressive subsample → refine for large datasets.

Strategy
--------
For datasets too large to grid-search over exhaustively (>>100k cells):

  1. **Stratified subsample** via KMeans++ per batch, preserving rare clusters.
  2. Run the **full grid search** on the subsample (via modality-agnostic
     ``full_grid_fn`` callable).
  3. **Rank candidates** by multi-metric composite score.
  4. Take **top‑K** candidates and **re‑validate** on the full dataset.
  5. Return the **best entry** from full re‑validation.

The ``full_grid_fn`` callable decouples funnel logic from the modality‑specific
step script.  The caller wraps ``grid_search_clustering`` + ``enrich_grid_results``
into a single ``(adata, cfg) → list[dict]`` function.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any

import numpy as np
import scanpy as sc
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from core.cluster.evaluation import (
    DEFAULT_MULTI_METRIC_WEIGHTS,
    enrich_grid_results,
)

logger = logging.getLogger(__name__)

# Priority order for auto-detecting the embedding slot
_FALLBACK_REP_ORDER = ("X_integrated", "X_pca", "X_scvi")

# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------


def _resolve_use_rep(adata: Any) -> str:
    """Return the first available embedding key from ``_FALLBACK_REP_ORDER``."""
    for key in _FALLBACK_REP_ORDER:
        if key in adata.obsm:
            return key
    raise KeyError(
        f"No suitable embedding in adata.obsm (tried {_FALLBACK_REP_ORDER}). "
        f"Available: {list(adata.obsm.keys())}"
    )


def _resolve_batch_key(adata: Any) -> str | None:
    """Return the canonical batch / sample column or ``None``."""
    for key in ("batch", "sample"):
        if key in adata.obs:
            return key
    return None


def _resolve_cell_type_key(adata: Any) -> str | None:
    """Return the canonical cell-type column or ``None``."""
    for key in ("true_cell_type", "cell_type", "annotation"):
        if key in adata.obs:
            return key
    return None


def _compute_composite_scores(
    results: list[dict[str, Any]],
    cfg: Any,
    log: logging.Logger | None = None,
) -> None:
    """Add ``composite_score`` to each entry via rank-based multi-metric scoring.

    Mirrors the logic in ``evaluation._select_multi_metric`` so the same weights
    and degradation rules apply.  Operates **in place** on *results*.
    """
    _log = log or logger
    if not results:
        return

    # --- Weights ---
    weights = dict(
        getattr(cfg.clustering, "multi_metric_weights", None) or DEFAULT_MULTI_METRIC_WEIGHTS
    )

    # Detect available metrics
    has_coherence = any("cluster_coherence" in r for r in results)
    has_split = any("splitting_gain" in r for r in results)
    has_kb = any("kb_annotatable_rate" in r for r in results)

    if not has_coherence:
        weights.pop("cluster_coherence", None)
    if not has_split:
        weights.pop("splitting_gain", None)
    if not has_kb:
        weights.pop("kb_annotatable_rate", None)

    # Fallback when everything degraded
    if not weights:
        weights = {"silhouette": 1.0}

    # Renormalise to sum = 1
    total = sum(weights.values())
    if total > 0:
        for k in weights:
            weights[k] /= total

    n = len(results)

    def _normalize(arr: np.ndarray) -> np.ndarray:
        """Rank‑based normalisation to [0, 1] (same as evaluation.py)."""
        if n <= 1:
            return np.array([0.5])
        # NaN guard
        arr_clean = np.where(np.isfinite(arr), arr, 0.0)
        # Low-variance guard
        rng = float(np.max(arr_clean) - np.min(arr_clean))
        if rng < 1e-8:
            return np.full(n, 0.5)
        ranks = stats.rankdata(arr_clean, method="average")
        return (ranks - 1) / (n - 1) if n > 1 else np.array([0.5])

    composite = np.zeros(n)
    for metric, w in weights.items():
        raw = np.array([r.get(metric, 0.0) or 0.0 for r in results], dtype=float)
        composite += w * _normalize(raw)

    for i, r in enumerate(results):
        r["composite_score"] = float(composite[i])


# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------


def subsample_stratified(
    adata: Any,
    size: int = 50000,
    use_rep: str | None = None,
    log: logging.Logger | None = None,
) -> tuple[Any, np.ndarray]:
    """Stratified subsample preserving rare clusters via KMeans++ per batch.

    Parameters
    ----------
    adata : AnnData
    size : int
        Target subsample size.  When ``adata.n_obs <= size`` the full data
        is returned unchanged.
    use_rep : str or None
        Key in ``adata.obsm`` for the embedding (PCA / integrated).
        Auto‑detected when ``None``.
    log : Logger or None

    Returns
    -------
    sub_adata : AnnData
        Subsetted **copy** of ``adata``.
    subsample_idx : np.ndarray
        Integer indices of selected cells in the original ``adata``.
    """
    _log = log or logger

    if use_rep is None:
        use_rep = _resolve_use_rep(adata)

    if use_rep not in adata.obsm:
        raise KeyError(f"Embedding {use_rep!r} not found in adata.obsm")

    n_obs = adata.n_obs
    if n_obs <= size:
        _log.info(
            "subsample_stratified: n_obs (%d) <= size (%d) — returning full data",
            n_obs,
            size,
        )
        return adata, np.arange(n_obs)

    # ---- 1. Batch partition ----
    batch_key = _resolve_batch_key(adata)
    if batch_key is None:
        # Treat everything as one batch
        batch_values: list[Any] = ["_all"]
        batch_counts: dict[Any, int] = {"_all": n_obs}
        batch_indices: dict[Any, np.ndarray] = {"_all": np.arange(n_obs)}
    else:
        cats = adata.obs[batch_key].values
        batch_values = list(np.unique(cats))
        batch_counts = {}
        batch_indices = {}
        for bv in batch_values:
            mask = cats == bv
            batch_indices[bv] = np.where(mask)[0]
            batch_counts[bv] = int(mask.sum())

    # ---- 2. KMeans++ per batch (proportional allocation) ----
    all_selected: set[int] = set()
    forced: set[int] = set()

    for bv in batch_values:
        idx = batch_indices[bv]
        n_batch = batch_counts[bv]
        proportion = n_batch / n_obs
        n_samples = max(1, round(proportion * size))

        if n_batch <= n_samples:
            # Take every cell from tiny batches
            all_selected.update(idx.tolist())
            continue

        emb = adata.obsm[use_rep][idx]

        kmeans = KMeans(
            n_clusters=n_samples,
            init="k-means++",
            random_state=42,
            n_init=3,
        )
        kmeans.fit(emb)
        centroids = kmeans.cluster_centers_

        # Nearest cell to each centroid → farthest‑point coverage
        batch_selected: list[int] = []
        for centroid in centroids:
            dists = np.linalg.norm(emb - centroid, axis=1)
            nearest = int(np.argmin(dists))
            batch_selected.append(int(idx[nearest]))
        all_selected.update(batch_selected)

    # ---- 3. Rare‑cluster preservation ----
    ct_key = _resolve_cell_type_key(adata)
    if ct_key is not None:
        ct_values = adata.obs[ct_key].values
        unique_cts = np.unique(ct_values)
        for ct in unique_cts:
            ct_mask = ct_values == ct
            ct_indices = np.where(ct_mask)[0]
            ct_count = len(ct_indices)
            ct_sampled = sum(1 for i in ct_indices if i in all_selected)

            if ct_count >= 10 and ct_sampled < 5:
                missing = 5 - ct_sampled
                available = [int(i) for i in ct_indices if i not in all_selected]
                if available:
                    rng = np.random.RandomState(42)
                    to_add = rng.choice(
                        np.array(available), min(missing, len(available)), replace=False
                    ).tolist()
                    all_selected.update(to_add)
                    forced.update(to_add)

    subsample_idx = np.array(sorted(all_selected), dtype=int)
    _log.info(
        "subsample_stratified: %d → %d cells (forced %d rare-cluster cells)",
        n_obs,
        len(subsample_idx),
        len(forced),
    )

    sub_adata = adata[subsample_idx].copy()
    return sub_adata, subsample_idx


def _validate_on_full(
    adata: Any,
    top_k_entries: list[dict[str, Any]],
    cfg: Any,
    use_rep: str | None = None,
    log: logging.Logger | None = None,
) -> list[dict[str, Any]]:
    """Re‑validate top‑K candidates on the full dataset.

    For each unique ``n_neighbors`` value, the KNN graph is built once and
    shared across all resolutions in that group.  Enrichment (stability,
    cluster coherence, splitting gain, KB rate) is delegated to
    :func:`enrich_grid_results`.
    """
    _log = log or logger
    if use_rep is None:
        use_rep = _resolve_use_rep(adata)

    # Group by n_neighbors so KNN is built once per group
    by_n: dict[int, list[dict[str, Any]]] = {}
    for entry in top_k_entries:
        n_val = entry["n_neighbors"]
        by_n.setdefault(n_val, []).append(entry)

    full_results: list[dict[str, Any]] = []

    for n_val, group in sorted(by_n.items()):
        _log.info(
            "Full validation: building KNN (n_neighbors=%d) on %d cells …",
            n_val,
            adata.n_obs,
        )
        sc.pp.neighbors(
            adata,
            n_neighbors=n_val,
            n_pcs=min(cfg.pca.n_pcs_use, adata.obsm[use_rep].shape[1]),
            use_rep=use_rep,
            random_state=cfg.execution.random_seed,
        )

        for entry in group:
            resolution = entry["resolution"]
            cluster_key = f"funnel_{n_val}_{resolution}"

            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=cluster_key,
                random_state=cfg.execution.random_seed,
                flavor=getattr(cfg.clustering, "leiden_flavor", "igraph"),
                directed=False,
            )

            labels = adata.obs[cluster_key].values
            n_clusters = int(adata.obs[cluster_key].nunique())

            # Silhouette score (sampled on large data to bound cost)
            n_pcs_use = min(cfg.pca.n_pcs_use, adata.obsm[use_rep].shape[1])
            if adata.n_obs > 10000:
                rng = np.random.RandomState(cfg.execution.random_seed)
                idx = rng.choice(adata.n_obs, 10000, replace=False)
                sil = float(silhouette_score(adata.obsm[use_rep][idx, :n_pcs_use], labels[idx]))
            else:
                sil = float(silhouette_score(adata.obsm[use_rep][:, :n_pcs_use], labels))

            full_results.append(
                {
                    "n_neighbors": n_val,
                    "resolution": resolution,
                    "cluster_key": cluster_key,
                    "n_clusters": n_clusters,
                    "silhouette_score": sil,
                }
            )

    # Enrich with stability, coherence, splitting_gain, kb_annotatable_rate
    enrich_grid_results(adata, full_results, cfg, log=_log, use_rep=use_rep)

    # Compute composite scores from enriched metrics
    _compute_composite_scores(full_results, cfg, log=_log)

    # Compute UMAP on full data (was missing — BugFix)
    # KNN is already built above with the per-group n_neighbors
    _log.info("Full validation: computing UMAP on %d cells …", adata.n_obs)
    try:
        sc.tl.umap(
            adata,
            min_dist=getattr(cfg.clustering, "umap_min_dist", 0.3),
            spread=getattr(cfg.clustering, "umap_spread", 1.0),
            random_state=cfg.execution.random_seed,
        )
    except Exception as e:
        _log.warning("UMAP on full data failed (non-fatal): %s", e)

    return full_results


def run_funnel_grid_search(
    adata: Any,
    cfg: Any,
    full_grid_fn: Callable[[Any, Any], list[dict[str, Any]]],
    log: logging.Logger | None = None,
) -> dict[str, Any]:
    """Progressive funnel grid search — subsample → rank → re‑validate.

    Parameters
    ----------
    adata : AnnData
    cfg : Config
        Must expose ``cfg.clustering`` (ClusteringSettings).
        Optional funnel fields (accessed via ``getattr`` with defaults):

        - ``funnel_subsample_size`` (default 50000)
        - ``funnel_top_k`` (default 3)

    full_grid_fn : Callable[[adata, cfg], list[dict]]
        Modality‑agnostic wrapper that runs the **full** grid search **plus**
        enrichment on the given (possibly subsampled) ``adata`` and returns an
        enriched ``results_summary`` — a list of dicts with keys including
        ``n_neighbors``, ``resolution``, ``cluster_key``, ``n_clusters``,
        ``silhouette_score``, ``stability_score``, ``cluster_coherence``,
        ``splitting_gain``, and ``kb_annotatable_rate``.

        The caller is expected to compose::

            def _full_grid_fn(adata, cfg):
                results = grid_search_clustering(adata, param_grid, ...)
                # rename score → silhouette_score
                for r in results:
                    if "score" in r:
                        r["silhouette_score"] = r.pop("score")
                enrich_grid_results(adata, results, cfg, use_rep=...)
                return results

    log : Logger or None

    Returns
    -------
    dict[str, Any]
        Best entry from full re‑validation (keys: *n_neighbors*, *resolution*,
        *cluster_key*, *n_clusters*, *silhouette_score*, *composite_score*, …).

        The funnel lineage is also stored at ``adata.uns["funnel_lineage"]``.
    """
    _log = log or logger
    use_rep = _resolve_use_rep(adata)

    # -- Step 1: Stratified subsample --
    subsample_size = getattr(cfg.clustering, "funnel_subsample_size", 50000)
    _log.info(
        "Funnel: subsampling (target=%d) via stratified KMeans++ …",
        subsample_size,
    )
    sub_adata, subsample_idx = subsample_stratified(
        adata,
        size=subsample_size,
        use_rep=use_rep,
        log=_log,
    )

    # -- Step 2: Run full grid on subsample --
    _log.info(
        "Funnel: running full grid search on subsample (%d cells) …",
        sub_adata.n_obs,
    )
    sub_results = full_grid_fn(sub_adata, cfg)

    if not sub_results:
        raise ValueError("Funnel: grid search on subsample returned no results")

    # -- Step 3: Rank by composite score, keep top-K --
    _compute_composite_scores(sub_results, cfg, log=_log)

    ranked = sorted(sub_results, key=lambda r: r.get("composite_score", 0.0), reverse=True)
    top_k_count = getattr(cfg.clustering, "funnel_top_k", 3)
    top_k = ranked[:top_k_count]

    _log.info(
        "Funnel: top-%d candidates (subsample):\n  %s",
        len(top_k),
        "\n  ".join(
            f"n_nei={e['n_neighbors']} r={e['resolution']:.2f} "
            f"composite={e.get('composite_score', 0):.4f}"
            for e in top_k
        ),
    )

    # -- Step 4: Re-validate top-K on full data --
    _log.info("Funnel: re‑validating top‑%d on full data …", len(top_k))
    full_results = _validate_on_full(
        adata,
        top_k,
        cfg,
        use_rep=use_rep,
        log=_log,
    )

    if not full_results:
        raise ValueError("Funnel: full re‑validation returned no results")

    # -- Step 5: Pick the best from full validation --
    best = max(full_results, key=lambda r: r.get("composite_score", 0.0))

    # -- Step 6: Record lineage in adata.uns --
    adata.uns["funnel_lineage"] = {
        "subsample_n": int(len(subsample_idx)),
        "subsample_idx": [int(i) for i in subsample_idx],
        "subsample_top_k_composites": [float(e.get("composite_score", 0.0)) for e in top_k],
        "full_top_k_composites": [float(e.get("composite_score", 0.0)) for e in full_results],
        "best_key": str(best.get("cluster_key", "")),
        "composite_delta_full_minus_sub": float(
            best.get("composite_score", 0.0)
            - (top_k[0].get("composite_score", 0.0) if top_k else 0.0)
        ),
    }

    _log.info(
        "Funnel: best from full validation → n_nei=%d r=%.2f composite=%.4f",
        best.get("n_neighbors"),
        best.get("resolution"),
    )

    # Clean up temp columns EXCEPT the best one
    best_key = best.get("cluster_key", "")
    for r in full_results:
        ck = r.get("cluster_key", "")
        if ck != best_key and ck in adata.obs.columns:
            del adata.obs[ck]

    return best
