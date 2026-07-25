#!/usr/bin/env python3
"""
Binary-search resolution targeting a desired cluster count.

Provides ``find_resolution_for_target_k`` which uses binary search over
resolution (with 3-run median to handle Leiden non-monotonicity) to find a
resolution that produces approximately *target_k* clusters.

Exports
-------
find_resolution_for_target_k(adata, target_k, n_neighbors, cfg, log=None)
    -> (resolution, actual_k)
target_grid_search(adata, cfg, log=None)
    -> list[dict]  (enriched results_summary)
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------


def find_resolution_for_target_k(
    adata: Any,
    target_k: int,
    n_neighbors: int,
    cfg: Any,
    log: logging.Logger | None = None,
) -> tuple[float, int]:
    """Binary-search Leiden resolution to hit a target cluster count.

    Builds the KNN graph once at *n_neighbors*, then performs a binary search
    over resolution in ``[0.1, 5.0]``.  At each candidate resolution Leiden is
    run 3 times with different seeds (42, 123, 456) and the median cluster
    count is used — this guards against non-monotonicity from seed variance.

    Convergence is reached when ``|actual_k - target_k| ≤ 1`` or when
    *max_iters* is exhausted.  On non-convergence the entry closest to
    *target_k* is returned with a warning.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.  Must already have a PCA or other embedding
        in ``adata.obsm`` referenced by ``cfg.clustering.use_rep``.
    target_k : int
        Desired number of clusters.
    n_neighbors : int
        Number of neighbours for the KNN graph.
    cfg : Config
        Pipeline config object.  Reads:

        - ``cfg.clustering.use_rep`` (str, default ``'X_pca'``)
        - ``cfg.clustering.target_search_max_iters`` (int, default ``10``)
    log : Logger or None
        Optional logger.  Uses module-level logger when *None*.

    Returns
    -------
    resolution : float
        Best resolution found (closest to *target_k*).
    actual_k : int
        Number of clusters at that resolution (median of 3 seeds).
    """
    _log = log or logger
    lo, hi = 0.1, 5.0

    best_res: float = lo
    best_k: int = 0
    best_dist: int | float = float("inf")
    convergence_curve: list[tuple[float, int, list[int]]] = []

    # Resolve config values with safe fallbacks
    max_iters: int = getattr(cfg.clustering, "target_search_max_iters", 10)
    use_rep: str = getattr(cfg.clustering, "use_rep", "X_pca")

    # --- Build KNN once for this n_neighbors ---
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    for iteration in range(max_iters):
        mid = (lo + hi) / 2.0

        # 3-run median for stability against seed variance
        ks: list[int] = []
        for seed in (42, 123, 456):
            sc.tl.leiden(
                adata,
                resolution=mid,
                random_state=seed,
                key_added="_target_tmp",
            )
            ks.append(int(adata.obs["_target_tmp"].nunique()))
            del adata.obs["_target_tmp"]

        actual_k = int(np.median(ks))
        convergence_curve.append((mid, actual_k, ks))

        dist = abs(actual_k - target_k)
        if dist < best_dist:
            best_res = mid
            best_k = actual_k
            best_dist = dist

        if dist <= 1:
            _log.info(
                "Target converged: resolution=%.3f -> k=%d (seeds=%s)",
                mid,
                actual_k,
                ks,
            )
            return (mid, actual_k)

        # Binary-search direction:
        #   More clusters → lower resolution to reduce clusters
        #   Fewer clusters → higher resolution to increase clusters
        if actual_k > target_k:
            hi = mid
        else:
            lo = mid

    # --- Non-convergence ---
    _log.warning(
        "Target search did not converge to k=%d in %d iters. "
        "Best: k=%d at res=%.3f (dist=%d). "
        "Curve: %s",
        target_k,
        max_iters,
        best_k,
        best_res,
        best_dist,
        convergence_curve,
    )
    return (best_res, best_k)


def target_grid_search(
    adata: Any,
    cfg: Any,
    log: logging.Logger | None = None,
) -> list[dict[str, Any]]:
    """Grid search over *n_neighbors* using target-K binary-search resolution.

    For each value in ``cfg.clustering.param_grid_n_neighbors``, calls
    :func:`find_resolution_for_target_k` to find a resolution that hits
    ``cfg.clustering.target_n_clusters`` clusters.  The resulting
    (n_neighbors, resolution) pairs are then enriched with multi-metric
    scores via :func:`core.cluster.evaluation.enrich_grid_results`.

    Parameters
    ----------
    adata : AnnData
    cfg : Config
        Pipeline config.  Reads:

        - ``cfg.clustering.param_grid_n_neighbors`` (list[int])
        - ``cfg.clustering.target_n_clusters`` (int)
    log : Logger or None

    Returns
    -------
    list[dict]
        Enriched results_summary with the same structure as
        :func:`core.cluster.evaluation.enrich_grid_results` output.
        Each entry includes *n_neighbors*, *resolution*, *n_clusters*,
        *cluster_key*, and all multi-metric scores.
    """
    _log = log or logger

    param_grid_n_neighbors: list[int] = getattr(
        cfg.clustering, "param_grid_n_neighbors", [15, 20, 30]
    )
    target_k: int = getattr(cfg.clustering, "target_n_clusters", 20)

    # --- De-duplicate neighbours (in case param_grid_n_neighbors has duplicates) ---
    n_vals = sorted(set(param_grid_n_neighbors))
    _log.info(
        "Target grid search: target_k=%d, n_neighbors=%s",
        target_k,
        n_vals,
    )

    results_summary: list[dict[str, Any]] = []
    for n_val in n_vals:
        resolution, actual_k = find_resolution_for_target_k(
            adata,
            target_k=target_k,
            n_neighbors=n_val,
            cfg=cfg,
            log=_log,
        )
        cluster_key = f"leiden_{n_val}_{resolution:.3f}"
        results_summary.append(
            {
                "n_neighbors": n_val,
                "resolution": resolution,
                "n_clusters": actual_k,
                "cluster_key": cluster_key,
            }
        )

    # --- Enrich with multi-metric scores ---
    # Lazy import to avoid circular dependency
    from core.cluster.evaluation import enrich_grid_results

    use_rep: str = getattr(cfg.clustering, "use_rep", "X_pca")
    enrich_grid_results(
        adata,
        results_summary,
        cfg=cfg,
        log=_log,
        use_rep=use_rep,
    )

    return results_summary
