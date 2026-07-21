#!/usr/bin/env python3
"""
Step 04: Multi-param Leiden + UMAP
=====================================
  - Grid search over n_neighbors × resolutions
  - Each KNN graph shared across all resolutions (no redundant recomputation)
  - Silhouette score for quality evaluation (sampled for large datasets)
  - Stores only the selected best combination in obsm / obs

Input:  03_processed.h5ad
Output: 04_clustered.h5ad
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import snapatac2 as snap
from joblib import Parallel, delayed
from sklearn.metrics import silhouette_score

from core.cluster.grid_search import grid_search_clustering, select_best_params
from core.config.schema import SILHOUETTE_SAMPLE_THRESHOLD
from core.utils import resolve_config, safe_write, setup_logger

# ---------------------------------------------------------------------------
#  ATAC-specific callables for the shared grid_search_clustering interface
# ---------------------------------------------------------------------------


def _atac_clusterer(adata, resolution=None, n_neighbors=None, random_seed=42, **kwargs):
    """Run snapatac2 leiden and return the observation column key.

    The key is constructed from n_neighbors and resolution so each
    combination gets a unique column (e.g. ``leiden_15_0.5``).
    """
    assert n_neighbors is not None, "n_neighbors must be provided"
    assert resolution is not None, "resolution must be provided"
    key = f"leiden_{n_neighbors}_{resolution}"
    snap.tl.leiden(adata, resolution=resolution, key_added=key, random_state=random_seed)
    return key


def _atac_neighbor_fn(adata, n_neighbors=None, **kwargs):
    """Compute the KNN graph with snapatac2."""
    assert n_neighbors is not None, "n_neighbors must be provided"
    snap.pp.knn(adata, n_neighbors=n_neighbors)


def _atac_umap_fn(adata, random_seed=42, **kwargs):
    """Compute the UMAP embedding with snapatac2."""
    snap.tl.umap(adata, random_state=random_seed)


def _atac_evaluation_fn(adata, cluster_key, random_seed=42, **kwargs):
    """Silhouette score on the spectral embedding (sampled for large datasets)."""
    x_spec = adata.obsm["X_spectral"]
    n_use = min(30, x_spec.shape[1])
    if adata.n_obs > SILHOUETTE_SAMPLE_THRESHOLD:
        rng = np.random.RandomState(random_seed)
        idx = rng.choice(adata.n_obs, SILHOUETTE_SAMPLE_THRESHOLD, replace=False)
        return float(silhouette_score(x_spec[idx, :n_use], adata.obs[cluster_key].values[idx]))
    return float(silhouette_score(x_spec[:, :n_use], adata.obs[cluster_key].values))


def _evaluate_n_neighbor_atac(data, n, resolutions, cfg, log):
    """Worker for parallel grid search (ATAC): evaluate one n_neighbors value.

    Delegates the inner grid (KNN → resolutions × Leiden + silhouette) to the
    shared ``grid_search_clustering`` on an in-memory copy of the data.

    Returns (n, results_summary_rows, umap_coords, leiden_cols_dict)
    or None on failure.
    """
    local = data.copy()
    try:
        results = grid_search_clustering(
            local,
            param_grid={"n_neighbors": [n], "resolution": list(resolutions)},
            clusterer=_atac_clusterer,
            neighbor_fn=_atac_neighbor_fn,
            umap_fn=_atac_umap_fn,
            evaluation_fn=_atac_evaluation_fn,
            group_key="n_neighbors",
            random_seed=cfg.execution.random_seed,
        )
    except Exception as e:
        log.error("Grid search failed (n_neighbors=%d): %s", n, e)
        return None

    if not results:
        return None

    # Extract UMAP coordinates from the copy
    umap_coords = local.obsm.get("X_umap").copy() if "X_umap" in local.obsm else None

    # Extract leiden labels from the copy before it goes out of scope
    leiden_cols = {}
    for r in results:
        ck = r.get("cluster_key")
        if ck and ck in local.obs:
            leiden_cols[ck] = local.obs[ck].values.copy()

    # Build summary rows in the format expected by select_best_params
    summary_rows = []
    for r in results:
        sil = r.get("score")
        sil_str = f", sil={sil:.4f}" if sil is not None else ""
        log.info(
            "  n=%d r=%.1f -> %d clusters%s",
            r["n_neighbors"],
            r["resolution"],
            r["n_clusters"],
            sil_str,
        )
        summary_rows.append(
            {
                "n_neighbors": r["n_neighbors"],
                "resolution": r["resolution"],
                "n_clusters": r["n_clusters"],
                "silhouette_score": sil,
            }
        )

    return (n, summary_rows, umap_coords, leiden_cols)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(cfg.log_dir, "04_cluster.log"))
    log.info("Step 04: Multi-param Leiden + UMAP")

    if os.path.exists(cfg.clustered_h5ad):
        log.info("Skip: %s exists.", cfg.clustered_h5ad)
        return

    # Read in backed mode then materialize to memory
    data = snap.read(cfg.processed_h5ad)
    if data.isbacked:
        data = data.to_memory()
    log.info("Loaded: %d cells, vars: %d", data.n_obs, data.n_vars)

    # Clean up any stray leiden columns from previous partial runs
    for col in list(data.obs.columns):
        if col.startswith("leiden_"):
            del data.obs[col]

    nns = getattr(cfg.clustering, "param_grid_n_neighbors", [15, 20, 30])
    resolutions = getattr(cfg.clustering, "param_grid_resolutions", [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])

    results_summary = []

    # ── Parallel outer loop over n_neighbors ──
    n_jobs = min(getattr(cfg.execution, "n_jobs", 4) or os.cpu_count() or 1, len(nns))
    log.info("Evaluating %d n_neighbors values with n_jobs=%d", len(nns), n_jobs)
    parallel_results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_evaluate_n_neighbor_atac)(data, n, resolutions, cfg, log) for n in nns
    )

    # ── Collect results back into main AnnData ──
    for r in parallel_results:
        if r is None:
            continue
        n, summary_rows, umap_coords, leiden_cols = r
        results_summary.extend(summary_rows)
        # Store UMAP coords per n_neighbors (overwrites for last, OK — kept only for grid summary)
        if umap_coords is not None:
            data.obsm[f"X_umap_{n}"] = umap_coords
        for key, labels in leiden_cols.items():
            data.obs[key] = labels

    if not results_summary:
        log.critical("All parameter combinations failed.")
        sys.exit(1)

    # ── Keep only the best combination in obsm/obs (reduces saved file size) ──
    method = getattr(cfg.clustering, "cluster_selection_method", "pareto_elbow")

    if method is not None and (
        getattr(cfg.clustering, "best_resolution", 1.0) != 1.0
        or getattr(cfg.clustering, "best_n_neighbors", 0) != 0
    ):
        log.warning(
            "best_resolution=%.1f / best_n_neighbors=%d are set but cluster_selection_method=%r will ignore them. "
            "Set cluster_selection_method=None to use manual mode.",
            cfg.clustering.best_resolution,
            getattr(cfg.clustering, "best_n_neighbors", 0),
            method,
        )

    best_n, best_r, method_name, reason = select_best_params(
        results_summary,
        method=method,  # type: ignore[reportArgumentType]
        best_resolution=cfg.clustering.best_resolution if method is None else None,
        best_n_neighbors=getattr(cfg.clustering, "best_n_neighbors", 0) if method is None else 0,
        log=log,
    )

    log.info(
        "Selected best params via %s: n_neighbors=%d, resolution=%.1f (%s)",
        method_name,
        best_n,
        best_r,
        reason,
    )

    best_key = f"leiden_{best_n}_{best_r}"
    best_umap_key = f"X_umap_{best_n}"
    if best_key in data.obs:
        data.obs["leiden"] = data.obs[best_key]
        # Set X_umap from the stored per-n copy
        if best_umap_key in data.obsm:
            data.obsm["X_umap"] = data.obsm[best_umap_key]
        # Clean up grid search columns — keep only the best
        for col in list(data.obs.columns):
            if col.startswith("leiden_") and col != "leiden":
                del data.obs[col]
        # Clean up per-n UMAP keys
        for key in list(data.obsm.keys()):
            if key.startswith("X_umap_") and key != "X_umap":
                del data.obsm[key]
        # Keep only the best UMAP in obsm
        data.uns["cluster_params"] = {
            "n_neighbors": best_n,
            "resolution": best_r,
            "method": method_name,
            "reason": reason,
        }

    safe_write(data, cfg.clustered_h5ad, cfg=cfg, compression_override=None)
    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
