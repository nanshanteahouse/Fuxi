#!/usr/bin/env python3
"""
Step 03: Multi-param Leiden + UMAP
=====================================
  - Grid search over n_neighbors × resolutions
  - Each KNN graph shared across all resolutions (no redundant recomputation)
  - Silhouette score for quality evaluation (sampled for large datasets)
  - Stores only the selected best combination in obsm / obs

Input:  02_processed.h5ad
Output: 03_clustered.h5ad
"""

import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, validate_adata
import numpy as np
import snapatac2 as snap
from sklearn.metrics import silhouette_score
from joblib import Parallel, delayed
from rna.utils.cluster_evaluation import select_best_params


def _evaluate_n_neighbor_atac(data, n, resolutions, CFG, log):
    """Worker for parallel grid search (ATAC): evaluate one n_neighbors value.

    Runs KNN → UMAP → serial Leiden over all resolutions on a copy.
    Returns (n, results_summary_rows, umap_coords, leiden_cols_dict)
    or None on failure.
    """
    local = data.copy()
    try:
        snap.pp.knn(local, n_neighbors=n)
    except Exception as e:
        log.error("KNN failed (n_neighbors=%d): %s", n, e)
        return None

    try:
        snap.tl.umap(local, random_state=CFG.random_seed)
        umap_coords = local.obsm['X_umap'].copy()
    except Exception as e:
        log.warning("UMAP failed for n=%d, skipping", n)
        return None

    X_spec = local.obsm['X_spectral']
    n_use = min(30, X_spec.shape[1])
    summary_rows = []
    leiden_cols = {}

    for res in resolutions:
        key = f'leiden_{n}_{res}'
        try:
            snap.tl.leiden(local, resolution=res, key_added=key,
                           random_state=CFG.random_seed)
            n_cl = int(local.obs[key].nunique())
            sil = None
            try:
                if local.n_obs > 10000:
                    rng = np.random.RandomState(CFG.random_seed)
                    idx = rng.choice(local.n_obs, 10000, replace=False)
                    sil = float(silhouette_score(X_spec[idx, :n_use],
                                                 local.obs[key].values[idx]))
                else:
                    sil = float(silhouette_score(X_spec[:, :n_use],
                                                 local.obs[key].values))
            except Exception:
                pass
            sil_str = f", sil={sil:.4f}" if sil is not None else ""
            log.info("  n=%d r=%.1f -> %d clusters%s", n, res, n_cl, sil_str)
            summary_rows.append({
                'n_neighbors': n, 'resolution': res,
                'n_clusters': n_cl, 'silhouette_score': sil,
            })
            leiden_cols[key] = local.obs[key].values.copy()
        except Exception as e:
            log.warning("  Leiden failed (n=%d, r=%.1f): %s", n, res, e)

    return (n, summary_rows, umap_coords, leiden_cols)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("03_cluster", os.path.join(CFG.log_dir, "03_cluster.log"))
    log.info("Step 03: Multi-param Leiden + UMAP")

    if os.path.exists(CFG.clustered_h5ad):
        log.info("Skip: %s exists.", CFG.clustered_h5ad)
        return

    # Read in backed mode then materialize to memory
    data = snap.read(CFG.processed_h5ad)
    if data.isbacked:
        data = data.to_memory()
    log.info("Loaded: %d cells, vars: %d", data.n_obs, data.n_vars)

    # Clean up any stray leiden columns from previous partial runs
    for col in list(data.obs.columns):
        if col.startswith('leiden_'):
            del data.obs[col]

    nns = getattr(CFG, 'param_grid_n_neighbors', [15, 20, 30])
    resolutions = getattr(CFG, 'param_grid_resolutions', [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])

    results_summary = []

    # ── Parallel outer loop over n_neighbors ──
    n_jobs = min(getattr(CFG, 'n_jobs', 4) or os.cpu_count() or 1, len(nns))
    log.info("Evaluating %d n_neighbors values with n_jobs=%d", len(nns), n_jobs)
    parallel_results = Parallel(n_jobs=n_jobs, prefer='threads')(
        delayed(_evaluate_n_neighbor_atac)(data, n, resolutions, CFG, log)
        for n in nns
    )

    # ── Collect results back into main AnnData ──
    for r in parallel_results:
        if r is None:
            continue
        n, summary_rows, umap_coords, leiden_cols = r
        results_summary.extend(summary_rows)
        # Store UMAP coords per n_neighbors (overwrites for last, OK — kept only for grid summary)
        data.obsm[f'X_umap_{n}'] = umap_coords
        for key, labels in leiden_cols.items():
            data.obs[key] = labels

    if not results_summary:
        log.critical("All parameter combinations failed.")
        sys.exit(1)

    # ── Keep only the best combination in obsm/obs (reduces saved file size) ──
    method = getattr(CFG, 'cluster_selection_method', 'pareto_elbow')

    if method is not None and (getattr(CFG, 'best_resolution', 1.0) != 1.0 or getattr(CFG, 'best_n_neighbors', 0) != 0):
        log.warning(
            "best_resolution=%.1f / best_n_neighbors=%d are set but cluster_selection_method=%r will ignore them. "
            "Set cluster_selection_method=None to use manual mode.",
            CFG.best_resolution, getattr(CFG, 'best_n_neighbors', 0), method,
        )

    best_n, best_r, method_name, reason = select_best_params(
        results_summary,
        method=method,
        best_resolution=CFG.best_resolution if method is None else None,
        best_n_neighbors=getattr(CFG, 'best_n_neighbors', 0) if method is None else 0,
    )

    log.info("Selected best params via %s: n_neighbors=%d, resolution=%.1f (%s)",
             method_name, best_n, best_r, reason)

    best_key = f"leiden_{best_n}_{best_r}"
    best_umap_key = f"X_umap_{best_n}"
    if best_key in data.obs:
        data.obs['leiden'] = data.obs[best_key]
        # Set X_umap from the stored per-n copy
        if best_umap_key in data.obsm:
            data.obsm['X_umap'] = data.obsm[best_umap_key]
        # Clean up grid search columns — keep only the best
        for col in list(data.obs.columns):
            if col.startswith('leiden_') and col != 'leiden':
                del data.obs[col]
        # Clean up per-n UMAP keys
        for key in list(data.obsm.keys()):
            if key.startswith('X_umap_') and key != 'X_umap':
                del data.obsm[key]
        # Keep only the best UMAP in obsm
        data.uns['cluster_params'] = {
            'n_neighbors': best_n,
            'resolution': best_r,
            'method': method_name,
            'reason': reason,
        }

    safe_write(data, CFG.clustered_h5ad, cfg=CFG, compression_override=None)
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
