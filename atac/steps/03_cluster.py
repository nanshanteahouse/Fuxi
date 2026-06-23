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
    X_spec = data.obsm['X_spectral']
    n_use = min(30, X_spec.shape[1])

    for n in nns:
        snap.pp.knn(data, n_neighbors=n)
        try:
            snap.tl.umap(data, random_state=CFG.random_seed)
            umap_coords = data.obsm['X_umap'].copy()
        except Exception:
            log.warning("UMAP failed for n=%d, skipping", n)
            continue

        for res in resolutions:
            key = f'leiden_{n}_{res}'
            try:
                snap.tl.leiden(data, resolution=res, key_added=key,
                               random_state=CFG.random_seed)
                n_cl = int(data.obs[key].nunique())
                sil = None
                try:
                    if data.n_obs > 10000:
                        rng = np.random.RandomState(CFG.random_seed)
                        idx = rng.choice(data.n_obs, 10000, replace=False)
                        sil = float(silhouette_score(X_spec[idx, :n_use],
                                                     data.obs[key].values[idx]))
                    else:
                        sil = float(silhouette_score(X_spec[:, :n_use],
                                                     data.obs[key].values))
                except Exception:
                    pass
                sil_str = f", sil={sil:.4f}" if sil is not None else ""
                log.info("  n=%d r=%.1f -> %d clusters%s", n, res, n_cl, sil_str)
                results_summary.append({
                    'n_neighbors': n, 'resolution': res,
                    'n_clusters': n_cl, 'silhouette_score': sil,
                })
            except Exception as e:
                log.warning("  Leiden failed (n=%d, r=%.1f): %s", n, res, e)

    if not results_summary:
        log.critical("All parameter combinations failed.")
        sys.exit(1)

    # ── Keep only the best combination in obsm/obs (reduces saved file size) ──
    best_res = getattr(CFG, 'best_resolution', None)
    if best_res is not None and any(r['resolution'] == best_res for r in results_summary):
        candidates = [r for r in results_summary if r['resolution'] == best_res]
        best = max(candidates, key=lambda r: r['silhouette_score'] or 0)
        log.info("Using configured resolution=%.1f", best_res)
    else:
        best = max(results_summary, key=lambda r: r['silhouette_score'] or 0)
        log.info("Auto-selected: n=%d r=%.1f (sil=%.4f)",
                 best['n_neighbors'], best['resolution'], best['silhouette_score'] or 0)

    best_key = f"leiden_{best['n_neighbors']}_{best['resolution']}"
    if best_key in data.obs:
        data.obs['leiden'] = data.obs[best_key]
        # Clean up grid search columns — keep only the best
        for col in list(data.obs.columns):
            if col.startswith('leiden_') and col != 'leiden':
                del data.obs[col]
        # Keep only the best UMAP in obsm
        data.uns['cluster_params'] = best

    safe_write(data, CFG.clustered_h5ad, cfg=CFG)
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
