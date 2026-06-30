#!/usr/bin/env python3
"""
Step 04: Neighbors + UMAP + Leiden clustering
================================================
  - Build PCA neighbor graph
  - Multi-param grid scan (n_neighbors × resolution)
  - Auto-select best params (silhouette score)
  - Generate UMAP visualizations

Input:  03_processed.h5ad
Output: 04_clustered.h5ad
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from rna.utils.cluster_evaluation import select_best_params, select_best_umap_params


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(CFG.log_dir, "04_cluster.log"))
    log.info("Step 04: Neighbors + UMAP + Leiden clustering")

    # ── Checkpoint ──────────────────────────────────────────────────────'
    output_path = os.path.join(CFG.h5ad_dir, "04_clustered.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    # ── Input ───────────────────────────────────────────────────────────
    input_path = os.path.join(CFG.h5ad_dir, "03_processed.h5ad")
    if not os.path.exists(input_path):
        log.error("Input not found: %s", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    use_rep = 'X_pca'
    if 'X_pca' not in adata.obsm:
        log.error("No PCA found in obsm. Run Step 03 first.")
        sys.exit(1)

    log.info("Using PCA representation: %s (%d PCs)", use_rep, CFG.n_pcs_use)

    # ── Parameter grid ──────────────────────────────────────────────────
    n_neighbors_grid = getattr(CFG, 'param_grid_n_neighbors', [15, 20, 30])
    resolutions_grid = getattr(CFG, 'param_grid_resolutions', [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    log.info("Grid: n_neighbors=%s, resolutions=%s", n_neighbors_grid, resolutions_grid)

    fig_dir = os.path.join(CFG.figure_dir, '04_cluster')
    os.makedirs(fig_dir, exist_ok=True)

    results_summary = []

    for n in n_neighbors_grid:
        # ── Neighbor graph ──
        log.info("Computing neighbors (n_neighbors=%d)...", n)
        try:
            sc.pp.neighbors(
                adata, n_neighbors=n,
                n_pcs=CFG.n_pcs_use, use_rep=use_rep,
                random_state=CFG.random_seed,
            )
        except Exception as e:
            log.error("Neighbor computation failed (n_neighbors=%d): %s", n, e)
            continue

        # ── UMAP (once per n_neighbors) ──
        log.info("  Computing UMAP...")
        try:
            sc.tl.umap(adata, min_dist=0.3, spread=1.0,
                       random_state=CFG.random_seed)
            umap_coords = adata.obsm['X_umap'].copy()
        except Exception as e:
            log.error("  UMAP failed: %s", e)
            continue

        # ── Leiden per resolution ──
        log.info("  Running %d resolutions...", len(resolutions_grid))
        for res in resolutions_grid:
            umap_key = f'umap_{n}_{res}'
            leiden_key = f'leiden_{n}_{res}'

            adata.obsm[umap_key] = umap_coords

            try:
                sc.tl.leiden(
                    adata, resolution=res, key_added=leiden_key,
                    random_state=CFG.random_seed,
                    flavor=getattr(CFG, 'leiden_flavor', 'igraph'),
                )
            except Exception as e:
                log.error("  Leiden failed (n=%d, r=%.1f): %s", n, res, e)
                continue

            n_clusters = int(adata.obs[leiden_key].nunique())

            # Silhouette score
            sil_score = None
            try:
                labels = adata.obs[leiden_key].values
                if adata.n_obs > 10000:
                    rng = np.random.RandomState(CFG.random_seed)
                    idx = rng.choice(adata.n_obs, 10000, replace=False)
                    sil_score = float(silhouette_score(
                        adata.obsm[use_rep][idx, :CFG.n_pcs_use],
                        labels[idx],
                    ))
                else:
                    sil_score = float(silhouette_score(
                        adata.obsm[use_rep][:, :CFG.n_pcs_use],
                        labels,
                    ))
            except Exception:
                pass

            score_str = f", silhouette={sil_score:.4f}" if sil_score is not None else ""
            log.info("  n=%d, r=%.1f -> %d clusters%s", n, res, n_clusters, score_str)

            results_summary.append({
                'n_neighbors': n,
                'resolution': res,
                'n_clusters': n_clusters,
                'silhouette_score': sil_score,
            })

        # ── Generate per-param UMAP plots ──
        for res in resolutions_grid:
            umap_key = f'umap_{n}_{res}'
            leiden_key = f'leiden_{n}_{res}'
            if umap_key not in adata.obsm or leiden_key not in adata.obs:
                continue
            saved = adata.obsm.get('X_umap')
            adata.obsm['X_umap'] = adata.obsm[umap_key].copy()
            try:
                safe_plot(sc.pl.umap, adata, color=leiden_key, show=False,
                          title=f'UMAP (n_neighbors={n}, resolution={res})')
                plt.savefig(
                    os.path.join(fig_dir,
                                 f'umap_grid_n{n}_r{res}.png'),
                    dpi=150, bbox_inches='tight')
                plt.close()
            except Exception as e:
                log.warning("    Plot save failed: %s", e)
            finally:
                if saved is not None:
                    adata.obsm['X_umap'] = saved

    if not results_summary:
        log.critical("All clustering computations failed")
        sys.exit(1)

    # ── Save parameter grid summary ──
    df_summary = pd.DataFrame(results_summary)
    csv_path = os.path.join(CFG.table_dir, 'param_grid_summary.csv')
    df_summary.to_csv(csv_path, index=False)
    log.info("Grid summary saved: %s", csv_path)

    # ── 自动选择最佳参数并生成最终 checkpoint ──
    df_summary = pd.DataFrame(results_summary)

    if not results_summary:
        log.critical("All neighbor/cluster computations failed — no parameter combination succeeded")
        sys.exit(1)

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
    leiden_col = f'leiden_{best_n}_{best_r}'
    umap_col = f'umap_{best_n}_{best_r}'

    if leiden_col in adata.obs and umap_col in adata.obsm:
        adata.obs['leiden'] = adata.obs[leiden_col].copy()
        adata.obsm['X_umap'] = adata.obsm[umap_col].copy()
        safe_write(adata, output_path, cfg=CFG)
        log.info("Final checkpoint saved: %s (resolution=%.1f)", output_path, best_r)
    else:
        log.warning("Best param combination not found in results, skipping auto-lock")

    # ── UMAP 参数扫描 (min_dist × spread) ──────────────────────────────────
    min_dist_grid = getattr(CFG, 'param_grid_min_dist', [0.3])
    spread_grid = getattr(CFG, 'param_grid_spread', [1.0])
    umap_method = getattr(CFG, 'umap_selection_method', 'convex_hull')
    best_md, best_sp, umap_method_label, sweep_results = select_best_umap_params(
        adata, best_n, min_dist_grid, spread_grid, umap_method, CFG, use_rep, log)

    # Rebuild UMAP with selected params and re-save checkpoint
    log.info("Rebuilding UMAP with selected params (min_dist=%.2f, spread=%.1f) [%s]...",
             best_md, best_sp, umap_method_label)
    try:
        sc.tl.umap(adata, min_dist=best_md, spread=best_sp,
                   random_state=CFG.random_seed)
        safe_write(adata, output_path, cfg=CFG)
        log.info("Checkpoint saved with final UMAP: %s", output_path)
    except Exception as e:
        log.warning("Final UMAP rebuild failed: %s", e)

    if sweep_results:
        # Summary CSV
        try:
            sweep_csv = os.path.join(CFG.table_dir, 'umap_min_dist_sweep_summary.csv')
            pd.DataFrame(sweep_results).to_csv(sweep_csv, index=False)
            log.info("UMAP sweep summary saved: %s", sweep_csv)
        except Exception as e:
            log.warning("UMAP sweep CSV save failed: %s", e)

        # Comparison figure
        try:
            n_md = len(min_dist_grid) if min_dist_grid else 1
            n_sp = len(spread_grid) if spread_grid else 1
            n_total = n_md * n_sp
            n_cols = min(3, n_total)
            n_rows = int(np.ceil(n_total / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols,
                                     figsize=(6 * n_cols, 5 * n_rows))
            axes_flat = axes.ravel() if n_total > 1 else [axes]
            for idx, r in enumerate(sweep_results):
                ax = axes_flat[idx]
                md = r['min_dist']
                sp = r['spread']
                try:
                    sc.tl.umap(adata, min_dist=md, spread=sp,
                               random_state=CFG.random_seed)
                    sc.pl.umap(adata, color='leiden', ax=ax, show=False,
                               legend_loc='on data', legend_fontsize=5,
                               title=f'min_dist={md}, spread={sp}')
                except Exception:
                    ax.text(0.5, 0.5, 'Error', ha='center', va='center',
                            transform=ax.transAxes)
            for j in range(len(sweep_results), len(axes_flat)):
                axes_flat[j].axis('off')
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, 'umap_min_dist_comparison.png'),
                        dpi=150, bbox_inches='tight')
            plt.close(fig)
            log.info("UMAP min_dist comparison plot saved")
        except Exception as e:
            log.warning("UMAP comparison plot generation failed: %s", e)

    # ── Grid summary plot ──
    try:
        n_n = len(n_neighbors_grid)
        n_r = len(resolutions_grid)
        fig, axes = plt.subplots(n_n, n_r, figsize=(5 * n_r + 2, 4 * n_n + 1), squeeze=False)
        for i, n in enumerate(n_neighbors_grid):
            for j, res in enumerate(resolutions_grid):
                ax = axes[i, j]
                umap_key = f'umap_{n}_{res}'
                leiden_key = f'leiden_{n}_{res}'
                if umap_key in adata.obsm and leiden_key in adata.obs:
                    saved_umap = adata.obsm['X_umap'].copy()
                    try:
                        adata.obsm['X_umap'] = adata.obsm[umap_key].copy()
                        sc.pl.umap(adata, color=leiden_key, ax=ax,
                                   show=False, legend_loc='on data',
                                   legend_fontsize=5, title=f'n={n}, r={res}')
                    except Exception:
                        ax.text(0.5, 0.5, 'Error', ha='center', va='center',
                                transform=ax.transAxes)
                    finally:
                        adata.obsm['X_umap'] = saved_umap
                else:
                    ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                            transform=ax.transAxes, fontsize=12)
                    ax.set_title(f'n={n}, r={res}')
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'umap_param_grid_summary.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        log.warning("Grid summary plot failed: %s", e)

    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
