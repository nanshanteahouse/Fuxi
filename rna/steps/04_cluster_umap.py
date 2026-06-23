#!/usr/bin/env python3
"""
Step 04: 邻居图 + UMAP + 多参数网格 Leiden 聚类
==================================================
  - 在 Harmony 校正后的 PCA 上建图
  - 多参数网格扫描 (n_neighbors × resolution)
  - 保存所有组合结果用于交互比较

输入: 03_integrated.h5ad
输出: 04_grid_results.h5ad (含所有参数组合的邻居图、UMAP、Leiden 标签)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(CFG.log_dir, "04_cluster_umap.log"))
    log.info("Step 04: Neighbors + UMAP + multi-param grid Leiden clustering")

    # ── 输入 ──
    input_path = CFG.integrated_h5ad
    log.info("Loaded: %s", input_path)
    adata = sc.read(input_path)
    log.info("  shape: %s", adata.shape)

    use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    log.info("use_rep: %s", use_rep)

    # ── 参数网格 ──
    n_neighbors_grid = getattr(CFG, 'param_grid_n_neighbors', [15, 20, 30])
    resolutions_grid = getattr(CFG, 'param_grid_resolutions', [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    log.info("Parameter grid: n_neighbors=%s, resolutions=%s", n_neighbors_grid, resolutions_grid)

    results_summary = []

    for n in n_neighbors_grid:
        # 邻居图
        log.info("Computing neighbors (n_neighbors=%d, use_rep=%s)...", n, use_rep)
        try:
            sc.pp.neighbors(
                adata, n_neighbors=n,
                n_pcs=CFG.n_pcs_use, use_rep=use_rep,
                random_state=CFG.random_seed,
            )
        except Exception as e:
            log.error("Neighbor computation failed (n_neighbors=%d): %s", n, e)
            continue

        # ── 计算 UMAP（每个 n_neighbors 只算一次，所有 resolution 共享）──
        log.info("  Computing UMAP (n_neighbors=%d)...", n)
        try:
            sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=CFG.random_seed)
            umap_coords = adata.obsm['X_umap'].copy()
        except Exception as e:
            log.error("  UMAP computation failed (n_neighbors=%d): %s", n, e)
            continue

        # ── 逐 resolution 计算 Leiden + Silhouette（串行，无 deepcopy）──
        log.info("  Serial computation of %d resolutions...", len(resolutions_grid))
        for res in resolutions_grid:
            umap_key = f'umap_{n}_{res}'
            leiden_key = f'leiden_{n}_{res}'

            # 所有 resolution 共享同一份 UMAP 坐标
            adata.obsm[umap_key] = umap_coords

            try:
                sc.tl.leiden(adata, resolution=res, key_added=leiden_key,
                             random_state=CFG.random_seed, flavor=CFG.leiden_flavor)
            except Exception as e:
                log.error("  Leiden failed (n=%d, r=%.1f): %s", n, res, e)
                continue

            leiden_labels = adata.obs[leiden_key]
            n_clusters = int(leiden_labels.nunique())

            # Silhouette score（PCA 空间，大数据集采样 10K 细胞）
            sil_score = None
            try:
                if adata.n_obs > 10000:
                    rng = np.random.RandomState(CFG.random_seed)
                    idx = rng.choice(adata.n_obs, 10000, replace=False)
                    sil_score = float(silhouette_score(
                        adata.obsm[use_rep][idx, :CFG.n_pcs_use],
                        leiden_labels.values[idx],
                    ))
                else:
                    sil_score = float(silhouette_score(
                        adata.obsm[use_rep][:, :CFG.n_pcs_use],
                        leiden_labels.values,
                    ))
            except Exception:
                pass

            score_str = f", silhouette={sil_score:.4f}" if sil_score is not None else ""
            log.info("  n=%d, r=%.1f → %d clusters%s", n, res, n_clusters, score_str)

            results_summary.append({
                'n_neighbors': n,
                'resolution': res,
                'n_clusters': n_clusters,
                'silhouette_score': sil_score,
            })

        # ── 单参数组合 UMAP 图 (逐 resolution，数据已就绪) ──
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
                    os.path.join(CFG.figure_dir,
                                 f'umap_grid_n{n}_r{res}.png'),
                    dpi=150, bbox_inches='tight')
                plt.close()
                log.info("    Plot saved: umap_grid_n%d_r%.1f.png", n, res)
            except Exception as e:
                log.warning("    Single-param UMAP plot save failed: %s", e)
            finally:
                if saved is not None:
                    adata.obsm['X_umap'] = saved

    # ── 汇总 CSV ──
    df_summary = pd.DataFrame(results_summary)
    csv_path = os.path.join(CFG.table_dir, 'param_grid_summary.csv')
    try:
        df_summary.to_csv(csv_path, index=False)
        log.info("Parameter grid summary saved: %s", csv_path)
        log.info("\n%s", df_summary.to_string())
    except Exception as e:
        log.warning("Summary CSV save failed: %s", e)

    if not results_summary:
        log.critical("All neighbor/cluster computations failed — no parameter combination succeeded")
        sys.exit(1)

    # ── 自动选择最佳参数并生成最终 checkpoint ──
    df_summary = pd.DataFrame(results_summary)

    if CFG.best_resolution is not None and any(r['resolution'] == CFG.best_resolution for r in results_summary):
        best = [r for r in results_summary if r['resolution'] == CFG.best_resolution][0]
        log.info("Using configured resolution=%.1f", CFG.best_resolution)
    else:
        # Auto-select: highest silhouette score
        best = max(results_summary, key=lambda r: r['silhouette_score'] or 0)
        log.info("Auto-selected params: n_neighbors=%d, resolution=%.1f (silhouette=%.4f)",
                 best['n_neighbors'], best['resolution'], best['silhouette_score'] or 0)

    best_n = best['n_neighbors']
    best_r = best['resolution']
    leiden_col = f'leiden_{best_n}_{best_r}'
    umap_col = f'umap_{best_n}_{best_r}'

    if leiden_col in adata.obs and umap_col in adata.obsm:
        adata.obs['leiden'] = adata.obs[leiden_col].copy()
        adata.obsm['X_umap'] = adata.obsm[umap_col].copy()
        safe_write(adata, CFG.cluster_h5ad, cfg=CFG)
        log.info("Final checkpoint saved: %s (resolution=%.1f)", CFG.cluster_h5ad, best_r)
    else:
        log.warning("Selected param combination (%s, %s) not in results, skipping auto-lock", leiden_col, umap_col)

    # ── 网格汇总图: 所有参数组合对比 ──
    n_n = len(n_neighbors_grid)
    n_r = len(resolutions_grid)
    try:
        fig, axes = plt.subplots(n_n, n_r,
                                 figsize=(5 * n_r + 2, 4 * n_n + 1),
                                 squeeze=False)
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
                                   legend_fontsize=5,
                                   title=f'n={n}, r={res}')
                    except Exception as e_sub:
                        log.warning("  Subplot failed (n=%d, r=%.1f): %s",
                                    n, res, e_sub)
                        ax.text(0.5, 0.5, 'Error', ha='center',
                                va='center', transform=ax.transAxes)
                    finally:
                        adata.obsm['X_umap'] = saved_umap
                else:
                    ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                            transform=ax.transAxes, fontsize=12)
                    ax.set_title(f'n={n}, r={res}')
        fig.tight_layout()
        fig.savefig(os.path.join(CFG.figure_dir, 'umap_param_grid_summary.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info("Parameter grid summary plot saved")
    except Exception as e:
        log.warning("Grid summary plot generation failed: %s", e)

    # ── 按 n_neighbors 分组的多分辨率对比图 ──
    for n in n_neighbors_grid:
        res_keys = [f'leiden_{n}_{r}' for r in resolutions_grid
                    if f'leiden_{n}_{r}' in adata.obs]
        n_res = len(res_keys)
        if n_res > 0:
            try:
                n_cols = min(3, n_res)
                n_rows = int(np.ceil(n_res / n_cols))
                fig, axes = plt.subplots(n_rows, n_cols,
                                         figsize=(6 * n_cols, 5 * n_rows))
                axes = axes.ravel() if n_res > 1 else [axes]
                for i, key in enumerate(res_keys):
                    sc.pl.umap(adata, color=key, ax=axes[i], show=False,
                               legend_loc='on data', legend_fontsize=6,
                               title=key)
                for j in range(len(res_keys), len(axes)):
                    axes[j].axis('off')
                fig.tight_layout()
                fig.savefig(
                    os.path.join(CFG.figure_dir,
                                 f'umap_leiden_n{n}_all_resolutions.pdf'),
                    dpi=150, bbox_inches='tight')
                plt.close(fig)
                log.info("  Multi-resolution UMAP plot (n=%d) saved", n)
            except Exception as e:
                log.warning("  Multi-resolution comparison plot (n=%d) failed: %s", n, e)

    # ── 保存临时 h5ad (非最终 checkpoint) ──
    temp_path = os.path.join(CFG.h5ad_dir, "04_grid_results.h5ad")
    try:
        safe_write(adata, temp_path, cfg=CFG)
        log.info("Temporary h5ad saved: %s", temp_path)
    except Exception as e:
        log.error("Temporary h5ad save failed: %s", e)

    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
