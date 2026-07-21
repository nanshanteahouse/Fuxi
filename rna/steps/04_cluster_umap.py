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

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score

from core.cluster.evaluation import select_best_umap_params
from core.cluster.grid_search import grid_search_clustering, select_best_params
from core.config.schema import SILHOUETTE_SAMPLE_THRESHOLD
from core.utils import resolve_config, safe_plot, safe_write, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(cfg.log_dir, "04_cluster_umap.log"))
    log.info("Step 04: Neighbors + UMAP + multi-param grid Leiden clustering")

    # ── 输入 ──
    input_path = cfg.integrated_h5ad
    log.info("Loaded: %s", input_path)
    adata = sc.read(input_path)
    log.info("  shape: %s", adata.shape)

    use_rep = "X_integrated" if "X_integrated" in adata.obsm else "X_pca"
    log.info("use_rep: %s", use_rep)

    # ── 参数网格 ──
    n_neighbors_grid = getattr(cfg.clustering, "param_grid_n_neighbors", [15, 20, 30])
    resolutions_grid = getattr(
        cfg.clustering, "param_grid_resolutions", [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    )
    # ── Adaptive resolution expansion for small datasets ──
    if getattr(cfg.clustering, "multi_metric_adaptive_resolution", False):
        if adata.n_obs < 3000:
            extra = [r for r in [3.0, 5.0] if r not in resolutions_grid]
            if extra:
                resolutions_grid = sorted(resolutions_grid + extra)
                log.info(
                    "Adaptive resolution expansion (n_cells=%d < 3000): added %s → %s",
                    adata.n_obs,
                    extra,
                    resolutions_grid,
                )
    log.info("Parameter grid: n_neighbors=%s, resolutions=%s", n_neighbors_grid, resolutions_grid)

    fig_dir = os.path.join(cfg.figure_dir, "04_cluster")
    os.makedirs(fig_dir, exist_ok=True)

    umap_min_dist = getattr(cfg.clustering, "umap_min_dist", 0.3)
    umap_spread = getattr(cfg.clustering, "umap_spread", 1.0)

    # ── Grid search via shared core ──

    # Scanpy-specific callables (close over CFG and use_rep)
    def _neighbors_fn(adata, n_neighbors=15, **kwargs):
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=cfg.pca.n_pcs_use,
            use_rep=use_rep,
            random_state=cfg.execution.random_seed,
        )

    def _umap_fn(adata, **kwargs):
        sc.tl.umap(
            adata,
            min_dist=umap_min_dist,
            spread=umap_spread,
            random_state=cfg.execution.random_seed,
        )

    def _clusterer_fn(adata, resolution=1.0, n_neighbors=15, **kwargs):
        leiden_key = f"leiden_{n_neighbors}_{resolution}"
        umap_key = f"umap_{n_neighbors}_{resolution}"
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=leiden_key,
            random_state=cfg.execution.random_seed,
            flavor=cfg.clustering.leiden_flavor,
            directed=False,
            n_iterations=2,
        )
        adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        return leiden_key

    def _evaluation_fn(adata, cluster_key, **kwargs):
        labels = adata.obs[cluster_key].values
        if adata.n_obs > SILHOUETTE_SAMPLE_THRESHOLD:
            rng = np.random.RandomState(cfg.execution.random_seed)
            idx = rng.choice(adata.n_obs, SILHOUETTE_SAMPLE_THRESHOLD, replace=False)
            return float(
                silhouette_score(
                    adata.obsm[use_rep][idx, : cfg.pca.n_pcs_use],
                    labels[idx],
                )
            )
        else:
            return float(
                silhouette_score(
                    adata.obsm[use_rep][:, : cfg.pca.n_pcs_use],
                    labels,
                )
            )

    results_summary = grid_search_clustering(
        adata,
        param_grid={
            "n_neighbors": n_neighbors_grid,
            "resolution": resolutions_grid,
        },
        clusterer=_clusterer_fn,
        neighbor_fn=_neighbors_fn,
        umap_fn=_umap_fn,
        evaluation_fn=_evaluation_fn,
        group_key="n_neighbors",
        n_jobs=cfg.execution.n_jobs,
        random_seed=cfg.execution.random_seed,
    )

    # Rename score → silhouette_score for select_best_params compatibility
    for r in results_summary:
        if "score" in r:
            r["silhouette_score"] = r.pop("score")
    # ── Granularity detection ──
    from core.cluster.evaluation import _detect_granularity

    granularity = _detect_granularity(results_summary)
    log.info("Granularity classification: %s", granularity)
    _de_gated_selected = False

    if granularity == "subtype":
        _de_gated_selected = True
        from core.cluster.evaluation import _select_de_gated

        log.info(
            "Granularity=subtype — using DE-gated resolution selection (bypassing enrichment)"
        )
        n_clusters, resolution, cluster_key, reason_str = _select_de_gated(
            results_summary,
            adata,
            de_gate_threshold=getattr(cfg.clustering, "multi_metric_de_gate_threshold", 25),
        )
        # Extract n_neighbors from selected entry
        best_n = None
        for r in results_summary:
            if r.get("resolution") == resolution and r.get("cluster_key") == cluster_key:
                best_n = r["n_neighbors"]
                break
        if best_n is None:
            best_n = results_summary[0]["n_neighbors"]
            log.warning(
                "DE-gated: resolution=%.2f entry not found, using n_neighbors=%d",
                resolution,
                best_n,
            )
        best_r = resolution
        method_name = "de_gated"
        reason = reason_str
        log.info(
            "DE-gated selection: n_neighbors=%d, resolution=%.2f (%s)", best_n, best_r, reason
        )
    else:
        # ── Multi-metric enrichment (for multi_metric selection method) ──
        import logging as _logging

        from core.cluster.evaluation import (
            _compute_cluster_coherence,
            _compute_splitting_gain,
            _compute_stability,
        )

        _log_enrich = _logging.getLogger(__name__)

        marker_dict = getattr(cfg.marker, "marker_dict", None) or {}
        has_markers = bool(marker_dict)
        n_stab_seeds = getattr(cfg.clustering, "multi_metric_n_stability_seeds", 5)
        dominance_threshold = getattr(cfg.clustering, "multi_metric_coverage_ratio_threshold", 1.5)

        # Group results by n_neighbors
        by_n = {}
        for r in results_summary:
            n = r.get("n_neighbors")
            by_n.setdefault(n, []).append(r)

        for n_val, group in by_n.items():
            # Rebuild KNN graph for this n_neighbors group
            try:
                sc.pp.neighbors(
                    adata,
                    n_neighbors=n_val,
                    n_pcs=cfg.pca.n_pcs_use,
                    use_rep=use_rep,
                    random_state=cfg.execution.random_seed,
                )
            except Exception as e:
                _log_enrich.warning(
                    "KNN rebuild failed for n_neighbors=%d: %s — skipping group", n_val, e
                )
                continue

            # Pre-compute per_cell_scores once per group (only if markers available AND adata.raw exists)
            per_cell_scores = {}
            if has_markers and adata.raw is not None:
                from anndata import utils as anndata_utils

                adata.raw._var.index = anndata_utils.make_index_unique(
                    adata.raw._var.index, join="-"
                )
                try:
                    for ct, genes in marker_dict.items():
                        valid_genes = [g for g in genes if g in adata.raw.var_names]
                        if valid_genes:
                            sc.tl.score_genes(
                                adata, gene_list=valid_genes, score_name=f"_score_{ct}"
                            )
                            per_cell_scores[ct] = adata.obs[f"_score_{ct}"].values.copy()
                    # Clean up temporary score columns
                    for col in list(adata.obs.columns):
                        if col.startswith("_score_") and col in adata.obs.columns:
                            adata.obs.drop(columns=[col], inplace=True)
                except Exception as e:
                    _log_enrich.warning(
                        "Marker score pre-computation failed: %s — falling back to no markers", e
                    )
                    has_markers = False
            elif has_markers and adata.raw is None:
                _log_enrich.warning(
                    "adata.raw is None — cannot compute marker coverage. Degrading to silhouette+stability only."
                )
                has_markers = False

            # Compute stability + marker_coverage for each combo in this group
            for entry in group:
                try:
                    resolution = entry["resolution"]
                    ck = entry["cluster_key"]
                    entry["stability_score"] = _compute_stability(
                        adata,
                        resolution=resolution,
                        leiden_flavor=cfg.clustering.leiden_flavor,
                        n_seeds=n_stab_seeds,
                    )
                    if has_markers and per_cell_scores:
                        entry["cluster_coherence"] = _compute_cluster_coherence(
                            adata,
                            ck,
                            per_cell_scores,
                            dominance_threshold=dominance_threshold,
                        )
                    # ── KB annotatable rate ──
                    if getattr(cfg, "tissue_kb", None) and per_cell_scores:
                        labels = adata.obs[ck].values
                        unique_clusters = np.unique(labels)
                        n_total = len(unique_clusters)
                        n_annotatable = 0
                        for cl in unique_clusters:
                            mask = labels == cl
                            best_score = 0.0
                            for ct in per_cell_scores:
                                scores = per_cell_scores[ct]
                                if scores is not None and len(scores) == len(labels):
                                    mean_val = float(np.mean(scores[mask]))
                                    if mean_val > best_score:
                                        best_score = mean_val
                            if best_score > 0.5:
                                n_annotatable += 1
                        rate = n_annotatable / n_total if n_total > 0 else 0.0
                        entry["kb_annotatable_rate"] = rate
                        _log_enrich.info(f"KB annotatable rate: {rate:.3f}")

                except Exception as e:
                    _log_enrich.warning(
                        "Enrichment failed for n_neighbors=%d, resolution=%.1f: %s",
                        entry.get("n_neighbors"),
                        entry.get("resolution"),
                        e,
                    )
                    entry["stability_score"] = None
                    entry["cluster_coherence"] = None
                    entry["kb_annotatable_rate"] = None

            # ── Compute splitting_gain for this n_neighbors group ──
            if len(group) >= 2:
                group_sorted = sorted(group, key=lambda e: e.get("resolution", 0.0))
                gains = _compute_splitting_gain(group_sorted)
                for entry in group:
                    entry["splitting_gain"] = gains.get(entry["resolution"], 0.0)
    # ── Single-param UMAP plots ──
    for n in n_neighbors_grid:
        for res in resolutions_grid:
            umap_key = f"umap_{n}_{res}"
            leiden_key = f"leiden_{n}_{res}"
            if umap_key not in adata.obsm or leiden_key not in adata.obs:
                continue
            saved = adata.obsm.get("X_umap")
            adata.obsm["X_umap"] = adata.obsm[umap_key].copy()
            try:
                safe_plot(
                    sc.pl.umap,
                    adata,
                    color=leiden_key,
                    show=False,
                    title=f"UMAP (n_neighbors={n}, resolution={res})",
                )
                plt.savefig(
                    os.path.join(fig_dir, f"umap_grid_n{n}_r{res}.png"),
                    dpi=150,
                    bbox_inches="tight",
                )
                plt.close()
                log.info("    Plot saved: umap_grid_n%d_r%.1f.png", n, res)
            except Exception as e:
                log.warning("    Single-param UMAP plot save failed: %s", e)
            finally:
                if saved is not None:
                    adata.obsm["X_umap"] = saved

    # ── 汇总 CSV ──
    df_summary = pd.DataFrame(results_summary)
    csv_path = os.path.join(cfg.table_dir, "param_grid_summary.csv")
    try:
        df_summary.to_csv(csv_path, index=False)
        log.info("Parameter grid summary saved: %s", csv_path)
        log.info("\n%s", df_summary.to_string())
    except Exception as e:
        log.warning("Summary CSV save failed: %s", e)

    if not results_summary:
        log.critical(
            "All neighbor/cluster computations failed — no parameter combination succeeded"
        )
        sys.exit(1)

    # ── 自动选择最佳参数并生成最终 checkpoint ──
    df_summary = pd.DataFrame(results_summary)

    method = getattr(cfg.clustering, "cluster_selection_method", "pareto_elbow")

    # Warn if best_resolution is explicitly set to non-default but will be ignored
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

    if not _de_gated_selected:
        best_n, best_r, method_name, reason = select_best_params(
            results_summary,
            method=method,
            best_resolution=cfg.clustering.best_resolution if method is None else None,
            best_n_neighbors=getattr(cfg.clustering, "best_n_neighbors", 0)
            if method is None
            else 0,
            multi_metric_weights=getattr(cfg.clustering, "multi_metric_weights", None),
            log=log,
        )

    log.info(
        "Selected best params via %s: n_neighbors=%d, resolution=%.1f (%s)",
        method_name,
        best_n,
        best_r,
        reason,
    )

    adata.uns["best_resolution"] = best_r
    adata.uns["best_n_neighbors"] = best_n
    adata.uns["cluster_selection_method"] = method_name

    leiden_col = f"leiden_{best_n}_{best_r}"
    umap_col = f"umap_{best_n}_{best_r}"

    if leiden_col in adata.obs and umap_col in adata.obsm:
        adata.obs["leiden"] = adata.obs[leiden_col].copy()
        adata.obsm["X_umap"] = adata.obsm[umap_col].copy()
        safe_write(adata, cfg.cluster_h5ad, cfg=cfg)
        log.info("Final checkpoint saved: %s (resolution=%.1f)", cfg.cluster_h5ad, best_r)
    else:
        log.warning(
            "Selected param combination (%s, %s) not in results, skipping auto-lock",
            leiden_col,
            umap_col,
        )

    # ── UMAP 参数扫描 (min_dist × spread) ──────────────────────────────────
    min_dist_grid = getattr(cfg.clustering, "param_grid_min_dist", [0.3])
    spread_grid = getattr(cfg.clustering, "param_grid_spread", [1.0])
    umap_method = getattr(cfg.clustering, "umap_selection_method", "convex_hull")
    use_paga = getattr(cfg.clustering, "umap_paga_init", False)

    # If PAGA init is enabled, compute PAGA backbone first
    if use_paga:
        log.info("Computing PAGA backbone for UMAP initialization...")
        sc.tl.paga(adata, groups="leiden")
        sc.pl.paga(adata, show=False)
        plt.savefig(os.path.join(fig_dir, "paga_backbone.png"), dpi=150, bbox_inches="tight")
        plt.close()
        log.info("  PAGA backbone computed and saved")

    best_md, best_sp, umap_method_label, sweep_results = select_best_umap_params(
        adata, best_n, min_dist_grid, spread_grid, umap_method, cfg, use_rep, log
    )

    # Rebuild UMAP with selected params and re-save checkpoint
    log.info(
        "Rebuilding UMAP with selected params (min_dist=%.2f, spread=%.1f) [%s]...",
        best_md,
        best_sp,
        umap_method_label,
    )
    try:
        sc.tl.umap(
            adata,
            min_dist=best_md,
            spread=best_sp,
            init_pos="paga" if use_paga else "spectral",
            random_state=cfg.execution.random_seed,
        )
        safe_write(adata, cfg.cluster_h5ad, cfg=cfg)
        log.info("Checkpoint saved with final UMAP: %s", cfg.cluster_h5ad)
    except Exception as e:
        log.warning("Final UMAP rebuild failed: %s", e)

    if sweep_results:
        # Summary CSV
        try:
            sweep_csv = os.path.join(cfg.table_dir, "umap_min_dist_sweep_summary.csv")
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
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
            axes_flat = axes.ravel() if n_total > 1 else [axes]
            for idx, r in enumerate(sweep_results):
                ax = axes_flat[idx]
                md = r["min_dist"]
                sp = r["spread"]
                try:
                    sc.tl.umap(
                        adata, min_dist=md, spread=sp, random_state=cfg.execution.random_seed
                    )
                    sc.pl.umap(
                        adata,
                        color="leiden",
                        ax=ax,
                        show=False,
                        legend_fontsize=8,
                        title=f"min_dist={md}, spread={sp}",
                    )
                except Exception:
                    ax.text(0.5, 0.5, "Error", ha="center", va="center", transform=ax.transAxes)
            for j in range(len(sweep_results), len(axes_flat)):
                axes_flat[j].axis("off")
            fig.tight_layout()
            fig.savefig(
                os.path.join(fig_dir, "umap_min_dist_comparison.png"), dpi=150, bbox_inches="tight"
            )
            plt.close(fig)
            log.info("UMAP min_dist comparison plot saved")
        except Exception as e:
            log.warning("UMAP comparison plot generation failed: %s", e)

    # ── 网格汇总图: 所有参数组合对比 ──
    n_n = len(n_neighbors_grid)
    n_r = len(resolutions_grid)
    try:
        fig, axes = plt.subplots(n_n, n_r, figsize=(5 * n_r + 2, 4 * n_n + 1), squeeze=False)
        for i, n in enumerate(n_neighbors_grid):
            for j, res in enumerate(resolutions_grid):
                ax = axes[i, j]
                umap_key = f"umap_{n}_{res}"
                leiden_key = f"leiden_{n}_{res}"
                if umap_key in adata.obsm and leiden_key in adata.obs:
                    saved_umap = adata.obsm["X_umap"].copy()
                    try:
                        adata.obsm["X_umap"] = adata.obsm[umap_key].copy()
                        sc.pl.umap(
                            adata,
                            color=leiden_key,
                            ax=ax,
                            legend_fontsize=5,
                            title=f"n={n}, r={res}",
                        )
                    except Exception as e_sub:
                        log.warning("  Subplot failed (n=%d, r=%.1f): %s", n, res, e_sub)
                        ax.text(
                            0.5, 0.5, "Error", ha="center", va="center", transform=ax.transAxes
                        )
                    finally:
                        adata.obsm["X_umap"] = saved_umap
                else:
                    ax.text(
                        0.5,
                        0.5,
                        "N/A",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=12,
                    )
                    ax.set_title(f"n={n}, r={res}")
        fig.tight_layout()
        fig.savefig(
            os.path.join(fig_dir, "umap_param_grid_summary.png"), dpi=150, bbox_inches="tight"
        )
        plt.close(fig)
        log.info("Parameter grid summary plot saved")
    except Exception as e:
        log.warning("Grid summary plot generation failed: %s", e)

    # ── 按 n_neighbors 分组的多分辨率对比图 ──
    for n in n_neighbors_grid:
        res_keys = [f"leiden_{n}_{r}" for r in resolutions_grid if f"leiden_{n}_{r}" in adata.obs]
        n_res = len(res_keys)
        if n_res > 0:
            try:
                n_cols = min(3, n_res)
                n_rows = int(np.ceil(n_res / n_cols))
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
                axes = axes.ravel() if n_res > 1 else [axes]
                for i, key in enumerate(res_keys):
                    sc.pl.umap(
                        adata, color=key, ax=axes[i], show=False, legend_fontsize=8, title=key
                    )
                for j in range(len(res_keys), len(axes)):
                    axes[j].axis("off")
                fig.tight_layout()
                fig.savefig(
                    os.path.join(fig_dir, f"umap_leiden_n{n}_all_resolutions.pdf"),
                    dpi=150,
                    bbox_inches="tight",
                )
                plt.close(fig)
                log.info("  Multi-resolution UMAP plot (n=%d) saved", n)
            except Exception as e:
                log.warning("  Multi-resolution comparison plot (n=%d) failed: %s", n, e)

    # ── 保存临时 h5ad (非最终 checkpoint) ──
    temp_path = os.path.join(cfg.h5ad_dir, "04_grid_results.h5ad")
    try:
        safe_write(adata, temp_path, cfg=cfg)
        log.info("Temporary h5ad saved: %s", temp_path)
    except Exception as e:
        log.error("Temporary h5ad save failed: %s", e)

    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
