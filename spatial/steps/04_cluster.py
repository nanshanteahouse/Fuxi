#!/usr/bin/env python3
"""
Step 04: Neighbors + UMAP + multi-param grid Leiden clustering
================================================
  - Build PCA neighbor graph
  - Multi-param grid scan (n_neighbors × resolution)
  - Multi-metric grid search: silhouette + stability + cluster_coherence + splitting_gain + kb_annotatable_rate
  - Adaptive resolution expansion for small datasets
  - Granularity detection (tissue vs subtype)
  - DE-gated resolution selection for subtype granularity
  - UMAP min_dist × spread sweep with select_best_umap_params
  - Auto-select best params via composite scoring
  - Generate UMAP visualizations

Input:  03_processed.h5ad
Output: 04_clustered.h5ad
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
    log = setup_logger("04_cluster", os.path.join(cfg.log_dir, "04_cluster.log"))
    log.info("Step 04: Neighbors + UMAP + multi-param grid Leiden clustering")

    # ── Checkpoint ──
    output_path = os.path.join(cfg.h5ad_dir, "04_clustered.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    # ── Input ──
    input_path = os.path.join(cfg.h5ad_dir, "03_processed.h5ad")
    if not os.path.exists(input_path):
        log.error("Input not found: %s", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    use_rep = "X_pca"
    if "X_pca" not in adata.obsm:
        log.error("No PCA found in obsm. Run Step 03 first.")
        sys.exit(1)

    log.info("Using PCA representation: %s (%d PCs)", use_rep, cfg.pca.n_pcs_use)

    # ── Parameter grid ──
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
                    "Adaptive resolution expansion (n_cells=%d < 3000): added %s \u2192 %s",
                    adata.n_obs,
                    extra,
                    resolutions_grid,
                )

    log.info("Grid: n_neighbors=%s, resolutions=%s", n_neighbors_grid, resolutions_grid)

    fig_dir = os.path.join(cfg.figure_dir, "04_cluster")
    os.makedirs(fig_dir, exist_ok=True)

    umap_min_dist = getattr(cfg.clustering, "umap_min_dist", 0.3)
    umap_spread = getattr(cfg.clustering, "umap_spread", 1.0)

    # ── Scanpy-specific callables for grid search ──
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
            flavor=getattr(cfg.clustering, "leiden_flavor", "igraph"),
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
        random_seed=cfg.execution.random_seed,
    )

    # Rename score \u2192 silhouette_score for select_best_params compatibility
    for r in results_summary:
        if "score" in r:
            r["silhouette_score"] = r.pop("score")

    # ═══ Multi-metric enrichment ═══
    from core.cluster.evaluation import _detect_granularity

    granularity = _detect_granularity(results_summary)
    log.info("Granularity classification: %s", granularity)

    _de_gated_selected = False

    if granularity == "subtype" and adata.raw is not None:
        _de_gated_selected = True
        from core.cluster.evaluation import _select_de_gated

        log.info(
            "Granularity=subtype \u2014 using DE-gated resolution selection (bypassing enrichment)"
        )
        n_clusters, resolution, cluster_key, reason_str = _select_de_gated(
            results_summary,
            adata,
            de_gate_threshold=getattr(cfg.clustering, "multi_metric_de_gate_threshold", 25),
        )
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
        if granularity == "subtype":
            log.warning(
                "adata.raw is None \u2014 cannot run DE-gated selection. "
                "Falling back to multi-metric enrichment."
            )

        # ── Multi-metric enrichment (for multi_metric selection method) ──
        from core.cluster.evaluation import enrich_grid_results

        enrich_grid_results(
            adata,
            results_summary,
            cfg,
            log=log,
            use_rep=use_rep,
        )

    # ── Generate per-param UMAP plots ──
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
                    cfg=cfg,
                )
                plt.savefig(
                    os.path.join(fig_dir, f"umap_grid_n{n}_r{res}.png"),
                    dpi=cfg.plot.figure_dpi,
                    bbox_inches="tight",
                )
                plt.close()
            except Exception as e:
                log.warning("    Plot save failed: %s", e)
            finally:
                if saved is not None:
                    adata.obsm["X_umap"] = saved

    # ── Summary CSV ──
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
            "All clustering computations failed \u2014 no parameter combination succeeded"
        )
        sys.exit(1)

    # ── Auto-select best params and generate final checkpoint ──
    df_summary = pd.DataFrame(results_summary)

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
        safe_write(adata, output_path, cfg=cfg)
        log.info("Final checkpoint saved: %s (resolution=%.1f)", output_path, best_r)
    else:
        log.warning(
            "Selected param combination (%s, %s) not in results, skipping auto-lock",
            leiden_col,
            umap_col,
        )

    # ── UMAP parameter sweep (min_dist \u00d7 spread) ──
    min_dist_grid = getattr(cfg.clustering, "param_grid_min_dist", [0.3])
    spread_grid = getattr(cfg.clustering, "param_grid_spread", [1.0])
    umap_method = getattr(cfg.clustering, "umap_selection_method", "convex_hull")
    use_paga = getattr(cfg.clustering, "umap_paga_init", False)

    if use_paga:
        log.warning(
            "umap_paga_init=True \u2014 PAGA is designed for single-cell resolution, not spots. "
            "Proceeding with PAGA init anyway per config."
        )

    best_md, best_sp, umap_method_label, sweep_results = select_best_umap_params(
        adata,
        best_n,
        min_dist_grid,
        spread_grid,
        umap_method,
        cfg,
        use_rep,
        log,
        metric=getattr(cfg.clustering, "umap_selection_metric", "trustworthiness"),
    )
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
        safe_write(adata, output_path, cfg=cfg)
        log.info("Checkpoint saved with final UMAP: %s", output_path)
    except Exception as e:
        log.warning("Final UMAP rebuild failed: %s", e)

    if sweep_results:
        try:
            sweep_csv = os.path.join(cfg.table_dir, "umap_min_dist_sweep_summary.csv")
            pd.DataFrame(sweep_results).to_csv(sweep_csv, index=False)
            log.info("UMAP sweep summary saved: %s", sweep_csv)
        except Exception as e:
            log.warning("UMAP sweep CSV save failed: %s", e)

        try:
            n_md = len(min_dist_grid) if min_dist_grid else 1
            n_sp = len(spread_grid) if spread_grid else 1
            n_total = n_md * n_sp
            n_cols = min(3, n_total)
            n_rows = int(np.ceil(n_total / n_cols))
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=(
                    cfg.plot.umap_panel_size[0] * n_cols,
                    cfg.plot.umap_panel_size[1] * n_rows,
                ),
            )
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
                except Exception as e:
                    log.debug("Plot annotation failed: %s", e)
                    ax.text(0.5, 0.5, "Error", ha="center", va="center", transform=ax.transAxes)
            for j in range(len(sweep_results), len(axes_flat)):
                axes_flat[j].axis("off")
            fig.tight_layout()
            fig.savefig(
                os.path.join(fig_dir, "umap_min_dist_comparison.png"),
                dpi=cfg.plot.figure_dpi,
                bbox_inches="tight",
            )
            plt.close(fig)
            log.info("UMAP min_dist comparison plot saved")
        except Exception as e:
            log.warning("UMAP comparison plot generation failed: %s", e)

    # ── Grid summary plot: all parameter combinations ──
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
                    except Exception as e:
                        log.debug("Plot annotation failed: %s", e)
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
            os.path.join(fig_dir, "umap_param_grid_summary.png"),
            dpi=cfg.plot.figure_dpi,
            bbox_inches="tight",
        )
        plt.close(fig)
        log.info("Parameter grid summary plot saved")
    except Exception as e:
        log.warning("Grid summary plot generation failed: %s", e)

    # ── Multi-resolution comparison plots per n_neighbors ──
    for n in n_neighbors_grid:
        res_keys = [f"leiden_{n}_{r}" for r in resolutions_grid if f"leiden_{n}_{r}" in adata.obs]
        n_res = len(res_keys)
        if n_res > 0:
            try:
                n_cols = min(3, n_res)
                n_rows = int(np.ceil(n_res / n_cols))
                fig, axes = plt.subplots(
                    n_rows,
                    n_cols,
                    figsize=(
                        cfg.plot.umap_panel_size[0] * n_cols,
                        cfg.plot.umap_panel_size[1] * n_rows,
                    ),
                )
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
                    dpi=cfg.plot.figure_dpi,
                    bbox_inches="tight",
                )
                plt.close(fig)
                log.info("  Multi-resolution UMAP plot (n=%d) saved", n)
            except Exception as e:
                log.warning("  Multi-resolution comparison plot (n=%d) failed: %s", n, e)

    # ── Save temporary h5ad (non-final checkpoint) ──
    temp_path = os.path.join(cfg.h5ad_dir, "04_grid_results.h5ad")
    try:
        safe_write(adata, temp_path, cfg=cfg)
        log.info("Temporary h5ad saved: %s", temp_path)
    except Exception as e:
        log.error("Temporary h5ad save failed: %s", e)

    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
