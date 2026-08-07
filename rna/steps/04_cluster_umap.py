#!/usr/bin/env python3
"""
Step 04: 邻居图 + UMAP + 多参数网格 Leiden 聚类
==================================================
  - 在 Harmony 校正后的 PCA 上建图
  - 多参数网格扫描 (n_neighbors × resolution)
  - 保存所有组合结果用于交互比较

输入: 03_integrated.h5ad
输出: 04_clustered.h5ad (所有参数组合的邻居图、UMAP、Leiden 标签)
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
from core.utils import (
    gpu_leiden,
    gpu_neighbors,
    gpu_umap,
    resolve_config,
    safe_write,
    save_figure,
    setup_logger,
    timed_substep,
)

# Plot-capped stratum sizes for UMAP scatter subsampling.
# Rare clusters (<= cap_small cells) are kept in full; medium clusters are
# capped at cap_medium; very large clusters at cap_large — the visual keeps
# rare populations visible instead of diluting them under big clusters.
# Configurable via cfg.clustering.umap_plot_cap_{small,medium,large} and
# umap_plot_stratum_large.
_PLOT_CAP_SMALL = 500
_PLOT_CAP_MEDIUM = 500
_PLOT_CAP_LARGE = 1000
_PLOT_STRATUM_LARGE = 50_000


def _plot_caps(cfg):
    c = cfg.clustering
    return (
        getattr(c, "umap_plot_cap_small", _PLOT_CAP_SMALL),
        getattr(c, "umap_plot_cap_medium", _PLOT_CAP_MEDIUM),
        getattr(c, "umap_plot_cap_large", _PLOT_CAP_LARGE),
        getattr(c, "umap_plot_stratum_large", _PLOT_STRATUM_LARGE),
    )


def _plot_subsample_idx(labels, max_cells, rng, caps=None):
    """Stratified subsample index for UMAP scatter plots.

    Tiers: label size <= cap_small kept in full; cap_small-50k capped at
    cap_medium; > stratum_large capped at cap_large. Then truncated to
    ``max_cells`` as a hard ceiling. ``caps`` is the tuple from
    ``_plot_caps(cfg)``; defaults to module constants when None.
    """
    cap_small, cap_medium, cap_large, stratum_large = caps or (
        _PLOT_CAP_SMALL,
        _PLOT_CAP_MEDIUM,
        _PLOT_CAP_LARGE,
        _PLOT_STRATUM_LARGE,
    )
    codes = pd.Categorical(labels).codes
    parts = []
    for code in np.unique(codes):
        sel = np.where(codes == code)[0]
        n = len(sel)
        if n <= cap_small:
            parts.append(sel)
        else:
            cap = cap_medium if n <= stratum_large else cap_large
            parts.append(rng.choice(sel, min(cap, n), replace=False))
    idx = np.concatenate(parts)
    if len(idx) > max_cells:
        idx = rng.choice(idx, max_cells, replace=False)
    return idx


def plot_step04_figures(h5ad_path: str, cfg, log):
    """Re-render Step 04 summary figures from a saved checkpoint (plot-only).

    Deterministic renderer of ``04_clustered.h5ad`` -> ``figures/04_cluster/``:
    parameter-grid summary, per-n_neighbors multi-resolution panels, and batch
    diagnostics. Uses a light AnnData (``X_umap`` + obs columns only) so memory
    stays ~100MB even for 1M+ cell datasets. Safe to re-run any time (after a
    plotting failure, or a dpi/colour/subsample change) without recomputing
    clustering.
    """
    import anndata
    from scipy import sparse as _sp

    fig_dir = os.path.join(cfg.figure_dir, "04_cluster")
    os.makedirs(fig_dir, exist_ok=True)
    backed = sc.read_h5ad(h5ad_path, backed="r")
    try:
        light = anndata.AnnData(
            X=_sp.csr_matrix((backed.n_obs, 1), dtype="float32"),
            obs=backed.obs,
            obsm={
                k: backed.obsm[k]
                for k in backed.obsm.keys()
                if k.startswith("X_umap") or k.startswith("umap_")
            },
        )
    finally:
        backed.file.close()

    # Recover the (n_neighbors, resolution) grid from obs leiden_* columns.
    pairs: list[tuple[int, float]] = []
    for c in light.obs.columns:
        parts = c.split("_")
        if len(parts) == 3 and parts[0] == "leiden":
            try:
                pairs.append((int(parts[1]), float(parts[2])))
            except ValueError:
                pass
    n_grid = sorted({n for n, _ in pairs})
    r_grid = sorted({r for _, r in pairs})
    if not n_grid:
        log.warning("plot_step04_figures: no leiden_* columns in checkpoint, nothing to draw")
        return

    # 1) Parameter grid summary: n_neighbors x resolution UMAP panels
    try:
        n_cols = len(r_grid)
        n_rows = len(n_grid)
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(5 * n_cols + 2, 4 * n_rows + 1), squeeze=False
        )
        for i, n in enumerate(n_grid):
            for j, r in enumerate(r_grid):
                ax = axes[i, j]
                umap_key = f"umap_{n}_{r}"
                leiden_key = f"leiden_{n}_{r}"
                if umap_key in light.obsm and leiden_key in light.obs:
                    saved = light.obsm["X_umap"]
                    light.obsm["X_umap"] = light.obsm[umap_key]
                    try:
                        _smart_plot_umap(
                            light,
                            color=leiden_key,
                            ax=ax,
                            title=f"n={n}, r={r}",
                            cfg=cfg,
                            log=log,
                            legend_fontsize=5,
                        )
                    finally:
                        light.obsm["X_umap"] = saved
                else:
                    ax.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax.transAxes)
                    ax.set_title(f"n={n}, r={r}")
        fig.tight_layout()
        save_figure(
            fig,
            os.path.join(fig_dir, "param_grid_summary"),
            cfg=cfg,
            dpi=cfg.plot.figure_dpi,
            bbox_inches="tight",
        )
        plt.close(fig)
        log.info("plot-only: Parameter grid summary saved")
    except Exception as e:
        log.warning("plot-only: Grid summary plot failed: %s", e)

    # 2) Multi-resolution panels per n_neighbors (final X_umap)
    for n in n_grid:
        keys = [f"leiden_{n}_{r}" for r in r_grid if f"leiden_{n}_{r}" in light.obs]
        if not keys:
            continue
        try:
            n_cols = min(3, len(keys))
            n_rows = int(np.ceil(len(keys) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
            axes = axes.ravel() if len(keys) > 1 else [axes]
            for i, key in enumerate(keys):
                _smart_plot_umap(
                    light, color=key, ax=axes[i], title=key, cfg=cfg, log=log, legend_fontsize=8
                )
            for j in range(len(keys), len(axes)):
                axes[j].axis("off")
            fig.tight_layout()
            save_figure(
                fig,
                os.path.join(fig_dir, f"leiden_multires_n{n}"),
                cfg=cfg,
                dpi=cfg.plot.figure_dpi,
                bbox_inches="tight",
            )
            plt.close(fig)
            log.info("plot-only: Multi-resolution Leiden plot (n=%d) saved", n)
        except Exception as e:
            log.warning("plot-only: Multi-resolution plot (n=%d) failed: %s", n, e)

    # 3) Batch diagnostics
    if getattr(cfg.clustering, "umap_color_by_batch", False):
        batch_key = (
            getattr(cfg.clustering, "batch_key_override", None) or cfg.integration.batch_key
        )
        if batch_key in light.obs:
            batches = light.obs[batch_key].unique()
            n_batches = len(batches)
            if n_batches < 2:
                log.info("plot-only: n_batches=%d < 2, skipping batch UMAP", n_batches)
            else:
                try:
                    fig_b, ax_b = plt.subplots(figsize=(6, 5))
                    if _smart_plot_umap(
                        light, batch_key, ax_b, f"UMAP colored by {batch_key}", cfg, log
                    ):
                        save_figure(
                            fig_b,
                            os.path.join(fig_dir, "batch_colored_umap"),
                            cfg=cfg,
                            dpi=cfg.plot.figure_dpi,
                            bbox_inches="tight",
                        )
                    plt.close(fig_b)
                    log.info("plot-only: Batch-colored UMAP saved")
                except Exception as e:
                    log.warning("plot-only: Batch-colored UMAP failed: %s", e)
                if n_batches <= 12:
                    try:
                        n_cols = min(4, n_batches)
                        n_rows = int(np.ceil(n_batches / n_cols))
                        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
                        axes = axes.ravel() if n_batches > 1 else [axes]
                        rng = np.random.RandomState(getattr(cfg.execution, "random_seed", 42))
                        bg_cap = min(50000, light.n_obs)
                        bg_idx = rng.choice(light.n_obs, bg_cap, replace=False)
                        for i, bv in enumerate(batches):
                            ax = axes[i]
                            mask = light.obs[batch_key] == bv
                            ax.scatter(
                                light.obsm["X_umap"][bg_idx, 0],
                                light.obsm["X_umap"][bg_idx, 1],
                                c="lightgray",
                                s=1,
                                alpha=0.3,
                                rasterized=True,
                            )
                            fg_all = np.where(mask)[0]
                            fg_sub = fg_all[
                                _plot_subsample_idx(
                                    light.obs.loc[mask, "leiden"].values,
                                    50000,
                                    rng,
                                    _plot_caps(cfg),
                                )
                            ]
                            ax.scatter(
                                light.obsm["X_umap"][fg_sub, 0],
                                light.obsm["X_umap"][fg_sub, 1],
                                c=light.obs.loc[fg_sub, "leiden"].astype("category").cat.codes,
                                cmap=cfg.plot.palette.categorical,
                                s=3,
                                alpha=0.8,
                                rasterized=True,
                            )
                            ax.set_title(f"{batch_key}={bv} ({mask.sum()} cells)", fontsize=9)
                            ax.set_xticks([])
                            ax.set_yticks([])
                        for j in range(n_batches, len(axes)):
                            axes[j].axis("off")
                        fig.tight_layout()
                        save_figure(
                            fig,
                            os.path.join(fig_dir, "batch_faceted_umap"),
                            cfg=cfg,
                            dpi=cfg.plot.figure_dpi,
                            bbox_inches="tight",
                        )
                        plt.close(fig)
                        log.info("plot-only: Split-by-batch faceted UMAP saved")
                    except Exception as e:
                        log.warning("plot-only: Faceted UMAP failed: %s", e)

    log.info("plot_step04_figures: done")


def _smart_plot_umap(adata, color, ax, title, cfg, log, legend_fontsize=8):
    """Render a UMAP scatter respecting cfg.clustering.umap_plot_mode.

    Modes (auto resolves at runtime based on adata.n_obs vs umap_plot_max_cells):
      * 'full'      — call sc.pl.umap directly (best visual fidelity, slow on >100k cells)
      * 'subsample' — random subsample to umap_plot_max_cells, then plain ax.scatter
      * 'skip'      — do nothing, return False
      * 'auto'      — 'full' when n_obs <= max_cells, else 'subsample'

    Returns True if something was drawn, False if skipped. Caller must still
    save the figure (plt.savefig) on success — this helper only draws onto *ax*.
    """

    mode = getattr(cfg.clustering, "umap_plot_mode", "auto")
    max_cells = getattr(cfg.clustering, "umap_plot_max_cells", 50000)
    if mode == "auto":
        mode = "subsample" if adata.n_obs > max_cells else "full"

    if mode == "skip":
        ax.text(
            0.5,
            0.5,
            "skipped\n(umap_plot_mode=skip)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10,
        )
        return False

    if mode == "subsample" and adata.n_obs > max_cells:
        rng = np.random.RandomState(getattr(cfg.execution, "random_seed", 42))
        idx = _plot_subsample_idx(adata.obs[color], max_cells, rng, _plot_caps(cfg))
        coords = adata.obsm["X_umap"][idx]
        # Map categorical labels to integer codes for scatter coloring.
        cat = pd.Categorical(adata.obs[color])
        codes = cat.codes[idx]
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=codes,
            cmap=getattr(cfg.plot.palette, "categorical", "tab20"),
            s=3,
            alpha=0.6,
            rasterized=True,
        )
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        log.info("    [plot] stratified subsample: %d/%d cells drawn", len(idx), adata.n_obs)
        return True

    # 'full' (or 'subsample' with n_obs <= max_cells)
    sc.pl.umap(
        adata,
        color=color,
        ax=ax,
        show=False,
        title=title,
        legend_fontsize=legend_fontsize,
    )
    return True


def _plot_umap_from_coords(coords, adata, color, ax, title, cfg, log, legend_fontsize=8):
    """Plot UMAP scatter directly from saved coords (no UMAP recomputation).

    Used by sweep comparison figure to avoid the O(47min × N combos) cost of
    re-running ``sc.tl.umap`` for cosmetic plotting. Reuses the per-combo
    ``coords`` saved in ``sweep_results`` by ``select_best_umap_params``.

    Honors ``cfg.clustering.umap_plot_mode`` (auto/full/subsample/skip) the
    same way ``_smart_plot_umap`` does — ensures consistent visual treatment
    between single-param plots and the sweep comparison figure.

    Args:
        coords: ``(n_obs, 2)`` ndarray of UMAP coordinates for this combo.
        adata: AnnData for obs lookups (color labels).
        color: column in ``adata.obs`` to color by (typically ``'leiden'``).
        ax: matplotlib axes to draw on.
        title: subplot title.
        cfg, log, legend_fontsize: same as ``_smart_plot_umap``.

    Returns:
        True if a scatter was drawn; False if ``skip`` mode suppressed it.
    """
    mode = getattr(cfg.clustering, "umap_plot_mode", "auto")
    max_cells = getattr(cfg.clustering, "umap_plot_max_cells", 50000)
    if mode == "auto":
        mode = "subsample" if adata.n_obs > max_cells else "full"
    if mode == "skip":
        ax.text(0.5, 0.5, f"{title}\n(skipped)", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return False
    if mode == "subsample" and adata.n_obs > max_cells:
        rng = np.random.RandomState(getattr(cfg.execution, "random_seed", 42))
        idx = _plot_subsample_idx(adata.obs[color], max_cells, rng, _plot_caps(cfg))
        sub_coords = coords[idx]
        cat = pd.Categorical(adata.obs[color])
        codes = cat.codes[idx]
        ax.scatter(
            sub_coords[:, 0],
            sub_coords[:, 1],
            c=codes,
            cmap=getattr(cfg.plot.palette, "categorical", "tab20"),
            s=3,
            alpha=0.6,
            rasterized=True,
        )
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        log.info("[plot] %s: stratified subsample %d/%d cells", title, len(idx), adata.n_obs)
        return True
    # full mode
    cat = pd.Categorical(adata.obs[color])
    codes = cat.codes
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=codes,
        cmap=getattr(cfg.plot.palette, "categorical", "tab20"),
        s=3,
        alpha=0.6,
        rasterized=True,
    )
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    return True


# ── Incremental write (copy+append, mode A) ─────────────────────────────────
# Step 04 owns exactly these h5ad keys; everything else (X, layers, raw,
# X_pca / X_integrated, 03-era uns like pca variance_ratio / integration
# results, 03-era obs columns) is preserved by the copy of 03 and never
# rewritten. Grid-derived keys (leiden_{n}_{r}, funnel_{n}_{r},
# umap_{n}_{r}) are matched by prefix — the grid is config-driven.
OWNED_OBSM = {"X_umap"}
OWNED_OBSP = {"connectivities", "distances"}
OWNED_UNS = {
    "neighbors",
    "umap",
    "grid_scan_mode",
    "best_resolution",
    "best_n_neighbors",
    "cluster_selection_method",
    "funnel_lineage",
    "paga",
    "leiden_sizes",
}


def _owned_obsm(adata):
    return {k: adata.obsm[k] for k in adata.obsm if k in OWNED_OBSM or k.startswith("umap_")}


def _owned_obsp(adata):
    return {k: adata.obsp[k] for k in OWNED_OBSP if k in adata.obsp}


def _owned_uns(adata):
    return {
        k: adata.uns[k]
        for k in adata.uns
        if k in OWNED_UNS or k.startswith("leiden_") or k.startswith("funnel_")
    }


def _owned_obs(adata):
    cols = [
        c
        for c in adata.obs.columns
        if c == "leiden" or c.startswith("leiden_") or c.startswith("funnel_")
    ]
    return adata.obs[cols] if cols else None


def _write_cluster_h5ad(adata, cfg, log):
    """Write 04_clustered.h5ad — copy+append (incremental_io) or full safe_write.

    Mode A: the target starts as a copy of 03_integrated.h5ad (X, layers, raw,
    03-era obs/uns all preserved), then only the keys step 04 owns are
    appended/overwritten. A failed append deletes the corrupt copy — 03 stays
    pristine, so re-running the step recovers. ``incremental_io=False`` (old
    configs default True via getattr) falls back to the full safe_write path.
    """
    if getattr(cfg, "incremental_io", True):
        import shutil

        shutil.copy2(cfg.integrated_h5ad, cfg.cluster_h5ad)
        try:
            from core.utils import write_h5ad_incremental

            write_h5ad_incremental(
                cfg.cluster_h5ad,
                obsm=_owned_obsm(adata),
                obsp=_owned_obsp(adata),
                uns=_owned_uns(adata),
                obs=_owned_obs(adata),
                logger=log,
            )
        except Exception:
            if os.path.exists(cfg.cluster_h5ad):
                os.remove(cfg.cluster_h5ad)
            raise
    else:
        safe_write(adata, cfg.cluster_h5ad, cfg=cfg)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args_parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Re-render Step 04 summary figures from the saved 04_clustered.h5ad checkpoint (no clustering).",
    )
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("04_cluster", os.path.join(cfg.log_dir, "04_cluster_umap.log"))

    if args.plot_only:
        log.info("Step 04: plot-only mode (re-render figures from checkpoint)")
        plot_step04_figures(cfg.cluster_h5ad, cfg, log)
        return
    log.info("Step 04: Neighbors + UMAP + multi-param grid Leiden clustering")

    # ── Memory guard: estimate step 04 peak vs budget before the full load ──
    # n_cells prefers the step-00 load_meta (persisted by the runner into
    # results/perf_report.json — filtering only shrinks, so it is a
    # conservative upper bound); falls back to a zero-copy h5py shape probe.
    from core.utils import check_memory_guard, estimate_step_peak, resolve_memory_settings

    input_path = cfg.integrated_h5ad
    _mem_policy, _mem_budget, _mem_guard = resolve_memory_settings(cfg)
    _n_cells = 0
    try:
        import json as _json
        import os as _os

        _pr_path = _os.path.join(cfg.results_dir, "perf_report.json")
        if _os.path.isfile(_pr_path):
            with open(_pr_path) as _f:
                _lm = (_json.load(_f).get("pipeline", {}) or {}).get("load_meta")
            if _lm and _lm.get("n_cells"):
                _n_cells = int(_lm["n_cells"])
    except Exception:
        pass
    if _n_cells <= 0:
        try:
            import h5py

            with h5py.File(input_path, "r") as _h5:
                _n_cells = int(_h5["X"].shape[0])
        except Exception:
            pass
    if _n_cells > 0:
        _est = {
            4: estimate_step_peak(4, _n_cells, 4000, policy=_mem_policy, budget_bytes=_mem_budget)
        }
        if _mem_budget > 0:
            log.info("[memory-guard] estimated step 04 peak: ~%.0f GB", _est[4])
        check_memory_guard(_est, _mem_budget, _mem_guard, logger_obj=log)

    # ── 输入 ──
    input_path = cfg.integrated_h5ad
    log.info("Loaded: %s", input_path)
    adata = sc.read(input_path)
    log.info("  shape: %s", adata.shape)

    use_rep = "X_integrated" if "X_integrated" in adata.obsm else "X_pca"
    log.info("use_rep: %s", use_rep)

    # ── 参数网格 ──
    n_neighbors_grid = getattr(cfg.clustering, "param_grid_n_neighbors", [15, 20, 30])
    # Adaptive n_neighbors by cohort size (param_grid_n_neighbors_adaptive=True):
    # n_neighbors is nearly flat for quality (measured silhouette range 7x smaller than
    # resolution) — a single scale-appropriate value replaces the 3-way grid sweep.
    # Calibrated to Multiome Academy guidance + our 15-dataset grid analysis.
    if getattr(cfg.clustering, "param_grid_n_neighbors_adaptive", False):
        _n = adata.n_obs
        if _n < 5000:
            n_neighbors_grid = [10]
        elif _n < 50000:
            n_neighbors_grid = [15]
        elif _n < 500000:
            n_neighbors_grid = [20]
        else:
            n_neighbors_grid = [30]
        log.info("Adaptive n_neighbors (n_obs=%d): %s", _n, n_neighbors_grid)
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
    table_dir = os.path.join(cfg.table_dir, "04_cluster")
    os.makedirs(table_dir, exist_ok=True)

    umap_min_dist = getattr(cfg.clustering, "umap_min_dist", 0.3)
    umap_spread = getattr(cfg.clustering, "umap_spread", 1.0)

    # ── KB suggestion for target_n_clusters ──
    if getattr(cfg.clustering, "target_n_clusters", None) is None and getattr(cfg, "tissue", None):
        from core.kb import suggest_target_n_clusters

        try:
            suggested = suggest_target_n_clusters(cfg.tissue, getattr(cfg, "species", None))
            if suggested is not None:
                log.info(
                    "KB suggests target_n_clusters=%d for tissue=%s species=%s. "
                    "Set clustering.target_n_clusters=%d in config to activate target mode.",
                    suggested,
                    cfg.tissue,
                    getattr(cfg, "species", None),
                    suggested,
                )
        except Exception:
            log.debug(
                "KB suggestion failed (tissue=%s, species=%s)",
                cfg.tissue,
                getattr(cfg, "species", None),
            )
    # ── Grid search via shared core ──

    # Scanpy-specific callables (close over CFG and use_rep)
    def _neighbors_fn(adata, n_neighbors=15, **kwargs):
        actual_dims = adata.obsm[use_rep].shape[1]
        gpu_neighbors(
            adata,
            log=log,
            device=cfg.execution.device,
            n_neighbors=n_neighbors,
            n_pcs=min(cfg.pca.n_pcs_use, actual_dims),
            use_rep=use_rep,
            random_state=cfg.execution.random_seed,
        )

    def _umap_fn(adata, **kwargs):
        gpu_umap(
            adata,
            log=log,
            device=cfg.execution.device,
            min_dist=umap_min_dist,
            spread=umap_spread,
            maxiter=cfg.clustering.umap_maxiter,
            n_epochs=cfg.clustering.umap_n_epochs,
            random_state=cfg.execution.random_seed,
        )

    def _clusterer_fn(adata, resolution=1.0, n_neighbors=15, **kwargs):
        leiden_key = f"leiden_{n_neighbors}_{resolution}"
        umap_key = f"umap_{n_neighbors}_{resolution}"
        gpu_leiden(
            adata,
            log=log,
            device=(
                _leiden_device := cfg.execution.device
                if cfg.execution.device == "cpu"
                or adata.n_obs >= getattr(cfg.clustering, "leiden_gpu_min_cells", 20_000)
                else "cpu"
            ),
            resolution=resolution,
            key_added=leiden_key,
            random_state=cfg.execution.random_seed,
            flavor=cfg.clustering.leiden_flavor,
            directed=False,
            n_iterations=cfg.clustering.leiden_n_iterations,
        )
        if "X_umap" in adata.obsm:
            adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        return leiden_key

    def _evaluation_fn(adata, cluster_key, **kwargs):
        labels = adata.obs[cluster_key].values
        _rep_full = adata.obsm[use_rep]
        if adata.n_obs > SILHOUETTE_SAMPLE_THRESHOLD:
            rng = np.random.RandomState(cfg.execution.random_seed)
            idx = rng.choice(adata.n_obs, SILHOUETTE_SAMPLE_THRESHOLD, replace=False)
            _rep = _rep_full[idx, : min(cfg.pca.n_pcs_use, _rep_full.shape[1])]
            _labels = labels[idx]
        else:
            _rep = _rep_full[:, : min(cfg.pca.n_pcs_use, _rep_full.shape[1])]
            _labels = labels
        if hasattr(_rep, "get") and not isinstance(_rep, np.ndarray):
            _rep = _rep.get()
        _rep = np.asarray(_rep, dtype=np.float64)
        return float(silhouette_score(_rep, _labels))

    with timed_substep("Grid search (Leiden)", log=log):
        # ── Three-mode dispatch ──
        if getattr(cfg.clustering, "target_n_clusters", None) is not None:
            log.warning("TARGET mode: target_n_clusters=%d", cfg.clustering.target_n_clusters)
            adata.uns["grid_scan_mode"] = "target"
            from core.cluster.target import target_grid_search

            results_summary = target_grid_search(adata, cfg, log=log)

        elif getattr(cfg.clustering, "funnel_enabled", False) and adata.n_obs > getattr(
            cfg.clustering, "funnel_threshold", 100000
        ):
            log.warning(
                "FUNNEL mode: n_obs=%d > threshold=%d. Subsample size=%d, top_k=%d",
                adata.n_obs,
                getattr(cfg.clustering, "funnel_threshold", 100000),
                getattr(cfg.clustering, "funnel_subsample_size", 50000),
                getattr(cfg.clustering, "funnel_top_k", 3),
            )
            adata.uns["grid_scan_mode"] = "funnel"

            def _funnel_full_grid_fn(sub_adata, cfg):
                from core.cluster.evaluation import enrich_grid_results

                sub_results = grid_search_clustering(
                    sub_adata,
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
                    log=log,
                )
                for r in sub_results:
                    if "score" in r:
                        r["silhouette_score"] = r.pop("score")
                enrich_grid_results(sub_adata, sub_results, cfg, log=log, use_rep=use_rep)
                return sub_results

            from core.cluster.funnel import run_funnel_grid_search

            best_entry = run_funnel_grid_search(adata, cfg, _funnel_full_grid_fn, log=log)
            results_summary = [best_entry]

        else:
            log.info("FULL GRID mode: n_obs=%d", adata.n_obs)
            adata.uns["grid_scan_mode"] = "full"
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
                log=log,
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
            pairwise_max_clusters=getattr(cfg.clustering, "de_pairwise_max_clusters", 30),
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
        # ── Sync GPU→CPU to prevent VRAM/CPU double memory pressure ──
        # Grid search left adata.X, raw, obsm, obsp on GPU via
        # rsc.get.anndata_to_GPU(). Enrichment builds CPU neighbor graphs
        # and runs marker scoring on raw.X. Without this sync, VRAM holds
        # grid search artifacts while CPU allocates KNN graphs — easily 10+ GB.
        try:
            import gc

            import cupy as cp

            # obsm: dense arrays (PCA, UMAP coords) → just .get()
            for key in list(adata.obsm.keys()):
                arr = adata.obsm[key]
                if hasattr(arr, "get"):
                    adata.obsm[key] = arr.get()
            # obsp: sparse neighbor graphs → csr_matrix
            for key in list(adata.obsp.keys()):
                mat = adata.obsp[key]
                if hasattr(mat, "get"):
                    from scipy.sparse import csr_matrix

                    adata.obsp[key] = csr_matrix(mat.get())
            # adata.X: HVG-filtered dense data → just .get(), NOT csr_matrix
            # (dense→CSR doubles memory; HVG data is effectively dense)
            if hasattr(adata.X, "get"):
                adata.X = adata.X.get()
            # adata.raw.X: backed sparse dataset, typically not on GPU
            if adata.raw is not None and hasattr(adata.raw.X, "get"):
                adata.raw.X = adata.raw.X.get()
            gc.collect()
            cp.get_default_memory_pool().free_all_blocks()
            log.info("GPU arrays synced to CPU — VRAM released before enrichment")
        except Exception:
            pass  # cupy not available or sync non-critical

        # ── Multi-metric enrichment (for multi_metric selection method) ──
        from core.cluster.evaluation import enrich_grid_results

        enrich_grid_results(
            adata,
            results_summary,
            cfg,
            log=log,
            use_rep=use_rep,
            stability_top_k=getattr(cfg.clustering, "multi_metric_stability_top_k", None),
        )

    # This loop is the largest matplotlib cost in step 04 on big datasets
    # (Li2026: 18 plots × 1M cells ≈ 2h). Config clustering.plot_per_combo
    # controls whether per-combo plots are generated. Disable to save ~80%
    # plot wall time at 1M scale. Summary/aggregate plots still generate.
    if cfg.clustering.plot_per_combo:
        for n in n_neighbors_grid:
            for res in resolutions_grid:
                umap_key = f"umap_{n}_{res}"
                leiden_key = f"leiden_{n}_{res}"
                if umap_key not in adata.obsm or leiden_key not in adata.obs:
                    continue
                saved = adata.obsm.get("X_umap")
                adata.obsm["X_umap"] = adata.obsm[umap_key].copy()
                try:
                    fig, ax = plt.subplots(figsize=(6, 5))
                    drawn = _smart_plot_umap(
                        adata,
                        color=leiden_key,
                        ax=ax,
                        title=f"UMAP (n_neighbors={n}, resolution={res})",
                        cfg=cfg,
                        log=log,
                    )
                    if drawn:
                        save_figure(
                            fig,
                            os.path.join(fig_dir, f"leiden_grid_n{n}_r{res}"),
                            cfg=cfg,
                            dpi=cfg.plot.figure_dpi,
                            bbox_inches="tight",
                        )
                        log.info("    Plot saved: umap_grid_n%d_r%.1f", n, res)
                    plt.close(fig)
                except Exception as e:
                    log.warning("    Single-param UMAP plot save failed: %s", e)
                finally:
                    if saved is not None:
                        adata.obsm["X_umap"] = saved
    else:
        log.info("Per-combo UMAP plots skipped (plot_per_combo=False)")

    # ── 汇总 CSV ──
    df_summary = pd.DataFrame(results_summary)
    csv_path = os.path.join(table_dir, "param_grid_summary.csv")
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

    method = getattr(cfg.clustering, "cluster_selection_method", "multi_metric")

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

    # Resolve cluster_key from results (funnel uses funnel_{n}_{r}, full grid uses leiden_{n}_{r})
    leiden_col = None
    umap_col = f"umap_{best_n}_{best_r}"
    for r in results_summary:
        if r.get("n_neighbors") == best_n and r.get("resolution") == best_r:
            leiden_col = r.get("cluster_key")
            break
    if leiden_col is None:
        leiden_col = f"leiden_{best_n}_{best_r}"

    if leiden_col in adata.obs:
        adata.obs["leiden"] = adata.obs[leiden_col].copy()
        # Funnel mode: UMAP already in X_umap from re-validation; full grid: copy per-combo UMAP
        if umap_col in adata.obsm:
            adata.obsm["X_umap"] = adata.obsm[umap_col].copy()
        with timed_substep("Save checkpoint (post-selection)", log=log):
            _write_cluster_h5ad(adata, cfg, log)
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

    # Guard: multi-value min_dist sweep needs an O(n^2) trustworthiness matrix
    # (76k cells ≈ 92GB). Silhouette is computed on the UMAP embedding, so a
    # sweep also makes silhouette non-comparable across grid combos. Enforce
    # single-value grids on large datasets; degrade to the first value.
    if len(min_dist_grid) > 1 and adata.n_obs > 30_000:
        log.warning(
            "UMAP min_dist sweep disabled for n_obs=%d > 30k (O(n^2) trustworthiness OOM + "
            "silhouette non-comparability). Using min_dist=%s only.",
            adata.n_obs,
            min_dist_grid[0],
        )
        min_dist_grid = min_dist_grid[:1]

    # If PAGA init is enabled, compute PAGA backbone first
    if use_paga:
        log.info("Computing PAGA backbone for UMAP initialization...")
        sc.tl.paga(adata, groups="leiden")
        sc.pl.paga(adata, show=False)
        save_figure(
            None,
            os.path.join(fig_dir, "paga_backbone"),
            cfg=cfg,
            dpi=cfg.plot.figure_dpi,
            bbox_inches="tight",
        )
        plt.close()
        log.info("  PAGA backbone computed and saved")

    with timed_substep("UMAP sweep (min_dist × spread)", log=log):
        best_md, best_sp, umap_method_label, sweep_results = select_best_umap_params(
            adata,
            best_n,
            min_dist_grid,
            spread_grid,
            umap_method,
            cfg,
            use_rep,
            log,
            device=cfg.execution.device,
            metric=getattr(cfg.clustering, "umap_selection_metric", "trustworthiness"),
        )

    # Rebuild UMAP with selected params and re-save checkpoint.
    # Warm-start from the existing X_umap (left by select_best_umap_params sweep)
    # instead of paying for spectral initialization again. Major win on
    # 1M-cell datasets where spectral init alone is ~10-15 min.
    _final_init = adata.obsm.get("X_umap")
    _final_init_source = "existing X_umap (warm start)"
    if _final_init is None:
        _final_init = "paga" if use_paga else "spectral"
        _final_init_source = str(_final_init)
    log.info(
        "Rebuilding UMAP with selected params (min_dist=%.2f, spread=%.1f) [%s], init=%s...",
        best_md,
        best_sp,
        umap_method_label,
        _final_init_source,
    )
    try:
        gpu_umap(
            adata,
            log=log,
            device=cfg.execution.device,
            min_dist=best_md,
            spread=best_sp,
            init_pos=_final_init,
            maxiter=cfg.clustering.umap_maxiter,
            n_epochs=cfg.clustering.umap_n_epochs,
            random_state=cfg.execution.random_seed,
        )
    except Exception as e:
        # UMAP recompute failure is non-fatal (checkpoint was already written
        # post-selection above) — warn and keep the earlier 04_clustered.h5ad.
        log.warning("Final UMAP rebuild failed: %s", e)
    else:
        # Write must NOT be swallowed here: mode A needs the exception to
        # propagate so _write_cluster_h5ad's corrupt-copy deletion is honored.
        _write_cluster_h5ad(adata, cfg, log)
        log.info("Checkpoint saved with final UMAP: %s", cfg.cluster_h5ad)

    if sweep_results:
        # Summary CSV
        try:
            sweep_csv = os.path.join(table_dir, "umap_min_dist_sweep_summary.csv")
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
                    # Reuse sweep coords (saved by select_best_umap_params) — avoids
                    # re-running sc.tl.umap (cold spectral init costs ~47min/combo on 1M cells).
                    sweep_coords = r.get("coords")
                    title = f"min_dist={md}, spread={sp}"
                    if sweep_coords is None:
                        # Fallback (defensive): no coords saved — skip rather than recompute
                        ax.text(
                            0.5,
                            0.5,
                            f"{title}\n(no coords saved)",
                            ha="center",
                            va="center",
                            transform=ax.transAxes,
                        )
                    else:
                        _plot_umap_from_coords(
                            sweep_coords,
                            adata,
                            color="leiden",
                            ax=ax,
                            title=title,
                            cfg=cfg,
                            log=log,
                            legend_fontsize=8,
                        )
                except Exception:
                    ax.text(0.5, 0.5, "Error", ha="center", va="center", transform=ax.transAxes)
            for j in range(len(sweep_results), len(axes_flat)):
                axes_flat[j].axis("off")
            fig.tight_layout()
            save_figure(
                fig,
                os.path.join(fig_dir, "min_dist_sweep"),
                cfg=cfg,
                dpi=cfg.plot.figure_dpi,
                bbox_inches="tight",
            )
            plt.close(fig)
            log.info("UMAP min_dist comparison plot saved")
        except Exception as e:
            log.warning("UMAP comparison plot generation failed: %s", e)

    # ── Batch UMAP 诊断 ──
    if getattr(cfg.clustering, "umap_color_by_batch", False):
        batch_key = (
            getattr(cfg.clustering, "batch_key_override", None) or cfg.integration.batch_key
        )
        if batch_key not in adata.obs:
            log.warning(
                "batch_key '%s' not found in adata.obs, skipping batch UMAP diagnostics", batch_key
            )
        else:
            batches = adata.obs[batch_key].unique()
            n_batches = len(batches)
            log.info("Batch UMAP diagnostics: key=%s, n_batches=%d", batch_key, n_batches)
            # Guard: skip if fewer than 2 batches (single-batch = meaningless)
            if n_batches < 2:
                log.info(
                    "  n_batches=%d < 2, skipping batch UMAP diagnostics (single-batch data)",
                    n_batches,
                )
            else:
                # 1) Colored UMAP — subsample-aware (stratified) renderer
                try:
                    fig_b, ax_b = plt.subplots(figsize=(6, 5))
                    drawn = _smart_plot_umap(
                        adata, batch_key, ax_b, f"UMAP colored by {batch_key}", cfg, log
                    )
                    if drawn:
                        save_figure(
                            None,
                            os.path.join(fig_dir, "batch_colored_umap"),
                            cfg=cfg,
                            dpi=cfg.plot.figure_dpi,
                            bbox_inches="tight",
                        )
                        log.info("  Batch-colored UMAP saved")
                    plt.close(fig_b)
                except Exception as e:
                    log.warning("  Batch-colored UMAP failed: %s", e)
                # 2) Split-by-batch faceted UMAP (degrade if >12 batches)
                if n_batches <= 12:
                    try:
                        n_cols = min(4, n_batches)
                        n_rows = int(np.ceil(n_batches / n_cols))
                        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
                        axes_flat = axes.ravel() if n_batches > 1 else [axes]
                        for i, batch_val in enumerate(batches):
                            ax = axes_flat[i]
                            mask = adata.obs[batch_key] == batch_val
                            rng = np.random.RandomState(getattr(cfg.execution, "random_seed", 42))
                            bg_cap = min(50000, adata.n_obs)
                            bg_idx = rng.choice(adata.n_obs, bg_cap, replace=False)
                            ax.scatter(
                                adata.obsm["X_umap"][bg_idx, 0],
                                adata.obsm["X_umap"][bg_idx, 1],
                                c="lightgray",
                                s=1,
                                alpha=0.3,
                                rasterized=True,
                            )
                            fg_all = np.where(mask)[0]
                            fg_sub = fg_all[
                                _plot_subsample_idx(
                                    adata.obs.loc[mask, "leiden"].values,
                                    50000,
                                    rng,
                                    _plot_caps(cfg),
                                )
                            ]
                            ax.scatter(
                                adata.obsm["X_umap"][fg_sub, 0],
                                adata.obsm["X_umap"][fg_sub, 1],
                                c=adata.obs.loc[fg_sub, "leiden"].astype("category").cat.codes,
                                cmap=cfg.plot.palette.categorical,
                                s=3,
                                alpha=0.8,
                                rasterized=True,
                            )
                            ax.set_title(
                                f"{batch_key}={batch_val} ({mask.sum()} cells)", fontsize=9
                            )
                            ax.set_xticks([])
                            ax.set_yticks([])
                        for j in range(n_batches, len(axes_flat)):
                            axes_flat[j].axis("off")
                        fig.tight_layout()
                        save_figure(
                            fig,
                            os.path.join(fig_dir, "batch_faceted_umap"),
                            cfg=cfg,
                            dpi=cfg.plot.figure_dpi,
                            bbox_inches="tight",
                        )
                        log.info("  Split-by-batch faceted UMAP saved")
                    except Exception as e:
                        log.warning("  Split-by-batch UMAP failed: %s", e)
                    finally:
                        try:
                            plt.close(fig)
                        except Exception:
                            pass
                else:
                    log.info("  n_batches=%d > 12, skipping faceted UMAP", n_batches)
                # 3) Leiden × batch crosstab heatmap (pure matplotlib — zero new deps)
                try:
                    ct = pd.crosstab(adata.obs["leiden"], adata.obs[batch_key])
                    fig, ax = plt.subplots(
                        figsize=(max(6, n_batches * 0.8), max(5, ct.shape[0] * 0.3))
                    )
                    im = ax.imshow(ct.values, aspect="auto", cmap=cfg.plot.palette.dotplot_fill)
                    for i in range(ct.shape[0]):
                        for j in range(ct.shape[1]):
                            val = ct.values[i, j]
                            ax.text(
                                j,
                                i,
                                str(val),
                                ha="center",
                                va="center",
                                fontsize=8,
                                color="white" if val > ct.values.max() / 2 else "black",
                            )
                    ax.set_xticks(range(ct.shape[1]))
                    ax.set_xticklabels(ct.columns, rotation=45, ha="right")
                    ax.set_yticks(range(ct.shape[0]))
                    ax.set_yticklabels(ct.index)
                    ax.set_xlabel(batch_key)
                    ax.set_ylabel("Leiden cluster")
                    ax.set_title(f"Cluster \u00d7 batch mixing ({batch_key})")
                    plt.colorbar(im, ax=ax, label="cell count")
                    fig.tight_layout()
                    save_figure(
                        fig,
                        os.path.join(fig_dir, "batch_mixing_heatmap"),
                        cfg=cfg,
                        dpi=cfg.plot.figure_dpi,
                        bbox_inches="tight",
                    )
                    log.info("  Batch mixing heatmap saved")
                except Exception as e:
                    log.warning("  Batch mixing heatmap failed: %s", e)
                finally:
                    try:
                        plt.close(fig)
                    except Exception:
                        pass

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
                        _smart_plot_umap(
                            adata,
                            color=leiden_key,
                            ax=ax,
                            title=f"n={n}, r={res}",
                            cfg=cfg,
                            log=log,
                            legend_fontsize=5,
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
        save_figure(
            fig,
            os.path.join(fig_dir, "param_grid_summary"),
            cfg=cfg,
            dpi=cfg.plot.figure_dpi,
            bbox_inches="tight",
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
                    _smart_plot_umap(
                        adata,
                        color=key,
                        ax=axes[i],
                        title=key,
                        cfg=cfg,
                        log=log,
                        legend_fontsize=8,
                    )
                for j in range(len(res_keys), len(axes)):
                    axes[j].axis("off")
                fig.tight_layout()
                save_figure(
                    fig,
                    os.path.join(fig_dir, f"leiden_multires_n{n}"),
                    cfg=cfg,
                    dpi=cfg.plot.figure_dpi,
                    bbox_inches="tight",
                )
                plt.close(fig)
                log.info("  Multi-resolution Leiden plot (n=%d) saved", n)
            except Exception as e:
                log.warning("  Multi-resolution comparison plot (n=%d) failed: %s", n, e)

    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
