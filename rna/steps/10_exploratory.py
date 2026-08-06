#!/usr/bin/env python3
"""
Step 10: 探索性分析
=======================
  1. 细胞组成随发育阶段的变化 (堆叠柱状图)
  2. UMAP 上的 QC 指标检查
  3. 已知标记基因的 UMAP 表达
  4. 聚类大小统计

输入: 04_clustered.h5ad (需要 Stage 05 运行后以获得 cell_type 注释)
输出: CSV 表格 + PNG 图片 (不修改 h5ad)
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

from core.utils import resolve_config, safe_plot, save_figure, setup_logger
from core.utils._optional import require_sccoda


def plot_composition(adata, group_col, stage_col, stage_order, fig_dir, table_dir, log, cfg=None):
    """绘制细胞类型随发育阶段的组成变化堆积图"""
    if group_col not in adata.obs or stage_col not in adata.obs:
        log.warning("Missing %s or %s, skipping composition plot", group_col, stage_col)
        return
    ct_counts = (
        adata.obs.groupby([stage_col, group_col], observed=True).size().reset_index(name="count")
    )
    ct_pivot = ct_counts.pivot_table(
        index=stage_col, columns=group_col, values="count", fill_value=0
    )
    avail_stages = [s for s in stage_order if s in ct_pivot.index]
    if not avail_stages:
        avail_stages = list(ct_pivot.index)
    ct_pivot = ct_pivot.reindex(avail_stages)
    ct_pivot = ct_pivot.div(ct_pivot.sum(axis=1), axis=0)

    n_types = ct_pivot.shape[1]
    _cmap = cfg.plot.palette.categorical if cfg else "tab20"
    colors = plt.get_cmap(_cmap)(np.linspace(0, 1, min(n_types, 20)))
    if n_types > 20:
        n_tile = int(np.ceil(n_types / 20))
        colors = np.vstack([colors] * n_tile)[:n_types]
    colors = [tuple(c) for c in colors]

    fig, ax = plt.subplots(figsize=(max(10, len(avail_stages) * 1.5), 6))
    ct_pivot.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.8)
    ax.set_xlabel("Developmental stage")
    ax.set_ylabel("Fraction of cells")
    ax.set_title(f"Cluster composition by stage ({group_col})")
    ax.legend(
        title=group_col, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title_fontsize=9
    )
    fig.tight_layout()
    _dpi = cfg.plot.figure_dpi if cfg else 150
    save_figure(
        fig,
        os.path.join(fig_dir, f"composition_by_stage_{group_col}"),
        cfg=cfg,
        dpi=_dpi,
        bbox_inches="tight",
    )
    plt.close(fig)
    log.info("  Composition plot saved: composition_by_stage_%s.png", group_col)

    # 导出 CSV
    ct_pivot.to_csv(os.path.join(table_dir, f"composition_by_stage_{group_col}.csv"))
    log.info("  Composition table exported")


def run_sccoda_composition(adata, cfg, log, fig_dir=None, table_dir=None):
    """Run scCODA compositional analysis for differential cell-type abundance.

    Requires the ``sccoda`` package. This function sets up the data and
    triggers the composition test. The actual MCMC sampling is handled
    by the scCODA backend.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cell-type labels and sample metadata.
    cfg : Config
        Pipeline configuration (uses ``cfg.exploratory.sccoda``).
    log : logging.Logger
        Logger instance.
    """
    require_sccoda()
    log.info("scCODA composition analysis requires manual MCMC execution")
    log.info("  -- configure and run via sccoda CLI or Python API")
    log.info(
        "  sample_col=%s  condition_col=%s  ref_cell_type=%s",
        cfg.exploratory.sccoda.sample_col,
        cfg.exploratory.sccoda.condition_col,
        cfg.exploratory.sccoda.reference_cell_type,
    )


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("10_exploratory", os.path.join(cfg.log_dir, "10_exploratory.log"))
    log.info("Step 10: Exploratory analysis")

    # Prefer annotated h5ad (has cell_type), fall back to clustered
    annotated_path = os.path.join(cfg.h5ad_dir, "05_annotated.h5ad")
    input_path = annotated_path if os.path.exists(annotated_path) else cfg.cluster_h5ad
    if not os.path.exists(annotated_path):
        log.warning("05_annotated.h5ad not found, falling back to: %s", cfg.cluster_h5ad)
    adata = sc.read(input_path)
    log.info("Loaded: %s — %d cells", input_path, adata.n_obs)

    # 轻量绘图对象（参照 04_cluster 的 plot_step04_figures 模式）：
    # 分层降采样保留稀有簇可见性；raw 只留 marker 基因，避免全基因副本
    marker_genes = []
    for genes in cfg.marker.marker_dict.values():
        marker_genes.extend([g for g in genes if g in adata.raw.var_names][:2])
    marker_genes = list(dict.fromkeys(marker_genes))

    if adata.n_obs > 20_000:
        _rng = np.random.default_rng(0)
        _strat = "cell_type" if "cell_type" in adata.obs else "leiden"
        _codes = pd.Categorical(adata.obs[_strat]).codes
        _parts = []
        for _code in np.unique(_codes):
            _sel = np.where(_codes == _code)[0]
            _n = len(_sel)
            if _n <= 500:
                _parts.append(_sel)
            else:
                _cap = 500 if _n <= 50_000 else 1000
                _parts.append(_rng.choice(_sel, min(_cap, _n), replace=False))
        _idx = np.concatenate(_parts)
        if len(_idx) > 20_000:
            _idx = _rng.choice(_idx, 20_000, replace=False)
        import anndata as _ad

        plot_adata = _ad.AnnData(
            X=adata[_idx].X,
            obs=adata.obs.iloc[_idx],
            var=adata.var.copy(),
            obsm={"X_umap": adata.obsm["X_umap"][_idx]},
        )
        if marker_genes:
            _raw_sel = adata.raw[_idx][:, marker_genes]
            plot_adata.raw = _ad.AnnData(X=_raw_sel.X, var=_raw_sel.var)
        log.info("  Stratified downsample to %d cells for UMAP plots", plot_adata.n_obs)
    else:
        plot_adata = adata

    fig_dir = os.path.join(cfg.figure_dir, "10_exploratory")
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    sc.settings.autoshow = False
    csv_dir = os.path.join(cfg.table_dir, "10_exploratory")
    os.makedirs(csv_dir, exist_ok=True)

    # 1. 细胞组成
    group_by = ["cell_type", "cell_type_sub", "leiden"]
    for g in group_by:
        if g in adata.obs:
            plot_composition(
                adata,
                g,
                "stage" if "stage" in adata.obs else "sample",
                cfg.sample_meta.stage_order,
                fig_dir,
                csv_dir,
                log,
                cfg=cfg,
            )

    # 1b. scCODA composition test (optional)
    if cfg.exploratory.composition_test == "sccoda":
        run_sccoda_composition(adata, cfg, log, fig_dir=fig_dir, table_dir=csv_dir)

    # 2. UMAP: QC 指标
    qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    qc_metrics = [m for m in qc_metrics if m in adata.obs]
    if qc_metrics:
        safe_plot(
            sc.pl.umap,
            plot_adata,
            color=qc_metrics,
            show=False,
            save="qc_umap",
            vmax="p99",
            ncols=3,
            cfg=cfg,
        )
        plt.close("all")

    # 3. UMAP: 标记基因
    all_markers = marker_genes
    if all_markers:
        n_markers = len(all_markers)
        batch_size = 12
        for batch_start in range(0, n_markers, batch_size):
            batch = all_markers[batch_start : batch_start + batch_size]
            safe_plot(
                sc.pl.umap,
                plot_adata,
                color=batch,
                use_raw=True,
                show=False,
                save=f"marker_umap_batch{batch_start}",
                vmax="p99",
                ncols=4,
                cfg=cfg,
            )
            plt.close("all")

    # 4. 标记基因 dotplot
    if all_markers:
        group_col = "cell_type" if "cell_type" in adata.obs else "leiden"
        safe_plot(
            sc.pl.dotplot,
            adata,
            var_names=all_markers,
            groupby=group_col,
            show=False,
            save="marker_dotplot",
            cfg=cfg,
        )
        plt.close("all")

    # 5. 聚类大小统计
    for group_col in ["cell_type", "leiden"]:
        if group_col not in adata.obs:
            continue
        sizes = adata.obs[group_col].value_counts().sort_index()
        log.info("  %s size distribution:", group_col)
        for label, cnt in sizes.items():
            log.info("    %s: %d cells (%.1f%%)", label, cnt, 100 * cnt / adata.n_obs)
        sizes.to_csv(os.path.join(csv_dir, f"{group_col}_sizes.csv"), header=["n_cells"])

    # 6. 额外元数据分组可视化
    # 来源 A: 用户显式配置的 meta_columns
    extra_cols = set()
    for obs_col in getattr(cfg.sample_meta, "meta_columns", {}).values():
        if obs_col and obs_col in adata.obs:
            extra_cols.add(obs_col)
    # 来源 B: 用户手动指定的 step10_groupby 覆盖
    for obs_col in getattr(cfg.marker, "step10_groupby", []):
        if obs_col in adata.obs:
            extra_cols.add(obs_col)
    # 来源 C: 自动发现非 pipeline 分类列 (meta_columns 覆盖不到的)
    pipeline_prefixes = (
        "n_genes",
        "log1p_",
        "total_counts",
        "pct_counts",
        "leiden_",
        "doublet_",
        "predicted_",
        "annot_",
        "marker_",
        "scater_",
        "X_",
        "umap_",
    )
    auto_skip = {
        "Barcode",
        "sample",
        "stage",
        "batch",
        "leiden",
        "cell_type",
        "cell_type_sub",
        "cell_state",
        "log_genes_per_umi",
    }
    for col in adata.obs.columns:
        if col in extra_cols or col in auto_skip:
            continue
        if col.startswith(pipeline_prefixes):
            continue
        n_unique = adata.obs[col].nunique()
        if 2 <= n_unique <= 50:
            extra_cols.add(col)
    # 系统自动生成的列
    if "predicted_sex" in adata.obs:
        extra_cols.add("predicted_sex")

    for col in sorted(extra_cols):
        if col not in adata.obs:
            continue
        n_unique = adata.obs[col].nunique()
        if n_unique > 50:
            continue

        # Distribution CSV
        sizes = adata.obs[col].value_counts().sort_index()
        log.info("  %s distribution:", col)
        for label, cnt in sizes.items():
            log.info("    %s: %d cells (%.1f%%)", label, cnt, 100 * cnt / adata.n_obs)
        sizes.to_csv(os.path.join(csv_dir, f"{col}_sizes.csv"), header=["n_cells"])

        # UMAP
        safe_plot(
            sc.pl.umap,
            plot_adata,
            color=col,
            show=False,
            legend_loc="on data" if len(sizes) < 30 else "right margin",
            save=f"umap_{col}",
            cfg=cfg,
        )
        plt.close("all")

    log.info("Step 10 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
