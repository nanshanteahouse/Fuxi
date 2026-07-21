#!/usr/bin/env python3
"""
Step 12: Gene Regulatory Network (GRN) analysis for spatial transcriptomics
=============================================================================
Pseudobulk aggregation per cell_type across spots → CollecTRI regulons → ULM activity inference → heatmap + tables.

Input:  05_annotated.h5ad (requires cell_type column)
Output:
  {table_dir}/grn/tf_activity_per_cell_type.csv   — TF activity scores
  {table_dir}/grn/tf_activity_pvals.csv           — associated p-values
  {figure_dir}/grn/tf_activity_heatmap.png        — dendrogram + heatmap
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import dendrogram, linkage

from core.pipeline.grn import compute_tf_relevance
from core.utils import resolve_config, setup_logger


def build_pseudobulk(adata, group_col: str, use_raw: bool = True, log=None) -> pd.DataFrame:
    import scipy.sparse as sp

    if group_col not in adata.obs:
        if log:
            log.warning("%s not in adata.obs - using 'leiden'", group_col)
        group_col = "leiden"

    groups = adata.obs[group_col].values
    src = adata.raw if use_raw and adata.raw else adata
    x = src.X
    var_names = src.var_names
    is_sparse = sp.issparse(x)

    if log:
        log.info(
            "Pseudobulk: %d spots -> %d groups",
            adata.n_obs,
            len(adata.obs[group_col].cat.categories),
        )

    unique_groups = adata.obs[group_col].cat.categories
    n_groups = len(unique_groups)
    n_genes = len(var_names)

    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    group_indices = np.array([group_to_idx[g] for g in groups])

    pseudo = np.zeros((n_groups, n_genes), dtype=np.float64)
    for g_idx in range(n_groups):
        mask = group_indices == g_idx
        if mask.any():
            subset = x[mask]
            if is_sparse:
                pseudo[g_idx] = subset.mean(axis=0).A1
            else:
                pseudo[g_idx] = subset.mean(axis=0)

    pseudo = np.log1p(pseudo)

    df = pd.DataFrame(pseudo, index=unique_groups, columns=var_names)
    if log:
        log.info("  Pseudobulk matrix: %d x %d", n_groups, n_genes)
    return df


def filter_regulon_net(net: pd.DataFrame, min_genes: int = 5, log=None) -> pd.DataFrame:
    n_before = net["source"].nunique()
    gene_counts = net.groupby("source")["target"].nunique()
    keep = gene_counts[gene_counts >= min_genes].index
    net_filt = net[net["source"].isin(keep)].copy()
    if log:
        log.info("Regulon filter (>=%d targets): %d -> %d TFs", min_genes, n_before, len(keep))
    return net_filt


def run_grn(pseudo_df: pd.DataFrame, net: pd.DataFrame, log) -> tuple:
    import decoupler as dc

    log.info(
        "Running ULM enrichment on %d cell groups x %d genes",
        pseudo_df.shape[0],
        pseudo_df.shape[1],
    )

    avail_genes = set(pseudo_df.columns)
    net = net[net["target"].isin(avail_genes)].copy()
    log.info("  Regulon edges covering available genes: %d", len(net))

    estimates, pvals = dc.mt.ulm(pseudo_df, net, verbose=False)

    est_df = pd.DataFrame(estimates, index=pseudo_df.index, columns=estimates.columns)
    pval_df = pd.DataFrame(pvals, index=pseudo_df.index, columns=estimates.columns)

    log.info("  Activity matrix: %d cell types x %d TFs", est_df.shape[0], est_df.shape[1])
    return est_df, pval_df, net


def top_variable_tfs(
    estimates_df: pd.DataFrame,
    n_top: int,
    log,
    mode: str = "off",
    tf_annotation: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if mode == "off":
        var = estimates_df.var(axis=0)
        top = var.nlargest(n_top).index
        log.info("Top %d TFs by variance: %s", n_top, ", ".join(top[:20].tolist()))
        return estimates_df[top]

    if mode == "soft":
        var = estimates_df.var(axis=0)
        top = var.nlargest(n_top).index
        log.info("Top %d TFs by variance: %s", n_top, ", ".join(top[:20].tolist()))

        if tf_annotation is not None:
            top_tf_set = set(top)
            overlap_count = tf_annotation[
                tf_annotation["tf"].isin(top_tf_set) & (tf_annotation["kb_overlap_ratio"] > 0)
            ].shape[0]
            log.info("soft mode: %d/%d top TFs have KB marker overlap", overlap_count, len(top))
        else:
            log.warning("soft mode: tf_annotation not provided — skipping KB overlap logging")

        return estimates_df[top]

    if mode == "hard":
        if tf_annotation is None:
            log.warning("tf_annotation not provided for hard mode — falling back to off mode")
            var = estimates_df.var(axis=0)
            top = var.nlargest(n_top).index
            log.info("Top %d TFs by variance: %s", n_top, ", ".join(top[:20].tolist()))
            return estimates_df[top]

        var = estimates_df.var(axis=0)
        var_rank = var.rank(ascending=False)

        tf_ratio_map = dict(zip(tf_annotation["tf"], tf_annotation["kb_overlap_ratio"]))
        kb_metric = pd.Series(
            {
                tf: tf_ratio_map.get(tf, 0.0) * estimates_df[tf].abs().mean()
                for tf in estimates_df.columns
            }
        )
        kb_rank = kb_metric.rank(ascending=False)

        combined = var_rank + kb_rank
        top = combined.nsmallest(n_top).index

        log.info("hard mode: selecting top %d TFs by combined variance+KB rank", n_top)
        log.info("Top %d TFs (hard mode): %s", n_top, ", ".join(top[:20].tolist()))

        return estimates_df[top]

    log.warning("Unknown mode '%s' — falling back to off mode", mode)
    var = estimates_df.var(axis=0)
    top = var.nlargest(n_top).index
    return estimates_df[top]


def export_results(estimates_df, top_df, pvals_df, net_top, cfg, log, kb_markers=None):
    table_dir = os.path.join(cfg.table_dir, "grn")
    os.makedirs(table_dir, exist_ok=True)

    if kb_markers is not None:
        _, tf_ann = compute_tf_relevance(estimates_df, net_top, kb_markers, log)
        ann_path = os.path.join(table_dir, "tf_annotation_table.csv")
        tf_ann.to_csv(ann_path, index=False)
        log.info("Exported: %s (%d TFs)", ann_path, len(tf_ann))

    path = os.path.join(table_dir, "tf_activity_per_cell_type.csv")
    estimates_df.to_csv(path)
    log.info("Exported: %s", path)

    path = os.path.join(table_dir, "tf_activity_pvals.csv")
    pvals_df.to_csv(path)
    log.info("Exported: %s", path)

    path = os.path.join(table_dir, "tf_target_edges.csv")
    net_top.to_csv(path, index=False)
    log.info("Exported: %s (%d edges)", path, len(net_top))

    target_counts = (
        net_top.groupby("source")["target"]
        .nunique()
        .reset_index()
        .rename(columns={"source": "tf", "target": "n_targets"})
        .sort_values("n_targets", ascending=False)
    )
    path = os.path.join(table_dir, "tf_target_counts.csv")
    target_counts.to_csv(path, index=False)
    log.info("Exported: %s (%d TFs)", path, len(target_counts))

    top_regulon_path = os.path.join(table_dir, "tf_activity_top_regulons.csv")
    top_df.to_csv(top_regulon_path)
    log.info("Exported: %s (%d cell types x %d TFs)", top_regulon_path, *top_df.shape)

    if net_top.empty:
        log.warning("Top-TF edge list is empty — no edges to export")


def plot_heatmap(top_df, cfg, log):
    fig_dir = os.path.join(cfg.figure_dir, "grn")
    os.makedirs(fig_dir, exist_ok=True)

    n_tfs = top_df.shape[1]
    n_cts = top_df.shape[0]
    if n_tfs < 2 or n_cts < 2:
        log.warning("Too few TFs or cell types for heatmap - skipping")
        return

    data = top_df.values.T
    z_rows = linkage(data, method="ward")
    z_cols = linkage(data.T, method="ward")

    longest_tf = max(len(n) for n in top_df.columns)
    tf_label_w = longest_tf * 0.075
    left_pct = max(0.18, min(0.40, longest_tf * 0.008))

    tf_yunit = min(0.30, max(0.18, 6.0 / max(1, n_tfs)))
    tf_hm = max(0.28, tf_yunit * 0.7)
    ct_xunit = 0.45

    top_margin = 0.8
    right_pad = tf_label_w + 0.55
    bottom_margin = 0.6

    left_margin = left_pct * max(5, n_cts * ct_xunit)
    heatmap_w = n_cts * ct_xunit
    heatmap_h = n_tfs * tf_hm

    fig_w = left_margin + heatmap_w + right_pad
    fig_h = top_margin + 0.5 + heatmap_h + bottom_margin

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[left_margin, heatmap_w],
        height_ratios=[0.5, heatmap_h],
        hspace=0.0,
        wspace=0.0,
        left=left_margin / fig_w,
        right=(left_margin + heatmap_w) / fig_w,
        top=(top_margin + 0.5 + heatmap_h) / fig_h,
        bottom=bottom_margin / fig_h,
    )

    ax_row = fig.add_subplot(gs[1, 0])
    d_rows = dendrogram(
        z_rows,
        ax=ax_row,
        orientation="left",
        link_color_func=lambda k: "#555555",
        above_threshold_color="#bbbbbb",
        no_labels=True,
    )
    row_idx = d_rows["leaves"]
    ax_row.invert_xaxis()
    ax_row.set_xticks([])
    ax_row.set_yticks([])
    for s in ax_row.spines.values():
        s.set_visible(False)

    ax_col = fig.add_subplot(gs[0, 1])
    d_cols = dendrogram(
        z_cols,
        ax=ax_col,
        orientation="top",
        link_color_func=lambda k: "#555555",
        above_threshold_color="#bbbbbb",
        no_labels=True,
    )
    col_idx = d_cols["leaves"]
    ax_col.set_yticks([])
    ax_col.set_xticks([])
    for s in ax_col.spines.values():
        s.set_visible(False)

    data_clust = data[row_idx, :][:, col_idx]
    tf_labels = top_df.columns[row_idx]
    ct_labels = top_df.index[col_idx]

    ax_hm = fig.add_subplot(gs[1, 1])
    vabs = np.percentile(np.abs(data_clust), 90)
    im = ax_hm.imshow(
        data_clust, aspect="auto", cmap="RdBu_r", vmin=-vabs, vmax=vabs, interpolation="nearest"
    )

    ax_hm.set_xticks(range(n_cts))
    ax_hm.set_xticklabels(ct_labels, rotation=35, ha="right", fontsize=8.5)
    ax_hm.set_yticks(range(n_tfs))
    ax_hm.set_yticklabels(tf_labels, fontsize=7.0)
    ax_hm.yaxis.tick_right()
    ax_hm.tick_params(length=0, pad=3)

    cax = fig.add_axes(
        [
            (left_margin + heatmap_w + tf_label_w + 0.12) / fig_w,
            (bottom_margin + heatmap_h * 0.15) / fig_h,
            0.012,
            (heatmap_h * 0.45) / fig_h,
        ]
    )
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.set_label("Activity score", fontsize=8)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle(
        f"TF Activity (ULM) - Top {n_tfs} Regulons",
        fontsize=11,
        fontweight="bold",
        x=(left_margin + heatmap_w / 2) / fig_w,
        y=(top_margin + 0.5 + heatmap_h + 0.15) / fig_h,
        ha="center",
        va="bottom",
    )

    path = os.path.join(fig_dir, "tf_activity_heatmap.png")
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info("Heatmap saved: %s", path)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("12_grn", os.path.join(cfg.log_dir, "12_grn.log"))
    log.info("Step 12: GRN regulatory network analysis (spatial)")

    if not getattr(cfg.grn, "run", True):
        log.info("GRN disabled - skipping")
        return

    input_path = os.path.join(cfg.h5ad_dir, "05_annotated.h5ad")
    adata = sc.read(input_path)
    log.info("Loaded: %s - %d spots, %d genes", input_path, adata.n_obs, adata.n_vars)

    assert adata.raw is not None, "Raw counts missing — re-run from step 03"

    group_col = "cell_type" if "cell_type" in adata.obs else "leiden"
    log.info("Group column: %s (%d categories)", group_col, adata.obs[group_col].nunique())

    pseudo_df = build_pseudobulk(adata, group_col, use_raw=True, log=log)

    species = getattr(cfg.grn, "species", "human")
    log.info("Regulon: CollecTRI (%s)", species)
    net = dc.op.collectri(organism=species)
    net = net[net["weight"] > 0].copy()

    min_size = getattr(cfg.grn, "min_regulon_size", 5)
    net = filter_regulon_net(net, min_genes=min_size, log=log)

    est_df, pval_df, net_filtered = run_grn(pseudo_df, net, log)

    kb_markers = None
    n_top = min(getattr(cfg.grn, "n_top_regulons", 50), est_df.shape[1])

    if getattr(cfg.grn, "use_kb_relevance", False):
        tissue = getattr(cfg, "tissue", "") or ""
        if tissue and tissue != "unknown":
            from core.kb import load_all_kb_markers

            kb_markers = load_all_kb_markers(tissue)
            if kb_markers:
                log.info("Loaded %d KB markers for tissue '%s'", len(kb_markers), tissue)
            else:
                log.info("No KB markers found for tissue '%s'", tissue)
        else:
            log.info("No tissue configured, skipping KB marker loading")

    grn_mode = getattr(cfg.grn, "tissue_mode", "off")
    if grn_mode not in {"off", "soft", "hard"}:
        raise ValueError(f"Invalid grn_tissue_mode: '{grn_mode}'. Must be off, soft, or hard.")

    if grn_mode in ("soft", "hard") and kb_markers is not None:
        annotated_activity_df, tf_ann = compute_tf_relevance(est_df, net_filtered, kb_markers, log)
    else:
        tf_ann = None

    top_df = top_variable_tfs(est_df, n_top, log, mode=grn_mode, tf_annotation=tf_ann)

    top_tfs = set(top_df.columns)
    net_top = net_filtered[net_filtered["source"].isin(top_tfs)].copy()
    log.info("Top-TF edges: %d (from %d total filtered edges)", len(net_top), len(net_filtered))

    if grn_mode == "hard" and tf_ann is not None and getattr(cfg.grn, "export_filtered", False):
        relevant_tfs = set(tf_ann[tf_ann["kb_overlap_ratio"] > 0]["tf"])
        top_df = top_df[[c for c in top_df.columns if c in relevant_tfs]]
        log.info("hard+export_filtered: reduced to %d tissue-relevant TFs", top_df.shape[1])

    export_results(est_df, top_df, pval_df, net_top, cfg, log, kb_markers=kb_markers)
    plot_heatmap(top_df, cfg, log)

    elapsed = time.time() - t0
    log.info("Step 12 complete (took %.1fs).", elapsed)


if __name__ == "__main__":
    main()
