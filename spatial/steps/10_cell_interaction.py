#!/usr/bin/env python3
"""
Step 10: Spatial Cell-Cell Interaction (CCI) analysis via LIANA+ bivariate metrics
====================================================================================
Spatial bivariate ligand-receptor co-expression analysis using local spatial metrics
(Cosine, Jaccard, Pearson, Spearman) and global Moran's R with permutation testing.

Method:
  * liana.method.bivariate — local + global spatial bivariate metrics

Input:  05_annotated.h5ad (requires spatial_connectivities in obsp + cell_type column)
Output:
  {table_dir}/10_cell_interaction/cci_spatial_interactions.csv    — full interaction scores
  {table_dir}/10_cell_interaction/cci_spatial_top.csv             — top N significant pairs
  {figure_dir}/10_cell_interaction/cci_spatial_heatmap.png        — spatial interaction heatmap
  {figure_dir}/10_cell_interaction/cci_spatial_dotplot.png        — top interactions dotplot

Dependencies: pip install liana>=1.0.0
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


# ═══════════════════════════════════════════════════════════════════════
#  Export
# ═══════════════════════════════════════════════════════════════════════

def export_spatial_results(lr_res, top_df, CFG, log):
    """Save spatial CCI interaction tables to CSV."""
    table_dir = os.path.join(CFG.table_dir, "10_cell_interaction")
    os.makedirs(table_dir, exist_ok=True)

    path = os.path.join(table_dir, "cci_spatial_interactions.csv")
    lr_res.to_csv(path, index=False)
    log.info("Exported: %s (%d rows)", path, len(lr_res))

    path = os.path.join(table_dir, "cci_spatial_top.csv")
    top_df.to_csv(path, index=False)
    log.info("Exported: %s (%d rows)", path, len(top_df))


# ═══════════════════════════════════════════════════════════════════════
#  Plots
# ═══════════════════════════════════════════════════════════════════════

def plot_spatial_heatmap(top_df, CFG, log):
    """Heatmap showing Moran's R scores for top LR interactions."""
    fig_dir = os.path.join(CFG.figure_dir, "10_cell_interaction")
    os.makedirs(fig_dir, exist_ok=True)

    # Bivariate output: ligand, receptor, morans -- no source/target columns
    if "ligand" not in top_df.columns or "receptor" not in top_df.columns:
        log.warning("Missing ligand/receptor columns -- skipping spatial heatmap")
        return

    score_col = "morans" if "morans" in top_df.columns else top_df.columns[0]
    pivot = top_df.pivot_table(
        index="ligand", columns="receptor",
        values=score_col, aggfunc="mean",
    )
    pivot = pivot.fillna(0)

    if pivot.empty or pivot.shape[0] < 2 or pivot.shape[1] < 2:
        log.warning("Not enough LR pairs for heatmap -- skipping")
        return

    n_rows, n_cols = pivot.shape
    fig_w = max(6, n_cols * 0.5 + 2.5)
    fig_h = max(4, n_rows * 0.4 + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xlabel("Receptor")
    ax.set_ylabel("Ligand")
    ax.set_title(
        f"Top {len(top_df)} Spatial CCI ({score_col})"
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("N interactions")

    fig.tight_layout()
    path = os.path.join(fig_dir, "cci_spatial_heatmap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", path)


def plot_spatial_dotplot(top_df, CFG, log):
    """Dotplot of top LR pairs with global Moran's R or similarity scores."""
    fig_dir = os.path.join(CFG.figure_dir, "10_cell_interaction")
    os.makedirs(fig_dir, exist_ok=True)

    # Determine value column for color
    score_candidates = ["morans", "morans_r", "global_morans", "magnitude", "magnitude_rank"]
    score_col = None
    for candidate in score_candidates:
        if candidate in top_df.columns:
            score_col = candidate
            break
    if score_col is None:
        for col in top_df.columns:
            if "rank" in col.lower() or "score" in col.lower() or "moran" in col.lower():
                score_col = col
                break
    if score_col is None:
        log.warning("No scoring column found -- skipping dotplot")
        return

    # Build interaction label
    top_df = top_df.copy()
    if "ligand_complex" in top_df.columns and "receptor_complex" in top_df.columns:
        top_df["interaction"] = (
            top_df["ligand_complex"].astype(str) + "_"
            + top_df["receptor_complex"].astype(str)
        )
    elif "ligand" in top_df.columns and "receptor" in top_df.columns:
        top_df["interaction"] = (
            top_df["ligand"].astype(str) + "_" + top_df["receptor"].astype(str)
        )
    else:
        log.warning("No ligand/receptor columns -- skipping dotplot")
        return

    # Bivariate has no source/target -- bar chart of top LR pairs
    n_interactions = min(len(top_df), 20)
    plot_df = top_df.head(n_interactions).iloc[::-1]  # reverse for bottom-up bars

    fig_w = max(8, n_interactions * 0.35 + 3.0)
    fig_h = max(4, n_interactions * 0.3 + 1.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    colors = plt.cm.RdYlBu_r((plot_df[score_col].values - plot_df[score_col].min())
                              / (plot_df[score_col].max() - plot_df[score_col].min() + 0.001))
    ax.barh(range(n_interactions), plot_df[score_col].values, color=colors)
    ax.set_yticks(range(n_interactions))
    ax.set_yticklabels(plot_df["interaction"].values, fontsize=7)
    ax.set_xlabel(score_col)
    ax.set_title(f"Top {n_interactions} Spatial CCI ({score_col})")

    fig.tight_layout()
    path = os.path.join(fig_dir, "cci_spatial_dotplot.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", path)


# ═══════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("10_cci", os.path.join(CFG.log_dir, "10_cci.log"))
    log.info("Step 10: Spatial Cell-Cell Interaction (CCI) via LIANA+ bivariate")

    # ── Gate check ──────────────────────────────────────────────────────
    if not getattr(CFG, "run_cci", True):
        log.info("run_cci=False — skipping")
        return

    # ── Load input ──────────────────────────────────────────────────────
    adata_path = CFG.annotated_h5ad
    if not os.path.exists(adata_path):
        log.error("Input not found: %s (run Step 05 first)", adata_path)
        sys.exit(1)

    log.info("Loading: %s", adata_path)
    adata = sc.read(adata_path)

    # ── Validate spatial prerequisites ──────────────────────────────────
    if "spatial_connectivities" not in adata.obsp:
        log.error(
            "'spatial_connectivities' not found in adata.obsp. "
            "Run spatial neighbors construction (Step 03) first."
        )
        sys.exit(1)

    n_spots = adata.n_obs
    avg_degree = adata.obsp["spatial_connectivities"].sum() / n_spots
    log.info("Spatial connectivities: %d spots, avg degree %.1f", n_spots, avg_degree)

    # Optionally rebuild spatial neighbors with custom radius
    cci_distance = getattr(CFG, "cci_spatial_distance", 0.0)
    if cci_distance > 0:
        import squidpy as sq
        log.info("Rebuilding spatial neighbors with radius=%.0f px", cci_distance)
        sq.gr.spatial_neighbors(
            adata, radius=cci_distance, key_added="spatial_connectivities",
        )
        new_avg = adata.obsp["spatial_connectivities"].sum() / n_spots
        log.info("  Updated avg degree: %.1f", new_avg)

    # ── Resolve grouping column for source/target labels ────────────────
    group_col = "cell_type" if "cell_type" in adata.obs.columns else "leiden"
    log.info("Cell type labels: %s (%d groups)", group_col, adata.obs[group_col].nunique())

    # LIANA bivariate needs obs labels for source/target annotation.
    # Ensure the grouping column is set as the obs label key.
    if group_col in adata.obs.columns:
        adata.obs["cci_label"] = adata.obs[group_col].astype(str)

    # ── Run CCI spatial analysis ────────────────────────────────────────
    from rna.utils.cell_interaction import (
        ensure_gene_symbols,
        run_cci_spatial,
        format_cci_results,
    )

    # ── Ensure gene symbols (LIANA uses HGNC symbols) ───────────────────
    adata = ensure_gene_symbols(adata, log=log)

    lr_res = run_cci_spatial(
        adata,
        resource_name=CFG.cci_lr_database,
        connectivity_key="spatial_connectivities",
        local_name="cosine",
        global_name="morans",
        n_perms=CFG.cci_permutations,
        log=log,
    )

    # ── Format & export ─────────────────────────────────────────────────
    # LIANA bivariate returns morans/morans_pvals columns;
    # sort by Moran's R descending (stronger spatial co-expression first)
    sort_col = "morans" if "morans" in lr_res.columns else lr_res.columns[0]
    top_df = format_cci_results(
        lr_res,
        n_top=CFG.cci_n_top_interactions,
        pval_col=sort_col,
        log=log,
    )

    export_spatial_results(lr_res, top_df, CFG, log)

    # ── Plot ────────────────────────────────────────────────────────────
    plot_spatial_heatmap(top_df, CFG, log)
    plot_spatial_dotplot(top_df, CFG, log)

    log.info("Step 10 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
