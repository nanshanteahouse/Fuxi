#!/usr/bin/env python3
"""
Step 12: Cell-Cell Interaction (CCI) analysis via LIANA+ permutation testing
==============================================================================
Permutation-based ligand-receptor interaction scoring across cell type groups.

Method:
  * rank_aggregate — consensus of multiple scoring methods (CellPhoneDB-like,
    logFC, NATMI, SingleCellSignalR, Connectome, CellChat, geometric_mean)

Input:  05_annotated.h5ad (requires cell_type column; uses .raw for counts)
Output:
  {table_dir}/12_cell_interaction/cci_interactions.csv         — full interaction scores
  {table_dir}/12_cell_interaction/cci_top_interactions.csv     — top N significant pairs
  {figure_dir}/12_cell_interaction/cci_heatmap.png             — LR interaction heatmap
  {figure_dir}/12_cell_interaction/cci_dotplot.png             — top interactions dotplot

Dependencies: pip install liana>=1.0.0
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

# Agg backend for headless environments
matplotlib.use("Agg")


# ═══════════════════════════════════════════════════════════════════════
#  Export
# ═══════════════════════════════════════════════════════════════════════

def export_results(lr_res, top_df, CFG, log):
    """Save CCI interaction tables to CSV."""
    table_dir = os.path.join(CFG.table_dir, "12_cell_interaction")
    os.makedirs(table_dir, exist_ok=True)

    path = os.path.join(table_dir, "cci_interactions.csv")
    lr_res.to_csv(path, index=False)
    log.info("Exported: %s (%d rows)", path, len(lr_res))

    path = os.path.join(table_dir, "cci_top_interactions.csv")
    top_df.to_csv(path, index=False)
    log.info("Exported: %s (%d rows)", path, len(top_df))


# ═══════════════════════════════════════════════════════════════════════
#  Plots
# ═══════════════════════════════════════════════════════════════════════

def plot_heatmap(top_df, CFG, log):
    """Heatmap of top interaction scores (source→target cell type pairs)."""
    from rna.utils.cell_interaction import format_cci_results

    fig_dir = os.path.join(CFG.figure_dir, "12_cell_interaction")
    os.makedirs(fig_dir, exist_ok=True)

    # Pivot to (source x target) matrix using magnitude_rank or similar
    if "interaction" not in top_df.columns:
        top_df = top_df.copy()
        top_df["interaction"] = (
            top_df["source"].astype(str) + "→" + top_df["target"].astype(str)
        )

    # Count interactions per source→target pair
    st_counts = (
        top_df.groupby(["source", "target"]).size()
        .reset_index(name="n_interactions")
    )
    pivot = st_counts.pivot(index="source", columns="target", values="n_interactions")
    pivot = pivot.fillna(0)

    if pivot.empty or pivot.shape[0] < 2:
        log.warning("Not enough source→target pairs for heatmap — skipping")
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
    ax.set_xlabel("Target cell type")
    ax.set_ylabel("Source cell type")
    ax.set_title(f"Top {len(top_df)} CCI Interactions (LIANA+ rank_aggregate)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("N interactions")

    fig.tight_layout()
    path = os.path.join(fig_dir, "cci_heatmap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", path)


def plot_dotplot(top_df, CFG, log):
    """Dotplot of top ligand-receptor pairs across source→target pairs."""
    fig_dir = os.path.join(CFG.figure_dir, "12_cell_interaction")
    os.makedirs(fig_dir, exist_ok=True)

    if "interaction" not in top_df.columns:
        top_df = top_df.copy()
        if "ligand" in top_df.columns and "receptor" in top_df.columns:
            top_df["interaction"] = (
                top_df["ligand"].astype(str) + "_" + top_df["receptor"].astype(str)
            )
        else:
            log.warning("No interaction label column — skipping dotplot")
            return

    # Build a pivot: rows = LR pairs, columns = source→target pairs
    top_df = top_df.copy()
    top_df["source_target"] = (
        top_df["source"].astype(str) + "→" + top_df["target"].astype(str)
    )

    # Determine value column for color
    score_col = None
    for candidate in ["magnitude_rank", "specificity_rank", "lrscore"]:
        if candidate in top_df.columns:
            score_col = candidate
            break

    if score_col is None:
        log.warning("No scoring column found — skipping dotplot")
        return

    pivot = top_df.pivot_table(
        index="interaction", columns="source_target",
        values=score_col, aggfunc="mean",
    )
    pivot = pivot.fillna(0)

    if pivot.empty:
        log.warning("Empty pivot for dotplot — skipping")
        return

    n_rows, n_cols = pivot.shape
    fig_w = max(8, n_cols * 0.5 + 4.0)
    fig_h = max(4, n_rows * 0.3 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Color-wash: normalize scores for consistent coloring
    vmin, vmax = pivot.values.min(), pivot.values.max()
    if vmin == vmax:
        vmin, vmax = vmin - 0.1, vmax + 0.1

    im = ax.imshow(pivot.values, aspect="auto", cmap="RdYlBu_r",
                   vmin=vmin, vmax=vmax)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_xlabel("Source → Target cell type")
    ax.set_ylabel("Ligand_Receptor pair")
    ax.set_title(f"Top {len(top_df)} CCI Interactions ({score_col})")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(score_col)

    fig.tight_layout()
    path = os.path.join(fig_dir, "cci_dotplot.png")
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
    log = setup_logger("12_cci", os.path.join(CFG.log_dir, "12_cci.log"))
    log.info("Step 12: Cell-Cell Interaction (CCI) analysis via LIANA+")

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

    # Resolve grouping column
    group_col = "cell_type" if "cell_type" in adata.obs.columns else "leiden"
    log.info("Grouping by: %s (%d groups)", group_col, adata.obs[group_col].nunique())

    use_raw = adata.raw is not None
    if use_raw:
        log.info("Using adata.raw for expression (raw counts)")
    else:
        log.warning("adata.raw is not available — using adata.X; "
                     "ensure it contains raw/normalized counts suitable for LIANA")

    n_jobs = getattr(CFG, "n_jobs", 1) or 1
    if n_jobs == 0:
        n_jobs = os.cpu_count() or 4
        log.info("n_jobs=0 → auto-detected %d cores", n_jobs)

    # ── Run CCI permutation testing ─────────────────────────────────────
    from rna.utils.cell_interaction import (
        ensure_gene_symbols,
        run_cci_permutation,
        format_cci_results,
    )

    # ── Ensure gene symbols (LIANA uses HGNC symbols) ───────────────────
    adata = ensure_gene_symbols(adata, log=log)

    lr_res = run_cci_permutation(
        adata,
        groupby_col=group_col,
        resource_name=CFG.cci_lr_database,
        n_perms=CFG.cci_permutations,
        use_raw=use_raw,
        n_jobs=n_jobs,
        log=log,
    )

    # ── Format & export ─────────────────────────────────────────────────
    top_df = format_cci_results(
        lr_res,
        n_top=CFG.cci_n_top_interactions,
        log=log,
    )

    export_results(lr_res, top_df, CFG, log)

    # ── Plot ────────────────────────────────────────────────────────────
    plot_heatmap(top_df, CFG, log)
    plot_dotplot(top_df, CFG, log)

    log.info("Step 12 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
