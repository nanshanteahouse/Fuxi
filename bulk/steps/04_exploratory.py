#!/usr/bin/env python3
"""
Step 04: Exploratory analysis and visualization for bulk RNA-seq.

Generates 4 figures:
  - PCA plot (samples colored by contrast group)
  - Sample distance heatmap (Euclidean distance on normalized counts)
  - Top DEG heatmap (Z-score normalized, top 50 significant genes)
  - Top gene boxplots (expression of top 10 DEGs per sample group)

Input:  02_de.h5ad (normalized counts) from CFG.h5ad_dir
        tables/02_de_significant.csv
Output: figures/04_pca.png
        figures/04_sample_heatmap.png
        figures/04_de_heatmap.png
        figures/04_top_genes_boxplot.png
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("04_exploratory", os.path.join(CFG.log_dir, "04_exploratory.log"))
    log.info("Step 04: Exploratory analysis / visualization")

    # ── Resume check ──────────────────────────────────────────────────
    pca_fig = os.path.join(CFG.figure_dir, "04_pca.png")
    if os.path.exists(pca_fig):
        log.info("Skip: %s already exists. Delete figures to force rerun.", pca_fig)
        return

    # ── Ensure figure dir ─────────────────────────────────────────────
    os.makedirs(CFG.figure_dir, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────
    de_h5ad = os.path.join(CFG.h5ad_dir, "02_de.h5ad")
    log.info("Loading %s", de_h5ad)
    adata = sc.read(de_h5ad)
    log.info("Loaded: %d samples x %d genes", adata.n_obs, adata.n_vars)

    contrast_col = getattr(CFG.bulk, "contrast_column", "condition")
    log.info("Contrast column: %s", contrast_col)

    # ── 1. PCA Plot ───────────────────────────────────────────────────
    _pca_plot(adata, contrast_col, CFG.figure_dir, log)

    # ── 2. Sample distance heatmap ────────────────────────────────────
    _sample_distance_heatmap(adata, contrast_col, CFG.figure_dir, log)

    # ── Load DEG table for DE-specific plots ──────────────────────────
    sig_path = os.path.join(CFG.table_dir, "02_de_significant.csv")
    if os.path.isfile(sig_path):
        sig_df = pd.read_csv(sig_path)
        n_sig = len(sig_df)
        log.info("Loaded %d significant DEGs from %s", n_sig, sig_path)
    else:
        sig_df = pd.DataFrame()
        n_sig = 0
        log.warning("DEG file not found: %s", sig_path)

    # ── 3. Top DEG heatmap ────────────────────────────────────────────
    if n_sig == 0:
        log.warning("No significant DEGs found — skipping DE heatmap")
    else:
        _de_heatmap(adata, sig_df, contrast_col, CFG.figure_dir, log)

    # ── 4. Top gene boxplots ──────────────────────────────────────────
    if n_sig == 0:
        log.warning("No significant DEGs found — skipping top gene boxplots")
    else:
        _top_genes_boxplot(adata, sig_df, contrast_col, CFG.figure_dir, log)

    log.info("Step 04 complete, took %.1fs", time.time() - t0)


# ═══════════════════════════════════════════════════════════════════════
#  Plot helpers
# ═══════════════════════════════════════════════════════════════════════


def _pca_plot(adata, contrast_col, figure_dir, log):
    """PCA plot of samples colored by contrast group."""
    try:
        if adata.n_obs < 2:
            log.warning("Skipping PCA: need at least 2 samples (got %d)", adata.n_obs)
            return

        n_comps = min(50, adata.n_obs - 1)
        log.info("Running PCA with %d components...", n_comps)

        # Use normalized counts from .X
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X_log = np.log1p(X)

        sc.pp.pca(adata, n_comps=n_comps)
        fig, ax = plt.subplots(figsize=(8, 6))

        if contrast_col in adata.obs.columns:
            sc.pl.pca(adata, color=contrast_col, ax=ax, show=False)
        else:
            # Fallback: no contrast column, just plot PC1 vs PC2
            sc.pl.pca(adata, ax=ax, show=False)
            log.warning("contrast_column '%s' not in adata.obs — PCA without color", contrast_col)

        fig.savefig(os.path.join(figure_dir, "04_pca.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("PCA plot saved: 04_pca.png")

    except Exception as e:
        log.warning("PCA plot failed: %s", e)


def _sample_distance_heatmap(adata, contrast_col, figure_dir, log):
    """Euclidean distance heatmap between samples (log1p normalized counts)."""
    try:
        log.info("Computing sample distance matrix...")

        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X_log = np.log1p(X)

        # Euclidean distance between samples
        dist_matrix = squareform(pdist(X_log, metric="euclidean"))
        sample_names = adata.obs_names.tolist()
        dist_df = pd.DataFrame(dist_matrix, index=sample_names, columns=sample_names)

        # Build annotation DataFrame for row colors
        row_colors_df = None
        if contrast_col in adata.obs.columns:
            groups = adata.obs[contrast_col].astype("category")
            palette = sns.color_palette("Set2", n_colors=len(groups.cat.categories))
            color_map = dict(zip(groups.cat.categories, palette))
            # Build color mapping without map() to avoid pandas MultiIndex issue
            cat_to_color = pd.Series(color_map)
            color_vals = cat_to_color[groups.values].values
            row_colors_df = pd.DataFrame({"condition": color_vals}, index=sample_names)

        n = len(sample_names)
        figsize = (max(8, n * 0.8), max(6, n * 0.6))
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            dist_df,
            cmap="viridis",
            xticklabels=sample_names,
            yticklabels=sample_names,
            linewidths=0.5,
            linecolor="gray",
            ax=ax,
        )
        ax.set_title("Sample Distance (Euclidean, log1p counts)")
        fig.tight_layout()
        fig.savefig(os.path.join(figure_dir, "04_sample_heatmap.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("Sample distance heatmap saved: 04_sample_heatmap.png")

    except Exception as e:
        log.warning("Sample distance heatmap failed: %s", e)


def _de_heatmap(adata, sig_df, contrast_col, figure_dir, log):
    """Z-score heatmap of top 50 significant DEGs across samples."""
    try:
        log.info("Generating top DEG heatmap...")

        # Get top 50 genes by padj
        if "padj" in sig_df.columns:
            top_genes = sig_df.nsmallest(50, "padj")["gene"].tolist()
        else:
            top_genes = sig_df.head(50)["gene"].tolist()
        n_top = len(top_genes)

        if n_top == 0:
            log.warning("No DEGs to plot — skipping DE heatmap")
            return

        log.info("Top %d DEGs for heatmap", n_top)

        # Get expression values
        gene_mask = [g in adata.var_names for g in top_genes]
        available = [g for g, m in zip(top_genes, gene_mask) if m]
        if not available:
            log.warning("No top DEGs found in adata.var_names — skipping DE heatmap")
            return

        X = adata[:, available].X
        if hasattr(X, "toarray"):
            X = X.toarray()

        # Z-score normalize (row-wise = per gene)
        X_z = zscore(X, axis=1)

        heatmap_df = pd.DataFrame(
            X_z,
            index=adata.obs_names,
            columns=available,
        ).T

        # Build annotation colors as plain DataFrame
        col_colors_df = None
        if contrast_col in adata.obs.columns:
            groups = adata.obs[contrast_col].astype("category")
            palette = sns.color_palette("Set2", n_colors=len(groups.cat.categories))
            color_map = dict(zip(groups.cat.categories, palette))
            # Build annotation without map() to avoid pandas MultiIndex issue
            cat_to_color = pd.Series(color_map)
            color_vals = cat_to_color[groups.values].values
            col_colors_df = pd.DataFrame({"condition": color_vals}, index=adata.obs_names)

        n_samp = adata.n_obs
        figsize = (max(8, n_samp * 0.5), max(10, n_top * 0.3))
        g = sns.clustermap(
            heatmap_df,
            z_score=None,  # already z-scored
            cmap="RdBu_r",
            center=0,
            vmin=-3, vmax=3,
            figsize=figsize,
            col_colors=col_colors_df,
            linewidths=0.3,
            xticklabels=True,
            yticklabels=True,
            dendrogram_ratio=(0.1, 0.2),
        )
        g.fig.suptitle(f"Top {n_top} DEGs — Z-score (log1p counts)", y=1.02)
        g.savefig(os.path.join(figure_dir, "04_de_heatmap.png"), dpi=150, bbox_inches="tight")
        plt.close(g.fig)
        log.info("DE heatmap saved: 04_de_heatmap.png")

    except Exception as e:
        log.warning("DE heatmap failed: %s", e)


def _top_genes_boxplot(adata, sig_df, contrast_col, figure_dir, log):
    """Expression boxplots for top 10 significant DEGs across sample groups."""
    try:
        log.info("Generating top gene boxplots...")

        # Get top 10 genes by padj
        if "padj" in sig_df.columns:
            top_genes = sig_df.nsmallest(10, "padj")["gene"].tolist()
        else:
            top_genes = sig_df.head(10)["gene"].tolist()
        n_top = len(top_genes)

        if n_top == 0:
            log.warning("No DEGs to plot — skipping boxplots")
            return

        log.info("Top %d DEGs for boxplots", n_top)

        # Filter to genes present in adata
        available = [g for g in top_genes if g in adata.var_names]
        if not available:
            log.warning("No top DEGs found in adata.var_names — skipping boxplots")
            return

        X = adata[:, available].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        expr_df = pd.DataFrame(X, index=adata.obs_names, columns=available)

        # Add group info
        if contrast_col in adata.obs.columns:
            expr_df["group"] = adata.obs[contrast_col].values
        else:
            expr_df["group"] = "all"

        # Melt to long form for seaborn
        plot_df = expr_df.reset_index().melt(
            id_vars=["index", "group"],
            var_name="gene",
            value_name="expression",
        )

        ncols = min(5, n_top)
        nrows = int(np.ceil(n_top / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.5),
                                 squeeze=False)

        for idx, gene in enumerate(available):
            row = idx // ncols
            col = idx % ncols
            ax = axes[row, col]

            gene_df = plot_df[plot_df["gene"] == gene]
            sns.boxplot(data=gene_df, x="group", y="expression", ax=ax, palette="Set2")
            sns.stripplot(data=gene_df, x="group", y="expression", ax=ax,
                          color="black", size=4, alpha=0.6, jitter=True)

            ax.set_title(gene, fontsize=10)
            ax.set_xlabel("")
            ax.set_ylabel("Normalized count")

        # Hide any unused subplots
        for idx in range(len(available), nrows * ncols):
            row = idx // ncols
            col = idx % ncols
            axes[row, col].set_visible(False)

        fig.suptitle("Top DEG Expression by Group", y=1.02, fontsize=12)
        fig.tight_layout()

        fig.savefig(os.path.join(figure_dir, "04_top_genes_boxplot.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("Top gene boxplots saved: 04_top_genes_boxplot.png")

    except Exception as e:
        log.warning("Top gene boxplots failed: %s", e)


if __name__ == "__main__":
    main()
