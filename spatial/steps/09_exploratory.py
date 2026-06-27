#!/usr/bin/env python3
"""
Step 09: Exploratory analysis for spatial transcriptomics
=============================================================
  1. Spatial scatter plots (cell types on tissue)
  2. Gene expression spatial maps (top markers + SVGs)
  3. Spot composition statistics
  4. Interactive spatial viewer (napari, if available)

Input:  05_annotated.h5ad + 06_svg.h5ad (if available)
Output: Figures + CSV tables
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_plot
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt


def spatial_cell_type_plot(adata, CFG, log):
    """Plot cell types on spatial coordinates."""
    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    if group_col not in adata.obs:
        log.warning("No '%s' column — skipping spatial cell type plot", group_col)
        return

    fig_dir = os.path.join(CFG.figure_dir, '09_exploratory')
    os.makedirs(fig_dir, exist_ok=True)

    try:
        safe_plot(sq.pl.spatial_scatter,
                  adata, color=group_col,
                  shape=None, size=1.5, show=False,
                  save='_09_spatial_celltype.png')
        plt.savefig(os.path.join(fig_dir, 'spatial_cell_type.png'),
                    dpi=200, bbox_inches='tight')
        plt.close()
        log.info("Spatial cell type plot saved")
    except Exception as e:
        log.warning("Spatial cell type plot failed: %s", e)


def spatial_gene_plots(adata, CFG, log):
    """Plot top marker genes on spatial coordinates."""
    fig_dir = os.path.join(CFG.figure_dir, '09_exploratory')
    os.makedirs(fig_dir, exist_ok=True)

    # Priority: SVGs > DE top markers > configured markers
    gene_candidates = []

    if 'spatially_variable' in adata.var:
        svg_genes = adata.var_names[adata.var['spatially_variable']].tolist()
        gene_candidates.extend(svg_genes[:4])
        log.info("  SVG candidates: %d", len(svg_genes))

    # Supplement with top DE markers from marker CSV if available
    marker_csv = os.path.join(CFG.table_dir, 'marker_genes_per_group.csv')
    if os.path.exists(marker_csv) and len(gene_candidates) < 8:
        marker_df = pd.read_csv(marker_csv)
        de_top = (marker_df
                  .sort_values('pvals_adj')
                  .groupby('group')
                  .head(2)['names']
                  .unique()
                  .tolist())
        gene_candidates.extend([g for g in de_top if g not in gene_candidates])
        gene_candidates = gene_candidates[:8]

    # Fallback to configured markers
    if not gene_candidates and CFG.marker_dict:
        for genes in CFG.marker_dict.values():
            gene_candidates.extend(genes[:2])
        gene_candidates = list(dict.fromkeys(gene_candidates))[:8]

    # Filter to available genes
    genes = [g for g in gene_candidates if g in adata.var_names][:9]
    if not genes:
        log.warning("No marker genes available for spatial plot")
        return

    log.info("Plotting %d genes on spatial coordinates: %s", len(genes), genes[:5])
    try:
        fig, axes = plt.subplots(
            max(1, (len(genes) + 2) // 3), min(3, len(genes)),
            figsize=(5 * min(3, len(genes)), 5 * max(1, (len(genes) + 2) // 3)),
            squeeze=False,
        )
        for i, gene in enumerate(genes):
            ax = axes[i // 3, i % 3]
            try:
                safe_plot(sq.pl.spatial_scatter, adata, color=gene,
                          shape=None, size=1.5, ax=ax, show=False)
            except Exception:
                ax.text(0.5, 0.5, gene, ha='center', va='center')

        for j in range(len(genes), axes.size):
            axes[j // 3, j % 3].axis('off')

        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'spatial_marker_genes.png'),
                    dpi=200, bbox_inches='tight')
        plt.close(fig)
        log.info("Spatial gene expression plot saved")
    except Exception as e:
        log.warning("Spatial gene plot failed: %s", e)


def composition_stats(adata, CFG, log):
    """Cluster/cell type composition statistics."""
    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'

    sizes = adata.obs[group_col].value_counts().sort_index()
    log.info("  %s size distribution:", group_col)
    for label, cnt in sizes.items():
        log.info("    %s: %d spots (%.1f%%)", label, cnt, 100 * cnt / adata.n_obs)

    sizes.to_csv(os.path.join(CFG.table_dir, f'{group_col}_sizes.csv'),
                 header=['n_spots'])

    # In-tissue distribution if available
    if 'in_tissue' in adata.obs:
        tissue_counts = adata.obs.groupby(group_col)['in_tissue'].sum()
        log.info("  Tissue spot counts per %s:", group_col)
        for label, cnt in tissue_counts.items():
            log.info("    %s: %d tissue spots", label, cnt)


def spatial_neighbors_summary(adata, CFG, log):
    """Summary statistics about spatial neighborhood graph."""
    if 'spatial_connectivities' not in adata.obsp:
        log.info("  No spatial connectivity graph — skip neighborhood stats")
        return

    try:
        connectivities = adata.obsp['spatial_connectivities']
        if scipy.sparse.issparse(connectivities):
            n_edges = connectivities.nnz
            avg_degree = n_edges / max(adata.n_obs, 1)
            log.info("  Spatial graph: %d edges, avg degree %.1f", n_edges, avg_degree)
    except Exception:
        pass


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("09_exploratory", os.path.join(CFG.log_dir, "09_exploratory.log"))
    log.info("Step 09: Exploratory spatial analysis")

    # Load the most complete data available
    input_path = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_path):
        log.error("Input not found: %s. Run Steps 00–05 first.", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    # Also check for SVG data
    svg_path = os.path.join(CFG.h5ad_dir, "06_svg.h5ad")
    if os.path.exists(svg_path):
        svg_adata = sc.read(svg_path)
        if 'spatially_variable' in svg_adata.var:
            adata.var['spatially_variable'] = svg_adata.var['spatially_variable']
            log.info("  SVG annotations loaded from %s", svg_path)

    # ── 1. Spatial cell type plots ──
    spatial_cell_type_plot(adata, CFG, log)

    # ── 2. Spatial gene expression plots ──
    spatial_gene_plots(adata, CFG, log)

    # ── 3. Composition statistics ──
    composition_stats(adata, CFG, log)

    # ── 4. Spatial neighborhood summary ──
    spatial_neighbors_summary(adata, CFG, log)

    # ── 5. UMAP summary plots ──
    sc.settings.figdir = os.path.join(CFG.figure_dir, '09_exploratory')
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False

    for col in ['cell_type', 'leiden', 'total_counts', 'n_genes_by_counts']:
        if col in adata.obs:
            try:
                safe_plot(sc.pl.umap, adata, color=col, show=False,
                          save=f'_09_{col}.pdf')
            except Exception:
                pass

    log.info("Step 09 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
