#!/usr/bin/env python3
"""
Step 06: DE + SVG + nhood enrichment + co-occurrence
========================================================================
  1. Per-cluster DE (Wilcoxon rank-sum) — marker_genes_per_group.csv
  2. Spatially variable genes via Moran's I (sq.gr.spatial_autocorr)
  3. Export SVG rankings
  4. Nhood enrichment (sq.gr.nhood_enrichment) — heatmap + CSV
  5. Co-occurrence (sq.gr.co_occurrence) — distance plot
  6. Interaction matrix (sq.gr.interaction_matrix) — heatmap

Input:  05_annotated.h5ad
Output: marker_genes_per_group.csv, svg_rankings.csv, nhood_enrichment_zscore.csv, co_occurrence_plot.png
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run_de_per_cluster(adata, CFG, log):
    """Wilcoxon rank-sum per cluster (annotated cell types)."""
    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    log.info("Differential expression (groupby='%s')...", group_col)

    sc.tl.rank_genes_groups(
        adata, groupby=group_col,
        method=CFG.de.method,
        n_genes=CFG.de.n_genes,
    )

    # Collect results
    rows = []
    for group in sorted(adata.obs[group_col].unique(), key=str):
        try:
            df = sc.get.rank_genes_groups_df(adata, group=str(group))
            df['group'] = group
            rows.append(df)
        except (KeyError, ValueError) as e:
            log.warning("  Group '%s': %s", group, e)

    if not rows:
        log.error("No DE results produced")
        return None

    marker_df = pd.concat(rows, ignore_index=True)
    csv_path = os.path.join(CFG.table_dir, 'marker_genes_per_group.csv')
    marker_df.to_csv(csv_path, index=False)
    log.info("DE results: %d rows saved → %s", len(marker_df), csv_path)
    return marker_df


def run_spatial_autocorr(adata, CFG, log):
    """Compute Moran's I spatial autocorrelation for SVGs."""
    log.info("Spatial autocorrelation (Moran's I)...")

    # Ensure spatial connectivity graph exists
    if 'spatial_connectivities' not in adata.obsp:
        log.warning("spatial_connectivities missing — rebuilding spatial graph")
        if CFG.spatial.neighbors_radius > 0:
            sq.gr.spatial_neighbors(
                adata, radius=CFG.spatial.neighbors_radius,
                coord_type='generic',
            )
        else:
            sq.gr.spatial_neighbors(
                adata, n_neighs=CFG.spatial.neighbors_n,
                coord_type='generic',
            )

    if 'spatial_connectivities' not in adata.obsp:
        log.error("Cannot compute SVGs — no spatial graph available")
        return None

    try:
        sq.gr.spatial_autocorr(
            adata,
            mode='moran',
            n_perms=100,
            n_jobs=CFG.execution.n_jobs,
        )
    except Exception as e:
        log.warning("Spatial autocorrelation failed: %s", e)
        return None

    # Extract Moran's I results
    moran_key = 'moranI'
    if moran_key in adata.uns:
        moran_df = adata.uns[moran_key].copy()
        if isinstance(moran_df, pd.DataFrame):
            moran_df = moran_df.sort_values('I', ascending=False)
            svg_csv = os.path.join(CFG.table_dir, 'svg_rankings.csv')
            moran_df.to_csv(svg_csv)
            log.info("SVG rankings saved: %s (%d genes)", svg_csv, len(moran_df))

            # Top SVGs
            n_sig = (moran_df['pval_sim'] < 0.05).sum() if 'pval_sim' in moran_df.columns else 0
            log.info("  Moran's I: %d genes, %d significant (p<0.05)", len(moran_df), n_sig)
            if n_sig > 0:
                top_svg = moran_df.head(5).index.tolist()
                log.info("  Top 5 SVGs: %s", top_svg)
            return moran_df
    else:
        log.warning("No 'moranI' key in adata.uns — spatial_autocorr may not have run")

    return None


def plot_top_svg(adata, moran_df, CFG, log):
    """Spatial plot of top spatially variable genes."""
    top_n = min(6, len(moran_df) if moran_df is not None else 0)
    if top_n == 0 or 'spatial' not in adata.obsm:
        return

    top_genes = [g for g in moran_df.head(top_n).index.tolist() if g in adata.var_names]
    if not top_genes:
        return

    fig_dir = os.path.join(CFG.figure_dir, '06_spatial_de')
    os.makedirs(fig_dir, exist_ok=True)

    try:
        fig, axes = plt.subplots(
            max(1, (len(top_genes) + 2) // 3), min(3, len(top_genes)),
            figsize=(5 * min(3, len(top_genes)), 5 * max(1, (len(top_genes) + 2) // 3)),
            squeeze=False,
        )
        for i, gene in enumerate(top_genes):
            ax = axes[i // 3, i % 3]
            safe_plot(sq.pl.spatial_scatter, adata, color=gene, ax=ax,
                      show=False, title=gene)

        # Hide unused subplots
        for j in range(len(top_genes), axes.size):
            axes[j // 3, j % 3].axis('off')

        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'top_svg_spatial.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info("Top SVG spatial plot saved")
    except Exception as e:
        log.warning("SVG spatial plot failed: %s", e)


def run_nhood_enrichment(adata, CFG, log):
    """Nhood enrichment analysis — cluster co-localization z-scores."""
    if 'spatial_connectivities' not in adata.obsp:
        log.warning("spatial_connectivities missing — skipping nhood enrichment")
        return None

    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    log.info("Nhood enrichment (groupby='%s')...", group_col)
    try:
        sq.gr.nhood_enrichment(adata, cluster_key=group_col)
        sq.pl.nhood_enrichment(adata, cluster_key=group_col, show=False)
        fig_path = os.path.join(CFG.figure_dir, '06_spatial_de', 'nhood_enrichment_heatmap.png')
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()

        enrich_key = 'nhood_enrichment'
        if enrich_key in adata.uns:
            enrich_df = pd.DataFrame(adata.uns[enrich_key])
            csv_path = os.path.join(CFG.table_dir, 'nhood_enrichment_zscore.csv')
            enrich_df.to_csv(csv_path)
            log.info("Nhood enrichment z-scores saved: %s", csv_path)
            return enrich_df
    except Exception as e:
        log.warning("Nhood enrichment failed: %s", e)
    return None


def run_co_occurrence(adata, CFG, log):
    """Co-occurrence analysis — cluster spatial co-occurrence vs distance."""
    if 'spatial_connectivities' not in adata.obsp:
        log.warning("spatial_connectivities missing — skipping co-occurrence")
        return

    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    log.info("Co-occurrence analysis (groupby='%s')...", group_col)
    try:
        sq.gr.co_occurrence(adata, cluster_key=group_col)
        sq.pl.co_occurrence(adata, cluster_key=group_col, show=False)
        fig_path = os.path.join(CFG.figure_dir, '06_spatial_de', 'co_occurrence_plot.png')
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        log.info("Co-occurrence plot saved: %s", fig_path)
    except Exception as e:
        log.warning("Co-occurrence failed: %s", e)


def run_interaction_matrix(adata, CFG, log):
    """Interaction matrix — cluster pairwise interaction counts."""
    if 'spatial_connectivities' not in adata.obsp:
        log.warning("spatial_connectivities missing — skipping interaction matrix")
        return

    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    log.info("Interaction matrix (groupby='%s')...", group_col)
    try:
        sq.gr.interaction_matrix(adata, cluster_key=group_col)
        sq.pl.interaction_matrix(adata, cluster_key=group_col, show=False)
        fig_path = os.path.join(CFG.figure_dir, '06_spatial_de', 'interaction_matrix_heatmap.png')
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        log.info("Interaction matrix heatmap saved: %s", fig_path)
    except Exception as e:
        log.warning("Interaction matrix failed: %s", e)



def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("06_spatial_de", os.path.join(CFG.log_dir, "06_spatial_de.log"))
    log.info("Step 06: DE + SVG + nhood enrichment + co-occurrence")

    input_path = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_path):
        log.error("Input not found: %s. Run Step 05 first.", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    # ── 1. Per-cluster DE ──
    marker_df = run_de_per_cluster(adata, CFG, log)

    # ── 2. Spatially variable genes ──
    if CFG.spatial.run_autocorr:
        moran_df = run_spatial_autocorr(adata, CFG, log)
        if moran_df is not None:
            plot_top_svg(adata, moran_df, CFG, log)

            # Subset to top SVGs for downstream analysis
            n_top = min(CFG.spatial.svg_n_top, len(moran_df))
            top_svg_genes = moran_df.head(n_top).index.tolist()
            valid_svg = [g for g in top_svg_genes if g in adata.var_names]
            if valid_svg:
                adata.var['spatially_variable'] = adata.var_names.isin(valid_svg)
                log.info("  Marked %d top SVGs in adata.var['spatially_variable']",
                         adata.var['spatially_variable'].sum())

    # ── 3. Nhood enrichment ──
    run_nhood_enrichment(adata, CFG, log)

    # ── 4. Co-occurrence ──
    run_co_occurrence(adata, CFG, log)

    # ── 5. Interaction matrix ──
    run_interaction_matrix(adata, CFG, log)
    # ── Plot top DE marker per cell type ──
    if marker_df is not None and 'cell_type' in adata.obs:
        group_col = 'cell_type'
        try:
            top_per_group = (
                marker_df[marker_df['logfoldchanges'] > 0.5]
                .sort_values('pvals_adj')
                .groupby('group')
                .head(3)
            )
            top_genes = [g for g in top_per_group['names'].unique().tolist()
                         if g in adata.var_names][:9]
            if top_genes:
                fig_dir = os.path.join(CFG.figure_dir, '06_spatial_de')
                os.makedirs(fig_dir, exist_ok=True)
                fig, axes = plt.subplots(
                    max(1, (len(top_genes) + 2) // 3), min(3, len(top_genes)),
                    figsize=(5 * min(3, len(top_genes)), 5 * max(1, (len(top_genes) + 2) // 3)),
                    squeeze=False,
                )
                for i, gene in enumerate(top_genes):
                    ax = axes[i // 3, i % 3]
                    try:
                        safe_plot(sc.pl.umap, adata, color=gene, ax=ax,
                                  show=False, title=gene, use_raw=True)
                    except Exception as e:
                        log.warning("Top marker UMAP annotation failed: %s", e)
                        ax.text(0.5, 0.5, 'Error', ha='center', va='center')
                for j in range(len(top_genes), axes.size):
                    axes[j // 3, j % 3].axis('off')
                fig.tight_layout()
                fig.savefig(os.path.join(fig_dir, 'top_markers_umap.png'),
                            dpi=150, bbox_inches='tight')
                plt.close(fig)
                log.info("Top marker UMAP plot saved")
        except Exception as e:
            log.warning("Top marker plot failed: %s", e)

    # ── Save h5ad (only if we modified it with SVG info) ──
    if CFG.spatial.run_autocorr:
        # Don't overwrite the annotated h5ad; save SVG stuff separately
        svg_h5ad = os.path.join(CFG.h5ad_dir, "06_svg.h5ad")
        safe_write(adata, svg_h5ad, cfg=CFG)
        log.info("SVG h5ad saved: %s", svg_h5ad)

    log.info("Step 06 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
