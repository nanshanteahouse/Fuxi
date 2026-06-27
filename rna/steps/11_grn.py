#!/usr/bin/env python3
"""
Step 11: Gene Regulatory Network (GRN) analysis via decoupler
==============================================================
Pseudobulk aggregation per cell_type → TF activity inference → heatmap + table.

Methods:
  * CollecTRI (default) — literature-based signed TF regulons, ~1,185 TFs
  * DoRothEA — curated TF regulons w/ confidence-level weighting (~429 TFs)

输入: 05_annotated.h5ad (requires cell_type column)
输出:
  {table_dir}/11_grn/tf_activity_per_cell_type.csv   — TF activity scores
  {table_dir}/11_grn/tf_activity_pvals.csv           — associated p-values
  {figure_dir}/11_grn/tf_activity_heatmap.png        — dendrogram + heatmap
  {h5ad_dir}/11_grn.h5ad                             — pseudobulk AnnData + obsm['X_tf_activity']
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram


def _as_dense(adata, use_raw: bool = False):
    """Return a dense (obs x var) expression matrix."""
    src = adata.raw if use_raw else adata
    X = src.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    elif hasattr(X, 'todense'):
        X = np.asarray(X.todense())
    return X


def build_pseudobulk(adata, group_col: str, use_raw: bool = True, log: object = None) -> pd.DataFrame:
    """Per-group mean expression -> (n_groups x n_genes) DataFrame.

    If data are stored as log1p-transformed counts in adata.X and raw counts
    in adata.raw, we work on raw counts and then add a pseudocount + log1p
    before returning.  Otherwise we use .X directly.
    """
    if group_col not in adata.obs:
        if log:
            log.warning("%s not in adata.obs - using 'leiden'", group_col)
        group_col = 'leiden'

    groups = adata.obs[group_col].values
    var_names = adata.raw.var_names if use_raw and adata.raw else adata.var_names

    if log:
        log.info("Pseudobulk: %d cells -> %d groups", adata.n_obs, len(adata.obs[group_col].cat.categories))

    unique_groups = adata.obs[group_col].cat.categories
    n_groups = len(unique_groups)
    n_genes = len(var_names)

    X = _as_dense(adata, use_raw=use_raw)

    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    group_indices = np.array([group_to_idx[g] for g in groups])

    pseudo = np.zeros((n_groups, n_genes), dtype=np.float64)
    for g_idx in range(n_groups):
        mask = group_indices == g_idx
        if mask.any():
            pseudo[g_idx] = X[mask].mean(axis=0)

    pseudo = np.log1p(pseudo)

    df = pd.DataFrame(pseudo, index=unique_groups, columns=var_names)
    if log:
        log.info("  Pseudobulk matrix: %d x %d", n_groups, n_genes)
    return df


def filter_regulon_net(net: pd.DataFrame, min_genes: int = 5, log: object = None) -> pd.DataFrame:
    """Remove regulons with too few target genes present."""
    n_before = net['source'].nunique()
    gene_counts = net.groupby('source')['target'].nunique()
    keep = gene_counts[gene_counts >= min_genes].index
    net_filt = net[net['source'].isin(keep)].copy()
    if log:
        log.info("Regulon filter (>=%d targets): %d -> %d TFs", min_genes, n_before, len(keep))
    return net_filt


def run_grn(pseudo_df: pd.DataFrame, net: pd.DataFrame, log: object) -> tuple:
    """Run ULM enrichment on pseudobulk data.

    Returns (estimates_df, pvals_df, filtered_net) where rows = groups, cols = TFs.
    filtered_net is the net filtered to available genes.
    """
    import decoupler as dc

    log.info("Running ULM enrichment on %d cells x %d genes", pseudo_df.shape[0], pseudo_df.shape[1])

    avail_genes = set(pseudo_df.columns)
    net = net[net['target'].isin(avail_genes)].copy()
    log.info("  Regulon edges covering available genes: %d", len(net))

    estimates, pvals = dc.mt.ulm(pseudo_df, net, verbose=False)

    est_df = pd.DataFrame(estimates, index=pseudo_df.index, columns=estimates.columns)
    pval_df = pd.DataFrame(pvals, index=pseudo_df.index, columns=estimates.columns)

    log.info("  Activity matrix: %d cell types x %d TFs", est_df.shape[0], est_df.shape[1])
    return est_df, pval_df, net


def top_variable_tfs(estimates_df: pd.DataFrame, n_top: int, log: object) -> pd.DataFrame:
    """Select top N TFs by activity variance across cell types."""
    var = estimates_df.var(axis=0)
    top = var.nlargest(n_top).index
    log.info("Top %d TFs by variance: %s", n_top, ', '.join(top[:20].tolist()))
    return estimates_df[top]


def export_results(estimates_df, top_df, pvals_df, net_top, CFG, log):
    """Save tables and checkpoint AnnData."""
    table_dir = os.path.join(CFG.table_dir, "11_grn")
    os.makedirs(table_dir, exist_ok=True)

    path = os.path.join(table_dir, "tf_activity_per_cell_type.csv")
    estimates_df.to_csv(path)
    log.info("Exported: %s", path)

    path = os.path.join(table_dir, "tf_activity_pvals.csv")
    pvals_df.to_csv(path)
    log.info("Exported: %s", path)

    # Export TF-to-target edges for top-variance TFs
    path = os.path.join(table_dir, "tf_target_edges.csv")
    net_top.to_csv(path, index=False)
    log.info("Exported: %s (%d edges)", path, len(net_top))

    # Export per-TF target gene count summary
    target_counts = (
        net_top.groupby('source')['target']
        .nunique()
        .reset_index()
        .rename(columns={'source': 'tf', 'target': 'n_targets'})
        .sort_values('n_targets', ascending=False)
    )
    path = os.path.join(table_dir, "tf_target_counts.csv")
    target_counts.to_csv(path, index=False)
    log.info("Exported: %s (%d TFs)", path, len(target_counts))

    if net_top.empty:
        log.warning("Top-TF edge list is empty — no edges to export")

    adata_pb = sc.AnnData(
        X=np.expm1(top_df.values),
        obs=pd.DataFrame(index=top_df.index),
        var=pd.DataFrame(index=top_df.columns),
    )
    adata_pb.obsm['X_tf_activity'] = estimates_df.values
    adata_pb.uns['tf_activity_cell_types'] = estimates_df.index.tolist()
    adata_pb.uns['tf_activity_tfs'] = estimates_df.columns.tolist()
    safe_write(adata_pb, CFG.grn_h5ad, cfg=CFG)
    log.info("Checkpoint saved: %s", CFG.grn_h5ad)


def plot_heatmap(top_df, CFG, log):
    """Clustered heatmap + row/col dendrograms for top TF activities."""
    fig_dir = os.path.join(CFG.figure_dir, "11_grn")
    os.makedirs(fig_dir, exist_ok=True)

    n_tfs = top_df.shape[1]
    n_cts = top_df.shape[0]
    if n_tfs < 2 or n_cts < 2:
        log.warning("Too few TFs or cell types for heatmap - skipping")
        return

    data = top_df.values.T  # (TFs x cell_types)
    z_rows = linkage(data, method='ward')
    z_cols = linkage(data.T, method='ward')

    # ---------- layout ----------
    longest_tf = max(len(n) for n in top_df.columns)
    tf_label_w = longest_tf * 0.075               # estimated inches for TF labels

    left_pct = max(0.18, min(0.40, longest_tf * 0.008))

    tf_yunit = min(0.30, max(0.18, 6.0 / max(1, n_tfs)))
    tf_hm = max(0.28, tf_yunit * 0.7)
    ct_xunit = 0.45

    top_margin = 0.8
    right_pad = tf_label_w + 0.55                 # TF labels + colorbar gap
    bottom_margin = 0.6

    left_margin = left_pct * max(5, n_cts * ct_xunit)
    heatmap_w = n_cts * ct_xunit
    heatmap_h = n_tfs * tf_hm

    fig_w = left_margin + heatmap_w + right_pad
    fig_h = top_margin + 0.5 + heatmap_h + bottom_margin

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(2, 2,
                          width_ratios=[left_margin, heatmap_w],
                          height_ratios=[0.5, heatmap_h],
                          hspace=0.0, wspace=0.0,
                          left=left_margin / fig_w,
                          right=(left_margin + heatmap_w) / fig_w,
                          top=(top_margin + 0.5 + heatmap_h) / fig_h,
                          bottom=bottom_margin / fig_h)

    # ---------- row dendrogram (left) ----------
    ax_row = fig.add_subplot(gs[1, 0])
    d_rows = dendrogram(z_rows, ax=ax_row, orientation='left',
                        link_color_func=lambda k: '#555555',
                        above_threshold_color='#bbbbbb',
                        no_labels=True)
    row_idx = d_rows['leaves']
    ax_row.invert_xaxis()
    ax_row.set_xticks([])
    ax_row.set_yticks([])
    for s in ax_row.spines.values():
        s.set_visible(False)

    # ---------- col dendrogram (top) ----------
    ax_col = fig.add_subplot(gs[0, 1])
    d_cols = dendrogram(z_cols, ax=ax_col, orientation='top',
                        link_color_func=lambda k: '#555555',
                        above_threshold_color='#bbbbbb',
                        no_labels=True)
    col_idx = d_cols['leaves']
    ax_col.set_yticks([])
    ax_col.set_xticks([])
    for s in ax_col.spines.values():
        s.set_visible(False)

    # ---------- heatmap ----------
    data_clust = data[row_idx, :][:, col_idx]
    tf_labels = top_df.columns[row_idx]
    ct_labels = top_df.index[col_idx]

    ax_hm = fig.add_subplot(gs[1, 1])
    vabs = np.percentile(np.abs(data_clust), 90)
    im = ax_hm.imshow(data_clust, aspect='auto', cmap='RdBu_r',
                       vmin=-vabs, vmax=vabs, interpolation='nearest')

    ax_hm.set_xticks(range(n_cts))
    ax_hm.set_xticklabels(ct_labels, rotation=35, ha='right', fontsize=8.5)
    ax_hm.set_yticks(range(n_tfs))
    ax_hm.set_yticklabels(tf_labels, fontsize=7.0)
    ax_hm.yaxis.tick_right()
    ax_hm.tick_params(length=0, pad=3)

    # ---------- colorbar ----------
    cax = fig.add_axes([
        (left_margin + heatmap_w + tf_label_w + 0.12) / fig_w,
        (bottom_margin + heatmap_h * 0.15) / fig_h,
        0.012,
        (heatmap_h * 0.45) / fig_h,
    ])
    cb = fig.colorbar(im, cax=cax, orientation='vertical')
    cb.set_label('Activity score', fontsize=8)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle(f'TF Activity (ULM) - Top {n_tfs} Regulons',
                 fontsize=11, fontweight='bold',
                 x=(left_margin + heatmap_w / 2) / fig_w,
                 y=(top_margin + 0.5 + heatmap_h + 0.15) / fig_h,
                 ha='center', va='bottom')

    path = os.path.join(fig_dir, "tf_activity_heatmap.png")
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log.info("Heatmap saved: %s", path)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("11_grn", os.path.join(CFG.log_dir, "11_grn.log"))
    log.info("Step 11: GRN regulatory network analysis (decoupler)")

    if not getattr(CFG, 'run_grn', True):
        log.info("run_grn=False - skipping")
        return

    # ---------- Load annotated data ----------
    adata = sc.read(CFG.annotated_h5ad)
    log.info("Loaded: %s - %d cells, %d genes", CFG.annotated_h5ad, adata.n_obs, adata.n_vars)

    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    log.info("Group column: %s (%d categories)", group_col, adata.obs[group_col].nunique())

    # ---------- Pseudobulk ----------
    use_raw = adata.raw is not None
    pseudo_df = build_pseudobulk(adata, group_col, use_raw=use_raw, log=log)

    # ---------- Fetch regulon network ----------
    import decoupler as dc

    species = getattr(CFG, 'grn_species', 'human')
    log.info("Regulon: CollecTRI (%s)", species)
    net = dc.op.collectri(organism=species)
    net = net[net['weight'] > 0].copy()

    min_size = getattr(CFG, 'grn_min_regulon_size', 5)
    net = filter_regulon_net(net, min_genes=min_size, log=log)

    # ---------- Run activity inference ----------
    est_df, pval_df, net_filtered = run_grn(pseudo_df, net, log)

    # ---------- Select top TFs ----------
    n_top = min(getattr(CFG, 'grn_n_top_regulons', 50), est_df.shape[1])
    top_df = top_variable_tfs(est_df, n_top, log)

    # ---------- Filter edge list to top-variance TFs ----------
    top_tfs = set(top_df.columns)
    net_top = net_filtered[net_filtered['source'].isin(top_tfs)].copy()
    log.info("Top-TF edges: %d (from %d total filtered edges)", len(net_top), len(net_filtered))

    # ---------- Export ----------
    export_results(est_df, top_df, pval_df, net_top, CFG, log)
    plot_heatmap(top_df, CFG, log)

    elapsed = time.time() - t0
    log.info("Step 11 complete (took %.1fs).", elapsed)


if __name__ == '__main__':
    main()
