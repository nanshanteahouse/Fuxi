#!/usr/bin/env python3
"""
Step 08: 轨迹分析 — PAGA + 扩散伪时间 + 分支分析
===================================================
继承 GSE169109 的深度轨迹分析:
  1. PAGA (在子聚类或聚类级别上)
  2. 根细胞自动识别 (ROI 类型或标记基因)
  3. 扩散伪时间 (DPT)
  4. 分支间差异表达
  5. 基因沿伪时间表达趋势

输入: 04_clustered.h5ad (需要 Stage 05 注释结果)
输出: 05_final.h5ad (含 PAGA, DPT, 分支结果) + tables + figures
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

def recompute_neighbors(adata, CFG, log):
    """确保邻居图和 UMAP 存在"""
    if 'neighbors' not in adata.uns:
        log.info("Recomputing neighbors...")
        use_rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
        sc.pp.neighbors(adata, n_pcs=CFG.n_pcs_use,
                        n_neighbors=CFG.n_neighbors,
                        use_rep=use_rep, random_state=CFG.random_seed)
    if 'X_umap' not in adata.obsm:
        log.info("Recomputing UMAP...")
        sc.tl.umap(adata, random_state=CFG.random_seed)
    log.info("Neighbors + UMAP ready")

def run_paga(adata, CFG, log):
    """PAGA 轨迹拓扑"""
    group_col = 'cell_type_sub' if 'cell_type_sub' in adata.obs else \
                ('cell_type' if 'cell_type' in adata.obs else 'leiden')
    log.info("PAGA (groupby=%s)...", group_col)
    sc.tl.paga(adata, groups=group_col)
    n_edges = np.sum(adata.uns['paga']['connectivities'].data > 0)
    log.info("  PAGA edges: %d", n_edges)
    safe_plot(sc.pl.paga, adata, color=group_col, show=False,
              save='_08_paga_graph.pdf', title='PAGA trajectory')
    safe_plot(sc.pl.paga_compare, adata, basis='umap', color=group_col,
              show=False, save='_08_paga_umap.pdf',
              edge_width_scale=0.5, title='PAGA on UMAP')

def find_root_cells(adata, CFG, log):
    """自动识别根细胞"""
    # 方法 1: 指定根细胞类型
    if CFG.root_cell_types:
        log.info("Root cells: type %s + earliest stage", CFG.root_cell_types)
        if 'stage' in adata.obs and CFG.stage_order:
            root_mask = (
                adata.obs['cell_type'].isin(CFG.root_cell_types)
                & (adata.obs['stage'] == CFG.stage_order[0])
            )
        else:
            root_mask = adata.obs['cell_type'].isin(CFG.root_cell_types)
        if root_mask.sum() > 0:
            log.info("  Root cells: %d", root_mask.sum())
            return root_mask.values
        else:
            log.warning("  Root cells of specified type not found, trying marker gene method")

    # 方法 2: 标记基因自动检测
    if CFG.root_markers:
        log.info("Root cells: marker gene method %s", CFG.root_markers)
        markers_present = [g for g in CFG.root_markers if g in adata.raw.var_names]
        if markers_present:
            group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
            cluster_scores = []
            for cl in adata.obs[group_col].cat.categories:
                mask = adata.obs[group_col] == cl
                sub = adata.raw[mask]
                gene_indices = [list(adata.raw.var_names).index(g) for g in markers_present]
                gene_exprs = sub.X[:, gene_indices]
                if issparse(gene_exprs):
                    scores = gene_exprs.mean(axis=0).A1.tolist()
                else:
                    scores = gene_exprs.mean(axis=0).tolist()
                cluster_scores.append((cl, np.mean(scores)))
            cluster_scores.sort(key=lambda x: -x[1])
            best_cl = cluster_scores[0][0]
            root_mask = adata.obs[group_col] == best_cl
            log.info("  High-score cluster: %s (score=%.4f)", best_cl, cluster_scores[0][1])
            log.info("  Root cells: %d", root_mask.sum())
            return root_mask.values

    # 方法 3: 回退到最早阶段的细胞
    log.warning("  Cannot auto-determine root, using earliest stage cells.")
    if 'stage' in adata.obs and CFG.stage_order:
        root_mask = adata.obs['stage'] == CFG.stage_order[0]
        log.info("  Root cells: %d (earliest stage %s)", root_mask.sum(), CFG.stage_order[0])
        return root_mask.values
    # 最终回退: 第一个细胞
    log.warning("  Final fallback: using first cell as root.")
    root_mask = np.zeros(adata.n_obs, dtype=bool)
    root_mask[0] = True
    return root_mask

def compute_dpt(adata, root_mask, CFG, log):
    """扩散图 + 扩散伪时间"""
    log.info("Diffusion map (n_comps=%d)...", CFG.n_diffmap_comps)
    sc.tl.diffmap(adata, n_comps=CFG.n_diffmap_comps)

    log.info("Diffusion pseudotime...")
    adata.uns['iroot'] = np.flatnonzero(root_mask)[0]
    for nb in [CFG.n_branchings, 1, 0]:
        try:
            sc.tl.dpt(adata, n_branchings=nb)
            log.info("  DPT complete (n_branchings=%d)", nb)
            break
        except ValueError:
            log.warning("DPT n_branchings=%d failed, trying n_branchings=%d", nb, 1 if nb > 1 else 0)
            continue
    log.info("  DPT range: %.3f – %.3f",
             adata.obs['dpt_pseudotime'].min(),
             adata.obs['dpt_pseudotime'].max())

    safe_plot(sc.pl.umap, adata, color='dpt_pseudotime', show=False,
              save='_08_pseudotime.pdf', cmap='plasma')
    safe_plot(sc.pl.diffmap, adata, color='dpt_pseudotime', show=False,
              save='_08_pseudotime_diffmap.pdf', cmap='plasma')

def branch_analysis(adata, CFG, log):
    """分支间差异表达 (GSE169109 策略)"""
    if 'cell_type' not in adata.obs:
        log.info("No cell_type annotation, skipping branch analysis.")
        return

    if hasattr(CFG, 'trajectory_branches') and CFG.trajectory_branches:
        branches = CFG.trajectory_branches
    else:
        # Auto-detect: use cell type pairs from PAGA graph
        if 'cell_type' in adata.obs:
            avail_types = list(adata.obs['cell_type'].cat.categories)
            branches = []
            for i in range(len(avail_types) - 1):
                branches.append((avail_types[i], avail_types[i + 1]))
        else:
            branches = []
    # 仅保留数据中存在的分支
    avail_types = set(adata.obs['cell_type'].cat.categories)
    branches = [(p, c) for p, c in branches if p in avail_types and c in avail_types]

    if not branches:
        log.info("No matching branch pairs, skipping.")
        return

    log.info("Branch differential expression analysis...")
    branch_results = []
    for parent, child in branches:
        mask = adata.obs['cell_type'].isin([parent, child])
        sub = adata[mask].copy()
        if sub.obs['cell_type'].value_counts().min() < 10:
            log.info("  %s → %s: insufficient cells", parent, child)
            continue
        try:
            sc.tl.rank_genes_groups(
                sub, groupby='cell_type', groups=[child], reference=parent,
                method='wilcoxon', n_genes=50, use_raw=True,
                random_state=CFG.random_seed,
            )
            de_df = sc.get.rank_genes_groups_df(sub, group=child)
            if CFG.de_pval_cutoff is not None:
                de_df = de_df[de_df['pvals_adj'] < CFG.de_pval_cutoff].copy()
            de_df['branch'] = f'{child}_vs_{parent}'
            branch_results.append(de_df)
            n_up = (de_df['logfoldchanges'] > 0).sum()
            n_down = (de_df['logfoldchanges'] < 0).sum()
            log.info("  %s → %s: %d DEGs (%d up, %d down)",
                     parent, child, len(de_df), n_up, n_down)
        except Exception as e:
            log.debug("  %s → %s failed: %s", parent, child, e)

    if branch_results:
        combined = pd.concat(branch_results, ignore_index=True)
        out_path = os.path.join(CFG.table_dir, 'branch_deg.csv')
        combined.to_csv(out_path, index=False)
        log.info("  Branch DEG exported: %s (%d rows)", out_path, len(combined))

def gene_trends(adata, CFG, log):
    """基因表达沿伪时间趋势"""
    if 'dpt_pseudotime' not in adata.obs:
        log.info("No DPT, skipping gene trends.")
        return

    dev_genes = ['SOX2', 'HES1', 'MKI67', 'PAX6', 'DCX', 'STMN2',
                 'RBFOX3', 'NEUROD1', 'GFAP', 'PDGFRA', 'MBP']
    dev_genes = [g for g in dev_genes if g in adata.raw.var_names]
    log.info("Gene trends along pseudotime (%d genes)...", len(dev_genes))

    # 前 6 个: 散点图
    for gene in dev_genes[:6]:
        safe_plot(sc.pl.scatter, adata, x='dpt_pseudotime', y=gene,
                  use_raw=True, show=False, save=f'_08_trend_{gene}.pdf')

    # 热图
    if len(dev_genes) >= 5:
        safe_plot(sc.pl.heatmap,
                  adata[adata.obs['dpt_pseudotime'].notna()].copy(),
                  var_names=dev_genes, groupby='dpt_pseudotime',
                  show=False, save='_08_dev_genes_heatmap.pdf')

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("08_trajectory", os.path.join(CFG.log_dir, "08_trajectory.log"))
    log.info("Step 08: Trajectory analysis")

    input_path = CFG.annotated_h5ad if os.path.exists(CFG.annotated_h5ad) else CFG.cluster_h5ad
    adata = sc.read(input_path)
    log.info("Loaded: %s — %d cells", input_path, adata.n_obs)

    # 设置图输出目录（必须在 plot 调用之前，否则 scanpy save= 默认写到 ./figures/）
    sc.settings.figdir = CFG.figure_dir

    recompute_neighbors(adata, CFG, log)
    run_paga(adata, CFG, log)
    root_mask = find_root_cells(adata, CFG, log)
    compute_dpt(adata, root_mask, CFG, log)
    branch_analysis(adata, CFG, log)
    gene_trends(adata, CFG, log)

    # 最终可视化 (figdir 已在上面设置)
    for color in ['stage', 'cell_type', 'cell_type_sub', 'dpt_pseudotime']:
        if color in adata.obs or color in adata.obsm:
            safe_plot(sc.pl.umap, adata, color=color, show=False,
                      save=f'_08_final_{color}.pdf')

    safe_write(adata, CFG.final_h5ad, cfg=CFG)
    log.info("Step 08 complete, took %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
