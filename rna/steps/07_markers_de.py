#!/usr/bin/env python3
"""
Step 07: 标记基因 + 差异表达分析
=====================================
三层分析:
  Layer 1: 每组 vs 其他 — Wilcoxon rank-sum (多注释层级)
  Layer 2: 相邻发育阶段配对比较 — per cell type
  Layer 3: 发育时间趋势基因 — Spearman 相关

输入: 05_annotated.h5ad (fallback: 04_clustered.h5ad)
输出: tables/*.csv + figures
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_plot
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata
import scipy.sparse as sp
from joblib import Parallel, delayed
import gc

def layer1_markers(adata, CFG, log, group_col):
    """全细胞类型标记基因 (Wilcoxon, each vs rest) for a given annotation column"""
    log.info("[Layer 1] Marker gene detection: groupby=%s", group_col)
    sc.tl.rank_genes_groups(
        adata, groupby=group_col, method='wilcoxon',
        n_genes=CFG.de_n_genes * 2, use_raw=True, pts=True,
        random_state=CFG.random_seed,
    )
    result = sc.get.rank_genes_groups_df(adata, group=None)
    if CFG.de_pval_cutoff is not None:
        result = result[result['pvals_adj'] < CFG.de_pval_cutoff]
    out_path = os.path.join(CFG.table_dir, f'marker_genes_per_group_{group_col}.csv')
    result.to_csv(out_path, index=False)
    log.info("  Exported: %s (%d rows)", out_path, len(result))

    for group in adata.obs[group_col].cat.categories:
        top5 = result[result['group'] == group].head(5)
        if len(top5) > 0:
            log.info("  %s top5: %s", group, ', '.join(top5['names'].values))
    return result


def _layer2_one_pair(ct, s1, s2, adata, ct_col, CFG, log):
    """Worker for parallel Layer 2 paired DE (one cell type, one stage pair)."""
    ct_mask = adata.obs[ct_col] == ct
    stage_mask = adata.obs['stage'].isin([s1, s2])
    sub = adata[ct_mask & stage_mask].copy()
    if sub.n_obs < 20:
        return None
    min_group_size = sub.obs['stage'].value_counts().min()
    if min_group_size < 10:
        log.warning("  %s %s vs %s: only %d cells in smallest group, skipping", ct, s2, s1, min_group_size)
        return None
    try:
        sc.tl.rank_genes_groups(
            sub, groupby='stage', groups=[s2], reference=s1,
            method='t-test', n_genes=CFG.de_n_genes,
            use_raw=True, random_state=CFG.random_seed,
        )
        de_df = sc.get.rank_genes_groups_df(sub, group=s2)
        if CFG.de_pval_cutoff is not None:
            de_df = de_df[de_df['pvals_adj'] < CFG.de_pval_cutoff].copy()
        de_df['cell_type'] = ct
        de_df['comparison'] = f'{s2}_vs_{s1}'
        result = (f'{ct}_{s2}_vs_{s1}', de_df)
        del sub
        gc.collect()
        return result
    except Exception as e:
        log.debug("  %s %s vs %s failed: %s", ct, s2, s1, e)
        del sub
        gc.collect()
        return None

def layer2_pairwise_de(adata, CFG, log, primary_col=None):
    """相邻发育阶段配对差异表达"""
    if 'stage' not in adata.obs or not CFG.stage_order:
        log.info("[Layer 2] No stage info, skipping.")
        return {}
    if not getattr(CFG, 'de_stage_pairwise', True):
        log.info("[Layer 2] de_stage_pairwise=False, skipping.")
        return {}
    stage_pairs = list(zip(CFG.stage_order[:-1], CFG.stage_order[1:]))
    ct_col = primary_col if primary_col else ('cell_type' if 'cell_type' in adata.obs else 'leiden')
    all_results = {}
    log.info("[Layer 2] Adjacent stage pairwise DE (%d pairs, %d types)...",
             len(stage_pairs), adata.obs[ct_col].nunique())

    tasks = [
        (ct, s1, s2)
        for ct in adata.obs[ct_col].cat.categories
        for s1, s2 in stage_pairs
    ]

    if tasks:
        n_jobs = min(getattr(CFG, 'n_jobs', 4) or os.cpu_count() or 1, len(tasks))
        results = Parallel(n_jobs=n_jobs, prefer='threads', require='sharedmem')(
            delayed(_layer2_one_pair)(ct, s1, s2, adata, ct_col, CFG, log)
            for ct, s1, s2 in tasks
        )
        for r in results:
            if r is not None:
                key, de_df = r
                all_results[key] = de_df

    if all_results:
        combined = pd.concat(all_results.values(), ignore_index=True)
        out_path = os.path.join(CFG.table_dir, 'pairwise_stage_de.csv')
        combined.to_csv(out_path, index=False)
        log.info("  Exported: %s (%d rows)", out_path, len(combined))
    return all_results

def layer3_temporal_trends(adata, CFG, log, primary_col=None):
    """发育时间趋势基因 (Spearman 相关 vs 发育顺序)"""
    if 'stage' not in adata.obs or not CFG.stage_order:
        log.info("[Layer 3] No stage info, skipping.")
        return pd.DataFrame()

    stage_numeric = {s: i for i, s in enumerate(CFG.stage_order)}
    ct_col = primary_col if primary_col else ('cell_type' if 'cell_type' in adata.obs else 'leiden')
    log.info("[Layer 3] Temporal trend analysis (per %s)...", ct_col)
    results = []

    for ct in adata.obs[ct_col].cat.categories:
        ct_mask = adata.obs[ct_col] == ct
        n_ct = ct_mask.sum()
        if n_ct < 50:
            continue
        stages = adata.obs.loc[ct_mask, 'stage']
        # 至少 3 个阶段且每阶段 >= 5 细胞
        valid_stages = [s for s in CFG.stage_order
                        if s in stages.values and (stages == s).sum() >= 5]
        if len(valid_stages) < 3:
            continue

        stage_means = {}
        for s in valid_stages:
            s_mask = (stages == s).values
            s_idx = np.flatnonzero(ct_mask.values)[s_mask]
            sub_X = adata.raw[s_idx].X
            mean_expr = sub_X.mean(axis=0).A1 if sp.issparse(sub_X) else sub_X.mean(axis=0)
            stage_means[s] = mean_expr

        stage_nums = np.array([stage_numeric[s] for s in valid_stages])
        mean_matrix = np.stack([stage_means[s] for s in valid_stages], axis=1)
        gene_names = adata.raw.var_names

        # Vectorized Spearman: rank each gene across stages, then Pearson = Spearman
        ranked_genes = np.apply_along_axis(rankdata, 1, mean_matrix)
        ranked_stages = rankdata(stage_nums)
        combined = np.vstack([ranked_genes, ranked_stages.reshape(1, -1)])
        corr_matrix = np.corrcoef(combined)
        corr = corr_matrix[:-1, -1]
        corr_idx = np.argsort(corr)[::-1]
        n_top = min(20, len(corr))

        for i in range(n_top):
            idx = corr_idx[i]
            results.append({
                'cell_type': ct, 'gene': gene_names[idx],
                'spearman_r': corr[idx], 'direction': 'up',
            })
        for i in range(n_top):
            idx = corr_idx[-1 - i]
            results.append({
                'cell_type': ct, 'gene': gene_names[idx],
                'spearman_r': corr[idx], 'direction': 'down',
            })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        out_path = os.path.join(CFG.table_dir, 'temporal_trend_genes.csv')
        results_df.to_csv(out_path, index=False)
        log.info("  Exported: %s (%d rows)", out_path, len(results_df))
    return results_df

def generate_figures(adata, markers_df, CFG, log, primary_col=None):
    sc.settings.figdir = os.path.join(CFG.figure_dir, '07_markers')
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False
    group_col = primary_col if primary_col else ('cell_type' if 'cell_type' in adata.obs else 'leiden')

    # Heatmap: 每个类型 top5 标记
    top5_per_group = (
        markers_df[markers_df['pvals_adj'] < 0.01]
        .groupby('group', observed=True)
        .apply(lambda x: x.nsmallest(5, 'pvals_adj'))
        .reset_index(drop=True)
    )
    top_genes = top5_per_group['names'].unique().tolist()
    top_genes = [g for g in top_genes if g in adata.raw.var_names][:30]
    if len(top_genes) >= 5:
        safe_plot(sc.pl.heatmap, adata, var_names=top_genes,
                  groupby=group_col, show=False, save='_07_marker_heatmap.pdf')

    # 关键标记基因 dotplot
    if CFG.marker_dict:
        all_markers = []
        for genes in CFG.marker_dict.values():
            all_markers.extend([g for g in genes if g in adata.raw.var_names][:2])
        all_markers = list(dict.fromkeys(all_markers))
        if all_markers:
            safe_plot(sc.pl.dotplot, adata, var_names=all_markers,
                      groupby=group_col, show=False, save='_07_dotplot.pdf')

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("07_de", os.path.join(CFG.log_dir, "07_markers_de.log"))
    log.info("Step 07: Marker genes + differential expression analysis")

    # 优先加载 05_annotated.h5ad，回退到 cluster_h5ad
    input_h5ad = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_h5ad):
        input_h5ad = CFG.cluster_h5ad
        log.warning("05_annotated.h5ad not found, falling back to: %s", input_h5ad)
    adata = sc.read(input_h5ad)
    log.info("Loaded: %s — %d cells", input_h5ad, adata.n_obs)

    # Quality awareness (v3.1.0+): check marker_validation PASS rate
    _pass_rate = None
    if 'marker_validation' in adata.obs and adata.n_obs > 0:
        _pass_cells = (adata.obs['marker_validation'] == 'PASS').sum()
        _pass_rate = _pass_cells / adata.n_obs
        log.info("marker_validation PASS rate: %.1f%%", _pass_rate * 100)
        _pass_rate_min = getattr(CFG, 'marker_validation_pass_rate_min', 0.1)
        if _pass_rate < _pass_rate_min:
            log.warning(
                "⚠  marker_validation PASS rate %.1f%% (<%.0f%%) — "
                "DE genes are computed on potentially unreliable cell_type "
                "labels. Results should be interpreted with caution.",
                _pass_rate * 100, _pass_rate_min * 100,
            )
    else:
        log.info("No marker_validation column — skipping quality check.")

    # 自动检测注释层级列
    annotation_cols = []
    for col in ['cell_type_sub', 'cell_type', 'leiden']:
        if col in adata.obs:
            annotation_cols.append(col)
    if not annotation_cols:
        log.error("No annotation columns found in adata.obs")
        sys.exit(1)
    log.info("Detected annotation columns: %s", annotation_cols)
    primary_col = annotation_cols[0]
    log.info("Primary annotation column: %s", primary_col)

    # Layer 1: 遍历所有注释层级进行标记基因检测 (parallel across annotation columns)
    all_markers = {}
    if len(annotation_cols) > 1:
        n_jobs = min(getattr(CFG, 'n_jobs', 4) or os.cpu_count() or 1, len(annotation_cols))
        log.info("Layer 1: parallel marker detection across %d annotation cols (n_jobs=%d)",
                 len(annotation_cols), n_jobs)
        parallel_layer1 = Parallel(n_jobs=n_jobs, prefer='threads')(
            delayed(layer1_markers)(adata.copy(), CFG, log, col)
            for col in annotation_cols
        )
        for col, result_df in zip(annotation_cols, parallel_layer1):
            all_markers[col] = result_df
    else:
        col = annotation_cols[0]
        all_markers[col] = layer1_markers(adata, CFG, log, group_col=col)

    # 导出兼容文件 (使用主注释列)
    combined_path = os.path.join(CFG.table_dir, 'marker_genes_per_group.csv')
    all_markers[primary_col].to_csv(combined_path, index=False)
    log.info("  Exported (compat): %s (%d rows)", combined_path, len(all_markers[primary_col]))

    # Layer 2 & 3: 使用主注释列
    layer2_pairwise_de(adata, CFG, log, primary_col=primary_col)
    layer3_temporal_trends(adata, CFG, log, primary_col=primary_col)
    generate_figures(adata, all_markers[primary_col], CFG, log, primary_col=primary_col)

    log.info("Step 07 complete, took %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
