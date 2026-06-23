#!/usr/bin/env python3
"""
Step 02: QC filtering (doublets already removed in Step 01)
============================================================
继承了 GSE169109 的最佳实践:
  1. QC 指标 (mito%, ribo%, 复杂度)
  2. 过滤 predicted_doublet 细胞 (由 Step 01 产生)
  3. 自适应 MAD 或全局阈值过滤

输入: 01_doublet.h5ad (含 doublet_scores, predicted_doublet 列)
输出: 02_qc.h5ad (过滤后的细胞 + QC 指标)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import scanpy as sc

def compute_qc_metrics(adata, cfg, log):
    log.info("Computing QC metrics...")
    mt_mask = adata.var_names.str.startswith(cfg.mt_gene_pattern)
    if cfg.mt_gene_list:
        mt_mask = mt_mask | adata.var_names.isin(cfg.mt_gene_list)
    adata.var['mt'] = mt_mask
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo'],
        percent_top=[20], log1p=True, inplace=True,
    )
    adata.obs['log_genes_per_umi'] = (
        np.log10(adata.obs['n_genes_by_counts'])
        / np.log10(adata.obs['total_counts'])
    ).replace([np.inf, -np.inf], np.nan)
    log.info("  Median genes/cell: %.0f", adata.obs['n_genes_by_counts'].median())
    log.info("  Median UMIs/cell: %.0f", adata.obs['total_counts'].median())
    log.info("  Median mito%%:    %.2f%%", adata.obs['pct_counts_mt'].median())
    log.info("  Median complexity: %.3f", adata.obs['log_genes_per_umi'].median())



def filter_cells(adata, cfg, log):
    n_before = adata.n_obs

    f_doublet = adata.obs['predicted_doublet']
    n_doublet = f_doublet.sum()
    log.info("Doublet filtering: removing %d predicted doublets (%.1f%%)",
             n_doublet, 100 * n_doublet / n_before if n_before else 0)

    log.info("Applying QC filtering...")
    min_g = cfg.min_genes
    max_g = cfg.max_genes
    max_m = cfg.max_pct_mito
    min_cpx = cfg.min_genes_per_umi

    f_genes_low  = adata.obs['n_genes_by_counts'] < min_g
    f_genes_high = adata.obs['n_genes_by_counts'] > max_g
    f_mito       = adata.obs['pct_counts_mt'] > max_m
    f_cpx        = adata.obs['log_genes_per_umi'] < min_cpx
    f_any = f_genes_low | f_genes_high | f_mito | f_cpx

    log.info("  Filtering breakdown:")
    log.info("    n_genes < %d:     %6d (%.1f%%)", min_g, f_genes_low.sum(), 100*f_genes_low.mean())
    log.info("    n_genes > %d:     %6d (%.1f%%)", max_g, f_genes_high.sum(), 100*f_genes_high.mean())
    log.info("    mito > %d%%:      %6d (%.1f%%)", max_m, f_mito.sum(), 100*f_mito.mean())
    log.info("    complexity < %.2f: %6d (%.1f%%)", min_cpx, f_cpx.sum(), 100*f_cpx.mean())
    log.info("    Total (dedup):    %6d (%.1f%%)", f_any.sum(), 100*f_any.mean())

    mask = ~f_doublet & ~f_any
    adata = adata[mask].copy()
    log.info("  After QC filtering: %d cells", adata.n_obs)
    return adata

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("02_qc", os.path.join(CFG.log_dir, "02_qc.log"))
    log.info("Step 02: QC filtering (doublets already removed in Step 01)")

    input_path = os.path.join(CFG.h5ad_dir, "01_doublet.h5ad")
    adata = sc.read(input_path)
    log.info("Loaded: %s — %d cells × %d genes",
             input_path, adata.n_obs, adata.n_vars)

    compute_qc_metrics(adata, CFG, log)
    adata = filter_cells(adata, CFG, log)
    sc.pp.filter_genes(adata, min_cells=CFG.min_cells_per_gene)
    log.info("After gene filtering: %d genes", adata.n_vars)

    safe_write(adata, CFG.qc_h5ad, cfg=CFG)
    log.info("Step 02 complete, took %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
