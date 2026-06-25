#!/usr/bin/env python3
"""
Step 01: QC filtering for spatial transcriptomics
====================================================
  1. Compute QC metrics (counts, genes, mito%)
  2. Filter spots by QC thresholds
  3. Filter genes by min_cells
  4. Tissue spot detection (if not already present)

Input:  00_raw.h5ad
Output: 01_qc.h5ad
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import scanpy as sc


def compute_qc_metrics(adata, cfg, log):
    """Compute per-spot QC metrics: counts, genes, mito%, ribo%."""
    log.info("Computing QC metrics...")

    # Mitochondrial genes
    mt_mask = adata.var_names.str.startswith(cfg.mt_gene_pattern)
    if cfg.mt_gene_list:
        mt_mask = mt_mask | adata.var_names.isin(cfg.mt_gene_list)
    adata.var['mt'] = mt_mask
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo'],
        percent_top=[20], log1p=True, inplace=True,
    )

    # Complexity metric
    adata.obs['log_genes_per_umi'] = (
        np.log10(adata.obs['n_genes_by_counts'])
        / np.log10(adata.obs['total_counts'])
    ).replace([np.inf, -np.inf], np.nan)

    log.info("  Median counts/spot: %.0f", adata.obs['total_counts'].median())
    log.info("  Median genes/spot:  %.0f", adata.obs['n_genes_by_counts'].median())
    log.info("  Median mito%%:       %.2f%%", adata.obs['pct_counts_mt'].median())
    log.info("  Median complexity:   %.3f", adata.obs['log_genes_per_umi'].median())


def filter_spots(adata, cfg, log):
    """Filter spots by QC thresholds and in_tissue flag."""
    n_before = adata.n_obs

    # Tissue spot filtering (Visium: under-tissue spots only)
    if 'in_tissue' in adata.obs:
        n_tissue = adata.obs['in_tissue'].sum()
        log.info("Tissue spots: %d / %d (%.1f%%)",
                 n_tissue, n_before, 100 * n_tissue / n_before if n_before else 0)

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
    log.info("    n_genes < %d:     %6d (%.1f%%)", min_g, f_genes_low.sum(), 100 * f_genes_low.mean())
    log.info("    n_genes > %d:     %6d (%.1f%%)", max_g, f_genes_high.sum(), 100 * f_genes_high.mean())
    log.info("    mito > %.0f%%:     %6d (%.1f%%)", max_m, f_mito.sum(), 100 * f_mito.mean())
    log.info("    complexity < %.2f: %6d (%.1f%%)", min_cpx, f_cpx.sum(), 100 * f_cpx.mean())
    log.info("    Total (dedup):    %6d (%.1f%%)", f_any.sum(), 100 * f_any.mean())

    mask = ~f_any
    if 'in_tissue' in adata.obs:
        mask = mask & adata.obs['in_tissue'].astype(bool)

    adata = adata[mask].copy()
    log.info("  After QC filtering: %d spots", adata.n_obs)
    return adata


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("01_qc", os.path.join(CFG.log_dir, "01_qc.log"))
    log.info("Step 01: QC filtering for spatial transcriptomics")

    input_path = CFG.raw_h5ad
    if not os.path.exists(input_path):
        log.error("Input not found: %s. Run Step 00 first.", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    compute_qc_metrics(adata, CFG, log)
    adata = filter_spots(adata, CFG, log)

    # Gene filtering
    sc.pp.filter_genes(adata, min_cells=CFG.min_cells_per_gene)
    log.info("After gene filtering: %d genes", adata.n_vars)

    # Save
    qc_out = os.path.join(CFG.h5ad_dir, "01_qc.h5ad")
    safe_write(adata, qc_out, cfg=CFG)
    log.info("Step 01 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
