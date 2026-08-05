#!/usr/bin/env python3
"""
Step 01: Sample-level quality control for bulk RNA-seq

Performs:
  - Library size and gene detection rate per sample
  - Flag samples below thresholds
  - Sample correlation heatmap
  - QC summary table

Output: 01_qc.h5ad + tables/01_qc_summary.csv + figures/01_sample_correlation.png
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

from core.utils import resolve_config, save_figure, setup_logger

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("01_qc", os.path.join(cfg.log_dir, "01_qc.log"))
    log.info("Step 01: Sample-level QC")

    if os.path.exists(os.path.join(cfg.h5ad_dir, "01_qc.h5ad")):
        log.info("Skip: 01_qc.h5ad already exists. Delete it to force rerun.")
        return

    # --- Load raw data ---
    log.info("Loading %s", cfg.raw_h5ad)
    adata = sc.read(cfg.raw_h5ad)
    log.info("Loaded: %d samples x %d genes", adata.n_obs, adata.n_vars)

    # --- Calculate QC metrics ---
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True, log1p=False)
    # Rename metrics for sample-level clarity
    adata.obs["total_counts"] = adata.obs["total_counts"]
    adata.obs["n_genes_by_counts"] = adata.obs["n_genes_by_counts"]

    log.info("QC metrics computed:")
    log.info(
        "  Total counts:  min=%.0f  median=%.0f  max=%.0f",
        adata.obs["total_counts"].min(),
        adata.obs["total_counts"].median(),
        adata.obs["total_counts"].max(),
    )
    log.info(
        "  Genes detected: min=%d  median=%d  max=%d",
        int(adata.obs["n_genes_by_counts"].min()),
        int(adata.obs["n_genes_by_counts"].median()),
        int(adata.obs["n_genes_by_counts"].max()),
    )

    # --- Flag low-quality samples ---
    min_counts = getattr(cfg.bulk, "min_counts_per_sample", int(1e6))
    min_genes = (
        getattr(cfg.bulk, "min_genes_per_sample", 5000)
        if not hasattr(cfg.bulk, "min_genes_per_sample")
        else 5000
    )

    # Use qc thresholds if bulk-specific ones are not set
    if not hasattr(cfg.bulk, "min_counts_per_sample"):
        min_counts = int(1e6)
    if not hasattr(cfg.bulk, "min_genes_per_sample"):
        qc_min_genes = getattr(cfg.qc, "min_genes", 500) * 10  # 10x for bulk
        min_genes = qc_min_genes

    flagged = (adata.obs["total_counts"] < min_counts) | (
        adata.obs["n_genes_by_counts"] < min_genes
    )
    adata.obs["qc_flagged"] = flagged
    n_flagged = flagged.sum()

    if n_flagged > 0:
        log.warning("Flagged %d/%d samples below QC thresholds:", n_flagged, adata.n_obs)
        for name in adata.obs_names[flagged]:
            log.warning(
                "  %s: %d counts, %d genes",
                name,
                int(adata.obs.loc[name, "total_counts"]),
                int(adata.obs.loc[name, "n_genes_by_counts"]),
            )

    # --- Sample correlation heatmap ---
    log.info("Generating sample correlation heatmap...")
    expr_log = np.log1p(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X)
    corr = np.corrcoef(expr_log)
    sample_names = adata.obs_names.tolist()

    fig, ax = plt.subplots(
        figsize=(max(8, len(sample_names) * 0.8), max(6, len(sample_names) * 0.6))
    )
    sns.heatmap(
        corr,
        xticklabels=sample_names,
        yticklabels=sample_names,
        annot=True,
        fmt=".2f",
        cmap=cfg.plot.palette.heatmap,
        vmin=-1,
        vmax=1,
        ax=ax,
        square=True,
    )
    ax.set_title("Sample Correlation (log1p counts)")
    fig.tight_layout()
    corr_path = os.path.join(cfg.figure_dir, "01_sample_correlation.png")
    os.makedirs(cfg.figure_dir, exist_ok=True)
    save_figure(fig, corr_path, cfg=cfg, bbox_inches="tight")
    plt.close(fig)
    log.info("  Saved: %s", corr_path)

    # --- QC summary table ---
    summary = pd.DataFrame(
        {
            "sample": adata.obs_names,
            "total_counts": adata.obs["total_counts"].values,
            "n_genes": adata.obs["n_genes_by_counts"].values,
            "flagged": flagged.values,
        }
    )
    summary_path = os.path.join(cfg.table_dir, "01_qc_summary.csv")
    os.makedirs(cfg.table_dir, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    log.info("  Saved: %s", summary_path)

    # --- Optional: remove flagged samples ---
    remove_failed = getattr(cfg.bulk, "remove_failed_qc", False)
    if remove_failed and n_flagged > 0:
        log.info("Removing %d flagged samples (remove_failed_qc=True)", n_flagged)
        adata = adata[~adata.obs["qc_flagged"]].copy()
        log.info("After filtering: %d samples x %d genes", adata.n_obs, adata.n_vars)

    # --- Save ---
    qc_h5ad = os.path.join(cfg.h5ad_dir, "01_qc.h5ad")
    log.info("Saving to %s...", qc_h5ad)
    from core.utils import safe_write

    safe_write(adata, qc_h5ad, cfg=cfg)
    log.info("Step 01 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
