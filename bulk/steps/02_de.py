#!/usr/bin/env python3
"""
Step 02: Differential expression analysis via PyDESeq2 for bulk RNA-seq.

Performs:
  - Design formula validation (terms must exist in adata.obs)
  - Contrast validation (treatment/baseline levels must exist)
  - DESeq2 normalization (median of ratios) + dispersion estimation + Wald test
  - Optional LFC shrinkage (apeglm)
  - DEG identification (padj < alpha)
  - Volcano plot with top 20 significant genes labeled
  - MA plot (mean of normalized counts vs log2 fold change)

Input:  01_qc.h5ad (from CFG.h5ad_dir)
Output: 02_de.h5ad           — normalized counts stored in .X
        tables/02_de_results.csv
        tables/02_de_significant.csv
        figures/02_volcano.png
        figures/02_ma_plot.png
"""
import sys
import os
import time
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core.utils import setup_logger, resolve_config, safe_write


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("02_de", os.path.join(CFG.log_dir, "02_de.log"))
    log.info("Step 02: DESeq2 differential expression")

    de_h5ad_path = os.path.join(CFG.h5ad_dir, "02_de.h5ad")

    # ── Resume check ──────────────────────────────────────────────────
    if os.path.exists(de_h5ad_path):
        log.info("Skip: %s already exists. Delete it to force rerun.", de_h5ad_path)
        return

    # ── Load input ────────────────────────────────────────────────────
    qc_h5ad = os.path.join(CFG.h5ad_dir, "01_qc.h5ad")
    log.info("Loading %s", qc_h5ad)
    adata = sc.read(qc_h5ad)
    log.info("Loaded: %d samples x %d genes", adata.n_obs, adata.n_vars)

    # ── Validate design formula ───────────────────────────────────────
    bulk = CFG.bulk
    design = bulk.design

    # Parse R-style formula "~term1 + term2" → extract term names
    # Strip the tilde and split on '+'
    design_raw = design.lstrip("~").strip()
    design_terms = [t.strip() for t in design_raw.split("+") if t.strip()]
    log.info("Design formula: %s → parsed terms: %s", design, design_terms)

    missing_terms = [t for t in design_terms if t not in adata.obs.columns]
    if missing_terms:
        log.error(
            "Design terms %s not found in adata.obs. Available columns: %s",
            missing_terms, list(adata.obs.columns),
        )
        sys.exit(1)
    log.info("All design terms present in adata.obs")

    # ── Validate contrast ─────────────────────────────────────────────
    contrast_col = bulk.contrast_column
    if contrast_col not in adata.obs.columns:
        log.error(
            "contrast_column '%s' not found in adata.obs. Available: %s",
            contrast_col, list(adata.obs.columns),
        )
        sys.exit(1)

    treatment = bulk.contrast_treatment
    baseline = bulk.contrast_baseline
    if not treatment or not baseline:
        log.error(
            "contrast_treatment and contrast_baseline must be non-empty "
            "(got: '%s', '%s')",
            treatment, baseline,
        )
        sys.exit(1)

    valid_levels = set(adata.obs[contrast_col].unique())
    for val_name, val in [("contrast_treatment", treatment),
                          ("contrast_baseline", baseline)]:
        if val not in valid_levels:
            log.error(
                "bulk.%s='%s' not found in adata.obs['%s']. "
                "Valid values: %s",
                val_name, val, contrast_col, sorted(valid_levels),
            )
            sys.exit(1)
    log.info(
        "Contrast: %s = '%s' (treatment) vs '%s' (baseline)",
        contrast_col, treatment, baseline,
    )

    # ── Extract raw integer count matrix ──────────────────────────────
    log.info("Extracting count matrix...")
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    # Warn if data appears non-integer
    if not np.issubdtype(X.dtype, np.integer) and not np.allclose(X, X.astype(int)):
        log.warning(
            "adata.X contains non-integer values (dtype=%s). "
            "Rounding to integer for DESeq2.",
            X.dtype,
        )

    X_int = np.round(X).astype(int)

    counts_df = pd.DataFrame(
        X_int,
        index=adata.obs_names,
        columns=adata.var_names,
    )
    log.info("Counts matrix: %d samples x %d genes", *counts_df.shape)

    # ── Build metadata DataFrame ──────────────────────────────────────
    # Include all design terms + the contrast column
    meta_cols = list(set(design_terms + [contrast_col]))
    metadata_df = adata.obs[meta_cols].copy()
    # Convert non-numeric columns to categorical (factors for PyDESeq2)
    for col in metadata_df.columns:
        if not pd.api.types.is_numeric_dtype(metadata_df[col].dtype):
            metadata_df[col] = metadata_df[col].astype("category")
    log.info("Metadata columns: %s", list(metadata_df.columns))

    # ── Filter low-count genes ────────────────────────────────────────
    min_counts = bulk.min_counts_per_gene
    gene_totals = counts_df.sum(axis=0)
    genes_to_keep = gene_totals >= min_counts
    n_removed = (~genes_to_keep).sum()
    if n_removed > 0:
        log.info(
            "Removing %d/%d genes with total count < %d",
            n_removed, counts_df.shape[1], min_counts,
        )
    else:
        log.info("All %d genes pass the min_counts_per_gene (%d) threshold",
                  counts_df.shape[1], min_counts)

    counts_df = counts_df.loc[:, genes_to_keep]
    log.info("After filtering: %d genes remain", counts_df.shape[1])

    # ── DESeq2 ────────────────────────────────────────────────────────
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        log.error("pydeseq2 not installed — run: pip install pydeseq2")
        sys.exit(1)

    log.info("Creating DeseqDataSet...")
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design=design,
        n_cpus=bulk.n_jobs if bulk.n_jobs > 0 else None,
    )

    log.info("Running DESeq2 (size factor estimation + dispersion + Wald test)...")
    dds.deseq2()
    log.info("DESeq2 complete")

    # ── Extract results via DeseqStats ─────────────────────────────────
    log.info("Extracting DESeq2 results...")
    stat = DeseqStats(
        dds,
        contrast=[bulk.contrast_column, bulk.contrast_treatment, bulk.contrast_baseline],
        alpha=bulk.alpha,
    )

    stat.summary()

    # ── Optional LFC shrinkage (call after summary so SE is populated) ─
    if bulk.lfc_shrink:
        log.info("Applying LFC shrinkage...")
        # Determine the LFC coefficient name (non-intercept contrast column)
        lfc_coeffs = [
            c for c in dds.varm["LFC"].columns
            if c != "Intercept" and c.startswith(bulk.contrast_column)
        ]
        coeff = lfc_coeffs[0] if lfc_coeffs else None
        if coeff:
            stat.lfc_shrink(coeff=coeff)
            log.info("LFC shrinkage applied on coefficient '%s'", coeff)
        else:
            log.warning(
                "Could not determine LFC coefficient from varm['LFC'] columns: %s. "
                "Skipping LFC shrinkage.",
                list(dds.varm["LFC"].columns),
            )

    results_df = stat.results_df.copy()

    # Move gene index to a column for CSV export
    results_df = results_df.reset_index().rename(columns={"index": "gene"})
    log.info("Raw results: %d genes", len(results_df))

    # ── Filter significant ────────────────────────────────────────────
    alpha = bulk.alpha
    sig_mask = results_df["padj"].fillna(1) < alpha
    sig_df = results_df.loc[sig_mask].copy()
    n_sig = len(sig_df)
    n_total = len(results_df)
    log.info(
        "Significant DEGs (padj < %s): %d / %d (%.1f%%)",
        alpha, n_sig, n_total,
        100.0 * n_sig / n_total if n_total > 0 else 0.0,
    )

    # ── Export CSVs ───────────────────────────────────────────────────
    os.makedirs(CFG.table_dir, exist_ok=True)

    results_path = os.path.join(CFG.table_dir, "02_de_results.csv")
    results_df.to_csv(results_path, index=False)
    log.info("Exported: %s (%d rows)", results_path, n_total)

    sig_path = os.path.join(CFG.table_dir, "02_de_significant.csv")
    sig_df.to_csv(sig_path, index=False)
    log.info("Exported: %s (%d rows)", sig_path, n_sig)

    # ── Volcano plot ──────────────────────────────────────────────────
    os.makedirs(CFG.figure_dir, exist_ok=True)
    _volcano_plot(results_df, alpha, treatment, baseline, CFG.figure_dir, log)

    # ── MA plot ───────────────────────────────────────────────────────
    _ma_plot(results_df, alpha, treatment, baseline, CFG.figure_dir, log)

    # ── Store normalized counts and save ──────────────────────────────
    log.info("Storing normalized counts...")

    # After deseq2(), normalized counts are available in dds.layers
    if hasattr(dds, "layers") and "normed_counts" in dds.layers:
        normalized = dds.layers["normed_counts"]  # samples x genes
        log.info("Using dds.layers['normed_counts']")
    else:
        # Fallback: use dds.X which is the normalized count matrix in recent pydeseq2
        normalized = dds.X
        log.info("Falling back to dds.X for normalized counts")

    # Create output AnnData with only genes that passed filtering
    kept_genes = counts_df.columns.tolist()
    out_adata = adata[:, kept_genes].copy()
    out_adata.X = normalized.astype(np.float32)

    # Carry forward the contrast information in uns
    out_adata.uns["de"] = {
        "contrast_column": bulk.contrast_column,
        "contrast_treatment": bulk.contrast_treatment,
        "contrast_baseline": bulk.contrast_baseline,
        "design": bulk.design,
        "alpha": bulk.alpha,
        "lfc_shrink": bulk.lfc_shrink,
        "n_significant": n_sig,
        "n_tested": n_total,
    }

    log.info("Saving to %s...", de_h5ad_path)
    safe_write(out_adata, de_h5ad_path, cfg=CFG)
    log.info("Step 02 complete, took %.1fs", time.time() - t0)


# ── Plot helpers ──────────────────────────────────────────────────────────

def _volcano_plot(results_df, alpha, treatment, baseline, figure_dir, log):
    """Generate volcano plot: -log10(padj) vs log2FoldChange."""
    try:
        plot_df = results_df.dropna(subset=["padj"]).copy()
        # Clip to avoid -inf in -log10(0)
        plot_df = plot_df[plot_df["padj"] > 0]

        if len(plot_df) == 0:
            log.warning("No valid padj values — skipping volcano plot")
            return

        plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"])

        fig, ax = plt.subplots(figsize=(10, 8))

        sig_mask = plot_df["padj"] < alpha

        # Non-significant
        ax.scatter(
            plot_df.loc[~sig_mask, "log2FoldChange"],
            plot_df.loc[~sig_mask, "neg_log10_padj"],
            s=3, alpha=0.3, color="grey", label="NS",
        )
        # Significant
        ax.scatter(
            plot_df.loc[sig_mask, "log2FoldChange"],
            plot_df.loc[sig_mask, "neg_log10_padj"],
            s=5, alpha=0.6, color="red",
            label=f"padj<{alpha}",
        )

        # Threshold lines
        ax.axhline(-np.log10(alpha), color="blue", linestyle="--", linewidth=0.8)
        ax.axvline(-1, color="gray", linestyle="--", linewidth=0.5)
        ax.axvline(1, color="gray", linestyle="--", linewidth=0.5)

        # Label top 20 significant genes by padj
        sig_plot = plot_df.loc[sig_mask]
        if len(sig_plot) > 0:
            top = sig_plot.nsmallest(20, "padj")
            for _, row in top.iterrows():
                label = row.get("gene", row.name)
                ax.annotate(
                    label,
                    (row["log2FoldChange"], row["neg_log10_padj"]),
                    fontsize=7,
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.3),
                )

        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10 adjusted p-value")
        ax.set_title(f"Volcano Plot: {treatment} vs {baseline}")
        ax.legend(loc="upper right")

        vol_path = os.path.join(figure_dir, "02_volcano.png")
        fig.savefig(vol_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("Volcano plot saved: %s", vol_path)

    except Exception as e:
        log.warning("Volcano plot failed: %s", e)


def _ma_plot(results_df, alpha, treatment, baseline, figure_dir, log):
    """Generate MA plot: mean of normalized counts vs log2 fold change."""
    try:
        plot_df = results_df.dropna(subset=["baseMean", "log2FoldChange"]).copy()

        if len(plot_df) == 0:
            log.warning("No valid baseMean/log2FoldChange values — skipping MA plot")
            return

        fig, ax = plt.subplots(figsize=(10, 8))

        sig_ma = None
        if "padj" in plot_df.columns:
            sig_ma = plot_df["padj"].fillna(1) < alpha

        # Non-significant
        ns_ma = plot_df if sig_ma is None else plot_df[~sig_ma]
        ax.scatter(
            ns_ma["baseMean"], ns_ma["log2FoldChange"],
            s=3, alpha=0.3, color="grey",
        )

        # Significant
        if sig_ma is not None and sig_ma.any():
            s_ma = plot_df[sig_ma]
            ax.scatter(
                s_ma["baseMean"], s_ma["log2FoldChange"],
                s=5, alpha=0.6, color="red",
                label=f"padj<{alpha}",
            )
            ax.legend(loc="upper right")

        # Base line at logFC=0
        ax.axhline(0, color="black", linestyle="-", linewidth=0.5)

        ax.set_xlabel("Mean of normalized counts")
        ax.set_ylabel("log2 Fold Change")
        ax.set_xscale("log")
        ax.set_title(f"MA Plot: {treatment} vs {baseline}")

        ma_path = os.path.join(figure_dir, "02_ma_plot.png")
        fig.savefig(ma_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("MA plot saved: %s", ma_path)

    except Exception as e:
        log.warning("MA plot failed: %s", e)


if __name__ == "__main__":
    main()
