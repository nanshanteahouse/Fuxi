#!/usr/bin/env python3
"""
Step 05: (Optional) Batch correction for bulk RNA-seq via pycombat.

When enabled, removes unwanted batch effects while preserving biological
contrast signal by treating ``CFG.bulk.contrast_column`` as a continuous
covariate.

Input:  01_qc.h5ad  (from CFG.h5ad_dir)
Output: 05_batch_corrected.h5ad
        figures/05_pca_comparison.png

Graceful degradation:
  - If ``CFG.bulk.batch_correct`` is False → skip with info log
  - If pycombat is not installed → log warning and return (no crash)
  - If batch_column not found in obs → log warning and skip
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

import matplotlib
import numpy as np
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core.utils import resolve_config, safe_write, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("05_batch", os.path.join(cfg.log_dir, "05_batch.log"))
    log.info("Step 05: Batch correction")

    # ── Check if batch correction is enabled ────────────────────────────
    if not cfg.bulk.batch_correct:
        log.info("Skipping batch correction (bulk.batch_correct=False)")
        return

    out_path = os.path.join(cfg.h5ad_dir, "05_batch_corrected.h5ad")
    if os.path.exists(out_path):
        log.info("Skip: %s already exists. Delete it to force rerun.", out_path)
        return

    # ── Graceful pycombat import ──────────────────────────────────────
    try:
        import pycombat
    except ImportError:
        log.warning("pycombat not installed. Install with: pip install pycombat")
        log.warning("Skipping batch correction step.")
        return

    # ── Load input ─────────────────────────────────────────────────────
    qc_h5ad = os.path.join(cfg.h5ad_dir, "01_qc.h5ad")
    if not os.path.exists(qc_h5ad):
        log.warning("Skip: %s not found. Run step 01 first.", qc_h5ad)
        return

    log.info("Loading %s", qc_h5ad)
    adata = sc.read(qc_h5ad)
    log.info("Loaded: %d samples x %d genes", adata.n_obs, adata.n_vars)

    # ── Validate batch column ──────────────────────────────────────────
    batch_col = cfg.bulk.batch_column
    if batch_col not in adata.obs.columns:
        log.warning(
            "Batch column '%s' not found in adata.obs. "
            "Available columns: %s. Skipping batch correction.",
            batch_col,
            list(adata.obs.columns),
        )
        return

    # ── Extract variables ──────────────────────────────────────────────
    contrast_col = getattr(cfg.bulk, "contrast_column", "condition")

    # Warn but continue if contrast column is absent (purely technical correction)
    if contrast_col not in adata.obs.columns:
        log.warning(
            "Contrast column '%s' not in adata.obs — correcting batches "
            "without preserving biological contrast.",
            contrast_col,
        )
        # Fallback: treat all samples as one biological group
        categorical_cols = [batch_col]
        contrast_available = False
    else:
        categorical_cols = [batch_col]
        contrast_available = True

    log.info("Batch column: %s | Contrast column: %s", batch_col, contrast_col)
    log.info("Batch levels: %s", sorted(adata.obs[batch_col].unique()))

    # ── Prepare count matrix ──────────────────────────────────────────
    log.info("Preparing count matrix for pycombat...")
    x = adata.X
    if hasattr(x, "toarray"):
        x = x.toarray()
    counts = np.asarray(x, dtype=np.float64)

    # pycombat operates on samples-as-columns (genes x samples)
    # adata.X is samples x genes → transpose to genes x samples
    counts_t = counts.T  # shape: (n_genes, n_samples)

    # ── Build covariate DataFrame ──────────────────────────────────────
    sample_info = (
        adata.obs[[batch_col, contrast_col]].copy()
        if contrast_available
        else adata.obs[[batch_col]].copy()
    )
    # Ensure all columns are string-typed for pycombat
    for col in sample_info.columns:
        sample_info[col] = sample_info[col].astype(str)

    log.info("Running pycombat... (this may take a moment)")
    try:
        corrected_t = pycombat.pycombat(
            counts_t,
            sample_info,
            categorical_cols=categorical_cols,
        )
    except Exception as e:
        log.warning("pycombat failed: %s", e)
        log.warning("Skipping batch correction.")
        return

    # Transpose back: samples x genes
    corrected = corrected_t.T.astype(np.float32)
    log.info("Batch correction complete: corrected matrix shape %s", corrected.shape)

    # ── Store in output AnnData ────────────────────────────────────────
    out_adata = adata.copy()
    out_adata.X = corrected
    out_adata.uns["batch_correction"] = {
        "method": "pycombat",
        "batch_column": batch_col,
        "contrast_column": contrast_col if contrast_available else None,
    }

    # ── Before / After PCA plots ───────────────────────────────────────
    os.makedirs(cfg.figure_dir, exist_ok=True)
    _pca_comparison(adata.X, out_adata.X, adata.obs, batch_col, cfg.figure_dir, log)

    # ── Save ───────────────────────────────────────────────────────────
    log.info("Saving to %s...", out_path)
    safe_write(out_adata, out_path, cfg=cfg)
    log.info("Step 05 complete, took %.1fs", time.time() - t0)


# ── Plot helper ──────────────────────────────────────────────────────────


def _pca_comparison(x_before, x_after, obs, batch_col, figure_dir, log):
    """Side-by-side PCA plots: before and after batch correction."""
    try:
        n_comps = min(50, x_before.shape[0] - 1)
        if n_comps < 2:
            log.warning("Skipping PCA comparison: need ≥2 samples (got %d)", x_before.shape[0])
            return

        # Before
        log.info("Computing PCA before correction...")
        adata_before = sc.AnnData(X=x_before, obs=obs)
        sc.pp.pca(adata_before, n_comps=n_comps)
        pca_before = adata_before.obsm["X_pca"]

        # After
        log.info("Computing PCA after correction...")
        adata_after = sc.AnnData(X=x_after, obs=obs)
        sc.pp.pca(adata_after, n_comps=n_comps)
        pca_after = adata_after.obsm["X_pca"]

        # ── Plot ───────────────────────────────────────────────────────
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        colors = _get_colors(obs, batch_col)
        for label in obs[batch_col].unique():
            mask = obs[batch_col] == label
            ax1.scatter(
                pca_before[mask, 0],
                pca_before[mask, 1],
                c=[colors[label]],
                label=label,
                edgecolors="black",
                linewidths=0.5,
                s=60,
            )
            ax2.scatter(
                pca_after[mask, 0],
                pca_after[mask, 1],
                c=[colors[label]],
                label=label,
                edgecolors="black",
                linewidths=0.5,
                s=60,
            )

        ax1.set_title("Before batch correction")
        ax1.set_xlabel(f"PC1 ({_var_explained(adata_before, 0):.1f}%)")
        ax1.set_ylabel(f"PC2 ({_var_explained(adata_before, 1):.1f}%)")
        ax1.legend(fontsize=8)

        ax2.set_title("After batch correction")
        ax2.set_xlabel(f"PC1 ({_var_explained(adata_after, 0):.1f}%)")
        ax2.set_ylabel(f"PC2 ({_var_explained(adata_after, 1):.1f}%)")
        ax2.legend(fontsize=8)

        fig.tight_layout()

        comparison_path = os.path.join(figure_dir, "05_pca_comparison.png")
        fig.savefig(comparison_path, dpi=150, bbox_inches="tight")
        log.info("Saved: %s", comparison_path)

        plt.close(fig)

    except Exception as e:
        log.warning("PCA comparison plot failed: %s", e)


def _get_colors(obs, col):
    """Return a color map for unique values of an obs column."""
    unique = obs[col].unique()
    palette = plt.cm.tab10(np.linspace(0, 1, len(unique)))
    return dict(zip(unique, palette))


def _var_explained(adata, pc_idx):
    """Return variance explained percentage for a given PC."""
    if "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        return adata.uns["pca"]["variance_ratio"][pc_idx] * 100
    return 0.0


if __name__ == "__main__":
    main()
