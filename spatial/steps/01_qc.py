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

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

from core.utils import resolve_config, safe_write, save_figure, setup_logger

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, median_abs_deviation


def compute_qc_metrics(adata, cfg, log):
    """Compute per-spot QC metrics: counts, genes, mito%, ribo%."""
    log.info("Computing QC metrics...")

    # Mitochondrial genes
    mt_mask = adata.var_names.str.startswith(cfg.qc.mt_gene_pattern)
    if cfg.qc.mt_gene_list:
        mt_mask = mt_mask | adata.var_names.isin(cfg.qc.mt_gene_list)
    adata.var["mt"] = mt_mask
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=[20],
        log1p=True,
        inplace=True,
    )

    # Complexity metric
    adata.obs["log_genes_per_umi"] = (
        np.log10(adata.obs["n_genes_by_counts"]) / np.log10(adata.obs["total_counts"])
    ).replace([np.inf, -np.inf], np.nan)

    log.info("  Median counts/spot: %.0f", adata.obs["total_counts"].median())
    log.info("  Median genes/spot:  %.0f", adata.obs["n_genes_by_counts"].median())
    log.info("  Median mito%%:       %.2f%%", adata.obs["pct_counts_mt"].median())
    if cfg.qc.is_nuclei:
        log.info(
            "  [snRNA-seq mode] Mitochondrial reads reflect cytoplasmic residue, not cell stress."
        )
    log.info("  Median complexity:   %.3f", adata.obs["log_genes_per_umi"].median())


# ══════════════════════════════════════════════════════════════════════════════
#  KDE Peak Detection
# ══════════════════════════════════════════════════════════════════════════════


def _detect_peaks(vals, log):
    """Detect peaks in nFeature distribution via KDE + peak finding.

    Spatial-conservative: min_samples=50 (vs RNA's 100).

    Returns: (n_peaks, peaks_x, peaks_y, x_range, density, is_multimodal)
    """
    vals = vals[np.isfinite(vals)]
    if len(vals) < 50:
        return 0, [], [], np.array([]), np.array([]), False
    try:
        kde = gaussian_kde(vals)
        x_range = np.linspace(vals.min(), vals.max(), 500)
        density = kde(x_range)
        peaks, props = find_peaks(density, prominence=density.max() * 0.02, distance=20)
        n_peaks = len(peaks)
        if n_peaks >= 2:
            peak_genes = x_range[peaks]
            is_multimodal = np.max(np.diff(np.sort(peak_genes))) > 500
        else:
            is_multimodal = False
        return (
            n_peaks,
            x_range[peaks].tolist(),
            density[peaks].tolist(),
            x_range,
            density,
            is_multimodal,
        )
    except Exception as e:
        log.warning("Peak detection failed: %s", e)
        return 0, [], [], np.array([]), np.array([]), False


# ══════════════════════════════════════════════════════════════════════════════
#  MAD Adaptive Thresholds
# ══════════════════════════════════════════════════════════════════════════════


def _mad_thresholds(adata, cfg, log):
    """Compute MAD-based adaptive thresholds for QC metrics.

    Returns dict with keys:
        n_genes_by_counts: (lo, hi)
        total_counts:      (None, hi)
        pct_counts_mt:     (None, hi)
        log_genes_per_umi: (lo, None)
    """
    thresholds = {}

    # ── Distribution reliability check ──
    _nf = adata.obs["n_genes_by_counts"].values.astype(np.float64)
    _nf = _nf[np.isfinite(_nf)]
    _med = np.median(_nf)
    _, _, _, _, _, _multimodal = _detect_peaks(_nf, log)
    _p10 = np.percentile(_nf, 10)
    _p90 = np.percentile(_nf, 90)

    _fallback = False
    _reason = ""
    if _multimodal and _med < 2000:
        _fallback = True
        _reason = f"bimodal distribution (median={_med:.0f} < 2000)"
    elif (_p90 - _p10) > 2000 and _med < 1500:
        _fallback = True
        _reason = f"wide unimodal (P10={_p10:.0f}, P90={_p90:.0f}, median={_med:.0f} < 1500)"

    if _fallback:
        log.warning("MAD reliability check FAILED: %s", _reason)
        log.warning(
            "Falling back to hard thresholds — median too low for MAD normality assumption."
        )
        return _hard_thresholds(cfg, log, adata=adata)

    # ---- n_genes_by_counts ----
    vals = adata.obs["n_genes_by_counts"].values.astype(np.float64)
    med = np.median(vals)
    mad = median_abs_deviation(vals, scale="normal")
    _maturity_mad = {"developing": 5.0}.get(cfg.tissue_maturity, cfg.qc.mad_n_mads)
    lo_mad = max(med - _maturity_mad * mad, 0)
    hi_mad = med + _maturity_mad * mad
    lo = max(lo_mad, cfg.qc.min_genes)
    if lo > cfg.qc.min_genes:
        lo_orig = lo
        lo = cfg.qc.min_genes
        log.warning(
            "  MAD lower bound (%.0f) exceeds min_genes (%.0f) — clamping to min_genes.",
            lo_orig,
            cfg.qc.min_genes,
        )
    if lo > med * 0.80:
        lo_orig = lo
        lo = med * 0.80
        log.warning(
            "  SAFETY CAP: min_genes=%.0f exceeds 80%% of median (median=%.0f, cap=%.0f) — clamping.",
            lo_orig,
            med,
            med * 0.80,
            cfg.tissue_maturity,
        )
    _safe_floor = (
        cfg.qc.min_mad_upper_genes_nuclei if cfg.qc.is_nuclei else cfg.qc.min_mad_upper_genes
    )
    hi = min(max(hi_mad, _safe_floor), cfg.qc.max_genes)
    thresholds["n_genes_by_counts"] = (lo, hi)
    log.info(
        "  n_genes_by_counts: median=%.0f, MAD=%.0f  ->  (lo=%.0f, hi=%.0f)  [adaptive x%.1f, %s]",
        med,
        mad,
        lo,
        hi,
        _maturity_mad,
        cfg.tissue_maturity,
    )

    # ── Distribution report ──
    try:
        n_peaks, _, _, _, _, is_multimodal = _detect_peaks(
            adata.obs["n_genes_by_counts"].values.astype(np.float64), log
        )
        p10 = np.percentile(_nf, 10)
        p90 = np.percentile(_nf, 90)
        log.info(
            "  nFeature distribution: %d KDE peak%s (P10=%.0f, P90=%.0f)%s",
            n_peaks,
            "" if n_peaks == 1 else "s",
            p10,
            p90,
            " — MULTIMODAL" if is_multimodal else "",
        )
    except Exception as e:
        log.warning("Reliability check failed: %s", e)

    # ---- total_counts ----
    if cfg.expression_type == "raw_counts":
        vals = adata.obs["total_counts"].values.astype(np.float64)
        med = np.median(vals)
        mad = median_abs_deviation(vals, scale="normal")
        hi = med + cfg.qc.ncount_max_mad * mad
        thresholds["total_counts"] = (None, hi)
        log.info(
            "  total_counts:       median=%.0f, MAD=%.0f  ->  hi=%.0f  [adaptive x%.1f]",
            med,
            mad,
            hi,
            cfg.qc.ncount_max_mad,
        )
    else:
        thresholds["total_counts"] = (None, None)
        log.info("  total_counts:       (skipped — expression_type=%s)", cfg.expression_type)

    # ---- pct_counts_mt ----
    if cfg.qc.is_nuclei:
        hi = cfg.qc.max_pct_mito_nuclei
        log.info("  pct_counts_mt:      hi=%.2f%%  [snRNA-seq: fixed threshold, MAD skipped]", hi)
    else:
        vals = adata.obs["pct_counts_mt"].values.astype(np.float64)
        med = np.median(vals)
        mad = median_abs_deviation(vals, scale="normal")
        hi_mad = med + cfg.qc.mad_n_mads * mad
        hi = min(hi_mad, cfg.qc.max_pct_mito)
        log.info(
            "  pct_counts_mt:      median=%.2f%%, MAD=%.2f%% ->  hi=%.2f%%  [adaptive, factor=%.1f]",
            med,
            mad,
            hi,
            cfg.qc.mad_n_mads,
        )
    thresholds["pct_counts_mt"] = (None, hi)

    # ---- log_genes_per_umi ----
    if cfg.expression_type == "raw_counts":
        vals = adata.obs["log_genes_per_umi"].values.astype(np.float64)
        finite = vals[np.isfinite(vals)]
        med = np.median(finite)
        mad = median_abs_deviation(finite, scale="normal")
        lo_mad = max(med - _maturity_mad * mad, 0)
        lo = max(lo_mad, cfg.qc.min_genes_per_umi)
        if lo > cfg.qc.min_genes_per_umi:
            lo_orig = lo
            lo = cfg.qc.min_genes_per_umi
            log.warning(
                "  MAD lower bound (%.4f) exceeds min_genes_per_umi (%.4f) — clamping.",
                lo_orig,
                cfg.qc.min_genes_per_umi,
            )
        thresholds["log_genes_per_umi"] = (lo, None)
        log.info(
            "  log_genes_per_umi:  median=%.4f, MAD=%.4f ->  lo=%.4f  [adaptive x%.1f, %s]",
            med,
            mad,
            lo,
            _maturity_mad,
            cfg.tissue_maturity,
        )
    else:
        thresholds["log_genes_per_umi"] = (None, None)
        log.info("  log_genes_per_umi:  (skipped — expression_type=%s)", cfg.expression_type)

    return thresholds


def _hard_thresholds(cfg, log, adata=None):
    """Build hard threshold dict from config values."""
    is_native = cfg.expression_type == "raw_counts"
    thresholds = {
        "n_genes_by_counts": (cfg.qc.min_genes, cfg.qc.max_genes),
        "total_counts": (None, None),
        "pct_counts_mt": (
            None,
            cfg.qc.max_pct_mito_nuclei if cfg.qc.is_nuclei else cfg.qc.max_pct_mito,
        ),
        "log_genes_per_umi": (cfg.qc.min_genes_per_umi, None) if is_native else (None, None),
    }
    log.info("  Using hard thresholds from config:")
    log.info("    n_genes_by_counts:  lo=%d, hi=%d", cfg.qc.min_genes, cfg.qc.max_genes)
    log.info("    total_counts:       (no limit)")
    log.info("    pct_counts_mt:      hi=%.1f%%", cfg.qc.max_pct_mito)
    if is_native:
        log.info("    log_genes_per_umi:  lo=%.4f", cfg.qc.min_genes_per_umi)
    else:
        log.info("    log_genes_per_umi:  (skipped — expression_type=%s)", cfg.expression_type)

    # ── Safety cap ──
    if adata is not None:
        nf = adata.obs["n_genes_by_counts"].values.astype(np.float64)
        nf = nf[np.isfinite(nf)]
        if len(nf) >= 10:
            med = np.median(nf)
            lo, hi = thresholds["n_genes_by_counts"]
            if lo > med * 0.80:
                thresholds["n_genes_by_counts"] = (med * 0.80, hi)
                log.warning(
                    "  SAFETY CAP: min_genes=%.0f exceeds 80%% of median (%.0f) — clamping to %.0f.",
                    lo,
                    med,
                    med * 0.80,
                )
    return thresholds


# ══════════════════════════════════════════════════════════════════════════════
#  Diagnostic Plots
# ══════════════════════════════════════════════════════════════════════════════


def _plot_qc_diagnostics(adata, thresholds, output_dir, log, cfg=None):
    """Save 3 diagnostic plots with threshold lines.

    Panels:
      A: nFeature distribution histogram
      B: nCount vs nFeature scatter (colored by % Mito)
      C: pct_mito distribution histogram
    """
    os.makedirs(output_dir, exist_ok=True)
    _dpi = cfg.plot.figure_dpi if cfg else 150
    _fmt = cfg.plot.figure_format if cfg else "png"
    _figsize = cfg.plot.qc_figure_size if cfg else [8, 5]
    _qchist = cfg.plot.palette.qc_hist if cfg else "steelblue"
    _qcthresh = cfg.plot.palette.qc_threshold if cfg else "red"
    _qcsecond = cfg.plot.palette.qc_second if cfg else "indianred"
    _pseudotime = cfg.plot.palette.pseudotime if cfg else "viridis"
    # ---- Panel A: nFeature distribution ----
    try:
        fig, ax = plt.subplots(figsize=_figsize)
        vals = adata.obs["n_genes_by_counts"].values
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, bins=100, color=_qchist, edgecolor="white", alpha=0.85)
        lo, hi = thresholds["n_genes_by_counts"]
        if lo is not None:
            ax.axvline(
                lo,
                color=_qcthresh,
                linestyle="--",
                linewidth=1.2,
                label=f"lo={lo:.0f}",
            )
        if hi is not None:
            ax.axvline(
                hi,
                color=_qcthresh,
                linestyle="--",
                linewidth=1.2,
                label=f"hi={hi:.0f}",
            )
        ax.set_xlabel("n_genes_by_counts (nFeature)")
        ax.set_ylabel("Number of spots")
        ax.set_title(f"nFeature distribution (N={adata.n_obs}, median={np.median(vals):.0f})")
        ax.legend(fontsize=9)
        fig.tight_layout()
        save_figure(
            fig, os.path.join(output_dir, "nFeature_distribution"), cfg=cfg, fmt=_fmt, dpi=_dpi
        )
        plt.close(fig)
        log.info("  Plot saved: nFeature_distribution.png")
    except Exception as e:
        log.warning("nFeature distribution plot failed: %s", e)

    # ---- Panel B: nCount vs nFeature scatter ----
    try:
        fig, ax = plt.subplots(figsize=(_figsize[0], _figsize[1] + 1))
        x = adata.obs["total_counts"].values
        y = adata.obs["n_genes_by_counts"].values
        c = adata.obs["pct_counts_mt"].values
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
        xp, yp, cp = x[finite], y[finite], c[finite]
        vmax = np.nanpercentile(c, 99) if np.isfinite(c).any() else None
        scat = ax.scatter(
            xp,
            yp,
            c=cp,
            cmap=_pseudotime,
            s=2,
            alpha=0.6,
            rasterized=True,
            vmax=vmax,
        )
        fig.colorbar(scat, ax=ax).set_label("% Mito")
        nfeat_lo, nfeat_hi = thresholds["n_genes_by_counts"]
        _, ncount_hi = thresholds["total_counts"]
        if nfeat_lo is not None:
            ax.axhline(nfeat_lo, color=_qcthresh, linestyle="--", linewidth=1.0)
        if nfeat_hi is not None:
            ax.axhline(nfeat_hi, color=_qcthresh, linestyle="--", linewidth=1.0)
        if ncount_hi is not None:
            ax.axvline(
                ncount_hi,
                color="orange",
                linestyle="--",
                linewidth=1.0,
                label=f"nCount hi={ncount_hi:.0f}",
            )
        ax.set_xlabel("total_counts (nCount)")
        ax.set_ylabel("n_genes_by_counts (nFeature)")
        ax.set_title(f"nCount vs nFeature (N={finite.sum()})")
        if ncount_hi is not None:
            ax.legend(fontsize=9)
        fig.tight_layout()
        save_figure(
            fig, os.path.join(output_dir, "nCount_vs_nFeature"), cfg=cfg, fmt=_fmt, dpi=_dpi
        )
        plt.close(fig)
        log.info("  Plot saved: nCount_vs_nFeature.png")
    except Exception as e:
        log.warning("nCount vs nFeature scatter plot failed: %s", e)

    # ---- Panel C: pct_mito distribution ----
    try:
        fig, ax = plt.subplots(figsize=_figsize)
        vals = adata.obs["pct_counts_mt"].values
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, bins=100, color=_qcsecond, edgecolor="white", alpha=0.85)
        _, hi = thresholds["pct_counts_mt"]
        if hi is not None:
            ax.axvline(
                hi,
                color=_qcthresh,
                linestyle="--",
                linewidth=1.2,
                label=f"hi={hi:.2f}%",
            )
        ax.set_xlabel("pct_counts_mt (% Mito)")
        ax.set_ylabel("Number of spots")
        ax.set_title(f"% Mito distribution (N={len(vals)}, median={np.median(vals):.2f}%)")
        if hi is not None:
            ax.legend(fontsize=9)
        fig.tight_layout()
        save_figure(
            fig, os.path.join(output_dir, "pct_mito_distribution"), cfg=cfg, fmt=_fmt, dpi=_dpi
        )
        plt.close(fig)
        log.info("  Plot saved: pct_mito_distribution.png")
    except Exception as e:
        log.warning("pct_mito distribution plot failed: %s", e)


# ══════════════════════════════════════════════════════════════════════════════
#  Threshold-dict filtering
# ══════════════════════════════════════════════════════════════════════════════


def _filter_cells(adata, thresholds, cfg, log):
    """Apply threshold dict for spot filtering (adaptive path)."""

    gf_lo, gf_hi = thresholds["n_genes_by_counts"]
    _, tc_hi = thresholds["total_counts"]
    _, mt_hi = thresholds["pct_counts_mt"]
    cpx_lo, _ = thresholds["log_genes_per_umi"]

    f_genes_low = (
        adata.obs["n_genes_by_counts"] < gf_lo
        if gf_lo is not None
        else pd.Series(False, index=adata.obs_names)
    )
    f_genes_high = (
        adata.obs["n_genes_by_counts"] > gf_hi
        if gf_hi is not None
        else pd.Series(False, index=adata.obs_names)
    )
    f_mito = (
        adata.obs["pct_counts_mt"] > mt_hi
        if mt_hi is not None
        else pd.Series(False, index=adata.obs_names)
    )
    f_count_hi = (
        adata.obs["total_counts"] > tc_hi
        if tc_hi is not None
        else pd.Series(False, index=adata.obs_names)
    )
    f_cpx = (
        adata.obs["log_genes_per_umi"] < cpx_lo
        if cpx_lo is not None
        else pd.Series(False, index=adata.obs_names)
    )
    f_any = f_genes_low | f_genes_high | f_mito | f_count_hi | f_cpx

    log.info("  Filtering breakdown:")
    if gf_lo is not None:
        log.info(
            "    n_genes < %.0f:       %6d (%.1f%%)",
            gf_lo,
            f_genes_low.sum(),
            100 * f_genes_low.mean(),
        )
    if gf_hi is not None:
        log.info(
            "    n_genes > %.0f:       %6d (%.1f%%)",
            gf_hi,
            f_genes_high.sum(),
            100 * f_genes_high.mean(),
        )
    if mt_hi is not None:
        log.info(
            "    mito > %.2f%%:        %6d (%.1f%%)", mt_hi, f_mito.sum(), 100 * f_mito.mean()
        )
    if tc_hi is not None:
        log.info(
            "    nCount > %.0f:        %6d (%.1f%%)",
            tc_hi,
            f_count_hi.sum(),
            100 * f_count_hi.mean(),
        )
    if cpx_lo is not None:
        log.info("    complexity < %.4f:   %6d (%.1f%%)", cpx_lo, f_cpx.sum(), 100 * f_cpx.mean())
    log.info("    Total (dedup):        %6d (%.1f%%)", f_any.sum(), 100 * f_any.mean())

    mask = ~f_any
    if "in_tissue" in adata.obs:
        mask = mask & adata.obs["in_tissue"].astype(bool)

    adata = adata[mask].copy()
    log.info("  After QC filtering: %d spots", adata.n_obs)
    return adata


def filter_spots(adata, cfg, log):
    """Filter spots by QC thresholds and in_tissue flag."""
    n_before = adata.n_obs

    # Tissue spot filtering (Visium: under-tissue spots only)
    if "in_tissue" in adata.obs:
        n_tissue = adata.obs["in_tissue"].sum()
        log.info(
            "Tissue spots: %d / %d (%.1f%%)",
            n_tissue,
            n_before,
            100 * n_tissue / n_before if n_before else 0,
        )

    log.info("Applying QC filtering...")
    min_g = cfg.qc.min_genes
    max_g = cfg.qc.max_genes
    max_m = cfg.qc.max_pct_mito_nuclei if cfg.qc.is_nuclei else cfg.qc.max_pct_mito
    min_cpx = cfg.qc.min_genes_per_umi

    f_genes_low = adata.obs["n_genes_by_counts"] < min_g
    f_genes_high = adata.obs["n_genes_by_counts"] > max_g
    f_mito = adata.obs["pct_counts_mt"] > max_m
    f_cpx = adata.obs["log_genes_per_umi"] < min_cpx
    f_any = f_genes_low | f_genes_high | f_mito | f_cpx

    log.info("  Filtering breakdown:")
    log.info(
        "    n_genes < %d:     %6d (%.1f%%)", min_g, f_genes_low.sum(), 100 * f_genes_low.mean()
    )
    log.info(
        "    n_genes > %d:     %6d (%.1f%%)", max_g, f_genes_high.sum(), 100 * f_genes_high.mean()
    )
    log.info("    mito > %.0f%%:     %6d (%.1f%%)", max_m, f_mito.sum(), 100 * f_mito.mean())
    log.info("    complexity < %.2f: %6d (%.1f%%)", min_cpx, f_cpx.sum(), 100 * f_cpx.mean())
    log.info("    Total (dedup):    %6d (%.1f%%)", f_any.sum(), 100 * f_any.mean())

    mask = ~f_any
    if "in_tissue" in adata.obs:
        mask = mask & adata.obs["in_tissue"].astype(bool)

    adata = adata[mask].copy()
    log.info("  After QC filtering: %d spots", adata.n_obs)
    return adata


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("01_qc", os.path.join(cfg.log_dir, "01_qc.log"))
    log.info("Step 01: QC filtering for spatial transcriptomics")

    input_path = cfg.raw_h5ad
    if not os.path.exists(input_path):
        log.error("Input not found: %s. Run Step 00 first.", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    compute_qc_metrics(adata, cfg, log)

    if cfg.qc.use_adaptive_thresholds:
        log.info(
            "Mode: adaptive (MAD × %.1f, nCount MAD × %.1f)",
            cfg.qc.mad_n_mads,
            cfg.qc.ncount_max_mad,
        )
        thresholds = _mad_thresholds(adata, cfg, log)
        fig_dir = os.path.join(cfg.figure_dir, "01_qc")
        _plot_qc_diagnostics(adata, thresholds, fig_dir, log, cfg=cfg)
        adata = _filter_cells(adata, thresholds, cfg, log)

        if adata.n_obs == 0:
            log.warning(
                "ZERO survivors with MAD thresholds — reloading raw data with hard thresholds."
            )
            adata = sc.read(input_path)
            log.info("Reloaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)
            compute_qc_metrics(adata, cfg, log)
            adata = filter_spots(adata, cfg, log)
    else:
        log.info("Mode: hard thresholds")
        adata = filter_spots(adata, cfg, log)

    sc.pp.filter_genes(adata, min_cells=cfg.qc.min_cells_per_gene)
    log.info("After gene filtering: %d genes", adata.n_vars)

    qc_out = os.path.join(cfg.h5ad_dir, "01_qc.h5ad")
    safe_write(adata, qc_out, cfg=cfg)
    log.info("Step 01 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
