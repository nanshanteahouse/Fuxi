#!/usr/bin/env python3
"""
Step 02: QC filtering (doublets already removed in Step 01)
============================================================
Best practices:
  1. QC metrics (mito%, ribo%, complexity)
  2. Filter predicted_doublet cells (from Step 01)
  3. Adaptive MAD or global threshold filtering
  4. Diagnostic plots (always generated for audit trail)

输入: 01_doublet.h5ad (含 doublet_scores, predicted_doublet 列)
输出: 02_qc.h5ad (过滤后的细胞 + QC 指标)
      {figure_dir}/02_qc/ (3 张诊断图)
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

from core.utils import resolve_config, setup_logger

try:
    import hdf5plugin  # zstd 滤镜；写 zstd 压缩文件需要
except ImportError:
    hdf5plugin = None
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, median_abs_deviation


def _detect_peaks(vals):
    """Detect peaks in nFeature distribution via KDE + peak finding.

    Returns: (n_peaks, peaks_x, peaks_y, x_range, density, is_multimodal)
    """
    vals = vals[np.isfinite(vals)]
    if len(vals) < 100:
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
    except Exception:
        return 0, [], [], np.array([]), np.array([]), False


# ══════════════════════════════════════════════════════════════════════════════
#  MAD 自适应阈值
# ══════════════════════════════════════════════════════════════════════════════


def _mad_thresholds(adata, cfg, log):
    """用 MAD (Median Absolute Deviation) 为每个 QC 指标计算自适应上下界。

    返回:
        dict: {'n_genes_by_counts': (lo, hi),
               'total_counts': (None, hi),       # nCount 通常仅上限
               'pct_counts_mt': (None, hi),
               'log_genes_per_umi': (lo, None)}
    """
    thresholds = {}

    # ── 分布可靠性检测：MAD 假设单峰、中位数足够的分布 ──
    _nf = adata.obs["n_genes_by_counts"].values.astype(np.float64)
    _nf = _nf[np.isfinite(_nf)]
    _med = np.median(_nf)
    _, _, _, _, _, _multimodal = _detect_peaks(_nf)
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

    # ---- n_genes_by_counts (nFeature_RNA) ----
    vals = adata.obs["n_genes_by_counts"].values.astype(np.float64)
    med = np.median(vals)
    mad = median_abs_deviation(vals, scale="normal")
    _maturity_mad = {"developing": 5.0}.get(cfg.tissue_maturity, cfg.qc.mad_n_mads)
    lo_mad = max(med - _maturity_mad * mad, 0)
    hi_mad = med + _maturity_mad * mad
    # 硬阈值做地板/天花板
    lo = max(lo_mad, cfg.qc.min_genes)
    # Safety clamp: MAD 下界不应超过用户硬阈值（参考 ddqc, Subramanian 2022）
    if lo > cfg.qc.min_genes:
        lo_orig = lo
        lo = cfg.qc.min_genes
        log.warning(
            "  MAD lower bound (%.0f) exceeds min_genes (%.0f) — clamping to min_genes. MAD is unreliable for this distribution; consider use_adaptive_thresholds=False.",
            lo_orig,
            cfg.qc.min_genes,
        )
    # Safety cap: lower bound must not exceed 80% of median — prevents over-filtering
    # in low-gene-count tissues (e.g. developing retina, fetal samples).
    if lo > med * 0.80:
        lo_orig = lo
        lo = med * 0.80
        log.warning(
            "  SAFETY CAP: min_genes=%.0f exceeds 80%% of median (median=%.0f, cap=%.0f) — "
            "clamping to %.0f. Over-filtering risk in low-gene-count tissue "
            "(tissue_maturity=%s).",
            lo_orig,
            med,
            med * 0.80,
            lo,
            cfg.tissue_maturity,
        )

    _safe_floor = (
        cfg.qc.min_mad_upper_genes_nuclei if cfg.qc.is_nuclei else cfg.qc.min_mad_upper_genes
    )
    hi = min(max(hi_mad, _safe_floor), cfg.qc.max_genes)
    thresholds["n_genes_by_counts"] = (lo, hi)
    log.info(
        "  n_genes_by_counts: median=%.0f, MAD=%.0f  →  (lo=%.0f, hi=%.0f)  [adaptive ×%.1f, %s]",
        med,
        mad,
        lo,
        hi,
        _maturity_mad,
        cfg.tissue_maturity,
    )

    # ── 分布检测 ──
    try:
        vals_nf = adata.obs["n_genes_by_counts"].values.astype(np.float64)
        n_peaks, _, _, _, _, is_multimodal = _detect_peaks(vals_nf)
        # Also report P10/P90 for context
        p10 = np.percentile(vals_nf[np.isfinite(vals_nf)], 10)
        p90 = np.percentile(vals_nf[np.isfinite(vals_nf)], 90)
        log.info(
            "  nFeature distribution: %d KDE peak%s (P10=%.0f, P90=%.0f)%s",
            n_peaks,
            "" if n_peaks == 1 else "s",
            p10,
            p90,
            " — MULTIMODAL: MAD may over-filter the high-expression population"
            if is_multimodal
            else "",
        )
    except Exception:
        pass

    # ---- total_counts (nCount_RNA) ----
    # 非 raw_counts 数据不设 total_counts 上限
    if cfg.expression_type == "raw_counts":
        vals = adata.obs["total_counts"].values.astype(np.float64)
        med = np.median(vals)
        mad = median_abs_deviation(vals, scale="normal")
        hi = med + cfg.qc.ncount_max_mad * mad
        thresholds["total_counts"] = (None, hi)
        log.info(
            "  total_counts:       median=%.0f, MAD=%.0f  →  hi=%.0f  [adaptive ×%.1f]",
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

    # ---- log_genes_per_umi (complexity) ----
    # 非 raw_counts 数据下复杂度指标无解释力，跳过
    if cfg.expression_type == "raw_counts":
        vals = adata.obs["log_genes_per_umi"].values.astype(np.float64)
        finite = vals[np.isfinite(vals)]
        med = np.median(finite)
        mad = median_abs_deviation(finite, scale="normal")
        lo_mad = max(med - _maturity_mad * mad, 0)
        lo = max(lo_mad, cfg.qc.min_genes_per_umi)
        # Safety clamp: MAD 下界不应超过用户硬阈值
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
            "  log_genes_per_umi:  median=%.4f, MAD=%.4f →  lo=%.4f  [adaptive ×%.1f, %s]",
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
    """从 Config 构建硬阈值字典（现有行为）。"""
    is_native = cfg.expression_type == "raw_counts"
    thresholds = {
        "n_genes_by_counts": (cfg.qc.min_genes, cfg.qc.max_genes),
        "total_counts": (None, None),  # raw_counts 下也未启用硬上限
        "pct_counts_mt": (
            None,
            cfg.qc.max_pct_mito_nuclei if cfg.qc.is_nuclei else cfg.qc.max_pct_mito,
        ),
        "log_genes_per_umi": (cfg.qc.min_genes_per_umi, None) if is_native else (None, None),
    }
    log.info("  Using hard thresholds from config:")
    log.info("    n_genes_by_counts:  lo=%d, hi=%d", cfg.qc.min_genes, cfg.qc.max_genes)
    log.info("    total_counts:       (no limit)")
    if is_native:
        log.info("    pct_counts_mt:      hi=%.1f%%", cfg.qc.max_pct_mito)
        log.info("    log_genes_per_umi:  lo=%.4f", cfg.qc.min_genes_per_umi)
    else:
        log.info("    pct_counts_mt:      hi=%.1f%%", cfg.qc.max_pct_mito)
        log.info("    log_genes_per_umi:  (skipped — expression_type=%s)", cfg.expression_type)

    # ── 安全帽（与 MAD 模式一致）：防止 min_genes 过高 ──
    if adata is not None:
        nf = adata.obs["n_genes_by_counts"].values.astype(np.float64)
        nf = nf[np.isfinite(nf)]
        if len(nf) >= 10:
            med = np.median(nf)
            lo, hi = thresholds["n_genes_by_counts"]
            if lo > med * 0.80:
                thresholds["n_genes_by_counts"] = (med * 0.80, hi)
                log.warning(
                    "  SAFETY CAP: min_genes=%.0f exceeds 80%% of median (%.0f) — "
                    "clamping to %.0f. Over-filtering risk in low-gene-count tissue "
                    "(tissue_maturity=%s).",
                    lo,
                    med,
                    med * 0.80,
                    cfg.tissue_maturity,
                )
    return thresholds


# ══════════════════════════════════════════════════════════════════════════════
#  诊断图
# ══════════════════════════════════════════════════════════════════════════════


def _plot_qc_diagnostics(adata, thresholds, fig_dir, mode_label, cfg, log):
    """生成 3 张 QC 诊断图，标注当前使用的阈值线。

    参数:
        adata:      AnnData
        thresholds: _mad_thresholds() 或 _hard_thresholds() 返回的 dict
        fig_dir:    输出目录 (如 results/figures/02_qc)
        mode_label: "adaptive (MAD)" 或 "hard"
        cfg:        Config (用于 snRNA-seq 模式判断)
        log:        logger
    """
    os.makedirs(fig_dir, exist_ok=True)

    # ---- Panel A: nFeature 分布直方图 ----
    try:
        _fig, _ax = plt.subplots(figsize=cfg.plot.qc_figure_size)
        vals = adata.obs["n_genes_by_counts"].values
        vals = vals[np.isfinite(vals)]
        _ax.hist(vals, bins=100, color=cfg.plot.palette.qc_hist, edgecolor="white", alpha=0.85)
        lo, hi = thresholds["n_genes_by_counts"]
        if lo is not None:
            _ax.axvline(
                lo,
                color=cfg.plot.palette.qc_threshold,
                linestyle="--",
                linewidth=1.2,
                label=f"lo={lo:.0f}",
            )
        if hi is not None:
            _ax.axvline(
                hi,
                color=cfg.plot.palette.qc_threshold,
                linestyle="--",
                linewidth=1.2,
                label=f"hi={hi:.0f}",
            )
        _ax.set_xlabel("n_genes_by_counts (nFeature_RNA)")
        _ax.set_ylabel("Number of cells")
        _ax.set_title(
            f"nFeature distribution (N={adata.n_obs}, "
            f"median={np.median(vals):.0f}, N={len(vals)}, mode={mode_label})"
        )
        _ax.legend(fontsize=9)
        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, "nFeature_distribution.png"), dpi=cfg.plot.figure_dpi)
        plt.close(_fig)
        log.info("  Plot saved: nFeature_distribution.png")
    except Exception as e:
        log.warning("nFeature distribution plot failed: %s", e)

    # ---- Panel B: nCount vs nFeature 散点图 (按 mito% 着色) ----
    try:
        _fig, _ax = plt.subplots(
            figsize=(cfg.plot.qc_figure_size[0], cfg.plot.qc_figure_size[1] + 1)
        )
        x = adata.obs["total_counts"].values
        y = adata.obs["n_genes_by_counts"].values
        c = adata.obs["pct_counts_mt"].values
        # 过滤 NaN (某些细胞可能缺失 QC 指标)
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
        xp, yp, cp = x[finite], y[finite], c[finite]
        vmax = np.nanpercentile(c, 99) if np.isfinite(c).any() else None
        scat = _ax.scatter(
            xp,
            yp,
            c=cp,
            cmap=cfg.plot.palette.pseudotime,
            s=2,
            alpha=0.6,
            rasterized=True,
            vmax=vmax,
        )
        cbar = _fig.colorbar(scat, ax=_ax)
        cbar.set_label("% Mito")
        # 阈值线
        nfeat_lo, nfeat_hi = thresholds["n_genes_by_counts"]
        _, ncount_hi = thresholds["total_counts"]
        if nfeat_lo is not None:
            _ax.axhline(
                nfeat_lo, color=cfg.plot.palette.qc_threshold, linestyle="--", linewidth=1.0
            )
        if nfeat_hi is not None:
            _ax.axhline(
                nfeat_hi, color=cfg.plot.palette.qc_threshold, linestyle="--", linewidth=1.0
            )
        if ncount_hi is not None:
            _ax.axvline(
                ncount_hi,
                color="orange",
                linestyle="--",
                linewidth=1.0,
                label=f"nCount hi={ncount_hi:.0f}",
            )
        _ax.set_xlabel("total_counts (nCount_RNA)")
        _ax.set_ylabel("n_genes_by_counts (nFeature_RNA)")
        _ax.set_title(f"nCount vs nFeature (N={finite.sum()}/{len(x)}, mode={mode_label})")
        if ncount_hi is not None:
            _ax.legend(fontsize=9)
        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, "nCount_vs_nFeature.png"), dpi=cfg.plot.figure_dpi)
        plt.close(_fig)
        log.info("  Plot saved: nCount_vs_nFeature.png")
    except Exception as e:
        log.warning("nCount vs nFeature scatter plot failed: %s", e)

    # ---- Panel C: % Mito 分布直方图 ----
    try:
        _fig, _ax = plt.subplots(figsize=cfg.plot.qc_figure_size)
        vals = adata.obs["pct_counts_mt"].values
        vals = vals[np.isfinite(vals)]
        _ax.hist(vals, bins=100, color=cfg.plot.palette.qc_second, edgecolor="white", alpha=0.85)
        _, hi = thresholds["pct_counts_mt"]
        if hi is not None:
            _ax.axvline(
                hi,
                color=cfg.plot.palette.qc_threshold,
                linestyle="--",
                linewidth=1.2,
                label=f"hi={hi:.2f}%",
            )
        _ax.set_xlabel("pct_counts_mt (% Mito)")
        _ax.set_ylabel("Number of cells")
        _suffix = " (snRNA-seq: residual cytoplasm)" if cfg.qc.is_nuclei else ""
        _ax.set_title(
            f"% Mito distribution (N={len(vals)}, "
            f"median={np.median(vals):.2f}%, mode={mode_label}){_suffix}"
        )
        if hi is not None:
            _ax.legend(fontsize=9)
        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, "pct_mito_distribution.png"), dpi=cfg.plot.figure_dpi)
        plt.close(_fig)
        log.info("  Plot saved: pct_mito_distribution.png")
    except Exception as e:
        log.warning("pct_mito distribution plot failed: %s", e)


def _plot_nfeature_kde(adata, fig_dir, mode_label, cfg, log):
    """Panel D: nFeature KDE density with peak markers."""
    try:
        vals = adata.obs["n_genes_by_counts"].values.astype(np.float64)
        n_peaks, peaks_x, peaks_y, x_range, density, is_multimodal = _detect_peaks(vals)
        if len(x_range) == 0:
            log.info("  KDE plot: skipped (too few cells)")
            return

        _fig, _ax = plt.subplots(figsize=cfg.plot.qc_figure_size)
        _ax.plot(x_range, density, color=cfg.plot.palette.qc_hist, linewidth=1.5)
        _ax.fill_between(x_range, density, alpha=0.15, color=cfg.plot.palette.qc_hist)

        y_offset = density.max() * 0.03
        for px, py in zip(peaks_x, peaks_y):
            _ax.scatter(
                px,
                py + y_offset,
                marker="^",
                s=80,
                color=cfg.plot.palette.qc_third,
                edgecolors="black",
                linewidths=0.5,
                zorder=5,
                label=f"Peak at {px:.0f} genes",
            )

        # Threshold lines
        _ax.axvline(
            cfg.qc.min_genes,
            color=cfg.plot.palette.qc_threshold,
            linestyle="--",
            linewidth=1.0,
            label=f"lo={cfg.qc.min_genes:.0f}",
        )
        _ax.axvline(
            cfg.qc.max_genes,
            color=cfg.plot.palette.qc_threshold,
            linestyle="--",
            linewidth=1.0,
            label=f"hi={cfg.qc.max_genes:.0f}",
        )

        assessment = "BIMODAL" if is_multimodal else ("MULTIPEAK" if n_peaks >= 2 else "UNIMODAL")
        assess_color = "#c0392b" if is_multimodal else ("#e67e22" if n_peaks >= 2 else "#27ae60")

        _ax.set_xlabel("n_genes_by_counts (nFeature_RNA)")
        _ax.set_ylabel("Density")
        _peak_label = "s" if n_peaks != 1 else ""
        _ax.set_title(f"nFeature KDE density ({n_peaks} peak{_peak_label}, mode={mode_label})")
        _ax.legend(fontsize=8, loc="upper right")

        _ax.annotate(
            f"{assessment}",
            xy=(0.02, 0.08),
            xycoords="axes fraction",
            fontsize=12,
            fontweight="bold",
            color=assess_color,
            ha="left",
            va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor=cfg.plot.palette.grn_facecolor, alpha=0.8
            ),
        )

        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, "nFeature_KDE_density.png"), dpi=cfg.plot.figure_dpi)
        plt.close(_fig)
        log.info("  Plot saved: nFeature_KDE_density.png (%s)", assessment)
    except Exception as e:
        log.warning("nFeature KDE plot failed: %s", e)
        log.warning("pct_mito distribution plot failed: %s", e)


# ══════════════════════════════════════════════════════════════════════════════
#  QC 指标计算
# ══════════════════════════════════════════════════════════════════════════════


def compute_qc_metrics(adata, cfg, log):
    log.info("Computing QC metrics...")

    # Auto-detect mt_gene_pattern for non-human species
    mt_candidates = ["MT-", "mt-", "Mt-"]
    mt_pattern = cfg.qc.mt_gene_pattern
    if not any(adata.var_names.str.startswith(mt_pattern)):
        for alt in mt_candidates:
            if any(adata.var_names.str.startswith(alt)):
                log.info("Auto-switched mt_gene_pattern: '%s' -> '%s'", mt_pattern, alt)
                cfg.qc.mt_gene_pattern = alt
                break
        else:  # no candidate matched
            if not cfg.qc.mt_gene_list:
                log.warning(
                    "No MT genes detected (%s candidates failed: %s)",
                    mt_pattern,
                    mt_candidates,
                )

    mt_mask = adata.var_names.str.startswith(cfg.qc.mt_gene_pattern)
    if cfg.qc.mt_gene_list:
        mt_mask = mt_mask | adata.var_names.isin(cfg.qc.mt_gene_list)
    adata.var["mt"] = mt_mask
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # ── 流式 QC metrics（复刻 scanpy describe_obs + describe_var 语义）──
    from scanpy.pp._qc import top_segment_proportions_sparse_csr

    x_mat = adata.X
    n_cells, n_genes_t = x_mat.shape
    block_size = 200_000
    mt_idx = np.where(mt_mask)[0]
    ribo_idx = np.where(adata.var["ribo"].values)[0]
    x_dtype = x_mat.dtype
    # getnnz dtype 跟随磁盘 indptr dtype（全量内存读行为），backed 块读会 downcast int32
    obs_n_genes = np.empty(n_cells, dtype=x_mat._indptr.dtype)
    obs_total = np.empty(n_cells, dtype=x_dtype)
    obs_top20 = np.empty(n_cells, dtype=np.float64)
    obs_mt = np.empty(n_cells, dtype=x_dtype)
    obs_ribo = np.empty(n_cells, dtype=x_dtype)
    var_nnz = np.zeros(n_genes_t, dtype=np.int64)
    var_sum = np.zeros(n_genes_t, dtype=x_dtype)
    ns = np.array([20], dtype=np.int64)
    for i in range(0, n_cells, block_size):
        xb = x_mat[i : i + block_size]
        obs_n_genes[i : i + block_size] = xb.getnnz(axis=1)
        s = np.asarray(xb.sum(axis=1)).ravel()
        obs_total[i : i + block_size] = s
        props = top_segment_proportions_sparse_csr(xb.data, xb.indptr, ns)
        obs_top20[i : i + block_size] = props[:, 0] * 100
        obs_mt[i : i + block_size] = np.asarray(xb[:, mt_idx].sum(axis=1)).ravel()
        obs_ribo[i : i + block_size] = np.asarray(xb[:, ribo_idx].sum(axis=1)).ravel()
        var_nnz += xb.getnnz(axis=0)
        var_sum += np.asarray(xb.sum(axis=0)).ravel()
    obs = adata.obs
    obs["n_genes_by_counts"] = obs_n_genes
    obs["log1p_n_genes_by_counts"] = np.log1p(obs_n_genes)
    obs["total_counts"] = obs_total
    obs["log1p_total_counts"] = np.log1p(obs_total)
    obs["pct_counts_in_top_20_genes"] = obs_top20
    obs["total_counts_mt"] = obs_mt
    obs["log1p_total_counts_mt"] = np.log1p(obs_mt)
    obs["pct_counts_mt"] = obs_mt / obs_total * 100
    obs["total_counts_ribo"] = obs_ribo
    obs["log1p_total_counts_ribo"] = np.log1p(obs_ribo)
    obs["pct_counts_ribo"] = obs_ribo / obs_total * 100
    var = adata.var
    var_mean = var_sum / n_cells
    var["n_cells_by_counts"] = var_nnz
    var["mean_counts"] = var_mean
    var["log1p_mean_counts"] = np.log1p(var_mean)
    var["pct_dropout_by_counts"] = (1 - var_nnz / n_cells) * 100
    var["total_counts"] = var_sum
    var["log1p_total_counts"] = np.log1p(var_sum)
    adata.obs["log_genes_per_umi"] = (
        np.log10(adata.obs["n_genes_by_counts"]) / np.log10(adata.obs["total_counts"])
    ).replace([np.inf, -np.inf], np.nan)
    log.info("  Median genes/cell: %.0f", adata.obs["n_genes_by_counts"].median())
    log.info("  Median UMIs/cell: %.0f", adata.obs["total_counts"].median())
    log.info("  Median mito%%:    %.2f%%", adata.obs["pct_counts_mt"].median())
    if cfg.qc.is_nuclei:
        log.info(
            "  [snRNA-seq mode] Mitochondrial reads reflect cytoplasmic residue, not cell stress."
        )
    log.info("  Median complexity: %.3f", adata.obs["log_genes_per_umi"].median())


# ══════════════════════════════════════════════════════════════════════════════
#  过滤
# ══════════════════════════════════════════════════════════════════════════════


def _nonzero_col_counts(x_mat, mask=None, block_size=500_000):
    """每基因非零计数（块式，与 scanpy `(X > 0).sum(0)` 语义一致：NaN 不计）。
    mask 传入时仅在保留行上计数（等价于先过滤行再数）。"""
    n, m = x_mat.shape
    counts = np.zeros(m, dtype=np.int64)
    for i in range(0, n, block_size):
        xr = x_mat[i : i + block_size]
        if mask is not None:
            mm = np.asarray(mask[i : i + block_size])
            xr = xr[mm]
        if xr.nnz == 0:
            continue
        keep = xr.data > 0
        counts += np.bincount(xr.indices[keep], minlength=m)
    return counts


def _write_qc_h5ad(adata, mask_obs, vmask, n_cells_counts, cfg, log):
    """Filter-on-write: 过滤推迟到写入阶段，逐块行+列过滤直接写 h5ad。

    内存峰值 = X 全量（加载必需）+ 单块临时（~O(块)），零矩阵移动。
    压缩解析与 safe_write 一致（per_step_h5ad_compression > cfg.h5ad_compression > gzip），
    先写同目录隐藏 tmp 再原子 os.replace。"""
    import h5py
    from anndata._io.h5ad import write_elem

    x_mat = adata.X
    n, m = x_mat.shape
    n_keep_o = int(mask_obs.sum())
    n_keep_v = int(vmask.sum())
    vmask_np = np.asarray(vmask)
    block_size = 200_000

    # 02 专属默认：zstd 全面优于 gzip1（实测 GSE137398：写 11.8→5.9s, 读 3.2→1.6s, 文件 251→217MB）
    # hdf5plugin 不可用时回退 gzip；FUXI_QC_COMPR 环境变量可覆盖（如 gzip/lzf）
    compression = os.environ.get("FUXI_QC_COMPR", "")
    if not compression:
        compression = getattr(cfg, "h5ad_compression", "gzip")
        if compression == "gzip":
            compression = "zstd" if hdf5plugin is not None else "gzip"
    kwargs = {"compression": compression}
    if compression == "gzip":
        kwargs["compression_opts"] = 1
    elif compression == "zstd":
        if hdf5plugin is None:
            kwargs = {"compression": "gzip", "compression_opts": 1}
        else:
            # 实测（GSE137398 76k cells, 335.9M nnz）:
            #   gzip1 写 11.8s/251MB 读 3.2s, zstd1 写 5.9s/217MB 读 1.6s
            #   —— zstd 写读均快 ~2x 且文件小 13%，全面胜出
            kwargs = dict(hdf5plugin.Zstd(clevel=1))

    target = os.environ.get("FUXI_QC_OUT", cfg.qc_h5ad)
    target_dir = os.path.dirname(target) or "."
    os.makedirs(target_dir, exist_ok=True)
    tmp_path = os.path.join(target_dir, f".{os.path.basename(target)}.tmp.{os.getpid()}")

    t0 = time.time()
    with h5py.File(tmp_path, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        for key in ["layers", "obsm", "obsp", "varm", "varp"]:
            f.create_group(key)
        f.create_group("uns")
        obs = adata.obs.loc[mask_obs].copy()
        var = adata.var.loc[vmask].copy()
        var["n_cells"] = n_cells_counts[vmask].astype(np.int64)
        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        xg = f.create_group("X")
        xg.attrs["encoding-type"] = "csr_matrix"
        xg.attrs["encoding-version"] = "0.1.0"
        xg.attrs["shape"] = (n_keep_o, n_keep_v)
        d_data = xg.create_dataset(
            "data", (0,), maxshape=(None,), dtype=x_mat.dtype, chunks=(65536,), **kwargs
        )
        d_idx = xg.create_dataset(
            "indices", (0,), maxshape=(None,), dtype=np.int64, chunks=(65536,), **kwargs
        )
        d_iptr = xg.create_dataset(
            "indptr", (0,), maxshape=(None,), dtype=np.int64, chunks=(4096,), **kwargs
        )
        from concurrent.futures import ThreadPoolExecutor

        new_indptr = np.zeros(n_keep_o + 1, dtype=np.int64)
        prefetch = 2
        starts = list(range(0, n, block_size))
        pos = 0
        w = 0

        def _load(st):
            mi = np.asarray(mask_obs[st : st + block_size])
            if not mi.any():
                return st, None
            xr = x_mat[st : st + block_size][mi][:, vmask_np]
            if xr.nnz == 0:
                return st, None
            return st, xr

        with ThreadPoolExecutor(max_workers=2) as ex:
            pending = {}
            for st in starts[:prefetch]:
                pending[st] = ex.submit(_load, st)
            for st in starts:
                nxt = st + prefetch * block_size
                if nxt < n:
                    pending[nxt] = ex.submit(_load, nxt)
                fut = pending.pop(st)
                _, xr = fut.result()
                if xr is None:
                    continue
                k = xr.shape[0]
                new_indptr[pos + 1 : pos + k + 1] = new_indptr[pos] + xr.indptr[1:]
                pos += k
                n_w = xr.nnz
                d_data.resize(w + n_w, axis=0)
                d_idx.resize(w + n_w, axis=0)
                d_data[w : w + n_w] = xr.data
                d_idx[w : w + n_w] = xr.indices
                w += n_w
                del xr
        d_iptr.resize(n_keep_o + 1, axis=0)
        d_iptr[:] = new_indptr
    os.replace(tmp_path, target)
    log.info(
        "  Saved %s (%.1f MB, filter-on-write %.1fs)",
        os.path.basename(target),
        os.path.getsize(target) / 1e6,
        time.time() - t0,
    )


def filter_cells(adata, thresholds, cfg, log):
    """根据阈值字典过滤细胞。

    参数:
        adata:      AnnData
        thresholds: _mad_thresholds() 或 _hard_thresholds() 返回的 dict
        cfg:        Config
        log:        logger

    返回:
        np.ndarray (bool): 保留细胞的布尔 mask（行过滤推迟到写入阶段执行）
    """
    n_before = adata.n_obs

    # ---- doublet ----
    if cfg.scrublet.run:
        f_doublet = adata.obs["predicted_doublet"]
        n_doublet = f_doublet.sum()
        log.info(
            "Doublet filtering: removing %d predicted doublets (%.1f%%)",
            n_doublet,
            100 * n_doublet / n_before if n_before else 0,
        )

        # ── Doublet guard: warn if Scrublet ran but predicted zero doublets ──
        _n_doublet = f_doublet.sum()
        if _n_doublet == 0:
            log.warning(
                "Doublet detection enabled but 0 doublets predicted. "
                "Likely cause: input not raw_counts (Scrublet was skipped) or "
                "expected_doublet_rate too low. No doublets removed; "
                "interpret downstream with caution.",
            )
    else:
        log.info("Doublet filtering: skipped (scrublet disabled).")
        f_doublet = pd.Series(False, index=adata.obs_names)

    # ---- 从 thresholds 读取各指标边界 ----
    gf_lo, gf_hi = thresholds["n_genes_by_counts"]
    _, tc_hi = thresholds["total_counts"]
    _, mt_hi = thresholds["pct_counts_mt"]
    cpx_lo, _ = thresholds["log_genes_per_umi"]

    # ---- 构建过滤条件 ----
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

    mask = ~f_doublet & ~f_any
    n_after = int(mask.sum())
    log.info("  After QC filtering: %d cells", n_after)
    return mask


# ══════════════════════════════════════════════════════════════════════════════
#  自适应阈值建议（在过滤之前，基于数据分布输出诊断意见）
# ══════════════════════════════════════════════════════════════════════════════


def _suggest_qc_thresholds(adata, cfg, log):
    """检查当前配置阈值 vs 数据分布，输出诊断建议。

    在 compute_qc_metrics 之后、filter_cells 之前调用。
    零副作用 —— 仅输出日志，不改变任何阈值。
    """
    nf = adata.obs["n_genes_by_counts"].values.astype(np.float64)
    nf = nf[np.isfinite(nf)]
    if len(nf) < 10:
        log.info("  [qc-suggest] Too few cells (<10) for reliable distribution diagnosis.")
        return

    med = np.median(nf)
    mad = median_abs_deviation(nf, scale="normal")
    p10 = np.percentile(nf, 10)
    p90 = np.percentile(nf, 90)
    n_peaks, _, _, _, _, is_multimodal = _detect_peaks(nf)

    log.info(
        "  [qc-suggest] nFeature distribution: median=%.0f, MAD=%.0f, "
        "P10=%.0f, P90=%.0f, peaks=%d%s",
        med,
        mad,
        p10,
        p90,
        n_peaks,
        " (MULTIMODAL)" if is_multimodal else "",
    )

    # ── 建议 1: min_genes 是否合理 ──
    if cfg.qc.min_genes > med * 0.80:
        suggested = max(int(med * 0.80), 100)
        log.warning(
            "  [qc-suggest] ⚠️  OVER-FILTERING RISK: min_genes=%d > 80%% of median (%.0f). "
            "Will remove %.0f%% of cells. Suggested min_genes=%.0f. "
            "(Fix: set qc.min_genes to %.0f in config, or set qc.use_adaptive_thresholds=true)",
            cfg.qc.min_genes,
            med,
            100 * (nf < cfg.qc.min_genes).sum() / len(nf),
            suggested,
            suggested,
        )
    else:
        removed_by_lo = (nf < cfg.qc.min_genes).sum()
        pct = 100 * removed_by_lo / len(nf)
        if pct < 1:
            log.info(
                "  [qc-suggest] Only %.1f%% cells below min_genes=%d — "
                "check if QC is too permissive (consider qc.use_adaptive_thresholds=true).",
                pct,
                cfg.qc.min_genes,
            )
        else:
            log.info(
                "  [qc-suggest] min_genes=%d will remove %.1f%% of cells — OK.",
                cfg.qc.min_genes,
                pct,
            )

    # ── 建议 2: max_genes 覆盖范围 ──
    if cfg.qc.max_genes < p90:
        pct_above = 100 * (nf > cfg.qc.max_genes).sum() / len(nf)
        log.warning(
            "  [qc-suggest] ⚠️  max_genes=%.0f < P90 (%.0f) — will remove %.1f%% of cells. "
            "Consider max_genes=%.0f.",
            cfg.qc.max_genes,
            p90,
            pct_above,
            int(p90 * 1.2),
        )

    # ── 建议 3: 线粒体比例异常 ──
    mt = adata.obs["pct_counts_mt"].values
    mt = mt[np.isfinite(mt)]
    mt_med = np.median(mt)
    mt_p90 = np.percentile(mt, 90)

    if mt_med > 10:
        log.warning(
            "  [qc-suggest] ⚠️  Median pct_mito=%.1f%% is high (P90=%.1f%%) — "
            "check for stressed/dead cell population or tissue-specific biology. "
            "Consider qc.use_adaptive_thresholds=true for MAD-based mito threshold.",
            mt_med,
            mt_p90,
        )
    if mt_med == 0:
        log.info(
            "  [qc-suggest] Median pct_mito=0%% — may indicate incorrect mt_gene_pattern "
            "(current: '%s') for this species.",
            cfg.qc.mt_gene_pattern,
        )

    # ── 建议 4: MAD 推荐值 vs 当前配置 ──
    recommended_lo = max(med - 3 * mad, 100)
    med + 3 * mad
    if abs(cfg.qc.min_genes - recommended_lo) > max(recommended_lo * 0.5, 100):
        log.info(
            "  [qc-suggest] 💡 MAD-based suggestion: min_genes=%.0f (config has %d, "
            "diff=%.0f). Consider adjusting config or enabling use_adaptive_thresholds.",
            recommended_lo,
            cfg.qc.min_genes,
            abs(cfg.qc.min_genes - recommended_lo),
        )

    # ── 建议 5: 发育组织诊断 ──
    if cfg.tissue_maturity != "developing" and med < 800:
        log.info(
            "  [qc-suggest] Low median gene count (%.0f < 800) — "
            "consider setting tissue_maturity: 'developing' in config for wider MAD bounds.",
            med,
        )

    # ── 建议 6: max_pct_mito 合理性 ──
    pct_above_mito = 100 * (mt > cfg.qc.max_pct_mito).sum() / len(mt) if len(mt) > 0 else 0
    if pct_above_mito > 15:
        log.warning(
            "  [qc-suggest] ⚠️  max_pct_mito=%.0f%% will remove %.1f%% of cells. "
            "Consider raising max_pct_mito or checking tissue-specific biology.",
            cfg.qc.max_pct_mito,
            pct_above_mito,
        )
    elif pct_above_mito < 0.5:
        log.info(
            "  [qc-suggest] Only %.1f%% cells above max_pct_mito=%.0f%% — mito filter is nearly inactive.",
            pct_above_mito,
            cfg.qc.max_pct_mito,
        )


# ══════════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════════


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("02_qc", os.path.join(cfg.log_dir, "02_qc.log"))
    log.info("Step 02: QC filtering (doublets already removed in Step 01)")

    input_path = os.path.join(cfg.h5ad_dir, "01_doublet.h5ad")
    adata = sc.read_h5ad(input_path, backed="r")
    log.info("Loaded: %s — %d cells × %d genes", input_path, adata.n_obs, adata.n_vars)

    # 1. 计算 QC 指标
    compute_qc_metrics(adata, cfg, log)
    # 1a. 自适应阈值建议（诊断模式，零副作用）
    _suggest_qc_thresholds(adata, cfg, log)

    # 2. 确定阈值 (MAD 或硬阈值)
    if cfg.qc.use_adaptive_thresholds:
        log.info(
            "Mode: adaptive (MAD × %.1f, nCount MAD × %.1f)",
            cfg.qc.mad_n_mads,
            cfg.qc.ncount_max_mad,
        )
        thresholds = _mad_thresholds(adata, cfg, log)
        mode_label = "adaptive (MAD)"
    else:
        log.info("Mode: hard thresholds")
        thresholds = _hard_thresholds(cfg, log, adata=adata)
        mode_label = "hard"

    # 3. 生成诊断图 (在任何过滤之前，展示原始分布 + 阈值线)
    fig_dir = os.path.join(cfg.figure_dir, "02_qc")
    _plot_qc_diagnostics(adata, thresholds, fig_dir, mode_label, cfg, log)

    # 3a. 生成 nFeature KDE 密度峰图
    _fig_dir = os.path.join(cfg.figure_dir, "02_qc")
    _plot_nfeature_kde(adata, _fig_dir, mode_label, cfg, log)

    n_before = adata.n_obs
    mask_obs = filter_cells(adata, thresholds, cfg, log)
    min_cells = cfg.qc.min_cells_per_gene
    n_cells_counts = _nonzero_col_counts(adata.X, mask_obs)
    vmask = n_cells_counts >= min_cells
    log.info("After gene filtering: %d genes", int(vmask.sum()))
    # ── Checkpoint: filter-on-write（过滤与写入融合，零矩阵拷贝）──
    _write_qc_h5ad(adata, mask_obs, vmask, n_cells_counts, cfg, log)
    n_after = int(mask_obs.sum())
    # 5. QC SUMMARY — 事后评估过滤率
    n_after = adata.n_obs
    pct_removed = 100 * (n_before - n_after) / n_before if n_before else 0
    if pct_removed > 30:
        log.warning(
            "QC SUMMARY: %d → %d cells (%.1f%% removed) — ⚠️  high filtering rate, review thresholds",
            n_before,
            n_after,
            pct_removed,
        )
    elif pct_removed < 1:
        log.info(
            "QC SUMMARY: %d → %d cells (%.1f%% removed) — QC may be too permissive",
            n_before,
            n_after,
            pct_removed,
        )
    else:
        log.info(
            "QC SUMMARY: %d → %d cells (%.1f%% removed) — OK",
            n_before,
            n_after,
            pct_removed,
        )

    log.info("Step 02 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
