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
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation, gaussian_kde
from scipy.signal import find_peaks



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
        return n_peaks, x_range[peaks].tolist(), density[peaks].tolist(), x_range, density, is_multimodal
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
    _nf = adata.obs['n_genes_by_counts'].values.astype(np.float64)
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
        log.warning("Falling back to hard thresholds — median too low for MAD normality assumption.")
        return _hard_thresholds(cfg, log)


    # ---- n_genes_by_counts (nFeature_RNA) ----
    vals = adata.obs['n_genes_by_counts'].values.astype(np.float64)
    med = np.median(vals)
    mad = median_abs_deviation(vals, scale='normal')
    _maturity_mad = {"developing": 5.0}.get(cfg.tissue_maturity, cfg.mad_n_mads)
    lo_mad = max(med - _maturity_mad * mad, 0)
    hi_mad = med + _maturity_mad * mad
    # 硬阈值做地板/天花板
    lo = max(lo_mad, cfg.min_genes)
    # Safety clamp: MAD 下界不应超过用户硬阈值（参考 ddqc, Subramanian 2022）
    if lo > cfg.min_genes:
        lo_orig = lo
        lo = cfg.min_genes
        log.warning(
            '  MAD lower bound (%.0f) exceeds min_genes (%.0f) — clamping to min_genes. MAD is unreliable for this distribution; consider use_adaptive_thresholds=False.',
            lo_orig, cfg.min_genes,
        )
    _safe_floor = cfg.min_mad_upper_genes_nuclei if cfg.is_nuclei else cfg.min_mad_upper_genes
    hi = min(max(hi_mad, _safe_floor), cfg.max_genes)
    thresholds['n_genes_by_counts'] = (lo, hi)
    log.info("  n_genes_by_counts: median=%.0f, MAD=%.0f  →  (lo=%.0f, hi=%.0f)  [adaptive ×%.1f, %s]",
             med, mad, lo, hi, _maturity_mad, cfg.tissue_maturity)

    # ── 分布检测 ──
    try:
        vals_nf = adata.obs['n_genes_by_counts'].values.astype(np.float64)
        n_peaks, _, _, _, _, is_multimodal = _detect_peaks(vals_nf)
        # Also report P10/P90 for context
        p10 = np.percentile(vals_nf[np.isfinite(vals_nf)], 10)
        p90 = np.percentile(vals_nf[np.isfinite(vals_nf)], 90)
        log.info("  nFeature distribution: %d KDE peak%s (P10=%.0f, P90=%.0f)%s",
                 n_peaks, "" if n_peaks == 1 else "s", p10, p90,
                 " — MULTIMODAL: MAD may over-filter the high-expression population" if is_multimodal else "")
    except Exception:
        pass

    # ---- total_counts (nCount_RNA) ----
    # 非 raw_counts 数据不设 total_counts 上限
    if cfg.expression_type == "raw_counts":
        vals = adata.obs['total_counts'].values.astype(np.float64)
        med = np.median(vals)
        mad = median_abs_deviation(vals, scale='normal')
        hi = med + cfg.qc_ncount_max_mad * mad
        thresholds['total_counts'] = (None, hi)
        log.info("  total_counts:       median=%.0f, MAD=%.0f  →  hi=%.0f  [adaptive ×%.1f]",
                 med, mad, hi, cfg.qc_ncount_max_mad)
    else:
        thresholds['total_counts'] = (None, None)
        log.info("  total_counts:       (skipped — expression_type=%s)", cfg.expression_type)

    # ---- pct_counts_mt ----
    if cfg.is_nuclei:
        hi = cfg.max_pct_mito_nuclei
        log.info("  pct_counts_mt:      hi=%.2f%%  [snRNA-seq: fixed threshold, MAD skipped]", hi)
    else:
        vals = adata.obs['pct_counts_mt'].values.astype(np.float64)
        med = np.median(vals)
        mad = median_abs_deviation(vals, scale='normal')
        hi_mad = med + cfg.mad_n_mads * mad
        hi = min(hi_mad, cfg.max_pct_mito)
        log.info("  pct_counts_mt:      median=%.2f%%, MAD=%.2f%% ->  hi=%.2f%%  [adaptive, factor=%.1f]",
                 med, mad, hi, cfg.mad_n_mads)
    thresholds['pct_counts_mt'] = (None, hi)

    # ---- log_genes_per_umi (complexity) ----
    # 非 raw_counts 数据下复杂度指标无解释力，跳过
    if cfg.expression_type == "raw_counts":
        vals = adata.obs['log_genes_per_umi'].values.astype(np.float64)
        finite = vals[np.isfinite(vals)]
        med = np.median(finite)
        mad = median_abs_deviation(finite, scale='normal')
        lo_mad = max(med - _maturity_mad * mad, 0)
        lo = max(lo_mad, cfg.min_genes_per_umi)
        # Safety clamp: MAD 下界不应超过用户硬阈值
        if lo > cfg.min_genes_per_umi:
            lo_orig = lo
            lo = cfg.min_genes_per_umi
            log.warning(
                '  MAD lower bound (%.4f) exceeds min_genes_per_umi (%.4f) — clamping.',
                lo_orig, cfg.min_genes_per_umi,
            )
        thresholds['log_genes_per_umi'] = (lo, None)
        log.info("  log_genes_per_umi:  median=%.4f, MAD=%.4f →  lo=%.4f  [adaptive ×%.1f, %s]",
                 med, mad, lo, _maturity_mad, cfg.tissue_maturity)
    else:
        thresholds['log_genes_per_umi'] = (None, None)
        log.info("  log_genes_per_umi:  (skipped — expression_type=%s)", cfg.expression_type)

    return thresholds


def _hard_thresholds(cfg, log):
    """从 Config 构建硬阈值字典（现有行为）。"""
    is_native = (cfg.expression_type == "raw_counts")
    thresholds = {
        'n_genes_by_counts': (cfg.min_genes, cfg.max_genes),
        'total_counts':      (None, None),                     # raw_counts 下也未启用硬上限
        'pct_counts_mt':     (None, cfg.max_pct_mito_nuclei if cfg.is_nuclei else cfg.max_pct_mito),
        'log_genes_per_umi': (cfg.min_genes_per_umi, None) if is_native else (None, None),
    }
    log.info("  Using hard thresholds from config:")
    log.info("    n_genes_by_counts:  lo=%d, hi=%d", cfg.min_genes, cfg.max_genes)
    log.info("    total_counts:       (no limit)")
    if is_native:
        log.info("    pct_counts_mt:      hi=%.1f%%", cfg.max_pct_mito)
        log.info("    log_genes_per_umi:  lo=%.4f", cfg.min_genes_per_umi)
    else:
        log.info("    pct_counts_mt:      hi=%.1f%%", cfg.max_pct_mito)
        log.info("    log_genes_per_umi:  (skipped — expression_type=%s)", cfg.expression_type)
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
        _fig, _ax = plt.subplots(figsize=(8, 5))
        vals = adata.obs['n_genes_by_counts'].values
        vals = vals[np.isfinite(vals)]
        _ax.hist(vals, bins=100, color='steelblue', edgecolor='white', alpha=0.85)
        lo, hi = thresholds['n_genes_by_counts']
        if lo is not None:
            _ax.axvline(lo, color='red', linestyle='--', linewidth=1.2,
                       label=f'lo={lo:.0f}')
        if hi is not None:
            _ax.axvline(hi, color='red', linestyle='--', linewidth=1.2,
                       label=f'hi={hi:.0f}')
        _ax.set_xlabel('n_genes_by_counts (nFeature_RNA)')
        _ax.set_ylabel('Number of cells')
        _ax.set_title(f'nFeature distribution (N={adata.n_obs}, '
                     f'median={np.median(vals):.0f}, N={len(vals)}, mode={mode_label})')
        _ax.legend(fontsize=9)
        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, 'nFeature_distribution.png'), dpi=150)
        plt.close(_fig)
        log.info("  Plot saved: nFeature_distribution.png")
    except Exception as e:
        log.warning("nFeature distribution plot failed: %s", e)

    # ---- Panel B: nCount vs nFeature 散点图 (按 mito% 着色) ----
    try:
        _fig, _ax = plt.subplots(figsize=(8, 6))
        x = adata.obs['total_counts'].values
        y = adata.obs['n_genes_by_counts'].values
        c = adata.obs['pct_counts_mt'].values
        # 过滤 NaN (某些细胞可能缺失 QC 指标)
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
        xp, yp, cp = x[finite], y[finite], c[finite]
        vmax = np.nanpercentile(c, 99) if np.isfinite(c).any() else None
        scat = _ax.scatter(xp, yp, c=cp, cmap='viridis', s=2, alpha=0.6,
                          rasterized=True, vmax=vmax)
        cbar = _fig.colorbar(scat, ax=_ax)
        cbar.set_label('% Mito')
        # 阈值线
        nfeat_lo, nfeat_hi = thresholds['n_genes_by_counts']
        _, ncount_hi = thresholds['total_counts']
        if nfeat_lo is not None:
            _ax.axhline(nfeat_lo, color='red', linestyle='--', linewidth=1.0)
        if nfeat_hi is not None:
            _ax.axhline(nfeat_hi, color='red', linestyle='--', linewidth=1.0)
        if ncount_hi is not None:
            _ax.axvline(ncount_hi, color='orange', linestyle='--', linewidth=1.0,
                       label=f'nCount hi={ncount_hi:.0f}')
        _ax.set_xlabel('total_counts (nCount_RNA)')
        _ax.set_ylabel('n_genes_by_counts (nFeature_RNA)')
        _ax.set_title(f'nCount vs nFeature (N={finite.sum()}/{len(x)}, mode={mode_label})')
        if ncount_hi is not None:
            _ax.legend(fontsize=9)
        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, 'nCount_vs_nFeature.png'), dpi=150)
        plt.close(_fig)
        log.info("  Plot saved: nCount_vs_nFeature.png")
    except Exception as e:
        log.warning("nCount vs nFeature scatter plot failed: %s", e)

    # ---- Panel C: % Mito 分布直方图 ----
    try:
        _fig, _ax = plt.subplots(figsize=(8, 5))
        vals = adata.obs['pct_counts_mt'].values
        vals = vals[np.isfinite(vals)]
        _ax.hist(vals, bins=100, color='indianred', edgecolor='white', alpha=0.85)
        _, hi = thresholds['pct_counts_mt']
        if hi is not None:
            _ax.axvline(hi, color='red', linestyle='--', linewidth=1.2,
                       label=f'hi={hi:.2f}%')
        _ax.set_xlabel('pct_counts_mt (% Mito)')
        _ax.set_ylabel('Number of cells')
        _suffix = " (snRNA-seq: residual cytoplasm)" if cfg.is_nuclei else ""
        _ax.set_title(f'% Mito distribution (N={len(vals)}, '
                     f'median={np.median(vals):.2f}%, mode={mode_label}){_suffix}')
        if hi is not None:
            _ax.legend(fontsize=9)
        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, 'pct_mito_distribution.png'), dpi=150)
        plt.close(_fig)
        log.info("  Plot saved: pct_mito_distribution.png")
    except Exception as e:
        log.warning("pct_mito distribution plot failed: %s", e)
def _plot_nfeature_kde(adata, fig_dir, mode_label, cfg, log):
    """Panel D: nFeature KDE density with peak markers."""
    try:
        vals = adata.obs['n_genes_by_counts'].values.astype(np.float64)
        n_peaks, peaks_x, peaks_y, x_range, density, is_multimodal = _detect_peaks(vals)
        if len(x_range) == 0:
            log.info("  KDE plot: skipped (too few cells)")
            return

        _fig, _ax = plt.subplots(figsize=(8, 5))
        _ax.plot(x_range, density, color='steelblue', linewidth=1.5)
        _ax.fill_between(x_range, density, alpha=0.15, color='steelblue')

        for px, py in zip(peaks_x, peaks_y):
            _ax.scatter(px, py, marker='^', s=80, color='darkorange',
                       edgecolors='black', linewidths=0.5, zorder=5,
                       label=f'Peak at {px:.0f} genes')

        # Threshold lines
        _ax.axvline(cfg.min_genes, color='red', linestyle='--', linewidth=1.0,
                   label=f'lo={cfg.min_genes:.0f}')
        _ax.axvline(cfg.max_genes, color='red', linestyle='--', linewidth=1.0,
                   label=f'hi={cfg.max_genes:.0f}')

        assessment = 'BIMODAL' if is_multimodal else ('MULTIPEAK' if n_peaks >= 2 else 'UNIMODAL')
        assess_color = '#c0392b' if is_multimodal else ('#e67e22' if n_peaks >= 2 else '#27ae60')

        _ax.set_xlabel('n_genes_by_counts (nFeature_RNA)')
        _ax.set_ylabel('Density')
        _peak_label = "s" if n_peaks != 1 else ""
        _ax.set_title(f'nFeature KDE density ({n_peaks} peak{_peak_label}, mode={mode_label})')
        _ax.legend(fontsize=8, loc='upper right')

        _ax.annotate(f'{assessment}', xy=(0.02, 0.08), xycoords='axes fraction',
                    fontsize=12, fontweight='bold', color=assess_color,
                    ha='left', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        _fig.tight_layout()
        _fig.savefig(os.path.join(fig_dir, 'nFeature_KDE_density.png'), dpi=150)
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
    if cfg.is_nuclei:
        log.info("  [snRNA-seq mode] Mitochondrial reads reflect cytoplasmic residue, not cell stress.")
    log.info("  Median complexity: %.3f", adata.obs['log_genes_per_umi'].median())


# ══════════════════════════════════════════════════════════════════════════════
#  过滤
# ══════════════════════════════════════════════════════════════════════════════

def filter_cells(adata, thresholds, cfg, log):
    """根据阈值字典过滤细胞。

    参数:
        adata:      AnnData
        thresholds: _mad_thresholds() 或 _hard_thresholds() 返回的 dict
        cfg:        Config
        log:        logger

    返回:
        AnnData (过滤后)
    """
    n_before = adata.n_obs

    # ---- doublet ----
    f_doublet = adata.obs['predicted_doublet']
    n_doublet = f_doublet.sum()
    log.info("Doublet filtering: removing %d predicted doublets (%.1f%%)",
             n_doublet, 100 * n_doublet / n_before if n_before else 0)

    # ---- 从 thresholds 读取各指标边界 ----
    gf_lo, gf_hi = thresholds['n_genes_by_counts']
    _,      tc_hi = thresholds['total_counts']
    _,      mt_hi = thresholds['pct_counts_mt']
    cpx_lo, _     = thresholds['log_genes_per_umi']

    # ---- 构建过滤条件 ----
    f_genes_low  = adata.obs['n_genes_by_counts'] < gf_lo if gf_lo is not None else pd.Series(False, index=adata.obs_names)
    f_genes_high = adata.obs['n_genes_by_counts'] > gf_hi if gf_hi is not None else pd.Series(False, index=adata.obs_names)
    f_mito       = adata.obs['pct_counts_mt'] > mt_hi if mt_hi is not None else pd.Series(False, index=adata.obs_names)
    f_count_hi   = adata.obs['total_counts'] > tc_hi if tc_hi is not None else pd.Series(False, index=adata.obs_names)
    f_cpx        = adata.obs['log_genes_per_umi'] < cpx_lo if cpx_lo is not None else pd.Series(False, index=adata.obs_names)
    f_any = f_genes_low | f_genes_high | f_mito | f_count_hi | f_cpx

    log.info("  Filtering breakdown:")
    if gf_lo is not None:
        log.info("    n_genes < %.0f:       %6d (%.1f%%)", gf_lo,
                 f_genes_low.sum(), 100 * f_genes_low.mean())
    if gf_hi is not None:
        log.info("    n_genes > %.0f:       %6d (%.1f%%)", gf_hi,
                 f_genes_high.sum(), 100 * f_genes_high.mean())
    if mt_hi is not None:
        log.info("    mito > %.2f%%:        %6d (%.1f%%)", mt_hi,
                 f_mito.sum(), 100 * f_mito.mean())
    if tc_hi is not None:
        log.info("    nCount > %.0f:        %6d (%.1f%%)", tc_hi,
                 f_count_hi.sum(), 100 * f_count_hi.mean())
    if cpx_lo is not None:
        log.info("    complexity < %.4f:   %6d (%.1f%%)", cpx_lo,
                 f_cpx.sum(), 100 * f_cpx.mean())
    log.info("    Total (dedup):        %6d (%.1f%%)", f_any.sum(), 100 * f_any.mean())

    mask = ~f_doublet & ~f_any
    adata = adata[mask].copy()
    log.info("  After QC filtering: %d cells", adata.n_obs)
    return adata


# ══════════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════════

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

    # 1. 计算 QC 指标
    compute_qc_metrics(adata, CFG, log)

    # 2. 确定阈值 (MAD 或硬阈值)
    if CFG.use_adaptive_thresholds:
        log.info("Mode: adaptive (MAD × %.1f, nCount MAD × %.1f)",
                 CFG.mad_n_mads, CFG.qc_ncount_max_mad)
        thresholds = _mad_thresholds(adata, CFG, log)
        mode_label = "adaptive (MAD)"
    else:
        log.info("Mode: hard thresholds")
        thresholds = _hard_thresholds(CFG, log)
        mode_label = "hard"

    # 3. 生成诊断图 (在任何过滤之前，展示原始分布 + 阈值线)
    fig_dir = os.path.join(CFG.figure_dir, '02_qc')
    _plot_qc_diagnostics(adata, thresholds, fig_dir, mode_label, CFG, log)

    # 3a. 生成 nFeature KDE 密度峰图
    _fig_dir = os.path.join(CFG.figure_dir, '02_qc')
    _plot_nfeature_kde(adata, _fig_dir, mode_label, CFG, log)

    # 4. 过滤
    adata = filter_cells(adata, thresholds, CFG, log)
    sc.pp.filter_genes(adata, min_cells=CFG.min_cells_per_gene)
    log.info("After gene filtering: %d genes", adata.n_vars)

    safe_write(adata, CFG.qc_h5ad, cfg=CFG)
    log.info("Step 02 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
