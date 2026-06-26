#!/usr/bin/env python3
"""
Step 01a: Scrublet 双细胞检测 (per sample, joblib 并行)
=========================================================
从 02_qc.py 中独立出的 Scrublet 步骤，不含 QC 指标计算或过滤。

输入: 00_raw.h5ad
输出: 01_doublet.h5ad (含 doublet_scores / predicted_doublet 列)
"""
import sys, os, time, argparse, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import numpy as np
import scanpy as sc
from joblib import Parallel, delayed


def _resolve_doublet_rate(cfg, n_cells: int) -> float:
    """返回用于 Scrublet 的 expected_doublet_rate。

    如果 config 显式设置了值，直接使用；否则按 10X 官方拟合：
      y (%) = 0.000759 * x + 0.053    (x = recovered cell 数)
    钳位在 [0.004, 0.15]。
    """
    if cfg.scrublet_expected_doublet_rate is not None:
        return cfg.scrublet_expected_doublet_rate
    rate = 0.00000759 * n_cells + 0.00000053
    return min(max(rate, 0.004), 0.15)


def run_scrublet_sample(adata_sub, sample_name, cfg):
    try:
        import scrublet as scr
        import scipy.sparse as sp
        expected_rate = _resolve_doublet_rate(cfg, adata_sub.n_obs)
        scrub = scr.Scrublet(
            adata_sub.X if isinstance(adata_sub.X, sp.spmatrix) else sp.csr_matrix(adata_sub.X),
            expected_doublet_rate=expected_rate,
            random_state=cfg.random_seed,
        )
        scores, predicted = scrub.scrub_doublets(
            min_counts=cfg.scrublet_min_counts,
            min_cells=cfg.scrublet_min_cells,
            min_gene_variability_pctl=cfg.scrublet_min_gene_var_pctl,
            n_prin_comps=cfg.scrublet_n_prin_comps,
        )
        if predicted is None:
            fallback = expected_rate
            warnings.warn(f"Scrublet auto-threshold failed for {sample_name}, "
                          f"falling back to manual threshold={fallback}")
            predicted = scrub.call_doublets(threshold=fallback)
        return scores, predicted
    except Exception as e:
        warnings.warn(f"Scrublet failed for {sample_name}: {e}")
        return np.zeros(adata_sub.n_obs), np.zeros(adata_sub.n_obs, dtype=bool)


def detect_doublets_parallel(adata, cfg, log):
    if not cfg.run_scrublet:
        log.info("Scrublet disabled, skipping doublet detection.")
        adata.obs['doublet_scores'] = 0.0
        adata.obs['predicted_doublet'] = False
        return

    # TPM/CPM/FPKM/log1p 数据 → Scrublet 的负二项分布假设不成立
    if cfg.expression_type != "raw_counts":
        log.warning(
            "Scrublet is designed for raw UMI counts. "
            "expression_type='%s' violates the negative-binomial assumption. "
            "Disabling Scrublet. Set run_scrublet=False to suppress this warning.",
            cfg.expression_type
        )
        adata.obs['doublet_scores'] = 0.0
        adata.obs['predicted_doublet'] = False
        return

    log.info("Running Scrublet (per sample, parallel)...")
    configured_key = cfg.scrublet_batch_key
    if configured_key in adata.obs:
        groupby_col = configured_key
        log.info("  Using configured batch column: %s", groupby_col)
    else:
        log.warning("Configured batch column '%s' not in adata.obs, falling back to 'sample'/'stage'",
                    configured_key)
        groupby_col = 'sample' if 'sample' in adata.obs else 'stage'
    if groupby_col not in adata.obs:
        log.warning("Group column (%s) not found, running Scrublet on all data.", groupby_col)
        scores, pred = run_scrublet_sample(adata, "all", cfg)
        adata.obs['doublet_scores'] = scores
        adata.obs['predicted_doublet'] = pred
        log.info("  Predicted doublets: %d / %d (%.1f%%)",
                 pred.sum(), adata.n_obs, 100 * pred.mean())
        return

    sample_groups = adata.obs.groupby(groupby_col, observed=True)

    # Memory-aware scheduling: large samples (>15000 cells) serially,
    # small samples (<=15000) in parallel to avoid OOM on big groups.
    MEMORY_THRESHOLD = 15000
    large_names, large_subsets, large_idxs = [], [], []
    small_names, small_subsets, small_idxs = [], [], []
    for name, idx in sample_groups.indices.items():
        if len(idx) > MEMORY_THRESHOLD:
            large_names.append(name)
            large_subsets.append(adata[idx])
            large_idxs.append(idx)
        else:
            small_names.append(name)
            small_subsets.append(adata[idx])
            small_idxs.append(idx)

    results = []

    if large_names:
        log.info("  Large groups (%s) — processing serially",
                 ", ".join(f"{n}({len(i)} cells)" for n, i in zip(large_names, large_idxs)))
    for sub, name in zip(large_subsets, large_names):
        results.append(run_scrublet_sample(sub, name, cfg))

    if small_subsets:
        n_jobs = min(cfg.n_jobs or os.cpu_count() or 1, len(small_names))
        log.info("  Small samples — processing %d groups in parallel (n_jobs=%d)",
                 len(small_names), n_jobs)
        small_results = Parallel(n_jobs=n_jobs)(
            delayed(run_scrublet_sample)(sub, name, cfg)
            for sub, name in zip(small_subsets, small_names)
        )
        results.extend(small_results)

    all_scores = np.zeros(adata.n_obs)
    all_pred = np.zeros(adata.n_obs, dtype=bool)
    all_names = large_names + small_names
    all_idxs = list(large_idxs) + list(small_idxs)
    for (scores, pred), name, idx in zip(results, all_names, all_idxs):
        all_scores[idx] = scores
        all_pred[idx] = pred
        used_rate = _resolve_doublet_rate(cfg, len(idx))
        log.info("  Sample %s: %d / %d doublets (%.1f%%) [expected_rate=%.4f]",
                 name, pred.sum(), len(idx), 100 * pred.mean(), used_rate)

    adata.obs['doublet_scores'] = all_scores
    adata.obs['predicted_doublet'] = all_pred
    log.info("  Total predicted doublets: %d / %d (%.1f%%)",
             all_pred.sum(), adata.n_obs, 100 * all_pred.mean())


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("01_doublet", os.path.join(CFG.log_dir, "01_doublet.log"))
    log.info("Step 01a: Scrublet doublet detection")

    adata = sc.read(CFG.raw_h5ad)
    log.info("Loaded: %s — %d cells × %d genes",
             CFG.raw_h5ad, adata.n_obs, adata.n_vars)

    detect_doublets_parallel(adata, CFG, log)

    out_path = os.path.join(CFG.h5ad_dir, "01_doublet.h5ad")
    safe_write(adata, out_path, cfg=CFG)
    log.info("Step 01a complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
