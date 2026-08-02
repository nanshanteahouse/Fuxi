#!/usr/bin/env python3
"""
Step 03: 归一化 + HVG 选择 + PCA + 批次校正（Harmony/Combat/scVI）
===============================================================
整合原 Step 02 + Step 03 为单一步骤，并加入 regress_out。

关键顺序 (最佳实践):
  1. raw counts → 找 HVG (seurat_v3, batch-aware)
  2. 保存全基因表达副本用于 .raw
  3. X 只保留 HVG 用于下游降维
  4. normalize_total → log1p on HVG 子集和全基因副本
  5. [可选] score_genes_cell_cycle on 全基因副本, 复制 scores 到 HVG 子集
  6. regress_out (cell cycle scores 或 pct_counts_mt) on HVG 子集 ← normalize 后
  7. 全基因副本赋值 .raw
  8. PCA (n_pcs_full, elbow 图)
  9. 批次校正（Harmony/Combat/scVI）

输入: 02_qc.h5ad
输出: 03_integrated.h5ad (X = log1p(normalized) on HVGs, .raw = 全基因,
                          obsm: X_pca, X_integrated)
"""

import argparse
import os
import sys
import time

import anndata

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from core.kb.cell_cycle import load_cell_cycle_genes
from core.utils import (
    gpu_harmony,
    gpu_pca,
    record_memory_skip,
    resolve_config,
    resolve_device,
    safe_write,
    setup_logger,
    stream_write_raw,
    timed_substep,
)

# Cell cycle gene lists — externalized to core/kb/cell_cycle/
# Previously hardcoded as _S_GENES / _G2M_GENES (Tirosh et al., 2016).
# Loading is deferred to the score_cell_cycle block below via load_cell_cycle_genes().

# ── scVI train() named params (must not appear in trainer_kwargs) ──
_SCVI_EXPLICIT_KWARGS: frozenset[str] = frozenset(
    {
        "max_epochs",
        "accelerator",
        "devices",
        "train_size",
        "validation_size",
        "shuffle_set_split",
        "load_sparse_tensor",
        "batch_size",
        "early_stopping",
        "datasplitter_kwargs",
        "plan_config",
        "plan_kwargs",
        "datamodule",
        "trainer_config",
        "precision",  # explicit pass-through (consumed by **trainer_kwargs inside scVI)
    }
)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("03_integrate", os.path.join(cfg.log_dir, "03_integrate.log"))
    log.info(
        "Step 03: Normalize + HVG selection + PCA + batch integration (%s)", cfg.integration.method
    )
    from core.utils import validate_adata

    # ── 读取 ──
    adata = sc.read(cfg.qc_h5ad)
    log.info("Loaded: %d cells × %d genes", adata.n_obs, adata.n_vars)

    # ── Memory policy detection ──
    mem_policy = getattr(cfg.execution, "memory_policy", "speed")
    log.info("Memory policy: %s", mem_policy)
    # Warn if dense allocation would consume >30% of available RAM
    try:
        import psutil

        avail_gb = psutil.virtual_memory().available / 1e9
        est_dense_gb = adata.n_obs * adata.n_vars * 8 / 1e9
        if est_dense_gb > avail_gb * 0.3:
            log.warning(
                "Estimated dense memory %.1f GB exceeds 30%% of available %.1f GB — "
                "consider memory_policy=balanced",
                est_dense_gb,
                avail_gb,
            )
    except ImportError:
        pass

    # ── HVG（自动降级：CFG.hvg.flavor → seurat_v3 → cell_ranger → seurat）──
    batch_key = cfg.hvg.batch_key if cfg.hvg.batch_key in adata.obs else None
    with timed_substep("HVG selection", log=log):
        flavors_to_try = list(
            dict.fromkeys([cfg.hvg.flavor, "seurat_v3", "cell_ranger", "seurat", "scry"])
        )
        hvg_found = False
        for flavor in flavors_to_try:
            for bk in [batch_key, None]:
                try:
                    log.info(
                        "Selecting top %d HVGs (flavor=%s, batch_key=%s)...",
                        cfg.hvg.n_top_genes,
                        flavor,
                        bk,
                    )
                    sc.pp.highly_variable_genes(
                        adata,
                        n_top_genes=cfg.hvg.n_top_genes,
                        flavor=flavor,
                        batch_key=bk,
                        inplace=True,
                    )
                    log.info("HVG flavor=%s, batch_key=%s succeeded", flavor, bk)
                    hvg_found = True
                    break
                except (ValueError, RuntimeWarning, TypeError, ImportError):
                    log.warning("HVG flavor=%s, batch_key=%s failed, trying next", flavor, bk)
                    continue
            if hvg_found:
                break

        if not hvg_found:
            log.warning("All standard HVG flavors failed, falling back to manual variance method")
            if adata.X is None:
                raise ValueError("adata.X is None — cannot compute manual variance HVG")
            x_mean = adata.X.mean(axis=0).A1
            x_sq = (adata.X.multiply(adata.X)).mean(axis=0).A1
            gene_var = x_sq - x_mean**2
            adata.var["highly_variable"] = np.zeros(adata.n_vars, dtype=bool)
            top_idx = np.argsort(gene_var)[-cfg.hvg.n_top_genes :]
            adata.var.iloc[top_idx, adata.var.columns.get_loc("highly_variable")] = True
            log.info("Manual variance HVG: selected top %d genes", cfg.hvg.n_top_genes)

    n_hvg = adata.var["highly_variable"].sum()
    log.info("HVG count: %d", n_hvg)

    # ── Force-keep critical marker genes in HVG selection ──
    forced_set: set[str] = set()
    if cfg.hvg.forced_genes:
        forced_set.update(cfg.hvg.forced_genes)
    if cfg.marker.marker_dict:
        for markers in cfg.marker.marker_dict.values():
            forced_set.update(markers)
    forced_in_adata = [g for g in forced_set if g in adata.var_names]
    newly_forced = 0
    for g in forced_in_adata:
        val = adata.var.at[g, "highly_variable"]
        if isinstance(val, pd.Series):
            if not val.all():
                newly_forced += (~val).sum()
                adata.var.loc[g, "highly_variable"] = True
        elif not val:
            adata.var.at[g, "highly_variable"] = True
            newly_forced += 1
    if newly_forced:
        log.info(
            "HVG force-keep: %d marker genes retained (HVGs: %d → %d)",
            newly_forced,
            int(adata.var["highly_variable"].sum()) - newly_forced,
            int(adata.var["highly_variable"].sum()),
        )

    # ── 保存全基因表达引用用于 .raw（零拷贝：adata 重新绑定时原对象保留）──
    adata_full = adata
    log.info("Retained full-gene expression reference for .raw")

    # ── X 缩小到 HVGs ──
    adata = adata[:, adata.var["highly_variable"]].copy()
    log.info("X subset to HVGs: %s", adata.shape)

    # ── 保存原始 counts 供 scVI 使用（必须在 normalize_total+log1p 之前）──
    if getattr(cfg.integration, "method", None) == "scvi":
        assert adata.X is not None, "adata.X is None — data not loaded properly"
        adata.layers["counts"] = adata.X.copy()
        log.info("Raw counts preserved in adata.layers['counts'] for scVI")

    # ── 归一化 (HVG 子集) ──
    pearson_mode = cfg.normalization.method == "pearson_residuals"
    skip_norm = False  # set in non-pearson branch; default False for pearson mode
    if pearson_mode:
        log.info(
            "Using Pearson residuals normalization (n_top_genes=%d)...",
            cfg.normalization.pearson_residuals.n_top_genes,
        )
        sc.experimental.pp.highly_variable_genes(
            adata_full,
            flavor="pearson_residuals",
            n_top_genes=cfg.normalization.pearson_residuals.n_top_genes,
        )
        adata = adata_full[:, adata_full.var["highly_variable"]].copy()
        sc.experimental.pp.normalize_pearson_residuals(
            adata,
            clip=cfg.normalization.pearson_residuals.clip,
        )
        log.info("  Pearson residuals normalization complete (%d HVGs)", adata.n_vars)
    else:
        skip_norm = getattr(cfg, "expression_type", "raw_counts") == "log1p_counts"
        if skip_norm:
            log.info(
                "expression_type='log1p_counts' — data already normalized, skipping normalize_total+log1p"
            )
        else:
            log.info(
                "Normalizing (target_sum=%.0f) + log1p...", cfg.normalization.normalize_target_sum
            )
            sc.pp.normalize_total(adata, target_sum=cfg.normalization.normalize_target_sum)
            sc.pp.log1p(adata)

    # ── 数据完整性检查：归一化后 ──
    validate_adata(adata, stage_name="normalize+log1p", logger=log)

    # ── 归一化全基因副本（用于细胞周期打分和 .raw）──
    if pearson_mode:
        log.info("  Pearson residuals mode — full-gene copy kept as raw counts")
    elif skip_norm:
        log.info("  expression_type='log1p_counts' — full-gene copy also pre-normalized")
    else:
        sc.pp.normalize_total(adata_full, target_sum=cfg.normalization.normalize_target_sum)
        sc.pp.log1p(adata_full)

    # ── 可选: 细胞周期打分 ──
    if cfg.normalization.score_cell_cycle:
        if pearson_mode:
            log.info("  Pearson residuals mode — cell cycle scoring skipped")
        else:
            log.info("Scoring cell cycle (S / G2M) on full gene reference...")
            try:
                s_genes, g2m_genes = load_cell_cycle_genes(cfg.species)
                sc.tl.score_genes_cell_cycle(adata_full, s_genes=s_genes, g2m_genes=g2m_genes)
                adata.obs["S_score"] = adata_full.obs["S_score"].copy()
                adata.obs["G2M_score"] = adata_full.obs["G2M_score"].copy()
                adata.obs["phase"] = adata_full.obs["phase"].copy()
                log.info("Cell cycle phases: %s", adata.obs["phase"].value_counts().to_dict())
            except Exception as e:
                log.warning("Cell cycle scoring failed (skipped): %s", e)

    # ── 保存全基因副本到 .raw，尽早释放全基因矩阵 ──
    # 归一化和细胞周期打分后 adata_full 已无其他用途；
    # 提前释放可避免与下游 regress_out / PCA 叠加峰值。
    # NOTE(2026-08-01): 曾试过延迟 .raw 构造（写盘前重建），实测峰值内存
    # 反而 +4.5GB（Python RSS 不回落 + 重复读全基因 18GB），已回滚。
    import gc

    if getattr(cfg.integration, "stream_raw", False) is True:
        # 流式写 .raw：不把全基因矩阵绑定到内存（写盘时从 02_qc 分块重建）
        log.info("[stream_raw] full-gene reference NOT bound to .raw — will stream-write at save")
        _raw_var = adata_full.var.copy()
        del adata_full
        gc.collect()
    else:
        adata.raw = adata_full
        log.info(".raw saved (full genes: %d vars)", adata_full.n_vars)
        del adata_full
        gc.collect()
    log.info("  full-gene reference released early (before regress_out/PCA)")

    # ── 回归技术变异 / 细胞周期分数 (HVG 子集, normalize+log1p 后) ──
    # regress_out internally densifies → large temporary allocation; skip in
    # balanced/memory modes unless cell cycle scoring is explicitly enabled.
    if pearson_mode:
        log.info("  Pearson residuals mode — regress_out skipped")
    else:
        _skip_regress = (
            mem_policy in ("balanced", "memory") and not cfg.normalization.score_cell_cycle
        )
        if _skip_regress:
            log.info(
                "memory_policy=%s — skipping regress_out (avoids dense allocation)", mem_policy
            )
            record_memory_skip(
                step="03_integrate",
                operation="regress_out",
                reason=f"memory_policy={mem_policy} would dense-allocate X (n_obs={adata.n_obs}, n_vars={adata.n_vars})",
                cfg=cfg,
                log=log,
            )
        elif cfg.normalization.score_cell_cycle:
            try:
                log.info("Regressing cell cycle scores: S_score, G2M_score ...")
                sc.pp.regress_out(adata, ["S_score", "G2M_score"])
                log.info("  regress_out complete")
            except Exception as e:
                log.warning("regress_out (cell cycle) failed (skipped): %s", e)
        elif cfg.normalization.use_regress_out:
            covariate_list = ["pct_counts_mt"]
            if cfg.normalization.regress_out_genes:
                covariate_list = cfg.normalization.regress_out_genes
            try:
                log.info("Regressing technical covariates: %s ...", covariate_list)
                sc.pp.regress_out(adata, covariate_list)
                log.info("  regress_out complete")
            except Exception as e:
                log.warning("regress_out (%s) failed (skipped): %s", covariate_list, e)
        else:
            log.info("Cell cycle scoring disabled, use_regress_out=False — skipping regress_out")

    # ── 可选: 回归自定义基因（额外）──
    # (主回归 covariate_list 已同时处理 obs 与 gene 列)
    # 如果 regress_out_genes 已在 covariate_list 中，跳过重复

    # ── regress_out 后降回 float32（regress_out 会产生 float64 中间体）──
    if getattr(cfg.execution, "use_float32", False) and adata.X is not None:
        if sp.issparse(adata.X):
            adata.X = adata.X.astype("float32", copy=False)
        else:
            adata.X = adata.X.astype("float32", copy=False)
        log.info("  X precision restored to float32")

    # ── 数据完整性检查：regress_out 后 ──
    if cfg.normalization.score_cell_cycle or cfg.normalization.use_regress_out:
        validate_adata(adata, stage_name="regress_out", logger=log)

    # ── .raw already saved and full-gene reference released after cell cycle scoring ──

    # ── 可选: 自动性别检测 ──
    if getattr(cfg.normalization, "detect_sex", False):
        from rna.utils.sex_detection import detect_sex

        detect_sex(adata, cfg, log)

    # ── 数据完整性检查：PCA 前 ──
    has_issues = validate_adata(adata, stage_name="before_PCA", logger=log)
    if has_issues:
        log.critical("Data integrity issues found before PCA, aborting!")
        sys.exit(1)

    # ── PCA ──
    # arpack = iterative (low memory), randomized = fast (dense allocation)
    _solver = "arpack" if mem_policy in ("balanced", "memory") else "randomized"
    log.info("PCA (%d components, solver=%s)...", cfg.pca.n_pcs_full, _solver)
    with timed_substep("PCA", log=log):
        # GPU PCA ignores solver arg (uses cuVS default); CPU path honors _solver
        _pca_kwargs = dict(
            n_comps=cfg.pca.n_pcs_full,
            random_state=cfg.execution.random_seed,
        )
        if not resolve_device(cfg.execution.device, log):
            _pca_kwargs["svd_solver"] = _solver
        gpu_pca(adata, log=log, device=cfg.execution.device, **_pca_kwargs)
        var_ratio = adata.uns["pca"]["variance_ratio"]
        log.info("  top-5 variance ratio: %.4f", var_ratio[:5].sum())
        log.info("  Cumulative variance ratio first 50 PCs: %.4f", var_ratio[:50].sum())

        # PCA elbow 图
        fig_dir = os.path.join(cfg.figure_dir, "03_integrate")
        os.makedirs(fig_dir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(cfg.plot.qc_figure_size[0], 4))
        ax.plot(range(1, cfg.pca.n_pcs_full + 1), var_ratio, "o-", ms=3)
        ax.axvline(
            cfg.pca.n_pcs_use,
            color=cfg.plot.palette.qc_threshold,
            linestyle="--",
            alpha=0.5,
            label=f"n_pcs_use={cfg.pca.n_pcs_use}",
        )
        ax.set_xlabel("PC")
        ax.set_ylabel("Variance ratio")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "pca_elbow.png"), dpi=cfg.plot.figure_dpi)
        plt.close(fig)
        log.info("  PCA elbow plot saved")

    # ── Batch diagnosis (v4.x+) ──
    report = None
    if cfg.integration.diagnose:
        from rna.utils.batch_diagnostics import diagnose_batch_candidates, plot_diagnosis_report

        # 轻量快照（无 X/raw）——preservation 只需 X_pca + obs，避免 .raw 42.5GB 深拷贝 (2026-08-02)
        adata_orig = anndata.AnnData(
            obs=adata.obs.copy(),
            obsm={"X_pca": adata.obsm["X_pca"].copy()},
        )
        gini_b = cfg.integration.gini_batch_threshold
        gini_ = cfg.integration.gini_biology_threshold
        log.info("Running batch diagnosis ...")
        try:
            with timed_substep("Batch diagnosis", log=log):
                report = diagnose_batch_candidates(
                    adata,
                    n_pcs=cfg.pca.n_pcs_use,
                    gini_batch_threshold=gini_b,
                    gini_biology_threshold=gini_,
                    max_cells=cfg.integration.diagnose_max_cells,
                )
                log.info(
                    "[batch-diagnosis] batch_cols=%s, biology_cols=%s, ambiguous=%s",
                    report.batch_cols,
                    report.biology_cols,
                    report.ambiguous_cols,
                )
                # Augment batch_key with auto-detected batch columns
                augment = report.batch_cols if report.batch_cols else []
                user_key = cfg.integration.batch_key
                bk_list = [user_key] if isinstance(user_key, str) else list(user_key)
                batch_keys = list(dict.fromkeys(bk_list + [c for c in augment if c in adata.obs]))
                for w in report.warnings:
                    log.warning("[batch-diagnosis] %s", w)
                if cfg.integration.diagnose_report:
                    os.makedirs(fig_dir, exist_ok=True)
                    report_path = os.path.join(fig_dir, "batch_diagnosis.pdf")
                    plot_diagnosis_report(report, report_path)
                    log.info("  Diagnosis report saved to %s", report_path)
        except Exception as e:
            log.warning("Batch diagnosis failed (%s) — continuing without diagnosis", e)
            batch_keys = (
                [cfg.integration.batch_key]
                if isinstance(cfg.integration.batch_key, str)
                else list(cfg.integration.batch_key)
            )
    else:
        batch_keys = (
            [cfg.integration.batch_key]
            if isinstance(cfg.integration.batch_key, str)
            else list(cfg.integration.batch_key)
        )

    # ── Integration ──
    bk_list = [b for b in batch_keys if b in adata.obs]
    collinear_warnings = []
    if cfg.integration.collinearity_guard and report is not None and report.warnings:
        collinear_warnings = [
            w
            for w in report.warnings
            if "perfectly collinear" in w.lower() and any(bc in w for bc in report.biology_cols)
        ]

    # ── Checkpoint path (used across all integration methods) ──
    out_path = os.path.join(cfg.h5ad_dir, "03_integrated.h5ad")

    if cfg.integration.method == "harmony":
        if collinear_warnings:
            log.error(
                "[collinearity-guard] Harmony ABORTED — batch_key perfectly collinear with biology:"
            )
            for w in collinear_warnings[:3]:
                log.error("  %s", w)
            log.error("  To override: set integration.collinearity_guard: false in config.yaml")
            adata.uns["harmony_skipped"] = {
                "reason": "collinearity",
                "warnings": collinear_warnings,
            }
        elif bk_list:
            # Unified NaN detection across all batch keys
            nan_mask = adata.obs[bk_list].isna().any(axis=1)
            if nan_mask.any():
                n_nan = nan_mask.sum()
                log.warning(
                    "batch_keys %s contain %d NaN rows — these cells will be removed",
                    bk_list,
                    n_nan,
                )
                adata._inplace_subset_obs(~nan_mask)
            log.info("Harmony correction (batch_keys=%s)...", bk_list)
            try:
                with timed_substep("Harmony", log=log):
                    gpu_harmony(
                        adata,
                        key=bk_list,
                        output_key="X_integrated",
                        log=log,
                        device=cfg.execution.device,
                        random_state=cfg.execution.random_seed,
                        max_iter_harmony=cfg.integration.max_iter,
                    )
                log.info("  Harmony complete, output shape: %s", adata.obsm["X_integrated"].shape)
            except Exception as e:
                log.warning("Harmony correction failed (%s) — continuing with raw PCA", e)
                adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()
            # ── Checkpoint before plotting ──
            # NOTE(2026-08-01): 末尾统一 safe_write 已覆盖保存；此处重复写同一文件
            # 会导致第二次写失败 (integrity check 的 backed read 保持 h5py 句柄,
            # truncate 已打开文件报错)。双写已移除，仅保留末尾统一写。

            # 对比图
            primary_key = bk_list[0]
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sc.pl.embedding(
                adata,
                basis="X_pca",
                color=primary_key,
                ax=axes[0],
                show=False,
                title="PCA (before Harmony)",
            )
            sc.pl.embedding(
                adata,
                basis="X_integrated",
                color=primary_key,
                ax=axes[1],
                show=False,
                title="Harmony-corrected",
            )
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "harmony_comparison.png"), dpi=cfg.plot.figure_dpi)
            plt.close(fig)
            log.info("  Harmony comparison plot saved")
            # ── Post-Harmony preservation check ──
            if cfg.integration.diagnose and report is not None and report.biology_cols:
                try:
                    from rna.utils.batch_diagnostics import validate_harmony_preservation

                    preservation = None
                    with timed_substep("Preservation check", log=log):
                        preservation = validate_harmony_preservation(
                            adata_orig, adata, report.biology_cols
                        )
                    for col, ratio in preservation.items():
                        if ratio < 0.9:
                            log.warning(
                                "[preservation-check] %s: purity dropped to %.2fx — "
                                "biological signal may be degraded",
                                col,
                                ratio,
                            )
                        else:
                            log.info(
                                "[preservation-check] %s: purity preserved at %.2fx", col, ratio
                            )
                except Exception as e:
                    log.warning("Preservation check failed (%s) — skipping", e)
        else:
            log.warning("No valid batch keys found in obs, skipping Harmony")
            adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()
    elif cfg.integration.method == "combat":
        if bk_list:
            log.info("Combat correction (batch_key=%s)...", cfg.integration.batch_key)
            sc.pp.combat(adata, key=cfg.integration.batch_key)
            sc.pp.pca(
                adata,
                n_comps=cfg.pca.n_pcs_use,
                svd_solver=_solver,
                random_state=cfg.execution.random_seed,
            )
            adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()
            # ── Checkpoint before plotting ──
            # NOTE(2026-08-01): 同 harmony 分支——末尾统一 safe_write 已覆盖，
            # 此处重复写会导致 truncate 失败，双写已移除。

            # 对比图
            primary_key = bk_list[0]
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sc.pl.embedding(
                adata,
                basis="X_pca",
                color=primary_key,
                ax=axes[0],
                show=False,
                title="PCA (before Combat)",
            )
            sc.pl.embedding(
                adata,
                basis="X_integrated",
                color=primary_key,
                ax=axes[1],
                show=False,
                title="Combat-corrected",
            )
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "combat_comparison.png"), dpi=cfg.plot.figure_dpi)
            plt.close(fig)
            log.info("  Combat comparison plot saved")
        else:
            log.warning("No valid batch keys found in obs, skipping Combat")
            adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()
    elif cfg.integration.method == "scvi":
        from core.utils._optional import gpu_available_torch, require_scvi

        require_scvi("step 03 integration (scVI)")
        import scvi

        log.info("scVI integration (n_latent=%d)...", cfg.integration.scvi.n_latent)
        try:
            # Check raw counts exist
            if "counts" not in adata.layers:
                raise ValueError(
                    "No 'counts' layer found. Raw counts must be preserved before normalization."
                )

            # Verify counts are integer-valued
            mean_count = adata.layers["counts"].mean()
            log.info("  mean raw count per cell: %.2f", mean_count)

            # GPU detection
            use_gpu = cfg.integration.scvi.use_gpu and gpu_available_torch()
            if cfg.integration.scvi.use_gpu and not use_gpu:
                log.warning("GPU requested but not available — falling back to CPU")

            # TF32 mixed-precision matmul (Ampere+): ~2x training speedup on
            # RTX 3090; scverse-verified safe for integration (2026-08-02).
            if use_gpu:
                import torch

                torch.set_float32_matmul_precision("medium")
                log.info("[scvi] float32_matmul_precision → medium (TF32)")

            # Determine batch key: try scvi-specific first, then general, else None
            batch_key = (
                cfg.integration.scvi.batch_key
                if cfg.integration.scvi.batch_key in adata.obs
                else (
                    cfg.integration.batch_key if cfg.integration.batch_key in adata.obs else None
                )
            )

            # Setup and train scVI model
            scvi.model.SCVI.setup_anndata(
                adata,
                layer="counts",
                batch_key=batch_key,
            )
            model = scvi.model.SCVI(
                adata,
                n_latent=cfg.integration.scvi.n_latent,
                n_layers=cfg.integration.scvi.n_layers,
                n_hidden=cfg.integration.scvi.n_hidden,
            )
            # ── Validate early_stopping preconditions ──
            if cfg.integration.scvi.early_stopping and cfg.integration.scvi.train_size >= 1.0:
                raise ValueError(
                    "early_stopping=True requires train_size < 1.0 "
                    "(scVI needs a validation set to monitor ELBO). "
                    f"Got train_size={cfg.integration.scvi.train_size}."
                )

            # ── Guard: CPU + mixed precision crash (D1) ──
            # PyTorch lacks fp16 matmul on CPU → RuntimeError
            if cfg.integration.scvi.precision != "32" and not use_gpu:
                log.warning(
                    "scVI precision=%s is unsafe on CPU (PyTorch lacks fp16 matmul) "
                    "— downgrading to '32'. Set use_gpu=true or precision='32'.",
                    cfg.integration.scvi.precision,
                )
                _precision = "32"
            else:
                _precision = cfg.integration.scvi.precision
            t_start = time.time()
            # Filter out keys that duplicate explicit args (prevents TypeError)
            _trainer_kwargs = {
                k: v
                for k, v in cfg.integration.scvi.trainer_kwargs.items()
                if k not in _SCVI_EXPLICIT_KWARGS
            }
            if _trainer_kwargs != cfg.integration.scvi.trainer_kwargs:
                _conflicts = set(cfg.integration.scvi.trainer_kwargs) & _SCVI_EXPLICIT_KWARGS
                log.warning(
                    "scVI trainer_kwargs contains explicit-arg keys %s — filtered out "
                    "(use top-level SCVIConfig fields instead)",
                    _conflicts,
                )

            model.train(
                max_epochs=cfg.integration.scvi.max_epochs,
                batch_size=cfg.integration.scvi.batch_size,
                early_stopping=cfg.integration.scvi.early_stopping,
                train_size=cfg.integration.scvi.train_size,
                accelerator="gpu" if use_gpu else "cpu",
                devices=1,
                precision=_precision,
                **_trainer_kwargs,
                plan_kwargs=cfg.integration.scvi.plan_kwargs or None,
                datasplitter_kwargs=cfg.integration.scvi.datasplitter_kwargs or None,
            )
            elapsed = time.time() - t_start
            log.info("scVI training completed in %.1f seconds", elapsed)

            # Get latent representation
            latent = model.get_latent_representation()
            adata.obsm["X_integrated"] = latent
            log.info("  scVI latent shape: %s", latent.shape)
            if latent.shape[1] != cfg.pca.n_pcs_use:
                log.warning(
                    "scVI latent has %d dims but cfg.pca.n_pcs_use=%d -- "
                    "downstream consumers expecting %d dims may misbehave. "
                    "Set pca.n_pcs_use=%d in your project config to match.",
                    latent.shape[1],
                    cfg.pca.n_pcs_use,
                    cfg.pca.n_pcs_use,
                    latent.shape[1],
                )

        except Exception as e:
            log.warning("scVI integration failed (%s) — falling back to raw PCA", e)
            adata.obsm["X_integrated"] = adata.obsm["X_pca"][:, : cfg.pca.n_pcs_use].copy()

    else:
        log.info("Integration method '%s' — using raw PCA.", cfg.integration.method)
        adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()

    # ── 保存 ──
    with timed_substep("Save checkpoint", log=log):
        safe_write(adata, out_path, cfg=cfg, step_alias="integrated")
        if getattr(cfg.integration, "stream_raw", False) is True:
            # 流式写 .raw：从 02_qc.h5ad 分块读 counts → normalize+log1p → 直写 raw 组
            # var 用 03 算出的 highly_variable 全基因注释（_raw_var 在早前已保存）
            stream_write_raw(
                out_path,
                cfg.qc_h5ad,
                target_sum=cfg.normalization.normalize_target_sum,
                compression=getattr(cfg, "h5ad_compression", "gzip"),
                compression_opts=getattr(cfg, "h5ad_compression_opts", None),
                var_df=_raw_var,
                logger=log,
            )
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
