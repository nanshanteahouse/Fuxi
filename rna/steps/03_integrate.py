#!/usr/bin/env python3
"""
Step 03: 归一化 + HVG 选择 + PCA + Harmony 批次校正（整合版）
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
  9. Harmony 批次校正

输入: 02_qc.h5ad
输出: 03_integrated.h5ad (X = log1p(normalized) on HVGs, .raw = 全基因,
                          obsm: X_pca, X_pca_harmony)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np

# Cell cycle gene lists (Tirosh et al., 2016) for sc.tl.score_genes_cell_cycle
_S_GENES = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP',
    'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','PIR51','MCM10',
    'RFWD3','FANCI','TK1','CDC45','CDC6','DSCC1','EXO1','TIPIN','E2F8',
    'GINS4','CASP8AP2','GMPS','BRIP1','CLSPN','HAT1','RRM2','RAD51',
    'RPA3','BRCA1',
]
_G2M_GENES = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
    'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2',
    'CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1',
    'KIF20B','HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1',
    'NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1',
    'ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA',
]


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("03_integrate", os.path.join(CFG.log_dir, "03_integrate.log"))
    log.info("Step 03: Normalize + HVG selection + PCA + Harmony (integrated)")
    from core.utils import validate_adata

    # ── 读取 ──
    adata = sc.read(CFG.qc_h5ad)
    log.info("Loaded: %d cells × %d genes", adata.n_obs, adata.n_vars)

    # ── HVG（自动降级：CFG.hvg_flavor → seurat_v3 → cell_ranger → seurat）──
    batch_key = CFG.hvg_batch_key if CFG.hvg_batch_key in adata.obs else None
    flavors_to_try = list(dict.fromkeys([CFG.hvg_flavor, 'seurat_v3', 'cell_ranger', 'seurat', 'scry']))
    hvg_found = False
    for flavor in flavors_to_try:
        for bk in [batch_key, None]:
            try:
                log.info("Selecting top %d HVGs (flavor=%s, batch_key=%s)...",
                         CFG.n_top_genes, flavor, bk)
                sc.pp.highly_variable_genes(
                    adata, n_top_genes=CFG.n_top_genes, flavor=flavor, batch_key=bk, inplace=True,
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
        X_mean = adata.X.mean(axis=0).A1
        X_sq = (adata.X.multiply(adata.X)).mean(axis=0).A1
        gene_var = X_sq - X_mean**2
        adata.var['highly_variable'] = np.zeros(adata.n_vars, dtype=bool)
        top_idx = np.argsort(gene_var)[-CFG.n_top_genes:]
        adata.var.iloc[top_idx, adata.var.columns.get_loc('highly_variable')] = True
        log.info("Manual variance HVG: selected top %d genes", CFG.n_top_genes)

    n_hvg = adata.var['highly_variable'].sum()
    log.info("HVG count: %d", n_hvg)

    # ── 保存全基因表达引用用于 .raw（零拷贝：adata 重新绑定时原对象保留）──
    adata_full = adata
    log.info("Retained full-gene expression reference for .raw")

    # ── X 缩小到 HVGs ──
    adata = adata[:, adata.var['highly_variable']].copy()
    log.info("X subset to HVGs: %s", adata.shape)

    # ── 归一化 (HVG 子集) ──
    skip_norm = getattr(CFG, 'expression_type', 'raw_counts') == 'log1p_counts'
    if skip_norm:
        log.info("expression_type='log1p_counts' — data already normalized, skipping normalize_total+log1p")
    else:
        log.info("Normalizing (target_sum=%.0f) + log1p...", CFG.normalize_target_sum)
        sc.pp.normalize_total(adata, target_sum=CFG.normalize_target_sum)
        sc.pp.log1p(adata)

    # ── 数据完整性检查：归一化后 ──
    validate_adata(adata, stage_name="normalize+log1p", logger=log)

    # ── 归一化全基因副本（用于细胞周期打分和 .raw）──
    if skip_norm:
        log.info("  expression_type='log1p_counts' — full-gene copy also pre-normalized")
    else:
        sc.pp.normalize_total(adata_full, target_sum=CFG.normalize_target_sum)
        sc.pp.log1p(adata_full)

    # ── 可选: 细胞周期打分 ──
    if CFG.score_cell_cycle:
        log.info("Scoring cell cycle (S / G2M) on full gene reference...")
        try:
            sc.tl.score_genes_cell_cycle(adata_full, s_genes=_S_GENES, g2m_genes=_G2M_GENES)
            adata.obs['S_score'] = adata_full.obs['S_score'].copy()
            adata.obs['G2M_score'] = adata_full.obs['G2M_score'].copy()
            adata.obs['phase'] = adata_full.obs['phase'].copy()
            log.info("Cell cycle phases: %s", adata.obs['phase'].value_counts().to_dict())
        except Exception as e:
            log.warning("Cell cycle scoring failed (skipped): %s", e)

    # ── 回归技术变异 / 细胞周期分数 (HVG 子集, normalize+log1p 后) ──
    if CFG.score_cell_cycle:
        try:
            log.info("Regressing cell cycle scores: S_score, G2M_score ...")
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
            log.info("  regress_out complete")
        except Exception as e:
            log.warning("regress_out (cell cycle) failed (skipped): %s", e)
    elif CFG.use_regress_out:
        try:
            log.info("Regressing technical covariates: pct_counts_mt ...")
            sc.pp.regress_out(adata, ['pct_counts_mt'])
            log.info("  regress_out complete")
        except Exception as e:
            log.warning("regress_out (pct_counts_mt) failed (skipped): %s", e)
    else:
        log.info("Cell cycle scoring disabled, use_regress_out=False — skipping regress_out")

    # ── 可选: 回归自定义基因 ──
    if CFG.regress_out_genes:
        valid_genes = [g for g in CFG.regress_out_genes if g in adata.var_names]
        if valid_genes:
            log.info("Regressing custom genes: %s ...", valid_genes)
            sc.pp.regress_out(adata, valid_genes)
            log.info("  regress_out (custom genes) complete")
        else:
            log.warning("regress_out_genes specified but none found in data: %s", CFG.regress_out_genes)

    # ── regress_out 后降回 float32（regress_out 会产生 float64 中间体）──
    if getattr(CFG, 'use_float32', False):
        if sp.issparse(adata.X):
            adata.X = adata.X.astype('float32', copy=False)
        else:
            adata.X = adata.X.astype('float32', copy=False)
        log.info("  X precision restored to float32")

    # ── 数据完整性检查：regress_out 后 ──
    if CFG.score_cell_cycle or CFG.use_regress_out:
        validate_adata(adata, stage_name="regress_out", logger=log)

    # ── 保存全基因副本到 .raw ──
    adata.raw = adata_full
    log.info(".raw saved (full genes: %d vars)", adata_full.n_vars)

    # ── 数据完整性检查：PCA 前 ──
    has_issues = validate_adata(adata, stage_name="before_PCA", logger=log)
    if has_issues:
        log.critical("Data integrity issues found before PCA, aborting!")
        sys.exit(1)

    # ── PCA ──
    log.info("PCA (%d components)...", CFG.n_pcs_full)
    sc.pp.pca(adata, n_comps=CFG.n_pcs_full,
              svd_solver='randomized', random_state=CFG.random_seed)
    var_ratio = adata.uns['pca']['variance_ratio']
    log.info("  top-5 variance ratio: %.4f", var_ratio[:5].sum())
    log.info("  Cumulative variance ratio first 50 PCs: %.4f", var_ratio[:50].sum())

    # PCA elbow 图
    fig_dir = os.path.join(CFG.figure_dir, '03_integrate')
    os.makedirs(fig_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(1, CFG.n_pcs_full + 1), var_ratio, 'o-', ms=3)
    ax.axvline(CFG.n_pcs_use, color='red', linestyle='--', alpha=0.5,
               label=f'n_pcs_use={CFG.n_pcs_use}')
    ax.set_xlabel('PC'); ax.set_ylabel('Variance ratio')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, 'pca_elbow.png'), dpi=150)
    plt.close(fig)
    log.info("  PCA elbow plot saved")

    # ── Harmony ──
    if CFG.use_harmony:
        from harmony import harmonize
        batch_key = CFG.harmony_batch_key
        if batch_key not in adata.obs:
            log.warning("Harmony batch_key '%s' not in obs, skipping correction", batch_key)
        else:
            # Check for NaN in batch_key column
            if bool(adata.obs[batch_key].isna().any()):
                n_nan = adata.obs[batch_key].isna().sum()
                log.warning("batch_key '%s' contains %d NaN — these cells will be removed", batch_key, n_nan)
                adata._inplace_subset_obs(~adata.obs[batch_key].isna())
            log.info("Harmony correction (batch_key=%s)...", batch_key)
            try:
                Z = harmonize(
                    adata.obsm['X_pca'][:, :CFG.n_pcs_use],
                    adata.obs,
                    batch_key=batch_key,
                    random_state=CFG.random_seed,
                    max_iter_harmony=CFG.harmony_max_iter,
                )
                adata.obsm['X_pca_harmony'] = Z
                log.info("  Harmony complete, output shape: %s", Z.shape)
            except Exception as e:
                log.warning("Harmony correction failed (%s) — continuing with raw PCA", e)
                adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()
            # 对比图
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sc.pl.embedding(adata, basis='X_pca', color=batch_key,
                            ax=axes[0], show=False, title='PCA (before Harmony)')
            sc.pl.embedding(adata, basis='X_pca_harmony', color=batch_key,
                            ax=axes[1], show=False, title='Harmony-corrected')
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, 'harmony_comparison.png'), dpi=150)
            plt.close(fig)
            log.info("  Harmony comparison plot saved")
    else:
        log.info("Harmony disabled, using raw PCA.")
        adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()

    # ── 保存 ──
    out_path = os.path.join(CFG.h5ad_dir, "03_integrated.h5ad")
    safe_write(adata, out_path, cfg=CFG)
    log.info("Step 03 complete, took %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
