#!/usr/bin/env python3
"""
Step 05: 细胞类型自动注释 (Major Lineage — AI + Score_genes 双模式)
=====================================================================
  双模式注释策略:
    1. AI 模式 (首选): 基于 marker 基因 + LLM 智能判断细胞类型
    2. Score_genes 模式 (回退): 基于已知 marker 基因打分
  输出主要细胞类型 (major lineage) 及亚型/状态/置信度信息。

输入: 04_clustered.h5ad
输出: 05_annotated.h5ad (新增 cell_type, cell_subtype, cell_state, annot_confidence, ... 列)
"""
import sys, os, time, argparse, json
# Add repo root so `from core.*` and `from rna.*` resolve correctly
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np
import logging

log: logging.Logger


def _warn_if_low_coverage(adata, log):
    """Emit WARNING if >50% of cells are annotated as Unknown."""
    if 'cell_type' not in adata.obs:
        return
    ct_counts = adata.obs['cell_type'].value_counts()
    n_unknown = ct_counts.get('Unknown', 0)
    if n_unknown > 0 and n_unknown / adata.n_obs > 0.5:
        log.warning("⚠️  %.0f%% of cells (%d/%d) annotated as 'Unknown' — KB may not cover this dataset",
                    n_unknown / adata.n_obs * 100, n_unknown, adata.n_obs)



# ═══════════════════════════════════════════════════════════════════════
#  旧有注释函数 (Score_genes 模式)
# ═══════════════════════════════════════════════════════════════════════

def run_annotation(adata, marker_dict, logger):
    """基于 marker 基因打分的细胞类型注释 (来自 05_annotate.py 原有逻辑)。"""
    if not marker_dict:
        logger.warning("marker_dict not configured, skipping annotation.")
        adata.obs['cell_type'] = adata.obs['leiden'].astype(str)
        return

    cell_types = list(marker_dict.keys())
    logger.info("Score-based annotation: %d candidate types", len(cell_types))

    for ct in cell_types:
        genes = marker_dict[ct]
        genes_present = [g for g in genes if g in adata.raw.var_names]
        if not genes_present:
            logger.warning("  %s: no marker genes found in data", ct)
            adata.obs[f'score_{ct}'] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=genes_present,
                          score_name=f'score_{ct}', random_state=42)

    # 每个聚类取最高分的类型
    score_cols = [f'score_{ct}' for ct in cell_types]
    cluster_scores = adata.obs.groupby('leiden', observed=True)[score_cols].mean()
    best_match = cluster_scores.idxmax(axis=1)
    best_ct = best_match.str.replace('score_', '')

    cluster_to_ct = dict(zip(best_ct.index, best_ct.values))
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_ct).astype('category')

    logger.info("Cluster → cell type mapping:")
    for label in sorted(adata.obs['leiden'].unique()):
        ct = cluster_to_ct[label]
        max_score = cluster_scores.loc[label, f'score_{ct}']
        logger.info("  Cluster %s → %s (score=%.3f)", label, ct, max_score)

    # 置信度: 最高分与次高分之差
    if len(cell_types) >= 2:
        sorted_scores = cluster_scores.apply(
            lambda row: row.sort_values(ascending=False).values, axis=1, result_type='expand'
        )
        confidence = sorted_scores.iloc[:, 0] - sorted_scores.iloc[:, 1]
        adata.obs['annotation_confidence'] = adata.obs['leiden'].map(confidence).astype(float).values
        low_conf = (adata.obs['annotation_confidence'] < 0.02).sum()
        if low_conf > 0:
            logger.info("  Low-confidence cells (<0.02): %d (%.1f%%)",
                        low_conf, 100 * low_conf / adata.n_obs)

    logger.info("Annotation complete: %d cell types", adata.obs['cell_type'].nunique())


def run_subclustering(adata, CFG, subcluster_types, resolution, min_cells, logger):
    """基于 parent cell_type 的子聚类 (来自 05_annotate.py 原有逻辑)。"""
    if not subcluster_types:
        logger.info("Subcluster types not configured, skipping.")
        adata.obs['cell_type_sub'] = adata.obs['cell_type'].astype(str)
        return

    logger.info("Subclustering: %s (resolution=%.1f)...", subcluster_types, resolution)
    adata.obs['cell_type_sub'] = adata.obs['cell_type'].astype(str)

    for parent_type in subcluster_types:
        mask = adata.obs['cell_type'] == parent_type
        n_cells = mask.sum()
        if n_cells < min_cells:
            logger.info("  %s: too few cells (%d < %d), skipping", parent_type, n_cells, min_cells)
            continue

        logger.info("  Subclustering %s (%d cells)...", parent_type, n_cells)
        sub = adata[mask].copy()
        sc.pp.neighbors(sub, n_pcs=50, use_rep='X_pca_harmony',
                        random_state=42)
        sc.tl.leiden(sub, resolution=resolution, key_added='subcluster',
                     random_state=42, flavor=CFG.clustering.leiden_flavor)
        labels = np.array(sub.obs['cell_type'].astype(str)
                          + '_' + sub.obs['subcluster'].astype(str))
        adata.obs.loc[mask, 'cell_type_sub'] = labels.tolist()

    adata.obs['cell_type_sub'] = adata.obs['cell_type_sub'].astype('category')
    n_sub = adata.obs['cell_type_sub'].nunique()
    logger.info("Subclustering complete: %d subtypes", n_sub)


# ═══════════════════════════════════════════════════════════════════════
#  AI 注释函数
# ═══════════════════════════════════════════════════════════════════════

def ai_annotate(adata, CFG, logger, std=None):
    """
    基于 LLM 的 AI 注释主流程。

    步骤:
      1. rank_genes_groups → 获取各聚类 marker 基因
      2. 保存 marker 基因 CSV
      3. 构建提示词 → 调用 LLM
      4. 解析 JSON 响应 → 映射注释到 adata.obs
      5. 生成 UMAP 可视化 & 导出注释表格

    返回:
        annotations dict (解析成功) 或 None (失败，触发回退)
    """
    # ── a. 计算 marker 基因 ───────────────────────────────────────────
    logger.info("Computing marker genes (Wilcoxon rank-sum)...")
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                            use_raw=True if adata.raw is not None else None)

    n_clusters = adata.obs['leiden'].nunique()
    # ── 自适应 max_tokens ───────────────────────────────────────────
    # 聚类数多时 JSON 注释响应会超过默认 4096 token 上限。线性放大
    # （300 token/聚类）但 floor=4096（避免小数据集浪费），cap=32768
    # （避免单次请求成本失控）。`max(...)` 保护用户已显式调高的预算。
    suggested_max_tokens = min(max(4096, n_clusters * 300), 32768)
    CFG.ai.max_tokens = max(
        getattr(CFG.ai, 'max_tokens', 4096), suggested_max_tokens
    )
    logger.info("Adaptive max_tokens: n_clusters=%d → max_tokens=%d",
                n_clusters, CFG.ai.max_tokens)

    compact = n_clusters > 20
    if compact:
        logger.info("n_clusters=%d (>20), using compact prompt mode", n_clusters)

    # ── b. 保存 marker 基因 CSV ───────────────────────────────────────
    marker_rows = []
    for cl in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        df = sc.get.rank_genes_groups_df(adata, group=str(cl))
        df['cluster'] = cl
        marker_rows.append(df)
    marker_df = pd.concat(marker_rows, ignore_index=True)
    marker_csv = os.path.join(CFG.table_dir, 'marker_genes_ai.csv')
    marker_df.to_csv(marker_csv, index=False)
    logger.info("Marker genes saved: %s", marker_csv)

    # ── c. 获取组织 & 物种 ────────────────────────────────────────────
    tissue = CFG.tissue
    species = CFG.species
    logger.info("Annotation context: tissue=%s, species=%s", tissue, species)

    # ── d. 构建提示词 ─────────────────────────────────────────────────
    from core.ai_prompts import build_annotation_prompt
    stages_present = sorted(adata.obs['stage'].unique().tolist()) if 'stage' in adata.obs else []
    extra_context = f"Developmental stages: {stages_present}" if stages_present else ""
    kb_candidates = std.get_candidates() if std else None
    sys_prompt, user_prompt = build_annotation_prompt(adata, tissue, species, precomputed_rank=True, extra_context=extra_context, compact=compact, kb_candidates=kb_candidates)

    # ── e. 调用 LLM ───────────────────────────────────────────────────
    from core.ai_caller import ai_query
    logger.info("Requesting cell type annotation from LLM (model=%s)...", CFG.ai.model)
    try:
        response = ai_query(sys_prompt, user_prompt, cfg=CFG.ai)
    except Exception as exc:
        logger.warning("LLM query failed: %s — falling back to score_genes method", exc)
        return None

    # ── f. 解析 JSON ──────────────────────────────────────────────────
    try:
        annotations = json.loads(response)
    except (json.JSONDecodeError, TypeError) as e:
        logger.warning("LLM response is not valid JSON (%s) — falling back to score_genes method", e)
        logger.warning("Raw response (first 500 chars): %s", response[:500])
        return None

    # 验证每聚类注释结构
    required_keys = {'cell_type', 'state', 'subtype', 'confidence', 'reasoning'}
    for cid, ann in annotations.items():
        if not isinstance(ann, dict):
            logger.warning("Cluster %s annotation is not a dict — falling back", cid)
            return None
        missing = required_keys - ann.keys()
        if missing:
            logger.warning("Cluster %s missing fields %s — falling back", cid, missing)
            return None

    logger.info("LLM annotation successful: %d clusters parsed", len(annotations))

    # ── Standardize cell types if standardizer is active ──────
    if std is not None:
        annotations = std.map_annotations(annotations)
        for cid in annotations:
            raw = annotations[cid].get('cell_type_raw', annotations[cid]['cell_type'])
            _, display_name, _ = std.standardize(raw)
            annotations[cid]['cell_type'] = display_name

    # ── g/h. 映射注释到 adata.obs ─────────────────────────────────────
    leiden_str = adata.obs['leiden'].astype(str)
    adata.obs['cell_type'] = leiden_str.map(
        {k: v['cell_type'] for k, v in annotations.items()}
    ).astype('category')
    adata.obs['cell_state'] = leiden_str.map(
        {k: v['state'] for k, v in annotations.items()}
    )
    adata.obs['cell_subtype'] = leiden_str.map(
        {k: v['subtype'] for k, v in annotations.items()}
    )
    adata.obs['annot_confidence'] = leiden_str.map(
        {k: v['confidence'] for k, v in annotations.items()}
    )
    adata.obs['annot_reasoning'] = leiden_str.map(
        {k: v['reasoning'] for k, v in annotations.items()}
    )
    if std is not None:
        adata.obs['cell_type_std'] = leiden_str.map(
            {k: v.get('cell_type_std', v['cell_type']) for k, v in annotations.items()}
        ).astype('category')
        adata.obs['cell_type_raw'] = leiden_str.map(
            {k: v.get('cell_type_raw', v['cell_type']) for k, v in annotations.items()}
        ).astype('category')
        adata.obs['marker_validation'] = leiden_str.map(
            {k: v.get('marker_validation', 'NO_ONTOLOGY') for k, v in annotations.items()}
        ).astype('category')

    # ── i. 保存注释 CSV ───────────────────────────────────────────────
    ann_records = []
    for cid in sorted(annotations.keys(), key=lambda x: int(x)):
        ann = annotations[cid]
        ann_records.append({
            'cluster': cid,
            'cell_type': ann['cell_type'],
            'cell_type_std': ann.get('cell_type_std', ann['cell_type']),
            'cell_type_raw': ann.get('cell_type_raw', ann['cell_type']),
            'marker_validation': ann.get('marker_validation', 'NO_ONTOLOGY'),
            'state': ann['state'],
            'subtype': ann['subtype'],
            'confidence': ann['confidence'],
            'reasoning': ann['reasoning'],
        })
    ann_df = pd.DataFrame(ann_records)
    ann_csv = os.path.join(CFG.table_dir, 'cell_type_annotations.csv')
    ann_df.to_csv(ann_csv, index=False)
    logger.info("Annotation table saved: %s", ann_csv)

    # 日志输出映射
    logger.info("Cluster → cell type mapping (AI):")
    for rec in ann_records:
        logger.info("  Cluster %s → %s (state=%s, subtype=%s, conf=%s)",
                    rec['cluster'], rec['cell_type'],
                    rec['state'], rec['subtype'], rec['confidence'])

    # ── j. UMAP 可视化 ────────────────────────────────────────────────
    sc.settings.figdir = os.path.join(CFG.figure_dir, '05_annotation')
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False

    # annot_label = cell_type (+ state 如果不为 N/A)
    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)
    state_not_na = adata.obs['cell_state'] != 'N/A'
    adata.obs.loc[state_not_na, 'annot_label'] = (
        adata.obs.loc[state_not_na, 'cell_type'].astype(str)
        + ' (' + adata.obs.loc[state_not_na, 'cell_state'].astype(str) + ')'
    )

    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype_ai.pdf')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label_ai.pdf')

    meta_dict = {
        'barcode': adata.obs_names,
        'UMAP_1': adata.obsm['X_umap'][:, 0],
        'UMAP_2': adata.obsm['X_umap'][:, 1],
        'cell_type': adata.obs['cell_type'].values,
        'cell_state': adata.obs['cell_state'].values,
        'cell_subtype': adata.obs['cell_subtype'].values,
        'annot_confidence': adata.obs['annot_confidence'].values,
    }
    if std is not None:
        meta_dict['cell_type_std'] = adata.obs['cell_type_std'].values
        meta_dict['cell_type_raw'] = adata.obs['cell_type_raw'].values
        meta_dict['marker_validation'] = adata.obs['marker_validation'].values
    meta_df = pd.DataFrame(meta_dict)
    meta_csv = os.path.join(CFG.table_dir, 'cell_metadata.csv')
    meta_df.to_csv(meta_csv, index=False)
    logger.info("Cell metadata exported: %s", meta_csv)

    return annotations



from rna.annotation_engine import run_unified_annotation as unified_annotate


# ═══════════════════════════════════════════════════════════════════════
#  Score_genes 模式 (回退)
# ═══════════════════════════════════════════════════════════════════════

def score_genes_mode(adata, CFG, logger):
    """Score_genes 回退模式 — 复用旧有 run_annotation + run_subclustering。"""
    logger.info("Score_genes mode — marker gene-based scoring annotation")

    run_annotation(adata, CFG.marker.marker_dict, logger)
    run_subclustering(adata, CFG, CFG.marker.subcluster_types,
                      CFG.marker.subcluster_resolution, CFG.marker.min_cells_subcluster, logger)

    # 统一列名: cell_type_sub → cell_subtype
    if 'cell_type_sub' in adata.obs:
        adata.obs['cell_subtype'] = adata.obs['cell_type_sub'].astype(str)

    # annot_label (这里仅为 cell_type，无 state 信息)
    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)

    # 置信度重命名
    if 'annotation_confidence' in adata.obs:
        adata.obs['annot_confidence'] = adata.obs['annotation_confidence']

    # 可视化
    sc.settings.figdir = os.path.join(CFG.figure_dir, '05_annotation')
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False
    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype.pdf')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label.pdf')
    if 'annotation_confidence' in adata.obs:
        safe_plot(sc.pl.umap, adata, color='annotation_confidence', show=False,
                  save='_05_confidence.pdf', cmap='viridis')

    # 导出细胞元数据
    meta_cols = ['barcode']
    if 'X_umap' in adata.obsm:
        meta_df = pd.DataFrame({
            'barcode': adata.obs_names,
            'UMAP_1': adata.obsm['X_umap'][:, 0],
            'UMAP_2': adata.obsm['X_umap'][:, 1],
        })
    else:
        meta_df = pd.DataFrame({'barcode': adata.obs_names})
    for col in ['cell_type', 'cell_subtype', 'cell_type_sub', 'annotation_confidence']:
        if col in adata.obs:
            meta_df[col] = adata.obs[col].values
    meta_csv = os.path.join(CFG.table_dir, 'cell_metadata.csv')
    meta_df.to_csv(meta_csv, index=False)
    logger.info("Cell metadata exported: %s", meta_csv)


# ═══════════════════════════════════════════════════════════════════════
#  主入口
# ═══════════════════════════════════════════════════════════════════════

def main():
    global log
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("05_annotate_major",
                        os.path.join(CFG.log_dir, "05_annotate_major.log"))

    # ── Annotation Standardizer ────────────────────────────────
    standardizer = getattr(CFG, 'tissue_ontology', None) or CFG.tissue_kb
    std = None
    if standardizer:
        from rna.annotation_standardizer import StandardOntology
        try:
            std = StandardOntology(standardizer)
            log.info("Annotation Standardizer active for tissue: %s", standardizer)
        except NotImplementedError:
            log.warning("Annotation Standardizer not available for tissue: %s", standardizer)
    log.info("Step 05: Cell type annotation (Major Lineage)")

    adata = sc.read(CFG.cluster_h5ad)
    log.info("Loaded: %s — %d cells, %d clusters",
             CFG.cluster_h5ad, adata.n_obs, adata.obs['leiden'].nunique())

    # ── 判断 AI 模式/Unified KB 模式是否可用 ────────────────────────────
    ai_enabled = getattr(CFG.ai, 'enabled', False)
    ai_annot_on = getattr(CFG.ai, 'ai_annotation', False)

    # ── Unified KB mode (if tissue_kb is set) ────────────────────────────
    if CFG.tissue_kb:
        log.info("Unified KB mode enabled \u2014 tissue_kb='%s'", CFG.tissue_kb)
        ann_result = unified_annotate(adata, CFG, log)
        if ann_result is not None:
            if std is not None:
                validation_results = std.validate(
                    adata,
                    top_n=CFG.marker.validation_n_top_genes,
                    min_overlap=CFG.marker.validation_min_overlap,
                    marginal_threshold=CFG.marker.validation_marginal_threshold,
                )
                log.info("Marker validation: %d/%d PASS",
                         sum(1 for r in validation_results if r['status'] == 'PASS'),
                         len(validation_results))
                validation_map = {r['cluster']: r['status'] for r in validation_results}
                adata.obs['marker_validation'] = adata.obs['leiden'].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
            _warn_if_low_coverage(adata, log)
            safe_write(adata, CFG.annotated_h5ad, cfg=CFG)
            log.info("Step 05 (Unified mode) complete, took %.1fs", time.time() - t0)
            return
        log.warning("Unified annotation failed, falling back to Score_genes mode")
    elif ai_enabled and ai_annot_on:
        log.info("AI mode enabled \u2014 using LLM for smart annotation")
        ann_result = ai_annotate(adata, CFG, log, std=std)
        if ann_result is not None:
            if std is not None:
                validation_results = std.validate(
                    adata,
                    top_n=CFG.marker.validation_n_top_genes,
                    min_overlap=CFG.marker.validation_min_overlap,
                    marginal_threshold=CFG.marker.validation_marginal_threshold,
                )
                log.info("Marker validation: %d/%d PASS",
                         sum(1 for r in validation_results if r['status'] == 'PASS'),
                         len(validation_results))
                validation_map = {r['cluster']: r['status'] for r in validation_results}
                adata.obs['marker_validation'] = adata.obs['leiden'].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
            _warn_if_low_coverage(adata, log)
            safe_write(adata, CFG.annotated_h5ad, cfg=CFG)
            log.info("Step 05 (AI mode) complete, took %.1fs", time.time() - t0)
            return
        log.warning("AI annotation failed, falling back to Score_genes mode")

    # ── Score_genes \u6a21\u5f0f (\u6240\u6709\u8def\u5f84\u56de\u9000) ─────────────────────────────────
    score_genes_mode(adata, CFG, log)
    if std is not None:
        validation_results = std.validate(
            adata,
            top_n=CFG.marker.validation_n_top_genes,
            min_overlap=CFG.marker.validation_min_overlap,
            marginal_threshold=CFG.marker.validation_marginal_threshold,
        )
        log.info("Marker validation: %d/%d PASS",
                 sum(1 for r in validation_results if r['status'] == 'PASS'),
                 len(validation_results))
        validation_map = {r['cluster']: r['status'] for r in validation_results}
        adata.obs['marker_validation'] = adata.obs['leiden'].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
        _warn_if_low_coverage(adata, log)
    safe_write(adata, CFG.annotated_h5ad, cfg=CFG)
    log.info("Step 05 (score_genes mode) complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
