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
# Also add rna/ so `from tissue_ontologies import load_kb` resolves
_rna_pkg = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
if _rna_pkg not in sys.path:
    sys.path.insert(0, _rna_pkg)
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np
import logging

log: logging.Logger


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
    groupby_kw = {'observed': True} if hasattr(pd.Categorical, 'observed') else {}
    cluster_scores = adata.obs.groupby('leiden', **groupby_kw)[score_cols].mean()
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
                     random_state=42, flavor=CFG.leiden_flavor)
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
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

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
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False

    # annot_label = cell_type (+ state 如果不为 N/A)
    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)
    state_not_na = adata.obs['cell_state'] != 'N/A'
    adata.obs.loc[state_not_na, 'annot_label'] = (
        adata.obs.loc[state_not_na, 'cell_type'].astype(str)
        + ' (' + adata.obs.loc[state_not_na, 'cell_state'].astype(str) + ')'
    )

    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype_ai.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label_ai.pdf', legend_loc='on data')

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



def unified_annotate(adata, CFG, logger):
    """
    KB-based unified annotation mode.

    Uses marker scoring + expert rules + evidence fusion, with AI fallback
    for low-confidence clusters.

    Returns
    -------
    dict or None
        Cluster -> FusionDecision mapping, or None on failure (triggers fallback).
    """
    # ── a. Compute marker genes ───────────────────────────────────────────
    logger.info("Computing marker genes (Wilcoxon rank-sum)...")
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

    n_clusters = adata.obs['leiden'].nunique()

    # ── b. Build per-cluster marker DataFrames ────────────────────────────
    marker_rows = []
    for cl in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        df = sc.get.rank_genes_groups_df(adata, group=str(cl))
        df['cluster'] = cl
        marker_rows.append(df)
    marker_df = pd.concat(marker_rows, ignore_index=True)
    marker_csv = os.path.join(CFG.table_dir, 'marker_genes_unified.csv')
    marker_df.to_csv(marker_csv, index=False)
    logger.info("Marker genes saved: %s", marker_csv)

    # ── c. Load KB ────────────────────────────────────────────────────────
    from tissue_ontologies import load_kb
    try:
        kb = load_kb(CFG.tissue_kb)
    except Exception as exc:
        logger.warning("Failed to load KB '%s': %s", CFG.tissue_kb, exc)
        return None

    n_types = sum(1 for k in kb if k != 'expert_rules')
    n_rules = len(kb.get('expert_rules', []))
    logger.info("Loaded KB: %s (%d cell types, %d rules)",
                CFG.tissue_kb, n_types, n_rules)

    # ── d. Full marker scoring + expert rules per cluster ─────────────────
    # Fix import shadowing: scripts/ at sys.path[0] shadows root/utils/ package
    import sys as _sys_module
    _scripts_dir = os.path.dirname(os.path.abspath(__file__))
    _orig_sys_path = list(_sys_module.path)
    _sys_module.path = [p for p in _sys_module.path if p != _scripts_dir]
    # Also clean any wrongly-resolved 'utils' from sys.modules
    _sys_module.modules.pop('utils', None)
    for _k in list(_sys_module.modules.keys()):
        if _k.startswith('utils.'):
            _sys_module.modules.pop(_k, None)
    try:
        from rna.utils.marker_scoring import (
            score_cluster_against_kb, apply_expert_rules,
            annotate_all_clusters, resolve_expert_rule_params,
            detect_low_quality_cluster,
        )
        from rna.utils.evidence_fusion import fuse_all_clusters
    finally:
        _sys_module.path = _orig_sys_path

    species = CFG.species

    # ── Resolve taxonomic class/order for phylogenetic weighting ──────
    # CFG.target_class/order take precedence; fall back to species lookup.
    from rna.ortholog import SPECIES_TO_CLASS
    target_class = CFG.target_class or SPECIES_TO_CLASS.get(species, "")
    target_order = CFG.target_order or ""
    if target_class:
        logger.info("Phylogenetic weighting: target_class=%s, target_order=%s",
                     target_class, target_order or "(none)")
    else:
        logger.info("Phylogenetic weighting: disabled (no target_class for '%s')",
                     species)

    # ── Resolve expert-rule constraint parameters ────────────────────
    rule_top_n, rule_pval = resolve_expert_rule_params(
        strictness=getattr(CFG, 'expert_rule_strictness', 'default'),
        top_n=getattr(CFG, 'expert_rule_top_n', 0),
        pval_cutoff=getattr(CFG, 'expert_rule_pval_cutoff', 0.0),
    )
    logger.info(
        "Expert rules: strictness=%s → top_n=%d, pval_cutoff=%.3f",
        getattr(CFG, 'expert_rule_strictness', 'default'),
        rule_top_n, rule_pval,
    )

    # ── Compute per-cluster marker scores and expert rules ───────────
    all_scores = {}
    all_rules = {}
    low_quality_clusters: dict[str, str] = {}  # cluster_str → reason
    clusters = sorted(
        marker_df['cluster'].unique(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )
    for cl in clusters:
        cl_str = str(cl)
        cl_mask = marker_df['cluster'] == cl
        cl_data = marker_df[cl_mask].copy()
        lfc_idx = cl_data['logfoldchanges'].argsort()[::-1]
        cl_data = cl_data.iloc[lfc_idx]

        # Path C: detect low-quality clusters (mito/ribo dominated)
        is_lq, lq_reason = detect_low_quality_cluster(cl_data)
        if is_lq:
            low_quality_clusters[cl_str] = lq_reason
            logger.info("Cluster %s flagged as low-quality: %s", cl_str, lq_reason)

        all_scores[cl_str] = score_cluster_against_kb(
            kb, cl_data, species=species,
            target_class=target_class, target_order=target_order,
            adaptive_top_n=True,
        )
        all_rules[cl_str] = apply_expert_rules(kb, cl_data,
                                                top_n=rule_top_n,
                                                pval_cutoff=rule_pval)

    decisions, fusion_quality = fuse_all_clusters(
        all_scores, all_rules, kb=kb, all_marker_dfs=marker_df,
        return_quality=True,
        low_quality_clusters=low_quality_clusters,
        unconstrained=getattr(CFG.ai, 'unconstrained_annotation', False),
    )
    logger.info("Evidence fusion: %d clusters processed", len(decisions))

    if not decisions:
        logger.warning("Evidence fusion produced no decisions — falling back")
        return None

    # Build cluster -> decision mapping (preserving fusion sort order)
    decision_clusters = sorted(
        all_scores.keys(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )
    decision_map = dict(zip(decision_clusters, decisions))

    # ── f. AI fallback for low-confidence clusters ────────────────────────
    ai_enabled = getattr(CFG.ai, 'enabled', False)
    ai_annot_on = getattr(CFG.ai, 'ai_annotation', False)

    low_conf_clusters = [
        d for d in decisions if d.confidence in ('low', 'unknown')
    ]
    ai_results = {}

    if low_conf_clusters and ai_enabled and ai_annot_on:
        logger.info(
            "AI fallback for %d low-confidence clusters", len(low_conf_clusters)
        )
        kb_candidates = sorted([k for k in kb if k != 'expert_rules' and not k.startswith('_')])
        tissue = CFG.tissue
        stages_present = (
            sorted(adata.obs['stage'].unique().tolist())
            if 'stage' in adata.obs else []
        )
        extra_context = (
            f"Developmental stages: {stages_present}" if stages_present else ""
        )

        # Unconstrained annotations require build_annotation_prompt import here
        from core.ai_prompts import build_annotation_prompt
        from core.ai_caller import ai_query

        unconstrained = getattr(CFG.ai, 'unconstrained_annotation', False)
        sys_prompt, user_prompt = build_annotation_prompt(
            adata, tissue, species, precomputed_rank=True,
            extra_context=extra_context,
            compact=n_clusters > 20,
            kb_candidates=kb_candidates,
            unconstrained=unconstrained,
        )

        try:
            response = ai_query(sys_prompt, user_prompt, cfg=CFG.ai)
            if response:
                ai_parsed = json.loads(response)
                for cid, ann in ai_parsed.items():
                    if isinstance(ann, dict) and 'cell_type' in ann:
                        ai_results[str(cid)] = ann['cell_type']
                logger.info(
                    "AI fallback: %d cluster suggestions received", len(ai_results)
                )
                # Re-run fusion with AI context
                decisions = fuse_all_clusters(
                    all_scores, all_rules, kb=kb,
                    all_marker_dfs=marker_df,
                    ai_results=ai_results,
                    low_quality_clusters=low_quality_clusters,
                    unconstrained=getattr(CFG.ai, 'unconstrained_annotation', False),
                )
                decision_map = dict(zip(decision_clusters, decisions))
        except Exception as exc:
            logger.warning(
                "AI fallback failed: %s \u2014 using pure KB results", exc
            )

    # ── g. Map decisions to adata.obs ─────────────────────────────────────
    leiden_str = adata.obs['leiden'].astype(str)

    # For low-quality clusters, downgrade: force Unknown + annotate reason.
    _forced_unknown = 0
    for cl_str, reason in low_quality_clusters.items():
        if cl_str in decision_map and decision_map[cl_str].confidence != 'unknown':
            old_ct = adata.obs['cell_type'].iloc[0] if cl_str == '0' else '—'
            decision_map[cl_str] = decision_map[cl_str]._replace(
                cell_type='Unknown',
                confidence='unknown',
                method='unknown',
                explanation=(
                    "Low-quality cluster ({}) — {}"
                    .format(reason, decision_map[cl_str].explanation[:120])
                ),
            )
            _forced_unknown += 1
    if _forced_unknown:
        logger.info(
            "Downgraded %d low-quality cluster(s) to Unknown: %s",
            _forced_unknown,
            ", ".join(
                "{} ({})".format(k, v) for k, v in low_quality_clusters.items()
                if k in decision_map
            ),
        )

    adata.obs['cell_type'] = leiden_str.map(
        {k: v.cell_type for k, v in decision_map.items()}
    ).astype('category')
    adata.obs['cell_state'] = leiden_str.map(
        {k: 'N/A' for k in decision_map}
    )
    adata.obs['cell_subtype'] = leiden_str.map(
        {k: 'N/A' for k in decision_map}
    )

    # annot_method: clean label from fusion method (+ AI suffix)
    def _clean_method(d):
        if d.method == 'expert_rule':
            return 'expert_rule'
        if d.method == 'unknown':
            return 'unknown'
        if d.ai_agreed or d.ai_suggested:
            return 'marker_scoring+ai'
        return 'marker_scoring'

    adata.obs['annot_method'] = leiden_str.map(
        {k: _clean_method(v) for k, v in decision_map.items()}
    )
    adata.obs['annot_confidence'] = leiden_str.map(
        {k: v.confidence for k, v in decision_map.items()}
    )
    adata.obs['annot_reasoning'] = leiden_str.map(
        {k: v.explanation for k, v in decision_map.items()}
    )
    adata.obs['annot_evidence'] = leiden_str.map(
        {k: json.dumps({
            'score': v.score,
            'method': v.method,
            'n_markers_found': v.n_markers_found,
            'ai_agreed': v.ai_agreed,
            'ai_suggested': v.ai_suggested,
            'diagnostic_category': v.diagnostic.category if v.diagnostic else None,
            'diagnostic_detail': v.diagnostic.detail if v.diagnostic else None,
            'top_competitors': v.diagnostic.top_competitors if v.diagnostic else [],
        }) for k, v in decision_map.items()}
    )

    # ── h. Save annotation CSV ────────────────────────────────────────────
    ann_records = []
    sort_key = lambda x: int(x) if str(x).isdigit() else str(x)
    for cl_name in sorted(decision_map.keys(), key=sort_key):
        d = decision_map[cl_name]
        ann_records.append({
            'cluster': cl_name,
            'cell_type': d.cell_type,
            'confidence': d.confidence,
            'method': _clean_method(d),
            'score': d.score,
            'n_markers_found': d.n_markers_found,
            'ai_agreed': d.ai_agreed,
            'ai_suggested': d.ai_suggested,
            'reasoning': d.explanation,
            'diagnostic_category': d.diagnostic.category if d.diagnostic else '',
        })
    ann_df = pd.DataFrame(ann_records)
    ann_csv = os.path.join(CFG.table_dir, 'cell_type_annotations.csv')
    ann_df.to_csv(ann_csv, index=False)
    logger.info("Annotation table saved: %s", ann_csv)

    logger.info("Cluster \u2192 cell type mapping (Unified):")
    for rec in ann_records:
        logger.info(
            "  Cluster %s \u2192 %s (conf=%s, method=%s)",
            rec['cluster'], rec['cell_type'],
            rec['confidence'], rec['method'],
        )

    # ── i. UMAP visualization ─────────────────────────────────────────────
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False

    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)

    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype_unified.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label_unified.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annot_confidence', show=False,
              save='_05_confidence_unified.pdf')

    # ── j. Cell metadata export ───────────────────────────────────────────
    meta_df = pd.DataFrame({
        'barcode': adata.obs_names,
        'UMAP_1': adata.obsm['X_umap'][:, 0],
        'UMAP_2': adata.obsm['X_umap'][:, 1],
        'cell_type': adata.obs['cell_type'].values,
        'cell_state': adata.obs['cell_state'].values,
        'cell_subtype': adata.obs['cell_subtype'].values,
        'annot_confidence': adata.obs['annot_confidence'].values,
        'annot_method': adata.obs['annot_method'].values,
    })
    meta_csv = os.path.join(CFG.table_dir, 'cell_metadata.csv')
    meta_df.to_csv(meta_csv, index=False)
    logger.info("Cell metadata exported: %s", meta_csv)

    # ── k. Annotation quality report ────────────────────────────────────────
    _write_quality_report(adata, ann_records, fusion_quality, CFG, logger)

    # ── l. Interactive review (--interactive flag) ──────────────────────────
    if getattr(CFG, 'interactive', False):
        _interactive_annotation_review(adata, fusion_quality, CFG, logger)

    return decision_map


def _write_quality_report(adata, ann_records, fusion_quality, CFG, logger):
    """Write 05_annotation_quality.json summarising annotation health."""
    import json as _json

    pass_cells = (
        (adata.obs['marker_validation'] == 'PASS').sum()
        if 'marker_validation' in adata.obs else 0
    )
    pass_rate = pass_cells / max(adata.n_obs, 1)

    ambiguity_clusters = []
    for rec in ann_records:
        reasoning = rec.get('reasoning', '')
        if 'also matched rules:' in reasoning:
            ambiguity_clusters.append(rec['cluster'])

    quality = {
        "pass_rate": round(pass_rate, 4),
        "total_clusters": len(ann_records),
        "annotated_by_rule": fusion_quality.get("annotated_by_rule", 0),
        "annotated_by_scoring": fusion_quality.get("annotated_by_scoring", 0),
        "unknown": fusion_quality.get("unknown", 0),
        "ambiguity_clusters": ambiguity_clusters,
        "ai_disagreement_rate": round(
            sum(1 for r in ann_records if not r.get('ai_agreed', True))
            / max(len(ann_records), 1), 4,
        ),
        "kb_blind_spot": pass_rate < 0.1,
        "recommended_strictness": (
            "relaxed" if pass_rate < 0.1 else
            "deep" if pass_rate < 0.3 else
            "default"
        ),
    }

    quality_path = os.path.join(CFG.table_dir, '05_annotation_quality.json')
    with open(quality_path, 'w', encoding='utf-8') as f:
        _json.dump(quality, f, indent=2)
    logger.info(
        "Annotation quality report: %s (pass_rate=%.1f%%)",
        quality_path, quality["pass_rate"] * 100,
    )


def _interactive_annotation_review(adata, fusion_quality, CFG, logger):
    """Present annotation quality summary and offer remediation choices.

    Only called when ``CFG.interactive`` is ``True``.
    """
    pass_cells = (
        (adata.obs['marker_validation'] == 'PASS').sum()
        if 'marker_validation' in adata.obs else 0
    )
    pass_rate = pass_cells / max(adata.n_obs, 1)
    n_total = fusion_quality.get("total", 0)
    n_rule = fusion_quality.get("annotated_by_rule", 0)
    n_scoring = fusion_quality.get("annotated_by_scoring", 0)
    n_unknown = fusion_quality.get("unknown", 0)
    n_ambiguity = fusion_quality.get("ambiguity", 0)

    print("\n" + "=" * 60)
    print("Annotation Quality Summary")
    print("=" * 60)
    print(f"  PASS rate:          {pass_rate * 100:.1f}%")
    print(f"  Annotated:          {n_rule} by rule, {n_scoring} by scoring")
    print(f"  Unknown:            {n_unknown}/{n_total}")
    if n_ambiguity > 0:
        print(f"  ⚠  High ambiguity:  {n_ambiguity} cluster(s) matched ≥3 rules")
    if pass_rate < 0.1:
        print(f"  ⚠  KB blind spot detected")
        rec = "relaxed" if pass_rate < 0.1 else "deep" if pass_rate < 0.3 else "default"
        print(f"  \U0001f4a1 Recommended:       strictness='{rec}'")
    print()

    if pass_rate < 0.1:
        rec = "relaxed" if pass_rate < 0.1 else "deep"
        try:
            choice = input(
                "KB coverage is very low on this dataset. Options:\n"
                f"  [r] Re-annotate with strictness='{rec}'\n"
                "  [s] Continue with score_genes fallback\n"
                "  [c] Continue with current labels (not recommended)\n"
                "  [a] Abort\n"
                "Choice> "
            ).strip().lower()
        except (EOFError, KeyboardInterrupt):
            logger.warning("Interactive input interrupted — continuing")
            return

        if choice == 'r':
            logger.info(
                "User chose: re-annotate with strictness='%s'", rec,
            )
            print(f"\nTo re-annotate, set in your config:\n"
                  f"  CFG.expert_rule_strictness = '{rec}'\n"
                  f"Or pass --config with the updated setting.\n")
        elif choice == 's':
            logger.info("User chose: score_genes fallback")
            print("Set CFG.tissue_kb = '' to use score_genes mode.\n")
        elif choice == 'a':
            raise SystemExit("Aborted by user.")
        else:
            logger.info("User chose: continue with current labels")


# ═══════════════════════════════════════════════════════════════════════
#  Score_genes 模式 (回退)
# ═══════════════════════════════════════════════════════════════════════

def score_genes_mode(adata, CFG, logger):
    """Score_genes 回退模式 — 复用旧有 run_annotation + run_subclustering。"""
    logger.info("Score_genes mode — marker gene-based scoring annotation")

    run_annotation(adata, CFG.marker_dict, logger)
    run_subclustering(adata, CFG, CFG.subcluster_types,
                      CFG.subcluster_resolution, CFG.min_cells_subcluster, logger)

    # 统一列名: cell_type_sub → cell_subtype
    if 'cell_type_sub' in adata.obs:
        adata.obs['cell_subtype'] = adata.obs['cell_type_sub'].astype(str)

    # annot_label (这里仅为 cell_type，无 state 信息)
    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)

    # 置信度重命名
    if 'annotation_confidence' in adata.obs:
        adata.obs['annot_confidence'] = adata.obs['annotation_confidence']

    # 可视化
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False
    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype.pdf', legend_loc='on data')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label.pdf', legend_loc='on data')
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
                    top_n=CFG.marker_validation_n_top_genes,
                    min_overlap=CFG.marker_validation_min_overlap,
                    marginal_threshold=CFG.marker_validation_marginal_threshold,
                )
                log.info("Marker validation: %d/%d PASS",
                         sum(1 for r in validation_results if r['status'] == 'PASS'),
                         len(validation_results))
                validation_map = {r['cluster']: r['status'] for r in validation_results}
                adata.obs['marker_validation'] = adata.obs['leiden'].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
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
                    top_n=CFG.marker_validation_n_top_genes,
                    min_overlap=CFG.marker_validation_min_overlap,
                    marginal_threshold=CFG.marker_validation_marginal_threshold,
                )
                log.info("Marker validation: %d/%d PASS",
                         sum(1 for r in validation_results if r['status'] == 'PASS'),
                         len(validation_results))
                validation_map = {r['cluster']: r['status'] for r in validation_results}
                adata.obs['marker_validation'] = adata.obs['leiden'].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
            safe_write(adata, CFG.annotated_h5ad, cfg=CFG)
            log.info("Step 05 (AI mode) complete, took %.1fs", time.time() - t0)
            return
        log.warning("AI annotation failed, falling back to Score_genes mode")

    # ── Score_genes \u6a21\u5f0f (\u6240\u6709\u8def\u5f84\u56de\u9000) ─────────────────────────────────
    score_genes_mode(adata, CFG, log)
    if std is not None:
        validation_results = std.validate(
            adata,
            top_n=CFG.marker_validation_n_top_genes,
            min_overlap=CFG.marker_validation_min_overlap,
            marginal_threshold=CFG.marker_validation_marginal_threshold,
        )
        log.info("Marker validation: %d/%d PASS",
                 sum(1 for r in validation_results if r['status'] == 'PASS'),
                 len(validation_results))
        validation_map = {r['cluster']: r['status'] for r in validation_results}
        adata.obs['marker_validation'] = adata.obs['leiden'].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
    safe_write(adata, CFG.annotated_h5ad)
    log.info("Step 05 (score_genes mode) complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
