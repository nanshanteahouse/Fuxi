"""
Unified annotation engine for cell type annotation.

Provides ``run_unified_annotation`` which performs KB-based unified annotation
with marker scoring, expert rules, evidence fusion, and optional AI fallback.

Extracted from ``rna/steps/05_annotate_major.py`` for cross-module reuse
(RNA, spatial, ATAC pipelines).
"""

import os
import json
import scanpy as sc
import pandas as pd
from core.utils import safe_plot
from collections import Counter




def _check_zero_scores_and_retry(
    kb, all_scores, marker_df, clusters, species,
    target_class, target_order, tissue_kb, logger,
):
    """Check for zero KB marker scores and attempt case-insensitive retry.

    Builds an uppercased deepcopy of the KB once and reuses it across two
    retry passes (the original KB is never mutated):

    1. Global retry — fires when *all* clusters have zero hits and there
       are >= 5 clusters.  Re-scores every cluster; replaces ``all_scores``
       if any cluster improves.

    2. Per-cluster retry — fires when any clusters remain zero-hit and
       there are >= 3 clusters.  Re-scores individual zero-hit clusters
       and updates them in place when their hit count improves.

    A diagnostic ERROR fires if zero hits persist across all clusters
    after retry.

    Parameters
    ----------
    kb : dict
        Original KB (never mutated).
    all_scores : dict
        ``{cluster_str: {type_key: Score}}`` from the initial scoring pass.
    marker_df : pd.DataFrame
        Concatenated marker DataFrames with ``names``, ``logfoldchanges``,
        and ``cluster`` columns.
    clusters : list
        Sorted cluster identifiers.
    species : str
        Species filter for KB scoring.
    target_class : str
        Taxonomic class for phylogenetic weighting.
    target_order : str
        Taxonomic order for phylogenetic weighting.
    tissue_kb : str
        KB identifier for diagnostic messages.
    logger : logging.Logger

    Returns
    -------
    tuple[dict, int, int]
        ``(all_scores, total_hits, n_clusters_total)`` — ``all_scores`` may
        be updated if the retry improved hits.
    """
    total_hits = sum(
        score.n_markers_found
        for cluster_scores in all_scores.values()
        for score in cluster_scores.values()
    )
    n_clusters_total = len(all_scores)

    # Detect retry eligibility up front so we can build the uppercased KB
    # copy at most once and reuse it across both retry paths.
    zero_hit_clusters = [
        cl_str for cl_str, cs in all_scores.items()
        if sum(s.n_markers_found for s in cs.values()) == 0
    ]
    needs_global_retry = total_hits == 0 and n_clusters_total >= 5
    needs_per_cluster_retry = bool(zero_hit_clusters) and n_clusters_total >= 3

    if needs_global_retry or needs_per_cluster_retry:
        import copy
        from rna.utils.marker_scoring import score_cluster_against_kb

        # Build uppercased KB deepcopy once; reused by both retry passes.
        kb_copy = copy.deepcopy(kb)
        for type_key in kb_copy:
            if type_key == "expert_rules" or type_key.startswith("_"):
                continue
            markers = kb_copy[type_key].get("markers", {})
            for tier in ("confirm", "add"):
                old_genes = markers.get(tier, {})
                if old_genes:
                    markers[tier] = {g.upper(): v for g, v in old_genes.items()}

        def _retry_cluster_score(cl_str):
            """Re-score a single cluster against the uppercased KB copy."""
            cl = int(cl_str) if cl_str.lstrip('-').isdigit() else cl_str
            cl_mask = marker_df['cluster'] == cl
            cl_data = marker_df[cl_mask].copy()
            lfc_idx = cl_data['logfoldchanges'].argsort()[::-1]
            cl_data = cl_data.iloc[lfc_idx]
            cl_data['names'] = cl_data['names'].str.upper()
            return score_cluster_against_kb(
                kb_copy, cl_data, species=species,
                target_class=target_class, target_order=target_order,
                adaptive_top_n=True,
            )

        # Global retry: all clusters zero-hit (n_clusters_total >= 5).
        # Re-score every cluster; replace all_scores if any improvement.
        if needs_global_retry:
            logger.warning(
                "Zero KB marker hits across all %d clusters — "
                "attempting case-insensitive normalization retry",
                n_clusters_total,
            )
            retry_scores = {str(cl): _retry_cluster_score(str(cl)) for cl in clusters}
            retry_total_hits = sum(
                s.n_markers_found
                for cluster_scores in retry_scores.values()
                for s in cluster_scores.values()
            )
            if retry_total_hits > total_hits:
                logger.info(
                    "Case-insensitive retry improved hits from %d to %d — "
                    "using normalized results",
                    total_hits, retry_total_hits,
                )
                all_scores = retry_scores
                total_hits = retry_total_hits
                # Recompute zero-hit list for the per-cluster pass below.
                zero_hit_clusters = [
                    cl_str for cl_str, cs in all_scores.items()
                    if sum(s.n_markers_found for s in cs.values()) == 0
                ]

        # Per-cluster retry: any remaining zero-hit clusters (n_clusters_total >= 3).
        # Updates individual clusters only when retry improves their hit count.
        if zero_hit_clusters and n_clusters_total >= 3:
            logger.warning(
                "%d/%d clusters have zero KB hits — retrying per-cluster "
                "with case-insensitive normalization",
                len(zero_hit_clusters), n_clusters_total,
            )
            per_cluster_fixed = 0
            for cl_str in zero_hit_clusters:
                retry_result = _retry_cluster_score(cl_str)
                retry_hits = sum(
                    s.n_markers_found for s in retry_result.values()
                )
                if retry_hits > 0:
                    all_scores[cl_str] = retry_result
                    per_cluster_fixed += 1
                    total_hits += retry_hits
            if per_cluster_fixed:
                logger.info(
                    "Per-cluster retry recovered %d/%d zero-hit clusters",
                    per_cluster_fixed, len(zero_hit_clusters),
                )
    # Zero-score diagnostic warning (fires regardless of retry outcome)
    if total_hits == 0 and n_clusters_total >= 5:
        logger.error(
            "Zero KB marker hits across all %d clusters for KB '%s' "
            "(species=%s) — likely species mismatch or missing cell types. "
            "Consider expanding KB coverage or using AI fallback.",
            n_clusters_total, tissue_kb, species,
        )

    return all_scores, total_hits, n_clusters_total

CATEGORY_PREFIX = "Broad_"


def _classify_broad_category(cluster_scores: dict, kb: dict) -> str:
    """Classify a cluster into a broad category using Fisher scores.

    Filters cluster_scores to only keys starting with CATEGORY_PREFIX
    (Broad_Progenitor, Broad_Neuron, Broad_Glia, Broad_Non-neural),
    then returns the category name with the highest Fisher score.

    Args:
        cluster_scores: Dict of {type_key: Score} from score_cluster_against_kb()
        kb: The full KB dict (unused; required for signature consistency)

    Returns:
        Category name string (e.g., "Broad_Neuron") or "" if no match.
    """
    from rna.utils.evidence_fusion import _resolve_score

    broad_entries = {
        k: v for k, v in cluster_scores.items()
        if k.startswith(CATEGORY_PREFIX)
    }
    if not broad_entries:
        return ""

    best_type = max(broad_entries, key=lambda k: _resolve_score(cluster_scores, k)[0])
    return best_type

def run_unified_annotation(adata, CFG, logger):
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
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                            use_raw=True if adata.raw is not None else None)
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
    from rna.tissue_ontologies import load_kb
    try:
        kb = load_kb(CFG.tissue_kb)
    except Exception as exc:
        logger.warning("Failed to load KB '%s': %s", CFG.tissue_kb, exc)
        return None

    n_types = sum(
        1 for k in kb
        if k != "expert_rules"
        and not k.startswith("_")
        and not k.startswith(CATEGORY_PREFIX)
    )
    n_rules = len(kb.get('expert_rules', []))
    logger.info("Loaded KB: %s (%d cell types, %d rules)",
                CFG.tissue_kb, n_types, n_rules)

    # ── d. Full marker scoring + expert rules per cluster ─────────────────
    from rna.utils.marker_scoring import (
        score_cluster_against_kb,
        annotate_all_clusters,
        detect_low_quality_cluster,
    )
    from rna.utils.marker_expert_rules import (
        apply_expert_rules, resolve_expert_rule_params,
    )
    from rna.utils.evidence_fusion import fuse_all_clusters

    species = CFG.species
    # ── Normalise species key (e.g. "danio_rerio" → "zebrafish") ────
    # _SPECIES_NORMALISE maps NCBI-style names to pipeline keys; this
    # ensures consistency with _SPECIES_SYNONYMS and SPECIES_TO_CLASS.
    from core.preprocess.format_detector import _SPECIES_NORMALISE
    species = _SPECIES_NORMALISE.get(species, species)

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
    _strictness = getattr(CFG.marker, 'expert_rule_strictness', 'default')
    _top_n = getattr(CFG.marker, 'expert_rule_top_n', 0)
    _pval = getattr(CFG.marker, 'expert_rule_pval_cutoff', 0.0)

    # ── developmental_mode: auto-relax constraints for developing tissues ──
    _dev_mode = getattr(CFG.marker, 'developmental_mode', False)
    if _dev_mode:
        if _strictness == 'default':
            _strictness = 'deep'
            logger.info("developmental_mode: auto-adjusted strictness from 'default' → 'deep'")
        logger.info("developmental_mode: enabled — relaxed constraints for developing tissue")

    rule_top_n, rule_pval = resolve_expert_rule_params(
        strictness=_strictness,
        top_n=_top_n,
        pval_cutoff=_pval,
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

    # ── e. Zero-score detection and case-insensitive retry ──────────────
    all_scores, total_hits, n_clusters_total = _check_zero_scores_and_retry(
        kb, all_scores, marker_df, clusters, species,
        target_class, target_order, CFG.tissue_kb, logger,
    )

    # Pre-filter: strip Broad_* synthetic types from scoring before fusion.
    # Fine-grained cell_type annotation must never see Broad_* entries.
    # The original all_scores (with Broad_* entries) is preserved for
    # _classify_broad_category() below.
    fine_scores = {}
    for cl_str, scores in all_scores.items():
        fine_scores[cl_str] = {
            k: v for k, v in scores.items()
            if not k.startswith(CATEGORY_PREFIX)
        }

    decisions, fusion_quality = fuse_all_clusters(
        fine_scores, all_rules, kb=kb, all_marker_dfs=marker_df,
        return_quality=True,
        low_quality_clusters=low_quality_clusters,
        unconstrained=getattr(CFG.ai, 'unconstrained_annotation', False),
    )
    logger.info("Evidence fusion: %d clusters processed", len(decisions))

    # Build cell_category_map using ORIGINAL all_scores (with Broad_* entries)
    cell_category_map = {}
    for cl_str in all_scores:
        category = _classify_broad_category(all_scores[cl_str], kb)
        cell_category_map[cl_str] = category


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
        from core.ai_prompts import build_annotation_prompt, build_hierarchical_annotation_prompt
        from core.ai_caller import ai_query

        unconstrained = getattr(CFG.ai, 'unconstrained_annotation', False)

        # Use hierarchical prompt if KB has _hierarchy section
        kb_hierarchy = kb.get("_hierarchy") if kb else None
        if kb_hierarchy:
            sys_prompt, user_prompt = build_hierarchical_annotation_prompt(
                adata=adata,
                tissue=tissue,
                species=species,
                kb_candidates=kb_candidates,
                kb_hierarchy=kb_hierarchy,
                precomputed_rank=True,
                extra_context=extra_context,
                compact=n_clusters > 20,
                unconstrained=unconstrained,
            )
        else:
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
                    fine_scores, all_rules, kb=kb,
                    all_marker_dfs=marker_df,
                    ai_results=ai_results,
                    low_quality_clusters=low_quality_clusters,
                    unconstrained=getattr(CFG.ai, 'unconstrained_annotation', False),
                )
                decision_map = dict(zip(decision_clusters, decisions))
        except Exception as exc:
            logger.warning(
                "AI fallback failed: %s — using pure KB results", exc
            )

    # ── g. Map decisions to adata.obs ─────────────────────────────────────
    leiden_str = adata.obs['leiden'].astype(str)

    # For low-quality clusters, downgrade: force Unknown + annotate reason.
    _forced_unknown = 0
    for cl_str, reason in low_quality_clusters.items():
        if cl_str in decision_map and decision_map[cl_str].confidence != 'unknown':
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

    # For low-quality clusters, set cell_category to empty string
    for cl_str in low_quality_clusters:
        if cl_str in cell_category_map:
            cell_category_map[cl_str] = ""


    adata.obs['cell_type'] = leiden_str.map(
        {k: v.cell_type for k, v in decision_map.items()}
    ).astype('category')
    adata.obs['cell_state'] = leiden_str.map(
        {k: 'N/A' for k in decision_map}
    )
    adata.obs['cell_subtype'] = leiden_str.map(
        {k: 'N/A' for k in decision_map}
    )
    adata.obs['cell_category'] = leiden_str.map(
        {k: cell_category_map.get(k, "") for k in decision_map}
    ).astype('category')


    # annot_method: clean label from fusion method (+ AI suffix)
    def _clean_method(d):
        if d.method == 'expert_rule':
            return 'expert_rule'
        if d.method == 'unknown':
            return 'unknown'
        if d.method == 'transition_state':
            return 'transition_state'
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
            'cell_category': cell_category_map.get(str(cl_name), ""),
        })
    ann_df = pd.DataFrame(ann_records)
    ann_csv = os.path.join(CFG.table_dir, 'cell_type_annotations.csv')
    ann_df.to_csv(ann_csv, index=False)
    logger.info("Annotation table saved: %s", ann_csv)

    logger.info("Cluster → cell type mapping (Unified):")
    for rec in ann_records:
        category_str = f" [{rec.get('cell_category', '')}]" if rec.get('cell_category') else ""
        logger.info(
            "  Cluster %s → %s%s (conf=%s, method=%s)",
            rec['cluster'], rec['cell_type'], category_str,
            rec['confidence'], rec['method'],
        )
    for rec in ann_records:
        logger.info(
            "  Cluster %s → %s (conf=%s, method=%s)",
            rec['cluster'], rec['cell_type'],
            rec['confidence'], rec['method'],
        )

    # ── i. UMAP visualization ─────────────────────────────────────────────
    sc.settings.figdir = os.path.join(CFG.figure_dir, '05_annotation')
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False

    adata.obs['annot_label'] = adata.obs['cell_type'].astype(str)

    safe_plot(sc.pl.umap, adata, color='cell_type', show=False,
              save='_05_celltype_unified.pdf')
    safe_plot(sc.pl.umap, adata, color='annot_label', show=False,
              save='_05_annot_label_unified.pdf')
    safe_plot(sc.pl.umap, adata, color='annot_confidence', show=False,
              save='_05_confidence_unified.pdf')
    safe_plot(sc.pl.umap, adata, color='cell_category', show=False,
              save='_05_cell_category_unified.png')

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
        'cell_category': adata.obs['cell_category'].values,
    })
    meta_csv = os.path.join(CFG.table_dir, 'cell_metadata.csv')
    meta_df.to_csv(meta_csv, index=False)
    logger.info("Cell metadata exported: %s", meta_csv)

    # ── k. Annotation quality report ────────────────────────────────────────
    _write_quality_report(adata, ann_records, fusion_quality, cell_category_map, decision_map, CFG, logger)

    # ── l. Interactive review (--interactive flag) ──────────────────────────
    if getattr(CFG, 'interactive', False):
        _interactive_annotation_review(adata, fusion_quality, CFG, logger)

    return decision_map


def _write_quality_report(adata, ann_records, fusion_quality, cell_category_map, decision_map, CFG, logger):
    """Write 05_annotation_quality.json summarising annotation health."""
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
        "categories_found": len(set(c for c in cell_category_map.values() if c)),
        "transition_clusters": sum(1 for d in decision_map.values() if d.method == "transition_state"),
        "category_distribution": {cat: count for cat, count in Counter(cell_category_map.values()).items() if cat},
    }

    quality_path = os.path.join(CFG.table_dir, '05_annotation_quality.json')
    with open(quality_path, 'w', encoding='utf-8') as f:
        json.dump(quality, f, indent=2)
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
        print(f"  💡 Recommended:       strictness='{rec}'")
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
                  f"  CFG.marker.expert_rule_strictness = '{rec}'\n"
                  f"Or pass --config with the updated setting.\n")
        elif choice == 's':
            logger.info("User chose: score_genes fallback")
            print("Set CFG.tissue_kb = '' to use score_genes mode.\n")
        elif choice == 'a':
            raise SystemExit("Aborted by user.")
        else:
            logger.info("User chose: continue with current labels")
