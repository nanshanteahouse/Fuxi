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
    from rna.tissue_ontologies import load_kb
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
    rule_top_n, rule_pval = resolve_expert_rule_params(
        strictness=getattr(CFG.marker, 'expert_rule_strictness', 'default'),
        top_n=getattr(CFG.marker, 'expert_rule_top_n', 0),
        pval_cutoff=getattr(CFG.marker, 'expert_rule_pval_cutoff', 0.0),
    )
    logger.info(
        "Expert rules: strictness=%s → top_n=%d, pval_cutoff=%.3f",
        getattr(CFG.marker, 'expert_rule_strictness', 'default'),
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

    logger.info("Cluster → cell type mapping (Unified):")
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
