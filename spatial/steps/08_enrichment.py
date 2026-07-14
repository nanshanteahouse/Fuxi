#!/usr/bin/env python3
"""
Step 08: GO/KEGG enrichment analysis
========================================
  Reads marker_genes_per_group.csv from Step 06.
  Runs ORA (Enrichr) and/or Pre-ranked GSEA via gseapy.

  Reuses core enrichment logic from the RNA pipeline.

Input:  marker_genes_per_group.csv (Step 06 output)
Output: 08_enrichment/ directory with CSVs + bubble plots
"""
import sys, os, time, argparse, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import pandas as pd
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("08_enrichment", os.path.join(CFG.log_dir, "08_enrichment.log"))
    log.info("Step 08: GO/KEGG enrichment analysis")

    if not CFG.enrichment.run:
        log.info("Enrichment analysis disabled (run_enrichment=False)")
        return

    # ── Read marker CSV ──────────────────────────────────────────────────
    marker_path = os.path.join(CFG.table_dir, "marker_genes_per_group.csv")
    if not os.path.exists(marker_path):
        log.error("Marker gene file not found: %s", marker_path)
        log.error("Run Step 06 (06_spatial_de.py) first.")
        sys.exit(1)

    marker_df = pd.read_csv(marker_path)
    log.info("Loaded marker genes: %d rows, %d groups",
             len(marker_df), marker_df['group'].nunique())

    log.info("Gene set libraries: %s", CFG.enrichment.gene_sets)
    log.info("Method: %s", CFG.enrichment.method)

    # ── Quality awareness (check annotation quality) ───────────────────
    quality_path = os.path.join(CFG.table_dir, '05_annotation_quality.json')
    if os.path.exists(quality_path):
        import json
        with open(quality_path, 'r') as f:
            q = json.load(f)
        pass_rate = q.get('pass_rate', 0)
        if pass_rate < getattr(CFG.marker, 'validation_pass_rate_min', 0.1):
            log.warning(
                "Annotation PASS rate %.1f%% < %.0f%% — enrichment results may be unreliable",
                pass_rate * 100, CFG.marker.validation_pass_rate_min * 100,
            )

    # ── Run enrichment for each gene set ──
    ora_results = {}
    prerank_results = {}

    for gs in CFG.enrichment.gene_sets:
        gs_name = gs.replace(' ', '_').replace('/', '_')

        if CFG.enrichment.method in ('ora', 'both'):
            log.info("[ORA] Gene set: %s", gs)
            try:
                from rna.steps import run_ora
                ora_df = run_ora(marker_df, gs, CFG, log)
                if ora_df is not None and not ora_df.empty:
                    ora_results[gs_name] = ora_df
            except ImportError:
                log.warning("ORA enrichment not available — rna/steps/09_enrichment.py may need to be runnable first")

        if CFG.enrichment.method in ('prerank', 'both'):
            log.info("[GSEA] Gene set: %s", gs)
            try:
                from rna.steps import run_prerank
                prerank_df = run_prerank(marker_df, gs, CFG, log)
                if prerank_df is not None and not prerank_df.empty:
                    prerank_results[gs_name] = prerank_df
            except ImportError:
                log.warning("GSEA enrichment not available")

    # ── Save results ──
    table_dir = os.path.join(CFG.table_dir, "08_enrichment")
    fig_dir = os.path.join(CFG.figure_dir, "08_enrichment")
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # ── Tissue-aware post-processing (v4.0+) ──
    tissue_mode = getattr(CFG.enrichment, 'tissue_mode', 'off')
    do_redundancy = getattr(CFG.enrichment, 'redundancy_cluster', False)
    do_kb = getattr(CFG.enrichment, 'use_kb_relevance', False)

    if tissue_mode != 'off' or do_redundancy or do_kb:
        from core.enrichment_tissue import (
            compute_pathway_relevance,
            cluster_redundant_pathways,
            filter_enrichment_by_tissue,
        )

        # 加载 KB markers
        kb_markers = None
        if do_kb and CFG.tissue:
            from rna.tissue_ontologies import load_kb
            kb = load_kb(CFG.tissue)
            kb_markers = set()
            for ct, entry in kb.items():
                if isinstance(entry, dict) and 'markers' in entry:
                    for tier in ('confirm', 'add'):
                        kb_markers.update(entry['markers'].get(tier, {}).keys())

        # 加载通路元数据
        pathway_whitelist = list(getattr(CFG.enrichment, 'tissue_pathways_whitelist', []))
        pathway_blacklist = list(getattr(CFG.enrichment, 'tissue_pathways_blacklist', []))
        if do_kb and CFG.tissue:
            from rna.tissue_ontologies import load_pathway_relevance
            pr = load_pathway_relevance(CFG.tissue)
            if pr:
                if not getattr(CFG.enrichment, 'tissue_pathways_whitelist', []):
                    pathway_whitelist = pr.get('key_pathways', [])
                if not getattr(CFG.enrichment, 'tissue_pathways_blacklist', []):
                    pathway_blacklist = pr.get('generic_pathways', [])

        for results_dict in [ora_results, prerank_results]:
            for gs_name, df in results_dict.items():
                if df.empty:
                    continue
                if do_kb and kb_markers:
                    df = compute_pathway_relevance(df, kb_markers, log)
                if do_redundancy and 'Overlap' in df.columns:
                    df = cluster_redundant_pathways(
                        df, 'Term', 'Overlap',
                        getattr(CFG.enrichment, 'redundancy_threshold', 0.6),
                    )
                if tissue_mode != 'off':
                    df = filter_enrichment_by_tissue(
                        df, tissue_mode,
                        pathway_whitelist, pathway_blacklist,
                        log=log,
                    )
                results_dict[gs_name] = df

    for gs_name, df in ora_results.items():
        path = os.path.join(table_dir, f"ora_{gs_name}_summary.csv")
        df.to_csv(path, index=False)
        log.info("  ORA exported: %s (%d rows)", path, len(df))

    for gs_name, df in prerank_results.items():
        path = os.path.join(table_dir, f"prerank_{gs_name}_summary.csv")
        df.to_csv(path, index=False)
        log.info("  GSEA exported: %s (%d rows)", path, len(df))

    total_ora = sum(len(df) for df in ora_results.values())
    total_gsea = sum(len(df) for df in prerank_results.values())
    log.info("Enrichment results: ORA %d rows, GSEA %d rows", total_ora, total_gsea)

    # ── Tissue-aware filtered copies (v4.0+) ──
    if tissue_mode != 'off':
        for gs_name, df in ora_results.items():
            if df.empty or 'tissue_relevant' not in df.columns:
                continue
            relevant = df[df['tissue_relevant'] == True]
            if not relevant.empty:
                path = os.path.join(table_dir, f"ora_{gs_name}_tissue_relevant.csv")
                relevant.to_csv(path, index=False)
                log.info("  Tissue-relevant ORA: %s (%d/%d rows)",
                         path, len(relevant), len(df))

        for gs_name, df in prerank_results.items():
            if df.empty or 'tissue_relevant' not in df.columns:
                continue
            relevant = df[df['tissue_relevant'] == True]
            if not relevant.empty:
                path = os.path.join(table_dir, f"prerank_{gs_name}_tissue_relevant.csv")
                relevant.to_csv(path, index=False)
                log.info("  Tissue-relevant GSEA: %s (%d/%d rows)",
                         path, len(relevant), len(df))

        total_relevant = 0
        total_all = 0
        for df in ora_results.values():
            if 'tissue_relevant' in df.columns:
                total_all += len(df)
                total_relevant += df['tissue_relevant'].sum()
        for df in prerank_results.values():
            if 'tissue_relevant' in df.columns:
                total_all += len(df)
                total_relevant += df['tissue_relevant'].sum()
        if total_all > 0:
            log.info("Tissue-aware enrichment: mode=%s, %d/%d (%.0f%%) pathways marked relevant",
                     tissue_mode, total_relevant, total_all,
                     total_relevant / total_all * 100)

    log.info("Step 08 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
