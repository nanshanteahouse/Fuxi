#!/usr/bin/env python3
"""
Step 08: GO/KEGG enrichment analysis
========================================
  Reads marker_genes_per_group.csv from Step 06.
  Runs ORA (Enrichr) and/or Pre-ranked GSEA via gseapy.

  Reuses core enrichment logic from the RNA pipeline.

Input:  marker_genes_per_group.csv (Step 06 output)
Output: enrichment/ directory with CSVs + bubble plots
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

    if not CFG.run_enrichment:
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

    log.info("Gene set libraries: %s", CFG.enrichment_gene_sets)
    log.info("Method: %s", CFG.enrichment_method)

    # ── Quality awareness (check annotation quality) ───────────────────
    quality_path = os.path.join(CFG.table_dir, '05_annotation_quality.json')
    if os.path.exists(quality_path):
        import json
        with open(quality_path, 'r') as f:
            q = json.load(f)
        pass_rate = q.get('pass_rate', 0)
        if pass_rate < getattr(CFG, 'marker_validation_pass_rate_min', 0.1):
            log.warning(
                "Annotation PASS rate %.1f%% < %.0f%% — enrichment results may be unreliable",
                pass_rate * 100, CFG.marker_validation_pass_rate_min * 100,
            )

    # ── Run enrichment for each gene set ──
    ora_results = {}
    prerank_results = {}

    for gs in CFG.enrichment_gene_sets:
        gs_name = gs.replace(' ', '_').replace('/', '_')

        if CFG.enrichment_method in ('ora', 'both'):
            log.info("[ORA] Gene set: %s", gs)
            try:
                from rna.steps._enrichment import run_ora
                ora_df = run_ora(marker_df, gs, CFG, log)
                if ora_df is not None and not ora_df.empty:
                    ora_results[gs_name] = ora_df
            except ImportError:
                log.warning("ORA enrichment not available — rna/steps/09_enrichment.py may need to be runnable first")

        if CFG.enrichment_method in ('prerank', 'both'):
            log.info("[GSEA] Gene set: %s", gs)
            try:
                from rna.steps._enrichment import run_prerank
                prerank_df = run_prerank(marker_df, gs, CFG, log)
                if prerank_df is not None and not prerank_df.empty:
                    prerank_results[gs_name] = prerank_df
            except ImportError:
                log.warning("GSEA enrichment not available")

    # ── Save results ──
    table_dir = os.path.join(CFG.table_dir, "enrichment")
    fig_dir = os.path.join(CFG.figure_dir, "enrichment")
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

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

    log.info("Step 08 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
