#!/usr/bin/env python3
"""
Step 05: Cell type annotation for spatial transcriptomics
=============================================================
  Reuses the RNA pipeline's annotation engine (AI + score_genes + KB modes)
  with spatial context. Three annotation modes:

    1. AI mode: LLM-based annotation (if CFG.ai.enabled + CFG.ai.ai_annotation)
    2. Score_genes mode: marker gene scoring fallback
    3. If CFG.tissue_kb is set: tissue Knowledge Base mode (full RNA engine)

Input:  04_clustered.h5ad
Output: 05_annotated.h5ad
"""

import argparse
import json
import os
import sys
import time

# Add repo root so `from core.*` and `from rna.*` resolve correctly
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import pandas as pd
import scanpy as sc

from core.utils import resolve_config, safe_plot, safe_write, setup_logger


def score_genes_mode(adata, cfg, log):
    """Marker gene scoring-based annotation fallback."""
    log.info("Score_genes mode — marker gene-based annotation")

    marker_dict = cfg.marker.marker_dict
    if not marker_dict:
        log.warning("marker_dict not configured, using leiden labels as cell_type")
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
        return

    cell_types = list(marker_dict.keys())
    log.info("Scoring %d candidate cell types", len(cell_types))

    for ct in cell_types:
        genes = marker_dict[ct]
        var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
        if adata.raw is None:
            log.warning("adata.raw is None; falling back to adata.var_names")
        genes_present = [g for g in genes if g in var_names]
        if not genes_present:
            log.warning("  %s: no marker genes found in data", ct)
            adata.obs[f"score_{ct}"] = 0.0
            continue
        sc.tl.score_genes(
            adata,
            gene_list=genes_present,
            score_name=f"score_{ct}",
            random_state=cfg.execution.random_seed,
        )

    # Assign best-scoring type per cluster
    score_cols = [f"score_{ct}" for ct in cell_types]
    groupby_kw = {"observed": True}
    cluster_scores = adata.obs.groupby("leiden", **groupby_kw)[score_cols].mean()
    best_match = cluster_scores.idxmax(axis=1)
    best_ct = best_match.str.replace("score_", "")

    cluster_to_ct = dict(zip(best_ct.index, best_ct.values))
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_ct).astype("category")

    log.info("Cluster -> cell_type mapping:")
    for label in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        ct = cluster_to_ct[label]
        max_score = cluster_scores.loc[label, f"score_{ct}"]
        log.info("  Cluster %s -> %s (score=%.3f)", label, ct, max_score)

    # Confidence: difference between top 2 scores
    if len(cell_types) >= 2:
        sorted_scores = cluster_scores.apply(
            lambda row: row.sort_values(ascending=False).values, axis=1, result_type="expand"
        )
        confidence = sorted_scores.iloc[:, 0] - sorted_scores.iloc[:, 1]
        adata.obs["annot_confidence"] = adata.obs["leiden"].map(confidence).astype(float).values


def ai_annotate(adata, cfg, log):
    """AI-based annotation via LLM."""
    log.info("AI mode — computing marker genes...")

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    n_clusters = adata.obs["leiden"].nunique()
    suggested_max_tokens = min(max(4096, n_clusters * 300), 32768)
    _saved_max_tokens = getattr(cfg.ai, "max_tokens", 4096)
    try:
        cfg.ai.max_tokens = max(_saved_max_tokens, suggested_max_tokens)
        log.info(
            "Adaptive max_tokens: n_clusters=%d -> max_tokens=%d", n_clusters, suggested_max_tokens
        )

        # Build marker CSV
        marker_rows = []
        for cl in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
            df = sc.get.rank_genes_groups_df(adata, group=str(cl))
            df["cluster"] = cl
            marker_rows.append(df)
        marker_df = pd.concat(marker_rows, ignore_index=True)
        marker_csv = os.path.join(cfg.table_dir, "marker_genes_ai.csv")
        marker_df.to_csv(marker_csv, index=False)
        log.info("Marker genes saved: %s", marker_csv)

        # Build prompt
        from core.ai.caller import ai_query
        from core.ai.prompts import build_annotation_prompt

        tissue = cfg.tissue
        species = cfg.species
        compact = n_clusters > 20

        # Build spatial context for AI prompt
        spatial_parts = [f"Spatial transcriptomics dataset — {adata.n_obs} spots"]
        platform = cfg.spatial.platform or ""
        if platform:
            spatial_parts.append(f"Platform: {platform}")
        if "spatial" in adata.obsm:
            coords = adata.obsm["spatial"]
            spatial_parts.append(
                f"Tissue coordinates: X [{coords[:, 0].min():.0f}, {coords[:, 0].max():.0f}], "
                f"Y [{coords[:, 1].min():.0f}, {coords[:, 1].max():.0f}]"
            )
        extra_context = " | ".join(spatial_parts)

        sys_prompt, user_prompt = build_annotation_prompt(
            adata,
            tissue,
            species,
            precomputed_rank=True,
            extra_context=extra_context,
            compact=compact,
        )

        # Call LLM
        log.info("Calling LLM (model=%s)...", cfg.ai.model)
        try:
            response = ai_query(sys_prompt, user_prompt, cfg=cfg.ai)
        except Exception as exc:
            log.warning("LLM query failed: %s", exc)
            return None

        # Parse JSON
        try:
            annotations = json.loads(response)
        except (json.JSONDecodeError, TypeError) as e:
            log.warning("LLM response is not valid JSON (%s)", e)
            return None

        # Validate
        required_keys = {"cell_type", "state", "subtype", "confidence", "reasoning"}
        for cid, ann in annotations.items():
            if not isinstance(ann, dict):
                return None
            if required_keys - ann.keys():
                log.warning("Cluster %s missing fields", cid)
                return None

        log.info("LLM annotation: %d clusters parsed", len(annotations))

        # Map to adata.obs
        leiden_str = adata.obs["leiden"].astype(str)
        adata.obs["cell_type"] = leiden_str.map(
            {k: v["cell_type"] for k, v in annotations.items()}
        ).astype("category")
        adata.obs["cell_state"] = leiden_str.map({k: v["state"] for k, v in annotations.items()})
        adata.obs["cell_subtype"] = leiden_str.map(
            {k: v["subtype"] for k, v in annotations.items()}
        )
        adata.obs["annot_confidence"] = leiden_str.map(
            {k: v["confidence"] for k, v in annotations.items()}
        )
        adata.obs["annot_reasoning"] = leiden_str.map(
            {k: v["reasoning"] for k, v in annotations.items()}
        )

        # Save annotation CSV
        ann_records = []
        for cid in sorted(annotations.keys(), key=lambda x: int(x)):
            ann_records.append(
                {
                    "cluster": cid,
                    "cell_type": annotations[cid]["cell_type"],
                    "subtype": annotations[cid]["subtype"],
                    "state": annotations[cid]["state"],
                    "confidence": annotations[cid]["confidence"],
                    "reasoning": annotations[cid]["reasoning"],
                }
            )
        ann_df = pd.DataFrame(ann_records)
        ann_csv = os.path.join(cfg.table_dir, "cell_type_annotations.csv")
        ann_df.to_csv(ann_csv, index=False)
        log.info("Annotation table saved: %s", ann_csv)

        return annotations
    finally:
        cfg.ai.max_tokens = _saved_max_tokens


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("05_annotate", os.path.join(cfg.log_dir, "05_annotate.log"))
    log.info("Step 05: Cell type annotation")

    output_path = os.path.join(cfg.h5ad_dir, "05_annotated.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    adata = sc.read(os.path.join(cfg.h5ad_dir, "04_clustered.h5ad"))
    log.info("Loaded: %d spots, %d clusters", adata.n_obs, adata.obs["leiden"].nunique())

    # ── Phase 1: scRNA marker-list transfer ──────────────────────────────
    # If rna_ref is configured, load scRNA-derived per-cell-type markers
    # and merge them into CFG.marker.marker_dict.  This enriches the score_genes
    # fallback without changing KB or AI mode behaviour.  User-configured
    # marker_dict entries take priority over auto-derived ones.
    _saved_marker_dict = cfg.marker.marker_dict
    try:
        if getattr(cfg, "rna_ref", ""):
            from core.utils import find_rna_marker_csv, load_scRNA_markers

            csv_path = find_rna_marker_csv(cfg=cfg, log=log)
            if csv_path and os.path.exists(csv_path):
                log.info("scRNA marker transfer: loading from %s", csv_path)
                try:
                    scrna_markers = load_scRNA_markers(
                        csv_path,
                        top_n=getattr(cfg, "rna_marker_top_n", 10),
                        pval_threshold=getattr(cfg, "rna_marker_pval_threshold", 0.05),
                        logfc_min=getattr(cfg, "rna_marker_logfc_min", 0.0),
                        log=log,
                    )
                    # Merge: scRNA markers as base, user markers override
                    merged = dict(scrna_markers)
                    merged.update(cfg.marker.marker_dict)
                    log.info(
                        "scRNA marker transfer: loaded %d cell types, "
                        "merged with %d user-configured types -> %d total",
                        len(scrna_markers),
                        len(cfg.marker.marker_dict),
                        len(merged),
                    )
                    cfg.marker.marker_dict = merged
                except Exception as e:
                    log.warning("scRNA marker transfer failed: %s - continuing", e)
            else:
                log.info(
                    "scRNA marker transfer: no marker CSV found for rna_ref='%s'", cfg.rna_ref
                )

        # ── Three annotation modes ──────────────────────────────────────────
        ai_enabled = getattr(cfg.ai, "enabled", False)
        ai_annot_on = getattr(cfg.ai, "ai_annotation", False)

        annot_result = None

        # Mode 1: Unified KB mode (if tissue_kb is set)
        if cfg.tissue_kb:
            log.info("Unified KB mode - tissue_kb='%s'", cfg.tissue_kb)
            try:
                from core.annotation.engine import run_unified_annotation as run_unified

                annot_result = run_unified(adata, cfg, log)
            except Exception as e:
                log.warning("Unified KB annotation failed: %s", e)

        # Mode 2: AI mode
        elif ai_enabled and ai_annot_on:
            log.info("AI mode - LLM-based annotation")
            annot_result = ai_annotate(adata, cfg, log)

        # Mode 3: Score_genes fallback
        if annot_result is None:
            log.info("Falling back to score_genes mode")
            score_genes_mode(adata, cfg, log)
    finally:
        cfg.marker.marker_dict = _saved_marker_dict

    # ── Spatial-aware UMAP visualization ──
    sc.settings.figdir = os.path.join(cfg.figure_dir, "05_annotation")
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False

    # Also compute spatial UMAPs
    if "cell_type" in adata.obs:
        try:
            safe_plot(sc.pl.umap, adata, color="cell_type", show=False, save="_05_celltype.png")
        except Exception as e:
            log.warning("UMAP cell_type plot failed: %s", e)

        # Spatial gene expression of top markers
        try:
            top_genes = []
            for ct in adata.obs["cell_type"].unique():
                if ct in cfg.marker.marker_dict:
                    top_genes.extend(cfg.marker.marker_dict[ct][:2])
            top_genes = list(dict.fromkeys(top_genes))[:8]
            for gene in top_genes:
                if gene in adata.var_names:
                    safe_plot(
                        sc.pl.umap,
                        adata,
                        color=gene,
                        show=False,
                        save=f"_05_marker_{gene}.png",
                        use_raw=True,
                    )
        except Exception as e:
            log.warning("Marker gene UMAP failed: %s", e)

    # ── Save cell metadata ──
    meta_cols = {
        "barcode": adata.obs_names,
        "UMAP_1": adata.obsm["X_umap"][:, 0],
        "UMAP_2": adata.obsm["X_umap"][:, 1],
    }
    if "cell_type" in adata.obs:
        meta_cols.update(
            {
                "cell_type": adata.obs["cell_type"].values,
                "leiden": adata.obs["leiden"].values,
            }
        )
    if "cell_state" in adata.obs:
        meta_cols["cell_state"] = adata.obs["cell_state"].values
    if "cell_subtype" in adata.obs:
        meta_cols["cell_subtype"] = adata.obs["cell_subtype"].values
    if "annot_confidence" in adata.obs:
        meta_cols["annot_confidence"] = adata.obs["annot_confidence"].values

    meta_df = pd.DataFrame(meta_cols)
    meta_csv = os.path.join(cfg.table_dir, "cell_metadata.csv")
    meta_df.to_csv(meta_csv, index=False)
    log.info("Cell metadata exported: %s", meta_csv)

    safe_write(adata, output_path, cfg=cfg)
    log.info("Step 05 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
