"""MCP tools for pipeline discovery and status (read-only).

Exposes: list_steps, pipeline_status
Does NOT import the heavy runner module — uses lightweight step definitions.
"""

from __future__ import annotations

import json
import os
from typing import Any

from mcp.server import MCPServer

# Lightweight step registry — mirrors runner.py MODALITY_MAP without heavy imports.
# Keep in sync with core/pipeline/runner.py when steps change.
_STEPS: dict[str, list[tuple[str, str, str]]] = {
    "rna": [
        ("00", "00_load.py", "Load raw data → 00_raw.h5ad"),
        ("01", "01_doublet.py", "Scrublet doublet detection (per sample) → 01_doublet.h5ad"),
        ("02", "02_qc.py", "QC filtering (doublets removed) → 02_qc.h5ad"),
        ("03", "03_integrate.py", "Normalize + HVG + PCA + Harmony → 03_integrated.h5ad"),
        ("04", "04_cluster_umap.py", "Multi-param UMAP + multi-resolution Leiden"),
        ("05", "05_annotate_major.py", "AI-assisted major cell type annotation (dual mode)"),
        ("06", "06_subcluster.py", "Interactive subtype analysis (requires --cell-type)"),
        ("07", "07_markers_de.py", "Differential expression (multi-layer)"),
        ("08", "08_trajectory.py", "PAGA + DPT trajectory analysis"),
        ("09", "09_enrichment.py", "GO/KEGG enrichment + AI interpretation"),
        ("10", "10_exploratory.py", "Exploratory analysis (composition/QC/marker)"),
        ("11", "11_grn.py", "GRN regulatory network analysis (decoupler) → 11_grn.h5ad"),
        ("12", "12_cell_interaction.py", "CCI cell-cell interaction (LIANA+) → tables + figures"),
    ],
    "atac": [
        ("00", "00_load.py", "Load fragments.tsv.gz → 00_raw.h5ad"),
        ("01", "01_doublet.py", "Scrublet doublet detection → 01_doublet.h5ad"),
        ("02", "02_qc.py", "QC filtering + TSS + peak calling + peak matrix → 02_filtered.h5ad"),
        (
            "03",
            "03_process.py",
            "Feature selection + spectral + Harmony + KNN → 03_processed.h5ad",
        ),
        ("04", "04_cluster.py", "Multi-param Leiden + UMAP → 04_clustered.h5ad"),
        ("05", "05_peaks.py", "Post-clustering peak calling → 05_peaks.h5ad"),
        ("06", "06_annotate.py", "AI-assisted chromatin state annotation → 05_annotated.h5ad"),
        ("07", "07_subcluster.py", "Subcluster analysis (placeholder)"),
        ("08", "08_marker_peaks.py", "Differential peak accessibility → marker_peaks.csv"),
        ("09", "09_motif.py", "Motif enrichment → motif_results.csv"),
        ("10", "10_trajectory.py", "ATAC pseudotime trajectory → 10_trajectory.h5ad"),
        (
            "11",
            "11_enrichment.py",
            "GO/KEGG enrichment on peak-associated genes → enrichment_*.csv",
        ),
        ("12", "12_exploratory.py", "Exploratory analysis (placeholder)"),
        ("13", "13_integrate.py", "RNA+ATAC integration via muon → 13_integrated.h5ad"),
    ],
    "spatial": [
        ("00", "00_load.py", "Load spatial data → 00_raw.h5ad (coords + image)"),
        ("01", "01_qc.py", "QC filtering (spots + tissue detection) → 01_qc.h5ad"),
        ("02", "02_image.py", "Image processing (sq.im.process) → 02_image.h5ad"),
        ("03", "03_normalize.py", "Normalize + HVG + spatial graph → 03_processed.h5ad"),
        ("04", "04_cluster.py", "PCA + UMAP + Leiden clustering → 04_clustered.h5ad"),
        ("05", "05_annotate.py", "Cell type annotation (AI / score_genes) → 05_annotated.h5ad"),
        (
            "06",
            "06_spatial_stats.py",
            "DE + SVG + nhood enrichment + co-occurrence → CSVs + figures",
        ),
        ("07", "07_trajectory.py", "Pseudotime analysis → 07_trajectory.h5ad"),
        ("08", "08_enrichment.py", "GO/KEGG enrichment → enrichment CSVs"),
        ("09", "09_exploratory.py", "Spatial visualization → figures + CSVs"),
        (
            "10",
            "10_cell_interaction.py",
            "CCI spatial cell-cell interaction (LIANA+) → tables + figures",
        ),
        ("11", "subcluster.py", "Conditional subclustering per cell type → 05_sub_{type}.h5ad"),
        ("12", "grn.py", "Conditional GRN analysis via decoupler → TF activity CSV + heatmap"),
    ],
    "bulk": [
        ("00", "00_load.py", "Load count matrix (CSV/TSV/h5ad) → 00_raw.h5ad"),
        ("01", "01_qc.py", "Sample QC (library size, gene detection) → 01_qc.h5ad"),
        ("02", "02_de.py", "DESeq2 normalization + DE → 02_de.h5ad + CSVs + figures"),
        ("03", "03_enrichment.py", "GO/KEGG enrichment (GSEApy) → tables/"),
        ("04", "04_exploratory.py", "PCA, heatmaps, volcano plots → figures/"),
        ("05", "05_batch.py", "Batch correction (optional, pycombat) → 05_batch_corrected.h5ad"),
    ],
}

_CHECKPOINTS: dict[str, list[str]] = {
    "rna": [
        "00_raw.h5ad",
        "01_doublet.h5ad",
        "02_qc.h5ad",
        "03_integrated.h5ad",
        "04_clustered.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "04_clustered.h5ad",
        "marker_genes_per_group.csv",
        "05_annotated.h5ad",
        "11_grn.h5ad",
        "05_annotated.h5ad",
    ],
    "atac": [
        "00_raw.h5ad",
        "01_doublet.h5ad",
        "02_filtered.h5ad",
        "03_processed.h5ad",
        "04_clustered.h5ad",
        "05_peaks.h5ad",
        "05_annotated.h5ad",
        "",
        "marker_peaks.csv",
        "motif_results.csv",
        "10_trajectory.h5ad",
        "enrichment_*.csv",
        "",
        "13_integrated.h5ad",
    ],
    "spatial": [
        "00_raw.h5ad",
        "01_qc.h5ad",
        "02_image.h5ad",
        "03_processed.h5ad",
        "04_clustered.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
        "05_annotated.h5ad",
    ],
    "bulk": [
        "00_raw.h5ad",
        "01_qc.h5ad",
        "02_de.h5ad",
        "",
        "",
        "05_batch_corrected.h5ad",
    ],
}
# Sentinel 完成度标记 — 原地写回 / 无产物步骤（mirrors runner.py 的
# *_{modality}_SENTINEL_FILES，Keep in sync）。这些步骤的锚定 checkpoint 在
# 更早步骤就已存在，完成度只看 sentinel 文件。
_SENTINELS: dict[str, dict[int, str]] = {
    "rna": {6: "05_annotated.h5ad.step06_done", 8: "05_final.h5ad.step08_done"},
}


def _check_checkpoint(h5ad_dir: str, ckpt: str) -> bool:
    """Check if a checkpoint file exists and is non-empty."""
    if not ckpt:
        return False
    path = os.path.join(h5ad_dir, ckpt)
    if "*" in path:
        import glob as glob_mod

        return bool(glob_mod.glob(path))
    return os.path.exists(path) and os.path.getsize(path) > 0


def _step_completed(h5ad_dir: str, ckpt: str, sentinel: str = "") -> bool | None:
    """单步完成判定。sentinel 步骤只看 sentinel（锚定文件不可作为完成标志）；
    否则按 checkpoint 文件。无任何 marker 时返回 None（= no checkpoint）。"""
    if sentinel:
        return _check_checkpoint(h5ad_dir, sentinel)
    if not ckpt:
        return None
    return _check_checkpoint(h5ad_dir, ckpt)


def register_pipeline_tools(server: MCPServer) -> None:
    """Register all pipeline-related MCP tools on the given server."""

    @server.tool()
    async def list_steps(modality: str = "rna") -> str:
        """List all available pipeline steps for a given modality.

        Use this to discover what steps are available before running them.

        Args:
            modality: One of "rna", "atac", "spatial", "bulk"

        Returns:
            JSON with step number, script name, and description for each step.
        """
        modality = modality.lower()
        if modality not in _STEPS:
            return json.dumps(
                {
                    "error": f"Unknown modality '{modality}'",
                    "valid": list(_STEPS.keys()),
                },
                indent=2,
            )

        steps_out = []
        for num, script, desc in _STEPS[modality]:
            steps_out.append(
                {
                    "step": int(num),
                    "script": script,
                    "description": desc,
                }
            )

        return json.dumps(
            {
                "modality": modality,
                "n_steps": len(steps_out),
                "steps": steps_out,
            },
            indent=2,
        )

    @server.tool()
    async def pipeline_status(config_path: str = "", modality: str = "rna", gse: str = "") -> str:
        """Check pipeline progress for a dataset — which steps are complete, which is next.

        Provide either config_path (e.g., "projects/rna/GSE12345/config_GSE12345.yaml")
        OR modality + gse to auto-discover the config and h5ad directory.

        Args:
            config_path: Path to the YAML config file
            modality: Modality (used if config_path not provided)
            gse: GSE accession ID (used if config_path not provided)

        Returns:
            JSON with per-step status (completed/pending), next step number, and config info.
        """
        h5ad_dir = ""
        actual_modality = modality.lower()

        if config_path:
            # Try to load config to get h5ad_dir
            try:
                from core.utils._config import resolve_config

                cfg = resolve_config(os.path.abspath(config_path))
                h5ad_dir = cfg.h5ad_dir or ""
                actual_modality = cfg.modality or actual_modality
            except Exception:
                # Config loading failed — try to infer h5ad_dir from path
                config_dir = os.path.dirname(os.path.abspath(config_path))
                h5ad_dir = config_dir
        elif gse:
            gse = gse.upper()
            # Auto-discover config directory
            config_dir = os.path.join("projects", actual_modality, gse)
            if os.path.isdir(config_dir):
                h5ad_dir = config_dir
                # Try to find exact h5ad output dir
                for sub in ["h5ad", "checkpoints", "output"]:
                    candidate = os.path.join(config_dir, sub)
                    if os.path.isdir(candidate):
                        h5ad_dir = candidate
                        break
        else:
            return json.dumps(
                {"error": "Must provide config_path OR (modality + gse)"},
                indent=2,
            )

        if actual_modality not in _STEPS:
            return json.dumps(
                {
                    "error": f"Unknown modality '{actual_modality}'",
                    "valid": list(_STEPS.keys()),
                },
                indent=2,
            )

        if not h5ad_dir:
            return json.dumps(
                {
                    "error": "Could not determine h5ad directory",
                    "suggestion": "Provide a config_path or ensure the dataset exists under projects/{modality}/{GSE}/",
                },
                indent=2,
            )

        steps = _STEPS[actual_modality]
        checkpoints = _CHECKPOINTS[actual_modality]

        result: dict[str, Any] = {
            "modality": actual_modality,
            "h5ad_dir": h5ad_dir,
        }

        if config_path:
            result["config_path"] = config_path

        if gse:
            result["gse"] = gse

        step_statuses = []
        first_incomplete = None

        mod_sentinels = _SENTINELS.get(actual_modality, {})
        for i, (num, _script, desc) in enumerate(steps):
            ckpt = checkpoints[i] if i < len(checkpoints) else ""
            sentinel = mod_sentinels.get(i)
            marker = sentinel or ckpt
            completed = _step_completed(h5ad_dir, ckpt, sentinel)
            step_statuses.append(
                {
                    "step": int(num),
                    "description": desc,
                    "checkpoint": marker or None,
                    "status": "completed"
                    if completed
                    else ("pending" if completed is False else "no_checkpoint"),
                }
            )
            if completed is False and first_incomplete is None:
                first_incomplete = int(num)

        result["step_statuses"] = step_statuses
        result["n_completed"] = sum(1 for s in step_statuses if s["status"] == "completed")
        result["n_total"] = len(step_statuses)

        if first_incomplete is not None:
            result["next_step"] = first_incomplete
            result["next_action"] = (
                f"python core/run_pipeline.py --modality {actual_modality} "
                f"--step {first_incomplete}" + (f" --config {config_path}" if config_path else "")
            )
        elif all(s["status"] in ("completed", "no_checkpoint") for s in step_statuses):
            result["next_step"] = None
            result["status"] = "all_complete"
        else:
            result["next_step"] = 0

        return json.dumps(result, indent=2, ensure_ascii=False)
