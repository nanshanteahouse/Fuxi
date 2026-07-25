#!/usr/bin/env python3
"""
preprocessor.py — Fuxi 预处理管线
====================================

从 GEO 下载目录到可运行管线的"一条龙"自动预处理。

工作流程:
  Phase 0 — 验证输入目录
  Phase 1 — 解压所有归档文件 (tar.gz, zip, gz, bz2)
  Phase 2 — SuperSeries 检测
  Phase 3 — 格式检测 & 模态推断
  Phase 4 — 生成 dataset.yaml
  Phase 5 — 生成 config_GSE_ID.yaml
  Phase 6 — 汇总报告

用法:
    python core/preprocess/preprocessor.py --gse GSE12345
    python core/preprocess/preprocessor.py --gse GSE12345 --dry-run
    python core/preprocess/preprocessor.py --gse GSE12345 --query-ncbi --verbose
"""

import argparse
import os
import shutil
import sys
import time
from datetime import datetime
from typing import Optional

# Add repo root to sys.path (consistent with all step scripts)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from core.config.dataset import (
    DatasetMeta,
    Meta,
    ModalityEntry,
    PipelineStatus,
    Resources,
    save_dataset,
)
from core.preprocess import archive_extractor as ae
from core.preprocess import format_detector as fd
from core.preprocess import superseries_detector as ssd
from core.preprocess.matrix_loader import (
    _detect_primary_format,
    _resolve_project_dir,
    generate_config,
)
from core.preprocess.metadata_parser import (
    generate_dataset_yaml,
)

# ═══════════════════════════════════════════════════════════════════════
#  Phase 3: Modality detection
# ═══════════════════════════════════════════════════════════════════════


def _infer_modality(classification: dict) -> str:
    """Infer the primary modality from classification results.

    Returns 'rna', 'atac', 'spatial', or 'multiome'.
    """
    has_rna = bool(
        classification.get("tenx_h5_dirs")
        or classification.get("tenx_mtx_dirs")
        or classification.get("h5ad_files")
        or classification.get("csv_files")
    )
    has_atac = bool(classification.get("fragment_dirs") or classification.get("tenx_peak_dirs"))
    if has_rna and has_atac:
        return "multiome"
    if has_atac:
        return "atac"
    # Default to RNA for h5ad / CSV cases
    return "rna"


# ═══════════════════════════════════════════════════════════════════════
#  Main entry point
# ═══════════════════════════════════════════════════════════════════════


def _group_files_by_accession(
    file_list: list[str], child_accessions: list[str]
) -> dict[str, list[str]]:
    """Group *file_list* into per-accession buckets using filename patterns.

    Delegates to ``superseries_detector.group_files_by_accession()``, which
    matches each file's basename against the given child GSE accessions.
    """
    return ssd.group_files_by_accession(file_list, child_accessions)


def _generate_parent_dataset_yaml(
    gse_id: str,
    superseries_info: dict,
    classification: dict,
    file_list: list[str],
    output_dir: str,
    data_root: Optional[str],
    input_dir_override: Optional[str] = None,
    dry_run: bool = False,
    force: bool = False,
    ncbi_assay_type: Optional[str] = None,
) -> Optional[str]:
    """Generate a parent-level dataset.yaml for a SuperSeries.

    This YAML serves as an index/placeholder — it records the SuperSeries
    metadata and lists the child accessions but does not reference pipeline
    runnables.  The actual per-child dataset.yaml files live under each
    child's project directory (flat layout).

    Returns the path to the generated file, or None.
    """
    species = fd.guess_species(file_list)
    if species == "unknown" and superseries_info:
        ncbi_species = superseries_info.get("species", "")
        if ncbi_species:
            species = fd._normalise_species(ncbi_species)
    species_key = species  # already normalised by guess_species()
    tissue = fd.guess_tissue(file_list)

    # Determine data format from classification
    data_format = _detect_primary_format(classification)

    # Determine base dir for path computation
    if input_dir_override:
        os.path.abspath(input_dir_override)
    elif data_root:
        os.path.join(data_root, gse_id)
    else:
        os.path.dirname(file_list[0]) if file_list else "."

    # Modality entry: placeholder
    modality_entry = ModalityEntry(
        name="SuperSeries",
        status="placeholder",
        format=data_format,
        file_count=len(file_list),
    )

    # Build subseries list
    subseries = []
    for child_acc in superseries_info.get("child_accessions", []):
        subseries.append(
            {
                "id": child_acc,
                "title": "",
                "modality": "scRNA-seq",
            }
        )

    ds = DatasetMeta(
        id=gse_id,
        type="SuperSeries",
        title=superseries_info.get("title", ""),
        species=species if species != "unknown" else "homo_sapiens",
        species_key=species_key if species_key != "unknown" else "human",
        tissue=tissue if tissue != "unknown" else None,
        modalities=[modality_entry],
        samples=[],
        subseries=subseries,
        comparisons=[],
        assay_type=ncbi_assay_type,
        resources=Resources(
            genome=fd.guess_genome(species),
        ),
        meta=Meta(
            created=datetime.now().isoformat(),
            generated_by="fuxi_preprocess",
            pipeline_status=PipelineStatus(),
        ),
    )

    os.makedirs(output_dir, exist_ok=True)
    yaml_path = os.path.join(output_dir, "dataset.yaml")

    if dry_run:
        print(f"    [DRY-RUN] Would write: {yaml_path}")
        return yaml_path

    if os.path.exists(yaml_path) and not force:
        print(f"    [SKIP] {yaml_path} — use --force to overwrite.")
        return yaml_path

    save_dataset(ds, yaml_path)
    print(f"    Written: {yaml_path}")
    return yaml_path


def run_preprocess(
    gse_id: Optional[str] = None,
    input_dir: Optional[str] = None,
    dataset_name: Optional[str] = None,
    data_root: Optional[str] = None,
    paper_context: Optional[dict] = None,
    query_ncbi: bool = False,
    download: bool = False,
    dry_run: bool = False,
    force: bool = False,
    no_extract: bool = False,
    modality: Optional[str] = None,
    output_dir: Optional[str] = None,
    verbose: bool = False,
    quiet: bool = False,
) -> int:
    """Run the full preprocessing pipeline on a dataset directory.

    Two calling conventions are supported:

    1.  GEO mode (--gse):  ``data_root`` / ``gse_id`` is the input dir.
    2.  Custom mode (--input-dir):  *input_dir* IS the input directory.

    When *gse_id* is set, it becomes the dataset id for naming output.
    In custom mode, *dataset_name* (or the basename of *input_dir*)
    supplies the dataset id.

    Args:
        gse_id:        GEO accession ID (GEO mode).
        input_dir:     Direct input directory path (custom mode).
        dataset_name:  Dataset identifier when *input_dir* is used
                       (defaults to basename of *input_dir*).
        data_root:     Root data directory (default: FUXI_DATA_ROOT env var).
        paper_context: Optional dict with paper-derived metadata values
                       (species, tissue, expression_type, genome, assay_type,
                       features, is_nuclei, etc.).
        query_ncbi:    Query NCBI for metadata (GEO mode with internet).
        download:     Auto-download GSE data from GEO before preprocessing
                       (GEO mode only; requires wget or curl).
        dry_run:       Report only, don't write files.
        force:         Overwrite existing files.
        no_extract:    Skip archive extraction.
        modality:      Force a modality override.
        output_dir:    Base output directory override.
        verbose:       Print detailed per-file detection info.
        quiet:         Minimal output.
    Returns:
        Exit code (0 = success, 1 = error).
    """
    # ── Phase 0: Resolve input directory & dataset id ──────────────────
    if input_dir is not None:
        # Custom mode — input_dir IS the data directory
        gse_dir = os.path.abspath(input_dir)
        if dataset_name:
            gse_id = dataset_name
        else:
            gse_id = os.path.basename(gse_dir.rstrip("/\\")) or "unknown_dataset"
    elif gse_id is not None:
        # GEO mode — gse_id under data_root
        if data_root is None:
            try:
                from core.utils import data_root as get_data_root

                data_root = get_data_root()
            except RuntimeError as e:
                print(f"[ERROR] {e}", file=sys.stderr)
                return 1
        gse_dir = os.path.join(data_root, gse_id)

        # ── Optional: auto-download from GEO ───────────────────────
        if download:
            from core.geo_downloader import download_gse

            if not quiet:
                print(f"\n[Phase 0] Auto-downloading {gse_id} from GEO...")
            dl_result = download_gse(
                gse_id=gse_id,
                dest_dir=gse_dir,
                dry_run=dry_run,
                skip_existing=not force,
                quiet=quiet,
            )
            # failed == -1 signals a hard listing/download failure
            if dl_result.get("failed", 0) < 0:
                print(
                    f"[ERROR] Failed to list/download files for {gse_id}.",
                    file=sys.stderr,
                )
                return 1
    else:
        print("[ERROR] Either --gse or --input-dir must be provided.", file=sys.stderr)
        return 1

    if not os.path.isdir(gse_dir):
        print(f"[ERROR] Directory not found: {gse_dir}", file=sys.stderr)
        if gse_id and not input_dir:
            print(
                "        Make sure FUXI_DATA_ROOT is set correctly and the dataset is downloaded.",
                file=sys.stderr,
            )
        return 1

    t_start = time.time()

    if not quiet:
        print(f"\n{'=' * 60}")
        print(f"Fuxi Preprocessing: {gse_id}")
        print(f"Data root: {data_root}")
        print(f"{'=' * 60}\n")

    # ── Phase 1: Archive extraction ────────────────────────────────────
    archive_results = []
    if not no_extract:
        if not quiet:
            print("[Phase 1] Scanning for archives...")
        all_files_before = ae.collect_file_tree(gse_dir)
        archives = fd.classify_files_by_format(all_files_before)["archives"]

        if archives:
            if not quiet:
                print(f"  Found {len(archives)} archive(s)")
            archive_results = ae.extract_all_archives(archives, verbose=verbose)
            errors = [r for r in archive_results if r["status"] == "error"]
            for err in errors:
                print(f"  [WARNING] Extraction failed: {err['error']}", file=sys.stderr)

            # Second pass: check for nested archives in extraction dirs
            second_pass_files = ae.collect_file_tree(gse_dir)
            second_archives = fd.classify_files_by_format(second_pass_files)["archives"]
            new_archives = [(p, f) for p, f in second_archives if (p, f) not in archives]
            if new_archives:
                if not quiet:
                    print(f"  Second pass: {len(new_archives)} nested archive(s) found")
                ae.extract_all_archives(new_archives, verbose=verbose)

            # Third pass (cap at 3)
            third_pass_files = ae.collect_file_tree(gse_dir)
            third_archives = fd.classify_files_by_format(third_pass_files)["archives"]
            newest = [
                (p, f)
                for p, f in third_archives
                if (p, f) not in archives and (p, f) not in new_archives
            ]
            if newest:
                if not quiet:
                    print(f"  Third pass: {len(newest)} deeply nested archive(s)")
                ae.extract_all_archives(newest, verbose=verbose)
        else:
            if not quiet:
                print("  No archives found.")

    # ── Collect all files after extraction ─────────────────────────────
    all_files = ae.collect_file_tree(gse_dir)
    if not all_files:
        print(f"[ERROR] No files found in {gse_dir} after extraction.", file=sys.stderr)
        return 1

    if not quiet:
        print(f"  Total files: {len(all_files)}")

    # ── Phase 2: SuperSeries detection ─────────────────────────────────
    if not quiet:
        print("\n[Phase 2] Checking for SuperSeries structure...")
    superseries_info = ssd.detect_superseries(
        root_dir=gse_dir,
        file_list=all_files,
        gse_id=gse_id,
        query_ncbi_flag=query_ncbi,
    )
    if superseries_info.get("is_superseries"):
        method = superseries_info.get("detected_by", "unknown")
        if not quiet:
            print(f"  SuperSeries detected (via {method}).")
            if superseries_info.get("subseries_dirs"):
                print(f"  Subseries dirs: {', '.join(superseries_info['subseries_dirs'])}")
            if superseries_info.get("child_accessions"):
                print(f"  Child accessions: {', '.join(superseries_info['child_accessions'])}")
    else:
        if not quiet:
            print("  Not a SuperSeries (single accession).")
    # ── Phase 3: Format detection ──────────────────────────────────────
    if not quiet:
        print("\n[Phase 3] Detecting file formats...")
    classification = fd.classify_files_by_format(all_files)

    if verbose:
        print(f"  Archives:   {len(classification['archives'])}")
        print(f"  10X MTX:    {len(classification['tenx_mtx_dirs'])} directories")
        print(f"  10X H5:     {len(classification['tenx_h5_dirs'])} directories")
        print(f"  Fragments:  {len(classification['fragment_dirs'])} directories")
        print(f"  h5ad:       {len(classification['h5ad_files'])} files")
        print(f"  CSV:        {len(classification['csv_files'])} files")
        print(f"  Metadata:   {len(classification['metadata_files'])} files")
        print(f"  Unmatched:  {len(classification['unmatched'])} files")

    # Determine modality
    if modality:
        detected_modality = modality
    else:
        detected_modality = _infer_modality(classification)
    if not quiet:
        print(f"  Inferred modality: {detected_modality}")

    # ── Extract assay_type from NCBI metadata ─────────────────────────
    ncbi_assay_type = superseries_info.get("assay_type")
    # Multiome special case: 10x Multiome is nucleus-based
    if ncbi_assay_type is None and detected_modality == "multiome":
        ncbi_assay_type = "snRNAseq"
        if not quiet:
            print("  [INFO] Multiome detected — auto-set assay_type=snRNAseq")

    # List unsupported files
    unsupported = classification.get("unsupported", [])
    if unsupported:
        print(f"  [WARNING] {len(unsupported)} unsupported format(s):", file=sys.stderr)
        for f in unsupported[:10]:
            print(f"    - {os.path.relpath(f, gse_dir)}", file=sys.stderr)

    # List unmatched files
    unmatched = classification.get("unmatched", [])
    if unmatched:
        if verbose:
            print(f"  Unmatched files ({len(unmatched)}):")
            for f in sorted(unmatched):
                print(f"    - {os.path.relpath(f, gse_dir)}")
        elif not quiet:
            print(f"  {len(unmatched)} unmatched file(s) (use -v for details)")

    # ── Resolve output modalities ────────────────────────────────────────
    if detected_modality == "multiome":
        modalities_out = ["rna", "atac"]
    else:
        modalities_out = [detected_modality]

    # ── Phase 4+5: Generate project files ────────────────────────────────
    if not quiet:
        print()

    if superseries_info.get("is_superseries"):
        # ── Flat layout: move files → sibling dirs, then generate ────────
        if not quiet:
            print("[Phase 4/5] SuperSeries: organising child accessions...")
        child_accs = superseries_info.get("child_accessions", [])

        # Group files by child accession using filename patterns
        file_groups = ssd.group_files_by_accession(all_files, child_accs)

        # Parent dir for sibling child dirs (= same level as gse_dir)
        data_parent = os.path.dirname(gse_dir.rstrip("/\\"))

        for child_gse in child_accs:
            child_files = file_groups.get(child_gse, [])
            if not child_files:
                if not quiet:
                    print(f"  [SKIP] {child_gse}: no matching files found")
                continue

            # Create sibling directory: e.g. /data/neurobiology/GSE123456/
            child_data_dir = os.path.join(data_parent, child_gse)
            os.makedirs(child_data_dir, exist_ok=True)

            if not quiet:
                print(f"  {child_gse}: {len(child_files)} file(s) → {child_data_dir}")

            # Move files from parent to sibling dir
            moved_files = []
            for f in child_files:
                dest = os.path.join(child_data_dir, os.path.basename(f))
                if os.path.abspath(f) == os.path.abspath(dest):
                    # File already in the right place
                    moved_files.append(dest)
                    continue
                if dry_run:
                    print(f"    [DRY-RUN] mv {f} → {dest}")
                    moved_files.append(dest)
                else:
                    if os.path.exists(dest) and not force:
                        print(f"    [SKIP] {dest} exists (use --force to overwrite)")
                        moved_files.append(dest)
                    else:
                        shutil.move(f, dest)
                        moved_files.append(dest)

            # Re-collect after move for unified file list
            if dry_run:
                # Files weren't actually moved — use original paths for detection
                child_abs_files = sorted(child_files)
            else:
                child_abs_files = sorted([f for f in moved_files if os.path.exists(f)])

            # Per-child format classification on the new location
            child_classification = fd.classify_files_by_format(child_abs_files)
            child_modality = modality or _infer_modality(child_classification)
            if child_modality == "multiome":
                child_modalities = ["rna", "atac"]
            else:
                child_modalities = [child_modality]

            for child_mod in child_modalities:
                # dataset.yaml IN the sibling data dir itself
                generate_dataset_yaml(
                    gse_id=child_gse,
                    modality=child_mod,
                    superseries_info={"is_superseries": False},
                    classification=child_classification,
                    file_list=child_abs_files,
                    output_dir=child_data_dir,
                    data_root=data_root,
                    input_dir_override=child_data_dir,
                    paper_context=paper_context,
                    dry_run=dry_run,
                    force=force,
                    ncbi_assay_type=ncbi_assay_type,
                )

                # config.py in projects/
                child_proj_dir = _resolve_project_dir(child_mod, child_gse, output_dir)
                generate_config(
                    gse_id=child_gse,
                    modality=child_mod,
                    classification=child_classification,
                    file_list=child_abs_files,
                    output_dir=child_proj_dir,
                    data_root=data_root,
                    input_dir_override=child_data_dir,
                    superseries_info=superseries_info,
                    paper_context=paper_context,
                    dry_run=dry_run,
                    force=force,
                )

        # Parent-level placeholder dataset.yaml — written to gse_dir itself
        if not quiet:
            print(f"  Parent: {gse_id} → {gse_dir}/dataset.yaml")
        _generate_parent_dataset_yaml(
            gse_id=gse_id,
            superseries_info=superseries_info,
            classification=classification,
            file_list=all_files,
            output_dir=gse_dir,
            data_root=data_root,
            input_dir_override=input_dir or gse_dir,
            dry_run=dry_run,
            force=force,
            ncbi_assay_type=ncbi_assay_type,
        )
    else:
        # ── Single-accession original flow ────────────────────────────────
        if not quiet:
            print("[Phase 4] Generating dataset.yaml...")
        for mod in modalities_out:
            proj_dir = _resolve_project_dir(mod, gse_id, output_dir)
            generate_dataset_yaml(
                gse_id=gse_id,
                modality=mod,
                superseries_info=superseries_info,
                classification=classification,
                file_list=all_files,
                output_dir=proj_dir,
                data_root=data_root,
                input_dir_override=input_dir,
                paper_context=paper_context,
                dry_run=dry_run,
                force=force,
                ncbi_assay_type=ncbi_assay_type,
            )

        if not quiet:
            print("\n[Phase 5] Generating config file...")
        for mod in modalities_out:
            proj_dir = _resolve_project_dir(mod, gse_id, output_dir)
            generate_config(
                gse_id=gse_id,
                modality=mod,
                classification=classification,
                file_list=all_files,
                output_dir=proj_dir,
                data_root=data_root,
                input_dir_override=input_dir,
                superseries_info=superseries_info,
                paper_context=paper_context,
                dry_run=dry_run,
                force=force,
            )

    # ── Phase 6: Summary ───────────────────────────────────────────────
    elapsed = time.time() - t_start
    if not quiet:
        print(f"\n{'=' * 60}")
        print(f"[Summary] {gse_id}")
        print(f"{'=' * 60}")
        print(
            f"  Type:         {'SuperSeries' if superseries_info.get('is_superseries') else 'SingleAccession'}"
        )
        print(f"  Modality:     {detected_modality}")
        print(f"  Data format:  {_detect_primary_format(classification, detected_modality)}")
        print(
            f"  Species:      {fd.guess_species(all_files) if fd.guess_species(all_files) != 'unknown' else fd._normalise_species(superseries_info.get('species', 'unknown'))}"
        )
        print(f"  Tissue:       {fd.guess_tissue(all_files)}")
        extracted = sum(1 for r in archive_results if r["status"] in ("extracted", "skipped"))
        if extracted:
            print(f"  Archives:     {extracted} processed")
        if unsupported:
            print(f"  Warnings:     {len(unsupported)} unsupported format(s)")
        if unmatched and not verbose:
            print(f"  Unmatched:    {len(unmatched)} file(s)")
        print(f"  Elapsed:      {elapsed:.1f}s")
        if not dry_run:
            print("\n  Generated:")
            if superseries_info.get("is_superseries"):
                data_parent = os.path.dirname(gse_dir.rstrip("/\\"))
                for child_gse in superseries_info.get("child_accessions", []):
                    child_dir = os.path.join(data_parent, child_gse)
                    print(
                        f"    {os.path.join(child_dir, 'dataset.yaml')}  (files moved to {child_dir}/)"
                    )
                    for mod in modalities_out:
                        out_dir = _resolve_project_dir(mod, child_gse, output_dir)
                        print(f"    {os.path.join(out_dir, 'config_' + child_gse + '.yaml')}")
                print(f"    {os.path.join(gse_dir, 'dataset.yaml')}  (parent index)")
            else:
                for mod in modalities_out:
                    out_dir = _resolve_project_dir(mod, gse_id, output_dir)
                    print(f"    {os.path.join(out_dir, 'dataset.yaml')}")
                    print(f"    {os.path.join(out_dir, 'config_' + gse_id + '.yaml')}")
        print("\n  Next steps:")
        if superseries_info.get("is_superseries"):
            print("    Each sub-series is an independent project. Review and edit")
            print("    the generated configs, then run each separately, e.g.:")
            for child_gse in superseries_info.get("child_accessions", []):
                for mod in modalities_out:
                    cfg_path = os.path.join(
                        _resolve_project_dir(mod, child_gse, output_dir),
                        f"config_{child_gse}.yaml",
                    )
                    print(
                        f"       python core/run_pipeline.py --modality {mod} --config {cfg_path}"
                    )
        else:
            print("    1. Review and edit the generated files")
            print("    2. Run the pipeline:")
            for mod in modalities_out:
                cfg_path = os.path.join(
                    _resolve_project_dir(mod, gse_id, output_dir), f"config_{gse_id}.yaml"
                )
                print(f"       python core/run_pipeline.py --modality {mod} --config {cfg_path}")
        print()

    return 0


# ═══════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description="Fuxi Pre-Processing Pipeline — Prepare downloaded GEO data "
        "for pipeline execution.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python core/preprocess/preprocessor.py --gse GSE12345
  python core/preprocess/preprocessor.py --gse GSE12345 --dry-run
  python core/preprocess/preprocessor.py --gse GSE12345 --query-ncbi --verbose
  python core/preprocess/preprocessor.py --gse GSE12345 --force
  python core/preprocess/preprocessor.py --gse GSE12345 --no-extract
  python core/preprocess/preprocessor.py --gse GSE12345 --modality atac
  python core/preprocess/preprocessor.py --input-dir /data/my_experiment --name my_label
  python core/preprocess/preprocessor.py --input-dir /data/my_experiment --dry-run
        """,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--gse",
        type=str,
        help="GEO accession ID (e.g., GSE12345). Requires FUXI_DATA_ROOT to be set.",
    )
    group.add_argument(
        "--input-dir",
        type=str,
        help="Direct path to the dataset directory. "
        "Use when there is no GEO accession (lab data, non-GEO source).",
    )
    parser.add_argument(
        "--name",
        type=str,
        default=None,
        help="Dataset identifier when using --input-dir (default: basename of --input-dir).",
    )
    parser.add_argument(
        "--data-root",
        type=str,
        default=None,
        help="Override FUXI_DATA_ROOT (default: from environment variable)",
    )
    parser.add_argument(
        "--query-ncbi",
        action="store_true",
        help="Query NCBI E-utilities for SuperSeries metadata (requires internet)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Detect and report only, no files written",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing dataset.yaml and config files",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Skip archive extraction (assume already extracted)",
    )
    parser.add_argument(
        "--modality",
        type=str,
        choices=["rna", "atac", "spatial", "multiome"],
        default=None,
        help="Force a specific modality (default: auto-detect)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Override output base directory "
        "(default: <repo>/projects/{modality}/{GSE_ID}/). "
        "Use to direct output to a temp location for review.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Auto-download GSE data from NCBI GEO before preprocessing "
        "(requires wget or curl; uses FUXI_DATA_ROOT).",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose output (show all detected files)",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Minimal output (errors only)",
    )

    args = parser.parse_args()

    exit_code = run_preprocess(
        gse_id=args.gse,
        input_dir=args.input_dir,
        dataset_name=args.name,
        data_root=args.data_root,
        query_ncbi=args.query_ncbi,
        download=args.download,
        dry_run=args.dry_run,
        force=args.force,
        no_extract=args.no_extract,
        modality=args.modality,
        output_dir=args.output_dir,
        verbose=args.verbose,
        quiet=args.quiet,
    )
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
