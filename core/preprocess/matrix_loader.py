#!/usr/bin/env python3
"""
matrix_loader.py — Phase 5 of the Fuxi preprocessing pipeline.

Detects primary matrix formats and generates ``config_GSE_ID.yaml``
files from config templates.
"""

import os
import sys
from typing import Optional

# Add repo root to sys.path (consistent with all step scripts)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from core.preprocess import format_detector as fd

# ── Template mapping: format → template file ──────────────────────────

TEMPLATE_MAP = {
    "10X_h5": "config_10X_h5.yaml",
    "10X_mtx": "config_10X_mtx.yaml",
    "csv_matrix": "config_csv_matrix.yaml",
    "h5ad": "config_10X_h5.yaml",  # reuse 10X_h5 template
    "10x_fragments": "config_fragments.yaml",
    "10x_peak_h5": "config_fragments.yaml",  # reuse ATAC template
    # ── Bulk entries ─────────────────────────────
    "count_matrix": "config_bulk.yaml",
    "tpm_matrix": "config_bulk.yaml",
    "bulk_h5ad": "config_bulk.yaml",
}


# ── Shared path helpers ──────────────────────────────────────────────


def _resolve_repo_root() -> str:
    """Return the absolute path to the repository root.

    matrix_loader.py is at: <repo>/core/preprocess/matrix_loader.py
    So __file__'s dirname → up 2 = core → up 1 = repo root.
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))  # .../core/preprocess
    return os.path.dirname(os.path.dirname(this_dir))  # .../ (repo root)


def _resolve_template_dir() -> str:
    """Return the path to the templates/config_templates/ directory."""
    return os.path.join(_resolve_repo_root(), "templates", "config_templates")


def _resolve_project_dir(modality: str, gse_id: str, output_dir: Optional[str] = None) -> str:
    """Return the output directory for a project config.

    Args:
        modality: 'rna', 'atac', 'spatial', or 'multiome'.
        gse_id:   GEO accession ID.
        output_dir: If set, use this as the base output directory instead
                    of the repo's projects/ tree.

    Returns:
        E.g.: <output_dir>/rna/GSE12345/ or projects/rna/GSE12345/
    """
    if output_dir:
        return os.path.join(output_dir, modality, gse_id)
    return os.path.join(_resolve_repo_root(), "projects", modality, gse_id)


# ═══════════════════════════════════════════════════════════════════════
#  Format / config helpers
# ═══════════════════════════════════════════════════════════════════════


def _detect_primary_format(classification: dict, modality: str = "") -> str:
    """Determine the primary data format from classification results.

    When *modality* is 'atac' and the classification matches both RNA and
    ATAC patterns, prefer the ATAC-specific format.
    """
    if modality == "atac":
        if classification.get("fragment_dirs"):
            return "10x_fragments"
        if classification.get("tenx_peak_dirs"):
            return "10x_peak_h5"

    if classification.get("tenx_h5_dirs"):
        return "10X_h5"
    if classification.get("tenx_mtx_dirs"):
        return "10X_mtx"
    if classification.get("fragment_dirs"):
        return "10x_fragments"
    if classification.get("tenx_peak_dirs"):
        return "10x_peak_h5"
    if classification.get("h5ad_files"):
        return "h5ad"
    if classification.get("csv_files"):
        return "csv_matrix"
    return "unknown"


def _fill_template(template_text: str, replacements: dict) -> str:
    """Replace {{KEY}} placeholders in *template_text* with values from *replacements*."""
    result = template_text
    for key, value in replacements.items():
        result = result.replace("{{" + key + "}}", str(value))
    return result


# ═══════════════════════════════════════════════════════════════════════
#  Post-process: inject paper-derived CFG fields
# ═══════════════════════════════════════════════════════════════════════


def _post_process_yaml(
    config_path: str, paper_context: dict, inject: Optional[dict] = None
) -> None:
    """Append paper-derived fields to a generated YAML config.

    Since the output is YAML (not Python), AST manipulation is no longer
    needed.  Instead we append additional YAML key-value lines.

    Args:
        config_path:   Path to the generated ``config_GSE_ID.yaml`` file.
        paper_context: Dict with optional keys ``features``, ``is_nuclei``,
                       ``tissue_kb``, ``tissue_ontology``.
        inject:        Optional dict of arbitrary key/value pairs to append.
    """
    if not paper_context and not inject:
        return

    lines_to_append: list[str] = []

    # -- paper_context (existing behaviour) --
    if paper_context:
        features = paper_context.get("features")
        if features is not None:
            genes_yaml = ", ".join(repr(g) for g in list(features))
            lines_to_append.append("marker:")
            lines_to_append.append(f"  marker_dict: {{extracted: [{genes_yaml}]}}")

        if paper_context.get("is_nuclei"):
            lines_to_append.append("qc:")
            lines_to_append.append("  is_nuclei: true")

        for key in ("tissue_kb", "tissue_ontology"):
            val = paper_context.get(key)
            if val is not None:
                lines_to_append.append(f"{key}: {repr(str(val))}")

    # -- inject dict (arbitrary key/value pairs) --
    if inject:
        for key, val in inject.items():
            lines_to_append.append(f"{key}: {repr(val)}")

    if not lines_to_append:
        return

    with open(config_path, "a", encoding="utf-8") as f:
        f.write("\n")
        f.write("\n".join(lines_to_append))
        f.write("\n")


# ═══════════════════════════════════════════════════════════════════════
#  Phase 5: Config generation
# ═══════════════════════════════════════════════════════════════════════


def generate_config(
    gse_id: str,
    modality: str,
    classification: dict,
    file_list: list[str],
    output_dir: str,
    data_root: Optional[str] = None,
    input_dir_override: Optional[str] = None,
    superseries_info: Optional[dict] = None,
    paper_context: Optional[dict] = None,
    dry_run: bool = False,
    force: bool = False,
) -> Optional[str]:
    """Generate a config_GSE_ID.py file.

    Args:
        gse_id:            GEO accession ID.
        modality:          'rna', 'atac', 'spatial', or 'multiome'.
        classification:    Format classification dict from format_detector.
        file_list:         List of file paths in the dataset.
        output_dir:        Output directory for the generated config.
        data_root:         Root data directory (default: FUXI_DATA_ROOT env var).
        input_dir_override: Override input directory path.
        superseries_info:  SuperSeries metadata dict (optional).
        paper_context:     Optional dict with paper-derived metadata values
                           (species, tissue, expression_type, genome,
                           assay_type, features, is_nuclei, etc.).
        dry_run:           Report only, don't write files.
        force:             Overwrite existing files.

    Returns the path to the generated config, or None.
    """
    data_format = _detect_primary_format(classification, modality)
    template_name = TEMPLATE_MAP.get(data_format)
    if not template_name:
        print(f"  [WARNING] No template for format '{data_format}' — skipping config generation")
        return None

    template_path = os.path.join(_resolve_template_dir(), template_name)
    if not os.path.exists(template_path):
        print(f"  [WARNING] Template not found: {template_path}")
        return None

    with open(template_path, "r", encoding="utf-8") as f:
        template_text = f.read()

    # Collect file paths relative to the input directory
    gse_dir = input_dir_override or os.path.join(data_root, gse_id)

    # Gather replacements
    species = fd.guess_species(file_list)
    if species == "unknown" and superseries_info:
        ncbi_species = superseries_info.get("species", "")
        if ncbi_species:
            species = fd._normalise_species(ncbi_species)
    tissue = fd.guess_tissue(file_list)
    genome = fd.guess_genome(species) or "hg38"

    # Detect primary file paths
    mtx_dir = ""
    mtx_prefix = ""
    h5_dir = ""
    matrix_file = ""
    barcodes_file = ""
    features_file = ""
    fragment_file = ""

    # 10X MTX
    for d, files in classification.get("tenx_mtx_dirs", {}).items():
        mtx_dir = os.path.relpath(d, gse_dir) if os.path.isabs(d) else d  # noqa: F841
        basenames = [os.path.basename(f) for f in files]
        # Find the common prefix before matrix/barcodes/features
        stripped = [fd.strip_known_suffix(b) for b in basenames]
        if stripped:
            prefix = os.path.commonprefix(stripped).rstrip("_.-")
            if prefix:
                mtx_prefix = prefix
        break

    # 10X H5
    for d, files in classification.get("tenx_h5_dirs", {}).items():
        h5_dir = os.path.relpath(d, gse_dir) if os.path.isabs(d) else d  # noqa: F841
        break

    # CSV
    for f in classification.get("csv_files", []):
        matrix_file = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
        break
    for f in classification.get("metadata_files", []):
        rf = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
        b = os.path.basename(rf).lower()
        if "barcode" in b:
            barcodes_file = rf
        elif "feature" in b or "gene" in b:
            features_file = rf
        elif not barcodes_file:
            barcodes_file = rf

    # ATAC fragments
    for d, files in classification.get("fragment_dirs", {}).items():
        for f in files:
            if "fragment" in os.path.basename(f).lower():
                fragment_file = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
                break
        break

    # Detect expression type
    expression_type = fd.detect_expression_type(classification, file_list)

    replacements = {
        "MTX_PREFIX": mtx_prefix,
        "MATRIX_FILE": matrix_file,
        "BARCODES_FILE": barcodes_file,
        "FEATURES_FILE": features_file,
        "FRAGMENT_FILE": fragment_file,
        "TISSUE": tissue,
        "SPECIES": species,
        "GENOME": genome,
        "EXPRESSION_TYPE": expression_type,
    }

    # Override heuristic values with paper_context where present
    if paper_context:
        for key in ("species", "tissue", "expression_type", "genome", "assay_type"):
            if key in paper_context and paper_context[key] is not None:
                replacements[key.upper()] = str(paper_context[key])

    filled = _fill_template(template_text, replacements)

    os.makedirs(output_dir, exist_ok=True)
    config_path = os.path.join(output_dir, f"config_{gse_id}.yaml")

    if dry_run:
        print(f"  [DRY-RUN] Would write: {config_path}")
        return config_path

    if os.path.exists(config_path) and not force:
        print(f"  [SKIP] Config already exists: {config_path}")
        print("         Use --force to overwrite.")
        return config_path

    with open(config_path, "w", encoding="utf-8") as f:
        f.write(filled)

    # Post-process: append paper-derived fields (marker_dict, is_nuclei, etc.)
    if paper_context:
        _post_process_yaml(config_path, paper_context)

    print(f"  Written: {config_path}")
    return config_path
