#!/usr/bin/env python3
"""
matrix_loader.py — Phase 5 of the Fuxi preprocessing pipeline.

Detects primary matrix formats and generates ``config_GSE_ID.py``
files from config templates.

Extracted from :mod:`core.preprocess.preprocessor` (Phase 5) to
keep the main preprocessor focused on orchestration.
"""

import ast
import os
import sys
from typing import Optional

# Add repo root to sys.path (consistent with all step scripts)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

from core.preprocess import format_detector as fd


# ── Template mapping: format → template file ──────────────────────────

TEMPLATE_MAP = {
    '10X_h5':       'config_10X_h5.py',
    '10X_mtx':      'config_10X_mtx.py',
    'csv_matrix':   'config_csv_matrix.py',
    'h5ad':         'config_10X_h5.py',       # reuse 10X_h5 template
    '10x_fragments': 'config_fragments.py',
    '10x_peak_h5':  'config_fragments.py',    # reuse ATAC template
}


# ── Shared path helpers ──────────────────────────────────────────────


def _resolve_repo_root() -> str:
    """Return the absolute path to the repository root.

    matrix_loader.py is at: <repo>/core/preprocess/matrix_loader.py
    So __file__'s dirname → up 2 = core → up 1 = repo root.
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))  # .../core/preprocess
    return os.path.dirname(os.path.dirname(this_dir))       # .../ (repo root)


def _resolve_template_dir() -> str:
    """Return the path to the templates/config_templates/ directory."""
    return os.path.join(_resolve_repo_root(), 'templates', 'config_templates')


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
    return os.path.join(_resolve_repo_root(), 'projects', modality, gse_id)


# ═══════════════════════════════════════════════════════════════════════
#  Format / config helpers
# ═══════════════════════════════════════════════════════════════════════


def _detect_primary_format(classification: dict, modality: str = '') -> str:
    """Determine the primary data format from classification results.

    When *modality* is 'atac' and the classification matches both RNA and
    ATAC patterns, prefer the ATAC-specific format.
    """
    if modality == 'atac':
        if classification.get('fragment_dirs'):
            return '10x_fragments'
        if classification.get('tenx_peak_dirs'):
            return '10x_peak_h5'

    if classification.get('tenx_h5_dirs'):
        return '10X_h5'
    if classification.get('tenx_mtx_dirs'):
        return '10X_mtx'
    if classification.get('fragment_dirs'):
        return '10x_fragments'
    if classification.get('tenx_peak_dirs'):
        return '10x_peak_h5'
    if classification.get('h5ad_files'):
        return 'h5ad'
    if classification.get('csv_files'):
        return 'csv_matrix'
    return 'unknown'


def _fill_template(template_text: str, replacements: dict) -> str:
    """Replace {{KEY}} placeholders in *template_text* with values from *replacements*."""
    result = template_text
    for key, value in replacements.items():
        result = result.replace('{{' + key + '}}', str(value))
    return result


# ═══════════════════════════════════════════════════════════════════════
#  Post-process: inject paper-derived CFG fields
# ═══════════════════════════════════════════════════════════════════════


def _post_process_config(config_path: str, paper_context: dict,
                      inject: Optional[dict] = None) -> None:
    """Inject paper-derived CFG fields using ast manipulation.

    Uses Python's ``ast`` module to parse the generated config, find
    existing ``CFG.*`` assignments, and inject or replace fields
    (``marker_dict``, ``is_nuclei``, ``tissue_kb``, ``tissue_ontology``),
    plus additional arbitrary ``CFG.*`` values from the *inject* dict.

    **Idempotent**: if a field already exists, its value is replaced
    in-place rather than duplicate lines being appended.

    Args:
        config_path:   Path to the generated ``config_GSE_ID.py`` file.
        paper_context: Dict with optional keys ``features``, ``is_nuclei``,
                       ``tissue_kb``, ``tissue_ontology``.
        inject:        Optional dict of arbitrary ``CFG.*`` key/value pairs
                       to inject (e.g. ``{'sample_keep': ['GSM1'], 'subset_suffix': '_test'}``).
                       Uses ``repr(val)`` to preserve Python literals.
    """
    if not paper_context and not inject:
        return

    with open(config_path, 'r', encoding='utf-8') as f:
        source = f.read()

    try:
        tree = ast.parse(source)
    except SyntaxError:
        # If ast parsing fails, skip gracefully
        return

    # Determine which fields to inject
    injections: dict[str, str] = {}

    # -- paper_context (existing behaviour) --
    if paper_context:
        features = paper_context.get('features')
        if features is not None:
            marker_dict_repr = repr({'extracted': list(features)})
            injections['marker_dict'] = f"CFG.marker_dict = {marker_dict_repr}"

        if paper_context.get('is_nuclei'):
            injections['is_nuclei'] = 'CFG.is_nuclei = True'

        for key in ('tissue_kb', 'tissue_ontology'):
            val = paper_context.get(key)
            if val is not None:
                injections[key] = f"CFG.{key} = {repr(str(val))}"

    # -- inject dict (arbitrary CFG.* value injection) --
    if inject:
        for key, val in inject.items():
            injections[key] = f"CFG.{key} = {repr(val)}"

    if not injections:
        return

    # Find all CFG assignments and the last CFG assignment end line
    existing_attrs: dict[str, tuple[int, int]] = {}
    last_cfg_end: int = 0

    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if (
                    isinstance(target, ast.Attribute)
                    and isinstance(target.value, ast.Name)
                    and target.value.id == 'CFG'
                ):
                    attr = target.attr
                    start_line = node.lineno
                    end_line = getattr(node, 'end_lineno', start_line)
                    existing_attrs[attr] = (start_line, end_line)
                    if end_line > last_cfg_end:
                        last_cfg_end = end_line

    lines = source.splitlines()

    # Phase 1: Replace existing fields (bottom-to-top to preserve indices)
    replace_ops = [(s, e, injections[a]) for a, (s, e) in existing_attrs.items() if a in injections]
    for start, end, new_line in sorted(replace_ops, key=lambda x: x[0], reverse=True):
        lines[start - 1:end] = [new_line]

    # Phase 2: Append new fields after the last CFG assignment
    new_fields = {k: v for k, v in injections.items() if k not in existing_attrs}
    if new_fields:
        # Re-parse modified source to find the new last CFG end
        modified_source = '\n'.join(lines)
        try:
            new_tree = ast.parse(modified_source)
            new_last_end: int = 0
            for node in ast.walk(new_tree):
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if (
                            isinstance(target, ast.Attribute)
                            and isinstance(target.value, ast.Name)
                            and target.value.id == 'CFG'
                        ):
                            end = getattr(node, 'end_lineno', node.lineno)
                            if end > new_last_end:
                                new_last_end = end
            insert_idx = new_last_end  # 1-indexed → line after last CFG assign
        except SyntaxError:
            insert_idx = len(lines)

        for new_line in new_fields.values():
            lines.insert(insert_idx, new_line)
            insert_idx += 1

    with open(config_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


# ═══════════════════════════════════════════════════════════════════════
#  Phase 5: Config generation
# ═══════════════════════════════════════════════════════════════════════


def generate_config(gse_id: str,
                    modality: str,
                    classification: dict,
                    file_list: list[str],
                    output_dir: str,
                    data_root: Optional[str] = None,
                    input_dir_override: Optional[str] = None,
                    superseries_info: Optional[dict] = None,
                    paper_context: Optional[dict] = None,
                    dry_run: bool = False,
                    force: bool = False) -> Optional[str]:
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

    with open(template_path, 'r', encoding='utf-8') as f:
        template_text = f.read()

    # Collect file paths relative to the input directory
    gse_dir = input_dir_override or os.path.join(data_root, gse_id)
    rel_files = [os.path.relpath(f, gse_dir) for f in file_list]

    # Gather replacements
    species = fd.guess_species(file_list)
    if species == 'unknown' and superseries_info:
        ncbi_species = superseries_info.get('species', '')
        if ncbi_species:
            species = fd._normalise_species(ncbi_species)
    tissue = fd.guess_tissue(file_list)
    genome = fd.guess_genome(species) or 'hg38'

    # Detect primary file paths
    mtx_dir = ''
    mtx_prefix = ''
    h5_dir = ''
    matrix_file = ''
    barcodes_file = ''
    features_file = ''
    fragment_file = ''

    # 10X MTX
    for d, files in classification.get('tenx_mtx_dirs', {}).items():
        mtx_dir = os.path.relpath(d, gse_dir) if os.path.isabs(d) else d
        # Heuristic: strip trailing directory separator + common suffix
        mtx_dir_norm = os.path.basename(mtx_dir.rstrip('/\\')) or mtx_dir
        basenames = [os.path.basename(f) for f in files]
        # Find the common prefix before matrix/barcodes/features
        stripped = [fd.strip_known_suffix(b) for b in basenames]
        if stripped:
            prefix = os.path.commonprefix(stripped).rstrip('_.-')
            if prefix:
                mtx_prefix = prefix
        break

    # 10X H5
    for d, files in classification.get('tenx_h5_dirs', {}).items():
        h5_dir = os.path.relpath(d, gse_dir) if os.path.isabs(d) else d
        break

    # CSV
    for f in classification.get('csv_files', []):
        matrix_file = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
        break
    for f in classification.get('metadata_files', []):
        rf = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
        b = os.path.basename(rf).lower()
        if 'barcode' in b:
            barcodes_file = rf
        elif 'feature' in b or 'gene' in b:
            features_file = rf
        elif not barcodes_file:
            barcodes_file = rf

    # ATAC fragments
    for d, files in classification.get('fragment_dirs', {}).items():
        for f in files:
            if 'fragment' in os.path.basename(f).lower():
                fragment_file = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
                break
        break

    # Detect expression type
    expression_type = fd.detect_expression_type(classification, file_list)

    replacements = {
        'MTX_PREFIX': mtx_prefix,
        'MATRIX_FILE': matrix_file,
        'BARCODES_FILE': barcodes_file,
        'FEATURES_FILE': features_file,
        'FRAGMENT_FILE': fragment_file,
        'TISSUE': tissue,
        'SPECIES': species,
        'GENOME': genome,
        'EXPRESSION_TYPE': expression_type,
    }

    # Override heuristic values with paper_context where present
    if paper_context:
        for key in ('species', 'tissue', 'expression_type', 'genome', 'assay_type'):
            if key in paper_context and paper_context[key] is not None:
                replacements[key.upper()] = str(paper_context[key])

    filled = _fill_template(template_text, replacements)

    os.makedirs(output_dir, exist_ok=True)
    config_path = os.path.join(output_dir, f'config_{gse_id}.py')

    if dry_run:
        print(f"  [DRY-RUN] Would write: {config_path}")
        return config_path

    if os.path.exists(config_path) and not force:
        print(f"  [SKIP] Config already exists: {config_path}")
        print(f"         Use --force to overwrite.")
        return config_path

    with open(config_path, 'w', encoding='utf-8') as f:
        f.write(filled)

    # Post-process: inject paper-derived CFG fields (marker_dict, is_nuclei, etc.)
    if paper_context:
        _post_process_config(config_path, paper_context)

    print(f"  Written: {config_path}")
    return config_path
