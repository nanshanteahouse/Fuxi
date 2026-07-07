#!/usr/bin/env python3
"""
metadata_parser.py — Phase 4 of the Fuxi preprocessing pipeline.

Extracts metadata from detected file structures and generates
``dataset.yaml`` files in the project output directory.

Extracted from :mod:`core.preprocess.preprocessor` (Phase 4) to
keep the main preprocessor focused on orchestration.
"""

import os
import sys
from datetime import datetime
from typing import Optional

# Add repo root to sys.path (consistent with all step scripts)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

from core.preprocess import format_detector as fd
from core.dataset_schema import (
    DatasetMeta, ModalityEntry, SampleEntry, FileEntry,
    Resources, PipelineStatus, Meta, save_dataset,
)
from core.preprocess.matrix_loader import (
    _detect_primary_format,
    _resolve_project_dir,
    _resolve_repo_root,
    _resolve_template_dir,
)


# ── Shared path helpers ──────────────────────────────────────────────


def _resolve_input_dir(gse_id: str, data_root: Optional[str], input_dir: Optional[str]) -> str:
    """Return the absolute input directory, handling GEO vs custom modes."""
    if input_dir is not None:
        return os.path.abspath(input_dir)
    if data_root is None:
        from core.utils import data_root as get_data_root
        data_root = get_data_root()
    return os.path.join(data_root, gse_id)


# ═══════════════════════════════════════════════════════════════════════
#  Phase 4: dataset.yaml generation
# ═══════════════════════════════════════════════════════════════════════


def generate_dataset_yaml(gse_id: str,
                          modality: str,
                          superseries_info: dict,
                          classification: dict,
                          file_list: list[str],
                          output_dir: str,
                          data_root: Optional[str],
                          input_dir_override: Optional[str] = None,
                          paper_context: Optional[dict] = None,
                          dry_run: bool = False,
                          force: bool = False,
                          ncbi_assay_type: Optional[str] = None) -> Optional[str]:
    """Generate a dataset.yaml metadata file.

    Returns the path to the generated YAML, or None.
    """
    species = fd.guess_species(file_list)
    if species == 'unknown' and superseries_info:
        ncbi_species = superseries_info.get('species', '')
        if ncbi_species:
            species = fd._normalise_species(ncbi_species)
    # Override with paper_context values where present
    if paper_context:
        if paper_context.get('species'):
            species = fd._normalise_species(str(paper_context['species']))
    species_key = species  # already normalised by species functions or paper_context
    tissue = fd.guess_tissue(file_list)
    if paper_context and paper_context.get('tissue'):
        tissue = paper_context['tissue']
    data_format = _detect_primary_format(classification, modality)

    # Determine the base directory for relative-path computation
    if input_dir_override:
        gse_dir = os.path.abspath(input_dir_override)
    elif data_root:
        gse_dir = os.path.join(data_root, gse_id)
    else:
        # Fallback: use the first file's directory
        gse_dir = os.path.dirname(file_list[0]) if file_list else '.'

    # ── Build ModalityEntry ──
    mod_name_map = {
        '10X_h5': 'scRNA-seq',
        '10X_mtx': 'scRNA-seq',
        'csv_matrix': 'scRNA-seq',
        'h5ad': 'scRNA-seq',
        '10x_fragments': 'scATAC-seq',
        '10x_peak_h5': 'scATAC-seq',
    }
    mod_name = mod_name_map.get(data_format, 'unknown')
    # When the caller forces a modality, align the entry name.
    forced_major = ''
    if modality in ('rna',):
        forced_major = 'scRNA-seq'
    elif modality in ('atac',):
        forced_major = 'scATAC-seq'
    elif modality in ('spatial',):
        forced_major = 'spatial_transcriptomics'
    if forced_major and mod_name not in (forced_major, 'unknown'):
        mod_name = forced_major

    modality_entry = ModalityEntry(
        name=mod_name,
        status='downloaded',
        format=data_format,
        file_count=len(file_list),
        total_size_gb=0.0,
        assay_type=ncbi_assay_type,
    )

    # ── Build SampleEntry list ──
    sample_groups = fd.group_files_by_sample(file_list)
    # Resolve the base data directory for path computations
    if input_dir_override:
        gse_dir = os.path.abspath(input_dir_override)
    elif data_root:
        gse_dir = os.path.join(data_root, gse_id)
    else:
        gse_dir = os.path.dirname(file_list[0]) if file_list else '.'
    samples = []
    for sample_id, files in sorted(sample_groups.items()):
        # Sanitize sample_id: use basename if it looks like a path, strip gse_dir prefix
        if os.path.isabs(sample_id) or '/' in sample_id or '\\' in sample_id:
            # It's a full path — derive a readable sample name
            sample_id = os.path.basename(sample_id.rstrip('/\\')) or sample_id
        # If sample_id is the gse_id itself, use 'all' as the sample name
        if sample_id.upper() == gse_id.upper():
            sample_id = 'all'
        # Classify files per sample
        rna_entries = []
        atac_entries = []
        for f in files:
            rel = os.path.relpath(f, gse_dir) if os.path.isabs(f) else f
            entry = FileEntry(file=rel, format='auto')
            b = os.path.basename(f).lower()
            if any(p in b for p in ('fragment', 'atac', 'peak', 'motif')):
                atac_entries.append(entry)
            else:
                rna_entries.append(entry)

        # Determine placement from modality name
        is_rna = mod_name == 'scRNA-seq'
        is_atac = mod_name == 'scATAC-seq'
        samples.append(SampleEntry(
            id=sample_id,
            label='',
            rna=rna_entries if is_rna else ([] if is_atac else rna_entries),
            atac=atac_entries if is_atac else ([] if is_rna else atac_entries),
            species=species if species != 'unknown' else None,
        ))

    # ── Build subseries list (if SuperSeries) ──
    subseries = []
    if superseries_info.get('is_superseries'):
        for child_acc in superseries_info.get('child_accessions', []):
            subseries.append({
                'id': child_acc,
                'title': '',
                'modality': mod_name,
            })
        # Also add directory-based subseries
        for dname in superseries_info.get('subseries_dirs', []):
            existing = {s.get('id') for s in subseries}
            if dname not in existing:
                subseries.append({
                    'id': dname,
                    'title': '',
                    'modality': 'unknown',
                })

    # ── Assemble DatasetMeta ──
    ds = DatasetMeta(
        id=gse_id,
        type='SuperSeries' if superseries_info.get('is_superseries') else 'SingleAccession',
        title=superseries_info.get('title', ''),
        species=species if species != 'unknown' else 'homo_sapiens',
        species_key=species_key if species_key != 'unknown' else 'human',
        tissue=tissue if tissue != 'unknown' else None,
        parent_superseries=None,
        modalities=[modality_entry],
        samples=samples,
        subseries=subseries,
        comparisons=[],
        assay_type=ncbi_assay_type,
        resources=Resources(
            genome=fd.guess_genome(species),
            technology='10x Genomics' if '10X' in data_format or '10x' in data_format else '',
        ),
        meta=Meta(
            created=datetime.now().isoformat(),
            generated_by='fuxi_preprocess',
            pipeline_status=PipelineStatus(),
        ),
    )

    os.makedirs(output_dir, exist_ok=True)
    yaml_path = os.path.join(output_dir, 'dataset.yaml')

    if dry_run:
        print(f"  [DRY-RUN] Would write: {yaml_path}")
        return yaml_path

    if os.path.exists(yaml_path) and not force:
        print(f"  [SKIP] dataset.yaml already exists: {yaml_path}")
        print(f"         Use --force to overwrite.")
        return yaml_path

    save_dataset(ds, yaml_path)
    print(f"  Written: {yaml_path}")
    return yaml_path
