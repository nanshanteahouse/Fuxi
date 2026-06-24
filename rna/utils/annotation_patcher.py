"""
utils/annotation_patcher.py — Safe annotation patch utility (v3.1.0+).

Provides :func:`apply_annotation_patches` which updates specific cluster
annotations on an already-annotated AnnData while preserving
``marker_validation`` metadata for unchanged clusters and recalculating it
for changed ones.

Usage from a one-liner patch script::

    import scanpy as sc
    from core.config import CFG
    from rna.annotation_standardizer import StandardOntology
    from rna.utils.annotation_patcher import apply_annotation_patches

    CFG.resolve_paths()
    adata = sc.read(CFG.annotated_h5ad)
    std = StandardOntology(CFG.tissue_ontology or CFG.tissue_kb)
    apply_annotation_patches(
        adata, {"5": "Rod Photoreceptor"}, cfg=CFG, std=std,
    )
    adata.write(CFG.annotated_h5ad)

This preserves pipeline metadata (marker_validation, annot_evidence, etc.)
that a bare ``.obs['cell_type']`` overwrite would lose.
"""

import json
import logging
import os
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


def apply_annotation_patches(
    adata,
    patches: dict,  # {cluster_str: new_cell_type_str}
    cfg=None,
    std: Optional[object] = None,
    log: Optional[logging.Logger] = None,
):
    """Update specific clusters' annotations while preserving pipeline metadata.

    Parameters
    ----------
    adata : AnnData
        Annotated data with ``.obs['cell_type']``, ``.obs['leiden']``, etc.
        Modified **in-place**.
    patches : dict
        Mapping of ``{cluster_id: new_cell_type_name}`` for clusters to update.
        Keys should be strings matching ``adata.obs['leiden']`` values.
    cfg : Config or None
        Pipeline config.  If provided, re-writes annotation CSV and quality
        report to ``cfg.table_dir``.
    std : StandardOntology or None
        If provided, recalculates ``marker_validation`` for patched clusters
        via :meth:`StandardOntology.validate`.
    log : Logger or None
        Logger instance.  Uses module-level logger when ``None``.

    Returns
    -------
    AnnData
        Modified *adata* (same object, in-place).
    """
    _log = log or logger

    if not patches:
        _log.info("No patches to apply.")
        return adata

    leiden_str = adata.obs['leiden'].astype(str)

    for cluster_id, new_ct in patches.items():
        mask = leiden_str == str(cluster_id)
        n_cells = mask.sum()
        if n_cells == 0:
            _log.warning(
                "Cluster '%s' not found in data, skipping.", cluster_id,
            )
            continue

        old_ct = adata.obs.loc[mask, 'cell_type'].iloc[0]
        _log.info(
            "Patching cluster %s: '%s' -> '%s' (%d cells)",
            cluster_id, old_ct, new_ct, n_cells,
        )

        # Update core annotation columns
        adata.obs.loc[mask, 'cell_type'] = new_ct
        if 'annot_confidence' in adata.obs:
            adata.obs.loc[mask, 'annot_confidence'] = 'patched'
        if 'annot_method' in adata.obs:
            adata.obs.loc[mask, 'annot_method'] = 'manual_patch'

        # Preserve old reasoning as context
        if 'annot_reasoning' in adata.obs:
            old_reasoning = adata.obs.loc[mask, 'annot_reasoning'].iloc[0]
            adata.obs.loc[mask, 'annot_reasoning'] = (
                f"[PATCHED from '{old_ct}'] {old_reasoning}"
            )

    # ── Recalculate marker_validation ───────────────────────────────────
    if std is not None:
        try:
            validation_results = std.validate(adata)
            validation_map = {
                r['cluster']: r['status'] for r in validation_results
            }
            adata.obs['marker_validation'] = leiden_str.map(
                lambda c: validation_map.get(c, "NO_ONTOLOGY")
            )
            _log.info(
                "marker_validation recalculated: %d/%d PASS",
                sum(1 for r in validation_results if r['status'] == 'PASS'),
                len(validation_results),
            )
        except Exception as exc:
            _log.warning(
                "marker_validation recalculation failed: %s — "
                "validation column may be stale", exc,
            )

    # ── Re-write annotation CSV ─────────────────────────────────────────
    if cfg is not None:
        _rewrite_annotation_csv(adata, cfg, _log)
        _rewrite_quality_report(adata, cfg, _log)

    return adata


# ═══════════════════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════════════════


def _rewrite_annotation_csv(adata, cfg, log):
    """Re-write ``cell_type_annotations.csv`` from current adata.obs."""
    leiden_str = adata.obs['leiden'].astype(str)
    cluster_ids = sorted(
        adata.obs['leiden'].unique(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )
    records = []
    for cl in cluster_ids:
        mask = leiden_str == str(cl)
        records.append({
            'cluster': str(cl),
            'cell_type': adata.obs.loc[mask, 'cell_type'].iloc[0],
            'confidence': (
                adata.obs.loc[mask, 'annot_confidence'].iloc[0]
                if 'annot_confidence' in adata.obs else 'N/A'
            ),
            'method': (
                adata.obs.loc[mask, 'annot_method'].iloc[0]
                if 'annot_method' in adata.obs else 'N/A'
            ),
            'reasoning': (
                adata.obs.loc[mask, 'annot_reasoning'].iloc[0]
                if 'annot_reasoning' in adata.obs else ''
            ),
        })
    ann_df = pd.DataFrame(records)
    ann_csv = os.path.join(cfg.table_dir, 'cell_type_annotations.csv')
    ann_df.to_csv(ann_csv, index=False)
    log.info("Annotation table re-written: %s", ann_csv)


def _rewrite_quality_report(adata, cfg, log):
    """Re-write ``05_annotation_quality.json`` from current adata.obs."""
    if 'marker_validation' not in adata.obs:
        log.info("No marker_validation column — skipping quality report.")
        return
    pass_cells = (adata.obs['marker_validation'] == 'PASS').sum()
    pass_rate = pass_cells / max(adata.n_obs, 1)
    quality = {
        "pass_rate": round(pass_rate, 4),
        "total_clusters": adata.obs['leiden'].nunique(),
        "note": "Generated by annotation_patcher — some clusters manually patched.",
    }
    q_path = os.path.join(cfg.table_dir, '05_annotation_quality.json')
    with open(q_path, 'w') as f:
        json.dump(quality, f, indent=2, ensure_ascii=False)
    log.info("Quality report re-written: %s", q_path)
