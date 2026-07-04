"""Cross-modality helpers — RNA result discovery for ATAC/spatial integration."""

import logging
import os
from typing import Dict, List, Optional

from core.utils._path import repo_root


def find_rna_h5ad(cfg=None, dataset_id: str = None, log=None) -> Optional[str]:
    """Auto-discover the RNA annotated h5ad for ATAC Step 09 integration.

    Invoked automatically by run_pipeline.py when --modality atac and
    cfg.rna_h5ad is empty.  Also callable from 09_integrate.py as a
    last-resort fallback.

    Search order (first existing file wins):
      1. cfg.rna_h5ad (explicitly set by user — always honoured first)
      2. projects/rna/{dataset_id}/results/h5ad/05_annotated.h5ad
      3. projects/rna/{dataset_id}/results/h5ad/04_clustered.h5ad
      4. projects/rna/{dataset_id}/results/h5ad/03_integrated.h5ad

    ``dataset_id`` is inferred from ``cfg.project_dir`` when omitted.
    """
    if log is None:
        log = logging.getLogger(__name__)

    # ── 1. Explicit setting always wins ──────────────────────────────
    if cfg is not None and getattr(cfg, 'rna_h5ad', ''):
        return cfg.rna_h5ad

    # ── 2. Derive dataset ID from config ─────────────────────────────
    if dataset_id is None and cfg is not None:
        proj = getattr(cfg, 'project_dir', '')
        if proj:
            dataset_id = os.path.basename(os.path.normpath(proj))
    if not dataset_id:
        return None

    # ── 3. Locate repo root ──────────────────────────────────────────
    repo = repo_root()  # e.g. /mnt/d/Projects/Fuxi or D:\Projects\Fuxi

    # Candidate h5ad files (most complete first)
    candidates = [
        "05_annotated.h5ad",
        "04_clustered.h5ad",
        "03_integrated.h5ad",
    ]

    rna_project_dir = os.path.join(repo, "projects", "rna", dataset_id)
    if not os.path.isdir(rna_project_dir):
        log.debug("find_rna_h5ad: no RNA project dir at %s", rna_project_dir)
        return None

    h5ad_dir = os.path.join(rna_project_dir, "results", "h5ad")
    if not os.path.isdir(h5ad_dir):
        # Fallback: try legacy structure where h5ad is in results/ directly
        h5ad_dir = os.path.join(rna_project_dir, "results")

    for fname in candidates:
        full = os.path.join(h5ad_dir, fname)
        if os.path.isfile(full) and os.path.getsize(full) > 0:
            log.info("Auto-discovered RNA h5ad: %s", full)
            return full

    log.debug("find_rna_h5ad: no RNA h5ad found in %s", h5ad_dir)
    return None


def find_rna_marker_csv(cfg=None, dataset_id: str = None, log=None) -> Optional[str]:
    """Auto-discover the scRNA marker CSV for spatial annotation (Phase 1).

    Invoked from ``spatial/steps/05_annotate.py`` when ``cfg.rna_ref``
    is set.  Returns the path to a ``marker_genes_per_group_cell_type.csv``
    file produced by RNA Step 07.

    Search order (first existing file wins):
      1. ``cfg.rna_ref`` — if path to an existing CSV file, use directly
      2. ``projects/rna/{id}/results/tables/marker_genes_per_group_cell_type.csv``
      3. ``projects/rna/{id}/results/tables/marker_genes_per_group.csv``  (fallback)

    ``dataset_id`` is inferred from ``cfg.rna_ref`` when that is a bare
    ID string, otherwise from ``cfg.project_dir``.
    """
    if log is None:
        log = logging.getLogger(__name__)

    # ── 1. Explicit path to CSV — use directly ────────────────────────
    if cfg is not None and getattr(cfg, 'rna_ref', ''):
        ref = cfg.rna_ref
        if ref.endswith('.csv') and os.path.isfile(ref):
            return ref
        # If rna_ref is a directory path, use its basename as dataset_id
        if os.path.isdir(ref):
            dataset_id = os.path.basename(os.path.normpath(ref))

    # ── 2. Resolve dataset_id from rna_ref string or project_dir ─────
    if dataset_id is None and cfg is not None:
        ref = getattr(cfg, 'rna_ref', '')
        if ref and not os.path.sep in ref and '/' not in ref:
            dataset_id = ref   # bare ID like "GSE235585"
    if dataset_id is None and cfg is not None:
        proj = getattr(cfg, 'project_dir', '')
        if proj:
            dataset_id = os.path.basename(os.path.normpath(proj))
    if not dataset_id:
        log.debug("find_rna_marker_csv: cannot determine dataset_id")
        return None

    # ── 3. Locate repo root ───────────────────────────────────────────
    repo = repo_root()

    # Candidate tables (cell-type markers preferred over leiden markers)
    candidates = [
        "marker_genes_per_group_cell_type.csv",
        "marker_genes_per_group.csv",
    ]

    rna_project_dir = os.path.join(repo, "projects", "rna", dataset_id)
    if not os.path.isdir(rna_project_dir):
        log.debug("find_rna_marker_csv: no RNA project dir at %s", rna_project_dir)
        return None

    tables_dir = os.path.join(rna_project_dir, "results", "tables")
    for fname in candidates:
        full = os.path.join(tables_dir, fname)
        if os.path.isfile(full) and os.path.getsize(full) > 0:
            log.info("Auto-discovered RNA marker CSV: %s", full)
            return full

    log.debug("find_rna_marker_csv: no marker CSV found in %s", tables_dir)
    return None


def load_scRNA_markers(
    csv_path: str,
    top_n: int = 10,
    pval_threshold: float = 0.05,
    logfc_min: float = 0.0,
    log=None,
) -> Dict[str, List[str]]:
    """Extract per-cell-type top-N marker genes from scRNA DE results.

    Reads a ``marker_genes_per_group_cell_type.csv`` file (produced by
    RNA Step 07) and converts it to ``CFG.marker_dict`` format.

    Parameters
    ----------
    csv_path : str
        Path to the RNA marker CSV (columns: ``group``, ``names``,
        ``scores``, ``logfoldchanges``, ``pvals``, ``pvals_adj``,
        ``pct_nz_group``, ``pct_nz_reference``).
    top_n : int
        Number of top marker genes to extract per cell type.
    pval_threshold : float
        Only consider genes with ``pvals_adj`` < this value.
    logfc_min : float
        Minimum ``logfoldchanges`` to include (0.0 = positive LFC only).

    Returns
    -------
    Dict[str, List[str]]
        ``{cell_type: [gene1, gene2, ..., geneN]}`` suitable for
        ``CFG.marker_dict``.
    """
    import pandas as pd

    if log is None:
        log = logging.getLogger(__name__)

    df = pd.read_csv(csv_path)
    log.info(
        "load_scRNA_markers: read %d rows from %s",
        len(df), os.path.basename(csv_path),
    )

    # Filter: significant + positively enriched
    mask = (
        (df['pvals_adj'] < pval_threshold)
        & (df['logfoldchanges'] > logfc_min)
    )
    df_filt = df[mask]
    log.info(
        "  after filtering (pval<%.0e, logfc>%.1f): %d rows retained",
        pval_threshold, logfc_min, len(df_filt),
    )

    # Group by cell type, sort by score descending, take top-N
    marker_dict: Dict[str, List[str]] = {}
    for group, group_df in df_filt.groupby('group', observed=True):
        top_genes = (
            group_df
            .sort_values('scores', ascending=False)
            .head(top_n)['names']
            .tolist()
        )
        if top_genes:
            marker_dict[group] = top_genes
            log.debug("  %s: %d markers → %s", group, len(top_genes),
                      top_genes[:5])

    log.info(
        "load_scRNA_markers: extracted markers for %d cell types "
        "(top_n=%d, pval<%.0e, logfc>%.1f)",
        len(marker_dict), top_n, pval_threshold, logfc_min,
    )
    return marker_dict
