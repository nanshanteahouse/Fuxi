"""Pipeline output parsers — read QC reports, marker genes, enrichment results.

All functions accept a ``config_path`` (path to a YAML config file), load
the :class:`~core.config.schema.Config` via :func:`load_yaml_config`, resolve
paths relative to the config location, and return structured data.

Missing files are handled gracefully — every function returns an empty result
(empty dict / empty list) instead of raising.
"""

from __future__ import annotations

import glob
import os
from typing import Any

import pandas as pd

from core.tui.backends.config import load_yaml_config


# ═══════════════════════════════════════════════════════════════════════
#  QC Report
# ═══════════════════════════════════════════════════════════════════════


def parse_qc_report(config_path: str) -> dict[str, Any]:
    """Read the QC summary CSV from *table_dir* and return a flat dict.

    Candidate filenames checked (first match wins):

    * ``qc_summary.csv``
    * ``qc_report.csv``

    Expected columns: ``n_cells``, ``median_genes``, ``median_umis``,
    ``median_pct_mito``, ``cells_before``, ``cells_after``, etc.

    Returns
    -------
    dict
        Empty dict when no QC file is found or the file is empty.
    """
    cfg = load_yaml_config(config_path)
    table_dir = cfg.table_dir

    candidates = ["qc_summary.csv", "qc_report.csv"]
    path = _first_existing(table_dir, candidates)
    if path is None:
        return {}

    try:
        df = pd.read_csv(path)
    except Exception:
        return {}

    if df.empty:
        return {}

    # If the CSV has one row, return it as a flat dict.
    if len(df) == 1:
        return df.iloc[0].dropna().to_dict()

    # Multiple rows — return as {"summary": [...], "columns": [...]}
    return {
        "summary": df.fillna("").to_dict(orient="records"),
        "columns": list(df.columns),
    }


# ═══════════════════════════════════════════════════════════════════════
#  Marker Genes
# ═══════════════════════════════════════════════════════════════════════


def parse_marker_genes(config_path: str) -> list[dict[str, Any]]:
    """Read the marker-gene-per-group CSV from *table_dir*.

    Candidates checked (first match wins):

    * ``marker_genes_per_group.csv``
    * ``marker_genes_per_group_filtered.csv``

    Expected columns: ``group``, ``names`` (gene symbol), ``scores``,
    ``logfoldchanges``, ``pvals``, ``pvals_adj``.

    Returns
    -------
    list[dict]
        Each dict is one marker-gene row.  Empty list when the file is
        missing or unreadable.
    """
    cfg = load_yaml_config(config_path)
    table_dir = cfg.table_dir

    candidates = [
        "marker_genes_per_group.csv",
        "marker_genes_per_group_filtered.csv",
    ]
    path = _first_existing(table_dir, candidates)
    if path is None:
        return []

    try:
        df = pd.read_csv(path)
    except Exception:
        return []

    if df.empty:
        return []

    # Normalise common column name variations
    _rename_cols(df)

    return df.fillna("").to_dict(orient="records")


# ═══════════════════════════════════════════════════════════════════════
#  Enrichment Results
# ═══════════════════════════════════════════════════════════════════════


def parse_enrichment(config_path: str) -> list[dict[str, Any]]:
    """Read enrichment CSV files from the pipeline output.

    Search locations (in order):

    1. ``{table_dir}/enrichment_results.csv`` — ATAC pipeline single-file
       output.
    2. ``{table_dir}/enrichment_*.csv`` — generic enrichment wildcard.
    3. ``{table_dir}/09_enrichment/ora_*.csv`` — RNA ORA per-gene-set.
    4. ``{table_dir}/09_enrichment/prerank_*.csv`` — RNA GSEA per-gene-set.
    5. ``{table_dir}/nhood_enrichment_zscore.csv`` — spatial
       neighbourhood enrichment.

    Each result dict includes a ``_source`` key indicating the filename
    (basename) so the caller can distinguish multiple enrichment outputs.

    Returns
    -------
    list[dict]
        Empty list when no enrichment files are found.
    """
    cfg = load_yaml_config(config_path)
    table_dir = cfg.table_dir

    results: list[dict[str, Any]] = []

    # 1. ATAC-style single enrichment_results.csv
    single_path = os.path.join(table_dir, "enrichment_results.csv")
    _append_csv(single_path, results, source_tag="enrichment_results.csv")

    # 2. Generic enrichment_*.csv
    for fp in sorted(glob.glob(os.path.join(table_dir, "enrichment_*.csv"))):
        if fp == single_path:
            continue  # already loaded above
        _append_csv(fp, results, source_tag=os.path.basename(fp))

    # 3. RNA ORA per gene-set
    ora_dir = os.path.join(table_dir, "09_enrichment")
    if os.path.isdir(ora_dir):
        for fp in sorted(glob.glob(os.path.join(ora_dir, "ora_*.csv"))):
            _append_csv(fp, results, source_tag=f"09_enrichment/{os.path.basename(fp)}")
        for fp in sorted(glob.glob(os.path.join(ora_dir, "prerank_*.csv"))):
            _append_csv(fp, results, source_tag=f"09_enrichment/{os.path.basename(fp)}")

    # 4. Spatial neighbourhood enrichment
    spatial_path = os.path.join(table_dir, "nhood_enrichment_zscore.csv")
    _append_csv(spatial_path, results, source_tag="nhood_enrichment_zscore.csv")

    return results


# ═══════════════════════════════════════════════════════════════════════
#  Report discovery
# ═══════════════════════════════════════════════════════════════════════


def list_available_reports(config_path: str) -> list[str]:
    """Return a list of recognised report filenames that exist on disk.

    Checks all the file locations that the other ``parse_*`` functions
    know about.

    Returns
    -------
    list[str]
        Human-readable labels, e.g. ``["QC report", "Marker genes",
        "Enrichment (enrichment_results.csv)"]``.
        Empty list if nothing is found.
    """
    cfg = load_yaml_config(config_path)
    table_dir = cfg.table_dir

    labels: list[str] = []

    # QC
    if _first_existing(table_dir, ["qc_summary.csv", "qc_report.csv"]):
        labels.append("QC report")

    # Marker genes
    if _first_existing(table_dir, [
        "marker_genes_per_group.csv",
        "marker_genes_per_group_filtered.csv",
    ]):
        labels.append("Marker genes")

    # Enrichment
    if os.path.isfile(os.path.join(table_dir, "enrichment_results.csv")):
        labels.append("Enrichment (enrichment_results.csv)")

    for fp in sorted(glob.glob(os.path.join(table_dir, "enrichment_*.csv"))):
        labels.append(f"Enrichment ({os.path.basename(fp)})")

    ora_dir = os.path.join(table_dir, "09_enrichment")
    if os.path.isdir(ora_dir):
        for fp in sorted(glob.glob(os.path.join(ora_dir, "ora_*.csv"))):
            labels.append(f"Enrichment ORA ({os.path.basename(fp)})")
        for fp in sorted(glob.glob(os.path.join(ora_dir, "prerank_*.csv"))):
            labels.append(f"Enrichment GSEA ({os.path.basename(fp)})")

    if os.path.isfile(os.path.join(table_dir, "nhood_enrichment_zscore.csv")):
        labels.append("Spatial neighbourhood enrichment")

    return labels


# ═══════════════════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════════════════


def _first_existing(directory: str, filenames: list[str]) -> str | None:
    """Return the full path of the first existing file in *directory*
    from *filenames*, or ``None`` if none exist."""
    for name in filenames:
        path = os.path.join(directory, name)
        if os.path.isfile(path):
            return path
    return None


def _append_csv(
    path: str,
    target: list[dict[str, Any]],
    source_tag: str,
) -> None:
    """Read a CSV and append its rows (as dicts) to *target*.

    Each row gets a ``_source`` key set to *source_tag*.
    Silently returns if the file does not exist or cannot be parsed.
    """
    if not os.path.isfile(path):
        return
    try:
        df = pd.read_csv(path)
    except Exception:
        return
    if df.empty:
        return
    _rename_cols(df)
    records = df.fillna("").to_dict(orient="records")
    for rec in records:
        rec["_source"] = source_tag
    target.extend(records)


def _rename_cols(df: pd.DataFrame) -> None:
    """Normalise common column-name variations in-place."""
    renames = {
        "names": "gene",
        "Term": "pathway",
        "Adjusted P-value": "pval_adj",
        "P-value": "pval",
        "Overlap": "overlap",
        "Genes": "genes",
    }
    # Only rename columns that actually exist
    existing = {k: v for k, v in renames.items() if k in df.columns}
    if existing:
        df.rename(columns=existing, inplace=True)
