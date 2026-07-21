#!/usr/bin/env python3
"""anatomy.py — Anatomical adjacency loading and CCI filtering across modalities

``core/anatomy.py`` is the cross-modality entry point for anatomical adjacency data.
RNA, Spatial, and future ATAC CCI steps use these two functions to load adjacency
and filter/annotate LIANA results.

Functions
---------
load_adjacency
    Load adjacency from a custom CSV file or a tissue knowledge base.
filter_cci_by_adjacency
    Annotate or filter CCI results against an adjacency table.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import pandas as pd

# ── Module-level logger ────────────────────────────────────────────────
_log = logging.getLogger(__name__)


# ======================================================================
#  Public API
# ======================================================================


def load_adjacency(
    tissue: str = "",
    custom_file: str = "",
    log: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Load anatomical adjacency data.

    Resolution priority:

        1. If **custom_file** is non-empty **and** the file exists → load as CSV.
        2. Else if **tissue** is non-empty → call
           ``rna.tissue_ontologies.load_adjacency(tissue)``.
        3. Else → return an empty DataFrame with columns
           ``["source", "target", "adjacency_type"]``.

    Parameters
    ----------
    tissue : str
        Tissue identifier (e.g. ``"retina"``).  Used as a fallback when
        *custom_file* is not provided or is missing.
    custom_file : str
        Path to a custom adjacency CSV file.  Expected columns:
        ``source``, ``target``, ``adjacency_type``.  Takes priority over
        *tissue*.
    log : logging.Logger or None
        Optional logger.  If ``None`` the module-level logger is used.

    Returns
    -------
    pd.DataFrame
        Adjacency table with columns ``["source", "target", "adjacency_type"]``.
        Returns an empty DataFrame (same columns) when no source is available.

    Examples
    --------
    >>> adj = load_adjacency(tissue="retina")
    >>> adj = load_adjacency(custom_file="adjacencies/retina.csv")
    >>> adj = load_adjacency()  # empty table
    """
    logger = log or _log

    # ── Priority 1: custom CSV file ──────────────────────────────────
    if custom_file:
        if os.path.isfile(custom_file):
            df = pd.read_csv(custom_file)
            required = {"source", "target", "adjacency_type"}
            missing = required - set(df.columns)
            if missing:
                logger.warning(
                    "Custom adjacency file '%s' missing columns: %s. Returning empty DataFrame.",
                    custom_file,
                    sorted(missing),
                )
                return pd.DataFrame(columns=["source", "target", "adjacency_type"])

            logger.info(
                "Loaded adjacency from custom file: %s (%d rows)",
                custom_file,
                len(df),
            )
            return df[["source", "target", "adjacency_type"]]
        else:
            logger.warning(
                "Custom adjacency file not found: '%s'. Falling through to tissue default.",
                custom_file,
            )

    # ── Priority 2: tissue knowledge base ────────────────────────────
    if tissue:
        try:
            # Lazy import to avoid circular dependency
            from core.kb import load_adjacency as _load_tissue_adj

            df = _load_tissue_adj(tissue)
            if df is None or df.empty:
                logger.info(
                    "No adjacency found for tissue '%s' (empty KB result).",
                    tissue,
                )
                return pd.DataFrame(columns=["source", "target", "adjacency_type"])

            logger.info(
                "Loaded adjacency for tissue '%s': %d rows",
                tissue,
                len(df),
            )
            return df
        except (ImportError, AttributeError):
            logger.warning(
                "core.kb.load_adjacency() is not available yet "
                "for tissue '%s'. Returning empty DataFrame.",
                tissue,
            )
        except Exception as exc:
            logger.warning(
                "Failed to load adjacency for tissue '%s': %s",
                tissue,
                exc,
            )

    # ── Priority 3: empty fallback ───────────────────────────────────
    logger.info("No adjacency source provided. Returning empty DataFrame.")
    return pd.DataFrame(columns=["source", "target", "adjacency_type"])


def filter_cci_by_adjacency(
    lr_res: pd.DataFrame,
    adjacency: pd.DataFrame,
    mode: str = "off",
    adjacency_types: Optional[list[str]] = None,
    log: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Annotate or filter LIANA CCI results by anatomical adjacency.

    Parameters
    ----------
    lr_res : pd.DataFrame
        LIANA interaction results.  Must contain columns ``source`` and
        ``target``.  May contain additional columns such as ``ligand``,
        ``receptor``, ``pvalue``, etc.
    adjacency : pd.DataFrame
        Anatomical adjacency table with columns ``source``, ``target``,
        and ``adjacency_type``.
    mode : {"off", "soft", "hard"}
        - ``"off"`` — return *lr_res* unchanged (default, no copy).
        - ``"soft"`` — add ``adjacent`` (bool) and ``adjacency_type`` (str)
          columns to a **copy** of *lr_res*.
        - ``"hard"`` — return a **filtered copy** of *lr_res* containing
          only pairs present in *adjacency* (bidirectional match).
    adjacency_types : list of str or None
        If non-empty, only adjacency rows whose ``adjacency_type`` value
        is in this list are considered.
    log : logging.Logger or None
        Optional logger.  If ``None`` the module-level logger is used.

    Returns
    -------
    pd.DataFrame
        - ``"off"`` mode: the original *lr_res* object (no copy).
        - ``"soft"`` mode: a copy of *lr_res* with ``adjacent`` and
          ``adjacency_type`` columns added.
        - ``"hard"`` mode: a filtered copy of *lr_res*.

    Raises
    ------
    ValueError
        If *mode* is not one of ``"off"``, ``"soft"``, ``"hard"``.

    Examples
    --------
    >>> adj = load_adjacency(tissue="retina")
    >>> result = filter_cci_by_adjacency(lr_res, adj, mode="soft")
    >>> filtered = filter_cci_by_adjacency(lr_res, adj, mode="hard")
    """
    logger = log or _log

    # ── Mode validation ──────────────────────────────────────────────
    if mode not in {"off", "soft", "hard"}:
        raise ValueError(f"Unknown mode: '{mode}'. Expected one of: 'off', 'soft', 'hard'.")

    if mode == "off":
        return lr_res

    # ── Validate required columns ────────────────────────────────────
    for df, name in [(lr_res, "lr_res"), (adjacency, "adjacency")]:
        for col in ("source", "target"):
            if col not in df.columns:
                logger.warning(
                    "Column '%s' missing from %s. Returning lr_res%s.",
                    col,
                    name,
                    " unchanged" if mode == "hard" else " copy",
                )
                if mode == "hard":
                    return lr_res
                result = lr_res.copy()
                result["adjacent"] = False
                result["adjacency_type"] = ""
                return result

    total: int = len(lr_res)

    # ── Filter adjacency by type (if requested) ──────────────────────
    adj: pd.DataFrame = adjacency.copy()
    if adjacency_types:
        adj = adj[adj["adjacency_type"].isin(adjacency_types)]

    if adj.empty:
        logger.info("Adjacency is empty after filtering. No pairs matched.")
        if mode == "soft":
            result = lr_res.copy()
            result["adjacent"] = False
            result["adjacency_type"] = ""
            return result
        # hard mode
        result = lr_res.iloc[0:0].copy()
        logger.info(
            "Adjacency filter: %d/%d pairs retained (0.0%%)",
            0,
            total,
        )
        return result

    # ── Build case-insensitive lookup map ────────────────────────────
    # Map: (lower_source, lower_target) -> adjacency_type
    adj_map: dict[tuple[str, str], str] = {}
    for _, row in adj.iterrows():
        key = (
            str(row["source"]).strip().lower(),
            str(row["target"]).strip().lower(),
        )
        if key not in adj_map:  # first match wins for duplicates
            adj_map[key] = str(row.get("adjacency_type", ""))

    # ── Apply matching ───────────────────────────────────────────────
    if mode == "soft":
        result = lr_res.copy()
        adjacent_flags: list[bool] = []
        type_values: list[str] = []

        for _, row in lr_res.iterrows():
            source = str(row["source"]).strip().lower()
            target = str(row["target"]).strip().lower()
            is_adj, adj_type = _match_bidirectional(source, target, adj_map)
            adjacent_flags.append(is_adj)
            type_values.append(adj_type)

        result["adjacent"] = adjacent_flags
        result["adjacency_type"] = type_values

        n_matched = sum(adjacent_flags)
        logger.info(
            "Adjacency annotation: %d/%d pairs adjacent",
            n_matched,
            total,
        )
        return result

    # mode == "hard"
    match_mask: list[bool] = []
    for _, row in lr_res.iterrows():
        source = str(row["source"]).strip().lower()
        target = str(row["target"]).strip().lower()
        is_adj, _ = _match_bidirectional(source, target, adj_map)
        match_mask.append(is_adj)

    result = lr_res.loc[match_mask].copy()
    n_kept = len(result)
    pct = (n_kept / total * 100.0) if total > 0 else 0.0
    logger.info(
        "Adjacency filter: %d/%d pairs retained (%.1f%%)",
        n_kept,
        total,
        pct,
    )
    return result


# ======================================================================
#  Internal helpers
# ======================================================================


def _match_bidirectional(
    source: str,
    target: str,
    adj_map: dict[tuple[str, str], str],
) -> tuple[bool, str]:
    """Check if a (source, target) pair matches any adjacency row.

    Matching is **bidirectional**: A → B and B → A are both considered.

    Parameters
    ----------
    source : str
        Lowercased source cell type.
    target : str
        Lowercased target cell type.
    adj_map : dict
        Mapping ``(lower_source, lower_target) → adjacency_type``.

    Returns
    -------
    tuple[bool, str]
        ``(is_adjacent, adjacency_type)``.  Returns ``("",)`` as the type
        when no match is found.
    """
    key = (source, target)
    if key in adj_map:
        return True, adj_map[key]

    rev_key = (target, source)
    if rev_key in adj_map:
        return True, adj_map[rev_key]

    return False, ""
