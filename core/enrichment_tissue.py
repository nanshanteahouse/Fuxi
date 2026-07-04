#!/usr/bin/env python3
"""enrichment_tissue.py — Tissue-aware enrichment post-processor for ORA/GSEA results

``core/enrichment_tissue.py`` is the cross-modality entry point for tissue-aware
enrichment processing.  RNA, Spatial, and future ATAC ORA/GSEA steps use
these three functions to compute knowledge-base relevance, cluster redundant
pathways, and filter/annotate enrichment results.

Functions
---------
compute_pathway_relevance
    Compute KB-overlap statistics and a combined relevance score for each
    pathway in an enrichment result table.
cluster_redundant_pathways
    Greedy cluster of redundant pathways by Jaccard gene-set similarity.
filter_enrichment_by_tissue
    Annotate or filter enrichment results by tissue relevance (whitelist,
    blacklist, and KB relevance score).
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd


# ── Module-level logger ────────────────────────────────────────────────
_log = logging.getLogger(__name__)


# ======================================================================
#  Public API
# ======================================================================


def compute_pathway_relevance(
    pathway_df: pd.DataFrame,
    kb_markers: Optional[set[str]] = None,
    log: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Annotate each pathway with its overlap against a set of KB marker genes.

    For every row in *pathway_df*, the ``Overlap`` column (comma-separated
    gene symbols) is intersected with *kb_markers*.  Four new columns are
    added:

    - ``kb_overlap_genes`` (list[str])  — genes found in both the pathway and KB
    - ``kb_overlap_count`` (int)        — size of the overlap
    - ``kb_overlap_ratio`` (float)      — overlap / total pathway genes
    - ``kb_relevance_score`` (float)    — combined significance + prior score

    When *kb_markers* is ``None`` or empty, ``kb_relevance_score`` degrades
    to ``-log10(Adjusted P-value + 1e-300)``.

    Parameters
    ----------
    pathway_df : pd.DataFrame
        Enrichr / GSEApy output.  Expected columns: ``Term``,
        ``Adjusted P-value`` (or ``P-value`` as fallback), and ``Overlap``
        (comma-separated gene list string).
    kb_markers : set of str or None
        Knowledge-base marker gene symbols (e.g. ``"RHO"``, ``"NRL"``).
        When ``None`` or empty the relevance score falls back to statistical
        significance only.
    log : logging.Logger or None
        Optional logger.  If ``None`` the module-level logger is used.

    Returns
    -------
    pd.DataFrame
        A **copy** of *pathway_df* with the four new columns added.

    Examples
    --------
    >>> markers = {"RHO", "NRL", "GNAT1", "VSX2"}
    >>> enriched = compute_pathway_relevance(ora_result, markers)
    """
    logger = log or _log
    df = pathway_df.copy()

    # ── Resolve the P-value column ─────────────────────────────────────
    pval_col: str = "Adjusted P-value"
    if pval_col not in df.columns:
        if "P-value" in df.columns:
            pval_col = "P-value"
            logger.debug("Column 'Adjusted P-value' not found; falling back to 'P-value'.")
        else:
            logger.warning(
                "Neither 'Adjusted P-value' nor 'P-value' column found. "
                "kb_relevance_score will be 0.0 for all rows."
            )
            pval_col = ""

    # ── Determine whether we have Overlap data ─────────────────────────
    has_overlap: bool = "Overlap" in df.columns
    if not has_overlap:
        logger.warning(
            "Column 'Overlap' not found. All kb_overlap columns will be 0."
        )

    # ── Normalise kb_markers ───────────────────────────────────────────
    markers: set[str] | None = kb_markers
    if markers is not None and len(markers) == 0:
        markers = None  # treat empty set as "not provided"

    # ── Compute per-row statistics ─────────────────────────────────────
    kb_overlap_genes: list[list[str]] = []
    kb_overlap_count: list[int] = []
    kb_overlap_ratio: list[float] = []

    for _, row in df.iterrows():
        if not has_overlap:
            kb_overlap_genes.append([])
            kb_overlap_count.append(0)
            kb_overlap_ratio.append(0.0)
            continue

        overlap_str = str(row.get("Overlap", ""))
        if not overlap_str:
            kb_overlap_genes.append([])
            kb_overlap_count.append(0)
            kb_overlap_ratio.append(0.0)
            continue

        genes = [g.strip() for g in overlap_str.split(",") if g.strip()]
        total_genes = len(genes)

        if markers and total_genes > 0:
            overlap_genes = [g for g in genes if g in markers]
            cnt = len(overlap_genes)
            ratio = cnt / total_genes
            kb_overlap_genes.append(overlap_genes)
            kb_overlap_count.append(cnt)
            kb_overlap_ratio.append(ratio)
        else:
            kb_overlap_genes.append([])
            kb_overlap_count.append(0)
            kb_overlap_ratio.append(0.0)

    # ── Compute combined relevance score ───────────────────────────────
    kb_relevance_score: list[float] = []
    for i in range(len(df)):
        if pval_col and pval_col in df.columns:
            raw_p = float(df.iloc[i].get(pval_col, 1.0))
            # Clamp non-positive values
            if raw_p <= 0.0:
                raw_p = 1e-300
            log_p = -np.log10(raw_p + 1e-300)
        else:
            log_p = 0.0

        if markers:
            score = log_p * (1.0 + kb_overlap_ratio[i])
        else:
            score = log_p  # degradation to pure statistical
        kb_relevance_score.append(score)

    # ── Assign output columns ──────────────────────────────────────────
    df = df.copy()  # ensure we own a fresh copy
    df["kb_overlap_genes"] = kb_overlap_genes
    df["kb_overlap_count"] = kb_overlap_count
    df["kb_overlap_ratio"] = kb_overlap_ratio
    df["kb_relevance_score"] = kb_relevance_score

    n_marked = sum(1 for c in kb_overlap_count if c > 0)
    logger.debug(
        "compute_pathway_relevance: %d/%d rows have KB overlap (markers=%s).",
        n_marked,
        len(df),
        "provided" if markers else "None/empty",
    )

    return df


def cluster_redundant_pathways(
    pathway_df: pd.DataFrame,
    _term_col: str = "Term",
    gene_col: str = "Overlap",
    similarity_threshold: float = 0.6,
) -> pd.DataFrame:
    """Greedy cluster of redundant pathways based on Jaccard gene-set similarity.

    Rows are sorted by ``Adjusted P-value`` ascending (most significant first).
    The first row starts a new cluster and becomes its **representative**.
    Each subsequent row is compared against all existing representatives; if
    any pairwise Jaccard similarity >= *similarity_threshold*, it joins that
    cluster.  Otherwise a new cluster is started.

    Two new columns are added:

    - ``redundant_cluster_id`` (int)  — 0-indexed cluster identifier
    - ``redundant_representative`` (bool) — ``True`` for the most-significant
      member of each cluster

    Parameters
    ----------
    pathway_df : pd.DataFrame
        Enrichment results.  Must contain columns *term_col*, *gene_col*,
        and ``Adjusted P-value``.
    term_col : str
        Column name for the pathway / term description (default ``"Term"``).
    gene_col : str
        Column name for the comma-separated gene list (default ``"Overlap"``).
        Each cell is split by ``","`` and stripped to form a set.
    similarity_threshold : float
        Minimum Jaccard similarity to consider two pathways redundant.
        Must be in ``[0.0, 1.0]`` (default ``0.6``).

    Returns
    -------
    pd.DataFrame
        A **copy** of *pathway_df* with the two new columns added.

    Raises
    ------
    ValueError
        If *similarity_threshold* is outside ``[0.0, 1.0]``.

    Examples
    --------
    >>> clustered = cluster_redundant_pathways(ora_result)
    """
    # ── Input validation ───────────────────────────────────────────────
    if not (0.0 <= similarity_threshold <= 1.0):
        raise ValueError(
            f"similarity_threshold must be in [0.0, 1.0], got {similarity_threshold}."
        )

    n: int = len(pathway_df)
    if n == 0:
        return pathway_df.copy()

    df = pathway_df.copy()

    # ── Edge case: single row ──────────────────────────────────────────
    if n == 1:
        df["redundant_cluster_id"] = 0
        df["redundant_representative"] = True
        return df

    # ── Sort by Adjusted P-value ascending (most significant first) ────
    sort_col = "Adjusted P-value"
    if sort_col not in df.columns:
        # Fall back to index order if the column is missing
        sort_col = None

    if sort_col:
        df = df.sort_values(by=sort_col, ascending=True).reset_index(drop=True)
    else:
        df = df.reset_index(drop=True)

    # ── Parse gene sets ────────────────────────────────────────────────
    gene_sets: list[frozenset[str]] = []
    for _, row in df.iterrows():
        raw = str(row.get(gene_col, ""))
        if raw:
            genes = frozenset(g.strip() for g in raw.split(",") if g.strip())
        else:
            genes = frozenset()
        gene_sets.append(genes)

    # ── Jaccard similarity ─────────────────────────────────────────────
    def _jaccard(a: frozenset[str], b: frozenset[str]) -> float:
        union = len(a | b)
        if union == 0:
            return 0.0
        return len(a & b) / union

    # ── Greedy clustering ──────────────────────────────────────────────
    cluster_ids: list[int] = [-1] * n
    representatives: list[tuple[int, frozenset[str]]] = []  # (idx, gene_set)

    for i in range(n):
        assigned = False
        for rep_idx, rep_genes in representatives:
            sim = _jaccard(gene_sets[i], rep_genes)
            if sim >= similarity_threshold:
                cluster_ids[i] = cluster_ids[rep_idx]
                assigned = True
                break

        if not assigned:
            cid = len(representatives)
            cluster_ids[i] = cid
            representatives.append((i, gene_sets[i]))

    # ── Determine representatives ──────────────────────────────────────
    is_rep: list[bool] = [False] * n
    for rep_idx, _ in representatives:
        is_rep[rep_idx] = True

    df["redundant_cluster_id"] = cluster_ids
    df["redundant_representative"] = is_rep

    return df


def filter_enrichment_by_tissue(
    result_df: pd.DataFrame,
    mode: str,
    pathway_whitelist: Optional[list[str]] = None,
    pathway_blacklist: Optional[list[str]] = None,
    log: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Annotate or filter enrichment results by tissue relevance.

    Multi-level decision priority for each row:

    1. If the term (case-insensitive partial match) appears in
       *pathway_whitelist* → ``tissue_relevant=True``, score = 1.0.
    2. Else if the term appears in *pathway_blacklist* →
       ``tissue_relevant=False``, score = -1.0.
    3. Else if the ``kb_relevance_score`` column exists **and** its value is
       greater than the median of all non-whitelist/non-blacklist rows →
       ``tissue_relevant=True``, score = normalised KB score.
    4. Else → ``tissue_relevant=False``, score = 0.0.

    Parameters
    ----------
    result_df : pd.DataFrame
        Enrichment result table.  Must contain a column ``Term``.
        May optionally contain ``kb_relevance_score``.
    mode : {"off", "soft", "hard"}
        - ``"off"`` — return *result_df* unchanged (no copy, no columns added).
        - ``"soft"`` — add ``tissue_relevant`` (bool) and
          ``tissue_relevance_score`` (float) columns; keep all rows.
        - ``"hard"`` — add the same columns; return **only** rows where
          ``tissue_relevant == True``.
    pathway_whitelist : list of str or None
        Pathway terms that are always marked relevant, regardless of score.
    pathway_blacklist : list of str or None
        Pathway terms that are always marked irrelevant.
    log : logging.Logger or None
        Optional logger.  If ``None`` the module-level logger is used.

    Returns
    -------
    pd.DataFrame
        - ``"off"`` mode: the original *result_df* object (no copy).
        - ``"soft"`` mode: a copy with ``tissue_relevant`` and
          ``tissue_relevance_score`` added.
        - ``"hard"`` mode: a filtered copy.

    Raises
    ------
    ValueError
        If *mode* is not one of ``"off"``, ``"soft"``, ``"hard"``.

    Examples
    --------
    >>> filtered = filter_enrichment_by_tissue(
    ...     result_df, mode="soft",
    ...     pathway_whitelist=["Phototransduction"],
    ...     pathway_blacklist=["Ribosome", "Oxidative Phosphorylation"],
    ... )
    """
    logger = log or _log

    # ── Mode validation ────────────────────────────────────────────────
    if mode not in {"off", "soft", "hard"}:
        raise ValueError(
            f"Unknown mode: '{mode}'. Expected one of: 'off', 'soft', 'hard'."
        )

    if mode == "off":
        return result_df

    # ── Prepare lookup sets ────────────────────────────────────────────
    whitelist = {t.lower() for t in (pathway_whitelist or [])}
    blacklist = {t.lower() for t in (pathway_blacklist or [])}
    has_kb_score = "kb_relevance_score" in result_df.columns

    # ── Compute median KB score across non-special rows ────────────────
    # Step 1: compute the initial decision (priority 1 & 2) to separate
    #         the "unbiased" rows whose KB scores define the threshold.
    n_total: int = len(result_df)
    terms_lower = [str(t).lower() for t in result_df["Term"]]

    # Phase-A decisions (priority 1 & 2)
    phase_a_relevant: list[bool] = []
    phase_a_unbiased: list[bool] = []  # rows not caught by whitelist/blacklist

    for t_lower in terms_lower:
        in_whitelist = any(w in t_lower for w in whitelist) if whitelist else False
        in_blacklist = any(b in t_lower for b in blacklist) if blacklist else False

        if in_whitelist:
            phase_a_relevant.append(True)
            phase_a_unbiased.append(False)
        elif in_blacklist:
            phase_a_relevant.append(False)
            phase_a_unbiased.append(False)
        else:
            phase_a_relevant.append(False)  # undecided yet
            phase_a_unbiased.append(True)

    # Determine median KB score from unbiased rows
    median_kb: float = 0.0
    if has_kb_score and any(phase_a_unbiased):
        unbiased_scores = [
            float(result_df.iloc[i]["kb_relevance_score"])
            for i in range(n_total)
            if phase_a_unbiased[i]
        ]
        if unbiased_scores:
            median_kb = float(np.median(unbiased_scores))
            logger.debug(
                "Median kb_relevance_score (unbiased rows): %.4f", median_kb
            )

    # ── Final decisions ────────────────────────────────────────────────
    tissue_relevant: list[bool] = []
    tissue_relevance_score: list[float] = []

    for i in range(n_total):
        # Priority 1: whitelist
        if any(w in terms_lower[i] for w in whitelist) if whitelist else False:
            tissue_relevant.append(True)
            tissue_relevance_score.append(1.0)
            continue

        # Priority 2: blacklist
        if any(b in terms_lower[i] for b in blacklist) if blacklist else False:
            tissue_relevant.append(False)
            tissue_relevance_score.append(-1.0)
            continue

        # Priority 3: KB relevance score above median
        if has_kb_score:
            kb_val = float(result_df.iloc[i].get("kb_relevance_score", 0.0))
            if kb_val > median_kb:
                # Normalize KB score to [0, 1] via min-max of all KB scores
                all_kb = result_df["kb_relevance_score"].astype(float)
                k_min, k_max = float(all_kb.min()), float(all_kb.max())
                if k_max > k_min:
                    norm_score = (kb_val - k_min) / (k_max - k_min)
                else:
                    norm_score = 0.5
                tissue_relevant.append(True)
                tissue_relevance_score.append(norm_score)
                continue

        # Priority 4: remainder
        tissue_relevant.append(False)
        tissue_relevance_score.append(0.0)

    # ── Build result ───────────────────────────────────────────────────
    if mode == "soft":
        result = result_df.copy()
        result["tissue_relevant"] = tissue_relevant
        result["tissue_relevance_score"] = tissue_relevance_score

        n_relevant = sum(tissue_relevant)
        logger.info(
            "Soft filter: %d/%d pathways marked tissue-relevant",
            n_relevant,
            n_total,
        )
        return result

    # mode == "hard"
    mask = tissue_relevant
    result = result_df.loc[mask].copy()
    result["tissue_relevant"] = [True] * len(result)
    result["tissue_relevance_score"] = [s for i, s in enumerate(tissue_relevance_score) if mask[i]]

    n_kept = len(result)
    pct = (n_kept / n_total * 100.0) if n_total > 0 else 0.0
    logger.info(
        "Hard filter: kept %d/%d pathways (%d%%)",
        n_kept,
        n_total,
        round(pct),
    )
    return result
