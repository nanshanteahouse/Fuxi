"""
Unified annotation engine for cell type annotation.

Provides ``run_unified_annotation`` which performs KB-based unified annotation
with marker scoring, expert rules, evidence fusion, and optional AI fallback.

Extracted from ``rna/steps/05_annotate_major.py`` for cross-module reuse
(RNA, spatial, ATAC pipelines).
"""

import json
import os
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from core.utils import safe_plot

_FUXI_FAST_RANKS_THRESHOLD = 200_000


def _patch_scanpy_fast_ranks() -> None:
    """Speed up ``sc.tl.rank_genes_groups`` (wilcoxon) on large sparse matrices.

    scanpy's ``_ranks`` slices gene-column chunks out of a CSR matrix, which is
    O(nnz) per chunk (scipy CSR column slices scan every row).  For big datasets
    (n_cells * n_genes > 1e10) this dominates runtime on a single core.  This
    patch replaces ``_ranks`` with an equivalent implementation that transposes
    the matrix once (CSR row slices are O(chunk nnz)) and then ranks chunks with
    scanpy's own multi-threaded ``rankdata``.  Output is bit-identical to the
    original (verified end-to-end on 208k cells x 31.7k genes).
    """
    import scanpy.tools._rank_genes_groups as rgg

    if getattr(rgg, "_fuxi_fast_ranks", False):
        return
    orig = rgg._ranks
    cache: dict = {}

    def fast_ranks(x, mask_obs=None, mask_obs_rest=None):
        if (
            mask_obs is not None
            or not sparse.issparse(x)
            or x.shape[0] < _FUXI_FAST_RANKS_THRESHOLD
        ):
            yield from orig(x, mask_obs, mask_obs_rest)
            return
        n_cells, n_genes = x.shape
        xt = cache.get("xt")
        if xt is None:
            xt = x.T.tocsr()
            cache["xt"] = xt
        # Large chunks keep per-chunk column count high so scanpy's
        # multi-threaded rankdata (one thread per column) saturates cores.
        # rankdata is column-independent, so chunk size never changes output.
        genes_per_mb = (2**20) // max(n_cells * 4, 1)  # float32 dense budget
        chunk = max(128, min(1024, genes_per_mb))
        for left in range(0, n_genes, chunk):
            right = min(left + chunk, n_genes)
            yield rgg.rankdata(xt[left:right].toarray().T), left, right

    rgg._ranks = fast_ranks
    rgg._fuxi_fast_ranks = True


def _patch_scanpy_wilcoxon() -> None:
    """Replace scanpy's wilcoxon rank-sum accumulation with a cumsum sweep.

    scanpy computes per-cluster rank sums as ``ranks[mask].sum(axis=0)`` for
    every cluster, re-walking each n_cells x chunk_genes ranks chunk once per
    cluster (70x on 70-cluster datasets, ~186 TB of memory traffic on the
    1.05M-cell retina atlas).  Here cells are sorted by cluster once, then each
    chunk is reduced with a single cumulative sum; cluster sums are the
    cumsum differences at cluster boundaries.  Ranks are exact integers, so
    results are bit-identical to the original.
    """
    import scanpy.tools._rank_genes_groups as rgg

    if getattr(rgg, "_fuxi_fast_wilcoxon", False):
        return
    orig_wilcoxon = rgg._RankGenes.wilcoxon

    def fast_wilcoxon(self, *, tie_correct):
        from scipy import stats

        self._basic_stats()
        n_genes = self.X.shape[1]
        if self.ireference is not None:
            yield from orig_wilcoxon(self, tie_correct=tie_correct)
            return
        n_groups = self.groups_masks_obs.shape[0]
        scores = np.zeros((n_groups, n_genes))
        n_cells = self.X.shape[0]
        if tie_correct:
            tc_coef = np.zeros((n_groups, n_genes))
        # One-time cluster ordering (cells are exactly one-hot in masks).
        labels = np.argmax(self.groups_masks_obs, axis=0)
        order = np.argsort(labels, kind="stable")
        bounds = np.searchsorted(labels[order], np.arange(n_groups + 1))
        starts = bounds[:-1]
        ends = bounds[1:]
        zero_row: np.ndarray | None = None
        for ranks, left, right in rgg._ranks(self.X):
            if tie_correct:
                tc_coef[:, left:right] = rgg._tiecorrect(ranks)
            rs = ranks[order]
            np.cumsum(rs, axis=0, out=rs)
            end_vals = rs[ends - 1, :]
            c = right - left
            if zero_row is None or zero_row.shape[1] != c:
                zero_row = np.zeros((1, c))
            start_vals = np.concatenate([zero_row, rs[:-1]])[starts]
            scores[:, left:right] = end_vals - start_vals
        for group_index, mask_obs in enumerate(self.groups_masks_obs):
            n_active = np.count_nonzero(mask_obs)
            coef = tc_coef[group_index] if tie_correct else 1
            std_dev = np.sqrt(coef * n_active * (n_cells - n_active) * (n_cells + 1) / 12.0)
            scores[group_index, :] = (
                scores[group_index, :] - (n_active * (n_cells + 1) / 2.0)
            ) / std_dev
            scores[np.isnan(scores)] = 0
            pvals = 2 * stats.distributions.norm.sf(np.abs(scores[group_index, :]))
            yield group_index, scores[group_index], pvals

    rgg._RankGenes.wilcoxon = fast_wilcoxon
    rgg._fuxi_fast_wilcoxon = True


def _ribo_fallback_pct_scores(adata, kb, cl_str, logger):
    """Re-score a ribo-dominated cluster via raw expression fractions (pct).

    Fisher scores saturate to 1.0 once ribosomal genes are filtered, so the
    hypergeometric rank cannot separate cell types.  This helper instead
    scores every KB cell type by the mean fraction of its positive markers
    expressed in the cluster (same discriminant used in external portability evaluation),
    and returns ``(top_type, top_score, top1_top2_gap)`` when the top type
    is decisive (score >= 0.25 and gap > 0.1), else ``None``.
    """
    if adata.raw is None:
        return None
    raw_vars = list(adata.raw.var_names)
    ridx = {g: i for i, g in enumerate(raw_vars)}
    mask = (adata.obs["leiden"].astype(str) == str(cl_str)).values
    if not mask.any():
        return None
    x_mat = adata.raw.X
    sub = x_mat[mask]

    best: list[tuple[str, float]] = []
    for tkey, tdata in kb.items():
        if tkey == "expert_rules" or tkey.startswith(("_", "Broad_")):
            continue
        markers = tdata.get("markers", {}) if isinstance(tdata, dict) else {}
        genes = []
        for tier in ("confirm", "add"):
            tm = markers.get(tier, {})
            genes += list(tm.keys()) if isinstance(tm, dict) else (tm or [])
        genes = [g for g in genes if g in ridx]
        if not genes:
            continue
        pcts = (sub[:, [ridx[g] for g in genes]] > 0).mean(axis=0)
        best.append((tkey, float(pcts.mean())))

    if len(best) < 2:
        return None
    best.sort(key=lambda kv: kv[1], reverse=True)
    top1, top2 = best[0], best[1]
    gap = top1[1] - top2[1]
    if top1[1] >= 0.25 and gap > 0.1:
        return (top1[0], top1[1], gap)
    return None


# ═══════════════════════════════════════════════════════════════════════
#  D3 — canonical-expression fallback + ai_only audit (plan §3 D3)
# ═══════════════════════════════════════════════════════════════════════
# ``consensus_levels`` maps a marker gene to its multi-source support level:
# "low" (1 source) | "medium" (2) | "high" (3-4) | "gold" (5+).
_CONSENSUS_LEVEL_RANK = {"low": 0, "medium": 2, "high": 3, "gold": 4}
_CANONICAL_TOP_N = 3


def _top_consensus_markers(kb, type_key, top_n=_CANONICAL_TOP_N):
    """Top-*top_n* confirm markers of *type_key* backed by >= 2 KB sources.

    Only markers whose ``consensus_levels`` entry ranks at least ``"medium"``
    (>= 2 independent sources) qualify, ordered gold > high > medium.  A
    single-source type (e.g. RGC_Alpha, only tran2019) has no qualifying
    marker, so the returned list is empty — the D3 canonical-expression
    check then degrades into a natural no-op.
    """
    entry = kb.get(type_key) or {}
    levels = entry.get("consensus_levels") or {}
    confirm = (entry.get("markers") or {}).get("confirm") or {}
    ranked = sorted(
        (g for g in confirm if _CONSENSUS_LEVEL_RANK.get(levels.get(g, "low"), 0) >= 2),
        key=lambda g: _CONSENSUS_LEVEL_RANK.get(levels.get(g, "low"), 0),
        reverse=True,
    )
    return ranked[:top_n]


def _cluster_marker_pcts(adata, cl_str, genes):
    """Fraction of cluster cells with raw counts > 0 for each *gene*.

    Reuses ``_ribo_fallback_pct_scores``'s raw-matrix access pattern — CSR
    row slice ``x_mat[mask]`` followed by column-index ``> 0`` fractions
    (no ``toarray()`` materialisation).

    Returns ``None`` when ``adata.raw`` is missing, the cluster is empty,
    or none of *genes* exist in the raw var index (the D3 canonical check
    then silently skips, mirroring ``_ribo_fallback_pct_scores``).
    """
    if adata.raw is None:
        return None
    raw_vars = list(adata.raw.var_names)
    ridx = {g: i for i, g in enumerate(raw_vars)}
    mask = (adata.obs["leiden"].astype(str) == str(cl_str)).values
    if not mask.any():
        return None
    present = [g for g in genes if g in ridx]
    if not present:
        return None
    x_mat = adata.raw.X
    sub = x_mat[mask]
    pcts = (sub[:, [ridx[g] for g in present]] > 0).mean(axis=0)
    return np.asarray(pcts, dtype=float).reshape(-1)


def _apply_canonical_expression_fallback(adata, kb, decision_map, CFG, logger):  # noqa: N803
    """D3 canonical-expression fallback for confident marker-scoring decisions.

    For ``marker_scoring``-family decisions with ``confidence in ("high",
    "medium")``, fetch the winning type's top-consensus confirm markers
    (>= 2 KB sources, top 3) and compute each marker's pct expression in the
    cluster (raw counts > 0).  When **all** of them fall below
    ``CFG.annotation.canonical_pct_floor`` (default 0.05), the confident
    label is unsupported by raw expression → downgrade ``confidence="low"``,
    record ``review_reason="no_canonical_expression"`` on the decision, and
    return ``{cluster_str: reason}`` as the review ledger.

    Silent skips (documented semantics):
    - ``adata.raw is None`` → expression statistics are undefined; the whole
      check is skipped (mirrors ``_ribo_fallback_pct_scores``).
    - winning type with no consensus>=2 markers (single-source types such as
      RGC_Alpha, only tran2019) → empty top-N set, naturally skipped.
    """
    if adata.raw is None:
        return {}
    floor = getattr(getattr(CFG, "annotation", None), "canonical_pct_floor", 0.05)
    reasons: dict[str, str] = {}
    for cl_str, d in decision_map.items():
        if d.confidence not in ("high", "medium"):
            continue
        if not d.method.startswith("marker_scoring"):
            continue
        markers = _top_consensus_markers(kb, d.cell_type)
        if not markers:
            continue
        pcts = _cluster_marker_pcts(adata, cl_str, markers)
        if pcts is None:
            continue
        if all(p < floor for p in pcts):
            logger.warning(
                "Cluster %s: all top-consensus markers of %s below pct floor %.3f"
                " (%s) — downgrading %s → low",
                cl_str,
                d.cell_type,
                floor,
                ", ".join(f"{g}={p:.3f}" for g, p in zip(markers, pcts)),
                d.confidence,
            )
            decision_map[cl_str] = d._replace(
                confidence="low",
                review_reason="no_canonical_expression",
                explanation=(
                    f"{d.explanation} | no_canonical_expression: top-consensus "
                    f"markers ({', '.join(markers)}) all below {floor:.2f} pct "
                    f"in cluster"
                ),
            )
            reasons[cl_str] = "no_canonical_expression"
    return reasons


def _flag_ai_only_decisions(decision_map, logger):
    """Flag AI-only decisions (D3) and force them down to low confidence.

    A decision is ``ai_only`` when AI was involved (``method == "ai_unconstrained"``
    or ``ai_agreed``/``ai_suggested``) and there is no meaningful KB marker
    support (``score < 0.25`` or ``n_markers_found == 0``).  Such decisions
    get ``review_reason="ai_only"`` and — when currently ``high``/``medium``
    — are forced down to ``confidence="low"``.

    Returns ``{cluster_str: "ai_only"}`` (the review ledger).
    """
    reasons: dict[str, str] = {}
    for cl_str, d in decision_map.items():
        ai_involved = d.method == "ai_unconstrained" or d.ai_agreed or bool(d.ai_suggested)
        if not ai_involved:
            continue
        if d.score >= 0.25 and d.n_markers_found > 0:
            continue  # KB marker support present — not ai_only
        if d.confidence in ("high", "medium"):
            logger.info(
                "Cluster %s: ai_only (no KB marker support) — forcing %s → low",
                cl_str,
                d.confidence,
            )
            decision_map[cl_str] = d._replace(
                confidence="low",
                review_reason="ai_only",
                explanation=f"{d.explanation} | ai_only: no KB marker support",
            )
        else:
            decision_map[cl_str] = d._replace(review_reason="ai_only")
        reasons[cl_str] = "ai_only"
    return reasons


def _check_zero_scores_and_retry(
    kb,
    all_scores,
    marker_df,
    clusters,
    species,
    target_class,
    target_order,
    tissue_kb,
    CFG,  # noqa: N803
    logger,
):
    """Check for zero KB marker scores and attempt case-insensitive retry.

    Builds an uppercased deepcopy of the KB once and reuses it across two
    retry passes (the original KB is never mutated):

    1. Global retry — fires when *all* clusters have zero hits and there
       are >= 5 clusters.  Re-scores every cluster; replaces ``all_scores``
       if any cluster improves.

    2. Per-cluster retry — fires when any clusters remain zero-hit and
       there are >= 3 clusters.  Re-scores individual zero-hit clusters
       and updates them in place when their hit count improves.

    A diagnostic ERROR fires if zero hits persist across all clusters
    after retry.

    Parameters
    ----------
    kb : dict
        Original KB (never mutated).
    all_scores : dict
        ``{cluster_str: {type_key: Score}}`` from the initial scoring pass.
    marker_df : pd.DataFrame
        Concatenated marker DataFrames with ``names``, ``logfoldchanges``,
        and ``cluster`` columns.
    clusters : list
        Sorted cluster identifiers.
    species : str
        Species filter for KB scoring.
    target_class : str
        Taxonomic class for phylogenetic weighting.
    target_order : str
        Taxonomic order for phylogenetic weighting.
    tissue_kb : str
        KB identifier for diagnostic messages.
    logger : logging.Logger

    Returns
    -------
    tuple[dict, int, int]
        ``(all_scores, total_hits, n_clusters_total)`` — ``all_scores`` may
        be updated if the retry improved hits.
    """
    total_hits = sum(
        score.n_markers_found
        for cluster_scores in all_scores.values()
        for score in cluster_scores.values()
    )
    n_clusters_total = len(all_scores)

    # Detect retry eligibility up front so we can build the uppercased KB
    # copy at most once and reuse it across both retry paths.
    _expand_steps = CFG.marker.candidate_pool_expand_steps
    zero_hit_clusters = [
        cl_str
        for cl_str, cs in all_scores.items()
        if sum(s.n_markers_found for s in cs.values()) == 0
    ]
    needs_global_retry = total_hits == 0 and n_clusters_total >= 5
    needs_per_cluster_retry = bool(zero_hit_clusters) and n_clusters_total >= 3

    if needs_global_retry or needs_per_cluster_retry:
        import copy

        from core.annotation.scoring import score_cluster_against_kb

        # Build uppercased KB deepcopy once; reused by both retry passes.
        kb_copy = copy.deepcopy(kb)
        for type_key in kb_copy:
            if type_key == "expert_rules" or type_key.startswith("_"):
                continue
            markers = kb_copy[type_key].get("markers", {})
            for tier in ("confirm", "add"):
                old_genes = markers.get(tier, {})
                if old_genes:
                    markers[tier] = {g.upper(): v for g, v in old_genes.items()}

        def _retry_cluster_score(cl_str, expand_steps=None):
            """Re-score a single cluster against the uppercased KB copy."""
            cl = int(cl_str) if cl_str.lstrip("-").isdigit() else cl_str
            cl_mask = marker_df["cluster"] == cl
            cl_data = marker_df[cl_mask].copy()
            lfc_idx = cl_data["logfoldchanges"].argsort()[::-1]
            cl_data = cl_data.iloc[lfc_idx]
            cl_data["names"] = cl_data["names"].str.upper()
            return score_cluster_against_kb(
                kb_copy,
                cl_data,
                species=species,
                target_class=target_class,
                target_order=target_order,
                adaptive_top_n=True,
                expand_steps=expand_steps,
            )

        # Global retry: all clusters zero-hit (n_clusters_total >= 5).
        # Re-score every cluster; replace all_scores if any improvement.
        if needs_global_retry:
            logger.warning(
                "Zero KB marker hits across all %d clusters — "
                "attempting case-insensitive normalization retry",
                n_clusters_total,
            )
            retry_scores = {}
            for cl in clusters:
                retry_scores[str(cl)] = _retry_cluster_score(str(cl), expand_steps=_expand_steps)
            retry_total_hits = sum(
                s.n_markers_found
                for cluster_scores in retry_scores.values()
                for s in cluster_scores.values()
            )
            if retry_total_hits > total_hits:
                logger.info(
                    "Case-insensitive retry improved hits from %d to %d — "
                    "using normalized results",
                    total_hits,
                    retry_total_hits,
                )
                all_scores = retry_scores
                total_hits = retry_total_hits
                # Recompute zero-hit list for the per-cluster pass below.
                zero_hit_clusters = [
                    cl_str
                    for cl_str, cs in all_scores.items()
                    if sum(s.n_markers_found for s in cs.values()) == 0
                ]

        # Per-cluster retry: any remaining zero-hit clusters (n_clusters_total >= 3).
        # Updates individual clusters only when retry improves their hit count.
        if zero_hit_clusters and n_clusters_total >= 3:
            logger.warning(
                "%d/%d clusters have zero KB hits — retrying per-cluster "
                "with case-insensitive normalization",
                len(zero_hit_clusters),
                n_clusters_total,
            )
            per_cluster_fixed = 0
            for cl_str in zero_hit_clusters:
                retry_result = _retry_cluster_score(cl_str, expand_steps=_expand_steps)
                retry_hits = sum(s.n_markers_found for s in retry_result.values())
                if retry_hits > 0:
                    all_scores[cl_str] = retry_result
                    per_cluster_fixed += 1
                    total_hits += retry_hits
            if per_cluster_fixed:
                logger.info(
                    "Per-cluster retry recovered %d/%d zero-hit clusters",
                    per_cluster_fixed,
                    len(zero_hit_clusters),
                )
    # Zero-score diagnostic warning (fires regardless of retry outcome)
    if total_hits == 0 and n_clusters_total >= 5:
        logger.error(
            "Zero KB marker hits across all %d clusters for KB '%s' "
            "(species=%s) — likely species mismatch or missing cell types. "
            "Consider expanding KB coverage or using AI fallback.",
            n_clusters_total,
            tissue_kb,
            species,
        )

    return all_scores, total_hits, n_clusters_total


CATEGORY_PREFIX = "Broad_"


def _classify_broad_category(
    cluster_scores: dict, kb: dict, allows_transitions: bool = False, transition_context: str = ""
) -> str:
    """Classify a cluster into a broad category using Fisher scores.

    Filters cluster_scores to only keys starting with CATEGORY_PREFIX
    (Broad_Progenitor, Broad_Neuron, Broad_Glia, Broad_Non-neural),
    then returns the category name with the highest Fisher score.

    Args:
        cluster_scores: Dict of {type_key: Score} from score_cluster_against_kb()
        kb: The full KB dict (unused; required for signature consistency)
    allows_transitions: bool
        When ``False``, skip Broad_Progenitor (adult tissue should never
        be classified into a progenitor category).
    transition_context: str
        Forward-looking context string for future tissue types (e.g., "tumor").

    Returns:
        Category name string (e.g., "Broad_Neuron") or "" if no match.
    """
    from rna.utils.evidence_fusion import _resolve_score

    broad_entries = {k: v for k, v in cluster_scores.items() if k.startswith(CATEGORY_PREFIX)}
    # Adult tissue: never assign Broad_Progenitor
    if not allows_transitions:
        broad_entries = {
            k: v for k, v in broad_entries.items() if k != f"{CATEGORY_PREFIX}Progenitor"
        }
    if not broad_entries:
        return ""

    best_type = max(broad_entries, key=lambda k: _resolve_score(cluster_scores, k)[0])
    return best_type


def _map_cell_state(decision, cell_category: str) -> str:
    """Map a FusionDecision to a ``cell_state`` (D5 six-class semantics).

    Rows match top-down; the first hit wins:

    ``method == "transition_state"`` → ``transient_transitional``
    ``method == "ambiguous"`` → ``N/A`` (downgraded; manual review)
    ``method == "developmental_potency"`` → ``differentiating`` (KADP;
    emitted before the Proliferating row so a KADP-named type that contains
    "Proliferating" still maps to ``differentiating``)
    ``"Proliferating" in cell_type`` → ``cycling``
    ``cell_category == "Broad_Progenitor"`` → ``committed_precursor``
    ``confidence in ("high", "medium")`` and ``cell_category in
    ("Broad_Neuron", "Broad_Glia", "Broad_Non-neural")`` → ``terminal``
    ``confidence == "low"`` → ``N/A``
    else (unknown etc.) → ``N/A``

    Here ``terminal`` = terminal node of the annotation hierarchy
    (non-precursor), NOT biological terminal differentiation — e.g. Müller
    Glia in developing tissue may retain plasticity, and its proliferating
    subset is captured by the ``cycling`` row ("Proliferating_MG" contains
    "Proliferating").  ``differentiating`` is emitted by the KADP potency
    path (method == "developmental_potency", plan annotation-kadp-metc).
    """
    if decision.method == "transition_state":
        return "transient_transitional"
    if decision.method == "ambiguous":
        return "N/A"
    # KADP (layer 3): differentiating precursor named by developmental
    # potency.  Row sits BEFORE the Proliferating row so a KADP type like
    # Proliferating_RPC still reads differentiating.
    if decision.method == "developmental_potency":
        return "differentiating"
    if "Proliferating" in decision.cell_type:
        return "cycling"
    if cell_category == "Broad_Progenitor":
        return "committed_precursor"
    if decision.confidence in ("high", "medium") and cell_category in (
        "Broad_Neuron",
        "Broad_Glia",
        "Broad_Non-neural",
    ):
        return "terminal"
    # D5 final rows: confidence == "low" and all remaining (unknown etc.) → N/A
    return "N/A"


def run_unified_annotation(adata, CFG, logger):  # noqa: N803
    """
    KB-based unified annotation mode.

    Uses marker scoring + expert rules + evidence fusion, with AI fallback
    for low-confidence clusters.

    Returns
    -------
    dict or None
        Cluster -> FusionDecision mapping, or None on failure (triggers fallback).
    """
    # ── a. Compute marker genes ───────────────────────────────────────────
    _n_genes = max(CFG.marker.candidate_pool_expand_steps)
    logger.info("Computing marker genes (Wilcoxon rank-sum, n_genes=%d)...", _n_genes)
    _patch_scanpy_fast_ranks()
    _patch_scanpy_wilcoxon()
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        n_genes=_n_genes,
        use_raw=True if adata.raw is not None else None,
    )
    n_clusters = adata.obs["leiden"].nunique()

    # ── b. Build per-cluster marker DataFrames ────────────────────────────
    marker_rows = []
    for cl in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        df = sc.get.rank_genes_groups_df(adata, group=str(cl))
        df["cluster"] = cl
        marker_rows.append(df)
    marker_df = pd.concat(marker_rows, ignore_index=True)
    marker_csv = os.path.join(CFG.table_dir, "marker_genes_unified.csv")
    marker_df.to_csv(marker_csv, index=False)
    logger.info("Marker genes saved: %s", marker_csv)

    # ── c. Load KB ────────────────────────────────────────────────────────
    from core.kb import load_kb

    try:
        kb = load_kb(CFG.tissue_kb)
    except Exception as exc:
        logger.warning("Failed to load KB '%s': %s", CFG.tissue_kb, exc)
        return None

    n_types = sum(
        1
        for k in kb
        if k != "expert_rules" and not k.startswith("_") and not k.startswith(CATEGORY_PREFIX)
    )
    n_rules = len(kb.get("expert_rules", []))
    logger.info("Loaded KB: %s (%d cell types, %d rules)", CFG.tissue_kb, n_types, n_rules)

    # ── d. Full marker scoring + expert rules per cluster ─────────────────
    from core.annotation.potency import KADPConfig
    from core.annotation.scoring import (
        detect_low_quality_cluster,
        score_cluster_against_kb,
    )
    from rna.utils.evidence_fusion import METCConfig, fuse_all_clusters
    from rna.utils.marker_expert_rules import (
        apply_expert_rules,
        resolve_expert_rule_params,
    )

    species = CFG.species
    # ── Normalise species key (e.g. "danio_rerio" → "zebrafish") ────
    # _SPECIES_NORMALISE maps NCBI-style names to pipeline keys; this
    # ensures consistency with _SPECIES_SYNONYMS and SPECIES_TO_CLASS.
    from core.preprocess.format_detector import _SPECIES_NORMALISE

    species = _SPECIES_NORMALISE.get(species, species)

    # ── Resolve taxonomic class/order for phylogenetic weighting ──────
    # CFG.target_class/order take precedence; fall back to species lookup.
    from rna.ortholog import SPECIES_TO_CLASS

    target_class = CFG.target_class or SPECIES_TO_CLASS.get(species, "")
    target_order = CFG.target_order or ""
    if target_class:
        logger.info(
            "Phylogenetic weighting: target_class=%s, target_order=%s",
            target_class,
            target_order or "(none)",
        )
    else:
        logger.info("Phylogenetic weighting: disabled (no target_class for '%s')", species)

    # ── Resolve expert-rule constraint parameters ────────────────────
    _strictness = getattr(CFG.marker, "expert_rule_strictness", "default")
    _top_n = getattr(CFG.marker, "expert_rule_top_n", 0)
    _pval = getattr(CFG.marker, "expert_rule_pval_cutoff", 0.0)

    # ── developmental_mode: auto-relax constraints for developing tissues ──
    _dev_mode = getattr(CFG.marker, "developmental_mode", False)
    if _dev_mode:
        if _strictness == "default":
            _strictness = "deep"
            logger.info("developmental_mode: auto-adjusted strictness from 'default' → 'deep'")
        logger.info("developmental_mode: enabled — relaxed constraints for developing tissue")

    _tissue_maturity = getattr(CFG, "tissue_maturity", "")
    # ── transition-aware gating ──────────────────────────────────────
    # Map tissue_maturity to a transition context string.
    # Future tissues (e.g. tumor) can add new contexts here.
    if _tissue_maturity == "developing":
        transition_context = "developmental"
    elif _tissue_maturity == "tumor":
        transition_context = "tumor"
    else:
        transition_context = ""
    allows_transitions = bool(transition_context) or _dev_mode

    # Extract incompatible transition pairs from tissue hierarchy
    _incompatible_transitions = kb.get("_hierarchy", {}).get("incompatible_transitions", [])

    rule_top_n, rule_pval = resolve_expert_rule_params(
        strictness=_strictness,
        top_n=_top_n,
        pval_cutoff=_pval,
    )

    # ── Compute per-cluster marker scores and expert rules ───────────
    all_scores = {}
    all_rules = {}
    all_top_genes: dict[str, list[str]] = {}
    low_quality_clusters: dict[str, str] = {}  # cluster_str → reason
    clusters = sorted(
        marker_df["cluster"].unique(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )
    for cl in clusters:
        cl_str = str(cl)
        cl_mask = marker_df["cluster"] == cl
        cl_data = marker_df[cl_mask].copy()
        lfc_idx = cl_data["logfoldchanges"].argsort()[::-1]
        cl_data = cl_data.iloc[lfc_idx]

        # Path C: detect low-quality clusters (mito/ribo dominated)
        is_lq, lq_reason = detect_low_quality_cluster(cl_data)
        if is_lq:
            low_quality_clusters[cl_str] = lq_reason
            logger.info("Cluster %s flagged as low-quality: %s", cl_str, lq_reason)

        all_scores[cl_str] = score_cluster_against_kb(
            kb,
            cl_data,
            species=species,
            target_class=target_class,
            target_order=target_order,
            adaptive_top_n=True,
            expand_steps=CFG.marker.candidate_pool_expand_steps,
        )

        # ribo_high fallback: scoring already filters RPL/RPS/MT- genes, so a
        # cluster flagged ribo-dominated may still carry a strong, unambiguous
        # cell-type signal once the ribosomal noise is removed.  Fisher scores
        # saturate to 1.0 when the filtered top markers overlap the KB heavily,
        # so when Fisher cannot separate top1/top2 we re-score with raw
        # expression fractions (pct of positive markers expressed), which is
        # the same discriminant used in external portability evaluation.  If the top
        # scored type is clearly separated from the runner-up, release the
        # low-quality hold so evidence fusion can annotate it (typically as a
        # low-confidence call) instead of forcing Unknown.  mito_high stays
        # conservative — mitochondrial dominance usually means dead cells.
        if is_lq and lq_reason.startswith("ribo_high"):
            _cl_scores = all_scores[cl_str]
            _ranked = sorted(_cl_scores.values(), key=lambda s: s.score, reverse=True)
            _fisher_separated = (
                len(_ranked) >= 2
                and _ranked[0].score >= 0.25
                and (_ranked[0].score - _ranked[1].score) > 0.1
            )
            _pct_pick: tuple[str, float, float] | None = None
            if not _fisher_separated:
                # Re-score via raw expression fractions (pct of markers).
                _pct_pick = _ribo_fallback_pct_scores(adata, kb, cl_str, logger)
                if _pct_pick is not None and _pct_pick[1] >= 0.25:
                    _fisher_separated = True
            if _fisher_separated:
                low_quality_clusters.pop(cl_str, None)
                if _pct_pick is not None:
                    logger.info(
                        "Cluster %s: %s but pct re-score decisive "
                        "(top1=%s %.3f, top1-top2=%.3f) — annotating",
                        cl_str,
                        lq_reason,
                        _pct_pick[0],
                        _pct_pick[1],
                        _pct_pick[2],
                    )
                else:
                    _top_key = max(_cl_scores, key=lambda k: _cl_scores[k].score)
                    logger.info(
                        "Cluster %s: %s but filtered Fisher decisive "
                        "(top1=%s %.3f, top1-top2=%.3f) — annotating",
                        cl_str,
                        lq_reason,
                        _top_key,
                        _ranked[0].score,
                        _ranked[0].score - _ranked[1].score,
                    )
        all_top_genes[cl_str] = cl_data["names"].head(20).tolist()
        all_rules[cl_str] = apply_expert_rules(
            kb, cl_data, top_n=rule_top_n, pval_cutoff=rule_pval
        )

    all_scores, total_hits, n_clusters_total = _check_zero_scores_and_retry(
        kb,
        all_scores,
        marker_df,
        clusters,
        species,
        target_class,
        target_order,
        CFG.tissue_kb,
        CFG,
        logger,
    )

    # Pre-filter: strip Broad_* synthetic types from scoring before fusion.
    # Fine-grained cell_type annotation must never see Broad_* entries.
    # The original all_scores (with Broad_* entries) is preserved for
    # _classify_broad_category() below.
    fine_scores = {}
    for cl_str, scores in all_scores.items():
        fine_scores[cl_str] = {
            k: v for k, v in scores.items() if not k.startswith(CATEGORY_PREFIX)
        }

    # ── CellTypist supplementary annotation ────────────────────────────
    celltypist_results: dict[str, str] = {}
    _celltypist_triggered = (
        getattr(CFG.annotation.celltypist, "enabled", False)
        or getattr(CFG.annotation, "method", "kb_unified") == "celltypist"
    )
    if _celltypist_triggered:
        _model_name = CFG.annotation.celltypist.model
        if not _model_name:
            logger.warning("CellTypist enabled but no model specified — skipping")
        else:
            try:
                from core.utils._optional import require_celltypist

                require_celltypist()
                import celltypist  # type: ignore[import-untyped]

                _model = celltypist.models.Model.load(model=_model_name)
                _mv = CFG.annotation.celltypist.majority_voting
                # celltypist >= 1.6 returns an AnnotationResult WITHOUT mutating
                # adata.obs: capture it and read per-cell labels from
                # ``_res.predicted_labels`` (a DataFrame whose column is
                # "majority_voting" when majority_voting=True else
                # "predicted_labels", index aligned with adata.obs_names;
                # reindex to obs_names defensively).
                _res = celltypist.annotate(adata, model=_model, majority_voting=_mv)
                _label_col = "majority_voting" if _mv else "predicted_labels"
                _labels = _res.predicted_labels
                if hasattr(_labels, "reindex"):
                    _labels = _labels.reindex(adata.obs_names)
                if _label_col in _labels.columns:
                    for cl in clusters:
                        cl_str = str(cl)
                        _mask = adata.obs["leiden"].astype(str) == cl_str
                        _types = _labels.loc[adata.obs_names[_mask], _label_col].mode()
                        if len(_types) > 0:
                            celltypist_results[cl_str] = str(_types[0])
                    logger.info(
                        "CellTypist: predicted %d/%d clusters via '%s'",
                        len(celltypist_results),
                        len(clusters),
                        _model_name,
                    )
                else:
                    logger.warning(
                        "CellTypist: '%s' column not found in AnnotationResult.predicted_labels — skipping",
                        _label_col,
                    )
            except Exception as exc:
                logger.warning("CellTypist prediction failed: %s — skipping", exc)

    # ── Layer-3 KADP developmental-potency config (plan annotation-kadp-metc)
    # Constructed ONCE from CFG.annotation and mirrored into BOTH
    # fuse_all_clusters calls (first pass + AI second pass) so a KADP naming
    # survives the AI-enhanced re-fusion (Oracle r1 BLOCKER 1).
    kadp_cfg = KADPConfig(
        enabled=getattr(CFG.annotation, "kadp_enabled", False),
        ratio_threshold=getattr(CFG.annotation, "kadp_ratio_threshold", 2.0),
        abs_threshold=getattr(CFG.annotation, "kadp_abs_threshold", 0.6),
        gap_threshold=getattr(CFG.annotation, "kadp_gap_threshold", 0.1),
        use_gap_criterion=getattr(CFG.annotation, "use_gap_criterion", False),
    )

    # ── Layer-4 METC multi-source voting config (plan annotation-kadp-metc) ──
    # Constructed ONCE from CFG.annotation and mirrored into BOTH
    # fuse_all_clusters calls (first pass + AI second pass) exactly like the
    # KADPConfig above (Oracle r1 BLOCKER 1 extended to METC by todo 10).
    metc_cfg = METCConfig(
        enabled=getattr(CFG.annotation, "metc_enabled", False),
        min_sources=getattr(CFG.annotation, "metc_min_sources", 3),
        min_distinct_transition=getattr(CFG.annotation, "metc_min_distinct_transition", 3),
    )

    # ── Label harmonization resources (todo 8) ────────────────────────
    # CellTypist and AI labels are resolved to canonical KB names through
    # the SHARED harmonize_label chain (parallel A/B evaluation, ambiguous
    # labels abstain).  Synonyms mirror into BOTH fuse_all_clusters calls.
    from core.kb import load_synonyms
    from rna.utils.evidence_fusion import harmonize_label

    synonyms = load_synonyms(CFG.tissue_kb)

    decisions, fusion_quality = fuse_all_clusters(
        fine_scores,
        all_rules,
        kb=kb,
        all_marker_dfs=marker_df,
        return_quality=True,
        low_quality_clusters=low_quality_clusters,
        unconstrained=getattr(CFG.ai, "unconstrained_annotation", False),
        allows_transitions=allows_transitions,
        incompatible_transitions=_incompatible_transitions,
        celltypist_results=celltypist_results,
        multi_peak_min_types=getattr(CFG.annotation, "multi_peak_min_types", 3),
        multi_peak_score_floor=getattr(CFG.annotation, "multi_peak_score_floor", 0.9),
        kadp_cfg=kadp_cfg,
        synonyms=synonyms,
        metc_cfg=metc_cfg,
    )
    logger.info("Evidence fusion: %d clusters processed", len(decisions))

    # Build cell_category_map using ORIGINAL all_scores (with Broad_* entries)
    cell_category_map = {}
    for cl_str in all_scores:
        category = _classify_broad_category(
            all_scores[cl_str], kb, allows_transitions=allows_transitions
        )
        cell_category_map[cl_str] = category

    if not decisions:
        logger.warning("Evidence fusion produced no decisions — falling back")
        return None

    # Build cluster -> decision mapping (preserving fusion sort order)
    decision_clusters = sorted(
        all_scores.keys(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )
    decision_map = dict(zip(decision_clusters, decisions))

    # ── Cell-type-aware category guard ─────────────────────────────────
    # Override the Fisher-score Broad_* winner with the canonical parent of
    # the chosen fine cell_type. Fixes systematic Broad_Neuron over-matching
    # for glial / immune / endothelial clusters (caused by weak Broad_Neuron
    # fallback markers TUBB3/ELAVL4/SLC17A6 hitting random housekeeping genes).
    _broad_parent_map: dict[str, str] = {}
    _hier = kb.get("_hierarchy") or {}
    for _cat_name, _cat_def in (_hier.get("categories") or {}).items():
        _broad_key = f"{CATEGORY_PREFIX}{_cat_name}"
        for _member in _cat_def.get("members") or []:
            _broad_parent_map[_member] = _broad_key

    if _broad_parent_map:
        _n_overrides = 0
        for _cl_str, _decision in decision_map.items():
            _fine_type = getattr(_decision, "cell_type", "") or ""
            _expected_broad = _broad_parent_map.get(_fine_type, "")
            _current_broad = cell_category_map.get(_cl_str, "")
            if _expected_broad and _current_broad != _expected_broad:
                # Respect the adult-tissue Broad_Progenitor exclusion:
                # only allow Progenitor override in developmental mode.
                if _expected_broad == f"{CATEGORY_PREFIX}Progenitor" and not allows_transitions:
                    continue
                logger.debug(
                    "Category guard: cluster %s (%s) %s → %s",
                    _cl_str,
                    _fine_type,
                    _current_broad,
                    _expected_broad,
                )
                cell_category_map[_cl_str] = _expected_broad
                _n_overrides += 1
        if _n_overrides:
            logger.info(
                "Category guard: overrode %d/%d cluster categories using fine cell_type parent",
                _n_overrides,
                len(decision_map),
            )

    # ── f. AI fallback for low-confidence clusters ────────────────────────
    ai_enabled = getattr(CFG.ai, "enabled", False)
    ai_annot_on = getattr(CFG.ai, "ai_annotation", False)

    # Two-segment selection: `_l1` is byte-identical to the baseline filter;
    # `_l2` additionally pulls ambiguous / transition_state candidates (incl.
    # confidence="transition") into the AI fallback when KADP or METC is
    # enabled.  With both flags off `_l2` is empty and the selection is
    # exactly the baseline.
    kadp_enabled = getattr(CFG.annotation, "kadp_enabled", False)
    metc_enabled = getattr(CFG.annotation, "metc_enabled", False)
    _l1 = [d for d in decisions if d.confidence in ("low", "unknown") and d.method != "ambiguous"]
    _l2 = (
        []
        if not (kadp_enabled or metc_enabled)
        else [
            d
            for d in decisions
            if d.confidence in ("low", "unknown", "transition")
            and d.method in ("ambiguous", "transition_state")
        ]
    )
    low_conf_clusters = _l1 + [d for d in _l2 if d not in _l1]
    ai_results = {}

    # Quality reported to _write_quality_report: the first-pass dict unless
    # the AI second pass re-fused, in which case the second pass's quality
    # wins (Oracle r3 MINOR 5) so 05_annotation_quality.json reflects the
    # AI-enhanced final decision_map.
    reported_quality = fusion_quality

    if low_conf_clusters and ai_enabled and ai_annot_on:
        logger.info("AI fallback for %d low-confidence clusters", len(low_conf_clusters))
        kb_candidates = sorted([k for k in kb if k != "expert_rules" and not k.startswith("_")])
        tissue = CFG.tissue
        stages_present = (
            sorted(adata.obs["stage"].unique().tolist()) if "stage" in adata.obs else []
        )
        extra_context = f"Developmental stages: {stages_present}" if stages_present else ""

        # Unconstrained annotations require build_annotation_prompt import here
        from core.ai.caller import ai_query
        from core.ai.prompts import build_annotation_prompt, build_hierarchical_annotation_prompt

        unconstrained = getattr(CFG.ai, "unconstrained_annotation", False)

        # Use hierarchical prompt if KB has _hierarchy section
        kb_hierarchy = kb.get("_hierarchy") if kb else None
        if kb_hierarchy:
            sys_prompt, user_prompt = build_hierarchical_annotation_prompt(
                adata=adata,
                tissue=tissue,
                species=species,
                kb_candidates=kb_candidates,
                kb_hierarchy=kb_hierarchy,
                precomputed_rank=True,
                extra_context=extra_context,
                compact=n_clusters > 20,
                unconstrained=unconstrained,
            )
        else:
            sys_prompt, user_prompt = build_annotation_prompt(
                adata,
                tissue,
                species,
                precomputed_rank=True,
                extra_context=extra_context,
                compact=n_clusters > 20,
                kb_candidates=kb_candidates,
                unconstrained=unconstrained,
            )

        try:
            response = ai_query(sys_prompt, user_prompt, cfg=CFG.ai)
            if response:
                ai_parsed = json.loads(response)
                for cid, ann in ai_parsed.items():
                    if isinstance(ann, dict) and "cell_type" in ann:
                        ai_results[str(cid)] = ann["cell_type"]
                logger.info("AI fallback: %d cluster suggestions received", len(ai_results))
                # AI labels share the SAME harmonization chain as CellTypist
                # (Oracle r2 MAJOR 1): each suggestion is pre-resolved to a
                # canonical KB name; a cluster whose label cannot be aligned
                # abstains (ai_suggestion blank — never inflates METC distinct).
                if ai_results:
                    _resolved_ai: dict = {}
                    for _cid, _label in ai_results.items():
                        _canon = harmonize_label(_label, kb, synonyms)
                        if _canon is not None:
                            _resolved_ai[_cid] = _canon
                    ai_results = _resolved_ai
                # Re-run fusion with AI context.  The second pass also
                # requests return_quality=True and its quality dict becomes
                # the reported one (Oracle r3 MINOR 5): 05_annotation_quality.json
                # must reflect the AI-enhanced final decision_map, not the
                # pre-AI first pass.
                decisions, reported_quality = fuse_all_clusters(
                    fine_scores,
                    all_rules,
                    kb=kb,
                    all_marker_dfs=marker_df,
                    ai_results=ai_results,
                    return_quality=True,
                    low_quality_clusters=low_quality_clusters,
                    unconstrained=getattr(CFG.ai, "unconstrained_annotation", False),
                    allows_transitions=allows_transitions,
                    incompatible_transitions=_incompatible_transitions,
                    celltypist_results=celltypist_results,
                    multi_peak_min_types=getattr(CFG.annotation, "multi_peak_min_types", 3),
                    multi_peak_score_floor=getattr(CFG.annotation, "multi_peak_score_floor", 0.9),
                    kadp_cfg=kadp_cfg,
                    synonyms=synonyms,
                    metc_cfg=metc_cfg,
                )
                decision_map = dict(zip(decision_clusters, decisions))
        except Exception as exc:
            logger.warning("AI fallback failed: %s — using pure KB results", exc)

    # ── g. Map decisions to adata.obs ─────────────────────────────────────
    leiden_str = adata.obs["leiden"].astype(str)

    # For low-quality clusters, downgrade: force Unknown + annotate reason.
    _forced_unknown = 0
    for cl_str, reason in low_quality_clusters.items():
        if cl_str in decision_map and decision_map[cl_str].confidence != "unknown":
            decision_map[cl_str] = decision_map[cl_str]._replace(
                cell_type="Unknown",
                confidence="unknown",
                method="unknown",
                explanation=(
                    "Low-quality cluster ({}) — {}".format(
                        reason, decision_map[cl_str].explanation[:120]
                    )
                ),
            )
            _forced_unknown += 1
    if _forced_unknown:
        logger.info(
            "Downgraded %d low-quality cluster(s) to Unknown: %s",
            _forced_unknown,
            ", ".join(
                "{} ({})".format(k, v)
                for k, v in low_quality_clusters.items()
                if k in decision_map
            ),
        )

    # For low-quality clusters, set cell_category to empty string
    for cl_str in low_quality_clusters:
        if cl_str in cell_category_map:
            cell_category_map[cl_str] = ""

    # ── Tiered annotation: enforce cell_type=L2, L3→cell_subtype ─────
    tiered_subtypes: dict[str, str] = {}
    tiered_candidates: dict[str, list] = {}
    if "_hierarchy" in kb and kb.get("_hierarchy", {}).get("categories"):
        from core.annotation.scoring import _build_kb_lookup
        from rna.utils.tiered_annotation import _build_subtype_map, resolve_tiered_label

        _kb_lookup = _build_kb_lookup(kb, species=species)
        _hierarchy = kb["_hierarchy"]
        _subtype_map = _build_subtype_map(_hierarchy)
        for cl_str in list(decision_map):
            # Tiered-block exemption (Oracle r3 MAJOR 2): a KADP decision
            # (method == "developmental_potency") skips the ENTIRE block —
            # tier / consensus / n_sources / subtype_resolution / L2-forcing
            # all keep their KADP defaults.  The cluster still gets a
            # tiered_subtypes default via the else branch below (N/A) and a
            # .get(k, []) fallback in tiered_candidates consumers.
            if decision_map[cl_str].method == "developmental_potency":
                continue
            _scores = all_scores.get(cl_str, {})
            _label, _ev = resolve_tiered_label(
                _scores,
                _hierarchy,
                _kb_lookup,
                all_top_genes.get(cl_str, []),
            )
            tiered_candidates[cl_str] = _ev.get("subtype_candidates", [])
            # Force fusion cell_type to its L2 ancestor.
            _ct = decision_map[cl_str].cell_type
            while _ct in _subtype_map:
                _ct = _subtype_map[_ct]
            # cell_subtype from tiered resolution.
            if _ev["tier"] == "L3" and _ev["subtype_resolution"] == "resolved":
                tiered_subtypes[cl_str] = _label
            elif _ev["subtype_resolution"] == "unresolved":
                tiered_subtypes[cl_str] = "unresolved"
            else:
                tiered_subtypes[cl_str] = "N/A"
            decision_map[cl_str] = decision_map[cl_str]._replace(
                cell_type=_ct,
                tier=("L3" if tiered_subtypes[cl_str] not in ("unresolved", "N/A") else "L2"),
                consensus=_ev["consensus"],
                n_sources=_ev["n_sources"],
                subtype_resolution=_ev["subtype_resolution"] or "na",
            )
    else:
        tiered_subtypes = {k: "N/A" for k in decision_map}

    # ── D3 evidence-strength audit: canonical-expression fallback + ai_only ──
    # Runs on the final (post-tiering) decisions, before any obs column is
    # written, so a downgrade propagates to cell_state / annot_confidence /
    # ann_records / review_queue downstream.  Both gates mutate decision_map
    # in place via NamedTuple._replace and record their reason on the
    # decision's ``review_reason`` field (the intermediate structure task 10
    # surfaces in the quality-report review_queue).
    review_reasons = _apply_canonical_expression_fallback(adata, kb, decision_map, CFG, logger)
    review_reasons.update(_flag_ai_only_decisions(decision_map, logger))
    if review_reasons:
        logger.info(
            "Evidence-strength audit: %d decision(s) flagged for review: %s",
            len(review_reasons),
            ", ".join(f"{cl}={reason}" for cl, reason in sorted(review_reasons.items())),
        )

    adata.obs["cell_type"] = leiden_str.map(
        {k: v.cell_type for k, v in decision_map.items()}
    ).astype("category")
    adata.obs["cell_state"] = leiden_str.map(
        {k: _map_cell_state(v, cell_category_map.get(k, "")) for k, v in decision_map.items()}
    )
    adata.obs["cell_subtype"] = leiden_str.map(
        {k: tiered_subtypes.get(k, "N/A") for k in decision_map}
    )
    adata.obs["cell_category"] = leiden_str.map(
        {k: cell_category_map.get(k, "") for k in decision_map}
    ).astype("category")

    # annot_method: clean label from fusion method (+ AI suffix)
    def _clean_method(d):
        if d.method == "expert_rule":
            return "expert_rule"
        if d.method == "unknown":
            return "unknown"
        if d.method == "ambiguous":
            return "ambiguous"
        if d.method == "transition_state":
            return "transition_state"
        # KADP (layer 3): must win over the marker_scoring+ai suffix — the
        # AI second pass may agree with the KADP name (ai_agreed=True), yet
        # the annot_method must stay developmental_potency.
        if d.method == "developmental_potency":
            return "developmental_potency"
        if d.ai_agreed or d.ai_suggested:
            return "marker_scoring+ai"
        return "marker_scoring"

    adata.obs["annot_method"] = leiden_str.map(
        {k: _clean_method(v) for k, v in decision_map.items()}
    )
    adata.obs["annot_confidence"] = leiden_str.map(
        {k: v.confidence for k, v in decision_map.items()}
    )
    adata.obs["annot_reasoning"] = leiden_str.map(
        {k: v.explanation for k, v in decision_map.items()}
    )
    # ── Persist tiered near-miss subtype candidates (additive, backward compatible) ──
    # Lazily imported here — a module-top import would create a static
    # core/ → rna/utils dependency (the JSON/CSV blocks are OUTSIDE the
    # tiered block where the other tiered helpers are imported).
    if tiered_candidates:
        from rna.utils.tiered_annotation import format_subtype_candidates

    adata.obs["annot_evidence"] = leiden_str.map(
        {
            k: json.dumps(
                {
                    "score": v.score,
                    "method": v.method,
                    "n_markers_found": v.n_markers_found,
                    "ai_agreed": v.ai_agreed,
                    "ai_suggested": v.ai_suggested,
                    "diagnostic_category": v.diagnostic.category if v.diagnostic else None,
                    "diagnostic_detail": v.diagnostic.detail if v.diagnostic else None,
                    "top_competitors": v.diagnostic.top_competitors if v.diagnostic else [],
                    "tier": v.tier,
                    "consensus": v.consensus,
                    "n_sources": v.n_sources,
                    "subtype_resolution": v.subtype_resolution,
                    "cell_subtype": tiered_subtypes.get(k, "N/A"),
                    "potency": v.potency,
                    # Layer-4 METC (todo 10): the per-source vote dict
                    # {marker/expert/ai/celltypist: label}, serialized inside
                    # this JSON string exactly like ``potency`` (None/omitted
                    # for non-METC decisions).  cell_metadata.csv has no
                    # source_votes column — annot_evidence only.
                    "source_votes": v.source_votes,
                }
                | (
                    {"subtype_candidates": tiered_candidates.get(k, [])}
                    if tiered_candidates
                    else {}
                )
            )
            for k, v in decision_map.items()
        }
    )

    # ── h. Save annotation CSV ────────────────────────────────────────────
    ann_records = []

    def sort_key(x):
        return int(x) if str(x).isdigit() else str(x)

    for cl_name in sorted(decision_map.keys(), key=sort_key):
        d = decision_map[cl_name]
        record = {
            "cluster": cl_name,
            "cell_type": d.cell_type,
            "confidence": d.confidence,
            "method": _clean_method(d),
            "score": d.score,
            "n_markers_found": d.n_markers_found,
            "ai_agreed": d.ai_agreed,
            "ai_suggested": d.ai_suggested,
            "reasoning": d.explanation,
            "diagnostic_category": d.diagnostic.category if d.diagnostic else "",
            "cell_category": cell_category_map.get(str(cl_name), ""),
            "cell_subtype": tiered_subtypes.get(str(cl_name), "N/A"),
            "tier": d.tier,
            "consensus": d.consensus,
            "n_sources": d.n_sources,
            "subtype_resolution": d.subtype_resolution,
            "potency": json.dumps(d.potency) if d.potency else "",
        }
        if tiered_candidates:
            record["subtype_candidates"] = format_subtype_candidates(
                tiered_candidates.get(str(cl_name), [])
            )
        ann_records.append(record)
    ann_df = pd.DataFrame(ann_records)
    ann_csv = os.path.join(CFG.table_dir, "cell_type_annotations.csv")
    ann_df.to_csv(ann_csv, index=False)
    logger.info("Annotation table saved: %s", ann_csv)

    logger.info("Cluster → cell type mapping (Unified):")
    for rec in ann_records:
        category_str = f" [{rec.get('cell_category', '')}]" if rec.get("cell_category") else ""
        logger.info(
            "  Cluster %s → %s%s (conf=%s, method=%s)",
            rec["cluster"],
            rec["cell_type"],
            category_str,
            rec["confidence"],
            rec["method"],
        )
    for rec in ann_records:
        logger.info(
            "  Cluster %s → %s (conf=%s, method=%s)",
            rec["cluster"],
            rec["cell_type"],
            rec["confidence"],
            rec["method"],
        )

    # ── i. UMAP visualization ─────────────────────────────────────────────
    sc.settings.figdir = os.path.join(CFG.figure_dir, "05_annotation")
    os.makedirs(sc.settings.figdir, exist_ok=True)
    sc.settings.autoshow = False

    adata.obs["annot_label"] = adata.obs["cell_type"].astype(str)

    safe_plot(sc.pl.umap, adata, cfg=CFG, color="cell_type", show=False, save="celltype_umap.pdf")
    safe_plot(
        sc.pl.umap, adata, cfg=CFG, color="annot_label", show=False, save="annot_label_umap.pdf"
    )
    safe_plot(
        sc.pl.umap,
        adata,
        cfg=CFG,
        color="annot_confidence",
        show=False,
        save="annot_confidence_umap.pdf",
    )
    safe_plot(
        sc.pl.umap,
        adata,
        cfg=CFG,
        color="cell_category",
        show=False,
        save="cell_category_umap.png",
    )

    # ── j. Cell metadata export ───────────────────────────────────────────
    # KADP potency (layer 3, plan annotation-kadp-metc): single column, JSON
    # string of the three-value dict, empty string when no potency present.
    potency_map = {
        k: (json.dumps(v.potency) if v.potency else "") for k, v in decision_map.items()
    }
    meta_df = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "UMAP_1": adata.obsm["X_umap"][:, 0],
            "UMAP_2": adata.obsm["X_umap"][:, 1],
            "cell_type": adata.obs["cell_type"].values,
            "cell_state": adata.obs["cell_state"].values,
            "cell_subtype": adata.obs["cell_subtype"].values,
            "annot_confidence": adata.obs["annot_confidence"].values,
            "annot_method": adata.obs["annot_method"].values,
            "cell_category": adata.obs["cell_category"].values,
            "potency": leiden_str.map(potency_map).values,
        }
    )
    meta_csv = os.path.join(CFG.table_dir, "cell_metadata.csv")
    meta_df.to_csv(meta_csv, index=False)
    logger.info("Cell metadata exported: %s", meta_csv)

    # ── k. Annotation quality report ────────────────────────────────────────
    _write_quality_report(
        adata,
        ann_records,
        reported_quality,
        cell_category_map,
        decision_map,
        CFG,
        logger,
        kb=kb,
    )

    # ── l. Interactive review (--interactive flag) ──────────────────────────
    if getattr(CFG, "interactive", False):
        _interactive_annotation_review(adata, reported_quality, CFG, logger)

    return decision_map


def _write_quality_report(
    adata,
    ann_records,
    fusion_quality,
    cell_category_map,
    decision_map,
    CFG,  # noqa: N803
    logger,
    kb: dict | None = None,
):
    """Write 05_annotation_quality.json summarising annotation health."""
    pass_cells = (
        (adata.obs["marker_validation"] == "PASS").sum() if "marker_validation" in adata.obs else 0
    )
    pass_rate = pass_cells / max(adata.n_obs, 1)

    ambiguity_clusters = []
    for rec in ann_records:
        reasoning = rec.get("reasoning", "")
        if "also matched rules:" in reasoning:
            ambiguity_clusters.append(rec["cluster"])

    # D8 — KB coverage: annotated types vs KB fine types vs ghost transition
    # endpoints.  F14: these ghost endpoints (transition-pair members that are
    # not annotated cell types) are REPORT-ONLY — they never enter KADP
    # precursor naming.  That is structurally enforced upstream by the
    # score > 0 filtering in core/annotation/potency.py (a ghost endpoint
    # has no marker scores, so it can never be a progenitor-pole argmax).
    annotated_types = {d.cell_type for d in decision_map.values()}
    transition_pairs = {
        tuple(d.cell_type[len("transitional: ") :].split("/"))
        for d in decision_map.values()
        if d.method == "transition_state"
    }
    ghost_endpoints = {t for p in transition_pairs for t in p} - annotated_types
    kb_fine_types = (
        {
            k
            for k in kb
            if not k.startswith("_") and k != "expert_rules" and not k.startswith("Broad_")
        }
        if kb
        else set()
    )
    kb_types_unannotated = kb_fine_types - annotated_types

    quality = {
        "pass_rate": round(pass_rate, 4),
        "total_clusters": len(ann_records),
        "annotated_by_rule": fusion_quality.get("annotated_by_rule", 0),
        "annotated_by_scoring": fusion_quality.get("annotated_by_scoring", 0),
        "unknown": fusion_quality.get("unknown", 0),
        # Layer-4 METC (todo 8/10): fraction of clusters with raw CellTypist
        # labels that aligned to a canonical KB name.  None when no raw labels
        # exist (Oracle r3 MINOR 6 — never a division by zero).
        "harmonization_rate": fusion_quality.get("harmonization_rate"),
        "ambiguity_clusters": ambiguity_clusters,
        "ai_disagreement_rate": round(
            sum(1 for r in ann_records if not r.get("ai_agreed", True)) / max(len(ann_records), 1),
            4,
        ),
        "kb_blind_spot": pass_rate < 0.1,
        "recommended_strictness": (
            "relaxed" if pass_rate < 0.1 else "deep" if pass_rate < 0.3 else "default"
        ),
        "categories_found": len(set(c for c in cell_category_map.values() if c)),
        "transition_clusters": [
            {"cluster": cl, "pair": d.cell_type[len("transitional: ") :]}
            for cl, d in decision_map.items()
            if d.method == "transition_state"
        ],
        "review_queue": [
            {
                "cluster": cl,
                # Task 10 (D6): every entry carries a review ``reason``.
                # decision.review_reason (the downgrade source) takes
                # precedence; a plain multi-peak ambiguity decision (no
                # review_reason) gets the canonical "ambiguous".
                "reason": getattr(d, "review_reason", "") or "ambiguous",
                # Tie detail is meaningful for multi-peak ambiguity AND for
                # METC entries (metc_divergent / metc_2way carry a ranked
                # top_competitors list of {cell_type, score} dicts — F13).
                # Kadp/consensus entries keep empty values.
                "n_tied_types": len(d.diagnostic.top_competitors)
                if d.diagnostic
                and (d.method == "ambiguous" or (d.review_reason or "").startswith("metc"))
                else 0,
                "top_types": [c["cell_type"] for c in d.diagnostic.top_competitors]
                if d.diagnostic
                and (d.method == "ambiguous" or (d.review_reason or "").startswith("metc"))
                else [],
            }
            for cl, d in decision_map.items()
            # D3: decisions carrying a review_reason (canonical-expression
            # fallback / ai_only / weak evidence) also enter the review queue.
            if d.method == "ambiguous" or getattr(d, "review_reason", "")
        ],
        "kb_coverage": {
            "annotated_types": sorted(annotated_types),
            "kb_types_unannotated": sorted(kb_types_unannotated),
            "ghost_endpoints": sorted(ghost_endpoints),
        },
        "category_distribution": {
            cat: count for cat, count in Counter(cell_category_map.values()).items() if cat
        },
    }

    quality_path = os.path.join(CFG.table_dir, "05_annotation_quality.json")
    with open(quality_path, "w", encoding="utf-8") as f:
        json.dump(quality, f, indent=2)
    logger.info(
        "Annotation quality report: %s (pass_rate=%.1f%%)",
        quality_path,
        quality["pass_rate"] * 100,
    )


def _interactive_annotation_review(adata, fusion_quality, CFG, logger):  # noqa: N803
    """Present annotation quality summary and offer remediation choices.

    Only called when ``CFG.interactive`` is ``True``.
    """
    pass_cells = (
        (adata.obs["marker_validation"] == "PASS").sum() if "marker_validation" in adata.obs else 0
    )
    pass_rate = pass_cells / max(adata.n_obs, 1)
    n_total = fusion_quality.get("total", 0)
    n_rule = fusion_quality.get("annotated_by_rule", 0)
    n_scoring = fusion_quality.get("annotated_by_scoring", 0)
    n_unknown = fusion_quality.get("unknown", 0)
    n_ambiguity = fusion_quality.get("ambiguity", 0)

    print("\n" + "=" * 60)
    print("Annotation Quality Summary")
    print("=" * 60)
    print(f"  PASS rate:          {pass_rate * 100:.1f}%")
    print(f"  Annotated:          {n_rule} by rule, {n_scoring} by scoring")
    print(f"  Unknown:            {n_unknown}/{n_total}")
    if n_ambiguity > 0:
        print(f"  ⚠  High ambiguity:  {n_ambiguity} cluster(s) matched ≥3 rules")
    if pass_rate < 0.1:
        print("  ⚠  KB blind spot detected")
        rec = "relaxed" if pass_rate < 0.1 else "deep" if pass_rate < 0.3 else "default"
        print(f"  💡 Recommended:       strictness='{rec}'")
    print()

    if pass_rate < 0.1:
        rec = "relaxed" if pass_rate < 0.1 else "deep"
        try:
            choice = (
                input(
                    "KB coverage is very low on this dataset. Options:\n"
                    f"  [r] Re-annotate with strictness='{rec}'\n"
                    "  [s] Continue with score_genes fallback\n"
                    "  [c] Continue with current labels (not recommended)\n"
                    "  [a] Abort\n"
                    "Choice> "
                )
                .strip()
                .lower()
            )
        except (EOFError, KeyboardInterrupt):
            logger.warning("Interactive input interrupted — continuing")
            return

        if choice == "r":
            logger.info(
                "User chose: re-annotate with strictness='%s'",
                rec,
            )
            print(
                f"\nTo re-annotate, set in your config:\n"
                f"  CFG.marker.expert_rule_strictness = '{rec}'\n"
                f"Or pass --config with the updated setting.\n"
            )
        elif choice == "s":
            logger.info("User chose: score_genes fallback")
            print("Set CFG.tissue_kb = '' to use score_genes mode.\n")
        elif choice == "a":
            raise SystemExit("Aborted by user.")
        else:
            logger.info("User chose: continue with current labels")
