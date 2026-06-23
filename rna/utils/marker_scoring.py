"""
utils/marker_scoring.py — Marker-based cell type scoring engine.

Purely computational annotation module (no LLM calls) that uses
hypergeometric enrichment tests and cosine similarity to match
cluster marker genes against an expert Knowledge Base (KB).

The KB is expected to have the following structure::

    kb = {
        "CellTypeA": {
            "markers": {
                "confirm": {"GENE1": ["PMID1", "PMID2"], "GENE2": ["PMID3"]},
                "add":    {"GENE3": ["PMID1"]},
                "refine": {"GENE4": ["PMID2"]},
            },
            "negative_markers": ["GENE5", "GENE6"],
            "species": ["human", "mouse"],
            "synonyms": ["TypeA", "Type A"],
            "parent": "LineageX",
        },
        "CellTypeB": ...,
        "expert_rules": [
            {
                "priority": 10,
                "condition": {
                    "markers_present": {"GENE1": 1.0},
                    "markers_absent": ["GENE5"],
                },
                "action": "CellTypeA",
            },
        ],
    }
"""

from typing import Any, Dict, List, NamedTuple, Optional
import logging
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
#  Data structures
# ═══════════════════════════════════════════════════════════════════════


class Score(NamedTuple):
    """Scoring result for one cluster–cell-type pair.

    Attributes
    ----------
    score : float
        Final composite score in [0, 1].
    p_value : float
        P-value from the hypergeometric (Fisher's exact) test.
    method : str
        Which scoring method contributed the max value:
        ``"hypergeometric"``, ``"cosine"``, or ``"none"``.
    n_markers_found : int
        Number of KB positive markers found in the cluster's top-20.
    negative_penalty : bool
        Whether a negative-marker penalty was applied.
    """
    score: float
    p_value: float
    method: str
    n_markers_found: int
    negative_penalty: bool


# ═══════════════════════════════════════════════════════════════════════
#  KB helpers
# ═══════════════════════════════════════════════════════════════════════


def _get_canonical_markers(kb: Dict[str, Any], type_key: str,
                           species: Optional[str] = None,
                           gene_names: Optional[List[str]] = None) -> List[str]:
    """Return the union of all *confirm* + *add* markers for *type_key*.

    If *species* is provided and the type's ``species`` list does not contain
    it, an empty list is returned (the type is considered irrelevant).

    However, if *gene_names* is provided and the majority of genes appear to
    already be human orthologs (i.e. ortholog conversion was applied), the
    species filter is **relaxed** — markers are returned so KB scoring can
    work on cross-species data that has been mapped to human gene symbols.
    """
    type_data = kb.get(type_key, {})
    if not type_data:
        return []

    markers_dict = type_data.get("markers", {})
    confirm = set(markers_dict.get("confirm", {}).keys())
    add = set(markers_dict.get("add", {}).keys())
    result = confirm | add

    if species:
        type_species = type_data.get("species", [])
        if species not in type_species:
            # Relax species filter if genes are mapped to human orthologs
            if gene_names and _looks_mapped_to_target(gene_names):
                logger.debug(
                    "Species '%s' not in KB for type '%s', but genes appear "
                    "mapped to target — keeping markers", species, type_key,
                )
                return list(result)
            return []

    return list(result)


def _looks_mapped_to_target(gene_names: List[str],
                            unmapped_prefix: str = "UNMAPPED_") -> bool:
    """Heuristic: check if gene names were already ortholog-mapped.

    Returns True if >= 75% of gene names look like standard symbols (not
    Ensembl IDs, not UNMAPPED_ prefixed), suggesting ortholog conversion
    was applied.
    """
    if not gene_names:
        return False
    sample = gene_names[:200]
    n_unmapped = sum(1 for g in sample if str(g).startswith(unmapped_prefix))
    n_ensembl = sum(1 for g in sample if re.match(r'^ENS[A-Z]{0,4}G\d{11}$', str(g)))
    ratio_mapped = 1.0 - (n_unmapped + n_ensembl) / len(sample)
    return ratio_mapped >= 0.75


def _negative_marker_penalty(kb: Dict[str, Any], type_key: str,
                             cluster_markers: pd.DataFrame) -> bool:
    """Return ``True`` if the cluster expresses >= 2 negative markers for *type_key*.

    Both raw KB (``negative_markers`` key) and lookup format (``negative``
    key) are accepted.
    """
    type_data = kb.get(type_key, {})
    if not type_data:
        return False

    neg_markers = (
        type_data.get("negative_markers")
        or type_data.get("negative")
        or []
    )
    if not neg_markers:
        return False

    top10 = set(cluster_markers.head(10)["names"].tolist())
    found = sum(1 for m in neg_markers if m in top10)
    return found >= 2


def _build_kb_lookup(kb: Dict[str, Any],
                     species: Optional[str] = None) -> Dict[str, Any]:
    """Convert a raw KB into a flat lookup dict for fast scoring.

    Returns
    -------
    dict
        ``{type_key: {"positive": [...], "negative": [...],
                       "species": [...], "synonyms": [...],
                       "parent": str, "marker_weights": {gene: count}}}``

    **Consensus merging**

    Each positive marker is associated with a list of supporting sources
    (e.g. PMIDs).  Markers appearing in **>= 2 sources** receive higher
    weight, which is stored in ``marker_weights``.  The ``positive`` list
    always contains the union of *confirm* + *add* markers regardless of
    source count; callers may use ``marker_weights`` for weighted scoring
    if desired.

    Duplicates are removed and the optional *species* filter is applied.
    """
    kb_all: Dict[str, Any] = {}

    for type_key, type_data in kb.items():
        if type_key == "expert_rules":
            continue

        positive = _get_canonical_markers(kb, type_key, species)
        negative = list(type_data.get("negative_markers", []))

        # Build source counts for consensus weighting.
        marker_weights: Dict[str, int] = {}
        for tier in ("confirm", "add", "refine"):
            tier_markers = type_data.get("markers", {}).get(tier, {})
            for gene, sources in tier_markers.items():
                n = len(sources) if isinstance(sources, (list, tuple)) else 1
                marker_weights[gene] = marker_weights.get(gene, 0) + n

        kb_all[type_key] = {
            "positive": positive,
            "negative": negative,
            "species": list(type_data.get("species", [])),
            "synonyms": list(type_data.get("synonyms", [])),
            "parent": type_data.get("parent", ""),
            "marker_weights": marker_weights,
        }

    return kb_all


# ═══════════════════════════════════════════════════════════════════════
#  Core scoring
# ═══════════════════════════════════════════════════════════════════════


def score_cluster_against_kb(kb: Dict[str, Any],
                             cluster_markers: pd.DataFrame,
                             species: Optional[str] = None
                             ) -> Dict[str, Score]:
    """Score one cluster against every cell type in the Knowledge Base.

    Parameters
    ----------
    kb : dict
        Either the raw KB (with ``markers``/``confirm``/``add``/``refine``
        structure) or a pre-built lookup from :func:`_build_kb_lookup`.
    cluster_markers : pd.DataFrame
        Top-N markers for a **single** cluster.  Must have columns
        ``names``, ``logfoldchanges``, ``pvals_adj``.  The first 20 rows
        are used as the cluster's marker set.
    species : str or None
        If set, only cell types matching this species are scored.

    Returns
    -------
    Dict[str, Score]
        Mapping from ``type_key`` → :class:`Score` for each cell type.
    """
    # ── Normalise input ─────────────────────────────────────────────
    sample_key = next((k for k in kb if k != "expert_rules"), None)
    is_raw = bool(sample_key and "markers" in kb[sample_key])

    if is_raw:
        kb_lookup = _build_kb_lookup(kb, species)
    else:
        kb_lookup = kb

    # Background = union of ALL positive markers across every type in KB.
    all_type_markers: set[str] = set()
    for entry in kb_lookup.values():
        all_type_markers.update(entry.get("positive", []))
    background_size = max(len(all_type_markers), 1)

    # Top-20 markers for this cluster.
    top_n = 20
    top_markers = cluster_markers.head(top_n).copy()
    top_gene_set = set(top_markers["names"].tolist())

    # Which top-20 genes fall within the KB background (union of all markers)?
    # Only these enter the Fisher contingency table — genes outside the
    # known marker universe are irrelevant for the enrichment test.
    top_in_bg = top_gene_set & all_type_markers
    n_top_in_bg = len(top_in_bg)

    results: Dict[str, Score] = {}

    for type_key, entry in kb_lookup.items():
        positive_markers = entry.get("positive", [])
        if not positive_markers:
            results[type_key] = Score(0.0, 1.0, "none", 0, False)
            continue

        positive_set = set(positive_markers)

        # ── 1. Hypergeometric (Fisher's exact) score ────────────────
        # Contingency table over the KB background genes:
        #                  | In top-20 | Not in top-20
        # -----------------|-----------|--------------
        # Is type marker   | a         | b
        # Not type marker  | c         | d
        a = len(positive_set & top_in_bg)   # type markers in top-20
        b = len(positive_set) - a           # type markers NOT in top-20
        c = n_top_in_bg - a                  # non-type KB markers in top-20
        d = max(background_size - a - b - c, 1)  # remaining KB markers
        if a > 0 and b >= 0 and c >= 0 and d > 0:
            table = [[a, b], [c, d]]
            _r = fisher_exact(table, alternative='greater')
            raw_p = float(str(_r[1]))
        else:
            raw_p = 1.0
        hypergeometric_score = 1.0 - raw_p

        # ── 2. Cosine similarity score ──────────────────────────────
        all_genes = list(set(top_markers["names"].tolist() + positive_markers))
        cluster_vec = np.array(
            [1.0 if g in top_gene_set else 0.0 for g in all_genes]
        )
        type_vec = np.array(
            [1.0 if g in positive_set else 0.0 for g in all_genes]
        )

        cluster_norm = np.linalg.norm(cluster_vec)
        type_norm = np.linalg.norm(type_vec)

        if cluster_norm > 0 and type_norm > 0:
            cos_sim = float(np.dot(cluster_vec, type_vec)
                            / (cluster_norm * type_norm))
        else:
            cos_sim = 0.0

        # ── 3. Confidence multiplier ────────────────────────────────
        n_type_markers = len(positive_markers)
        if n_type_markers > 5:
            conf_mult = 1.0
        elif n_type_markers >= 3:
            conf_mult = 0.8
        else:
            conf_mult = 0.5

        # ── 4. Combine ──────────────────────────────────────────────
        if hypergeometric_score >= cos_sim:
            base_method = "hypergeometric"
            base_score = hypergeometric_score
        else:
            base_method = "cosine"
            base_score = cos_sim

        final_score = base_score * conf_mult

        # ── 5. Negative-marker penalty ──────────────────────────────
        neg_penalty = _negative_marker_penalty(
            kb if is_raw else kb_lookup, type_key, cluster_markers
        )
        if neg_penalty:
            final_score *= 0.5

        results[type_key] = Score(
            score=final_score,
            p_value=float(raw_p),
            method=base_method,
            n_markers_found=a,
            negative_penalty=neg_penalty,
        )

    return results


def annotate_all_clusters(kb_all_markers: Dict[str, Any],
                          all_marker_dfs: pd.DataFrame,
                          species: str) -> pd.DataFrame:
    """Score every cluster and assign the best-matching cell type.

    Parameters
    ----------
    kb_all_markers : dict
        Pre-built lookup dict from :func:`_build_kb_lookup`.
    all_marker_dfs : pd.DataFrame
        Concatenated ``rank_genes_groups`` output for **all** clusters.
        Must have columns ``names``, ``logfoldchanges``, ``pvals_adj``,
        and **``cluster``** (the cluster identifier).
    species : str
        Species filter passed through to :func:`_build_kb_lookup`.

    Returns
    -------
    pd.DataFrame
        Columns: ``["cluster", "cell_type", "score", "p_value",
                     "method", "n_markers_found"]``

        Only clusters with a best score >= 0.25 are included.
    """
    records: List[Dict[str, Any]] = []

    clusters = sorted(
        all_marker_dfs["cluster"].unique(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )
    for cl in clusters:
        cl_mask = all_marker_dfs["cluster"] == cl
        cl_sort = all_marker_dfs[cl_mask].copy()
        # Sort by logfoldchanges descending.
        lfc_idx = cl_sort["logfoldchanges"].argsort()[::-1]
        cl_sort = cl_sort.iloc[lfc_idx]

        scores = score_cluster_against_kb(
            kb_all_markers, cl_sort, species=species
        )

        if not scores:
            continue

        # Pick the type with the highest score.
        best_type = max(scores, key=lambda k: scores[k].score)
        best_score = scores[best_type]

        if best_score.score >= 0.25:
            records.append({
                "cluster": cl,
                "cell_type": best_type,
                "score": best_score.score,
                "p_value": best_score.p_value,
                "method": best_score.method,
                "n_markers_found": best_score.n_markers_found,
            })

    return pd.DataFrame(records)


# ═══════════════════════════════════════════════════════════════════════
#  Expert-rules engine
# ═══════════════════════════════════════════════════════════════════════


def apply_expert_rules(kb: Dict[str, Any],
                       cluster_markers: pd.DataFrame) -> Optional[str]:
    """Try to deterministically match *cluster_markers* via expert rules.

    Iterates over ``kb["expert_rules"]`` (sorted by ``priority`` descending).
    The first rule whose conditions are all satisfied wins.

    Each rule has the structure::

        {
            "priority": int,
            "condition": {
                "markers_present": {"GENE": min_logFC, ...},
                "markers_absent": ["GENE", ...],   # optional
            },
            "action": "CellTypeKey",
        }

    Parameters
    ----------
    kb : dict
        Raw KB (must contain an ``"expert_rules"`` list).
    cluster_markers : pd.DataFrame
        Marker DataFrame for a single cluster (columns: ``names``,
        ``logfoldchanges``).

    Returns
    -------
    str or None
        The matched ``action`` cell-type key, or ``None`` if no rule fires.
    """
    rules = kb.get("expert_rules", [])
    if not rules:
        return None

    # Sort by priority descending (higher = more specific).
    sorted_rules = sorted(
        rules, key=lambda r: r.get("priority", 0), reverse=True
    )
    # Build a fast lookup: gene_name -> logfoldchanges.
    marker_map: Dict[str, float] = {}
    for _, row in cluster_markers.iterrows():
        gene_name = str(row["names"])
        marker_map[gene_name] = float(row["logfoldchanges"])

    cluster_genes = set(cluster_markers["names"].tolist())

    for rule in sorted_rules:
        condition = rule.get("condition", {})
        markers_present: Dict[str, float] = condition.get("markers_present", {})
        markers_absent: List[str] = condition.get("markers_absent", [])

        # All required markers must be in cluster markers at sufficient logFC.
        passed = True
        for gene, min_logfc in markers_present.items():
            if gene not in marker_map or marker_map[gene] < min_logfc:
                passed = False
                break

        if not passed:
            continue

        # No exclusion markers should be present.
        for gene in markers_absent:
            if gene in cluster_genes:
                passed = False
                break

        if passed:
            return rule.get("action")

    return None
