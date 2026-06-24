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
            "class": "Mammalia",          # taxonomic class (纲)
            "order": "Primates",          # taxonomic order (目)
            "classes": ["Mammalia"],      # all contributing classes
            "orders": ["Primates"],       # all contributing orders
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
        "_meta": {
            "total_sources": N,
            "classes": ["Mammalia"],
            "orders": ["Primates"],
        },
    }

**Phylogenetic weighting** (v3.0.0+): When ``target_class`` and/or
``target_order`` are passed, marker scores receive a multiplicative
phylogenetic weight:

    ========================  =====
    Source → Target            Weight
    ========================  =====
    Same class + same order    1.0
    Same class, different order 0.8
    Different class, multi-class marker  0.9
    Different class, single-class marker 0.6
    ========================  =====

This penalises markers that only appear in a distant taxonomic group
while preserving the power of markers that are **conserved across
multiple classes** (a strong signal of biological relevance).
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
            # Try common-name ↔ scientific-name normalisation.
            if _species_matches(species, type_species):
                # Species matches after normalisation — keep markers.
                return list(result)
            # Species does not match at all.
            # Relax filter if genes are mapped to human orthologs.
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


# Common-name → scientific-name mappings used by _species_matches().
# Extend this table as new species are added to the KB.
_SPECIES_SYNONYMS: Dict[str, str] = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "macaque": "Macaca fascicularis",
    "marmoset": "Callithrix jacchus",
    "zebrafish": "Danio rerio",
    "chicken": "Gallus gallus",
    "lamprey": "Petromyzon marinus",
    "frog": "Xenopus laevis",
    "pig": "Sus scrofa",
    "cow": "Bos taurus",
    "sheep": "Ovis aries",
    "ferret": "Mustela putorius furo",
    "squirrel": "Ictidomys tridecemlineatus",
    "opossum": "Didelphis marsupialis",
    "treeshrew": "Tupaia belangeri",
    "anolis": "Anolis sagrei",
    "deer_mouse": "Peromyscus maniculatus",
    "striped_mouse": "Rhabdomys pumilio",
}


def _species_matches(user_species: str, kb_species_list: list[str]) -> bool:
    """Return ``True`` if *user_species* matches any entry in *kb_species_list*.

    Normalises common-name aliases (e.g. ``"human"`` → ``"Homo sapiens"``)
    and performs case-insensitive comparison.
    """
    if not user_species or not kb_species_list:
        return False
    normalised_user = _SPECIES_SYNONYMS.get(
        user_species.strip().lower(), user_species.strip()
    )
    for ks in kb_species_list:
        normalised_ks = _SPECIES_SYNONYMS.get(
            ks.strip().lower(), ks.strip()
        )
        if normalised_user.lower() == normalised_ks.lower():
            return True
    return False


def _negative_marker_penalty(kb: Dict[str, Any], type_key: str,
                             cluster_markers: pd.DataFrame) -> bool:
    """Return ``True`` if the cluster expresses >= 2 negative markers for *type_key*.

    Both raw KB (``negative_markers`` key) and lookup format (``negative``
    key) are accepted.
    """
    type_data = _type_data(kb, type_key)
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


def _type_data(kb: Dict[str, Any], type_key: str) -> Dict[str, Any]:
    """Return the raw type entry from *kb*, regardless of raw vs lookup format.

    In raw KB the entry is ``kb[type_key]``.  In lookup format the entry may
    be flattened; we return the raw entry when available, otherwise an empty
    dict.
    """
    raw_entry = kb.get(type_key, {})
    # Raw KB entries have "markers" key; lookup entries have "positive".
    if "markers" in raw_entry or "positive" in raw_entry:
        return raw_entry
    return {}


def phylogenetic_weight(source_class: str,
                         target_class: str,
                         target_order: str = "",
                         source_order: str = "",
                         source_classes_contrib: list | None = None
                         ) -> float:
    """Return a multiplicative weight based on taxonomic distance.

    The weight reflects how relevant a KB cell-type marker set is for
    annotating data from a target taxonomic class (纲) / order (目).

    Parameters
    ----------
    source_class : str
        The primary class of the KB cell type (e.g. ``"Mammalia"``).
        Empty string means the class is unknown — weight falls through to
        the default (0.8) permissive value.
    target_class : str
        The desired class (e.g. ``"Mammalia"``).  An empty string means
        *no phylogenetic filtering* — weight is always 1.0.
    target_order : str
        The desired order (e.g. ``"Primates"``).  Only applied when
        *target_class* is also non-empty.
    source_order : str
        The primary order of the KB cell type.  Compared with
        *target_order* for intra-class fine-tuning.
    source_classes_contrib : list[str] or None
        All classes (across sources) that contributed markers to this type.
        If a target class appears here (e.g. because this type was observed
        across multiple classes in a cross-species study), it is treated as
        *same class* (weight 1.0).  When ``None``, only ``source_class``
        is used for the class check.

    Returns
    -------
    float
        Weight in [0.0, 1.0].

    Weight table
    ------------
    ================================  =====
    Source → Target                    Weight
    ================================  =====
    Same class + same order           1.0
    Same class + different order      0.8
    Different class, multi-class src  0.9
    Different class, single-class src 0.6
    Unknown source class              0.8
    Empty target_class (no filter)    1.0
    ================================  =====
    """
    # No phylogenetic target — use everything at full weight.
    if not target_class:
        return 1.0

    # Unknown source class — permissive fallback.
    if not source_class:
        return 0.8

    # Normalise for comparison.
    tc = target_class.strip().lower()
    sc = source_class.strip().lower()

    # Collect all contributing classes for breadth checks.
    contrib: set[str] = {sc}
    if source_classes_contrib:
        contrib.update(c.strip().lower() for c in source_classes_contrib if c)

    # Same class → start at 1.0, then apply order-level fine-tuning.
    if tc in contrib:
        if target_order and source_order:
            to = target_order.strip().lower()
            so = source_order.strip().lower()
            if to == so:
                return 1.0    # same class + same order
            else:
                return 0.8    # same class + different order
        return 1.0            # same class, no order filtering

    # Different class.  Check whether this marker set is multi-class
    # (conserved across classes) — that raises the weight.
    if len(contrib) >= 3:
        # Conserved across >=3 classes → strong biological signal.
        return 0.9
    elif len(contrib) >= 2:
        # Two classes → moderate conservation signal.
        return 0.8

    # Single distant class → significant penalty.
    return 0.6


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
        if type_key == "expert_rules" or type_key.startswith("_"):
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
                             species: Optional[str] = None,
                             target_class: str = "",
                             target_order: str = "",
                             adaptive_top_n: bool = False,
                             ) -> Dict[str, Score]:
    """Score one cluster against every cell type in the Knowledge Base.

    Parameters
    ----------
    kb : dict
        Either the raw KB (with ``markers``/``confirm``/``add``/``refine``
        structure) or a pre-built lookup from :func:`_build_kb_lookup`.
    cluster_markers : pd.DataFrame
        Top-N markers for a **single** cluster.  Must have columns
        ``names``, ``logfoldchanges``, ``pvals_adj``.  Non-protein-coding
        genes (lncRNA, MT-, RPL/RPS) are **automatically filtered out**
        before scoring.
    species : str or None
        If set, only cell types matching this species are scored.
    target_class : str
        Desired taxonomic class (纲), e.g. ``"Mammalia"``.  Empty string
        (default) disables phylogenetic weighting.
    target_order : str
        Desired taxonomic order (目), e.g. ``"Primates"``.  Empty string
        (default) disables order-level weighting.  Only has effect when
        *target_class* is also non-empty.
    adaptive_top_n : bool
        When ``True``, start with ``top_n=20`` but expand to 50 or 100
        if the KB-marker density in the top-20 is below 5%.  This helps
        on developmental / organoid data where known markers rank deep
        in the DE list.  Default ``False``.

    Returns
    -------
    Dict[str, Score]
        Mapping from ``type_key`` → :class:`Score` for each cell type.
    """
    # ── Pre-filter: keep only protein-coding / meaningful genes ─────
    from re import compile as _re_compile
    _RE_LNCRNA = _re_compile(
        r'^(LINC\d|AC\d|AL\d|AP\d|BX\d|FAM\d+[A-Z]|C\d+orf|'
        r'RP\d+-|CTC-|CTD-|RP11-|XXyac-|LLNLF-|WI2-|XXbac-)'
    )
    _RE_RIBO = _re_compile(r'^(RPL|RPS|MRPL|MRPS)\d*')
    # Filter: keep rows where names don't match noise patterns,
    # preserving logFC sort order (caller pre-sorts).
    _all_rows = cluster_markers.copy()
    _keep_mask = pd.Series(True, index=_all_rows.index)
    for _i, _row in _all_rows.iterrows():
        _g = str(_row["names"])
        if _g.startswith("MT-"):
            _keep_mask[_i] = False
        elif _RE_RIBO.match(_g):
            _keep_mask[_i] = False
        elif _RE_LNCRNA.match(_g):
            _keep_mask[_i] = False
    _filtered = _all_rows[_keep_mask]

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

    # ── Adaptive top-N ──────────────────────────────────────────────
    # Start at top_n=20.  If KB-marker density is <5%, expand to 50,
    # then 100, then 200.  This helps when known markers rank deep on
    # developmental / organoid data (KB self-audit Tier 0 finding).
    top_n = 20
    top_markers = _filtered.head(top_n).copy()
    top_gene_set = set(top_markers["names"].tolist())
    top_in_bg = top_gene_set & all_type_markers

    if adaptive_top_n:
        for _candidate_n in (50, 100, 200):
            if len(top_in_bg) / max(background_size, 1) >= 0.05:
                break
            _cand = _filtered.head(_candidate_n)
            _cand_set = set(_cand["names"].tolist())
            _cand_in_bg = _cand_set & all_type_markers
            if len(_cand_in_bg) > len(top_in_bg):
                top_n = _candidate_n
                top_markers = _cand
                top_gene_set = _cand_set
                top_in_bg = _cand_in_bg

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

        # ── 4b. Phylogenetic weight ─────────────────────────────────
        if target_class:
            # Read class/order/classes from the raw KB if available;
            # fall back to the lookup entry.
            type_data_raw = _type_data(kb, type_key)
            source_cls = type_data_raw.get("class", "") if type_data_raw else ""
            source_ord = type_data_raw.get("order", "") if type_data_raw else ""
            source_classes = type_data_raw.get("classes", []) if type_data_raw else []
            p_weight = phylogenetic_weight(
                source_cls, target_class,
                target_order=target_order,
                source_order=source_ord,
                source_classes_contrib=source_classes,
            )
            final_score *= p_weight

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
                          species: str,
                          target_class: str = "",
                          target_order: str = "",
                          adaptive_top_n: bool = False,
                          ) -> pd.DataFrame:
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
    target_class : str
        Taxonomic class for phylogenetic weighting (see
        :func:`score_cluster_against_kb`).
    target_order : str
        Taxonomic order for phylogenetic weighting.
    adaptive_top_n : bool
        Passed through to :func:`score_cluster_against_kb`.

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
            kb_all_markers, cl_sort, species=species,
            target_class=target_class, target_order=target_order,
            adaptive_top_n=adaptive_top_n,
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

_STRICTNESS_TEMPLATES: Dict[str, tuple] = {
    "strict":   (50,    0.01),
    "default":  (50,    0.05),
    "deep":     (200,   0.05),
    "wide":     (1000,  0.05),
    "relaxed":  (5000,  0.05),
    "manual":   (None,  None),
}


def resolve_expert_rule_params(
    strictness: str = "default",
    top_n: int = 0,
    pval_cutoff: float = 0.0,
) -> tuple[int, float]:
    """Resolve expert-rule constraints from a strictness template + overrides.

    Explicit *top_n* / *pval_cutoff* values take precedence over the
    template.  When ``strictness="manual"`` both must be set explicitly.

    Parameters
    ----------
    strictness : str
        One of ``"strict"``, ``"default"``, ``"deep"``, ``"wide"``,
        ``"relaxed"``, or ``"manual"``.
    top_n : int
        Explicit top-N override.  0 = use template value.
    pval_cutoff : float
        Explicit p-value override.  0.0 = use template value.

    Returns
    -------
    tuple[int, float]
        ``(resolved_top_n, resolved_pval_cutoff)``.

    Raises
    ------
    ValueError
        When ``strictness="manual"`` but *top_n* or *pval_cutoff* is unset.
    """
    template = _STRICTNESS_TEMPLATES.get(strictness)
    if template is None:
        logger.warning(
            "Unknown expert_rule_strictness '%s' — falling back to 'default'",
            strictness,
        )
        template = _STRICTNESS_TEMPLATES["default"]
    template_top_n, template_pval = template

    if strictness == "manual":
        if top_n <= 0 or pval_cutoff <= 0.0:
            raise ValueError(
                "expert_rule_strictness='manual' requires both "
                "expert_rule_top_n (>0) and expert_rule_pval_cutoff (>0.0)"
            )
        return top_n, pval_cutoff

    return (
        top_n if top_n > 0 else template_top_n,
        pval_cutoff if pval_cutoff > 0.0 else template_pval,
    )


def apply_expert_rules(kb: Dict[str, Any],
                       cluster_markers: pd.DataFrame,
                       top_n: int = 50,
                       pval_cutoff: float = 0.05,
                       ) -> tuple[Optional[str], list[Dict[str, Any]]]:
    """Try to deterministically match *cluster_markers* via expert rules.

    Iterates over ``kb["expert_rules"]`` (sorted by ``priority`` descending).
    All matching rules are collected; the highest-priority winner is returned
    as the first element.

    Only genes within the **top-N** DE genes and with ``pvals_adj <
    pval_cutoff`` are considered.  This prevents low-significance or
    deep-ranking genes from spuriously triggering rules (see KB self-audit
    Tier 0 findings).

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
        ``logfoldchanges``, and ideally ``pvals_adj``).
    top_n : int
        Only examine the top *top_n* DE genes.  Default 50.
    pval_cutoff : float
        Only consider genes with ``pvals_adj < pval_cutoff``.  Default 0.05.
        Silently ignored when the ``pvals_adj`` column is absent.

    Returns
    -------
    tuple[Optional[str], list[Dict[str, Any]]]
        ``(matched_action, all_matched_rules)``.

        *matched_action* — The winning rule's ``"action"`` key, or ``None``
        if no rule fired.
        *all_matched_rules* — Every rule that passed, in priority order
        (highest first).  Empty list when nothing matched.
    """
    rules = kb.get("expert_rules", [])
    if not rules:
        return None, []

    # ── Constrain to top-N statistically-significant DE genes ──────────
    de_subset = cluster_markers.head(top_n)
    if 'pvals_adj' in de_subset.columns:
        de_subset = de_subset[de_subset['pvals_adj'] < pval_cutoff]

    # Sort by priority descending (higher = more specific).
    sorted_rules = sorted(
        rules, key=lambda r: r.get("priority", 0), reverse=True
    )
    # Build a fast lookup: gene_name -> logfoldchanges.
    marker_map: Dict[str, float] = {}
    for _, row in de_subset.iterrows():
        gene_name = str(row["names"])
        marker_map[gene_name] = float(row["logfoldchanges"])

    cluster_genes = set(de_subset["names"].tolist())

    all_matched: list[Dict[str, Any]] = []

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
            all_matched.append(rule)

    if not all_matched:
        return None, []
    best = all_matched[0]
    return best.get("action"), all_matched


def detect_low_quality_cluster(cluster_markers: pd.DataFrame,
                                top_n: int = 20) -> tuple[bool, str]:
    """Detect clusters dominated by mitochondrial or ribosomal genes.

    Parameters
    ----------
    cluster_markers : pd.DataFrame
        DE gene DataFrame for a single cluster, sorted by logFC descending.
    top_n : int
        How many top genes to examine.  Default 20.

    Returns
    -------
    tuple[bool, str]
        ``(is_low_quality, reason)``.  *reason* is an empty string when
        the cluster is not flagged.
    """
    top = cluster_markers.head(top_n)
    genes = top["names"].tolist()

    n_mito = sum(1 for g in genes if str(g).startswith("MT-"))
    n_ribo = sum(
        1 for g in genes
        if str(g).startswith(("RPL", "RPS", "MRPL", "MRPS"))
    )

    if n_mito >= 3:
        return True, "mito_high ({} MT- genes in top-{})".format(n_mito, top_n)
    if n_ribo >= 5:
        return True, "ribo_high ({} ribosomal genes in top-{})".format(n_ribo, top_n)
    if n_mito + n_ribo >= 6:
        return True, "mito_ribo_mixed (MT={} RIBO={})".format(n_mito, n_ribo)

    return False, ""
