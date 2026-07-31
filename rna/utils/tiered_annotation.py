"""utils/tiered_annotation.py - Hierarchical (L1/L2/L3) cell-type label resolution.

Pure functions that sit *between* marker scoring and evidence fusion.
Given flat per-type ``Score`` results and the KB hierarchy, elects the best
label at the appropriate tier:

* **L1** (Broad_*) - never a ``cell_type``; only feeds ``cell_category``.
* **L2** (major type, e.g. ``RGC``) - default ``cell_type``.
* **L3** (subtype, e.g. ``RGC_Foxp2``) - elected as ``cell_subtype`` **only**
  when it passes all three gates (see :func:`resolve_tiered_label`).

The KB is intentionally asymmetric (RGC has 7 subtypes, Amacrine has none).
This module honours that asymmetry instead of forcing uniform granularity.
"""

from typing import Any

from core.annotation.scoring import Score

# ═══════════════════════════════════════════════════════════════════════
#  Thresholds (overridable via function args)
# ═══════════════════════════════════════════════════════════════════════

DEFAULT_MIN_SCORE = 0.25
DEFAULT_SUBTYPE_DELTA = 0.08
DEFAULT_CONSENSUS_FLOOR = "medium"

# Consensus level → comparable rank (higher = more reliable).
_CONSENSUS_RANK: dict[str, int] = {"gold": 4, "high": 3, "medium": 2, "low": 1}


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════


def _is_broad(key: str) -> bool:
    """L1 synthetic category keys start with ``Broad_``."""
    return key.startswith("Broad_")


def _build_subtype_map(hierarchy: dict[str, Any]) -> dict[str, str]:
    """``{subtype_key: l2_parent}`` from ``_hierarchy.categories[*].subtypes``."""
    mapping: dict[str, str] = {}
    for cat_def in hierarchy.get("categories", {}).values():
        for l2_key, l2_def in (cat_def.get("subtypes") or {}).items():
            for sub in l2_def.get("members", []):
                mapping[sub] = l2_key
    return mapping


def _classify_tier(type_key: str, subtype_map: dict[str, str]) -> str:
    """Classify a scored type_key as ``"L1"`` / ``"L2"`` / ``"L3"``."""
    if _is_broad(type_key):
        return "L1"
    if type_key in subtype_map:
        return "L3"
    return "L2"


def _private_markers_of(type_key: str, kb_lookup: dict[str, Any]) -> set[str]:
    return set(kb_lookup.get(type_key, {}).get("_private_markers", []))


def _consensus_levels_of(type_key: str, kb_lookup: dict[str, Any]) -> dict[str, str]:
    return kb_lookup.get(type_key, {}).get("consensus_levels", {})


def _consensus_of_hits(type_key: str, cluster_genes: set[str], kb_lookup: dict[str, Any]) -> str:
    """Best consensus level among this type's markers hit in the cluster."""
    consensus_levels = _consensus_levels_of(type_key, kb_lookup)
    hit = set(consensus_levels) & cluster_genes
    if not hit:
        return ""
    return max(
        (consensus_levels[g] for g in hit),
        key=lambda lv: _CONSENSUS_RANK.get(lv, 0),
        default="",
    )


def _best_consensus_among_private(
    type_key: str, cluster_genes: set[str], kb_lookup: dict[str, Any]
) -> str:
    """Best consensus among hit *private* markers (gate C input)."""
    private = _private_markers_of(type_key, kb_lookup)
    consensus_levels = _consensus_levels_of(type_key, kb_lookup)
    hit_private = private & cluster_genes
    if not hit_private:
        return ""
    return max(
        (consensus_levels.get(g, "low") for g in hit_private),
        key=lambda lv: _CONSENSUS_RANK.get(lv, 0),
        default="",
    )


def _build_evidence(
    type_key: str,
    tier: str,
    resolution: str,
    cluster_genes: set[str],
    kb_lookup: dict[str, Any],
    available_subtypes: list[str],
    subtype_candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    """Assemble the evidence metadata dict for the chosen label."""
    entry = kb_lookup.get(type_key, {})
    marker_weights = entry.get("marker_weights", {})
    # n_sources ≈ the most well-supported marker's source count for this type.
    n_sources = max(marker_weights.values(), default=0) if marker_weights else 0
    pm_hits = len(_private_markers_of(type_key, kb_lookup) & cluster_genes)
    return {
        "tier": tier,
        "subtype_resolution": resolution,
        "consensus": _consensus_of_hits(type_key, cluster_genes, kb_lookup),
        "n_sources": n_sources,
        "private_markers_hit": pm_hits,
        "available_subtypes": available_subtypes,
        "subtype_candidates": subtype_candidates,
    }


def format_subtype_candidates(candidates: list[dict[str, Any]]) -> str:
    """Render near-miss subtype candidates as ``name:score[gates]`` entries.

    Entries are ``; ``-joined with the score rounded to 2dp. Brackets are
    omitted entirely when ``failed_gates`` is empty (a winning subtype renders
    e.g. ``RGC_Foxp2:0.58``). An empty input list renders the empty string.
    """
    if not candidates:
        return ""
    parts: list[str] = []
    for cand in candidates:
        name = cand["type"]
        score = f"{cand['score']:.2f}"
        gates = cand.get("failed_gates", [])
        if gates:
            parts.append(f"{name}:{score}[{','.join(gates)}]")
        else:
            parts.append(f"{name}:{score}")
    return "; ".join(parts)


# ═══════════════════════════════════════════════════════════════════════
#  Core decision
# ═══════════════════════════════════════════════════════════════════════


def resolve_tiered_label(
    cluster_scores: dict[str, Score],
    hierarchy: dict[str, Any],
    kb_lookup: dict[str, Any],
    cluster_top_genes: list[str],
    min_score: float = DEFAULT_MIN_SCORE,
    subtype_delta: float = DEFAULT_SUBTYPE_DELTA,
    consensus_floor: str = DEFAULT_CONSENSUS_FLOOR,
) -> tuple[str, dict[str, Any]]:
    """Elect the best tiered label for a cluster.

    A subtype (L3) wins over its parent L2 **only when** it passes every gate:

    * **Gate A (score proximity)** - ``subtype.score >= best_L2.score - subtype_delta``
      *and* ``subtype.score >= min_score``. Rejects noise spikes.
    * **Gate B (private marker)** - at least one of the subtype's
      ``_private_markers`` is found in ``cluster_top_genes``. Rejects L2 signal
      bleeding through shared markers.
    * **Gate C (consensus)** - the best consensus level among the hit private
      markers is ``>= consensus_floor``. Rejects single-source markers.

    When no subtype passes, the L2 label is returned with
    ``subtype_resolution="unresolved"``. When the L2 has no subtypes defined in
    the KB, ``subtype_resolution="na"``.

    Parameters
    ----------
    cluster_scores : dict
        ``{type_key: Score}`` from ``score_cluster_against_kb()`` — **including**
        ``Broad_*`` entries (used only for tier classification, never as a label).
    hierarchy : dict
        ``kb["_hierarchy"]`` carrying ``categories[*].subtypes``.
    kb_lookup : dict
        ``{type_key: {...}}`` from ``_build_kb_lookup()``; must carry
        ``_private_markers``, ``consensus_levels``, and ``marker_weights``.
    cluster_top_genes : list
        The cluster's top-N DE genes (for private-marker gate B).
    min_score, subtype_delta, consensus_floor : float/str
        Gate thresholds (see module defaults).

    Returns
    -------
    (label, evidence) : tuple
        *label* is an L2 type (always) or an L3 subtype (when resolved).
        *evidence* carries ``tier``, ``subtype_resolution``, ``consensus``,
        ``n_sources``, ``private_markers_hit``, ``available_subtypes``,
        ``subtype_candidates``.
        Empty *label* signals the caller to mark the cluster Unknown.
    """
    cluster_genes = set(cluster_top_genes)
    subtype_map = _build_subtype_map(hierarchy)
    floor_rank = _CONSENSUS_RANK.get(consensus_floor, _CONSENSUS_RANK["medium"])

    # ── 1. Partition scored types by tier ───────────────────────────────
    l2_scores: dict[str, Score] = {}
    l3_scores_by_parent: dict[str, dict[str, Score]] = {}
    for key, sc in cluster_scores.items():
        tier = _classify_tier(key, subtype_map)
        if tier == "L2":
            l2_scores[key] = sc
        elif tier == "L3":
            l3_scores_by_parent.setdefault(subtype_map[key], {})[key] = sc
        # L1 (Broad_*) ignored for cell_type election.

    # ── 2. No L2 candidate → cannot resolve (caller marks Unknown) ──────
    if not l2_scores:
        return "", {
            "tier": "",
            "subtype_resolution": "na",
            "consensus": "",
            "n_sources": 0,
            "private_markers_hit": 0,
            "available_subtypes": [],
            "subtype_candidates": [],
        }

    # ── 3. Best L2 by score ─────────────────────────────────────────────
    best_l2_key = max(l2_scores, key=lambda k: l2_scores[k].score)
    best_l2_score = l2_scores[best_l2_key].score

    # Subtypes defined in the KB for this L2 (may be empty).
    available = [k for k, parent in subtype_map.items() if parent == best_l2_key]

    if not available:
        # L2 has no subtypes in the KB — "na" (no capacity to subdivide).
        return best_l2_key, _build_evidence(
            best_l2_key, "L2", "na", cluster_genes, kb_lookup, available, []
        )

    # ── 4. Gate each L3 subtype of best L2 ──────────────────────────────
    # Each subtype's gates A/B/C are evaluated INDEPENDENTLY and every failed
    # gate is recorded; an empty ``failed_gates`` list means all three passed.
    best_subtype: str | None = None
    best_subtype_score = -1.0
    subtype_candidates: list[dict[str, Any]] = []
    for sub_key, sub_sc in l3_scores_by_parent.get(best_l2_key, {}).items():
        failed_gates: list[str] = []
        # Gate A: score proximity + floor.
        if sub_sc.score < best_l2_score - subtype_delta or sub_sc.score < min_score:
            failed_gates.append("A")
        # Gate B: at least one private marker hit.
        private_markers = _private_markers_of(sub_key, kb_lookup)
        pm_hits = len(private_markers & cluster_genes)
        if not pm_hits:
            failed_gates.append("B")
        # Gate C: best consensus among hit private markers >= floor.
        pm_consensus = _best_consensus_among_private(sub_key, cluster_genes, kb_lookup)
        if _CONSENSUS_RANK.get(pm_consensus, 0) < floor_rank:
            failed_gates.append("C")
        subtype_candidates.append(
            {
                "type": sub_key,
                "score": sub_sc.score,
                "failed_gates": failed_gates,
                "private_markers_hit": pm_hits,
            }
        )
        # All gates passed — pick the highest-scoring qualifying subtype.
        if not failed_gates and sub_sc.score > best_subtype_score:
            best_subtype = sub_key
            best_subtype_score = sub_sc.score
    subtype_candidates.sort(key=lambda c: c["score"], reverse=True)

    # ── 5. Emit decision ────────────────────────────────────────────────
    if best_subtype is not None:
        return best_subtype, _build_evidence(
            best_subtype,
            "L3",
            "resolved",
            cluster_genes,
            kb_lookup,
            available,
            subtype_candidates,
        )
    # No subtype passed all gates — L2 unresolved.
    return best_l2_key, _build_evidence(
        best_l2_key,
        "L2",
        "unresolved",
        cluster_genes,
        kb_lookup,
        available,
        subtype_candidates,
    )
