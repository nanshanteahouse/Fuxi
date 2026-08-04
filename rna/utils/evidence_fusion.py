"""
utils/evidence_fusion.py — Evidence fusion decision engine.

Combines marker scoring results (:mod:`utils.marker_scoring`), expert rules,
and optional AI suggestions to produce a final cell type annotation decision
for each cluster.

Decision priority (hard-coded tiers):
    1. Expert rule (if triggered)
    2. Marker score >= 0.7       →  high confidence
    3. Marker score 0.5–0.7     →  medium confidence (AI-aware)
    4. Marker score 0.25–0.5    →  low confidence (AI-aware)
    5. All else                 →  Unknown
"""

import json
from dataclasses import dataclass
from typing import NamedTuple, Optional

import pandas as pd

from core.annotation.potency import (
    KADPConfig,
    compute_potency,
    derive_developmental_poles,
    evaluate_passes,
)


class DiagnosticInfo(NamedTuple):
    """Diagnostic context for Uncertain/Unknown clusters (v3.1.0+).

    Attributes
    ----------
    category : str
        One of ``'no_kb_match'`` | ``'low_quality_data'`` | ``'ambiguous'`` |
        ``'weak_signal'`` | ``'true_unknown'``.
    top_competitors : list
        Top-3 ``(cell_type, score)`` competitors, if any.
    detail : str
        Human-readable diagnostic detail.
    """

    category: str
    top_competitors: list
    detail: str


class FusionDecision(NamedTuple):
    """Final annotation decision for one cluster.

    Attributes
    ----------
    cell_type : str
        Final cell type name (from KB).
    confidence : str
        ``'high'`` | ``'medium'`` | ``'low'`` | ``'unknown'`` | ``'rule'`` | ``'transition'``.
    score : float
        The score that led to this decision.
    method : str
        Which tier produced the decision.
    n_markers_found : int
        How many KB markers matched in the cluster's top-20.
    ai_agreed : bool
        Did AI agree with the marker-based decision?
    ai_suggested : str
        What AI suggested (if AI was called).
    explanation : str
        Human-readable explanation.
    alternative_rules : list
        Other expert rules that also matched this cluster (if expert_rule
        was the winning method).  Empty list otherwise.
    diagnostic : DiagnosticInfo or None
        Diagnostic context for Unknown/Uncertain clusters (v3.1.0+).
        ``None`` for all non-Unknown decisions.
    cell_category : str
        Broad lineage category (e.g. ``"Broad_Neuron"``).
    tier : str
        Hierarchy tier of the chosen label: ``"L2"`` (major type) or
        ``"L3"`` (subtype). ``"L1"`` is never a final cell_type.
    consensus : str
        Best consensus level among hit markers for the chosen type.
    n_sources : int
        Number of KB sources supporting the chosen type.
    subtype_resolution : str
        ``"resolved"`` (L3 subtype elected), ``"unresolved"`` (subtypes
        existed in KB but none passed the gates), or ``"na"`` (KB has no
        subtypes for this L2 type).
    review_reason : str
        Non-empty when this decision should be surfaced in the engine-side
        review queue.  Currently set to ``"single_marker_rule"`` by the D5
        arbitration for uncorroborated single-marker expert-rule hits that
        carried no strong marker-scoring competitor.
    potency : Optional[dict] or None
        KADP developmental-potency payload (layer 3), set when this cluster
        was named by the developmental-potency path.  Always a three-value
        dict {'ratio': float, 'abs': float, 'gap': float} — never a bare
        float.  None for every non-KADP decision.
    source_votes : Optional[dict] or None
        Layer-4 METC multi-source vote payload (marker/expert/ai/celltypist
        votes), set when METC arbitrated this cluster.  METC assignments
        always carry a fresh dict; None otherwise.
    """

    cell_type: str
    confidence: str
    score: float
    method: str
    n_markers_found: int
    ai_agreed: bool
    ai_suggested: str
    explanation: str
    alternative_rules: list
    diagnostic: Optional[DiagnosticInfo] = None
    cell_category: str = ""
    tier: str = ""
    consensus: str = ""
    n_sources: int = 0
    subtype_resolution: str = ""
    review_reason: str = ""
    # -- Append-only tail fields (plan todo 1a).  Optional[dict] = None is
    #    deliberate: never a shared mutable default, and the potency payload
    #    must keep all three values (ratio/abs/gap) intact — no bare float.
    potency: Optional[dict] = None
    source_votes: Optional[dict] = None


@dataclass
class METCConfig:
    """Layer-4 multi-evidence transition consensus (METC) config (todo 1c).

    Controls the METC arbitration branch in :func:`fuse_evidence`.  Every
    field defaults to *disabled* so existing callers keep byte-identical
    behavior unless METC is explicitly enabled.
    """

    enabled: bool = False
    min_sources: int = 3
    min_distinct_transition: int = 3


def _is_transition_state(
    marker_scores: dict,
    kb: dict,
    delta_threshold: float = 0.15,
    min_score: float = 0.25,
    incompatible_transitions: Optional[list] = None,
) -> Optional[tuple]:
    """Detect transition between two types of same lineage.

    When Top-2 Fisher scores differ by < delta_threshold AND share
    the same parent (broad category), return (top1_type, top2_type).
    Otherwise return None.

    Parameters
    ----------
    marker_scores : dict
        {type_key: Score or float} with Broad_* keys already stripped.
    kb : dict
        Full KB dict with parent and marker information.
    delta_threshold : float
        Maximum score difference for transition consideration (default 0.15).
    min_score : float
        Minimum absolute score floor for the top candidate (default 0.25).
    incompatible_transitions : list or None
        List of [type_a, type_b] pairs that are explicitly forbidden as
        transitions despite meeting all other criteria.

    Returns
    -------
    tuple or None
        (top1_type, top2_type) if a valid transition is detected, else None.

    Notes
    -----
    CRITICAL: This function receives marker_scores that has ALREADY had
    Broad_* keys stripped (by the caller upstream). It must NOT operate on
    raw scores containing Broad_* entries.

    Two additional gates are applied after passing the basic checks:

    **Incompatible pairs gate:** If *incompatible_transitions* is provided,
    each forbidden pair is checked order-insensitively against the
    candidate (top1, top2). Any match returns None.

    **Private marker overlap gate (P0.4):** When both types share the same
    Broad_* parent, the parent\'s markers are subtracted from each type\'s
    marker set. If the remaining (type-specific) markers share no overlap,
    the Fisher proximity is attributed to shared Broad_* markers alone
    and the transition is rejected. When neither type has any marker
    beyond the parent\'s set, the gate falls back to the KB
    ``_private_markers`` of each type: the transition is accepted only
    when the two types share at least one private marker (both missing
    or disjoint → rejected).
    """
    if len(marker_scores) < 2:
        return None
    sorted_types = sorted(
        marker_scores.items(), key=lambda x: _resolve_score(marker_scores, x[0])[0], reverse=True
    )
    top1_key, _ = sorted_types[0]
    top2_key, _ = sorted_types[1]
    score1, _ = _resolve_score(marker_scores, top1_key)
    score2, _ = _resolve_score(marker_scores, top2_key)
    if score1 - score2 >= delta_threshold:
        return None
    # Absolute score floor: top score must be >= min_score (same cutoff that
    # separates marker_scoring_low from unknown in DECISION_TIERS).
    if score1 < min_score:
        return None
    parent1 = kb.get(top1_key, {}).get("parent", "")
    parent2 = kb.get(top2_key, {}).get("parent", "")
    if not parent1 or not parent2 or parent1 != parent2:
        return None

    # Gate: incompatible transition pairs
    if incompatible_transitions:
        for pair in incompatible_transitions:
            if {top1_key, top2_key} == set(pair):
                return None

    # Gate: type-specific markers must overlap (not just Broad_* shared markers)
    if parent1 and parent2 and parent1 == parent2:
        parent_markers = set(kb.get(parent1, {}).get("markers", {}).get("confirm", {}).keys())
        t1_confirm = set(kb.get(top1_key, {}).get("markers", {}).get("confirm", {}).keys())
        t2_confirm = set(kb.get(top2_key, {}).get("markers", {}).get("confirm", {}).keys())
        t1_add = set(kb.get(top1_key, {}).get("markers", {}).get("add", {}).keys())
        t2_add = set(kb.get(top2_key, {}).get("markers", {}).get("add", {}).keys())
        t1_specific = (t1_confirm | t1_add) - parent_markers
        t2_specific = (t2_confirm | t2_add) - parent_markers
        if t1_specific or t2_specific:
            specific_overlap = t1_specific & t2_specific
            if not specific_overlap:
                return None
        else:
            # P0.4: neither type has markers beyond the shared parent's set.
            # The Fisher proximity is then attributable to shared Broad_*
            # markers alone — fall back to each type's KB `_private_markers`
            # and reject unless the two types share at least one.
            private1 = set(kb.get(top1_key, {}).get("_private_markers", []))
            private2 = set(kb.get(top2_key, {}).get("_private_markers", []))
            if not (private1 & private2):
                return None

    return (top1_key, top2_key)


# Decision priority tiers — evaluated in order.
# Each tier is a (name, callable) where
# callable(score, expert_rule_result, ai_suggestion) → bool.
DECISION_TIERS = [
    ("expert_rule", lambda s, e, a: e is not None),  # Tier 0
    ("transition_state", lambda s, e, a: False),  # Tier 1 (detected pre-loop)
    ("marker_scoring_high", lambda s, e, a: s >= 0.7),  # Tier 1
    ("marker_scoring_medium", lambda s, e, a: 0.5 <= s < 0.7),  # Tier 2
    ("marker_scoring_low", lambda s, e, a: 0.25 <= s < 0.5),  # Tier 3
    ("unknown", lambda s, e, a: True),  # Tier 4
]

_CONFIDENCE_MAP = {
    "expert_rule": "rule",
    "marker_scoring_high": "high",
    "marker_scoring_medium": "medium",
    "marker_scoring_low": "low",
    "unknown": "unknown",
    "transition_state": "transition",
}


# WEAK_EVIDENCE — evidence-strength labels that must never produce a
# confident label (plan D1/D2).  A marker-scoring winner carrying one of
# these is capped at ``confidence="low"`` regardless of the tier mapping,
# and its cluster is flagged for the engine-side review queue via
# ``FusionDecision.review_reason = evidence_type``.
WEAK_EVIDENCE = frozenset({"single_marker", "window_padding", "weak_multi", "zero_evidence"})


# ── Label normalisation ──────────────────────────────────────────────
# AI-generated labels and KB cell-type keys often differ only by
# whitespace vs underscores (e.g. "Amacrine Cell" vs "Amacrine_Cell").
# Normalise both sides to a canonical form before comparison so that
# ai_agreed reflects genuine biological disagreement rather than
# formatting differences.


def _normalise_label(label: Optional[str]) -> str:
    """Canonicalise a cell-type label for fuzzy comparison.

    Collapses runs of non-alphanumeric characters into a single
    underscore, lowercases, and strips leading/trailing underscores.
    Non-ASCII characters (e.g. ``ü`` in ``Müller``) are decomposed
    to their ASCII base form via NFKD normalisation.
    """
    if not label:
        return ""
    import re
    import unicodedata

    # NFKD decomposes accents / umlauts: ü -> u + combining diaeresis
    nfkd = unicodedata.normalize("NFKD", label)
    # Drop combining chars and other non-ASCII
    ascii_label = nfkd.encode("ascii", "ignore").decode("ascii")
    # Replace all non-alphanumeric runs with a single underscore.
    normalised = re.sub(r"[^a-zA-Z0-9]+", "_", ascii_label)
    return normalised.strip("_").lower()


def _labels_match(a: Optional[str], b: Optional[str]) -> bool:
    """Return ``True`` if *a* and *b* refer to the same cell type.

    Uses :func:`_normalise_label` to ignore whitespace/underscore/
    hyphen/punctuation differences.
    """
    return _normalise_label(a) == _normalise_label(b)


# ═══════════════════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════════════════


def _resolve_score(marker_scores: dict, type_key: str) -> tuple:
    """Extract score and n_markers_found for *type_key* from *marker_scores*.

    Handles both :class:`~utils.marker_scoring.Score` objects and bare floats
    (useful for simplified contexts / testing).
    """
    entry = marker_scores.get(type_key)
    if entry is None:
        return 0.0, 0
    if isinstance(entry, (int, float)):
        return float(entry), 0
    return float(entry.score), int(entry.n_markers_found)


def _find_best_type(marker_scores: dict) -> tuple:
    """Return ``(best_type, best_score, n_markers_found)``.

    Returns ``(None, 0.0, 0)`` when *marker_scores* is empty.
    """
    if not marker_scores:
        return None, 0.0, 0

    best_type = max(marker_scores, key=lambda k: _resolve_score(marker_scores, k)[0])
    best_score, n_markers = _resolve_score(marker_scores, best_type)
    return best_type, best_score, n_markers


def _winning_rule_entry(
    alternative_rules: Optional[list],
    expert_rule_result: Optional[str],
) -> Optional[dict]:
    """Locate the winning rule's metadata entry among *alternative_rules*.

    ``apply_expert_rules`` returns the priority-sorted list of matched rules
    as its second element with the winner first, so the first entry whose
    ``"action"`` equals *expert_rule_result* is the winner; when no match
    is found the first entry is returned as a fallback.  ``None`` when
    *alternative_rules* is empty (callers that omit the metadata are treated
    as legacy: every rule hit is corroborated).
    """
    if not alternative_rules:
        return None
    for rule in alternative_rules:
        if rule.get("action") == expert_rule_result:
            return rule
    return alternative_rules[0]


def _rule_is_single_marker(rule: Optional[dict]) -> bool:
    """True when a rule requires exactly one ``markers_present`` gene.

    The D5 arbitration applies only to *single-marker* rules; multi-marker
    rules (two or more independent ``markers_present`` genes) keep the legacy
    full-confidence behavior even when uncorroborated.
    """
    if not rule:
        return False
    condition = rule.get("condition", {}) or {}
    return len(condition.get("markers_present", {}) or {}) == 1


def _has_strong_multi_marker_competitor(marker_scores: dict) -> bool:
    """Whether a top-2 marker-scoring type carries strong multi-marker evidence.

    D5 competitor check: returns ``True`` when any of the two highest-scoring
    types has ``evidence_type == "multi_marker"`` with ``consensus`` in
    ``("gold", "high")``.  Bare-float score entries (no :class:`~utils.marker_scoring.Score`
    metadata) never qualify, preserving simplified test contexts.
    """
    if not marker_scores:
        return False
    ranked = sorted(
        marker_scores.items(),
        key=lambda kv: _resolve_score(marker_scores, kv[0])[0],
        reverse=True,
    )
    for type_key, _ in ranked[:2]:
        entry = marker_scores[type_key]
        if isinstance(entry, (int, float)):
            continue
        if getattr(entry, "evidence_type", "") == "multi_marker" and getattr(
            entry, "consensus", ""
        ) in ("gold", "high"):
            return True
    return False


def _explain(
    cell_type: str,
    method: str,
    score: float,
    n_markers: int,
    best_type: Optional[str],
    ai_suggestion: Optional[str],
    ai_agreed: bool,
    confidence: str = "",
    alternative_rules: Optional[list] = None,
) -> str:
    """Build a human-readable explanation."""
    if method == "expert_rule":
        parts = [f"Expert rule matched: {cell_type}"]
        if alternative_rules and len(alternative_rules) > 1:
            alt_names = [r.get("action") for r in alternative_rules[1:]]
            parts.append(f"(also matched rules: {', '.join(alt_names)})")
        if score > 0:
            parts.append(f"(marker score: {score:.3f})")
    elif method == "unknown":
        parts = ["No cell type could be confidently assigned"]
        if best_type and score > 0:
            parts.append(f"(best match: {best_type}, score: {score:.3f})")
    elif method == "transition_state":
        alt_text = "\n".join(alternative_rules) if alternative_rules else ""
        return (
            f"Transitional state detected: {cell_type}\n"
            f"  Confidence: {confidence}  Score: {score:.4f}\n"
            f"  {alt_text}"
        )

    else:
        parts = [
            f"Marker scoring selected {cell_type} with score {score:.3f}",
        ]
        if n_markers > 0:
            parts.append(f"({n_markers} KB markers found in cluster top-20)")

    if ai_suggestion:
        if ai_agreed:
            parts.append("\u2014 AI agreed with this assignment")
        else:
            parts.append(
                f"\u2014 AI suggested '{ai_suggestion}' (different from marker-based result)"
            )

    return " ".join(parts)


def _build_diagnostic_summary(decisions: list) -> dict:
    """Build a diagnostic category count summary from fusion decisions."""
    summary = {}
    for d in decisions:
        if d.diagnostic and d.diagnostic.category:
            cat = d.diagnostic.category
            summary[cat] = summary.get(cat, 0) + 1
    return summary


def _classify_unknown(
    marker_scores: dict,
    low_quality_reason: str = "",
) -> DiagnosticInfo:
    """Classify an Unknown/Uncertain decision into a diagnostic category.

    Parameters
    ----------
    marker_scores : dict
        ``{type_key: Score or float}`` from marker scoring.
    low_quality_reason : str
        Non-empty if the cluster was flagged by
        :func:`~utils.marker_scoring.detect_low_quality_cluster`.

    Returns
    -------
    DiagnosticInfo
    """
    scored = [(k, _resolve_score(marker_scores, k)) for k in marker_scores]
    scored.sort(key=lambda x: -x[1][0])

    top3 = [(t, round(s, 4)) for t, (s, _) in scored[:3] if s > 0]

    if low_quality_reason:
        return DiagnosticInfo(
            category="low_quality_data",
            top_competitors=top3,
            detail=(f"Cluster flagged as low-quality: {low_quality_reason}"),
        )

    if scored and scored[0][1][0] >= 0.25:
        ambiguous_candidates = [(t, round(s, 4)) for t, (s, _) in scored if s >= 0.25]
        if len(ambiguous_candidates) >= 2:
            names = ", ".join(t for t, _ in ambiguous_candidates[:5])
            return DiagnosticInfo(
                category="ambiguous",
                top_competitors=top3,
                detail=(f"Multiple cell types with score >= 0.25: {names}"),
            )

    if scored and 0 < scored[0][1][0] < 0.25:
        return DiagnosticInfo(
            category="weak_signal",
            top_competitors=top3,
            detail=(
                f"Best score {scored[0][1][0]:.4f} below 0.25 threshold "
                f"(best type: {scored[0][0]})"
            ),
        )

    if not any(s > 0 for _, (s, _) in scored):
        return DiagnosticInfo(
            category="no_kb_match",
            top_competitors=[],
            detail="No KB cell type had any marker overlap with this cluster.",
        )

    return DiagnosticInfo(
        category="true_unknown",
        top_competitors=top3,
        detail="Could not determine cell type by any method.",
    )


# ═══════════════════════════════════════════════════════════════════════
def _kadp_name_candidate(
    marker_scores: dict,
    kb: dict,
    candidate: "FusionDecision",
    cfg: KADPConfig,
) -> Optional["FusionDecision"]:
    """Name an ambiguous multi-peak candidate via developmental potency (todo 4).

    When the Progenitor pole dominates the marker scores, the ambiguous
    candidate is replaced by a ``developmental_potency`` decision naming the
    argmax progenitor type.  The replacement is built with ``._replace`` so
    every other field (score / n_markers_found / ai_* / alternative_rules) is
    inherited from the candidate -- a fresh FusionDecision is never constructed.

    Returns the named decision, or ``None`` when potency does not pass -- the
    caller then returns the candidate unchanged (exit semantics: a KADP miss
    never falls into the tier loop).
    """
    result = compute_potency(
        marker_scores,
        derive_developmental_poles(kb),
        cfg,
    )
    if not evaluate_passes(result, cfg):
        return None
    best = result.best_progenitor_type
    if best is None:  # naming precondition: max_prog > 0 (defensive)
        return None
    potency = result.to_potency_dict()
    # The tied multi-peak list doubles as the competitor list -- it already
    # carries the {'cell_type', 'score'} dict shape the corpus uses.
    competitors = list(candidate.diagnostic.top_competitors) if candidate.diagnostic else []
    return candidate._replace(
        cell_type=best,
        confidence="medium",
        method="developmental_potency",
        explanation=(
            f"Developmental potency named '{best}' as differentiating precursor "
            f"(max progenitor {result.max_prog:.3f} vs max terminal "
            f"{result.max_term:.3f}, ratio {result.ratio:.3f}, "
            f"gap {result.gap:.3f})"
        ),
        diagnostic=DiagnosticInfo(
            category="developmental_potency",
            top_competitors=competitors,
            detail=(
                f"Developmental potency naming -- max progenitor "
                f"{result.max_prog:.3f} vs max terminal {result.max_term:.3f}; "
                f"{json.dumps(potency)}"
            ),
        ),
        review_reason="kadp_precursor",
        potency=potency,
        source_votes=None,
    )


#  Public API
# ═══════════════════════════════════════════════════════════════════════


def fuse_evidence(
    marker_scores: dict,
    expert_rule_result: Optional[str],
    kb: Optional[dict] = None,
    cluster_markers: Optional[pd.DataFrame] = None,
    ai_suggestion: Optional[str] = None,
    alternative_rules: Optional[list] = None,
    low_quality_reason: str = "",
    unconstrained: bool = False,
    allows_transitions: bool = False,
    incompatible_transitions: Optional[list] = None,
    multi_peak_min_types: int = 3,
    multi_peak_score_floor: float = 0.9,
    kadp_cfg: Optional[KADPConfig] = None,
    metc_cfg: Optional[METCConfig] = None,
    celltypist_suggestion: Optional[str] = None,
) -> "FusionDecision":
    """Combine marker scores, expert rules, and AI into one decision.

    Parameters
    ----------
    marker_scores : dict
        Output of :func:`utils.marker_scoring.score_cluster_against_kb`.
        Maps ``type_key → Score`` (or bare ``float`` in simplified contexts).
    expert_rule_result : str or None
        Output of :func:`utils.marker_expert_rules.apply_expert_rules`.
    kb : dict or None
        Full KB dict (reserved for explanation enrichment).
    cluster_markers : pd.DataFrame or None
        Marker DataFrame for this cluster (reserved for future use).
    ai_suggestion : str or None
        AI-proposed cell type, if available.
    alternative_rules : list or None
        Other expert rules that also matched (from
        :func:`apply_expert_rules`' second return element).
    low_quality_reason : str
        Non-empty if the cluster was flagged by
        :func:`~utils.marker_scoring.detect_low_quality_cluster` (v3.1.0+).
    allows_transitions : bool
        When ``True``, enables transition-state detection between closely
        related cell types (e.g. developing/tissue-remodelling contexts).
    incompatible_transitions : list or None
        List of [type_a, type_b] pairs explicitly forbidden as
        developmental transitions.  Passed through to
        :func:`_is_transition_state`.
    multi_peak_min_types : int
        Minimum number of types with score >= *multi_peak_score_floor*
        before the transition path downgrades to ``ambiguous`` (default 3).
    multi_peak_score_floor : float
        Score floor for counting a tied type in the multi-peak check
        (default 0.9).
    kadp_cfg : KADPConfig or None
        Layer-3 KADP developmental-potency config (todo 4).  When enabled and
        a multi-peak ``ambiguous`` candidate carries progenitor-dominant
        potency, the candidate is named as its argmax precursor type
        (``method="developmental_potency"``).  Default ``None`` keeps the
        baseline candidate-return behavior byte-identical.
    metc_cfg : METCConfig or None
        Layer-4 METC multi-source transition consensus config (todo 9).
        Reserved in this todo: enables the candidate exit gate so METC
        arbitration can slot in later; no arbitration happens yet.
    celltypist_suggestion : str or None
        CellTypist label for this cluster, pre-harmonized engine-side (todo 8).
        Passed through for METC arbitration (todo 9) — not consumed here.

    Returns
    -------
    FusionDecision

    Notes
    -----
    **D5 weak-rule arbitration.** A corroborated expert-rule hit keeps the
    legacy full-confidence early return (``confidence="rule"``).  An
    *uncorroborated single-marker* rule hit (``corroborated=False`` and
    exactly one ``markers_present`` gene) does not early-return: when a
    top-2 marker-scoring type carries ``evidence_type="multi_marker"`` with
    ``consensus`` in ``("gold", "high")`` the rule yields to marker scoring
    (its label is recorded in the explanation as "weak rule overridden");
    otherwise it returns ``confidence="low"`` with
    ``review_reason="single_marker_rule"`` for engine-side review-queue
    collection.

    **D2 weak-evidence cap.** A marker-scoring winner (tiers 1-3) whose
    ``evidence_type`` is in ``WEAK_EVIDENCE`` (single-marker /
    window-padding / weak-multi / zero-evidence) or is ``"ai_only"`` is capped
    at ``confidence="low"`` regardless of the tier mapping, and
    ``review_reason=evidence_type`` flags the cluster for the engine-side
    review queue.  ``ai_only`` is normally set engine-side (task 9); this cap
    only guards against a tier-level Score that already carries it.

    **M1 boundary.** In developmental RPC contexts the "correct" label is
    a marker-context-conditional question owned by **layer 3**; layers 1+2
    only guarantee that weak evidence never produces a confident label.
    """

    # ── Tier 0: expert rule (highest priority) ─────────────────────────
    # D5 weak-rule arbitration: a *corroborated* rule hit keeps the legacy
    # full-confidence early return; an *uncorroborated single-marker* rule
    # hit does not — it either yields to a strong marker-scoring competitor
    # or returns at low confidence flagged for engine-side review.
    #
    # Boundary (M1): in developmental RPC contexts the "correct" label is a
    # context-conditional question owned by layer 3.  Layers 1+2 only
    # guarantee that weak evidence never produces a confident label.
    weak_rule_override_label: Optional[str] = None
    if expert_rule_result is not None:
        rule_score, rule_n = _resolve_score(marker_scores, expert_rule_result)
        ai_agreed = _labels_match(ai_suggestion, expert_rule_result) if ai_suggestion else False

        # D5: read the winning rule's corroboration metadata (populated by
        # apply_expert_rules on its second return element).
        winning_rule = _winning_rule_entry(alternative_rules, expert_rule_result)
        corroborated = bool(winning_rule.get("corroborated", True)) if winning_rule else True
        uncorroborated_single = (not corroborated) and _rule_is_single_marker(winning_rule)

        if uncorroborated_single and _has_strong_multi_marker_competitor(marker_scores):
            # (a) Yield: a strong multi-marker candidate wins.  Run the
            # marker-scoring path; the rule label is recorded in the
            # explanation and the cluster is NOT labelled by the rule.
            weak_rule_override_label = expert_rule_result
            expert_rule_result = None
        elif uncorroborated_single:
            # (b) No strong competitor: low-confidence rule label + review.
            score_note = f" (marker score: {rule_score:.3f})" if rule_score > 0 else ""
            return FusionDecision(
                cell_type=expert_rule_result,
                confidence="low",
                score=rule_score,
                method="expert_rule",
                n_markers_found=rule_n,
                ai_agreed=ai_agreed,
                ai_suggested=ai_suggestion or "",
                explanation=(
                    f"Expert rule matched '{expert_rule_result}' but its "
                    f"corroborators were absent from the DE top-N subset. "
                    f"Uncorroborated single-marker rule — labelled at low "
                    f"confidence.{score_note}"
                ),
                alternative_rules=alternative_rules or [],
                review_reason="single_marker_rule",
            )
        else:
            # Corroborated hit (or non-single-marker rule): legacy behavior.
            # Quality gate (v3.1.0+): if Fisher scoring completely disagrees
            # with the expert rule (zero KB marker overlap), downgrade confidence
            # from 'rule' to 'low'.  This prevents noise-triggered rules (e.g.
            # a gene buried at rank 4000 in relaxed mode) from outranking well-
            # scored Fisher matches in downstream analysis.
            if rule_score < 0.25 and rule_n == 0:
                conf = "low"
                warning_note = (
                    f"Expert rule matched '{expert_rule_result}' but independent "
                    f"marker scoring found zero KB marker overlap (score={rule_score:.3f}, "
                    f"n_markers=0). Downgrading confidence from 'rule' to 'low'."
                )
            else:
                conf = "rule"
                warning_note = ""

            explanation_parts = []
            if warning_note:
                explanation_parts.append(warning_note)
            explanation_parts.append(
                _explain(
                    expert_rule_result,
                    "expert_rule",
                    rule_score,
                    rule_n,
                    expert_rule_result,
                    ai_suggestion,
                    ai_agreed,
                    alternative_rules=alternative_rules,
                )
            )

            return FusionDecision(
                cell_type=expert_rule_result,
                confidence=conf,
                score=rule_score,
                method="expert_rule",
                n_markers_found=rule_n,
                ai_agreed=ai_agreed,
                ai_suggested=ai_suggestion or "",
                explanation=" | ".join(explanation_parts),
                alternative_rules=alternative_rules or [],
            )

    # ── Unconstrained AI mode: accept AI suggestion directly ──────────
    if (
        unconstrained
        and ai_suggestion
        and (
            not marker_scores
            or max((_resolve_score(marker_scores, k)[0] for k in marker_scores), default=0) < 0.25
        )
    ):
        return FusionDecision(
            cell_type=ai_suggestion,
            confidence="medium",
            score=0.0,
            method="ai_unconstrained",
            n_markers_found=0,
            ai_agreed=True,
            ai_suggested=ai_suggestion,
            explanation=f"Unconstrained AI mode — accepted AI suggestion '{ai_suggestion}' (no KB match).",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category="weak_signal" if marker_scores else "no_kb_match",
                top_competitors=[],
                detail=f"AI assigned '{ai_suggestion}' in unconstrained mode.",
            ),
        )

    # ── No scores → early exit ─────────────────────────────────────────
    if not marker_scores:
        return FusionDecision(
            cell_type="Unknown",
            confidence="unknown",
            score=0.0,
            method="unknown",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested=ai_suggestion or "",
            explanation="No marker scores available for this cluster.",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category="true_unknown",
                top_competitors=[],
                detail="No marker scores calculated — empty or missing data.",
            ),
        )

    # ── Find the best-scoring cell type ─────────────────────────────────
    best_type, best_score, n_markers = _find_best_type(marker_scores)

    # ── Transition State Detection ──
    # Check if this cluster is in transition between two closely-scoring
    # types of the same lineage. Runs BEFORE the normal tier loop.
    # Broad_* keys are already stripped from marker_scores by the
    # caller (run_unified_annotation), so _find_best_type() and
    # _is_transition_state() only see fine-grained types.
    if kb is not None and allows_transitions:
        transition = _is_transition_state(
            marker_scores,
            kb,
            incompatible_transitions=incompatible_transitions,
        )
        if transition is not None:
            t1, t2 = transition
            # Multi-peak downgrade (D3): applies only to confirmed transition
            # candidates — whole-dataset Fisher saturation must not degrade
            # correctly annotated clusters (rerun: 14/14 clusters ≥4 types ≥0.9,
            # 9 degraded incl. former MG/RGC). Non-candidates keep their
            # marker-scoring labels; → METC for root-cause.
            ranked = sorted(
                marker_scores.items(),
                key=lambda kv: _resolve_score(marker_scores, kv[0])[0],
                reverse=True,
            )
            floor = multi_peak_score_floor
            tied = [
                (k, _resolve_score(marker_scores, k)[0])
                for k, _ in ranked
                if _resolve_score(marker_scores, k)[0] >= floor
            ]
            if len(tied) >= multi_peak_min_types:
                top3 = ", ".join(f"{k}={s:.3f}" for k, s in tied[:3])
                candidate: "FusionDecision" = FusionDecision(
                    cell_type="Unknown",
                    confidence="unknown",
                    score=best_score,
                    method="ambiguous",
                    n_markers_found=n_markers,
                    ai_agreed=False,
                    ai_suggested=ai_suggestion or "",
                    explanation=_explain(
                        "Unknown",
                        "ambiguous",
                        best_score,
                        n_markers,
                        best_type,
                        ai_suggestion,
                        False,
                    ),
                    alternative_rules=[],
                    diagnostic=DiagnosticInfo(
                        category="ambiguous",
                        top_competitors=[{"cell_type": k, "score": s} for k, s in tied],
                        detail=f"Multi-peak tie: {len(tied)} types >= {floor} (top: {top3})",
                    ),
                )
            else:
                delta = abs(
                    _resolve_score(marker_scores, t1)[0] - _resolve_score(marker_scores, t2)[0]
                )
                parent = kb.get(t1, {}).get("parent", "")
                explanation = _explain(
                    cell_type=f"transitional: {t1}/{t2}",
                    method="transition_state",
                    score=best_score,
                    n_markers=n_markers,
                    best_type=t1,
                    ai_suggestion=ai_suggestion,
                    ai_agreed=False,
                    confidence="transition",
                    alternative_rules=[
                        f"Top-2 scores within {delta:.3f}, shared lineage {parent}",
                        f"  {t1}: {_resolve_score(marker_scores, t1)[0]:.3f}",
                        f"  {t2}: {_resolve_score(marker_scores, t2)[0]:.3f}",
                    ],
                )
                candidate = FusionDecision(
                    cell_type=f"transitional: {t1}/{t2}",
                    confidence="transition",
                    score=best_score,
                    method="transition_state",
                    n_markers_found=n_markers,
                    ai_agreed=False,
                    ai_suggested=ai_suggestion if ai_suggestion else "",
                    explanation=explanation,
                    alternative_rules=[],
                    diagnostic=None,
                    cell_category="",
                )

            # ── Layer 3/4 candidate exit semantics (plan todo 4) ─────────
            # The candidate (ambiguous or transition_state) is the terminal
            # output for this cluster.  KADP naming / METC arbitration run
            # only when explicitly enabled; a KADP miss or an un-arbitrated
            # candidate is returned unchanged — byte-identical to the baseline
            # early returns above.  The tier loop below is reachable only for
            # non-ambiguous/transition_state candidates.
            if (kadp_cfg and kadp_cfg.enabled) or (metc_cfg and metc_cfg.enabled):
                if kadp_cfg and kadp_cfg.enabled and candidate.method == "ambiguous":
                    named = _kadp_name_candidate(
                        marker_scores=marker_scores,
                        kb=kb,
                        candidate=candidate,
                        cfg=kadp_cfg,
                    )
                    if named is not None:
                        return named
                # METC arbitration lands in todo 9; until then (and for any
                # non-arbitrated candidate incl. n_spoke < 3) the candidate is
                # returned unchanged — never falls into the tier loop.
                return candidate
            return candidate

    # ── Apply tiers 1–4 ────────────────────────────────────────────────
    for tier_name, tier_fn in DECISION_TIERS:
        if tier_name == "expert_rule":
            continue  # already handled above

        if not tier_fn(best_score, expert_rule_result, ai_suggestion):
            continue

        if tier_name == "unknown":
            diagnostic = _classify_unknown(
                marker_scores,
                low_quality_reason=low_quality_reason,
            )
            explanation = _explain(
                "Unknown",
                "unknown",
                best_score,
                n_markers,
                best_type,
                ai_suggestion,
                False,
            )
            if weak_rule_override_label is not None:
                explanation = (
                    f"Rule '{weak_rule_override_label}' overridden by strong "
                    f"marker scoring (weak rule overridden) | {explanation}"
                )
            return FusionDecision(
                cell_type="Unknown",
                confidence="unknown",
                score=best_score,
                method="unknown",
                n_markers_found=n_markers,
                ai_agreed=False,
                ai_suggested=ai_suggestion or "",
                explanation=explanation,
                alternative_rules=alternative_rules or [],
                diagnostic=diagnostic,
            )

        # Tiers 1-3: marker-scoring-based decisions.
        ai_agreed = _labels_match(ai_suggestion, best_type) if ai_suggestion else False

        # D2 weak-evidence cap: confidence is derived from *evidence strength*,
        # not from the tier alone.  A winner whose evidence_type is in
        # WEAK_EVIDENCE (single-marker / window-padding / weak-multi /
        # zero-evidence) or is "ai_only" is capped at confidence="low"
        # regardless of the tier mapping, and the cluster is flagged for the
        # engine-side review queue via review_reason=evidence_type.
        winner = marker_scores.get(best_type)
        winner_evidence = (
            getattr(winner, "evidence_type", "") if not isinstance(winner, (int, float)) else ""
        )
        if winner_evidence in WEAK_EVIDENCE or winner_evidence == "ai_only":
            confidence = "low"
            review_reason = winner_evidence
        else:
            confidence = _CONFIDENCE_MAP[tier_name]
            review_reason = ""

        explanation = _explain(
            best_type,
            tier_name,
            best_score,
            n_markers,
            best_type,
            ai_suggestion,
            ai_agreed,
        )
        if weak_rule_override_label is not None:
            explanation = (
                f"Rule '{weak_rule_override_label}' overridden by strong "
                f"marker scoring (weak rule overridden) | {explanation}"
            )
        return FusionDecision(
            cell_type=best_type,
            confidence=confidence,
            score=best_score,
            method=tier_name,
            n_markers_found=n_markers,
            ai_agreed=ai_agreed,
            ai_suggested=ai_suggestion or "",
            explanation=explanation,
            alternative_rules=[],
            review_reason=review_reason,
        )

    # Fallback (should never reach here — 'unknown' always matches)
    return FusionDecision(
        "Unknown", "unknown", 0.0, "unknown", 0, False, "", "Fallback: no tier matched.", []
    )


def fuse_all_clusters(
    all_scores: dict,
    all_rules: dict,
    kb: Optional[dict] = None,
    all_marker_dfs: Optional[pd.DataFrame] = None,
    ai_results: Optional[dict] = None,
    return_quality: bool = False,
    low_quality_clusters: Optional[dict] = None,
    unconstrained: bool = False,
    allows_transitions: bool = False,
    incompatible_transitions: Optional[list] = None,
    celltypist_results: Optional[dict] = None,
    multi_peak_min_types: int = 3,
    multi_peak_score_floor: float = 0.9,
) -> list | tuple[list, dict]:
    """Process all clusters and return a list of :class:`FusionDecision`.

    Parameters
    ----------
    all_scores : dict
        ``{cluster_id: {type_key: Score}}``.
    all_rules : dict
        ``{cluster_id: expert_rule_result_or_None}``.
    kb : dict or None
        Full KB dict (passed through to :func:`fuse_evidence`).
    all_marker_dfs : pd.DataFrame or None
        Concatenated ``rank_genes_groups`` output with a ``cluster`` column.
    ai_results : dict or None
        ``{cluster_id: AI-suggested cell type}``.
    return_quality : bool
        When ``True``, also return a quality metadata dict
        ``{annotated_by_rule, unknown, ambiguity, ai_agreed}``.
    low_quality_clusters : dict or None
        ``{cluster_id: reason_str}`` from
        :func:`~utils.marker_scoring.detect_low_quality_cluster` (v3.1.0+).
    allows_transitions : bool
        When ``True``, enables transition-state detection between closely
        related cell types.  Passed through to each
        :func:`fuse_evidence` call.
    incompatible_transitions : list or None
        List of [type_a, type_b] pairs explicitly forbidden as
        developmental transitions.  Passed through to each
        :func:`fuse_evidence` call.
    multi_peak_min_types : int
        Passed through to :func:`fuse_evidence` (default 3).
    multi_peak_score_floor : float
        Passed through to :func:`fuse_evidence` (default 0.9).

    Returns
    -------
    list[FusionDecision]  or  tuple[list[FusionDecision], dict]
        One decision per cluster, sorted by cluster id.
        When *return_quality* is ``True``, returns ``(decisions, quality)``.
    """
    if ai_results is None:
        ai_results = {}
    if low_quality_clusters is None:
        low_quality_clusters = {}

    decisions: list = []
    clusters = sorted(
        all_scores.keys(),
        key=lambda x: int(x) if str(x).isdigit() else str(x),
    )

    for cl in clusters:
        cl_markers = None
        if all_marker_dfs is not None and "cluster" in all_marker_dfs.columns:
            cl_mask = all_marker_dfs["cluster"] == cl
            cl_markers = all_marker_dfs[cl_mask].copy()

        rule_value = all_rules.get(cl)
        if isinstance(rule_value, tuple):
            rule_result, alt_rules = rule_value
        else:
            rule_result, alt_rules = rule_value, []

        decision = fuse_evidence(
            marker_scores=all_scores.get(cl, {}),
            expert_rule_result=rule_result,
            kb=kb,
            cluster_markers=cl_markers,
            ai_suggestion=ai_results.get(cl),
            alternative_rules=alt_rules,
            low_quality_reason=low_quality_clusters.get(str(cl), ""),
            unconstrained=unconstrained,
            allows_transitions=allows_transitions,
            incompatible_transitions=incompatible_transitions,
            multi_peak_min_types=multi_peak_min_types,
            multi_peak_score_floor=multi_peak_score_floor,
        )
        decisions.append(decision)

    if return_quality:
        quality = {
            "annotated_by_rule": sum(1 for d in decisions if d.method == "expert_rule"),
            "annotated_by_scoring": sum(
                1 for d in decisions if d.method.startswith("marker_scoring")
            ),
            "unknown": sum(1 for d in decisions if d.confidence == "unknown"),
            "ambiguity": sum(1 for d in decisions if len(d.alternative_rules) >= 3),
            "ai_agreed": sum(1 for d in decisions if d.ai_agreed),
            "total": len(decisions),
            "diagnostic_summary": _build_diagnostic_summary(decisions),
            "celltypist": bool(celltypist_results) if celltypist_results else False,
        }
        return decisions, quality
    return decisions
