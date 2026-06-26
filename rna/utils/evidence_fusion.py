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

from typing import NamedTuple, Optional

import pandas as pd


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
        ``'high'`` | ``'medium'`` | ``'low'`` | ``'unknown'`` | ``'rule'``.
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


# Decision priority tiers — evaluated in order.
# Each tier is a (name, callable) where
# callable(score, expert_rule_result, ai_suggestion) → bool.
DECISION_TIERS = [
    ('expert_rule',           lambda s, e, a: e is not None),          # Tier 0
    ('marker_scoring_high',   lambda s, e, a: s >= 0.7),               # Tier 1
    ('marker_scoring_medium', lambda s, e, a: 0.5 <= s < 0.7),         # Tier 2
    ('marker_scoring_low',    lambda s, e, a: 0.25 <= s < 0.5),        # Tier 3
    ('unknown',               lambda s, e, a: True),                    # Tier 4
]

_CONFIDENCE_MAP = {
    'expert_rule': 'rule',
    'marker_scoring_high': 'high',
    'marker_scoring_medium': 'medium',
    'marker_scoring_low': 'low',
    'unknown': 'unknown',
}


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
    nfkd = unicodedata.normalize('NFKD', label)
    # Drop combining chars and other non-ASCII
    ascii_label = nfkd.encode('ascii', 'ignore').decode('ascii')
    # Replace all non-alphanumeric runs with a single underscore.
    normalised = re.sub(r'[^a-zA-Z0-9]+', '_', ascii_label)
    return normalised.strip('_').lower()


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


def _explain(
    cell_type: str,
    method: str,
    score: float,
    n_markers: int,
    best_type: Optional[str],
    ai_suggestion: Optional[str],
    ai_agreed: bool,
    alternative_rules: Optional[list] = None,
) -> str:
    """Build a human-readable explanation."""
    if method == 'expert_rule':
        parts = [f"Expert rule matched: {cell_type}"]
        if alternative_rules and len(alternative_rules) > 1:
            alt_names = [r.get("action") for r in alternative_rules[1:]]
            parts.append(f"(also matched rules: {', '.join(alt_names)})")
        if score > 0:
            parts.append(f"(marker score: {score:.3f})")
    elif method == 'unknown':
        parts = ["No cell type could be confidently assigned"]
        if best_type and score > 0:
            parts.append(f"(best match: {best_type}, score: {score:.3f})")
    else:
        parts = [
            f"Marker scoring selected {cell_type} "
            f"with score {score:.3f}",
        ]
        if n_markers > 0:
            parts.append(f"({n_markers} KB markers found in cluster top-20)")

    if ai_suggestion:
        if ai_agreed:
            parts.append("\u2014 AI agreed with this assignment")
        else:
            parts.append(
                f"\u2014 AI suggested '{ai_suggestion}' "
                f"(different from marker-based result)"
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
    scored = [(k, _resolve_score(marker_scores, k))
              for k in marker_scores]
    scored.sort(key=lambda x: -x[1][0])

    top3 = [(t, round(s, 4)) for t, (s, _) in scored[:3] if s > 0]

    if low_quality_reason:
        return DiagnosticInfo(
            category='low_quality_data',
            top_competitors=top3,
            detail=(
                f"Cluster flagged as low-quality: {low_quality_reason}"
            ),
        )

    if scored and scored[0][1][0] >= 0.25:
        ambiguous_candidates = [(t, round(s, 4))
                                for t, (s, _) in scored if s >= 0.25]
        if len(ambiguous_candidates) >= 2:
            names = ", ".join(t for t, _ in ambiguous_candidates[:5])
            return DiagnosticInfo(
                category='ambiguous',
                top_competitors=top3,
                detail=(
                    f"Multiple cell types with score >= 0.25: {names}"
                ),
            )

    if scored and 0 < scored[0][1][0] < 0.25:
        return DiagnosticInfo(
            category='weak_signal',
            top_competitors=top3,
            detail=(
                f"Best score {scored[0][1][0]:.4f} below 0.25 threshold "
                f"(best type: {scored[0][0]})"
            ),
        )

    if not any(s > 0 for _, (s, _) in scored):
        return DiagnosticInfo(
            category='no_kb_match',
            top_competitors=[],
            detail="No KB cell type had any marker overlap with this cluster.",
        )

    return DiagnosticInfo(
        category='true_unknown',
        top_competitors=top3,
        detail="Could not determine cell type by any method.",
    )


# ═══════════════════════════════════════════════════════════════════════
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
) -> 'FusionDecision':
    """Combine marker scores, expert rules, and AI into one decision.

    Parameters
    ----------
    marker_scores : dict
        Output of :func:`utils.marker_scoring.score_cluster_against_kb`.
        Maps ``type_key → Score`` (or bare ``float`` in simplified contexts).
    expert_rule_result : str or None
        Output of :func:`utils.marker_scoring.apply_expert_rules`.
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

    Returns
    -------
    FusionDecision
    """
    # ── Tier 0: expert rule (highest priority) ─────────────────────────
    if expert_rule_result is not None:
        rule_score, rule_n = _resolve_score(marker_scores, expert_rule_result)
        ai_agreed = _labels_match(ai_suggestion, expert_rule_result) if ai_suggestion else False

        # Quality gate (v3.1.0+): if Fisher scoring completely disagrees
        # with the expert rule (zero KB marker overlap), downgrade confidence
        # from 'rule' to 'low'.  This prevents noise-triggered rules (e.g.
        # a gene buried at rank 4000 in relaxed mode) from outranking well-
        # scored Fisher matches in downstream analysis.
        if rule_score < 0.25 and rule_n == 0:
            conf = 'low'
            warning_note = (
                f"Expert rule matched '{expert_rule_result}' but independent "
                f"marker scoring found zero KB marker overlap (score={rule_score:.3f}, "
                f"n_markers=0). Downgrading confidence from 'rule' to 'low'."
            )
        else:
            conf = 'rule'
            warning_note = ""

        explanation_parts = []
        if warning_note:
            explanation_parts.append(warning_note)
        explanation_parts.append(_explain(
            expert_rule_result, 'expert_rule', rule_score, rule_n,
            expert_rule_result, ai_suggestion, ai_agreed,
            alternative_rules=alternative_rules,
        ))

        return FusionDecision(
            cell_type=expert_rule_result,
            confidence=conf,
            score=rule_score,
            method='expert_rule',
            n_markers_found=rule_n,
            ai_agreed=ai_agreed,
            ai_suggested=ai_suggestion or '',
            explanation=" | ".join(explanation_parts),
            alternative_rules=alternative_rules or [],
        )

    # ── Unconstrained AI mode: accept AI suggestion directly ──────────
    if unconstrained and ai_suggestion and (
        not marker_scores
        or max((_resolve_score(marker_scores, k)[0] for k in marker_scores), default=0) < 0.25
    ):
        return FusionDecision(
            cell_type=ai_suggestion,
            confidence='medium',
            score=0.0,
            method='ai_unconstrained',
            n_markers_found=0,
            ai_agreed=True,
            ai_suggested=ai_suggestion,
            explanation=f"Unconstrained AI mode — accepted AI suggestion '{ai_suggestion}' (no KB match).",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category='weak_signal' if marker_scores else 'no_kb_match',
                top_competitors=[],
                detail=f"AI assigned '{ai_suggestion}' in unconstrained mode.",
            ),
        )

    # ── No scores → early exit ─────────────────────────────────────────
    if not marker_scores:
        return FusionDecision(
            cell_type='Unknown',
            confidence='unknown',
            score=0.0,
            method='unknown',
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested=ai_suggestion or '',
            explanation="No marker scores available for this cluster.",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category='true_unknown',
                top_competitors=[],
                detail="No marker scores calculated — empty or missing data.",
            ),
        )

    # ── Find the best-scoring cell type ─────────────────────────────────
    best_type, best_score, n_markers = _find_best_type(marker_scores)

    # ── Apply tiers 1–4 ────────────────────────────────────────────────
    for tier_name, tier_fn in DECISION_TIERS:
        if tier_name == 'expert_rule':
            continue  # already handled above

        if not tier_fn(best_score, expert_rule_result, ai_suggestion):
            continue

        if tier_name == 'unknown':
            diagnostic = _classify_unknown(
                marker_scores, low_quality_reason=low_quality_reason,
            )
            return FusionDecision(
                cell_type='Unknown',
                confidence='unknown',
                score=best_score,
                method='unknown',
                n_markers_found=n_markers,
                ai_agreed=False,
                ai_suggested=ai_suggestion or '',
                explanation=_explain(
                    'Unknown', 'unknown', best_score, n_markers,
                    best_type, ai_suggestion, False,
                ),
                alternative_rules=alternative_rules or [],
                diagnostic=diagnostic,
            )

        # Tiers 1–3: marker-scoring-based decisions
        if tier_name == 'marker_scoring_medium':
            ai_agreed = _labels_match(ai_suggestion, best_type) if ai_suggestion else True
        elif tier_name == 'marker_scoring_low':
            ai_agreed = _labels_match(ai_suggestion, best_type) if ai_suggestion else False
        else:
            # marker_scoring_high — AI not required for agreement
            ai_agreed = _labels_match(ai_suggestion, best_type) if ai_suggestion else True

        return FusionDecision(
            cell_type=best_type,
            confidence=_CONFIDENCE_MAP[tier_name],
            score=best_score,
            method=tier_name,
            n_markers_found=n_markers,
            ai_agreed=ai_agreed,
            ai_suggested=ai_suggestion or '',
            explanation=_explain(
                best_type, tier_name, best_score, n_markers,
                best_type, ai_suggestion, ai_agreed,
            ),
            alternative_rules=[],
        )

    # Fallback (should never reach here — 'unknown' always matches)
    return FusionDecision('Unknown', 'unknown', 0.0, 'unknown', 0, False, '',
                          'Fallback: no tier matched.', [])


def fuse_all_clusters(
    all_scores: dict,
    all_rules: dict,
    kb: Optional[dict] = None,
    all_marker_dfs: Optional[pd.DataFrame] = None,
    ai_results: Optional[dict] = None,
    return_quality: bool = False,
    low_quality_clusters: Optional[dict] = None,
    unconstrained: bool = False,
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
        if all_marker_dfs is not None and 'cluster' in all_marker_dfs.columns:
            cl_mask = all_marker_dfs['cluster'] == cl
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
        )
        decisions.append(decision)

    if return_quality:
        quality = {
            "annotated_by_rule": sum(
                1 for d in decisions if d.method == "expert_rule"
            ),
            "annotated_by_scoring": sum(
                1 for d in decisions if d.method.startswith("marker_scoring")
            ),
            "unknown": sum(
                1 for d in decisions if d.confidence == "unknown"
            ),
            "ambiguity": sum(
                1 for d in decisions if len(d.alternative_rules) >= 3
            ),
            "ai_agreed": sum(1 for d in decisions if d.ai_agreed),
            "total": len(decisions),
            "diagnostic_summary": _build_diagnostic_summary(decisions),
        }
        return decisions, quality
    return decisions
