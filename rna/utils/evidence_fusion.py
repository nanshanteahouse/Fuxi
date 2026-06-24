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

    Returns
    -------
    FusionDecision
    """
    # ── Tier 0: expert rule (highest priority) ─────────────────────────
    if expert_rule_result is not None:
        rule_score, rule_n = _resolve_score(marker_scores, expert_rule_result)
        ai_agreed = (ai_suggestion == expert_rule_result) if ai_suggestion else False
        return FusionDecision(
            cell_type=expert_rule_result,
            confidence='rule',
            score=rule_score,
            method='expert_rule',
            n_markers_found=rule_n,
            ai_agreed=ai_agreed,
            ai_suggested=ai_suggestion or '',
            explanation=_explain(
                expert_rule_result, 'expert_rule', rule_score, rule_n,
                expert_rule_result, ai_suggestion, ai_agreed,
                alternative_rules=alternative_rules,
            ),
            alternative_rules=alternative_rules or [],
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
            )

        # Tiers 1–3: marker-scoring-based decisions
        if tier_name == 'marker_scoring_medium':
            ai_agreed = (ai_suggestion == best_type) if ai_suggestion else True
        elif tier_name == 'marker_scoring_low':
            ai_agreed = (ai_suggestion == best_type) if ai_suggestion else False
        else:
            # marker_scoring_high — AI not required for agreement
            ai_agreed = (ai_suggestion == best_type) if ai_suggestion else True

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

    Returns
    -------
    list[FusionDecision]  or  tuple[list[FusionDecision], dict]
        One decision per cluster, sorted by cluster id.
        When *return_quality* is ``True``, returns ``(decisions, quality)``.
    """
    if ai_results is None:
        ai_results = {}

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
        }
        return decisions, quality
    return decisions
