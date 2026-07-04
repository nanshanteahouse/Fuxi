"""
marker_expert_rules.py — Expert-rule cell-type assignment (extracted from marker_scoring)

Provides deterministic cell-type matching via expert-authored rules,
independent of the probabilistic marker scoring pipeline.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


# ── Strictness template constants ─────────────────────────────────

_STRICTNESS_TEMPLATES: Dict[str, tuple] = {
    "strict":   (50,    0.01),
    "default":  (50,    0.05),
    "deep":     (200,   0.05),
    "wide":     (1000,  0.05),
    "relaxed":  (5000,  0.05),
    "manual":   (None,  None),
}


# ── Public API ────────────────────────────────────────────────────

def resolve_expert_rule_params(
    strictness: str = "default",
    top_n: int = 0,
    pval_cutoff: float = 0.0,
) -> Tuple[int, float]:
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


def apply_expert_rules(
    kb: Dict[str, Any],
    cluster_markers: pd.DataFrame,
    top_n: int = 50,
    pval_cutoff: float = 0.05,
) -> Tuple[Optional[str], List[Dict[str, Any]]]:
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
