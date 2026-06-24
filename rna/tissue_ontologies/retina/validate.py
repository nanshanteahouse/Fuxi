"""
tissue_ontologies/retina/validate.py — KB validation routines.

Usage::

    from tissue_ontologies.retina import retina_expert_kb
    from tissue_ontologies.retina.validate import validate_kb

    is_valid, errors = validate_kb(retina_expert_kb)
    if not is_valid:
        for e in errors:
            print(f"  ERROR: {e}")
"""

import logging
from typing import Any, Dict, List, Set, Tuple

logger = logging.getLogger(__name__)


def validate_kb(kb: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """Validate a merged KB for structural and semantic correctness.

    Checks performed
    -----------------
    1. **Marker minimum** – every cell type must have >= 3 unique positive
       markers.
    2. **Source independence** – markers must come from >= 2 independent
       source files.
    3. **Negative-positive separation** – ``negative_markers`` must not
       contain any gene that also appears in ``confirm`` or ``add``.
    4. **Expert-rule sanity** – no two rules with the same condition map to
       different actions; every rule's ``action`` type must exist in the KB.
    5. **Species well-formedness** – all ``species`` entries are non-empty
       strings.
    6. **Consensus-level integrity** – every marker gene listed in
       ``confirm`` or ``add`` has a ``consensus_levels`` entry.

    Parameters
    ----------
    kb : dict
        The KB produced by :func:`tissue_ontologies.retina.merge.build_final_kb`.

    Returns
    -------
    (bool, list[str])
        ``(is_valid, [error_message, ...])``.
    """
    errors: List[str] = []

    # Collect all known type keys (excluding "expert_rules" and "_meta")
    type_keys: Set[str] = {
        k for k in kb
        if k not in ("expert_rules", "_meta")
        and not k.startswith("_")
        and isinstance(kb[k], dict)
    }

    # Known species across the entire KB
    known_species: Set[str] = set()

    # ── 1-3. Per-type checks ────────────────────────────────────────
    for type_key in sorted(type_keys):
        type_data = kb[type_key]
        markers = type_data.get("markers", {})
        confirm: Dict[str, Any] = markers.get("confirm", {})
        add: Dict[str, Any] = markers.get("add", {})
        refine: Dict[str, Any] = markers.get("refine", {})

        positive_genes: Set[str] = set(confirm.keys()) | set(add.keys())

        # 1. Minimum 3 markers
        if len(positive_genes) < 3:
            errors.append(
                f"'{type_key}': only {len(positive_genes)} positive marker(s) "
                f"(need >= 3)"
            )

        # 2. Sources >= 2
        all_sources: Set[str] = set()
        for gene, src_list in confirm.items():
            if isinstance(src_list, (list, tuple)):
                all_sources.update(src_list)
        for gene, src_list in add.items():
            if isinstance(src_list, (list, tuple)):
                all_sources.update(src_list)
        if len(all_sources) < 2:
            errors.append(
                f"'{type_key}': markers from only {len(all_sources)} "
                f"independent source(s) (need >= 2)"
            )

        # 3. Negative vs positive overlap
        neg: List[str] = type_data.get("negative_markers", [])
        neg_set: Set[str] = set(neg)
        overlap: Set[str] = neg_set & positive_genes
        if overlap:
            errors.append(
                f"'{type_key}': negative_markers overlap with positive "
                f"markers: {sorted(overlap)}"
            )

        # 4. Consensus-levels completeness
        consensus_levels: Dict[str, str] = type_data.get("consensus_levels", {})
        for gene in positive_genes:
            if gene not in consensus_levels:
                errors.append(
                    f"'{type_key}': marker '{gene}' missing from consensus_levels"
                )

        # Collect species
        species_list = type_data.get("species", [])
        known_species.update(species_list)

    # ── 4. Expert-rule checks ───────────────────────────────────────
    rules: List[Dict[str, Any]] = kb.get("expert_rules", [])
    seen_conditions: Dict[str, str] = {}  # condition_key → action

    for idx, rule in enumerate(rules):
        action = rule.get("action", "")
        if action not in type_keys:
            errors.append(
                f"expert_rules[{idx}]: action '{action}' does not match "
                f"any KB cell type"
            )

        condition = rule.get("condition", {})
        markers_present = condition.get("markers_present", {})
        if not markers_present:
            errors.append(f"expert_rules[{idx}]: empty condition.markers_present")

        # Build deterministic key for condition (include absent markers)
        markers_absent = condition.get("markers_absent", [])
        absent_str = "!" + ",".join(sorted(markers_absent)) if markers_absent else ""
        cond_key = ";".join(
            f"{g}:{markers_present[g]}" for g in sorted(markers_present.keys())
        ) + absent_str

        if cond_key in seen_conditions and seen_conditions[cond_key] != action:
            errors.append(
                f"expert_rules[{idx}]: contradictory rules – condition "
                f"'{cond_key}' maps to both '{seen_conditions[cond_key]}' "
                f"and '{action}'"
            )
        seen_conditions.setdefault(cond_key, action)

    # ── 5. Species well-formedness ──────────────────────────────────
    for sp in known_species:
        if not sp or not isinstance(sp, str):
            errors.append(f"Invalid species entry: {sp!r}")

    if not known_species:
        errors.append("No species found in any cell type entry")

    return len(errors) == 0, errors
