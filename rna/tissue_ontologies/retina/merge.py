"""
tissue_ontologies/retina/merge.py — Retina KB merge engine.

Loads 7 source files from ``sources/``, merges their markers with consensus
scoring, detects conflicts, and emits a unified KB dict consumable by
``marker_scoring.py``.

Usage::

    from tissue_ontologies.retina.merge import retina_expert_kb
"""

import importlib.util
import logging
import os
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════
#  Type-key synonym table for cross-source merging
# ═══════════════════════════════════════════════════════════════════════

TYPE_ALIASES: Dict[str, str] = {
    "Retinal_Ganglion_Cell": "RGC",
}

# ═══════════════════════════════════════════════════════════════════════
#  Source loading
# ═══════════════════════════════════════════════════════════════════════


def load_all_sources(sources_dir: Optional[str] = None) -> List[Dict[str, Any]]:
    """Auto-discover and import all ``.py`` source files from ``sources/``.

    Excludes files whose name starts with ``_`` (notably ``_TEMPLATE.py``).

    Parameters
    ----------
    sources_dir : str or None
        Path to the sources directory.  When ``None`` (default) it is
        inferred relative to this file's location.

    Returns
    -------
    list[dict]
        Each entry has keys ``meta``, ``markers``, ``novel_types``,
        ``expert_rules``, ``conflicts``, ``source_id``.
    """
    if sources_dir is None:
        sources_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "sources"
        )

    sources: List[Dict[str, Any]] = []
    entries = sorted(os.listdir(sources_dir))

    for entry in entries:
        if not entry.endswith(".py"):
            continue
        if entry.startswith("_"):
            continue  # skip _TEMPLATE.py etc.

        module_name = entry[:-3]
        filepath = os.path.join(sources_dir, entry)

        spec = importlib.util.spec_from_file_location(module_name, filepath)
        if spec is None or spec.loader is None:
            logger.warning("Cannot load source module: %s", entry)
            continue

        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except Exception:
            logger.exception("Error executing source module: %s", entry)
            continue

        meta = getattr(mod, "source_meta", {})
        sources.append(
            {
                "meta": meta,
                "markers": getattr(mod, "markers", {}),
                "novel_types": getattr(mod, "novel_types", []),
                "expert_rules": getattr(mod, "expert_rules", []),
                "conflicts": getattr(mod, "conflicts", []),
                "source_id": meta.get("id", module_name),
            }
        )

    logger.info("Loaded %d source(s) from %s", len(sources), sources_dir)
    return sources


# ═══════════════════════════════════════════════════════════════════════
#  Consensus level
# ═══════════════════════════════════════════════════════════════════════


def compute_consensus_level(source_count: int, n_sources: int = 7) -> str:
    """Map a marker's source-support count to a qualitative label.

    Thresholds (absolute counts across all available sources):

        * ``gold``   – 5+ sources
        * ``high``   – 3–4 sources
        * ``medium`` – 2 sources
        * ``low``    – 1 source

    Parameters
    ----------
    source_count : int
        Number of distinct sources that list this marker for the type.
    n_sources : int
        Total number of source files available (informational, 7 by default).

    Returns
    -------
    str
        One of ``"gold"``, ``"high"``, ``"medium"``, ``"low"``.
    """
    if source_count >= 5:
        return "gold"
    if source_count >= 3:
        return "high"
    if source_count >= 2:
        return "medium"
    return "low"


# ═══════════════════════════════════════════════════════════════════════
#  Type-key normalisation
# ═══════════════════════════════════════════════════════════════════════


def _normalize_type_key(key: str) -> str:
    """Map a source-internal type key to the canonical KB name.

    Currently handles:

        * ``Retinal_Ganglion_Cell`` → ``RGC``
    """
    return TYPE_ALIASES.get(key, key)


# ═══════════════════════════════════════════════════════════════════════
#  Marker merging
# ═══════════════════════════════════════════════════════════════════════


def merge_markers(sources: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate markers across sources, grouped by canonical cell type.

    Uses :data:`TYPE_ALIASES` to normalise type keys; logs a warning when
    a synonym match is made.

    Parameters
    ----------
    sources : list[dict]
        Source dicts as returned by :func:`load_all_sources`.

    Returns
    -------
    dict
        ``{canonical_type: {...}}`` with internal tracking of which sources
        contributed each marker (see source for full key list).
    """
    merged: Dict[str, Any] = {}

    for src in sources:
        src_id = src["source_id"]
        src_class = src["meta"].get("class", "")
        src_order = src["meta"].get("order", "")

        # ── Main markers dict ─────────────────────────────────────
        for raw_key, marker_data in src.get("markers", {}).items():
            canonical = _normalize_type_key(raw_key)

            if canonical not in merged:
                merged[canonical] = {
                    "confirm": {},
                    "add": {},
                    "refine": {},
                    "negative_markers": set(),
                    "species": set(),
                    "synonyms": set(),
                    "parent": "",
                    "source_ids": set(),
                    "classes": set(),
                    "orders": set(),
                }

            entry = merged[canonical]
            entry["source_ids"].add(src_id)
            if src_class:
                entry["classes"].add(src_class)
            if src_order:
                entry["orders"].add(src_order)

            if raw_key != canonical:
                entry["synonyms"].add(raw_key)

            # Confirm
            for gene in marker_data.get("confirm", {}):
                _register_marker(entry["confirm"], gene, src_id)

            # Add
            for gene in marker_data.get("add", {}):
                _register_marker(entry["add"], gene, src_id)

            # Refine
            for gene, refine_data in marker_data.get("refine", {}).items():
                entry["refine"].setdefault(gene, []).append(refine_data)

            # Species from source meta
            entry["species"].update(src["meta"].get("species", []))

        # ── Novel types ───────────────────────────────────────────
        for nt in src.get("novel_types", []):
            nt_name = nt.get("name", "")
            if not nt_name:
                continue

            canonical = _normalize_type_key(nt_name)

            if canonical not in merged:
                merged[canonical] = {
                    "confirm": {},
                    "add": {},
                    "refine": {},
                    "negative_markers": set(),
                    "species": set(),
                    "synonyms": set(),
                    "parent": "",
                    "source_ids": set(),
                    "classes": set(),
                    "orders": set(),
                }

            entry = merged[canonical]
            entry["source_ids"].add(src_id)
            if src_class:
                entry["classes"].add(src_class)
            if src_order:
                entry["orders"].add(src_order)
            if nt.get("parent"):
                entry["parent"] = nt["parent"]

            # Novel-type markers go into "add" (they are novel per definition)
            for gene in nt.get("markers", []):
                _register_marker(entry["add"], gene, src_id)

            for sp in nt.get("species", []):
                entry["species"].add(sp)

    return merged


def _register_marker(
    dest: Dict[str, Dict[str, Any]], gene: str, src_id: str
) -> None:
    """Register *src_id* as a source for *gene* in *dest*.

    *dest* is one of the ``"confirm"`` or ``"add"`` sub-dicts inside a
    merged-type entry.
    """
    if gene not in dest:
        dest[gene] = {"source_ids": [], "source_count": 0}
    if src_id not in dest[gene]["source_ids"]:
        dest[gene]["source_ids"].append(src_id)
        dest[gene]["source_count"] += 1


# ═══════════════════════════════════════════════════════════════════════
#  Conflict detection
# ═══════════════════════════════════════════════════════════════════════


def detect_conflicts(sources: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """Find annotation conflicts between sources.

    Two kinds of conflict are reported:

    1. **Explicit** – conflicts that are declared in the source data
       (the ``conflicts`` list in each source file).
    2. **Cross-type** – a gene that serves as a positive marker for two
       different cell types across different sources.  These *may* be
       legitimate (e.g., CALB1 marks both amacrine and horizontal cells),
       but flagging them helps downstream review.

    Returns
    -------
    dict
        ``{"explicit": [...], "cross_type": [...]}``
    """
    explicit: List[Dict[str, Any]] = []
    for src in sources:
        for c in src.get("conflicts", []):
            explicit.append(
                {
                    "type_a": c.get("type_a", {}),
                    "type_b": c.get("type_b", {}),
                    "notes": c.get("notes", ""),
                    "source": c.get("source", {}),
                }
            )

    # Cross-type detection: build gene → {type → set(source_ids)}
    gene_type_map: Dict[str, Dict[str, Set[str]]] = {}
    for src in sources:
        src_id = src["source_id"]
        for raw_key, marker_data in src.get("markers", {}).items():
            canonical = _normalize_type_key(raw_key)
            for tier in ("confirm", "add"):
                for gene in marker_data.get(tier, {}):
                    gene_type_map.setdefault(gene, {}).setdefault(canonical, set()).add(
                        src_id
                    )

    cross_type: List[Dict[str, Any]] = []
    for gene, type_map in gene_type_map.items():
        types = list(type_map.keys())
        if len(types) < 2:
            continue
        # Report all type pairs
        for i in range(len(types)):
            for j in range(i + 1, len(types)):
                t1, t2 = types[i], types[j]
                cross_type.append(
                    {
                        "type_a": {"cell_type": t1, "marker": gene},
                        "type_b": {"cell_type": t2, "marker": gene},
                        "notes": f"'{gene}' is a positive marker for both "
                        f"{t1} and {t2} across different sources",
                        "source": {"a": list(type_map[t1]), "b": list(type_map[t2])},
                    }
                )

    return {"explicit": explicit, "cross_type": cross_type}


# ═══════════════════════════════════════════════════════════════════════
#  Conflict resolution
# ═══════════════════════════════════════════════════════════════════════


def resolve_conflicts(
    conflicts: Dict[str, List[Dict[str, Any]]], sources: List[Dict[str, Any]]
) -> Dict[str, List[Dict[str, Any]]]:
    """Resolve detectable conflicts; flag the rest for manual review.

    **Auto-resolve** (cross-type only):
    If a gene is a positive marker for type A in >= 3 sources and for type B
    in strictly fewer sources, the conflict is auto-resolved in favour of
    type A (majority wins).

    **Flagged** (always):
        - All explicitly reported conflicts from source ``conflicts`` lists.
        - Cross-type conflicts where neither side has a clear majority.

    Returns
    -------
    dict
        ``{"resolved": [...], "flagged": [...]}``
    """
    resolved: List[Dict[str, Any]] = []
    flagged: List[Dict[str, Any]] = []

    total = len(sources)

    # Cross-type auto-resolution
    for c in conflicts.get("cross_type", []):
        src_a = set(c["source"].get("a", []))
        src_b = set(c["source"].get("b", []))
        n_a, n_b = len(src_a), len(src_b)

        if n_a >= 3 and n_a > n_b:
            resolved.append(
                {
                    **c,
                    "resolution": (
                        f"Prefer '{c['type_a']['cell_type']}' "
                        f"({n_a}/{total} sources vs {n_b})"
                    ),
                }
            )
        elif n_b >= 3 and n_b > n_a:
            resolved.append(
                {
                    **c,
                    "resolution": (
                        f"Prefer '{c['type_b']['cell_type']}' "
                        f"({n_b}/{total} sources vs {n_a})"
                    ),
                }
            )
        else:
            flagged.append(
                {
                    **c,
                    "reason": f"Cannot auto-resolve ({n_a} vs {n_b} sources)",
                }
            )

    # Explicit conflicts are always flagged
    for c in conflicts.get("explicit", []):
        flagged.append({**c, "reason": "Explicitly reported conflict"})

    return {"resolved": resolved, "flagged": flagged}


# ═══════════════════════════════════════════════════════════════════════
#  Expert rules merging
# ═══════════════════════════════════════════════════════════════════════


def _rule_dedup_key(rule: Dict[str, Any]) -> str:
    """Deterministic string key for a rule (condition + action).

    Normalises the action through ``TYPE_ALIASES`` so that rules with
    ``Retinal_Ganglion_Cell`` and ``RGC`` actions are treated as duplicates.
    Includes ``markers_absent`` in the key so that rules that differ only
    by exclusion markers are NOT treated as duplicates.
    """
    condition = rule.get("condition", {})
    action = _normalize_type_key(rule.get("action", ""))
    markers = condition.get("markers_present", {})
    absent = condition.get("markers_absent", [])
    sorted_genes = sorted(markers.keys())
    sorted_absent = sorted(absent) if absent else []
    marker_str = ",".join(f"{g}:{markers[g]}" for g in sorted_genes)
    absent_str = "!" + ",".join(sorted_absent) if sorted_absent else ""
    return f"{action}|{marker_str}{absent_str}"


def merge_rules(sources: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Merge expert rules from all sources, deduplicated by condition+action.

    When the same logical rule appears in multiple sources, the entry with
    the **lowest** priority number (highest priority) is kept.

    Returns
    -------
    list[dict]
        Rules sorted by ``priority`` ascending.
    """
    seen: Dict[str, Dict[str, Any]] = {}

    for src in sources:
        for rule in src.get("expert_rules", []):
            # Normalize action through type aliases
            rule_normalized = dict(rule)
            rule_normalized["action"] = _normalize_type_key(rule.get("action", ""))
            key = _rule_dedup_key(rule_normalized)
            existing_priority = seen.get(key, {}).get("priority", 999)
            if key not in seen or rule_normalized.get("priority", 999) < existing_priority:
                seen[key] = rule_normalized

    return sorted(seen.values(), key=lambda r: r.get("priority", 999))


# ═══════════════════════════════════════════════════════════════════════
#  KB assembly
# ═══════════════════════════════════════════════════════════════════════


def build_final_kb(
    merged_types: Dict[str, Any],
    merged_rules: List[Dict[str, Any]],
    sources: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Build the unified KB dict consumable by ``marker_scoring.py``.

    The output format matches what :func:`utils.marker_scoring._build_kb_lookup`
    expects::

        {
            "CellTypeA": {
                "markers": {
                    "confirm": {"GENE1": ["src1_id", "src2_id"], ...},
                    "add":    {"GENEn": ["srcN_id"], ...},
                    "refine": {"GENE": {"note": "...", "threshold": "...",
                                        "pmid": "..."}, ...},
                },
                "negative_markers": [...],
                "species": [...],
                "synonyms": [...],
                "parent": "...",
                "class": "...",       # str (most common class across sources)
                "order": "...",       # str (most common order across sources)
                "classes": [...],     # list[str] all contributing classes
                "orders": [...],      # list[str] all contributing orders
            },
            "expert_rules": [...],
            "_meta": {
                "total_sources": N,
                "classes": {"Mammalia", ...},
                "orders": {"Primates", ...},
            },
        }

    Parameters
    ----------
    merged_types : dict
        Output of :func:`merge_markers`.
    merged_rules : list[dict]
        Output of :func:`merge_rules`.
    sources : list[dict]
        Source list (used to determine ``total_sources`` for consensus).

    Returns
    -------
    dict
        Complete KB ready for consumption by ``marker_scoring.py``.
    """
    total_sources = len(sources)
    kb: Dict[str, Any] = {}

    # Collect all classes and orders for _meta
    all_classes: Set[str] = set()
    all_orders: Set[str] = set()

    for type_key, type_data in merged_types.items():
        confirm_out: Dict[str, List[str]] = {}
        add_out: Dict[str, List[str]] = {}
        refine_out: Dict[str, Dict[str, str]] = {}

        # Confirm
        for gene, info in type_data.get("confirm", {}).items():
            confirm_out[gene] = info["source_ids"]

        # Add
        for gene, info in type_data.get("add", {}).items():
            add_out[gene] = info["source_ids"]

        # Refine — merge multiple sources by concatenating notes
        for gene, refine_list in type_data.get("refine", {}).items():
            if not refine_list:
                continue
            first = refine_list[0]
            merged_refine: Dict[str, str] = {
                "note": first.get("note", ""),
                "threshold": first.get("threshold", ""),
                "pmid": first.get("pmid", ""),
            }
            n_extra = len(refine_list) - 1
            if n_extra > 0:
                merged_refine["note"] = (
                    first.get("note", "")
                    + f" (refined by {n_extra + 1} source(s))"
                )
            refine_out[gene] = merged_refine

        # Consensus levels
        consensus_levels: Dict[str, str] = {}
        gene_sources: Dict[str, Set[str]] = {}
        for gene, info in type_data.get("confirm", {}).items():
            gene_sources.setdefault(gene, set()).update(info["source_ids"])
        for gene, info in type_data.get("add", {}).items():
            gene_sources.setdefault(gene, set()).update(info["source_ids"])
        for gene, src_set in gene_sources.items():
            consensus_levels[gene] = compute_consensus_level(
                len(src_set), total_sources
            )

        # Resolve class/order — use most common; fall back to sorted list
        classes_list = sorted(type_data.get("classes", set()))
        orders_list = sorted(type_data.get("orders", set()))
        resolved_class = classes_list[0] if len(classes_list) == 1 else (
            ", ".join(classes_list) if classes_list else ""
        )
        resolved_order = orders_list[0] if len(orders_list) == 1 else (
            ", ".join(orders_list) if orders_list else ""
        )

        all_classes.update(classes_list)
        all_orders.update(orders_list)

        kb[type_key] = {
            "markers": {
                "confirm": confirm_out,
                "add": add_out,
                "refine": refine_out,
            },
            "negative_markers": sorted(type_data.get("negative_markers", set())),
            "species": sorted(type_data.get("species", set())),
            "synonyms": sorted(type_data.get("synonyms", set())),
            "parent": type_data.get("parent", ""),
            "consensus_levels": consensus_levels,
            "class": resolved_class,
            "order": resolved_order,
            "classes": classes_list,
            "orders": orders_list,
        }

    kb["expert_rules"] = merged_rules
    kb["_meta"] = {
        "total_sources": total_sources,
        "classes": sorted(all_classes),
        "orders": sorted(all_orders),
    }
    return kb


# ═══════════════════════════════════════════════════════════════════════
#  Module-level convenience: build the KB once at import time
# ═══════════════════════════════════════════════════════════════════════

_sources: List[Dict[str, Any]] = load_all_sources()
_merged_types: Dict[str, Any] = merge_markers(_sources)
_merged_rules: List[Dict[str, Any]] = merge_rules(_sources)
_conflicts: Dict[str, List[Dict[str, Any]]] = detect_conflicts(_sources)
_resolved: Dict[str, List[Dict[str, Any]]] = resolve_conflicts(_conflicts, _sources)
retina_expert_kb: Dict[str, Any] = build_final_kb(
    _merged_types, _merged_rules, _sources
)

# Log a summary of the built KB
logger.info(
    "Built retina_expert_kb: %d types, %d rules, %d conflicts flagged",
    sum(1 for k in retina_expert_kb if k != "expert_rules"),
    len(retina_expert_kb.get("expert_rules", [])),
    len(_resolved.get("flagged", [])),
)
