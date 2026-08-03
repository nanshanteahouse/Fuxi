"""
core/kb/merge.py — Tissue-agnostic KB merge engine.

Loads all source YAML files from a ``sources/`` directory, merges their markers
with consensus scoring, detects conflicts, resolves them, and emits a unified KB
dict consumable by ``marker_scoring.py``.

Usage::

    from core.kb.merge import build_tissue_kb
    kb = build_tissue_kb("core/kb/retina/sources",
                         type_aliases={"Retinal_Ganglion_Cell": "RGC"},
                         hierarchy_yaml_path="core/kb/retina/hierarchy.yaml")
"""

import logging
import os
from typing import Any, Dict, List, Optional, Set

import yaml

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
#  Source loading
# ═══════════════════════════════════════════════════════════════════════


def load_all_sources(sources_dir: str) -> List[Dict[str, Any]]:
    """Auto-discover and load all YAML source files from *sources_dir*.

    Excludes ``schema.yaml`` (documentation only) and files starting with ``_``.

    Parameters
    ----------
    sources_dir : str
        Path to the sources directory containing per-publication ``.yaml`` files.

    Returns
    -------
    list[dict]
        Each entry has keys ``meta``, ``markers``, ``novel_types``,
        ``expert_rules``, ``conflicts``, ``source_id``.
    """
    sources: List[Dict[str, Any]] = []
    entries = sorted(os.listdir(sources_dir))

    for entry in entries:
        if not entry.endswith(".yaml"):
            continue
        if entry.startswith("_"):
            continue  # skip _TEMPLATE etc.
        if entry == "schema.yaml":
            continue  # documentation only

        filepath = os.path.join(sources_dir, entry)

        with open(filepath, "r", encoding="utf-8") as fh:
            try:
                data = yaml.safe_load(fh)
            except Exception:
                logger.exception("Error parsing YAML source: %s", entry)
                continue

        if not isinstance(data, dict):
            logger.warning("Skipping non-dict YAML source: %s", entry)
            continue

        meta = data.get("source_meta", {})
        if not meta:
            logger.warning("Skipping source with empty source_meta: %s", entry)
            continue

        module_name = entry[:-5]  # strip '.yaml'
        sources.append(
            {
                "meta": meta,
                "markers": data.get("markers", {}),
                "novel_types": data.get("novel_types", []),
                "expert_rules": data.get("expert_rules", []),
                "conflicts": data.get("conflicts", []),
                "source_id": meta.get("id", module_name),
            }
        )

    logger.info("Loaded %d source(s) from %s", len(sources), sources_dir)
    return sources


# ═══════════════════════════════════════════════════════════════════════
#  Consensus level
# ═══════════════════════════════════════════════════════════════════════


def load_source_independence(yaml_path: str) -> Optional[Dict[str, Any]]:
    """Load the source-independence table (``_source_independence.yaml``).

    Per-source entries carry::

        roots: [GSE...]            # primary root datasets (independent votes)
        modalities: [rna, atac]    # union across roots (method diversity)
        injected_types: {Type: note}   # markers NOT this paper's own data
        validation_datasets: [GSE...]  # datasets where validation ran

    Returns ``None`` when the file does not exist (legacy behaviour).
    """
    if not os.path.isfile(yaml_path):
        return None
    with open(yaml_path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    return {
        "lambda_negative": float(data.get("lambda_negative", 0.5)),
        "sources": data.get("sources", {}),
    }


def compute_effective_source_count(
    gene: str,
    type_key: str,
    src_set: Set[str],
    sources: List[Dict[str, Any]],
    indep: Dict[str, Any],
) -> float:
    """Three-dimensional effective source count for one (gene, type).

    Replaces the naive "number of source files" vote with an evidence-
    independence aware consensus:

    1. **Independent roots** — each source votes via its root datasets;
       distinct root GSEs count once (shared-root re-publication does not
       multiply evidence).
    2. **Method diversity** — a root shared by sources using different
       modalities (e.g. RNA + ATAC) adds +0.5 (chapter-8 rule).
    3. **Negative coverage** — subtract ``lambda_negative`` per *independent*
       validation dataset where this gene was falsified (deduplicated via
       ``validation_datasets``), so a failure recorded by 6 re-copied
       sources still counts as 1 independent negative vote.

    Injected types (canonical copies, gene-score copies, absent from the
    paper's supplement) contribute 0 votes.  Sources without root data are
    treated as a single private root (conservative 1 vote).
    """
    indep_sources = indep.get("sources", {})
    # ── Dimension 1+2: source votes with shared-root discount ─────────
    # Vote unit = source paper (one vote per independent source).
    # A source's root datasets prove independence: unique roots → full
    # vote; roots fully covered by other sources (pure data re-publication)
    # → half vote (chapter-8 L1: other's data, own analysis).
    # Injected types → 0 votes.  No root info → conservative 1 vote.
    votes: Dict[str, float] = {}
    root_of: Dict[str, Set[str]] = {}
    for sid in src_set:
        info = indep_sources.get(sid, {})
        if type_key in info.get("injected_types", {}):
            votes[sid] = 0.0
            continue
        roots = set(info.get("roots") or [])
        if not roots:
            votes[sid] = 1.0
        else:
            votes[sid] = 1.0
            for r in roots:
                root_of.setdefault(r, set()).add(sid)

    shared_bonus = 0.0
    for _root, sharing in root_of.items():
        if len(sharing) < 2:
            continue
        mods: Set[str] = set()
        for sid in sharing:
            mods.update(indep_sources.get(sid, {}).get("modalities", []))
        if len(mods) >= 2:
            shared_bonus = 0.5

    for sid, roots in [(s, set(indep_sources.get(s, {}).get("roots") or [])) for s in votes]:
        if not roots or votes[sid] == 0.0:
            continue
        others = set(votes) - {sid}
        others_roots: Set[str] = set()
        for o in others:
            others_roots.update(indep_sources.get(o, {}).get("roots") or [])
        if roots <= others_roots:
            votes[sid] = 0.5

    positive = sum(votes.values()) + shared_bonus

    # ── Dimension 3: independent negative validation datasets ──────────
    neg_datasets: Set[str] = set()
    for src in sources:
        if src["source_id"] not in src_set:
            continue
        audit = (src.get("markers", {}).get(type_key, {}) or {}).get("audit", {})
        csv = audit.get("cross_species_validated")
        if isinstance(csv, dict) and csv.get(gene) is False:
            neg_datasets.update(audit.get("expression_validated") or [])

    lam = float(indep.get("lambda_negative", 0.5))
    return max(positive - lam * len(neg_datasets), 0.0)


def compute_consensus_level(source_count: int) -> str:
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


def _normalize_type_key(key: str, type_aliases: Optional[Dict[str, str]] = None) -> str:
    """Map a source-internal type key to the canonical KB name.

    Parameters
    ----------
    key : str
        Raw type key from a source file.
    type_aliases : dict or None
        Mapping of alternative names to canonical names, e.g.
        ``{"Retinal_Ganglion_Cell": "RGC"}``.  When ``None`` no aliasing
        is applied.

    Returns
    -------
    str
        Canonical type key.
    """
    if type_aliases is None:
        return key
    return type_aliases.get(key, key)


# ═══════════════════════════════════════════════════════════════════════
#  Marker merging
# ═══════════════════════════════════════════════════════════════════════


def merge_markers(
    sources: List[Dict[str, Any]],
    type_aliases: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Aggregate markers across sources, grouped by canonical cell type.

    Uses *type_aliases* to normalise type keys; logs a warning when
    a synonym match is made.

    Parameters
    ----------
    sources : list[dict]
        Source dicts as returned by :func:`load_all_sources`.
    type_aliases : dict or None
        Mapping for type-key normalisation.

    Returns
    -------
    dict
        ``{canonical_type: {...}}`` with internal tracking of which sources
        contributed each marker.
    """
    merged: Dict[str, Any] = {}

    for src in sources:
        src_id = src["source_id"]
        src_class = src["meta"].get("class", "")
        src_order = src["meta"].get("order", "")

        # ── Main markers dict ─────────────────────────────────────
        for raw_key, marker_data in src.get("markers", {}).items():
            canonical = _normalize_type_key(raw_key, type_aliases)

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
            for gene in marker_data.get("confirm") or {}:
                _register_marker(entry["confirm"], gene, src_id)

            # Add
            for gene in marker_data.get("add") or {}:
                _register_marker(entry["add"], gene, src_id)

            # Refine
            for gene, refine_data in (marker_data.get("refine") or {}).items():
                entry["refine"].setdefault(gene, []).append(refine_data)

            # Negative markers (union across sources)
            neg = marker_data.get("negative_markers") or []
            if isinstance(neg, list):
                entry["negative_markers"].update(neg)

            # Species from source meta
            entry["species"].update(src["meta"].get("species", []))

        # ── Novel types ───────────────────────────────────────────
        for nt in src.get("novel_types", []):
            if isinstance(nt, str):
                nt_name = nt
            else:
                nt_name = nt.get("name", "") if isinstance(nt, dict) else ""
            if not nt_name:
                continue

            canonical = _normalize_type_key(nt_name, type_aliases)

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
            if isinstance(nt, dict):
                if nt.get("parent"):
                    entry["parent"] = nt["parent"]
                # Novel-type markers go into 'add' (they are novel per definition)
                for gene in nt.get("markers", []):
                    _register_marker(entry["add"], gene, src_id)
                for sp in nt.get("species", []):
                    entry["species"].add(sp)

    return merged


def _register_marker(dest: Dict[str, Dict[str, Any]], gene: str, src_id: str) -> None:
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


def detect_conflicts(
    sources: List[Dict[str, Any]],
    type_aliases: Optional[Dict[str, str]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
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
            canonical = _normalize_type_key(raw_key, type_aliases)
            for tier in ("confirm", "add"):
                for gene in marker_data.get(tier) or {}:
                    gene_type_map.setdefault(gene, {}).setdefault(canonical, set()).add(src_id)

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
                        f"Prefer '{c['type_a']['cell_type']}' ({n_a}/{total} sources vs {n_b})"
                    ),
                }
            )
        elif n_b >= 3 and n_b > n_a:
            resolved.append(
                {
                    **c,
                    "resolution": (
                        f"Prefer '{c['type_b']['cell_type']}' ({n_b}/{total} sources vs {n_a})"
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


def _rule_dedup_key(rule: Dict[str, Any], type_aliases: Optional[Dict[str, str]] = None) -> str:
    """Deterministic string key for a rule (condition + action).

    Normalises the action through *type_aliases* so that rules with
    ``Retinal_Ganglion_Cell`` and ``RGC`` actions are treated as duplicates.
    Includes ``markers_absent`` in the key so that rules that differ only
    by exclusion markers are NOT treated as duplicates.
    """
    condition = rule.get("condition", {})
    action = _normalize_type_key(rule.get("action", ""), type_aliases)
    markers = condition.get("markers_present", {})
    absent = condition.get("markers_absent", [])
    sorted_genes = sorted(markers.keys())
    sorted_absent = sorted(absent) if absent else []
    marker_str = ",".join(f"{g}:{markers[g]}" for g in sorted_genes)
    absent_str = "!" + ",".join(sorted_absent) if sorted_absent else ""
    return f"{action}|{marker_str}{absent_str}"


def merge_rules(
    sources: List[Dict[str, Any]],
    type_aliases: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
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
            rule_normalized["action"] = _normalize_type_key(rule.get("action", ""), type_aliases)
            key = _rule_dedup_key(rule_normalized, type_aliases)
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
    source_independence: Optional[Dict[str, Any]] = None,
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
                    first.get("note", "") + f" (refined by {n_extra + 1} source(s))"
                )
            refine_out[gene] = merged_refine

        # Consensus levels
        consensus_levels: Dict[str, str] = {}
        consensus_effective: Dict[str, float] = {}
        gene_sources: Dict[str, Set[str]] = {}
        for gene, info in type_data.get("confirm", {}).items():
            gene_sources.setdefault(gene, set()).update(info["source_ids"])
        for gene, info in type_data.get("add", {}).items():
            gene_sources.setdefault(gene, set()).update(info["source_ids"])
        for gene, src_set in gene_sources.items():
            if source_independence is not None:
                effective = compute_effective_source_count(
                    gene, type_key, src_set, sources, source_independence
                )
                consensus_effective[gene] = round(effective, 2)
                consensus_levels[gene] = compute_consensus_level(int(round(effective)))
            else:
                consensus_levels[gene] = compute_consensus_level(len(src_set))

        # Resolve class/order — use most common; fall back to sorted list
        classes_list = sorted(type_data.get("classes", set()))
        orders_list = sorted(type_data.get("orders", set()))
        resolved_class = (
            classes_list[0]
            if len(classes_list) == 1
            else (", ".join(classes_list) if classes_list else "")
        )
        resolved_order = (
            orders_list[0]
            if len(orders_list) == 1
            else (", ".join(orders_list) if orders_list else "")
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
            "consensus_effective_counts": consensus_effective,
            "single_source_type": all(v <= 1.0 for v in consensus_effective.values())
            if consensus_effective
            else False,
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
#  Orchestrator: build a complete tissue KB in one call
# ═══════════════════════════════════════════════════════════════════════


def build_tissue_kb(
    sources_dir: str,
    type_aliases: Optional[Dict[str, str]] = None,
    hierarchy_yaml_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a complete Knowledge Base for a tissue from YAML sources.

    This is the primary entry point.  It loads all YAML sources from
    *sources_dir*, merges markers, detects and resolves conflicts,
    merges expert rules, assembles the final KB dict, and optionally
    builds a cell-type hierarchy.

    Parameters
    ----------
    sources_dir : str
        Path to the directory containing per-publication ``{name}.yaml`` files.
    type_aliases : dict or None
        Mapping of alternative type-key names → canonical names.
    hierarchy_yaml_path : str or None
        If provided, builds ``Broad_*`` category entries and computes
        private markers using :mod:`rna.utils.hierarchy`.

    Returns
    -------
    dict
        Complete KB with keys like ``"Rod_Photoreceptor"``,
        ``"expert_rules"``, and ``"_meta"``.
    """
    sources = load_all_sources(sources_dir)
    merged_types = merge_markers(sources, type_aliases)
    merged_rules = merge_rules(sources, type_aliases)
    conflicts = detect_conflicts(sources, type_aliases)
    _resolved = resolve_conflicts(conflicts, sources)

    # Evidence-independence aware consensus (optional).
    # Auto-loads ``_source_independence.yaml`` when present; None = legacy.
    source_independence = None
    if os.path.isdir(sources_dir):
        _indep_path = os.path.join(sources_dir, "_source_independence.yaml")
        source_independence = load_source_independence(_indep_path)
    kb = build_final_kb(merged_types, merged_rules, sources, source_independence)

    # Optionally build hierarchy
    if hierarchy_yaml_path and os.path.isfile(hierarchy_yaml_path):
        from rna.utils.hierarchy import (
            build_hierarchy,
            compute_private_markers,
            load_hierarchy_yaml,
        )

        cfg = load_hierarchy_yaml(hierarchy_yaml_path)
        build_hierarchy(kb, cfg)
        compute_private_markers(kb, cfg)

    # Log summary
    n_types = sum(1 for k in kb if k not in ("expert_rules", "_meta") and not k.startswith("_"))
    logger.info(
        "Built tissue KB: %d types, %d rules, %d conflicts flagged",
        n_types,
        len(kb.get("expert_rules", [])),
        len(_resolved.get("flagged", [])),
    )

    return kb
