"""
``rna/utils/hierarchy.py`` -- Tissue-agnostic framework for building cell-type hierarchy.

Provides functions to load a hierarchy YAML config, derive consensus markers
for broad categories, build synthetic ``Broad_*`` KB entries, backfill parent
fields, compute private markers, and extract incompatible transition pairs.

Usage
-----
>>> cfg = load_hierarchy_yaml("path/to/hierarchy.yaml")
>>> kb = build_hierarchy(kb, cfg)
>>> compute_private_markers(kb, cfg)
>>> pairs = get_incompatible_pairs(kb)

YAML config format
------------------
.. code-block:: yaml

    categories:
        Progenitor:
            label: Progenitor
            fallback_markers: [VSX2, SOX2, HES1]
            members:
                - RPC
                - Proliferating_RPC
                - ...

        Neuron:
            label: Neuron
            fallback_markers: [TUBB3, ELAVL4, SLC17A6]
            members:
                - RGC
                - Cone_Photoreceptor
                - ...

    incompatible_transitions:
        - - "Progenitor"
          - "Neuron"

The only module-level constant that may be tissue-specific is
:data:`CATEGORY_PREFIX` (``"Broad_"`` by default).  It can be overridden via
the config if needed by passing a ``category_prefix`` key.
"""

from __future__ import annotations

import os
from typing import Any

CATEGORY_PREFIX = "Broad_"


# ═══════════════════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════════════════


def load_hierarchy_yaml(path: str) -> dict[str, Any]:
    """Load and validate a hierarchy YAML config.

    Expected YAML structure
    -----------------------
    ``categories`` : dict
        Mapping of category keys to their definition.  Each value must have
        ``label`` (str) and ``members`` (list[str]).
        ``fallback_markers`` (list[str]) is optional.

    ``incompatible_transitions`` : list[list[str]], optional
        Pairs of category keys that are biologically incompatible.

    Parameters
    ----------
    path : str
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed and validated hierarchy config.

    Raises
    ------
    FileNotFoundError
        If the YAML file does not exist.
    ValueError
        If the YAML structure is invalid.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Hierarchy YAML not found: {path}")

    import yaml  # defer import for lazy loading

    with open(path, encoding="utf-8") as f:
        cfg: Any = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        raise ValueError("Hierarchy YAML must contain a top-level mapping")

    if "categories" not in cfg:
        raise ValueError("Hierarchy YAML must have a 'categories' key at the top level")

    categories = cfg["categories"]
    if not isinstance(categories, dict):
        raise ValueError("'categories' must be a dict")

    for cat_key, cat_def in categories.items():
        if not isinstance(cat_def, dict):
            raise ValueError(f"Category {cat_key!r} must be a dict")
        if "label" not in cat_def:
            raise ValueError(f"Category {cat_key!r} is missing required 'label' field")
        if "members" not in cat_def:
            raise ValueError(f"Category {cat_key!r} is missing required 'members' list")
        if not isinstance(cat_def["members"], list):
            raise ValueError(
                f"Category {cat_key!r} 'members' must be a list, "
                f"got {type(cat_def['members']).__name__}"
            )
        if "fallback_markers" in cat_def and not isinstance(cat_def["fallback_markers"], list):
            raise ValueError(f"Category {cat_key!r} 'fallback_markers' must be a list")

    if "incompatible_transitions" in cfg:
        it = cfg["incompatible_transitions"]
        if not isinstance(it, list):
            raise ValueError("'incompatible_transitions' must be a list")
        for i, pair in enumerate(it):
            if not isinstance(pair, list) or len(pair) != 2:
                raise ValueError(f"incompatible_transitions[{i}] must be a 2-element list")

    return cfg


def derive_category_markers(
    kb: dict[str, Any],
    members: list[str],
    fallback: list[str] | None = None,
    n_min: int = 3,
) -> list[str]:
    """Derive consensus markers from member types' confirm markers.

    Algorithm (mirrors the current ``merge.py`` approach):

    1. Count how many member types list each gene as a **confirm** marker.
    2. ``>= 50 %`` member coverage → take top-3.
    3. If fewer than ``n_min``, relax to ``>= 30 %`` → take top-5 → keep
       top-3 (or ``n_min``).
    4. If still fewer than ``n_min``, return the *fallback* list.

    Parameters
    ----------
    kb : dict
        The KB dict.  Each member type must exist as a key and have a
        ``markers.confirm`` sub-dict.
    members : list of str
        Type key strings whose confirm markers will be surveyed.
    fallback : list of str or None
        Marker genes to return when consensus derivation falls short.
    n_min : int
        Minimum number of markers to return (default 3).  Only the final
        threshold check uses this value; the step-2 / step-3 top-N are
        fixed at 3 and 5 respectively, matching the original algorithm.

    Returns
    -------
    list[str]
        Marker gene symbols.  May be empty if *fallback* is ``None`` and
        consensus fails.
    """
    gene_type_count: dict[str, int] = {}
    n_members = len(members)

    for type_key in members:
        entry = kb.get(type_key)
        if entry is None:
            continue
        confirm: dict[str, Any] = entry.get("markers", {}).get("confirm", {})
        for gene in confirm:
            gene_type_count[gene] = gene_type_count.get(gene, 0) + 1

    if not gene_type_count:
        return list(fallback) if fallback else []

    # Tier 1: >= 50 % coverage -> top-3
    half = n_members / 2.0
    candidates = sorted(
        [(cnt, gene) for gene, cnt in gene_type_count.items() if cnt >= half],
        reverse=True,
    )
    markers = [gene for _, gene in candidates[:3]]  # fixed: top-3 per original

    # Tier 2: relax to >= 30 % coverage, take top-5, keep n_min
    if len(markers) < n_min:
        third = n_members * 0.3
        candidates = sorted(
            [(cnt, gene) for gene, cnt in gene_type_count.items() if cnt >= third],
            reverse=True,
        )
        markers = [gene for _, gene in candidates[:5]][:n_min]

    # Tier 3: fallback
    if len(markers) < n_min:
        markers = list(fallback) if fallback else markers

    return markers


def compute_private_markers(kb: dict[str, Any], hierarchy: dict[str, Any]) -> None:
    """Compute and attach ``_private_markers`` to each fine-grained type in the KB.

    A gene is considered *private* to a type if it appears as a **confirm**
    marker in **≤ 2** types across the entire KB (excluding special keys).

    **Mutation**: Each fine-grained type gets a new ``_private_markers``
    entry: ``kb[type_key]["_private_markers"] = [gene, ...]``.

    Types that already have ``_private_markers`` are **not** recomputed
    (they are skipped).

    **Skip list**: ``expert_rules``, ``_meta``, ``_hierarchy``, and any key
    starting with :data:`CATEGORY_PREFIX` (``"Broad_"``).

    Parameters
    ----------
    kb : dict
        The KB dict (mutated in-place).
    hierarchy : dict
        Unused; accepted for signature consistency so callers can pass the
        config dict from :func:`load_hierarchy_yaml` without special-casing.
    """
    # Identify fine-grained types (skip special keys + Broad_*)
    fine_types = [
        k
        for k in kb
        if k not in ("expert_rules", "_meta", "_hierarchy") and not k.startswith(CATEGORY_PREFIX)
    ]

    # Only process types that have not already been computed
    todo = [k for k in fine_types if "_private_markers" not in kb.get(k, {})]
    if not todo:
        return

    # Count how many fine-grained types have each gene as a confirm marker
    gene_type_count: dict[str, int] = {}
    for type_key in fine_types:
        entry = kb.get(type_key)
        if entry is None:
            continue
        confirm: dict[str, Any] = entry.get("markers", {}).get("confirm", {})
        for gene in confirm:
            gene_type_count[gene] = gene_type_count.get(gene, 0) + 1

    # For each type, genes that appear in <= 2 types total are "private"
    for type_key in todo:
        entry = kb.get(type_key)
        if entry is None:
            continue
        confirm = entry.get("markers", {}).get("confirm", {})
        private = [gene for gene in confirm if gene_type_count.get(gene, 0) <= 2]
        entry["_private_markers"] = private


def build_hierarchy(kb: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    """Build ``_hierarchy`` dict, create ``Broad_*`` synthetic KB entries, and
    backfill parent fields.

    Three operations are performed:

    1. **Derive markers** for each category via :func:`derive_category_markers`
       and create a ``_hierarchy`` section in the KB.
    2. **Create Broad_* synthetic entries** in the KB (matching the format
       from the existing ``merge.py`` lines 756-774).
    3. **Backfill parent** fields on all fine-grained member types.

    Parameters
    ----------
    kb : dict
        The KB dict (**MUTATED** in-place -- Broad_* entries added, parents
        set, ``_hierarchy`` attached).
    config : dict
        Parsed hierarchy YAML (output of :func:`load_hierarchy_yaml`).
        Expected keys: ``categories`` (dict) and optionally
        ``incompatible_transitions`` (list of 2-element lists).

    Returns
    -------
    dict
        The same *kb* object (mutated in-place, returned for convenience).
    """
    categories: dict[str, Any] = config.get("categories", {})
    category_prefix = config.get("category_prefix", CATEGORY_PREFIX)

    # Force CATEGORY_PREFIX module-level if not overridden in config.
    # This lets consumers override with config["category_prefix"] while
    # keeping the module default for unmodified callers.
    prefix = category_prefix if category_prefix != CATEGORY_PREFIX else CATEGORY_PREFIX

    # ── 1. Build _hierarchy categories ──────────────────────────────────
    hierarchy_categories: dict[str, Any] = {}
    for cat_key, cat_def in categories.items():
        members: list[str] = cat_def.get("members", [])
        fallback: list[str] | None = cat_def.get("fallback_markers")
        markers = derive_category_markers(kb, members, fallback)

        hierarchy_categories[cat_key] = {
            "label": cat_def.get("label", cat_key),
            "markers": {
                "confirm": {gene: ["_hierarchy"] for gene in markers},
            },
            "members": members,
        }

    # ── 2. Create Broad_* synthetic entries ─────────────────────────────
    for cat_key, cat_def in categories.items():
        broad_key = f"{prefix}{cat_key}"
        members = cat_def.get("members", [])
        fallback = cat_def.get("fallback_markers")
        markers = derive_category_markers(kb, members, fallback)

        kb[broad_key] = {
            "markers": {
                "confirm": {gene: ["_hierarchy"] for gene in markers},
                "add": {},
                "refine": {},
            },
            "negative_markers": [],
            "species": ["Vertebrata"],
            "synonyms": [cat_def.get("label", cat_key)],
            "parent": "",
            "consensus_levels": {},
            "class": "",
            "order": "",
            "classes": [],
            "orders": [],
        }

    # ── 3. Attach _hierarchy ────────────────────────────────────────────
    incompatible = config.get("incompatible_transitions", [])
    _validate_pair_names(incompatible, kb, path="hierarchy.yaml")
    kb["_hierarchy"] = {
        "categories": hierarchy_categories,
        "incompatible_transitions": incompatible,
    }

    # ── 4. Backfill parent fields on member types ───────────────────────
    for type_key in list(kb.keys()):
        if type_key in ("expert_rules", "_meta", "_hierarchy"):
            continue
        if type_key.startswith(prefix):
            continue
        for cat_key, cat_def in categories.items():
            if type_key in cat_def.get("members", []):
                kb[type_key]["parent"] = f"{prefix}{cat_key}"
                break

    return kb


def get_incompatible_pairs(kb: dict[str, Any]) -> list[list[str]]:
    """Extract incompatible transition pairs from ``kb["_hierarchy"]``.

    Parameters
    ----------
    kb : dict
        The KB dict (expected to have a ``_hierarchy`` key, but gracefully
        returns ``[]`` if not).

    Returns
    -------
    list[list[str]]
        List of 2-element category-name pairs.  ``[]`` if no hierarchy or
        no incompatible transitions are defined.
    """
    hierarchy: dict[str, Any] | None = kb.get("_hierarchy")
    if not hierarchy:
        return []
    return hierarchy.get("incompatible_transitions", [])


def _validate_pair_names(
    pairs: list[list[str]], kb: dict[str, Any], path: str = "hierarchy.yaml"
) -> None:
    """Validate that all names in incompatible pairs exist in the KB.

    Raises ValueError with a clear diagnostic listing every name that
    is missing from the KB (with a hint to use canonical names, not
    display names or synonyms).

    Called inside :func:`build_hierarchy` so the mismatch is caught at
    KB build time, not silently ignored at pipeline runtime.
    """
    if not pairs:
        return
    missing: list[str] = []
    for pair in pairs:
        for name in pair:
            if name not in kb:
                missing.append(name)
    if missing:
        unique = sorted(set(missing))
        raise ValueError(
            f"Incompatible-transition names not found in KB ({path}): "
            f"{', '.join(unique)}.  "
            f"Use canonical KB keys (e.g. 'RGC', 'Bipolar_Cell'), "
            f"not display names or synonyms."
        )
