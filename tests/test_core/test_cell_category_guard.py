"""Tests for the cell-type-aware category guard in ``core.annotation.engine``.

Replicates the inline guard logic (lines 479–516 of engine.py) to verify:

- (×11) ``_broad_parent_map`` maps every known fine type to its correct
  ``Broad_<Category>`` parent.
- (×1)  Every fine type in the mock hierarchy is present in the built map.
- (×1)  ``allows_transitions=False`` rejects a ``Broad_Progenitor`` override.
- (×1)  ``allows_transitions=True``  permits a ``Broad_Progenitor`` override.
- (×1)  ``_n_overrides > 0`` writes a summary ``logger.info`` line.

Known limitation
----------------
The mock hierarchy only contains 4 categories with 13 fine types (the real
``retina/hierarchy.yaml`` has 7+ categories and 50+ types).  This is
intentional — the guard logic is category-count agnostic.
"""

from __future__ import annotations

import logging

import pytest

from core.annotation.engine import CATEGORY_PREFIX

# ── Mock hierarchy (must NOT read real retina/hierarchy.yaml) ──────────────

MOCK_HIERARCHY: dict = {
    "_hierarchy": {
        "categories": {
            "Neuron": {
                "members": [
                    "Bipolar_Cell",
                    "Rod_Photoreceptor",
                    "Cone_Photoreceptor",
                    "RGC",
                    "Amacrine_Cell",
                    "Horizontal_Cell",
                ],
            },
            "Glia": {
                "members": ["Muller_Glia", "Astrocyte"],
            },
            "Non-neural": {
                "members": ["Vascular_Endothelial", "Fetal_RPE", "Pericyte"],
            },
            "Progenitor": {
                "members": ["PRPC", "NRPC"],
            },
        },
    },
}

# Expected mappings for the 11 parametrised cases — excludes Progenitor members.
# fmt: off
FINE_TO_BROAD_CASES: list[tuple[str, str]] = [
    ("Bipolar_Cell",         "Broad_Neuron"),
    ("Rod_Photoreceptor",    "Broad_Neuron"),
    ("Cone_Photoreceptor",   "Broad_Neuron"),
    ("RGC",                  "Broad_Neuron"),
    ("Amacrine_Cell",        "Broad_Neuron"),
    ("Horizontal_Cell",      "Broad_Neuron"),
    ("Muller_Glia",          "Broad_Glia"),
    ("Astrocyte",            "Broad_Glia"),
    ("Vascular_Endothelial", "Broad_Non-neural"),
    ("Fetal_RPE",            "Broad_Non-neural"),
    ("Pericyte",             "Broad_Non-neural"),
]
# fmt: on

# All 13 fine types (Neuron=6 + Glia=2 + Non-neural=3 + Progenitor=2)
ALL_FINE_TYPES: set[str] = {
    "Bipolar_Cell",
    "Rod_Photoreceptor",
    "Cone_Photoreceptor",
    "RGC",
    "Amacrine_Cell",
    "Horizontal_Cell",
    "Muller_Glia",
    "Astrocyte",
    "Vascular_Endothelial",
    "Fetal_RPE",
    "Pericyte",
    "PRPC",
    "NRPC",
}


# ── Helpers that replicate the inline guard logic (engine.py lines 484–516) ─


def _build_broad_parent_map(hierarchy: dict) -> dict[str, str]:
    """Build ``_broad_parent_map`` — mirrors engine.py lines 484–489."""
    parent_map: dict[str, str] = {}
    _hier = hierarchy.get("_hierarchy") or {}
    for _cat_name, _cat_def in (_hier.get("categories") or {}).items():
        _broad_key = f"{CATEGORY_PREFIX}{_cat_name}"
        for _member in _cat_def.get("members") or []:
            parent_map[_member] = _broad_key
    return parent_map


def _apply_category_guard(
    cell_category_map: dict[str, str],
    decision_map: dict[str, object],
    broad_parent_map: dict[str, str],
    allows_transitions: bool = False,
    logger: logging.Logger | None = None,
) -> tuple[dict[str, str], int]:
    """Apply category overrides — mirrors engine.py lines 491–516.

    Returns the (mutated) *cell_category_map* and the override count.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    _n_overrides = 0
    for _cl_str, _decision in decision_map.items():
        _fine_type = getattr(_decision, "cell_type", "") or ""
        _expected_broad = broad_parent_map.get(_fine_type, "")
        _current_broad = cell_category_map.get(_cl_str, "")

        if _expected_broad and _current_broad != _expected_broad:
            # Adult-tissue Broad_Progenitor exclusion guard
            if _expected_broad == f"{CATEGORY_PREFIX}Progenitor" and not allows_transitions:
                continue

            logger.debug(
                "Category guard: cluster %s (%s) %s → %s",
                _cl_str,
                _fine_type,
                _current_broad,
                _expected_broad,
            )
            cell_category_map[_cl_str] = _expected_broad
            _n_overrides += 1

    if _n_overrides:
        logger.info(
            "Category guard: overrode %d/%d cluster categories using fine cell_type parent",
            _n_overrides,
            len(decision_map),
        )

    return cell_category_map, _n_overrides


# ── Mock decision that exposes a ``cell_type`` attribute ───────────────────


class _MockDecision:
    """Minimal stand-in for ``FusionDecision`` — only ``cell_type`` is needed."""

    def __init__(self, cell_type: str) -> None:
        self.cell_type = cell_type


# ═══════════════════════════════════════════════════════════════════════════
# 11 parametrised type-mapping tests
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("fine_type,expected_broad", FINE_TO_BROAD_CASES)
def test_fine_type_to_broad_mapping(fine_type: str, expected_broad: str) -> None:
    """Fine type is mapped to its canonical ``Broad_<Category>`` parent."""
    broad_parent_map = _build_broad_parent_map(MOCK_HIERARCHY)
    assert broad_parent_map[fine_type] == expected_broad


# ═══════════════════════════════════════════════════════════════════════════
# 4 boundary tests
# ═══════════════════════════════════════════════════════════════════════════


def test_broad_parent_map_built_from_hierarchy() -> None:
    """``_broad_parent_map`` contains every fine type from the hierarchy."""
    broad_parent_map = _build_broad_parent_map(MOCK_HIERARCHY)
    assert set(broad_parent_map.keys()) == ALL_FINE_TYPES
    assert len(broad_parent_map) == len(ALL_FINE_TYPES)


def test_adult_tissue_progenitor_excluded() -> None:
    """``allows_transitions=False`` → ``Broad_Progenitor`` override rejected.

    A cluster that the Fisher scorer labelled as ``Broad_Neuron`` but whose
    chosen ``cell_type`` is ``PRPC`` (Progenitor member) must **not** be
    overridden when ``allows_transitions`` is ``False``.
    """
    broad_parent_map = _build_broad_parent_map(MOCK_HIERARCHY)
    decision_map = {"0": _MockDecision("PRPC")}
    cell_category_map: dict[str, str] = {"0": "Broad_Neuron"}

    _, n_overrides = _apply_category_guard(
        cell_category_map=cell_category_map,
        decision_map=decision_map,
        broad_parent_map=broad_parent_map,
        allows_transitions=False,
    )

    assert cell_category_map["0"] == "Broad_Neuron", (
        "Progenitor override must be rejected in adult-tissue mode"
    )
    assert n_overrides == 0


def test_developing_tissue_progenitor_allowed() -> None:
    """``allows_transitions=True`` → ``Broad_Progenitor`` override allowed.

    The same scenario as above but with developmental mode enabled: the
    override should be applied.
    """
    broad_parent_map = _build_broad_parent_map(MOCK_HIERARCHY)
    decision_map = {"0": _MockDecision("PRPC")}
    cell_category_map: dict[str, str] = {"0": "Broad_Neuron"}

    _, n_overrides = _apply_category_guard(
        cell_category_map=cell_category_map,
        decision_map=decision_map,
        broad_parent_map=broad_parent_map,
        allows_transitions=True,
    )

    assert cell_category_map["0"] == "Broad_Progenitor", (
        "Progenitor override must be allowed in developmental mode"
    )
    assert n_overrides == 1


def test_override_count_logged(caplog: pytest.LogCaptureFixture) -> None:
    """``_n_overrides > 0`` triggers a summary ``logger.info`` line."""
    caplog.set_level(logging.INFO)

    broad_parent_map = _build_broad_parent_map(MOCK_HIERARCHY)
    decision_map = {
        "0": _MockDecision("Bipolar_Cell"),
        "1": _MockDecision("Muller_Glia"),
    }
    # Both are misclassified by the Fisher scorer
    cell_category_map: dict[str, str] = {
        "0": "Broad_Glia",  # should be Broad_Neuron
        "1": "Broad_Neuron",  # should be Broad_Glia
    }

    logger = logging.getLogger("test_cell_category_guard")

    _, n_overrides = _apply_category_guard(
        cell_category_map=cell_category_map,
        decision_map=decision_map,
        broad_parent_map=broad_parent_map,
        allows_transitions=False,
        logger=logger,
    )

    assert n_overrides == 2
    assert "overrode 2/2 cluster categories" in caplog.text
