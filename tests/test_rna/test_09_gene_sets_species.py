"""Test _adapt_gene_sets_for_species in the 09_enrichment step module."""

import importlib.util
import os
import sys

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 09_enrichment module via file path ───────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "09_enrichment.py")
_spec = importlib.util.spec_from_file_location(
    "rna.steps._09_enrichment_test",
    _STEP_PATH,
)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

_adapt = _mod._adapt_gene_sets_for_species  # function under test


# ═══════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════


def test_mouse_swaps_human_to_mouse():
    """species='mouse' replaces _Human with _Mouse in gene set names."""
    result = _adapt(["KEGG_2021_Human", "GO_Biological_Process_2023"], "mouse")
    assert result == ["KEGG_2021_Mouse", "GO_Biological_Process_2023"]


def test_human_unchanged():
    """species='human' leaves gene sets unchanged."""
    result = _adapt(["KEGG_2021_Human", "GO_Biological_Process_2023"], "human")
    assert result == ["KEGG_2021_Human", "GO_Biological_Process_2023"]


def test_macaque_keeps_human():
    """species='macaque' (no dedicated libraries) keeps _Human."""
    result = _adapt(["KEGG_2021_Human", "GO_Biological_Process_2023"], "macaque")
    assert result == ["KEGG_2021_Human", "GO_Biological_Process_2023"]


def test_case_insensitive_mouse():
    """species='Mouse' (capitalised) is treated the same as 'mouse'."""
    result = _adapt(["KEGG_2021_Human"], "Mouse")
    assert result == ["KEGG_2021_Mouse"]


def test_empty_gene_sets_returns_empty():
    """Empty input list returns an empty list (does not crash)."""
    result = _adapt([], "mouse")
    assert result == []
