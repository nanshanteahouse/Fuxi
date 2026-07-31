"""Tests for ``rna/steps/06_subcluster.py::build_subtype_candidates``.

Plan todo 5 (tiered-subtype-reuse): the helper returns the deduped, sorted
union of three candidate sources for a subcluster cell type:

  (a) L3 members of ``cell_type`` from ``kb["_hierarchy"]["categories"][*]["subtypes"][<L2>]["members"]``
  (b) distinct non-empty ``cell_subtype`` values in ``adata_subset.obs``
      that are NOT in {"unresolved", "N/A"} and NOT equal to ``cell_type``
  (c) distinct subtype names from each cell's parsed ``annot_evidence`` JSON
      ``subtype_candidates[*]["type"]`` (tolerates missing column / malformed JSON)

Must NOT raise when ``kb`` has no ``_hierarchy`` or ``annot_evidence`` is absent.
"""

import importlib.util
import json
import os
import sys

import numpy as np
from anndata import AnnData

# ── Ensure repo root is on sys.path (conftest.py also does this) ──────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 06_subcluster module via file path (step scripts are not
#    importable normally — they only run as subprocesses). ──────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "06_subcluster.py")
_spec = importlib.util.spec_from_file_location("rna.steps._06_subcluster_candidates", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_subtype_candidates = _mod.build_subtype_candidates


# ── Fixtures ────────────────────────────────────────────────────────────


def _kb_with_retina_hierarchy():
    """Minimal KB whose _hierarchy mirrors core/kb/retina/hierarchy.yaml shape."""
    return {
        "_hierarchy": {
            "categories": {
                "Neuron": {
                    "label": "Neuron",
                    "members": ["RGC", "Rod_Photoreceptor"],
                    "subtypes": {
                        "RGC": {
                            "members": ["RGC_Alpha", "RGC_Foxp2", "RGC_Neurod2"],
                        },
                        "Bipolar_Cell": {
                            "members": ["BC1A", "BC2"],
                        },
                    },
                },
            },
            "incompatible_transitions": [],
        }
    }


def _make_adata(n_cells=3, cell_subtype=None, annot_evidence=None):
    """Small synthetic AnnData; columns supplied only when requested."""
    rng = np.random.RandomState(0)
    adata = AnnData(rng.rand(n_cells, 5))
    if cell_subtype is not None:
        adata.obs["cell_subtype"] = cell_subtype
    if annot_evidence is not None:
        adata.obs["annot_evidence"] = annot_evidence
    return adata


def _evidence(*candidates):
    """JSON string for one cell's annot_evidence with the given candidates."""
    return json.dumps({"score": 0.7, "subtype_candidates": list(candidates)})


# ── Tests ───────────────────────────────────────────────────────────────


def test_union_of_all_three_sources_merged_and_deduped():
    """Sources (a) KB members + (b) resolved labels + (c) near-miss names."""
    kb = _kb_with_retina_hierarchy()
    adata = _make_adata(
        cell_subtype=["RGC_Alpha", "unresolved", "RGC_Tbr1"],
        annot_evidence=[
            _evidence({"type": "RGC_W3", "score": 0.5}),
            _evidence({"type": "RGC_Foxp2", "score": 0.4}),  # dup w/ KB member
            _evidence({"type": "RGC_ipRGC", "score": 0.3}),
        ],
    )
    result = build_subtype_candidates(adata, kb, "RGC")
    assert result == [
        "RGC_Alpha",
        "RGC_Foxp2",
        "RGC_Neurod2",
        "RGC_Tbr1",
        "RGC_W3",
        "RGC_ipRGC",
    ]


def test_missing_hierarchy_keeps_only_sources_b_and_c():
    """kb without _hierarchy → (a) drops out, (b)+(c) still present."""
    adata = _make_adata(
        n_cells=3,
        cell_subtype=["RGC_Alpha", "N/A", "RGC_Tbr1"],
        annot_evidence=[
            _evidence({"type": "RGC_W3", "score": 0.5}),
            _evidence({"type": "RGC_Foxp2", "score": 0.4}),
            _evidence({"type": "RGC_Alpha", "score": 0.6}),
        ],
    )
    result = build_subtype_candidates(adata, {"no": "hierarchy"}, "RGC")
    assert result == ["RGC_Alpha", "RGC_Foxp2", "RGC_Tbr1", "RGC_W3"]


def test_missing_annot_evidence_keeps_only_sources_a_and_b():
    """adata without annot_evidence column → (c) drops out, (a)+(b) present."""
    kb = _kb_with_retina_hierarchy()
    adata = _make_adata(cell_subtype=["RGC_Alpha", "unresolved", "RGC_Tbr1"])
    result = build_subtype_candidates(adata, kb, "RGC")
    assert result == ["RGC_Alpha", "RGC_Foxp2", "RGC_Neurod2", "RGC_Tbr1"]


def test_malformed_json_skipped_silently():
    """A non-JSON annot_evidence cell is skipped; no raise."""
    kb = _kb_with_retina_hierarchy()
    adata = _make_adata(
        n_cells=2,
        cell_subtype=["RGC_Alpha", "unresolved"],
        annot_evidence=[
            "{not json",
            _evidence({"type": "RGC_W3", "score": 0.5}),
        ],
    )
    result = build_subtype_candidates(adata, kb, "RGC")
    assert result == ["RGC_Alpha", "RGC_Foxp2", "RGC_Neurod2", "RGC_W3"]


def test_kb_members_included_even_when_resolved_candidates_empty():
    """Source (a) is independent of scoring: empty subtype_candidates in the
    annot_evidence JSON must NOT suppress KB-defined members."""
    kb = _kb_with_retina_hierarchy()
    adata = _make_adata(
        n_cells=2,
        cell_subtype=["RGC_Alpha", "unresolved"],
        annot_evidence=[_evidence(), _evidence()],
    )
    result = build_subtype_candidates(adata, kb, "RGC")
    assert result == ["RGC_Alpha", "RGC_Foxp2", "RGC_Neurod2"]


def test_cell_type_without_subtypes_and_no_other_sources_returns_empty():
    """cell_type absent from hierarchy subtypes + no resolved/near-miss → []."""
    kb = _kb_with_retina_hierarchy()
    adata = _make_adata(n_cells=2, cell_subtype=["Rod_Photoreceptor", "Rod_Photoreceptor"])
    result = build_subtype_candidates(adata, kb, "Rod_Photoreceptor")
    assert result == []
