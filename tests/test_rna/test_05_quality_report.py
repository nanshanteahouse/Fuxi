"""Tests for rna/steps/05_annotate_major.py — _update_quality_report_pass_rate.

Covers all three annotation paths (Unified KB, AI, Score_genes) and
the developing-tissue weighted pass-rate logic, plus filename correctness.
"""

import importlib.util
import json
import os
import sys
import tempfile
from unittest.mock import MagicMock

import numpy as np
import pytest
from anndata import AnnData

# ── Ensure repo root is on sys.path (conftest.py also does this) ──────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 05_annotate_major module via file path ───────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "05_annotate_major.py")
_spec = importlib.util.spec_from_file_location("rna.steps._05_quality_report_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
_update_quality_report_pass_rate = _mod._update_quality_report_pass_rate


# ── Shared helpers ────────────────────────────────────────────────────


def _make_synthetic_adata() -> AnnData:
    """Create a synthetic AnnData with leiden clusters and marker_validation column.

    Returns
    -------
    30 cells x 100 genes, 3 leiden clusters, marker_validation = 20 PASS +
    5 MARGINAL + 5 FAIL.  Optional .raw with 200 genes (RAW_0 … RAW_199).
    """
    rng = np.random.RandomState(42)
    n_cells = 30
    n_genes = 100

    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs["leiden"] = np.array(["0"] * 10 + ["1"] * 10 + ["2"] * 10)

    # 20 PASS + 5 MARGINAL + 5 FAIL (CRITICAL: function early-returns
    # if marker_validation column is missing — see lines 88–93)
    adata.obs["marker_validation"] = ["PASS"] * 20 + ["MARGINAL"] * 5 + ["FAIL"] * 5

    # Optional .raw
    n_raw = 200
    raw_x = rng.poisson(lam=1.0, size=(n_cells, n_raw)).astype(np.float32)
    raw = AnnData(raw_x)
    raw.var_names = [f"RAW_{i}" for i in range(n_raw)]
    adata.raw = raw

    return adata


def _make_cfg(tmp_dir: str, tissue_maturity: str = "") -> MagicMock:
    """Create a minimal Config mock suitable for _update_quality_report_pass_rate."""
    cfg = MagicMock()
    cfg.table_dir = tmp_dir
    cfg.tissue_maturity = tissue_maturity
    return cfg


# ── Tests ─────────────────────────────────────────────────────────────


def test_happy_path_unified() -> None:
    """Unified KB path: pre-existing JSON file is updated with correct pass_rate.

    Given: a quality report already exists (written by
           ``engine._write_quality_report`` with ``pass_rate=0``).
    When:  ``_update_quality_report_pass_rate`` is called.
    Then:  the original keys are preserved and the new pass-rate keys are
           computed from the ``marker_validation`` column.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = _make_synthetic_adata()
        cfg = _make_cfg(tmp_dir)
        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        # Pre-write a JSON file simulating what engine._write_quality_report wrote
        initial = {"pass_rate": 0, "source": "unified_kb"}
        with open(quality_path, "w") as f:
            json.dump(initial, f)

        _update_quality_report_pass_rate(adata, cfg)

        assert os.path.exists(quality_path), "Quality file should still exist"
        with open(quality_path) as f:
            result = json.load(f)

        # Original key preserved
        assert result["source"] == "unified_kb"

        # Newly computed keys
        assert result["strict_pass_rate"] == pytest.approx(0.6667)
        # Non-developing: pass_rate == strict_pass_rate
        assert result["pass_rate"] == pytest.approx(0.6667)
        assert result["pass_cells"] == 20
        assert result["marginal_cells"] == 5


def test_happy_path_ai() -> None:
    """AI path: no pre-existing file → function creates one from scratch.

    Given: no quality report file exists.
    When:  ``_update_quality_report_pass_rate`` is called.
    Then:  the JSON file is created at the expected path with correct values.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = _make_synthetic_adata()
        cfg = _make_cfg(tmp_dir)
        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        assert not os.path.exists(quality_path)

        _update_quality_report_pass_rate(adata, cfg)

        assert os.path.exists(quality_path), "Quality file should be created when none existed"
        with open(quality_path) as f:
            result = json.load(f)

        assert "pass_rate" in result
        assert "strict_pass_rate" in result
        assert result["pass_cells"] == 20
        assert result["marginal_cells"] == 5


def test_happy_path_score_genes() -> None:
    """Score_genes path: no pre-existing file, function computes pass_rate correctly.

    Given: no quality report file exists and all non-MARGINAL cells PASS.
    When:  ``_update_quality_report_pass_rate`` is called.
    Then:  ``strict_pass_rate == pass_rate == 0.8``.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        # 10 cells: 8 PASS + 2 FAIL (no MARGINAL)
        rng = np.random.RandomState(42)
        n_cells = 10
        n_genes = 50
        x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
        adata = AnnData(x)
        adata.obs["marker_validation"] = ["PASS"] * 8 + ["FAIL"] * 2

        cfg = _make_cfg(tmp_dir)
        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        _update_quality_report_pass_rate(adata, cfg)

        assert os.path.exists(quality_path)
        with open(quality_path) as f:
            result = json.load(f)

        # 8/10 = 0.8, and no MARGINAL → pass_rate == strict
        assert result["strict_pass_rate"] == pytest.approx(0.8)
        assert result["pass_rate"] == pytest.approx(0.8)
        assert result["pass_cells"] == 8
        assert result["marginal_cells"] == 0


def test_dev_mode_weighted() -> None:
    """Developing tissue: MARGINAL clusters count at half weight.

    Given: ``tissue_maturity='developing'``, 20 PASS + 5 MARGINAL (25 cells).
    When:  ``_update_quality_report_pass_rate`` computes weighted pass_rate.
    Then:  ``strict_pass_rate = 20/25 = 0.8``,
           ``pass_rate = (20 + 5×0.5)/25 = 22.5/25 = 0.9``.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        # 25 cells: 20 PASS + 5 MARGINAL (no FAIL) — values chosen so that
        # strict=20/25=0.8 and weighted=(20+5×0.5)/25=22.5/25=0.9.
        rng = np.random.RandomState(42)
        n_cells = 25
        n_genes = 50
        x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
        adata = AnnData(x)
        adata.obs["marker_validation"] = ["PASS"] * 20 + ["MARGINAL"] * 5

        cfg = _make_cfg(tmp_dir, tissue_maturity="developing")
        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        _update_quality_report_pass_rate(adata, cfg)

        with open(quality_path) as f:
            result = json.load(f)

        # strict_pass_rate = 20/25 = 0.8
        assert result["strict_pass_rate"] == pytest.approx(0.8)
        # weighted pass_rate = (20 + 5×0.5)/25 = 0.9
        assert result["pass_rate"] == pytest.approx(0.9)
        assert result["pass_cells"] == 20
        assert result["marginal_cells"] == 5


def test_failure_no_file_creates_one() -> None:
    """No existing JSON file → function creates a fresh quality dict.

    Given: no pre-existing quality report file.
    When:  ``_update_quality_report_pass_rate`` is called.
    Then:  all six expected keys are present in the freshly created dict.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = _make_synthetic_adata()
        cfg = _make_cfg(tmp_dir)
        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        assert not os.path.exists(quality_path)

        _update_quality_report_pass_rate(adata, cfg)

        assert os.path.exists(quality_path)
        with open(quality_path) as f:
            result = json.load(f)

        expected_keys = {
            "pass_rate",
            "strict_pass_rate",
            "pass_cells",
            "marginal_cells",
            "kb_blind_spot",
            "recommended_strictness",
        }
        assert set(result.keys()) == expected_keys, (
            f"Expected keys {expected_keys}, got {set(result.keys())}"
        )
        # Sanity-check a few derived fields
        assert result["kb_blind_spot"] is False  # 20/30 ≥ 0.1
        assert result["recommended_strictness"] == "default"  # 20/30 ≥ 0.3


def test_failure_filename_correctness() -> None:
    """Assert written filename is ``05_annotation_quality.json`` (NOT ``_step05.json``).

    Given: ``cfg.table_dir`` is a fresh temp directory.
    When:  ``_update_quality_report_pass_rate`` writes its report.
    Then:  only ``05_annotation_quality.json`` exists — the wrong name
           ``_step05.json`` must not be created.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = _make_synthetic_adata()
        cfg = _make_cfg(tmp_dir)

        _update_quality_report_pass_rate(adata, cfg)

        correct_path = os.path.join(tmp_dir, "05_annotation_quality.json")
        wrong_path = os.path.join(tmp_dir, "_step05.json")

        assert os.path.exists(correct_path), "Correct file 05_annotation_quality.json should exist"
        assert not os.path.exists(wrong_path), "Wrong file _step05.json should NOT be created"


def test_new_engine_fields_survive_pass_rate_update() -> None:
    """Engine-added fields (review_queue/transition_clusters/kb_coverage) survive pass_rate rewrite.



    Given: a quality report already written by ``engine._write_quality_report``

           carrying the D6/D7/D8 fields (review_queue list, transition_clusters

           detail list, kb_coverage dict) and ``pass_rate=0``.

    When:  ``_update_quality_report_pass_rate`` rewrites the pass-rate keys.

    Then:  the engine fields are still present, unchanged — the pass-rate

           update must not clobber them.

    """

    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = _make_synthetic_adata()

        cfg = _make_cfg(tmp_dir)

        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        initial = {
            "pass_rate": 0,
            "source": "unified_kb",
            "review_queue": [
                # old format (task 10 backward compat: no reason)
                {"cluster": "3", "n_tied_types": 4, "top_types": ["RGC", "Cone", "Rod", "NRPC"]},
                # new format: reason-carrying entry
                {"cluster": "5", "reason": "no_canonical_expression"},
            ],
            "transition_clusters": [{"cluster": "0", "pair": "RGC/Amacrine"}],
            "kb_coverage": {
                "annotated_types": ["RGC", "Unknown"],
                "kb_types_unannotated": ["Cone", "Rod"],
                "ghost_endpoints": ["Amacrine"],
            },
        }

        with open(quality_path, "w") as f:
            json.dump(initial, f)

        _update_quality_report_pass_rate(adata, cfg)

        with open(quality_path) as f:
            result = json.load(f)

        # pass-rate keys were still updated on top of the engine fields

        assert result["pass_rate"] == pytest.approx(0.6667)

        assert result["strict_pass_rate"] == pytest.approx(0.6667)

        # D6/D7/D8 fields preserved unchanged

        assert result["review_queue"] == initial["review_queue"]

        assert result["transition_clusters"] == initial["transition_clusters"]

        assert result["kb_coverage"] == initial["kb_coverage"]

        # and the original source key too

        assert result["source"] == "unified_kb"
