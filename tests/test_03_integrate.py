"""Tests for rna/steps/03_integrate.py — Harmony collinearity hard guard.

T3 (P0-CRITICAL) from cross-batch-critical-fixes plan:
  Prevent Harmony from erasing biology when the batch_key is perfectly
  collinear with a biology column (e.g., sample collinear with genotype).

Tests
-----
- test_T3_collinearity_guard_aborts_harmony — collinear batch_key+biology
  → harmony_skipped set, X_pca_harmony NOT created
- test_T3_collinearity_guard_disabled_harmony_runs — same collinear
  columns but collinearity_guard=False → Harmony runs (X_pca_harmony created)
- test_T3_collinearity_report_none_skips_guard — report is None
  (diagnose=False) → guard no-ops safely, Harmony runs
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from anndata import AnnData

# ── Ensure repo root is on sys.path ──────────────────────────────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 03_integrate module via file path ───────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "03_integrate.py")
_spec = importlib.util.spec_from_file_location(
    "rna.steps._03_integrate_test", _STEP_PATH,
)
assert _spec is not None and _spec.loader is not None, (
    f"Could not load {_STEP_PATH}"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═══════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════


def _make_adata(
    n_cells: int = 30,
    n_genes: int = 100,
    seed: int = 42,
) -> AnnData:
    """Create a minimal AnnData with X_pca pre-computed."""
    rng = np.random.RandomState(seed)
    X = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(X)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs["sample"] = rng.choice(["S1", "S2"], n_cells)
    adata.obs["genotype"] = ["WT"] * (n_cells // 2) + ["KO"] * (n_cells // 2)
    adata.obsm["X_pca"] = rng.randn(n_cells, 50)
    adata.var["highly_variable"] = [True] * n_genes
    adata.uns["pca"] = {"variance_ratio": np.zeros(100)}
    return adata


def _make_cfg(
    collinearity_guard: bool = True,
    diagnose: bool = True,
    use_harmony: bool = True,
) -> MagicMock:
    """Create a Config mock sufficient for the guard + Harmony section."""
    cfg = MagicMock()

    # Harmony settings
    cfg.harmony = MagicMock()
    cfg.harmony.diagnose = diagnose
    cfg.harmony.collinearity_guard = collinearity_guard
    cfg.harmony.use_harmony = use_harmony
    cfg.harmony.batch_key = "sample"
    cfg.harmony.max_iter = 20
    cfg.harmony.diagnose_report = False
    cfg.harmony.gini_batch_threshold = 0.3
    cfg.harmony.gini_biology_threshold = 0.6

    # PCA settings
    cfg.pca = MagicMock()
    cfg.pca.n_pcs_use = 50
    cfg.pca.n_pcs_full = 100

    # Execution settings
    cfg.execution = MagicMock()
    cfg.execution.random_seed = 42
    cfg.execution.use_float32 = False

    # HVG settings
    cfg.hvg = MagicMock()
    cfg.hvg.n_top_genes = 100
    cfg.hvg.flavor = "seurat_v3"
    cfg.hvg.batch_key = None

    # Normalization (disable all heavy work)
    cfg.normalization = MagicMock()
    cfg.normalization.score_cell_cycle = False
    cfg.normalization.use_regress_out = False
    cfg.normalization.regress_out_genes = []
    cfg.normalization.normalize_target_sum = 10000
    cfg.normalization.detect_sex = False

    # Other required fields
    cfg.expression_type = "raw_counts"
    cfg.figure_dir = "/tmp"
    cfg.h5ad_dir = "/tmp"
    cfg.log_dir = "/tmp"
    cfg.qc_h5ad = "/tmp/test.h5ad"
    cfg.modality = "rna"
    cfg.tissue = "test"
    cfg.species = "human"

    return cfg


def _make_report(
    warnings: list[str] | None = None,
) -> MagicMock:
    """Create a BatchDiagnosisReport-like mock with optional warnings."""
    report = MagicMock()
    report.batch_cols = ["sample"]
    report.biology_cols = ["genotype"]
    report.ambiguous_cols = []
    report.warnings = warnings or []
    return report


def _capture_adata_on_save(
    captured: list,
) -> callable:
    """Return a safe_write side-effect that captures the final adata."""
    def _side_effect(adata, *args, **kwargs):
        captured.append(adata)
    return _side_effect


# ═══════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════


def test_T3_collinearity_guard_aborts_harmony() -> None:
    """Collinear batch_key+biology → harmony_skipped set, X_pca_harmony NOT created.

    Given:  diagnose=True, collinearity_guard=True, report has collinearity warning.
    When:   main() runs through the guard section.
    Then:   adata.uns['harmony_skipped'] is set; X_pca_harmony is absent.
    """
    adata = _make_adata()
    cfg = _make_cfg(collinearity_guard=True, diagnose=True)
    report = _make_report(warnings=[
        "Column 'sample' is perfectly collinear with 'genotype' (V=1.0) — redundant column.",
    ])

    captured: list[AnnData] = []
    harmonize_mock = MagicMock()

    with (
        patch.object(
            _mod.argparse.ArgumentParser, "parse_args",
            return_value=argparse.Namespace(config="/tmp/test.yaml"),
        ),
        patch.object(_mod, "resolve_config", return_value=cfg),
        patch.object(_mod, "setup_logger", return_value=MagicMock()),
        patch.object(_mod.sc, "read", return_value=adata),
        patch.object(_mod.sc.pp, "highly_variable_genes"),
        patch.object(_mod.sc.pp, "normalize_total"),
        patch.object(_mod.sc.pp, "log1p"),
        patch.object(_mod.sc.pp, "pca"),
        patch("core.utils.validate_adata", return_value=False),
        patch(
            "rna.utils.batch_diagnostics.diagnose_batch_candidates",
            return_value=report,
        ),
        patch("harmony.harmonize", harmonize_mock),
        patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
    ):
        _mod.main()

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]

    # Guard must have fired
    assert "harmony_skipped" in result.uns, (
        "harmony_skipped should be set when guard fires"
    )
    assert result.uns["harmony_skipped"]["reason"] == "collinearity"
    assert len(result.uns["harmony_skipped"]["warnings"]) > 0

    # X_pca_harmony must NOT be created when guard fires
    assert "X_pca_harmony" not in result.obsm, (
        "X_pca_harmony should NOT be created when guard aborts Harmony"
    )

    # Harmony must NOT have been called
    harmonize_mock.assert_not_called()


def test_T3_collinearity_guard_disabled_harmony_runs() -> None:
    """Same collinear columns but collinearity_guard=False => Harmony runs.

    Given:  diagnose=True, collinearity_guard=False, report has collinearity warning.
    When:   main() runs.
    Then:   X_pca_harmony is created; harmony_skipped is absent.
    """
    adata = _make_adata()
    cfg = _make_cfg(collinearity_guard=False, diagnose=True)
    report = _make_report(warnings=[
        "Column 'sample' is perfectly collinear with 'genotype' (V=1.0) — redundant column.",
    ])

    captured: list[AnnData] = []

    def _fake_harmonize(Z, *args, **kwargs):
        return np.random.RandomState(42).randn(Z.shape[0], 50)

    with (
        patch.object(
            _mod.argparse.ArgumentParser, "parse_args",
            return_value=argparse.Namespace(config="/tmp/test.yaml"),
        ),
        patch.object(_mod, "resolve_config", return_value=cfg),
        patch.object(_mod, "setup_logger", return_value=MagicMock()),
        patch.object(_mod.sc, "read", return_value=adata),
        patch.object(_mod.sc.pp, "highly_variable_genes"),
        patch.object(_mod.sc.pp, "normalize_total"),
        patch.object(_mod.sc.pp, "log1p"),
        patch.object(_mod.sc.pp, "pca"),
        patch("core.utils.validate_adata", return_value=False),
        patch(
            "rna.utils.batch_diagnostics.diagnose_batch_candidates",
            return_value=report,
        ),
        patch("harmony.harmonize", side_effect=_fake_harmonize),
        patch.object(
            _mod.sc.pl, "embedding",
            return_value=None,
        ),
        patch(
            "rna.utils.batch_diagnostics.validate_harmony_preservation",
            return_value={"genotype": 0.95},
        ),
        patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
    ):
        _mod.main()

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]

    # Harmony must have run
    assert "X_pca_harmony" in result.obsm, (
        "X_pca_harmony should be created when guard is disabled"
    )

    # Guard must NOT have fired
    assert "harmony_skipped" not in result.uns, (
        "harmony_skipped should NOT be set when guard is disabled"
    )


def test_T3_collinearity_report_none_skips_guard() -> None:
    """report is None (diagnose=False) => guard no-ops safely, Harmony runs.

    Given:  diagnose=False (report stays None), collinearity_guard=True.
    When:   main() runs.
    Then:   Guard no-ops (no AttributeError); X_pca_harmony is created.
    """
    adata = _make_adata()
    cfg = _make_cfg(collinearity_guard=True, diagnose=False)

    captured: list[AnnData] = []

    def _fake_harmonize(Z, *args, **kwargs):
        return np.random.RandomState(42).randn(Z.shape[0], 50)

    with (
        patch.object(
            _mod.argparse.ArgumentParser, "parse_args",
            return_value=argparse.Namespace(config="/tmp/test.yaml"),
        ),
        patch.object(_mod, "resolve_config", return_value=cfg),
        patch.object(_mod, "setup_logger", return_value=MagicMock()),
        patch.object(_mod.sc, "read", return_value=adata),
        patch.object(_mod.sc.pp, "highly_variable_genes"),
        patch.object(_mod.sc.pp, "normalize_total"),
        patch.object(_mod.sc.pp, "log1p"),
        patch.object(_mod.sc.pp, "pca"),
        patch("core.utils.validate_adata", return_value=False),
        patch("harmony.harmonize", side_effect=_fake_harmonize),
        patch.object(
            _mod.sc.pl, "embedding",
            return_value=None,
        ),
        patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
    ):
        _mod.main()

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]

    # Guard must NOT have fired (report is None)
    assert "harmony_skipped" not in result.uns, (
        "harmony_skipped should NOT be set when report is None"
    )

    # Harmony must have run normally
    assert "X_pca_harmony" in result.obsm, (
        "X_pca_harmony should be created when guard no-ops"
    )
