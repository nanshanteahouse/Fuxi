"""Unit tests for ``rna/utils/sex_detection.py`` — :func:`detect_sex`.

Test scenarios
--------------
1. Explicit ``sex`` column exists → early return, no ``predicted_sex`` added.
2. ``gender`` column alias → same early return.
3. Both ``sex`` and ``gender`` → ``sex`` takes priority.
4. ``adata.raw`` is ``None`` → warning, early return.
5. No sex-linked genes in ``adata.raw.var_names`` → warning, early return.
6. Mixed-sex mouse dataset (Xist + Eif2s3y) → correct Female/Male/Unknown labels.
7. Ambiguous cells (both Xist and Eif2s3y positive) → ``Ambiguous`` label.
8. Human panel (XIST + RPS4Y1) → species autodetection.
9. Log messages verified via ``caplog``.
10. CSV report written to ``CFG.table_dir``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from rna.utils.sex_detection import detect_sex

# ── helpers ─────────────────────────────────────────────────────────────────


class _MockCFG:
    """Minimal CFG stub — only provides a writable ``table_dir``."""

    def __init__(self, table_dir: str | Path) -> None:
        self.table_dir = str(table_dir)


def _make_test_adata(
    n_cells: int = 100,
    n_genes: int = 1000,
    *,
    with_raw: bool = True,
    include_xist: bool = True,
    include_y: bool = True,
) -> sc.AnnData:
    """Build a synthetic AnnData for sex-detection tests.

    The ``.raw`` layer (if created) is **zero-filled** except for explicit
    positive cells, giving deterministic, repeatable results.

    Parameters
    ----------
    n_cells
        Number of cells (rows).
    n_genes
        Number of genes (columns).
    with_raw
        Whether to populate ``adata.raw``.
    include_xist
        Place ``Xist`` at var index 0 if ``True``.
    include_y
        Place ``Eif2s3y`` at var index 1 if ``True``.

    Returns
    -------
    AnnData
        A fresh AnnData object suitable for ``detect_sex``.
    """
    rng = np.random.RandomState(42)
    X = rng.poisson(0.5, size=(n_cells, n_genes)).astype(np.float32)
    var_names = [f"Gene_{i}" for i in range(n_genes)]

    if include_xist:
        var_names[0] = "Xist"
    if include_y:
        var_names[1] = "Eif2s3y"

    index = pd.Index(var_names)
    adata = sc.AnnData(X=X, var=pd.DataFrame(index=index))

    if with_raw:
        # Clean raw layer — all zeros except where we seed expression
        raw_X = np.zeros((n_cells, n_genes), dtype=np.float32)
        if include_xist:
            raw_X[:30, 0] = 1.0  # cells  0…29   Xist+
        if include_y:
            raw_X[50:80, 1] = 1.0  # cells 50…79   Eif2s3y+
        raw_adata = sc.AnnData(X=raw_X, var=pd.DataFrame(index=index))
        adata.raw = raw_adata

    return adata


# ── existing-column tests ──────────────────────────────────────────────────


class TestExistingColumns:
    """detect_sex should return early when 'sex' or 'gender' is present."""

    def test_explicit_sex_column(self, tmp_path: Path) -> None:
        """'sex' column exists → early return, no predicted_sex added."""
        adata = _make_test_adata()
        adata.obs["sex"] = ["Male"] * 50 + ["Female"] * 50
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "predicted_sex" not in adata.obs

    def test_gender_column_alias(self, tmp_path: Path) -> None:
        """'gender' column exists → same early return behaviour."""
        adata = _make_test_adata()
        adata.obs["gender"] = ["M"] * 50 + ["F"] * 50
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "predicted_sex" not in adata.obs

    def test_sex_preferred_over_gender(self, tmp_path: Path) -> None:
        """When both 'sex' and 'gender' exist, 'sex' takes priority (first
        match in the ``("sex", "gender")`` tuple)."""
        adata = _make_test_adata()
        adata.obs["sex"] = ["Male"] * 100
        adata.obs["gender"] = ["F"] * 100
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "predicted_sex" not in adata.obs


# ── missing-data tests ─────────────────────────────────────────────────────


class TestMissingData:
    """detect_sex should return early when prerequisites are not met."""

    def test_no_raw(self, tmp_path: Path) -> None:
        """adata.raw is None → warning log, return, no predicted_sex."""
        adata = _make_test_adata(with_raw=False)
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "predicted_sex" not in adata.obs

    def test_no_sex_genes(self, tmp_path: Path) -> None:
        """No sex-linked genes in raw → warning log, return, no predicted_sex."""
        adata = _make_test_adata(
            with_raw=True, include_xist=False, include_y=False
        )
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "predicted_sex" not in adata.obs


# ── prediction tests ──────────────────────────────────────────────────────


class TestPrediction:
    """detect_sex predictions with various expression patterns."""

    def test_mixed_sex_mouse(self, tmp_path: Path) -> None:
        """Mixed-sex dataset with Xist⁺ (female) and Eif2s3y⁺ (male) cells.

        Expected layout (100 cells):
          • cells  0…29  Xist⁺          → Female   (30)
          • cells 30…49  neither        → Unknown  (20)
          • cells 50…79  Eif2s3y⁺       → Male     (30)
          • cells 80…99  neither        → Unknown  (20)
        """
        adata = _make_test_adata(n_cells=100, with_raw=True)
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "predicted_sex" in adata.obs
        assert (adata.obs["predicted_sex"] == "Female").sum() == 30
        assert (adata.obs["predicted_sex"] == "Male").sum() == 30
        assert (adata.obs["predicted_sex"] == "Unknown").sum() == 40
        assert (adata.obs["predicted_sex"] == "Ambiguous").sum() == 0

    def test_ambiguous_cells(self, tmp_path: Path) -> None:
        """Cells expressing both female and male markers → ``Ambiguous``."""
        n = 100
        var_names = [f"Gene_{i}" for i in range(100)]
        var_names[0] = "Xist"
        var_names[1] = "Eif2s3y"
        index = pd.Index(var_names)

        raw_X = np.zeros((n, 100), dtype=np.float32)
        raw_X[:10, 0] = 1.0  # cells  0…9  Xist⁺       → Female
        raw_X[10:20, 1] = 1.0  # cells 10…19 Eif2s3y⁺   → Male
        raw_X[20:30, 0] = 1.0  # cells 20…29 both⁺      → Ambiguous
        raw_X[20:30, 1] = 1.0

        adata = sc.AnnData(
            X=np.zeros((n, 100), dtype=np.float32),
            var=pd.DataFrame(index=index),
        )
        adata.raw = sc.AnnData(X=raw_X, var=pd.DataFrame(index=index))

        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert (adata.obs["predicted_sex"] == "Female").sum() == 10
        assert (adata.obs["predicted_sex"] == "Male").sum() == 10
        assert (adata.obs["predicted_sex"] == "Ambiguous").sum() == 10
        assert (adata.obs["predicted_sex"] == "Unknown").sum() == 70

    def test_human_panel(self, tmp_path: Path) -> None:
        """Human gene names (XIST + RPS4Y1) detected and labeled correctly."""
        var_names = [f"Gene_{i}" for i in range(100)]
        var_names[0] = "XIST"
        var_names[1] = "RPS4Y1"
        index = pd.Index(var_names)

        raw_X = np.zeros((50, 100), dtype=np.float32)
        raw_X[:15, 0] = 1.0  # cells  0…14 XIST⁺    → Female
        raw_X[25:40, 1] = 1.0  # cells 25…39 RPS4Y1⁺  → Male

        adata = sc.AnnData(
            X=np.zeros((50, 100), dtype=np.float32),
            var=pd.DataFrame(index=index),
        )
        adata.raw = sc.AnnData(X=raw_X, var=pd.DataFrame(index=index))

        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert (adata.obs["predicted_sex"] == "Female").sum() == 15
        assert (adata.obs["predicted_sex"] == "Male").sum() == 15


# ── logging tests ──────────────────────────────────────────────────────────


class TestLogging:
    """Verify key log messages emitted by detect_sex."""

    def test_logs_existing_column(
        self, caplog: pytest.LogCaptureFixture, tmp_path: Path
    ) -> None:
        caplog.set_level(logging.INFO)
        adata = _make_test_adata()
        adata.obs["sex"] = ["Male"] * 50 + ["Female"] * 50
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "already present" in caplog.text

    def test_logs_no_raw_warning(
        self, caplog: pytest.LogCaptureFixture, tmp_path: Path
    ) -> None:
        caplog.set_level(logging.WARNING)
        adata = _make_test_adata(with_raw=False)
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "No raw data" in caplog.text

    def test_logs_no_sex_genes_warning(
        self, caplog: pytest.LogCaptureFixture, tmp_path: Path
    ) -> None:
        caplog.set_level(logging.WARNING)
        adata = _make_test_adata(
            with_raw=True, include_xist=False, include_y=False
        )
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "No sex-linked genes" in caplog.text

    def test_logs_mouse_panel(
        self, caplog: pytest.LogCaptureFixture, tmp_path: Path
    ) -> None:
        caplog.set_level(logging.INFO)
        adata = _make_test_adata()
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "mouse panel" in caplog.text

    def test_logs_human_panel(
        self, caplog: pytest.LogCaptureFixture, tmp_path: Path
    ) -> None:
        caplog.set_level(logging.INFO)
        var_names = [f"Gene_{i}" for i in range(50)]
        var_names[0] = "XIST"
        var_names[1] = "DDX3Y"
        index = pd.Index(var_names)
        raw_X = np.zeros((50, 50), dtype=np.float32)
        raw_X[:10, 0] = 1.0
        raw_X[20:30, 1] = 1.0
        adata = sc.AnnData(
            X=np.zeros((50, 50), dtype=np.float32),
            var=pd.DataFrame(index=index),
        )
        adata.raw = sc.AnnData(X=raw_X, var=pd.DataFrame(index=index))
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        assert "human panel" in caplog.text


# ── CSV report test ────────────────────────────────────────────────────────


class TestCSVReport:
    """The ``sex_report.csv`` should be written when predictions are made."""

    def test_report_saved(self, tmp_path: Path) -> None:
        adata = _make_test_adata()
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        report_path = tmp_path / "sex_report.csv"
        assert report_path.exists()
        df = pd.read_csv(report_path)
        assert "barcode" in df.columns
        assert "predicted_sex" in df.columns
        assert len(df) == adata.n_obs

    def test_no_report_when_sex_exists(self, tmp_path: Path) -> None:
        """Early return paths should NOT write the report."""
        adata = _make_test_adata()
        adata.obs["sex"] = ["Male"] * 50 + ["Female"] * 50
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        report_path = tmp_path / "sex_report.csv"
        assert not report_path.exists()

    def test_no_report_when_raw_missing(self, tmp_path: Path) -> None:
        adata = _make_test_adata(with_raw=False)
        detect_sex(adata, _MockCFG(tmp_path), logging.getLogger("test"))
        report_path = tmp_path / "sex_report.csv"
        assert not report_path.exists()
