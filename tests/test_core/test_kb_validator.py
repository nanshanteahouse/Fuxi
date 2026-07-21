"""Tests for core/kb_validator.py — Empirical KB marker validation layer.

Uses mock AnnData objects and a controlled KB to verify validation logic
without depending on real pipeline outputs or the full retina KB.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

# ═════════════════════════════════════════════════════════════════════════════
# Mock KB — minimal controlled subset for deterministic tests
# ═════════════════════════════════════════════════════════════════════════════

MOCK_KB = {
    "Rod_Photoreceptor": {
        "markers": {
            "confirm": {"RHO": ["src_1", "src_2"], "GNAT1": ["src_1"]},
            "add": {"SAG": ["src_3"]},
            "refine": {"NRL": {"note": "key TF", "threshold": "", "pmid": "123"}},
        }
    },
    "RGC": {
        "markers": {
            "confirm": {"POU4F1": ["src_1"], "POU4F2": ["src_2"]},
            "add": {"RBPMS": ["src_3"]},
            "refine": {},
        }
    },
}


# ═════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═════════════════════════════════════════════════════════════════════════════


def _make_mock_adata() -> AnnData:
    """Create a small AnnData with controlled expression for 2 cell types.

    Cell type assignments (10 cells total):
      - 5 Rod_Photoreceptor cells: RHO=5/5, GNAT1=3/5, SAG=1/5, NRL=0/5
      - 5 RGC cells: POU4F1=5/5, POU4F2=2/5, RBPMS=0/5

    All Rod markers are OFF in RGC cells and vice versa.
    """
    genes = ["RHO", "GNAT1", "SAG", "NRL", "POU4F1", "POU4F2", "RBPMS"]
    n_cells = 10
    x = np.zeros((n_cells, len(genes)), dtype=np.float32)

    # Rod cells (indices 0-4): RHO, GNAT1, SAG expressed; NRL not
    x[0:5, 0] = 10.0  # RHO   — all 5 express
    x[1:4, 1] = 5.0  # GNAT1 — cells 1,2,3 (3 of 5)
    x[2, 2] = 2.0  # SAG   — cell 2 only (1 of 5)

    # RGC cells (indices 5-9): POU4F1, POU4F2 expressed; RBPMS not
    x[5:10, 4] = 8.0  # POU4F1 — all 5 express
    x[6:8, 5] = 4.0  # POU4F2 — cells 6,7 (2 of 5)

    obs = pd.DataFrame(
        {"cell_type": ["Rod_Photoreceptor"] * 5 + ["RGC"] * 5},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=genes)

    adata = AnnData(x, obs=obs, var=var)
    return adata


@pytest.fixture
def mock_adata() -> AnnData:
    """Return a controlled mock AnnData for validation testing."""
    return _make_mock_adata()


@pytest.fixture
def mock_h5ad_path(mock_adata, tmp_path) -> Path:
    """Write mock AnnData to a temporary .h5ad file."""
    import scanpy as sc

    path = tmp_path / "test_mock.h5ad"
    sc.readwrite.write(str(path), mock_adata)
    return path


# ═════════════════════════════════════════════════════════════════════════════
# Shared helper: create KbValidator with mocked KB, no ontology (fast)
# ═════════════════════════════════════════════════════════════════════════════


def _make_validator(kb=None, use_ontology: bool = False):
    """Create KbValidator with mocked load_kb and controlled ontology flag."""
    if kb is None:
        kb = copy.deepcopy(MOCK_KB)
    with patch("core.kb_validator.load_kb", return_value=kb):
        from core.kb_validator import KbValidator

        return KbValidator(tissue="retina", use_ontology=use_ontology)
    """Create KbValidator with mocked load_kb and controlled ontology flag."""
    with patch("core.kb_validator.load_kb", return_value=kb):
        from core.kb_validator import KbValidator

        return KbValidator(tissue="retina", use_ontology=use_ontology)


# ═════════════════════════════════════════════════════════════════════════════
# Test 1: Core validation — pct_expressed < 0.3 flagged as not_validated
# ═════════════════════════════════════════════════════════════════════════════


class TestKbValidatorCore:
    """Core validation logic — ontology OFF, controlled mock KB."""

    def test_flags_low_expression_as_not_validated(self, mock_adata):
        """Genes with pct_expressed < 0.3 are marked not_validated."""
        validator = _make_validator()
        result = validator.validate(mock_adata, annotation_col="cell_type")

        # NRL: 0/5 Rod cells express → not_validated
        row = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "NRL")]
        assert len(row) == 1
        assert not row.iloc[0]["validated"]
        assert row.iloc[0]["pct_expressed"] == pytest.approx(0.0)

        # SAG: 1/5 = 20% → not_validated
        row = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "SAG")]
        assert len(row) == 1
        assert not row.iloc[0]["validated"]
        assert row.iloc[0]["pct_expressed"] == pytest.approx(0.2)

        # RBPMS: 0/5 RGC cells → not_validated
        row = result[(result.cell_type == "RGC") & (result.gene == "RBPMS")]
        assert len(row) == 1
        assert not row.iloc[0]["validated"]

    def test_flags_high_expression_as_validated(self, mock_adata):
        """Genes with pct_expressed >= 0.3 are marked validated."""
        validator = _make_validator()
        result = validator.validate(mock_adata, annotation_col="cell_type")

        # RHO: 5/5 = 100% → validated
        row = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "RHO")]
        assert len(row) == 1
        assert row.iloc[0]["validated"]
        assert row.iloc[0]["pct_expressed"] == pytest.approx(1.0)

        # GNAT1: 3/5 = 60% → validated
        row = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "GNAT1")]
        assert len(row) == 1
        assert row.iloc[0]["validated"]
        assert row.iloc[0]["pct_expressed"] == pytest.approx(0.6)

        # POU4F1: 5/5 = 100% → validated
        row = result[(result.cell_type == "RGC") & (result.gene == "POU4F1")]
        assert len(row) == 1
        assert row.iloc[0]["validated"]

        # POU4F2: 2/5 = 40% → validated
        row = result[(result.cell_type == "RGC") & (result.gene == "POU4F2")]
        assert len(row) == 1
        assert row.iloc[0]["validated"]
        assert row.iloc[0]["pct_expressed"] == pytest.approx(0.4)

    def test_result_dataframe_columns(self, mock_adata):
        """Result DataFrame has expected columns and correct row count."""
        validator = _make_validator()
        result = validator.validate(mock_adata, annotation_col="cell_type")

        expected_cols = {
            "cell_type",
            "gene",
            "tier",
            "validated",
            "pct_expressed",
            "mean_expression",
        }
        assert set(result.columns) >= expected_cols
        # 2 confirm + 1 add + 1 refine (Rod) + 2 confirm + 1 add (RGC) = 7
        assert len(result) == 7
        assert result["tier"].isin({"confirm", "add", "refine"}).all()

    def test_marker_tiers_propagate_correctly(self, mock_adata):
        """Each marker's tier comes from KB (confirm / add / refine)."""
        validator = _make_validator()
        result = validator.validate(mock_adata, annotation_col="cell_type")

        tier_map = dict(zip(result["gene"], result["tier"]))
        assert tier_map["RHO"] == "confirm"
        assert tier_map["GNAT1"] == "confirm"
        assert tier_map["SAG"] == "add"
        assert tier_map["NRL"] == "refine"
        assert tier_map["POU4F1"] == "confirm"
        assert tier_map["POU4F2"] == "confirm"
        assert tier_map["RBPMS"] == "add"


# ═════════════════════════════════════════════════════════════════════════════
# Test 2: Missing genes handled gracefully (no crash)
# ═════════════════════════════════════════════════════════════════════════════

KB_WITH_MISSING = copy.deepcopy(MOCK_KB)
KB_WITH_MISSING["Rod_Photoreceptor"]["markers"]["confirm"]["MISSING_GENE"] = ["src_x"]


class TestMissingGenes:
    """Genes in KB but absent from dataset should not crash."""

    def test_missing_gene_in_kb_no_crash(self, mock_adata):
        """A KB marker not in adata.var_names is recorded with 0 expression."""
        validator = _make_validator(kb=KB_WITH_MISSING)
        result = validator.validate(mock_adata, annotation_col="cell_type")

        row = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "MISSING_GENE")]
        assert len(row) == 1
        assert not row.iloc[0]["validated"]
        assert row.iloc[0]["pct_expressed"] == 0.0
        assert row.iloc[0]["mean_expression"] == 0.0


# ═════════════════════════════════════════════════════════════════════════════
# Test 3: Case-insensitive gene matching (cross-species)
# ═════════════════════════════════════════════════════════════════════════════


class TestCaseInsensitiveGeneMatching:
    """Cross-species: mouse (lowercase) genes match human (UPPERCASE) KB."""

    @staticmethod
    def _make_lowercase_adata() -> AnnData:
        """Same expression as mock but with lowercase var_names (mouse convention)."""
        genes_upper = ["RHO", "GNAT1", "SAG", "NRL", "POU4F1", "POU4F2", "RBPMS"]
        genes_lower = [g.lower() for g in genes_upper]
        n_cells = 10
        x = np.zeros((n_cells, len(genes_lower)), dtype=np.float32)

        x[0:5, 0] = 10.0  # rho — all Rod express
        x[1:4, 1] = 5.0  # gnat1
        x[2, 2] = 2.0  # sag
        x[5:10, 4] = 8.0  # pou4f1 — all RGC
        x[6:8, 5] = 4.0  # pou4f2

        obs = pd.DataFrame(
            {"cell_type": ["Rod_Photoreceptor"] * 5 + ["RGC"] * 5},
            index=[f"cell_{i}" for i in range(n_cells)],
        )
        var = pd.DataFrame(index=genes_lower)
        return AnnData(x, obs=obs, var=var)

    def test_case_insensitive_lowercase_mouse_genes(self):
        """Lowercase mouse gene names match uppercase KB markers."""
        adata = self._make_lowercase_adata()
        validator = _make_validator()
        result = validator.validate(adata, annotation_col="cell_type")

        # RHO (KB uppercase) should match "rho" (var lowercase)
        row = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "RHO")]
        assert len(row) == 1
        assert row.iloc[0]["validated"]  # 100% expression
        assert row.iloc[0]["pct_expressed"] == pytest.approx(1.0)

        # POU4F2 matches pou4f2
        row = result[(result.cell_type == "RGC") & (result.gene == "POU4F2")]
        assert len(row) == 1
        assert row.iloc[0]["validated"]  # 40% → validated
        assert row.iloc[0]["pct_expressed"] == pytest.approx(0.4)


# ═════════════════════════════════════════════════════════════════════════════
# Test 4: CLI entry point (fast: in-process, no subprocess)
# ═════════════════════════════════════════════════════════════════════════════


class TestCLI:
    """CLI: --h5ad, --annotation, --output, --no-ontology via main()."""

    def test_cli_outputs_csv(self, tmp_path):
        """CLI with --no-ontology and --output writes CSV."""
        import scanpy as sc

        adata = _make_mock_adata()
        h5ad_path = tmp_path / "test.h5ad"
        sc.readwrite.write(str(h5ad_path), adata)

        output_path = tmp_path / "out.csv"
        args = [
            "kb_validator.py",
            "--h5ad",
            str(h5ad_path),
            "--annotation",
            "cell_type",
            "--no-ontology",
            "--output",
            str(output_path),
        ]
        _orig = sys.argv[:]
        sys.argv = args
        try:
            import core.kb_validator as kbm

            with patch.object(kbm, "load_kb", return_value=copy.deepcopy(MOCK_KB)):
                kbm.main()
        finally:
            sys.argv = _orig

        assert output_path.exists()
        df = pd.read_csv(output_path)
        assert len(df) > 0
        assert "cell_type" in df.columns
        assert "gene" in df.columns
        assert "validated" in df.columns

    def test_cli_prints_summary(self, tmp_path, capsys):
        """CLI without --output prints summary."""
        import scanpy as sc

        adata = _make_mock_adata()
        h5ad_path = tmp_path / "test2.h5ad"
        sc.readwrite.write(str(h5ad_path), adata)

        args = [
            "kb_validator.py",
            "--h5ad",
            str(h5ad_path),
            "--annotation",
            "cell_type",
            "--no-ontology",
        ]
        _orig = sys.argv[:]
        sys.argv = args
        try:
            import core.kb_validator as kbm

            with patch.object(kbm, "load_kb", return_value=copy.deepcopy(MOCK_KB)):
                kbm.main()
        finally:
            sys.argv = _orig

        captured = capsys.readouterr()
        assert "confirm" in captured.out.lower() or "validated" in captured.out.lower()


# ═════════════════════════════════════════════════════════════════════════════
# Additional edge-case tests
# ═════════════════════════════════════════════════════════════════════════════


class TestEdgeCases:
    """Edge cases: empty adata, no matching KB types, missing annotation col."""

    def test_empty_adata_returns_empty_dataframe(self):
        """Empty AnnData should not crash."""
        validator = _make_validator()
        adata = AnnData(
            X=np.zeros((0, 3), dtype=np.float32),
            obs=pd.DataFrame({"cell_type": []}, index=pd.Index([], dtype=str)),
        )
        result = validator.validate(adata, annotation_col="cell_type")
        assert len(result) == 0
        assert isinstance(result, pd.DataFrame)
        assert set(result.columns) >= {"cell_type", "gene", "tier", "validated"}

    def test_no_matching_kb_types_returns_empty(self):
        """adata with cell types not in KB returns empty DataFrame."""
        validator = _make_validator()
        adata = AnnData(
            X=np.zeros((3, 2), dtype=np.float32),
            obs=pd.DataFrame(
                {"cell_type": ["Unknown_Type"] * 3}, index=pd.Index(["c1", "c2", "c3"])
            ),
        )
        result = validator.validate(adata, annotation_col="cell_type")
        assert len(result) == 0

    def test_missing_annotation_col_raises_keyerror(self, mock_adata):
        """Non-existent annotation column raises KeyError."""
        validator = _make_validator()
        with pytest.raises(KeyError):
            validator.validate(mock_adata, annotation_col="nonexistent_col")


# ═════════════════════════════════════════════════════════════════════════════
# Integration: validate against real retina KB (smoke test)
# ═════════════════════════════════════════════════════════════════════════════


class TestRealKBIntegration:
    """Minimum smoke test using the real retina KB (no network, full ontology)."""

    def test_real_kb_loads_and_validates_with_ontology(self, mock_adata):
        """KbValidator with real retina KB and StandardOntology runs successfully."""
        from core.kb_validator import KbValidator

        validator = KbValidator(tissue="retina", use_ontology=True)
        result = validator.validate(mock_adata, annotation_col="cell_type")
        assert len(result) > 0
        assert "validated" in result.columns

        # RHO (5/5 Rod cells) should validate with real KB
        rho_rows = result[(result.cell_type == "Rod_Photoreceptor") & (result.gene == "RHO")]
        if len(rho_rows) > 0:
            assert rho_rows.iloc[0]["validated"]
