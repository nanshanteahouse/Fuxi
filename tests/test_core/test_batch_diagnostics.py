"""Tests for rna/utils/batch_diagnostics.py — batch effect diagnosis module.

TDD (red phase): all tests fail until the module is implemented at
``rna/utils/batch_diagnostics.py``.

Synthetic AnnData fixtures provide controlled structure so diagnostics
can be validated against known ground truth.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA

# Ensure repo root is in path (redundant with conftest.py but explicit)
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Module under test (will raise ImportError until implemented in T2)
from rna.utils.batch_diagnostics import (  # noqa: E402
    BatchDiagnosisReport,
    ColumnDiagnosis,
    _compute_cramer_v,
    _compute_gini_criterion,
    _compute_purity_one_shot,
    diagnose_batch_candidates,
)

# ═════════════════════════════════════════════════════════════════════════════
# Helper
# ═════════════════════════════════════════════════════════════════════════════


def _make_adata(
    n_cells: int,
    n_genes: int,
    obs_dict: dict[str, pd.Categorical | pd.Series],
    n_pcs: int = 10,
    seed: int = 42,
) -> AnnData:
    """Build a synthetic AnnData with PCA."""
    rng = np.random.RandomState(seed)
    x_raw = rng.randn(n_cells, n_genes).astype(np.float32)
    adata = AnnData(
        X=csr_matrix(rng.negative_binomial(2, 0.5, size=(n_cells, n_genes)).astype(np.float32)),
        obs=pd.DataFrame(obs_dict, index=[f"cell_{i}" for i in range(n_cells)]),
    )
    pca = PCA(n_components=n_pcs, random_state=seed)
    adata.obsm["X_pca"] = pca.fit_transform(x_raw)
    return adata


# ═════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═════════════════════════════════════════════════════════════════════════════


@pytest.fixture
def adata_clean() -> AnnData:
    """500 cells, 20 genes, 4 obs columns.

    - ``batch``: 4 random groups (A/B/C/D) — *diffuse* PCA signal
      (no structured offset added for batch groups).
    - ``biology``: 2 groups (TissueA/TissueB) — *concentrated* PCA signal
      via structured PC1/PC2 offset in the first 5 genes.
    """
    np.random.seed(42)
    n_cells, n_genes = 500, 20

    # Base expression
    x_raw = np.random.randn(n_cells, n_genes).astype(np.float32)

    # Obs columns
    batch_labels = np.random.choice(["A", "B", "C", "D"], n_cells)
    biology_labels = np.random.choice(["TissueA", "TissueB"], n_cells)

    # Add structured biology signal in first 5 genes (concentrated → high Gini)
    x_raw[biology_labels == "TissueA", :5] += 2.0
    x_raw[biology_labels == "TissueB", :5] -= 1.0

    adata = AnnData(
        X=csr_matrix(
            np.random.negative_binomial(2, 0.5, size=(n_cells, n_genes)).astype(np.float32)
        ),
        obs=pd.DataFrame(
            {
                "batch": pd.Categorical(batch_labels),
                "biology": pd.Categorical(biology_labels),
                "sample_id": pd.Categorical(np.random.choice(["S1", "S2", "S3"], n_cells)),
                "n_counts": np.random.poisson(1000, n_cells).astype(float),
            },
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
    )

    # Compute PCA — biology structure is captured in first few components
    pca = PCA(n_components=10, random_state=42)
    adata.obsm["X_pca"] = pca.fit_transform(x_raw)

    return adata


@pytest.fixture
def adata_collinear() -> AnnData:
    """``batch`` and ``biology`` columns are identical → Cramer's V = 1.0."""
    np.random.seed(42)
    labels = np.random.choice(["X", "Y"], 100)
    x_raw = np.random.randn(100, 10).astype(np.float32)
    x_raw[labels == "X", :3] += 1.5
    x_raw[labels == "Y", :3] -= 1.0

    adata = AnnData(
        X=csr_matrix(np.random.negative_binomial(2, 0.5, size=(100, 10)).astype(np.float32)),
        obs=pd.DataFrame(
            {
                "batch": pd.Categorical(labels),
                "biology": pd.Categorical(labels.copy()),
            },
            index=[f"cell_{i}" for i in range(100)],
        ),
    )
    pca = PCA(n_components=5, random_state=42)
    adata.obsm["X_pca"] = pca.fit_transform(x_raw)
    return adata


@pytest.fixture
def adata_empty_obs() -> AnnData:
    """AnnData with zero categorical obs columns (only numeric)."""
    return _make_adata(
        n_cells=50,
        n_genes=10,
        obs_dict={"n_counts": np.random.poisson(1000, 50).astype(float)},
        n_pcs=5,
    )


@pytest.fixture
def adata_all_batch() -> AnnData:
    """All obs columns have random labels (diffuse — no biology signal)."""
    return _make_adata(
        n_cells=200,
        n_genes=15,
        obs_dict={
            "batch1": pd.Categorical(np.random.choice(["A", "B", "C"], 200)),
            "batch2": pd.Categorical(np.random.choice(["X", "Y", "Z"], 200)),
            "donor": pd.Categorical(np.random.choice(["D1", "D2", "D3", "D4", "D5"], 200)),
        },
        n_pcs=10,
    )


@pytest.fixture
def adata_single_col() -> AnnData:
    """Exactly one categorical obs column."""
    return _make_adata(
        n_cells=100,
        n_genes=10,
        obs_dict={
            "condition": pd.Categorical(np.random.choice(["control", "treated"], 100)),
        },
        n_pcs=5,
    )


# ═════════════════════════════════════════════════════════════════════════════
# Tests: diagnose_batch_candidates — global classification
# ═════════════════════════════════════════════════════════════════════════════


def test_diagnose_classifies_batch_correctly(adata_clean: AnnData) -> None:
    """batch column (diffuse signal) → low Gini, judged ``'batch'``."""
    report = diagnose_batch_candidates(adata_clean)
    diag = next(d for d in report.column_diagnoses if d.column == "batch")
    assert diag.gini_criterion < 0.3, f"Expected low Gini for batch, got {diag.gini_criterion:.4f}"
    assert diag.judgment == "batch", f"Expected judgment='batch', got {diag.judgment!r}"


def test_diagnose_classifies_biology_correctly(adata_clean: AnnData) -> None:
    """biology column (concentrated signal) → high Gini, judged ``'biology'``."""
    report = diagnose_batch_candidates(adata_clean)
    diag = next(d for d in report.column_diagnoses if d.column == "biology")
    assert diag.gini_criterion > 0.5, (
        f"Expected high Gini for biology, got {diag.gini_criterion:.4f}"
    )
    assert diag.judgment == "biology", f"Expected judgment='biology', got {diag.judgment!r}"


# ═════════════════════════════════════════════════════════════════════════════
# Tests: _compute_gini_criterion
# ═════════════════════════════════════════════════════════════════════════════


def test_gini_edge_cases() -> None:
    """All-equal R² → Gini ≈ 0; single-dominant → Gini ≈ 1."""
    # All equal → no concentration
    equal = np.ones(20)
    gini_eq = _compute_gini_criterion(equal)
    assert gini_eq < 0.05, f"Expected Gini≈0 for equal array, got {gini_eq:.6f}"

    # Single dominant R² value (all mass in one entry)
    dominant = np.zeros(20)
    dominant[-1] = 1.0  # largest value in ascending order
    gini_dom = _compute_gini_criterion(dominant)
    assert gini_dom > 0.8, f"Expected Gini≈1 for dominant array, got {gini_dom:.6f}"


# ═════════════════════════════════════════════════════════════════════════════
# Tests: _compute_cramer_v
# ═════════════════════════════════════════════════════════════════════════════


def test_cramer_v_perfect_collinearity(adata_collinear: AnnData) -> None:
    """Identical batch and biology columns → Cramer's V = 1.0."""
    v = _compute_cramer_v(
        adata_collinear.obs["batch"],
        adata_collinear.obs["biology"],
    )
    assert v == pytest.approx(1.0, abs=1e-6), (
        f"Expected Cramer's V=1.0 for identical columns, got {v:.6f}"
    )


def test_cramer_v_independent_columns(adata_clean: AnnData) -> None:
    """Independent batch and biology → Cramer's V < 1.0 (no perfect collinearity)."""
    v = _compute_cramer_v(
        adata_clean.obs["batch"],
        adata_clean.obs["biology"],
    )
    # Random 4×2 assignment should not be perfectly collinear
    assert v < 1.0, "Independent columns should not have V=1.0"


# ═════════════════════════════════════════════════════════════════════════════
# Tests: Edge-case data
# ═════════════════════════════════════════════════════════════════════════════


def test_empty_obs_graceful(adata_empty_obs: AnnData) -> None:
    """No categorical obs columns → empty report (no crash)."""
    report = diagnose_batch_candidates(adata_empty_obs)
    assert len(report.column_diagnoses) == 0, (
        "Expected zero diagnoses when no categorical columns exist"
    )
    assert len(report.batch_cols) == 0
    assert len(report.biology_cols) == 0


def test_all_batch_no_biology(adata_all_batch: AnnData) -> None:
    """All columns are random/diffuse → batch_cols populated, biology_cols empty."""
    report = diagnose_batch_candidates(adata_all_batch)
    assert len(report.batch_cols) > 0, "Expected at least one column classified as batch"
    assert len(report.biology_cols) == 0, (
        "Expected no biology columns when all signals are diffuse"
    )


def test_single_unique_value_skipped() -> None:
    """Column with exactly 1 unique value → skipped silently, not errored."""
    n_cells, n_genes = 50, 10
    rng = np.random.RandomState(42)
    x_raw = rng.randn(n_cells, n_genes).astype(np.float32)

    adata = AnnData(
        X=csr_matrix(x_raw.copy()),
        obs=pd.DataFrame(
            {"constant": pd.Categorical(["A"] * n_cells)},
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
    )
    pca = PCA(n_components=5, random_state=42)
    adata.obsm["X_pca"] = pca.fit_transform(x_raw)

    report = diagnose_batch_candidates(adata)
    # The single-value column should be skipped entirely
    col_names = [d.column for d in report.column_diagnoses]
    assert "constant" not in col_names, "Single-value column should be skipped"


def test_moderate_signal_ambiguous() -> None:
    """Column with moderate structured signal → Gini middle band → 'ambiguous'."""
    np.random.seed(123)
    n_cells, n_genes = 200, 15
    x = np.random.randn(n_cells, n_genes).astype(np.float32)

    # Create a column with moderate structured signal (less than biology)
    labels = np.array(["G1"] * (n_cells // 2) + ["G2"] * (n_cells - n_cells // 2))
    np.random.shuffle(labels)
    x[labels == "G1", :3] += 0.3  # moderate offset → Gini lands in the ambiguous band

    adata = AnnData(
        X=csr_matrix(x),
        obs=pd.DataFrame(
            {"moderate": pd.Categorical(labels)}, index=[f"cell_{i}" for i in range(n_cells)]
        ),
    )
    pca = PCA(n_components=5, random_state=123)
    adata.obsm["X_pca"] = pca.fit_transform(x)

    report = diagnose_batch_candidates(adata)
    assert len(report.column_diagnoses) == 1
    diag = report.column_diagnoses[0]
    # The permutation test was removed (a358bfc) — classification is now
    # Gini-threshold only. Moderate signal → ambiguous middle band.
    assert diag.judgment == "ambiguous"
    assert report.ambiguous_cols == ["moderate"]
    # Ambiguous columns take no part in the suggested batch key
    assert report.suggested_batch_key == []


# ═════════════════════════════════════════════════════════════════════════════
# Tests: _compute_purity_one_shot
# ═════════════════════════════════════════════════════════════════════════════


def test_purity_use_rep_parameter(adata_clean: AnnData) -> None:
    """``_compute_purity_one_shot`` with ``use_rep='X_pca'`` vs default."""
    purity_default = _compute_purity_one_shot(adata_clean, "batch")
    purity_pca = _compute_purity_one_shot(adata_clean, "batch", use_rep="X_pca")
    assert isinstance(purity_default, float)
    assert isinstance(purity_pca, float)
    assert 0.0 <= purity_default <= 1.0
    assert 0.0 <= purity_pca <= 1.0
    # Both should produce valid purity scores
    assert purity_default > 0.0, "Purity should be >0 for a valid column"
    assert purity_pca > 0.0, "Purity should be >0 for a valid column"


def test_pca_not_computed_raises() -> None:
    """AnnData without ``X_pca`` in ``.obsm`` → ``ValueError``."""
    np.random.seed(42)
    adata = AnnData(
        X=csr_matrix(np.random.randn(50, 10).astype(np.float32)),
    )
    # No PCA / no embedding available
    with pytest.raises(ValueError, match="X_pca|obsm"):
        _compute_purity_one_shot(adata, "batch")


def test_small_ncells_fallback_purity() -> None:
    """<3 cells → purity defaults to 1.0 (too few for neighbours/clustering)."""
    np.random.seed(42)
    x_raw = np.random.randn(2, 5).astype(np.float32)
    adata = AnnData(
        X=csr_matrix(x_raw.copy()),
        obs=pd.DataFrame(
            {"batch": pd.Categorical(["A", "B"])},
            index=[f"cell_{i}" for i in range(2)],
        ),
    )
    pca = PCA(n_components=2, random_state=42)
    adata.obsm["X_pca"] = pca.fit_transform(x_raw)

    purity = _compute_purity_one_shot(adata, "batch")
    assert purity == pytest.approx(1.0), f"Expected purity=1.0 for <3 cells, got {purity}"


def test_purity_readonly_contract(adata_clean: AnnData) -> None:
    """``_compute_purity_one_shot`` copies internally — original unchanged."""
    original_obsm_keys = set(adata_clean.obsm.keys())
    original_obs_keys = set(adata_clean.obs.columns)

    _compute_purity_one_shot(adata_clean, "batch")

    # No new keys leaked into the original object
    assert set(adata_clean.obsm.keys()) == original_obsm_keys, (
        "New obsm keys were added to the original AnnData"
    )
    assert set(adata_clean.obs.columns) == original_obs_keys, (
        "New obs columns were added to the original AnnData"
    )


# ═════════════════════════════════════════════════════════════════════════════
# Tests: report types and validation
# ═════════════════════════════════════════════════════════════════════════════


def test_diagnose_returns_named_tuple(adata_single_col: AnnData) -> None:
    """``diagnose_batch_candidates`` returns a ``BatchDiagnosisReport``."""
    report = diagnose_batch_candidates(adata_single_col)
    assert isinstance(report, BatchDiagnosisReport)
    # Must have the core fields
    assert hasattr(report, "column_diagnoses")
    assert hasattr(report, "batch_cols")
    assert hasattr(report, "biology_cols")


def test_column_diagnosis_has_expected_fields(adata_clean: AnnData) -> None:
    """Each ``ColumnDiagnosis`` exposes column, gini, judgment, cramer_v, n_unique."""
    report = diagnose_batch_candidates(adata_clean)
    for diag in report.column_diagnoses:
        assert isinstance(diag, ColumnDiagnosis)
        assert isinstance(diag.column, str)
        assert isinstance(diag.gini_criterion, float)
        assert isinstance(diag.judgment, str)
        assert isinstance(diag.n_unique, int) and diag.n_unique > 0
        assert isinstance(diag.cramer_v, dict)
        assert diag.judgment in ("batch", "biology", "ambiguous", "skip")
