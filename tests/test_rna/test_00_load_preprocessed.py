"""Tests for the preprocessed (TSV/CSV + metadata columns) load path.

The branch under test is ``rna/steps/00_load.py`` ``elif data_format ==
"preprocessed"``: auto meta/expr boundary detection + by-name column
assembly of per-sample expression tables into a cells×genes CSR.

Strategy (plan T4, metis G2/G6):
- ``_oracle_preprocessed`` reimplements the CURRENT pre-rework build —
  ``pd.concat(all_dfs, axis=0, ignore_index=True)`` → dense
  ``expr.values.astype(np.float32)`` → ``sp.csr_matrix`` — and the oracle
  tests pin its exact outputs on synthetic fixtures. These run against the
  unchanged module and define the reference behaviour.
- the reworked builder ``_build_preprocessed_sparse`` (chunked row reads,
  by-name column alignment, sparse block assembly — never a full dense
  frame) must reproduce the oracle output for identical fixtures: same var
  order, NaN-filled missing columns, ``_{i}`` barcode suffixes, duplicate
  barcode handling, dtype float32, exact-0.0 dropped, NaN preserved.
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_preprocessed_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Fixture builders / helpers
# ═════════════════════════════════════════════════════════════════════


def _make_frame(rows, cols_meta, cols_expr, data) -> pd.DataFrame:
    """Build a DataFrame from ``rows`` row-major value lists.

    ``data`` is a list of ``rows`` lists, each of length
    ``len(cols_meta) + len(cols_expr)``.
    """
    df = pd.DataFrame(data, columns=list(cols_meta) + list(cols_expr))
    assert len(df) == rows, f"fixture expects {rows} rows, got {len(df)}"
    return df


def _write_frame(path: Path, df: pd.DataFrame, sep: str) -> Path:
    df.to_csv(path, sep=sep, index=False)
    return path


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_preprocessed")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


def _new_build(file_list, sep, meta_cols):
    """Invoke the reworked chunked-sparse builder exactly like the branch does."""
    return _mod._build_preprocessed_sparse(file_list, sep, meta_cols, _make_logger())


def _assert_adata_equal(new: sc.AnnData, oracle: sc.AnnData) -> None:
    """New-path output must match the reference oracle exactly."""
    assert new.shape == oracle.shape
    assert list(new.var_names) == list(oracle.var_names), "var order must match oracle"
    assert list(new.obs_names) == list(oracle.obs_names), "obs_names order must match oracle"
    assert list(new.obs.columns) == list(oracle.obs.columns), "obs column order must match"
    for col in new.obs.columns:
        pd.testing.assert_series_equal(
            new.obs[col], oracle.obs[col], check_names=False, check_dtype=False
        )
    assert new.X.dtype == np.float32, "X dtype must be float32"
    assert np.array_equal(new.X.toarray(), oracle.X.toarray(), equal_nan=True)


GENES = ["G001", "G002", "G003", "G004"]

# Values use fractional floats so the 3-metric classifier sees them as
# expression (small-integer categorical columns classify as metadata).
SINGLE_ROWS = [
    ["BC1", "R", 1.0, 0.0, 5.0, 0.0],
    ["BC2", "M", 0.0, 3.0, 0.0, 2.0],
    ["BC3", "R", 2.0, 0.0, 0.0, 0.5],
    ["BC4", "M", 0.0, 4.0, 0.0, 0.0],
    ["BC5", "R", 0.0, 0.0, 6.0, 1.0],
]


# ═════════════════════════════════════════════════════════════════════
#  Reference oracle — the CURRENT (pre-rework) build
# ═════════════════════════════════════════════════════════════════════


def _oracle_preprocessed(file_list, sep, meta_cols):
    """Reimplements the current dense build: read all → concat → split → csr.

    Mirrors rna/steps/00_load.py L826-859 (pd.concat(axis=0, ignore_index=True),
    positional meta/expr split, duplicate-gene keep-first, barcode ``_{i}``
    suffixing, ``expr.values.astype(np.float32)`` → ``sp.csr_matrix``). The
    meta_columns rename / auto-categorize loops live downstream in ``main()``
    and are untouched by the rework, so they are not reproduced here.
    """
    import warnings

    with warnings.catch_warnings():
        # pandas 2.3 deprecation: concat of an empty frame with a data frame
        # emits FutureWarning; production runs without warnings-as-errors, so the
        # oracle suppresses it to pin the VALUES (the reworked path never concats
        # empty frames — it skips header-only files' empty blocks).
        warnings.simplefilter("ignore", FutureWarning)
        all_dfs = [pd.read_csv(f, sep=sep) for f in file_list]
        combined = pd.concat(all_dfs, axis=0, ignore_index=True)
    meta = combined.iloc[:, :meta_cols].copy()
    expr = combined.iloc[:, meta_cols:]

    if expr.columns.duplicated().any():
        expr = expr.loc[:, ~expr.columns.duplicated(keep="first")]

    barcodes = meta.iloc[:, 0].values.astype(str)
    if not pd.Index(barcodes).is_unique:
        barcodes = [f"{bc}_{i}" for i, bc in enumerate(barcodes)]

    adata = sc.AnnData(X=sp.csr_matrix(expr.values.astype(np.float32)))
    adata.obs_names = barcodes
    adata.var_names = expr.columns.astype(str)
    for col_idx in range(1, meta_cols):
        col_name = str(meta.columns[col_idx]).strip()
        adata.obs[col_name] = meta.iloc[:, col_idx].values
    return adata


class TestOracle:
    """Pin the current dense build on synthetic fixtures (pre-rework oracle)."""

    def test_single_file(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            _make_frame(5, ["barcode", "celltype"], GENES, SINGLE_ROWS),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1)], "\t", 2)

        assert adata.shape == (5, 4)
        assert list(adata.var_names) == GENES
        assert list(adata.obs_names) == ["BC1", "BC2", "BC3", "BC4", "BC5"]
        assert list(adata.obs.columns) == ["celltype"]
        assert list(adata.obs["celltype"].values) == ["R", "M", "R", "M", "R"]
        assert adata.X.dtype == np.float32
        expected = np.array(
            [
                [1.0, 0.0, 5.0, 0.0],
                [0.0, 3.0, 0.0, 2.0],
                [2.0, 0.0, 0.0, 0.5],
                [0.0, 4.0, 0.0, 0.0],
                [0.0, 0.0, 6.0, 1.0],
            ],
            dtype=np.float32,
        )
        assert np.array_equal(adata.X.toarray(), expected, equal_nan=True)

    def test_multi_file_aligned_columns(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            _make_frame(3, ["barcode", "celltype"], GENES, SINGLE_ROWS[:3]),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            _make_frame(2, ["barcode", "celltype"], GENES, SINGLE_ROWS[3:]),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1), str(f2)], "\t", 2)

        assert adata.shape == (5, 4)
        assert list(adata.var_names) == GENES
        assert list(adata.obs_names) == ["BC1", "BC2", "BC3", "BC4", "BC5"]

    def test_missing_columns_fill_nan(self, tmp_path: Path) -> None:
        # file1: 3 genes (no G004); file2: G001/G003/G004 (no G002).
        # pd.concat aligns BY NAME: G002 → NaN for file2 rows, G004 → NaN
        # for file1 rows (union appends G004 at the end).
        f1 = _write_frame(
            tmp_path / "a.tsv",
            pd.DataFrame(
                [
                    ["BC1", "R", 1.0, 0.0, 5.0],
                    ["BC2", "M", 0.0, 3.0, 0.0],
                ],
                columns=["barcode", "celltype", "G001", "G002", "G003"],
            ),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            pd.DataFrame(
                [
                    ["BC9", "M", 9.0, 9.0, 9.0],
                    ["BC10", "R", 8.0, 8.0, 8.0],
                ],
                columns=["barcode", "celltype", "G001", "G003", "G004"],
            ),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1), str(f2)], "\t", 2)

        assert adata.shape == (4, 4)
        assert list(adata.var_names) == GENES
        x = adata.X.toarray()
        assert np.isnan(x[0, 3]) and np.isnan(x[1, 3])  # G004 missing in file1
        assert np.isnan(x[2, 1]) and np.isnan(x[3, 1])  # G002 missing in file2
        assert x[2, 0] == 9.0 and x[3, 0] == 8.0  # G001 aligned by name
        assert x[2, 2] == 9.0 and x[3, 2] == 8.0  # G003 aligned by name
        assert x[2, 3] == 9.0 and x[3, 3] == 8.0  # G004 from file2

    def test_different_gene_order_across_files(self, tmp_path: Path) -> None:
        # Same genes, DIFFERENT order: concat keeps file-0 order, aligns by name.
        f1 = _write_frame(
            tmp_path / "a.tsv",
            _make_frame(2, ["barcode", "celltype"], GENES, SINGLE_ROWS[:2]),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            pd.DataFrame(
                [
                    ["BC9", "M", 7.0, 8.0, 9.0, 10.0],  # G004, G003, G002, G001
                    ["BC10", "R", 6.0, 5.0, 4.0, 3.0],
                ],
                columns=["barcode", "celltype", "G004", "G003", "G002", "G001"],
            ),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1), str(f2)], "\t", 2)

        assert list(adata.var_names) == ["G001", "G002", "G003", "G004"]
        x = adata.X.toarray()
        assert x[2, 0] == 10.0 and x[2, 1] == 9.0 and x[2, 2] == 8.0 and x[2, 3] == 7.0
        assert x[3, 0] == 3.0 and x[3, 1] == 4.0 and x[3, 2] == 5.0 and x[3, 3] == 6.0

    def test_duplicate_header_mangled_by_read_csv(self, tmp_path: Path) -> None:
        # pd.read_csv mangles within-file duplicate columns (G001 → G001.1),
        # so a duplicated gene header surfaces as TWO distinct vars. Both the
        # oracle and the reworked path must agree on this reality.
        f1 = _write_frame(
            tmp_path / "a.tsv",
            pd.DataFrame(
                [
                    ["BC1", "R", 1.0, 0.0, 5.0, 2.0],
                    ["BC2", "M", 0.0, 3.0, 0.0, 4.0],
                ],
                columns=["barcode", "celltype", "G001", "G002", "G003", "G001"],
            ),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1)], "\t", 2)

        assert adata.shape == (2, 4)
        assert list(adata.var_names) == ["G001", "G002", "G003", "G001.1"]
        x = adata.X.toarray()
        assert x[0, 3] == 2.0 and x[1, 3] == 4.0  # mangled G001.1 keeps second values

    def test_duplicate_barcodes_suffix(self, tmp_path: Path) -> None:
        rows = [
            ["BC1", "R", 1.0, 0.0, 5.0, 0.0],
            ["BC1", "M", 0.0, 3.0, 0.0, 2.0],
            ["BC2", "R", 2.0, 0.0, 0.0, 0.0],
            ["BC2", "M", 0.0, 4.0, 0.0, 0.0],
        ]
        f1 = _write_frame(
            tmp_path / "a.tsv", _make_frame(4, ["barcode", "celltype"], GENES, rows), "\t"
        )
        adata = _oracle_preprocessed([str(f1)], "\t", 2)

        # suffix enumeration runs over the COMBINED row order (ignore_index=True)
        assert list(adata.obs_names) == ["BC1_0", "BC1_1", "BC2_2", "BC2_3"]

    def test_duplicate_barcodes_across_files(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            _make_frame(1, ["barcode", "celltype"], GENES, [SINGLE_ROWS[0]]),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            _make_frame(1, ["barcode", "celltype"], GENES, [SINGLE_ROWS[0]]),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1), str(f2)], "\t", 2)

        assert list(adata.obs_names) == ["BC1_0", "BC1_1"]

    def test_zero_cells(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            pd.DataFrame(columns=["barcode", "celltype"] + GENES),
            "\t",
        )
        adata = _oracle_preprocessed([str(f1)], "\t", 2)
        assert adata.shape == (0, 4)
        assert adata.X.dtype == np.float32

    def test_nan_preserved_zero_dropped(self, tmp_path: Path) -> None:
        rows = [
            ["BC1", "R", np.nan, 0.0, 5.0, 0.0],
            ["BC2", "M", 0.0, np.nan, 0.0, 2.0],
        ]
        f1 = _write_frame(
            tmp_path / "a.tsv", _make_frame(2, ["barcode", "celltype"], GENES, rows), "\t"
        )
        adata = _oracle_preprocessed([str(f1)], "\t", 2)
        x = adata.X.toarray()
        assert np.isnan(x[0, 0]) and np.isnan(x[1, 1])
        assert x[0, 1] == 0.0 and x[1, 0] == 0.0 and x[1, 2] == 0.0  # exact zeros dropped


# ═════════════════════════════════════════════════════════════════════
#  Meta/expr boundary classification oracle
# ═════════════════════════════════════════════════════════════════════


def _oracle_classify(sample: pd.DataFrame) -> list[str]:
    """Reimplementation of the current 3-metric column classification."""
    n_sampled = len(sample)
    classifications = []
    for col in sample.columns:
        numeric = pd.to_numeric(sample[col], errors="coerce")
        numeric_ratio = numeric.notna().sum() / n_sampled
        if numeric_ratio < 0.5:
            classifications.append("M")
        else:
            non_na = numeric.dropna()
            zero_frac = (non_na == 0).sum() / len(non_na) if len(non_na) > 0 else 0
            if zero_frac > 0.8:
                classifications.append("E")
                continue
            is_small_int = False
            if len(non_na) > 0:
                if all(v == int(v) for v in non_na):
                    rng = non_na.max() - non_na.min()
                    if rng < 50:
                        is_small_int = True
            if is_small_int:
                classifications.append("M")
            else:
                unique_ratio = numeric.nunique() / n_sampled
                if unique_ratio < 0.5:
                    classifications.append("M")
                else:
                    classifications.append("E")
    return classifications


class TestClassifyOracle:
    def test_meta_and_expr(self) -> None:
        sample = pd.DataFrame(
            {
                "barcode": ["BC1", "BC2", "BC3", "BC4", "BC5"],
                "celltype": ["R", "M", "R", "M", "R"],
                "G001": [1.5, 0.0, 2.5, 0.0, 0.0],  # sparse numeric → E
                "G002": [5.1, 3.2, 8.7, 4.9, 9.3],  # high-cardinality numeric → E
            }
        )
        assert _oracle_classify(sample) == ["M", "M", "E", "E"]

    def test_sparse_expr_after_dense(self) -> None:
        sample = pd.DataFrame(
            {
                "barcode": ["BC1", "BC2", "BC3", "BC4", "BC5"],
                "G001": [1.5, 0.0, 2.5, 0.0, 0.0],
            }
        )
        assert _oracle_classify(sample) == ["M", "E"]

    def test_small_int_categorical_is_meta(self) -> None:
        sample = pd.DataFrame(
            {
                "cluster": [1, 2, 1, 2, 3],  # small-int categorical → M
                "G001": [1.0, 0.0, 2.5, 0.0, 0.0],
            }
        )
        assert _oracle_classify(sample) == ["M", "E"]


# ═════════════════════════════════════════════════════════════════════
#  New path (chunked sparse build) — must equal the oracle
# ═════════════════════════════════════════════════════════════════════


class TestNewPathEqualsOracle:
    def test_single_file(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv", _make_frame(5, ["barcode", "celltype"], GENES, SINGLE_ROWS), "\t"
        )
        files = [str(f1)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_multi_file_aligned_columns(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            _make_frame(3, ["barcode", "celltype"], GENES, SINGLE_ROWS[:3]),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            _make_frame(2, ["barcode", "celltype"], GENES, SINGLE_ROWS[3:]),
            "\t",
        )
        files = [str(f1), str(f2)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_missing_columns_fill_nan(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            pd.DataFrame(
                [
                    ["BC1", "R", 1.0, 0.0, 5.0],
                    ["BC2", "M", 0.0, 3.0, 0.0],
                ],
                columns=["barcode", "celltype", "G001", "G002", "G003"],
            ),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            pd.DataFrame(
                [
                    ["BC9", "M", 9.0, 9.0, 9.0],
                    ["BC10", "R", 8.0, 8.0, 8.0],
                ],
                columns=["barcode", "celltype", "G001", "G003", "G004"],
            ),
            "\t",
        )
        files = [str(f1), str(f2)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_different_gene_order_across_files(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            _make_frame(2, ["barcode", "celltype"], GENES, SINGLE_ROWS[:2]),
            "\t",
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            pd.DataFrame(
                [
                    ["BC9", "M", 7.0, 8.0, 9.0, 10.0],
                    ["BC10", "R", 6.0, 5.0, 4.0, 3.0],
                ],
                columns=["barcode", "celltype", "G004", "G003", "G002", "G001"],
            ),
            "\t",
        )
        files = [str(f1), str(f2)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_duplicate_header_mangled_by_read_csv(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            pd.DataFrame(
                [
                    ["BC1", "R", 1.0, 0.0, 5.0, 2.0],
                    ["BC2", "M", 0.0, 3.0, 0.0, 4.0],
                ],
                columns=["barcode", "celltype", "G001", "G002", "G003", "G001"],
            ),
            "\t",
        )
        files = [str(f1)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_duplicate_barcodes_suffix(self, tmp_path: Path) -> None:
        rows = [
            ["BC1", "R", 1.0, 0.0, 5.0, 0.0],
            ["BC1", "M", 0.0, 3.0, 0.0, 2.0],
            ["BC2", "R", 2.0, 0.0, 0.0, 0.0],
            ["BC2", "M", 0.0, 4.0, 0.0, 0.0],
        ]
        f1 = _write_frame(
            tmp_path / "a.tsv", _make_frame(4, ["barcode", "celltype"], GENES, rows), "\t"
        )
        files = [str(f1)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_zero_cells(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv", pd.DataFrame(columns=["barcode", "celltype"] + GENES), "\t"
        )
        f2 = _write_frame(
            tmp_path / "b.tsv", pd.DataFrame(columns=["barcode", "celltype"] + GENES), "\t"
        )
        files = [str(f1), str(f2)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_empty_then_data(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv", pd.DataFrame(columns=["barcode", "celltype"] + GENES), "\t"
        )
        f2 = _write_frame(
            tmp_path / "b.tsv",
            _make_frame(2, ["barcode", "celltype"], GENES, SINGLE_ROWS[1:3]),
            "\t",
        )
        files = [str(f1), str(f2)]
        _assert_adata_equal(_new_build(files, "\t", 2), _oracle_preprocessed(files, "\t", 2))

    def test_comma_separator(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.csv", _make_frame(5, ["barcode", "celltype"], GENES, SINGLE_ROWS), ","
        )
        files = [str(f1)]
        _assert_adata_equal(_new_build(files, ",", 2), _oracle_preprocessed(files, ",", 2))


class TestNewPathBoundaryDetection:
    """The extracted 3-metric classifier must match its oracle."""

    def test_classifier_matches_oracle(self) -> None:
        sample = pd.DataFrame(
            {
                "barcode": ["BC1", "BC2", "BC3", "BC4", "BC5"],
                "celltype": ["R", "M", "R", "M", "R"],
                "G001": [1.5, 0.0, 2.5, 0.0, 0.0],
                "G002": [5.1, 3.2, 8.7, 4.9, 9.3],
            }
        )
        assert _mod._classify_preprocessed_columns(sample) == _oracle_classify(sample)

    def test_detect_boundary_smoke(self, tmp_path: Path) -> None:
        rows = [
            ["BC1", "R", 1.5, 0.0, 5.0, 0.0],
            ["BC2", "M", 0.0, 3.0, 0.0, 2.0],
            ["BC3", "R", 2.5, 0.0, 0.0, 0.5],
            ["BC4", "M", 0.0, 4.0, 0.0, 0.0],
            ["BC5", "R", 0.0, 0.0, 6.0, 1.0],
        ]
        f1 = _write_frame(
            tmp_path / "a.tsv", _make_frame(5, ["barcode", "celltype"], GENES, rows), "\t"
        )
        meta_cols = _mod._detect_preprocessed_boundary([str(f1)], "\t", _make_logger())
        assert meta_cols == 2

    def test_detect_boundary_all_expr_errors(self, tmp_path: Path) -> None:
        f1 = _write_frame(
            tmp_path / "a.tsv",
            pd.DataFrame(
                {
                    "G001": [1.5, 0.0, 2.5, 0.0, 0.0],
                    "G002": [5.1, 3.2, 8.7, 4.9, 9.3],
                }
            ),
            "\t",
        )
        with pytest.raises(SystemExit):
            _mod._detect_preprocessed_boundary([str(f1)], "\t", _make_logger())


# ═════════════════════════════════════════════════════════════════════
#  Duplicate-gene dedup (keep-first) semantics
# ═════════════════════════════════════════════════════════════════════


class TestDuplicateGeneDedup:
    """keep-first dedup over the union gene columns (expr.columns.duplicated)."""

    def test_keep_first_wins(self) -> None:
        assert _mod._unique_columns_keep_first(["G1", "G2", "G1", "G3", "G2"]) == [
            "G1",
            "G2",
            "G3",
        ]

    def test_no_dups_passthrough(self) -> None:
        assert _mod._unique_columns_keep_first(["G1", "G2", "G3"]) == ["G1", "G2", "G3"]

    def test_logs_warning_on_dups(self) -> None:
        log = _make_logger()
        log.warning = MagicMock()
        _mod._unique_columns_keep_first(["G1", "G1"], log)
        assert log.warning.called


# ═════════════════════════════════════════════════════════════════════
#  Memory-boundedness guard
# ═════════════════════════════════════════════════════════════════════


class TestMemoryBounded:
    """The builder must never materialize the full dense expression frame.

    Primary check is structural (deterministic): the builder reads rows in
    ``chunksize`` blocks and assembles sparse CSR slices with ``sp.vstack``;
    the 3× dense-copy pattern of the old branch (``pd.concat(all_dfs)`` +
    whole-frame ``expr.values.astype``) must be absent from its body.
    """

    def test_structural_chunked_sparse(self) -> None:
        import inspect

        src = inspect.getsource(_mod._build_preprocessed_sparse)
        assert "chunksize" in src, "builder must read files in row-chunks"
        assert "pd.read_csv(" in src, "builder must use pd.read_csv"
        assert "sp.vstack" in src, "builder must assemble sparse CSR blocks"
        assert "pd.concat(all_dfs" not in src, "must not concat all expression frames"
        assert "combined = pd.concat(" not in src, "must not build the combined dense frame"
        # float32 conversion happens per chunk, never on a whole-file frame
        assert ".values.astype(np.float32)" in src

    def test_rss_smoke_wide_fixture(self, tmp_path: Path) -> None:
        """A wide synthetic case runs far under a dense-materialization budget.

        Density is low (~5%) so the sparse path holds nnz only; a full dense
        read of the same shape would need n_cells×n_genes×~20 B. Generous
        threshold to avoid flakiness — this guards against gross regressions
        (e.g. reading all files into one DataFrame), not exact accounting.
        """
        import resource

        n_cells, n_genes = 2000, 8000
        rng = np.random.RandomState(7)
        dense = rng.binomial(1, 0.05, size=(n_cells, n_genes)).astype(np.float32)
        expr = dense * rng.random_sample((n_cells, n_genes)).astype(np.float32)
        df = pd.DataFrame(expr)
        df.insert(0, "celltype", "R")
        df.insert(0, "barcode", [f"BC{i:05d}" for i in range(n_cells)])
        f1 = _write_frame(tmp_path / "wide.tsv", df, "\t")
        del df, expr, dense

        before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        adata = _new_build([str(f1)], "\t", 2)
        after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        delta = (after - before) * 1024  # ru_maxrss is in kB on Linux
        # Dense materialization of this shape would be ≳1.6 GiB; sparse is a
        # few hundred MB. 1.5 GiB is generous enough to never be flaky.
        assert delta < 1.5 * 1024**3, f"peak RSS delta {delta / 2**30:.2f} GiB too high"
        assert adata.shape == (n_cells, n_genes)
        assert adata.X.nnz <= (n_cells * n_genes * 0.05) + 1
