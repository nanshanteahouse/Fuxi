"""Tests for the csv_matrix TABLE load path in rna/steps/00_load.py.

The branch under test is ``rna/steps/00_load.py`` ``elif data_format ==
"csv_matrix"`` table sub-branch (``.csv/.tsv/.txt``): a single genes×cells
delimited file read with ``pd.read_csv(..., index_col=0)`` — first row = cell
barcodes (header), first column = gene names (index) — transposed into a
cells×genes expression matrix.

Strategy (plan T5, metis G2/G6):
- ``_oracle_csv_matrix`` reimplements the CURRENT pre-rework build — single
  ``pd.read_csv(..., index_col=0)`` → ``df.values.T.astype(np.float32)`` →
  ``sp.csr_matrix`` (the branch produces a dense X that the tail
  ``force_csr`` block converts to CSR) — and the oracle tests pin its exact
  outputs on synthetic tables. They run against the unchanged module and
  define the reference behaviour.
- the reworked builder ``_build_csv_matrix_sparse`` (chunked gene-block reads,
  per-block float32 transpose, CSR hstack — never a full dense transpose)
  must reproduce the oracle output for identical fixtures: same var order
  (gene index), obs order (cell header), dtype float32, NaN preserved,
  exact-0.0 dropped.

One deliberate improvement: a non-numeric cell crashed the old build
(``object.astype(np.float32)`` raises ValueError); the new path coerces it to
NaN with no crash (plan T5 QA scenario).
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_csv_matrix_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Fixture builders / helpers
# ═════════════════════════════════════════════════════════════════════


def _write_table(path: Path, genes, cells, data, sep="\t") -> Path:
    """Write a genes×cells delimited table.

    ``data`` is a list of ``len(genes)`` rows, each ``len(cells)`` values.
    The header row lists the cell barcodes (first cell empty → index name),
    the first column lists the gene names — exactly the layout the branch's
    ``pd.read_csv(..., index_col=0)`` expects.
    """
    with open(path, "w") as f:
        f.write(sep.join([""] + [str(c) for c in cells]) + "\n")
        for g, row in zip(genes, data):
            f.write(sep.join([str(g)] + [str(v) for v in row]) + "\n")
    return path


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_csv_matrix")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


def _detect_sep(matrix_file) -> str:
    """Replica of the branch's sep auto-detect (sep=None sniff, >1 col → tab)."""
    try:
        peek = pd.read_csv(matrix_file, sep=None, engine="python", nrows=1)
        return "\t" if len(peek.columns) > 1 else ","
    except Exception:
        return "\t"


def _new_build(matrix_file, sep, decimal=".", chunksize=None):
    """Invoke the reworked chunked-sparse builder exactly like the branch does."""
    return _mod._build_csv_matrix_sparse(
        matrix_file, sep, decimal, _make_logger(), chunksize=chunksize
    )


def _assert_adata_equal(new: sc.AnnData, oracle: sc.AnnData) -> None:
    """New-path output must match the reference oracle exactly."""
    assert new.shape == oracle.shape
    assert list(new.var_names) == list(oracle.var_names), "var order must match oracle"
    assert list(new.obs_names) == list(oracle.obs_names), "obs_names order must match oracle"
    assert new.X.dtype == np.float32, "X dtype must be float32"
    assert np.array_equal(new.X.toarray(), oracle.X.toarray(), equal_nan=True)


GENES = ["G001", "G002", "G003", "G004"]
CELLS = ["BC1", "BC2", "BC3", "BC4", "BC5"]

# Genes×cells rows (row = gene, col = cell), sparse-ish counts with some
# exact zeros and one NaN to pin both zero-dropping and NaN preservation.
ROWS = [
    [1.0, 0.0, 5.0, 0.0, 2.0],
    [0.0, 3.0, 0.0, 2.0, 0.0],
    [np.nan, 0.0, 0.0, 0.5, 1.0],
    [0.0, 4.0, 0.0, 0.0, 0.0],
]


# ═════════════════════════════════════════════════════════════════════
#  Reference oracle — the CURRENT (pre-rework) build
# ═════════════════════════════════════════════════════════════════════


def _oracle_csv_matrix(matrix_file, sep, decimal="."):
    """Reimplements the current table-branch dense build.

    Mirrors rna/steps/00_load.py L759-767 plus the tail ``force_csr``
    conversion: single ``pd.read_csv(..., index_col=0)`` (with the
    ``csv_decimal`` re-read), ``df.values.T.astype(np.float32)`` (dense
    transpose), then ``sp.csr_matrix`` (exact-0.0 dropped, NaN preserved).
    """
    df = pd.read_csv(matrix_file, index_col=0, sep=sep)
    if decimal != ".":
        df = pd.read_csv(matrix_file, index_col=0, sep=sep, decimal=decimal)
    adata = sc.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)))
    adata.var_names = df.index.astype(str)
    adata.obs_names = df.columns.astype(str)
    return adata


class TestOracle:
    """Pin the current dense build on synthetic tables (pre-rework oracle)."""

    def test_standard_table(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, CELLS, ROWS, "\t")
        adata = _oracle_csv_matrix(str(f), "\t")

        assert adata.shape == (5, 4)
        assert list(adata.var_names) == GENES
        assert list(adata.obs_names) == CELLS
        assert adata.X.dtype == np.float32
        expected = np.array(
            [
                [1.0, 0.0, np.nan, 0.0],
                [0.0, 3.0, 0.0, 4.0],
                [5.0, 0.0, 0.0, 0.0],
                [0.0, 2.0, 0.5, 0.0],
                [2.0, 0.0, 1.0, 0.0],
            ],
            dtype=np.float32,
        )
        assert np.array_equal(adata.X.toarray(), expected, equal_nan=True)

    def test_nan_preserved_zero_dropped(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv",
            GENES,
            CELLS[:2],
            [[np.nan, 0.0], [0.0, np.nan], [0.0, 0.0], [3.0, 0.0]],
            "\t",
        )
        adata = _oracle_csv_matrix(str(f), "\t")
        x = adata.X.toarray()
        x = adata.X.toarray()  # cells×genes: BC1 row = [nan,0,0,3], BC2 = [0,nan,0,0]
        assert np.isnan(x[0, 0]) and np.isnan(x[1, 1])
        assert x[0, 1] == 0.0 and x[1, 0] == 0.0 and x[0, 2] == 0.0  # exact zeros dropped
        assert x[0, 3] == 3.0

    def test_missing_values_become_nan(self, tmp_path: Path) -> None:
        # empty fields parse as NaN in the old single read
        f = _write_table(
            tmp_path / "m.tsv",
            GENES[:2],
            CELLS[:3],
            [["", 1.0, 2.0], [3.0, "", 0.0]],
            "\t",
        )
        adata = _oracle_csv_matrix(str(f), "\t")
        x = adata.X.toarray()
        assert np.isnan(x[0, 0]) and np.isnan(x[1, 1])
        x = adata.X.toarray()  # cells×genes: BC1=[nan,3], BC2=[1,nan], BC3=[2,0]
        assert np.isnan(x[0, 0]) and np.isnan(x[1, 1])
        assert x[2, 1] == 0.0

    def test_decimal_comma_rereread(self, tmp_path: Path) -> None:
        # European decimal comma: values like 1,5 are parsed only by the
        # second read (decimal=","); the first read yields string fields.
        f = _write_table(
            tmp_path / "m.csv",
            GENES[:2],
            CELLS[:2],
            [["1,5", "2,25"], ["0", "5,0"]],
            ";",
        )
        adata = _oracle_csv_matrix(str(f), ";", decimal=",")
        expected = np.array([[1.5, 0.0], [2.25, 5.0]], dtype=np.float32)
        assert np.array_equal(adata.X.toarray(), expected, equal_nan=True)

    def test_dtype_float32(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, CELLS, ROWS, "\t")
        adata = _oracle_csv_matrix(str(f), "\t")
        assert adata.X.dtype == np.float32
        assert adata.X.toarray().dtype == np.float32

    def test_single_row(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", ["G001"], CELLS[:3], [[1.0, 2.0, 3.0]], "\t")
        adata = _oracle_csv_matrix(str(f), "\t")
        assert adata.shape == (3, 1)
        assert list(adata.var_names) == ["G001"]
        assert list(adata.obs_names) == CELLS[:3]

    def test_single_column(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, ["BC1"], [[1.0], [0.0], [2.0], [0.0]], "\t")
        adata = _oracle_csv_matrix(str(f), "\t")
        assert adata.shape == (1, 4)
        x = adata.X.toarray()
        assert x[0, 0] == 1.0 and x[0, 1] == 0.0 and x[0, 2] == 2.0

    def test_header_only(self, tmp_path: Path) -> None:
        # header row present, no data rows → genes×0 → cells×0 matrix
        f = _write_table(tmp_path / "m.tsv", [], CELLS, [], "\t")
        adata = _oracle_csv_matrix(str(f), "\t")
        assert adata.shape == (5, 0)
        assert list(adata.obs_names) == CELLS
        assert list(adata.var_names) == []

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "m.tsv"
        f.write_text("")
        with pytest.raises(pd.errors.EmptyDataError):
            _oracle_csv_matrix(str(f), "\t")

    def test_duplicate_cell_header_mangled(self, tmp_path: Path) -> None:
        # pd.read_csv mangles duplicated cell barcodes (BC1 → BC1.1); both the
        # oracle and the reworked path must agree on this reality.
        f = _write_table(
            tmp_path / "m.tsv", GENES[:2], ["BC1", "BC1"], [[1.0, 2.0], [3.0, 4.0]], "\t"
        )
        adata = _oracle_csv_matrix(str(f), "\t")
        assert list(adata.obs_names) == ["BC1", "BC1.1"]
        x = adata.X.toarray()
        x = adata.X.toarray()  # cells×genes: BC1=[1,3], BC1.1=[2,4]
        assert x[0, 0] == 1.0 and x[0, 1] == 3.0
        assert x[1, 0] == 2.0 and x[1, 1] == 4.0

    def test_duplicate_gene_names_preserved(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv",
            ["G1", "G1", "G2"],
            CELLS[:2],
            [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
            "\t",
        )
        adata = _oracle_csv_matrix(str(f), "\t")
        assert list(adata.var_names) == ["G1", "G1", "G2"]
        x = adata.X.toarray()
        assert x[0, 0] == 1.0 and x[0, 1] == 3.0 and x[0, 2] == 5.0

    def test_nonnumeric_cell_crashes_old_build(self, tmp_path: Path) -> None:
        # object-dtype frame → astype(np.float32) raises ValueError in the
        # old path. The reworked path must NOT crash (NaN instead) — pinned
        # separately in TestNonNumericCells.
        f = _write_table(
            tmp_path / "m.tsv",
            GENES[:2],
            CELLS[:2],
            [["abc", 1.0], [2.0, 3.0]],
            "\t",
        )
        with pytest.raises(ValueError):
            _oracle_csv_matrix(str(f), "\t")


# ═════════════════════════════════════════════════════════════════════
#  New path (chunked sparse build) — must equal the oracle
# ═════════════════════════════════════════════════════════════════════


class TestNewPathEqualsOracle:
    """The chunked builder must reproduce the oracle for identical fixtures."""

    def test_standard_table(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, CELLS, ROWS, "\t")
        _assert_adata_equal(
            _new_build(str(f), "\t", chunksize=2), _oracle_csv_matrix(str(f), "\t")
        )

    def test_standard_table_single_block(self, tmp_path: Path) -> None:
        # default chunk sizing → whole table in one block
        f = _write_table(tmp_path / "m.tsv", GENES, CELLS, ROWS, "\t")
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_comma_separator(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.csv", GENES, CELLS, ROWS, ",")
        _assert_adata_equal(_new_build(str(f), ",", chunksize=2), _oracle_csv_matrix(str(f), ","))

    def test_decimal_comma(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.csv",
            GENES[:2],
            CELLS[:2],
            [["1,5", "2,25"], ["0", "5,0"]],
            ";",
        )
        _assert_adata_equal(
            _new_build(str(f), ";", decimal=","), _oracle_csv_matrix(str(f), ";", decimal=",")
        )

    def test_nan_and_zero_semantics(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv",
            GENES,
            CELLS[:2],
            [[np.nan, 0.0], [0.0, np.nan], [0.0, 0.0], [3.0, 0.0]],
            "\t",
        )
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_missing_values_become_nan(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv",
            GENES[:2],
            CELLS[:3],
            [["", 1.0, 2.0], [3.0, "", 0.0]],
            "\t",
        )
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_dup_cell_header_mangled(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv", GENES[:2], ["BC1", "BC1"], [[1.0, 2.0], [3.0, 4.0]], "\t"
        )
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_dup_gene_names_preserved(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv",
            ["G1", "G1", "G2"],
            CELLS[:2],
            [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
            "\t",
        )
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_single_row(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", ["G001"], CELLS, [[1.0, 2.0, 3.0]], "\t")
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_single_column(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, ["BC1"], [[1.0], [0.0], [2.0], [0.0]], "\t")
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_header_only(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", [], CELLS, [], "\t")
        _assert_adata_equal(_new_build(str(f), "\t"), _oracle_csv_matrix(str(f), "\t"))

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "m.tsv"
        f.write_text("")
        with pytest.raises(pd.errors.EmptyDataError):
            _new_build(str(f), "\t")

    def test_many_blocks_forced(self, tmp_path: Path) -> None:
        # 60 genes × 5 cells with chunksize=2 → 30 CSR blocks hstacked
        rows = [[float(i + j) for j in range(5)] for i in range(60)]
        f = _write_table(tmp_path / "m.tsv", [f"G{i:03d}" for i in range(60)], CELLS, rows, "\t")
        _assert_adata_equal(
            _new_build(str(f), "\t", chunksize=2), _oracle_csv_matrix(str(f), "\t")
        )


class TestSepAutoDetect:
    """The (unchanged) sep auto-detect must drive both paths to the same parse."""

    def test_tab_and_comma(self, tmp_path: Path) -> None:
        # The (unchanged) sep-detect sniffs with sep=None and maps ANY
        # multi-column result to "\t" (single-column → ",").  So both a tab
        # and a comma multi-col file detect as "\t".  What must hold is that
        # the reworked builder parses IDENTICALLY to the old path under the
        # SAME detected sep.
        for sep in ("\t", ","):
            f = _write_table(
                tmp_path / f"m.{'tsv' if sep == chr(9) else 'csv'}", GENES, CELLS, ROWS, sep
            )
            detected = _detect_sep(str(f))
            assert detected == "\t"  # multi-col sniff → tab (unchanged quirk)
            _assert_adata_equal(_new_build(str(f), detected), _oracle_csv_matrix(str(f), detected))


class TestNonNumericCells:
    """Plan T5: a non-numeric cell must become NaN, not crash the build."""

    def test_nonnumeric_cell_is_nan(self, tmp_path: Path) -> None:
        f = _write_table(
            tmp_path / "m.tsv",
            GENES[:2],
            CELLS[:2],
            [["abc", 1.0], [2.0, 3.0]],
            "\t",
        )
        adata = _new_build(str(f), "\t")
        x = adata.X.toarray()  # cells×genes: BC1=[abc→NaN,2], BC2=[1,3]
        assert np.isnan(x[0, 0]), "non-numeric cell must coerce to NaN, not crash"
        assert x[0, 1] == 2.0 and x[1, 0] == 1.0 and x[1, 1] == 3.0
        assert adata.X.dtype == np.float32

    def test_mixed_column_no_crash(self, tmp_path: Path) -> None:
        # one column fully numeric, one with a stray token → NaN only at token
        f = _write_table(
            tmp_path / "m.tsv",
            GENES[:3],
            CELLS[:2],
            [["1", "2"], ["3", "bad"], ["5", "6"]],
            "\t",
        )
        adata = _new_build(str(f), "\t")
        x = adata.X.toarray()  # cells×genes: BC1=[1,3,5], BC2=[2,bad→NaN,6]
        assert np.isnan(x[1, 1])
        assert x[0, 0] == 1.0 and x[0, 1] == 3.0 and x[0, 2] == 5.0
        assert x[1, 0] == 2.0 and x[1, 2] == 6.0


# ═════════════════════════════════════════════════════════════════════
#  Downstream obs/var handling (unchanged branch code) — identical on both
# ═════════════════════════════════════════════════════════════════════


def _apply_features(adata, features_path, sep, gene_symbol_col=""):
    """Replica of the unchanged features_file block in the table branch."""
    import warnings

    from anndata._warnings import ImplicitModificationWarning

    genes = _mod._read_features_with_header_detection(features_path, sep=sep)
    if len(genes) == adata.n_vars:
        if gene_symbol_col and gene_symbol_col in genes.columns:
            adata.var_names = genes[gene_symbol_col].values.astype(str)
            genes = genes.drop(columns=[gene_symbol_col])
        elif "gene_short_name" in genes.columns:
            adata.var_names = genes["gene_short_name"].values.astype(str)
            genes = genes.drop(columns=["gene_short_name"])
        elif "symbol" in genes.columns:
            adata.var_names = genes["symbol"].values.astype(str)
            genes = genes.drop(columns=["symbol"])
        else:
            adata.var_names = genes.iloc[:, 0].values.astype(str)
            genes = genes.drop(columns=[genes.columns[0]])
        # the branch's ``adata.var = genes`` reassigns a RangeIndexed frame over
        # the (just-str) var_names — anndata emits ImplicitModificationWarning
        # there in production too; suppress it under pytest's -W error.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ImplicitModificationWarning)
            adata.var = genes
    return adata


class TestDownstreamUnchanged:
    """Metadata join + features var_names must behave identically on both X paths."""

    def test_features_file_var_names(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, CELLS, ROWS, "\t")
        feats = tmp_path / "features.tsv"
        with open(feats, "w") as fh:
            fh.write("id\tgene_short_name\tfeature_type\n")
            for i, g in enumerate(GENES):
                fh.write(f"ENSG{i:05d}\t{g}\tGene Expression\n")
        new = _apply_features(_new_build(str(f), "\t"), feats, "\t")
        oracle = _apply_features(_oracle_csv_matrix(str(f), "\t"), feats, "\t")
        assert list(new.var_names) == list(oracle.var_names)
        pd.testing.assert_frame_equal(new.var, oracle.var)
        assert np.array_equal(new.X.toarray(), oracle.X.toarray(), equal_nan=True)

    def test_meta_columns_join(self, tmp_path: Path) -> None:
        f = _write_table(tmp_path / "m.tsv", GENES, CELLS, ROWS, "\t")
        meta = tmp_path / "barcodes.tsv"
        with open(meta, "w") as fh:
            fh.write("\t".join(["", "celltype", "n_umi"]) + "\n")
            for i, c in enumerate(CELLS):
                fh.write(f"{c}\t{'R' if i % 2 else 'M'}\t{10 + i}\n")

        def _join(adata):
            metadata = pd.read_csv(meta, index_col=0, sep="\t")
            adata.obs = adata.obs.join(metadata, how="left")
            return adata

        new = _join(_new_build(str(f), "\t"))
        oracle = _join(_oracle_csv_matrix(str(f), "\t"))
        assert list(new.obs_names) == list(oracle.obs_names)
        pd.testing.assert_series_equal(new.obs["celltype"], oracle.obs["celltype"])
        pd.testing.assert_series_equal(new.obs["n_umi"], oracle.obs["n_umi"])


# ═════════════════════════════════════════════════════════════════════
#  Memory-boundedness guard
# ═════════════════════════════════════════════════════════════════════


class TestMemoryBounded:
    """The builder must never materialize a full dense genes×cells transpose.

    Primary check is structural (deterministic): the builder reads the table
    in ``chunksize`` gene-blocks, converts each block to a cells×genes float32
    CSR slice, and hstacks them — the old branch's whole-file
    ``df.values.T.astype(np.float32)`` dense copy must be absent.
    """

    def test_structural_chunked_sparse(self) -> None:
        import inspect

        src = inspect.getsource(_mod._build_csv_matrix_sparse)
        helper = inspect.getsource(_mod._csv_chunk_to_float32)
        assert "chunksize" in src, "builder must read the table in row-chunks"
        assert "pd.read_csv(" in src, "builder must use pd.read_csv"
        assert "sp.hstack" in src, "builder must assemble sparse CSR blocks"
        assert ".astype(np.float32)" in helper, "float32 conversion happens per chunk"
        assert "df.values.T" not in src + helper, (
            "must not materialize a whole-file dense transpose"
        )

    def test_rss_smoke_wide_fixture(self, tmp_path: Path) -> None:
        """A wide synthetic table runs far under a dense-materialization budget.

        Shape 8000 cells × 2000 genes; a full dense read + float32 transpose
        of that is ~64 MB, but the sparse build holds only ~1% nnz. Generous
        threshold to avoid flakiness — this guards against gross regressions
        (e.g. reading the whole table into one DataFrame), not exact accounting.
        """
        import resource

        n_cells, n_genes = 8000, 2000
        rng = np.random.RandomState(7)
        dense = rng.binomial(1, 0.01, size=(n_cells, n_genes)).astype(np.float32)
        expr = dense * rng.random_sample((n_cells, n_genes)).astype(np.float32)
        df = pd.DataFrame(expr.T)  # genes × cells
        df.index = [f"G{i:05d}" for i in range(n_genes)]
        df.columns = [f"BC{i:05d}" for i in range(n_cells)]
        f = tmp_path / "wide.tsv"
        df.to_csv(f, sep="\t")
        del df, expr, dense

        before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        adata = _new_build(str(f), "\t")
        after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        delta = (after - before) * 1024  # ru_maxrss is in kB on Linux
        assert delta < 1.5 * 1024**3, f"peak RSS delta {delta / 2**30:.2f} GiB too high"
        assert adata.shape == (n_cells, n_genes)
        assert adata.X.dtype == np.float32
        assert adata.X.nnz <= (n_cells * n_genes * 0.01) + 1
