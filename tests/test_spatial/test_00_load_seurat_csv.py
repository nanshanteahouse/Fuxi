"""TDD red-phase tests for the ``seurat_csv`` load path in spatial/steps/00_load.py.

Branch under test: the future ``elif data_format == "seurat_csv"`` in
``spatial/steps/00_load.py`` — Seurat wide-format ``counts.csv.gz`` (rows =
ENSG gene IDs, columns = spot barcodes ``cell_N_M``) plus an ``md.csv.gz``
metadata table, assembled into an AnnData with X = spots×genes sparse matrix
(float32, int64 indptr), ``obsm['spatial']`` from the md ``pixel_x``/``pixel_y``
columns, obs carrying the md metadata columns, and ``var_names`` left as the raw
ENSG IDs (no symbol conversion).

TDD state: RED. The branch is not implemented yet (plan
``.omo/plans/spatial-pipeline-rewrite-phase2.md`` todo 3, implemented in todo 4).
``TestOracle`` pins the expected behaviour on synthetic GSE235583-shaped files
and is fully self-contained — it passes now.  ``TestNewPathEqualsOracle``
asserts that the not-yet-existing module helper ``_load_seurat_csv`` reproduces
the oracle, so it FAILS at an explicit ``hasattr`` guard (a clear AssertionError,
not an AttributeError) until todo 4 lands.

Synthetic files mirror the real GSE235583 layout (verified against the local
copy at ``$FUXI_DATA_ROOT`` — never used in the tests themselves):

- ``counts.csv.gz`` header: empty index-name field, then ``cell_1_1,cell_2_1,...``
- ``md.csv.gz`` header: empty index-name field, then
  ``orig.ident,nCount_originalexp,nFeature_originalexp,Sample,Barcode,Section,
  Spot_Y,Spot_X,Image_Y,Image_X,pixel_x,pixel_y,...`` — so file field 6 = Barcode
  and file field 7 = Section (1-based, including the leading rowname column).

Contract for the future implementation (todo 4): extract a module-level helper
mirroring ``rna/steps/00_load.py::_build_csv_matrix_sparse``::

    def _load_seurat_csv(counts_file: str, md_file: str, log: logging.Logger) -> sc.AnnData:
        ...

The optional reassembly of the ``GSE235583_RAW/`` Space-Ranger aux files into
``uns['spatial']`` is todo-4 scope and is not pinned here.
"""

from __future__ import annotations

import gzip
import importlib
import importlib.util
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "spatial", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location(
    "spatial.steps._00_load_seurat_csv_test", _STEP_PATH
)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Fixture builders / helpers
# ═════════════════════════════════════════════════════════════════════

GENES = [
    "ENSG00000243485",
    "ENSG00000237613",
    "ENSG00000186092",
    "ENSG00000238009",
    "ENSG00000239945",
    "ENSG00000241860",
    "ENSG00000241599",
    "ENSG00000235146",
]
N_SPOTS = 12
SPOTS = [f"cell_{i + 1}_1" for i in range(N_SPOTS)]


def _fmt_val(v: float) -> str:
    """Serialize one count value the way Seurat's wide CSV does (ints, ``nan``)."""
    v = float(v)
    if np.isnan(v):
        return "nan"
    return str(int(v)) if v.is_integer() else repr(v)


def _make_md(spots: list[str], seed: int) -> pd.DataFrame:
    """Build an md.csv table shaped like the real GSE235583 ``md.csv.gz``.

    Column order reproduces the real header for the first six named columns so
    that file field 6 = Barcode and file field 7 = Section (1-based, counting
    the leading rowname column), followed by the pixel coordinates and QC extras.
    """
    rng = np.random.RandomState(seed)
    n = len(spots)
    md = pd.DataFrame(index=spots)
    md["orig.ident"] = "AD3_D10"
    md["nCount_originalexp"] = rng.randint(300, 4000, size=n)
    md["nFeature_originalexp"] = rng.randint(200, 1800, size=n)
    md["Sample"] = "GSM7505850_AD3_D10_A1"
    md["Barcode"] = [f"AAACCCAAGAATTGTC-{i:02d}" for i in range(n)]
    md["Section"] = 1
    md["Spot_Y"] = rng.randint(0, 128, size=n)
    md["Spot_X"] = rng.randint(0, 128, size=n)
    md["Image_Y"] = rng.randint(0, 2000, size=n)
    md["Image_X"] = rng.randint(0, 2000, size=n)
    md["pixel_x"] = rng.uniform(50.0, 600.0, size=n)
    md["pixel_y"] = rng.uniform(50.0, 600.0, size=n)
    md["id"] = "AD3_D10_Retinal_Organoids_A1"
    md["day"] = 10
    md["tissue"] = "AD3"
    md["section"] = "A1"
    md["seurat_clusters"] = rng.randint(0, 3, size=n)
    return md


def _make_seurat_fixture(
    tmp_path: Path,
    genes: list[str] | None = None,
    spots: list[str] | None = None,
    seed: int = 0,
    with_nan: bool = False,
) -> tuple[Path, Path, np.ndarray, pd.DataFrame]:
    """Write a synthetic GSE235583-shaped ``counts.csv.gz`` + ``md.csv.gz``.

    Returns ``(counts_path, md_path, matrix, md)`` where ``matrix`` is the
    genes×spots float32 count table that was serialized and ``md`` the metadata
    DataFrame (index = spot IDs) that was serialized — so tests can assert
    against the exact written content.  Deterministic via ``RandomState(seed)``.
    """
    genes = genes if genes is not None else GENES
    spots = spots if spots is not None else SPOTS
    rng = np.random.RandomState(seed)
    matrix = rng.poisson(lam=0.2, size=(len(genes), len(spots))).astype(np.float32)
    if with_nan:
        matrix[0, 0] = np.nan
    md = _make_md(spots, seed + 1)
    counts_path = tmp_path / "counts.csv.gz"
    md_path = tmp_path / "md.csv.gz"
    with gzip.open(counts_path, "wt") as f:
        f.write("," + ",".join(spots) + "\n")
        for g, row in zip(genes, matrix):
            f.write(g + "," + ",".join(_fmt_val(v) for v in row) + "\n")
    md.to_csv(md_path, compression="gzip")
    return counts_path, md_path, matrix, md


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_load_seurat_csv")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


# ═════════════════════════════════════════════════════════════════════
#  Reference oracle — the expected seurat_csv build
# ═════════════════════════════════════════════════════════════════════


def _oracle_seurat_csv(counts_file: str, md_file: str) -> sc.AnnData:
    """Reference reimplementation of the expected ``seurat_csv`` load.

    Mirrors the target branch: one ``pd.read_csv(index_col=0)`` for the wide
    counts table (genes×spots), dense float32 transpose, ``sp.csr_matrix`` with
    the indptr widened to int64; then the md table for obs and
    ``obsm['spatial']`` from the ``pixel_x``/``pixel_y`` columns.
    """
    counts = pd.read_csv(counts_file, index_col=0)
    x_csr = sp.csr_matrix(counts.values.T.astype(np.float32))
    x_csr.indptr = x_csr.indptr.astype(np.int64)
    adata = sc.AnnData(X=x_csr)
    adata.var_names = counts.index.astype(str)
    adata.obs_names = counts.columns.astype(str)
    md = pd.read_csv(md_file, index_col=0)
    adata.obs = md
    adata.obsm["spatial"] = md[["pixel_x", "pixel_y"]].to_numpy()
    return adata


def _new_load(counts_file: str, md_file: str) -> sc.AnnData:
    """Invoke the future ``_load_seurat_csv`` builder in spatial/steps/00_load.py.

    The helper does not exist yet — that is the TDD red phase.  A clear
    AssertionError (rather than an AttributeError) points at the missing piece
    until plan todo 4 implements the branch.
    """
    assert hasattr(_mod, "_load_seurat_csv"), (
        "spatial 00_load seurat_csv loader not implemented yet — expected TDD red "
        "phase (plan todo 3); add _load_seurat_csv(counts_file, md_file, log) to "
        "spatial/steps/00_load.py (plan todo 4)"
    )
    return _mod._load_seurat_csv(counts_file, md_file, _make_logger())


def _assert_adata_equal(new: sc.AnnData, oracle: sc.AnnData) -> None:
    """New-path output must match the reference oracle exactly."""
    assert new.shape == oracle.shape
    assert list(new.var_names) == list(oracle.var_names), "var order must match oracle"
    assert list(new.obs_names) == list(oracle.obs_names), "obs_names order must match oracle"
    assert new.X.dtype == np.float32, "X dtype must be float32"
    assert sp.isspmatrix_csr(new.X), "X must be CSR"
    assert new.X.indptr.dtype == np.int64, "indptr must be int64"
    assert np.array_equal(new.X.toarray(), oracle.X.toarray(), equal_nan=True)
    pd.testing.assert_frame_equal(new.obs, oracle.obs)
    assert "spatial" in new.obsm
    np.testing.assert_allclose(new.obsm["spatial"], oracle.obsm["spatial"])


# ═════════════════════════════════════════════════════════════════════
#  Oracle — pins the expected seurat_csv behaviour (passes now)
# ═════════════════════════════════════════════════════════════════════


class TestOracle:
    """Pin the expected seurat_csv load behaviour on synthetic GSE235583 files."""

    def test_transposed_sparse_matrix(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))

        assert adata.shape == (N_SPOTS, len(GENES))
        assert list(adata.var_names) == GENES
        assert list(adata.obs_names) == SPOTS
        assert sp.isspmatrix_csr(adata.X)
        assert adata.X.dtype == np.float32
        assert adata.X.indptr.dtype == np.int64
        # row = spot, column = gene must equal the serialized matrix[gene, spot]
        assert np.array_equal(adata.X.toarray(), matrix.T, equal_nan=True)

    def test_obsm_spatial_from_pixel_coords(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))

        assert "spatial" in adata.obsm
        assert adata.obsm["spatial"].dtype == np.float64
        np.testing.assert_allclose(adata.obsm["spatial"], md[["pixel_x", "pixel_y"]].to_numpy())

    def test_obs_metadata_from_md(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))

        assert list(adata.obs.columns) == list(md.columns)
        pd.testing.assert_frame_equal(adata.obs, md)
        assert adata.obs["Barcode"].tolist() == md["Barcode"].tolist()
        assert adata.obs["Section"].tolist() == md["Section"].tolist()

    def test_md_column_positions(self, tmp_path: Path) -> None:
        # Real GSE235583 md.csv: file field 6 = Barcode, field 7 = Section
        # (1-based, counting the leading rowname column).
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        fields = list(pd.read_csv(md_path, nrows=0).columns)
        assert fields[5] == "Barcode"
        assert fields[6] == "Section"

    def test_gene_names_preserved_ensg(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))

        # var_names stay the raw ENSG IDs — unchanged, no symbol mapping
        assert all(str(g).startswith("ENSG") for g in adata.var_names)
        assert list(adata.var_names) == GENES

    def test_exact_zeros_dropped(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))

        assert adata.X.nnz == np.count_nonzero(matrix.T)
        assert adata.X.nnz < matrix.size  # Poisson(0.2) → mostly zeros

    def test_nan_preserved(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path, with_nan=True)
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))

        assert np.isnan(adata.X.toarray()).sum() == 1
        assert np.array_equal(adata.X.toarray(), matrix.T, equal_nan=True)

    def test_single_gene(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(
            tmp_path, genes=["ENSG00000243485"], spots=SPOTS[:3]
        )
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))
        assert adata.shape == (3, 1)
        assert list(adata.var_names) == ["ENSG00000243485"]
        assert list(adata.obs_names) == SPOTS[:3]

    def test_single_spot(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(
            tmp_path, genes=GENES[:3], spots=["cell_1_1"]
        )
        adata = _oracle_seurat_csv(str(counts_path), str(md_path))
        assert adata.shape == (1, 3)
        assert list(adata.obs_names) == ["cell_1_1"]
        assert list(adata.var_names) == GENES[:3]


# ═════════════════════════════════════════════════════════════════════
#  New path — must equal the oracle (RED until plan todo 4 lands)
# ═════════════════════════════════════════════════════════════════════


class TestNewPathEqualsOracle:
    """The future seurat_csv branch must reproduce the oracle on identical fixtures.

    RED (TDD): ``_mod._load_seurat_csv`` is not implemented yet, so every test
    here fails at the ``hasattr`` guard until plan todo 4 lands.
    """

    def test_standard_fixture(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path)
        _assert_adata_equal(
            _new_load(str(counts_path), str(md_path)),
            _oracle_seurat_csv(str(counts_path), str(md_path)),
        )

    def test_small_sparse(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(
            tmp_path, genes=GENES[:4], spots=SPOTS[:5]
        )
        _assert_adata_equal(
            _new_load(str(counts_path), str(md_path)),
            _oracle_seurat_csv(str(counts_path), str(md_path)),
        )

    def test_nan_semantics(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path, with_nan=True)
        _assert_adata_equal(
            _new_load(str(counts_path), str(md_path)),
            _oracle_seurat_csv(str(counts_path), str(md_path)),
        )

    def test_single_gene(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(
            tmp_path, genes=GENES[:1], spots=SPOTS
        )
        _assert_adata_equal(
            _new_load(str(counts_path), str(md_path)),
            _oracle_seurat_csv(str(counts_path), str(md_path)),
        )

    def test_single_spot(self, tmp_path: Path) -> None:
        counts_path, md_path, matrix, md = _make_seurat_fixture(
            tmp_path, genes=GENES, spots=SPOTS[:1]
        )
        _assert_adata_equal(
            _new_load(str(counts_path), str(md_path)),
            _oracle_seurat_csv(str(counts_path), str(md_path)),
        )

    def test_many_genes(self, tmp_path: Path) -> None:
        genes = [f"ENSG{i:011d}" for i in range(60)]
        counts_path, md_path, matrix, md = _make_seurat_fixture(tmp_path, genes=genes)
        _assert_adata_equal(
            _new_load(str(counts_path), str(md_path)),
            _oracle_seurat_csv(str(counts_path), str(md_path)),
        )
