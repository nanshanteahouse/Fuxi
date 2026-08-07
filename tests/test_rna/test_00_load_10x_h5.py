"""Tests for multi-file 10X HDF5 loading in rna/steps/00_load.py.

Covers
------
- reference oracle: reimplements the CURRENT (pre-rework) multi-file behavior
  (sequential ``sc.concat(join="outer")`` with per-file ``-{i}`` barcode
  suffixes and per-file ``var_names_make_unique``) and pins its output
  (n_obs / n_vars / var order / barcode suffixes / sample / x values)
- fast path: identical gene sets in EXACT order → sparse vstack; output must
  equal the oracle with NO var reordering
- fallback path: differing gene sets → batched outer-join concat; x equal to
  the oracle after reordering columns to the oracle's var order, var SET equal
- permuted gene order (same set, different order) → NOT vstack, falls to batched
- batch-boundary coverage (batch=2 with 5 files)
- fail-fast: corrupt/missing file → SystemExit naming the file

Fixtures are hand-built 10x v3 feature-barcode HDF5 files (scanpy 1.12.2 has
no ``write_10x_h5``): ``matrix/{data,indices,indptr,shape,features/{id,name,
feature_type,genome},barcodes}``.
"""

from __future__ import annotations

import importlib
import os
import sys
import warnings
from pathlib import Path
from unittest.mock import MagicMock

import h5py
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_h5_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Fixture builders (10x v3 h5, by hand — scanpy has no write_10x_h5)
# ═════════════════════════════════════════════════════════════════════


def _write_10x_h5(
    path: Path,
    gene_names: list[str],
    gene_ids: list[str],
    barcodes: list[str],
    x: np.ndarray,
) -> None:
    """Write a minimal 10x v3 feature-barcode matrix.

    ``x`` is a dense (n_cells, n_genes) array.  The file stores indptr over
    cells (AnnData orientation) while ``matrix/shape`` is labeled
    ``[n_genes, n_cells]`` (10x convention) — exactly what scanpy expects.
    """
    x = np.asarray(x, dtype=np.float32)
    csr = sp.csr_matrix(x)
    n_cells, n_genes = x.shape
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int32))
        g.create_dataset("data", data=csr.data.astype(np.int32))
        g.create_dataset("indices", data=csr.indices.astype(np.int64))
        g.create_dataset("indptr", data=csr.indptr.astype(np.int64))
        g.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes], dtype="S20"))
        feats = g.create_group("features")
        feats.create_dataset("name", data=np.array([g.encode() for g in gene_names], dtype="S20"))
        feats.create_dataset("id", data=np.array([i.encode() for i in gene_ids], dtype="S20"))
        feats.create_dataset(
            "feature_type",
            data=np.array([b"Gene Expression"] * len(gene_names), dtype="S20"),
        )
        feats.create_dataset("genome", data=np.array([b"GRCh38"] * len(gene_names), dtype="S20"))


def _rand_counts(seed: int, n_cells: int, n_genes: int) -> np.ndarray:
    return np.random.RandomState(seed).randint(0, 5, size=(n_cells, n_genes)).astype(np.float32)


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════


def _sample_name(fname: str, suffix: str) -> str:
    """Reimplement the module's sample-name extraction from the pattern suffix."""
    sample_name = os.path.basename(fname)
    if suffix and sample_name.endswith(suffix):
        sample_name = sample_name[: -len(suffix)].rstrip("_")
    elif suffix:
        alt = suffix.lstrip("_")
        if alt and sample_name.endswith(alt):
            sample_name = sample_name[: -len(alt)].rstrip("_")
    else:
        sample_name = os.path.splitext(sample_name)[0]
    return sample_name


def _oracle_multi_h5(h5_files: list[str], suffix: str):
    """Reimplement the CURRENT (pre-rework) multi-file 10X_h5 loop exactly.

    Sequential ``sc.concat(join="outer")`` with per-file ``-{i}`` barcode
    suffixes and per-file ``var_names_make_unique`` (dup handling).
    """
    import gc
    import warnings as _warn

    adata = sc.read_10x_h5(h5_files[0])
    if adata.var_names.duplicated().any():
        adata.var_names_make_unique()
    adata.obs_names = [f"{bc}-0" for bc in adata.obs_names]
    adata.obs["sample"] = _sample_name(h5_files[0], suffix)

    for i, h5_file in enumerate(h5_files[1:], start=1):
        with _warn.catch_warnings():
            _warn.simplefilter("ignore", UserWarning)
            a = sc.read_10x_h5(h5_file)
        if a.var_names.duplicated().any():
            a.var_names_make_unique()
        a.obs_names = [f"{bc}-{i}" for bc in a.obs_names]
        a.obs["sample"] = _sample_name(h5_file, suffix)
        adata = sc.concat([adata, a], join="outer")
        del a
        _ = gc.collect()
    return adata


def _make_cfg(h5_dir: str, pattern: str = "*.h5", mtx_concat_batch: int = 0) -> MagicMock:
    cfg = MagicMock()
    cfg.data_input.h5_dir = h5_dir
    cfg.data_input.h5_file_pattern = pattern
    cfg.data_input.mtx_concat_batch = mtx_concat_batch
    return cfg


def _make_capturing_log() -> MagicMock:
    log = MagicMock()
    return log


def _write_sample_files(
    tmp_path: Path,
    files: list[tuple[str, list[str], np.ndarray]],
) -> list[str]:
    """Write a set of h5 fixtures; return sorted glob results for ``*.h5``."""
    d = tmp_path / "h5"
    d.mkdir()
    for name, genes, x in files:
        # distinct barcode prefix per file: "sampleA" -> "A", "s0" -> "0"
        prefix = name[6] if name.startswith("sample") and len(name) > 6 else name[0]
        _write_10x_h5(
            d / name,
            genes,
            [f"ENSG{i:05d}" for i in range(len(genes))],
            [f"AAACCC{prefix}{i:03d}" for i in range(x.shape[0])],
            x,
        )
    return sorted(str(p) for p in d.glob("*.h5"))


GENES_A = [f"GENE_{i:03d}" for i in range(20)]
GENES_B = GENES_A + ["GENE_EXTRA"]  # superset of A


# ═════════════════════════════════════════════════════════════════════
#  Reference oracle — pins the CURRENT (pre-rework) behavior
# ═════════════════════════════════════════════════════════════════════


class TestOracle:
    def test_identical_gene_sets_pinned_values(self, tmp_path: Path) -> None:
        x0 = _rand_counts(1, 5, len(GENES_A))
        x1 = _rand_counts(2, 7, len(GENES_A))
        files = _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_A, x1)]
        )
        adata = _oracle_multi_h5(files, ".h5")

        assert adata.n_obs == 12
        assert adata.n_vars == len(GENES_A)
        # var order = first file's order (identical sets)
        assert list(adata.var_names) == GENES_A
        # barcode suffixes -0 / -1 in file order
        assert list(adata.obs_names) == [f"AAACCCA{i:03d}-0" for i in range(5)] + [
            f"AAACCCB{i:03d}-1" for i in range(7)
        ]
        # sample column from pattern suffix
        assert list(adata.obs["sample"]) == ["sampleA"] * 5 + ["sampleB"] * 7
        # x == vstack of per-file dense matrices (float32)
        expected = np.vstack([x0, x1]).astype(np.float32)
        assert np.array_equal(adata.X.toarray(), expected)

    def test_differing_gene_sets_pinned_values(self, tmp_path: Path) -> None:
        x0 = _rand_counts(3, 5, len(GENES_A))
        x1 = _rand_counts(4, 7, len(GENES_B))
        files = _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_B, x1)]
        )
        adata = _oracle_multi_h5(files, ".h5")

        assert adata.n_obs == 12
        assert adata.n_vars == len(GENES_B)  # union
        assert list(adata.var_names) == GENES_B  # A first, then EXTRA appended
        assert list(adata.obs_names) == [f"AAACCCA{i:03d}-0" for i in range(5)] + [
            f"AAACCCB{i:03d}-1" for i in range(7)
        ]
        # GENE_EXTRA column: sampleA block zero-filled, sampleB block its values
        extra_col = list(adata.var_names).index("GENE_EXTRA")
        mask_a = np.array([s == "sampleA" for s in adata.obs["sample"]])
        assert adata.X[mask_a, extra_col].toarray().sum() == 0
        assert np.allclose(adata.X[~mask_a, extra_col].toarray().ravel(), x1[:, -1])


# ═════════════════════════════════════════════════════════════════════
#  Reworked path — vstack fast path (identical gene sets, EXACT order)
# ═════════════════════════════════════════════════════════════════════


class TestVstackFastPath:
    def test_output_equals_oracle_exact_var_order(self, tmp_path: Path) -> None:
        x0 = _rand_counts(5, 5, len(GENES_A))
        x1 = _rand_counts(6, 7, len(GENES_A))
        files = _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_A, x1)]
        )
        oracle = _oracle_multi_h5(files, ".h5")
        log = _make_capturing_log()
        new = _mod._load_multi_sample_10x_h5(_make_cfg(str(tmp_path / "h5")), log)

        assert new.n_obs == oracle.n_obs
        assert new.n_vars == oracle.n_vars
        # EXACT same var order — no reorder allowed on this path
        assert list(new.var_names) == list(oracle.var_names)
        assert list(new.obs_names) == list(oracle.obs_names)
        assert list(new.obs["sample"]) == list(oracle.obs["sample"])
        assert np.array_equal(new.X.toarray(), oracle.X.toarray(), equal_nan=True)

        # fast path taken
        msgs = [str(a) for a in log.info.call_args_list]
        assert any("vstack" in m for m in msgs)

    def test_three_files_suffixes_pinned(self, tmp_path: Path) -> None:
        x0 = _rand_counts(7, 3, len(GENES_A))
        x1 = _rand_counts(8, 4, len(GENES_A))
        x2 = _rand_counts(9, 5, len(GENES_A))
        files = _write_sample_files(
            tmp_path,
            [
                ("sampleA.h5", GENES_A, x0),
                ("sampleB.h5", GENES_A, x1),
                ("sampleC.h5", GENES_A, x2),
            ],
        )
        oracle = _oracle_multi_h5(files, ".h5")
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        assert list(new.obs_names) == list(oracle.obs_names)
        assert list(new.obs_names[:3]) == [f"AAACCCA{i:03d}-0" for i in range(3)]
        assert list(new.obs_names[3:7]) == [f"AAACCCB{i:03d}-1" for i in range(4)]
        assert list(new.obs_names[7:]) == [f"AAACCCC{i:03d}-2" for i in range(5)]
        assert np.array_equal(new.X.toarray(), oracle.X.toarray(), equal_nan=True)


# ═════════════════════════════════════════════════════════════════════
#  Reworked path — fallback (permuted order / differing gene sets)
# ═════════════════════════════════════════════════════════════════════


class TestFallbackPath:
    def test_permuted_gene_order_not_vstack(self, tmp_path: Path) -> None:
        # same gene SET, but file 2's genes are permuted → ordered comparison
        # must reject the vstack fast path (set equality would silently misalign)
        genes_perm = list(reversed(GENES_A))
        x0 = _rand_counts(10, 5, len(GENES_A))
        x1 = _rand_counts(11, 7, len(GENES_A))
        files = _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", genes_perm, x1)]
        )
        oracle = _oracle_multi_h5(files, ".h5")
        log = _make_capturing_log()
        new = _mod._load_multi_sample_10x_h5(_make_cfg(str(tmp_path / "h5")), log)

        msgs = [str(a) for a in log.info.call_args_list]
        assert not any("vstack" in m for m in msgs)

        # fallback concat: reorder new columns to oracle var order then compare
        reorder = [list(new.var_names).index(g) for g in oracle.var_names]
        assert set(new.var_names) == set(oracle.var_names)
        assert np.array_equal(new.X.toarray()[:, reorder], oracle.X.toarray(), equal_nan=True)

    def test_differing_gene_sets_batched(self, tmp_path: Path) -> None:
        x0 = _rand_counts(12, 5, len(GENES_A))
        x1 = _rand_counts(13, 7, len(GENES_B))
        files = _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_B, x1)]
        )
        oracle = _oracle_multi_h5(files, ".h5")
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        assert new.n_obs == oracle.n_obs
        assert new.n_vars == oracle.n_vars
        assert set(new.var_names) == set(oracle.var_names)
        reorder = [list(new.var_names).index(g) for g in oracle.var_names]
        assert np.array_equal(new.X.toarray()[:, reorder], oracle.X.toarray(), equal_nan=True)

    def test_batch_boundary_5_files_batch2(self, tmp_path: Path) -> None:
        # 5 files, batch=2 → tree merge groups [0,1],[2,3],[4] → [[01],[23],[4]]
        # must equal the sequential oracle exactly
        x = [
            _rand_counts(20 + i, 4, len(GENES_A) if i % 2 == 0 else len(GENES_B)) for i in range(5)
        ]
        files = _write_sample_files(
            tmp_path,
            [(f"s{i}.h5", GENES_A if i % 2 == 0 else GENES_B, x[i]) for i in range(5)],
        )
        oracle = _oracle_multi_h5(files, ".h5")
        log = _make_capturing_log()
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5"), mtx_concat_batch=2), log
        )
        assert new.n_obs == 20
        assert new.n_vars == len(GENES_B)
        reorder = [list(new.var_names).index(g) for g in oracle.var_names]
        assert np.array_equal(new.X.toarray()[:, reorder], oracle.X.toarray(), equal_nan=True)


# ═════════════════════════════════════════════════════════════════════
#  Fail-fast
# ═════════════════════════════════════════════════════════════════════


class TestFailFast:
    def test_corrupt_file_exits_naming_file(self, tmp_path: Path) -> None:
        d = tmp_path / "h5"
        d.mkdir()
        # valid file first
        _write_10x_h5(
            d / "good.h5",
            GENES_A,
            [f"ENSG{i:05d}" for i in range(len(GENES_A))],
            ["AAACCCg001"],
            _rand_counts(30, 1, len(GENES_A)),
        )
        # corrupt file (not a valid h5)
        (d / "broken.h5").write_text("this is not an hdf5 file")

        log = _make_capturing_log()
        with pytest.raises(SystemExit):
            _mod._load_multi_sample_10x_h5(_make_cfg(str(d)), log)
        msg = " ".join(str(a) for a in log.error.call_args.args)
        assert "broken.h5" in msg

    def test_no_matching_files_exits(self, tmp_path: Path) -> None:
        d = tmp_path / "h5"
        d.mkdir()
        log = _make_capturing_log()
        with pytest.raises(SystemExit):
            _mod._load_multi_sample_10x_h5(_make_cfg(str(d), pattern="*.h5"), log)
        msg = " ".join(str(a) for a in log.error.call_args.args)
        assert "*.h5" in msg


# ═════════════════════════════════════════════════════════════════════
#  In-place CSR assembly (T1b) — identical-gene-set fast path
# ═════════════════════════════════════════════════════════════════════


def _oracle_vstack(xs):
    """Reference CSR for the fast path: per-file dense arrays, vstacked."""
    return sp.vstack([sp.csr_matrix(x.astype(np.float32)) for x in xs], format="csr")


class TestInPlaceAssembly:
    def test_multi_file_csr_equals_oracle_arrays(self, tmp_path: Path) -> None:
        x0 = _rand_counts(60, 5, len(GENES_A))
        x1 = _rand_counts(61, 7, len(GENES_A))
        x2 = _rand_counts(62, 3, len(GENES_A))
        _write_sample_files(
            tmp_path,
            [
                ("sampleA.h5", GENES_A, x0),
                ("sampleB.h5", GENES_A, x1),
                ("sampleC.h5", GENES_A, x2),
            ],
        )
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        ref = _oracle_vstack([x0, x1, x2])
        # direct CSR equality — per-row sorted, no re-sort needed
        assert np.array_equal(new.X.data, ref.data, equal_nan=True)
        assert np.array_equal(new.X.indices, ref.indices)
        assert np.array_equal(new.X.indptr, ref.indptr)
        assert np.array_equal(new.X.toarray(), ref.toarray(), equal_nan=True)

    def test_indptr_offsets_across_file_boundaries(self, tmp_path: Path) -> None:
        # crafted: file0 has an empty row; boundary must jump by file0 total nnz
        x0 = np.array(
            [
                [1, 0, 2, 0],  # row0 nnz 2
                [0, 0, 0, 0],  # row1 nnz 0 (empty)
                [3, 0, 0, 4],  # row2 nnz 2
            ],
            dtype=np.float32,
        )
        x1 = np.array([[5, 0, 0, 0], [0, 6, 0, 0]], dtype=np.float32)
        _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A[:4], x0), ("sampleB.h5", GENES_A[:4], x1)]
        )
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        x = new.X
        # merged 5 rows; file0 rows 0-2 then file1 rows 3-4
        assert list(x.indptr) == [0, 2, 2, 4, 5, 6]
        assert x.indptr[3] == 4  # boundary start == file0 total nnz
        ref = _oracle_vstack([x0, x1])
        assert np.array_equal(x.indptr, ref.indptr)

    def test_single_file_edge(self, tmp_path: Path) -> None:
        x0 = _rand_counts(70, 8, len(GENES_A))
        _write_sample_files(tmp_path, [("sampleA.h5", GENES_A, x0)])
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        ref = _oracle_vstack([x0])
        assert new.n_obs == x0.shape[0]
        assert new.n_vars == x0.shape[1]
        assert np.array_equal(new.X.data, ref.data, equal_nan=True)
        assert np.array_equal(new.X.indices, ref.indices)
        assert np.array_equal(new.X.indptr, ref.indptr)

    def test_empty_rows_edge(self, tmp_path: Path) -> None:
        # all rows of file1 are zero → indptr must still advance correctly
        x0 = _rand_counts(80, 3, len(GENES_A))
        x1 = np.zeros((4, len(GENES_A)), dtype=np.float32)
        _write_sample_files(tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_A, x1)])
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        ref = _oracle_vstack([x0, x1])
        assert np.array_equal(new.X.data, ref.data, equal_nan=True)
        assert np.array_equal(new.X.indices, ref.indices)
        assert np.array_equal(new.X.indptr, ref.indptr)
        assert list(new.X.indptr[3:]) == [ref.nnz] * 5  # all-zero block

    def test_dtypes_float32_int32(self, tmp_path: Path) -> None:
        x0 = _rand_counts(90, 5, len(GENES_A))
        x1 = _rand_counts(91, 7, len(GENES_A))
        _write_sample_files(tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_A, x1)])
        new = _mod._load_multi_sample_10x_h5(
            _make_cfg(str(tmp_path / "h5")), _make_capturing_log()
        )
        assert new.X.data.dtype == np.float32
        assert new.X.indices.dtype == np.int32
        assert new.X.indptr.dtype == np.int32

    def test_read_10x_h5_meta_parity_with_scanpy(self, tmp_path: Path) -> None:
        # duplicate gene to exercise var_names_make_unique replication
        genes = GENES_A[:5] + [GENES_A[2]]
        x = _rand_counts(92, 6, len(genes))
        d = tmp_path / "h5"
        d.mkdir()
        _write_10x_h5(
            d / "s0.h5",
            genes,
            [f"ENSG{i:05d}" for i in range(len(genes))],
            [f"BC{i:03d}" for i in range(6)],
            x,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            a = sc.read_10x_h5(str(d / "s0.h5"))
        if a.var_names.duplicated().any():
            a.var_names_make_unique()
        n_obs, nnz, var_names, obs_index, var_df = _mod._read_10x_h5_meta(str(d / "s0.h5"))
        assert n_obs == a.n_obs
        assert nnz == a.X.nnz
        assert list(var_names) == list(a.var_names)
        assert list(obs_index) == list(a.obs.index)
        assert list(var_df.index) == list(a.var_names)
        assert list(var_df.columns) == list(a.var.columns)
        for c in a.var.columns:
            assert list(var_df[c]) == list(a.var[c])

    def test_assemble_csr_inplace_equals_vstack(self, tmp_path: Path) -> None:
        x0 = _rand_counts(93, 4, len(GENES_A))
        x1 = _rand_counts(94, 6, len(GENES_A))
        files = _write_sample_files(
            tmp_path, [("sampleA.h5", GENES_A, x0), ("sampleB.h5", GENES_A, x1)]
        )
        meta = []
        for fp in files:
            n_obs, nnz, *_ = _mod._read_10x_h5_meta(fp)
            meta.append((fp, n_obs, nnz))
        out = _mod._assemble_csr_inplace(meta, len(GENES_A))
        ref = _oracle_vstack([x0, x1])
        assert np.array_equal(out.data, ref.data, equal_nan=True)
        assert np.array_equal(out.indices, ref.indices)
        assert np.array_equal(out.indptr, ref.indptr)
        assert out.data.dtype == np.float32
        assert out.indices.dtype == np.int32
        assert out.indptr.dtype == np.int32

    def test_indptr_add_handles_int32_overflow(self) -> None:
        # regression (T1b QA): per-file indptr is int32, but the global offset
        # exceeds int32 max at 2.8e9 nnz — the naive ``x.indptr[:-1] + off``
        # form raises OverflowError on numpy 2; the int64 cast must not.
        x = sp.csr_matrix(np.array([[1, 0], [0, 2]], dtype=np.float32))
        assert x.indptr.dtype == np.int32
        big_off = 2_150_585_452  # > int32 max
        vals = x.indptr[:-1].astype(np.int64) + big_off
        assert list(vals) == [big_off, big_off + 1]
        with pytest.raises(OverflowError):
            _ = x.indptr[:-1] + big_off  # the pre-fix form, as a guard
