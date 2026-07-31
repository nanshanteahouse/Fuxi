"""Tests for incremental h5ad I/O — h5ad-incremental-io Phase 1.

Covers: int8 codes overflow (1.0a), NaN→"nan" string corruption (1.0b),
write_elem migration of the lightweight obs writer (1.1), in-place writeback
backup path (1.2-B), and the write_h5ad_incremental engine + safe_write
delta_only routing (1.5).
"""

from __future__ import annotations

import os

import h5py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

from core.utils import write_obs_columns_lightweight

# ── fixtures ────────────────────────────────────────────────────────────


def _make_h5ad(tmp_path, n_obs=300, n_var=20, compression="gzip") -> tuple:
    """Create a base h5ad file; return (base_adata, path)."""
    rng = np.random.default_rng(0)
    adata = sc.AnnData(X=rng.random((n_obs, n_var)))
    adata.obs["batch"] = pd.Categorical(["a", "b", "c"] * (n_obs // 3))
    adata.obsm["X_pca"] = rng.random((n_obs, 10))
    path = str(tmp_path / "base.h5ad")
    adata.write(path, compression=compression)
    return adata, path


@pytest.fixture
def base_h5ad(tmp_path):
    return _make_h5ad(tmp_path)


@pytest.fixture
def base_h5ad_uncompressed(tmp_path):
    return _make_h5ad(tmp_path, compression=None)


# ── 1.0a / 1.1: categorical round-trip incl. >127 categories ────────────


def test_obs_categorical_over_127_roundtrip(base_h5ad):
    """>127 categories must round-trip without int8 wraparound (1.0a)."""
    _, path = base_h5ad
    n = 300
    cats = pd.Categorical([f"cat_{i % 200}" for i in range(n)])
    write_obs_columns_lightweight(path, pd.DataFrame({"big": cats}))

    adata = sc.read(path)
    assert (adata.obs["big"] == cats).all()
    assert adata.obs["big"].nunique() == 200

    # codes must not be int8 (would wrap at 127 categories)
    with h5py.File(path, "r") as f:
        codes_dtype = f["obs/big/codes"].dtype
    assert codes_dtype in (np.int16, np.int32)


def test_obs_categorical_small_roundtrip(base_h5ad):
    """Small category counts must keep values + categories intact."""
    _, path = base_h5ad
    n = 300
    cats = pd.Categorical(["alpha", "beta", "gamma"] * (n // 3))
    write_obs_columns_lightweight(path, pd.DataFrame({"ct": cats}))

    adata = sc.read(path)
    assert (adata.obs["ct"] == cats).all()
    assert list(adata.obs["ct"].cat.categories) == ["alpha", "beta", "gamma"]


# ── 1.0b / 1.1: string columns with NaN ────────────────────────────────


def test_obs_string_nan_roundtrip(base_h5ad):
    """NaN in a string column must stay missing, not become literal 'nan' (1.0b)."""
    _, path = base_h5ad
    n = 300
    s = pd.Series(["alpha", None, "beta"] * (n // 3))
    write_obs_columns_lightweight(path, pd.DataFrame({"lab": s}))

    adata = sc.read(path)
    assert adata.obs["lab"].isna().sum() == n // 3
    assert "nan" not in set(adata.obs["lab"].dropna().astype(str))


def test_obs_string_no_nan_roundtrip(base_h5ad):
    """Plain string column round-trips exactly."""
    _, path = base_h5ad
    n = 300
    s = pd.Series([f"sample_{i % 12}" for i in range(n)])
    write_obs_columns_lightweight(path, pd.DataFrame({"lab": s}))

    adata = sc.read(path)
    assert adata.obs["lab"].to_numpy().tolist() == s.tolist()


# ── 1.1: column-order & dtype preservation ─────────────────────────────


def test_obs_column_order_preserved(base_h5ad):
    """Existing columns keep order; new columns append at the end."""
    _, path = base_h5ad
    n = 300
    df = pd.DataFrame(
        {
            "aa": np.arange(n, dtype=np.float64),
            "bb": np.arange(n, dtype=np.int64),
        }
    )
    write_obs_columns_lightweight(path, df)

    adata = sc.read(path)
    assert list(adata.obs.columns) == ["batch", "aa", "bb"]


def test_obs_dtype_preserved(base_h5ad):
    """float64/int64/bool round-trip with their original dtypes."""
    _, path = base_h5ad
    n = 300
    df = pd.DataFrame(
        {
            "f64": np.arange(n, dtype=np.float64) + 0.5,
            "i64": np.arange(n, dtype=np.int64),
            "f32": np.arange(n, dtype=np.float32),
            "i32": np.arange(n, dtype=np.int32),
            "flag": np.array([True, False] * (n // 2), dtype=bool),
        }
    )
    write_obs_columns_lightweight(path, df)

    adata = sc.read(path)
    assert adata.obs["f64"].dtype == np.float64
    assert adata.obs["i64"].dtype == np.int64
    assert adata.obs["f32"].dtype == np.float32
    assert adata.obs["i32"].dtype == np.int32
    assert adata.obs["flag"].dtype == np.bool_


def test_obs_overwrite_existing_column(base_h5ad):
    """Re-writing an existing obs column replaces its values (overwrite semantics)."""
    _, path = base_h5ad
    n = 300
    write_obs_columns_lightweight(path, pd.DataFrame({"batch": np.arange(n)}))
    adata = sc.read(path)
    assert adata.obs["batch"].dtype == np.int64
    assert adata.obs["batch"].iloc[0] == 0


# ── 1.5: write_h5ad_incremental engine ───────────────────────────────


def test_engine_obs_obsm_obsp_uns_roundtrip(base_h5ad):
    """Engine writes obs/obsm/obsp/uns in one call; all round-trip."""
    from core.utils import write_h5ad_incremental

    _, path = base_h5ad
    n = 300
    rng = np.random.default_rng(1)
    obs_df = pd.DataFrame(
        {
            "leiden": pd.Categorical([f"cl{i % 8}" for i in range(n)]),
            "score": rng.random(n),
        }
    )
    obsm = {"X_umap": rng.random((n, 2))}
    obsp = {"connectivities": sp.csr_matrix(rng.random((n, n)) < 0.01, dtype=np.float32)}
    uns = {
        "neighbors": {
            "params": {"n_neighbors": 15, "method": "umap"},
            "distances_key": "distances",
        },
        "leiden": {"params": {"resolution": 0.5}},
        "scalar": 42,
    }

    write_h5ad_incremental(path, obs=obs_df, obsm=obsm, obsp=obsp, uns=uns)

    adata = sc.read(path)
    assert adata.obs["leiden"].to_numpy().tolist() == obs_df["leiden"].to_numpy().tolist()
    assert np.allclose(adata.obs["score"], obs_df["score"])
    assert np.allclose(adata.obsm["X_umap"], obsm["X_umap"])
    got = adata.obsp["connectivities"].tocsr()
    assert (got != obsp["connectivities"]).nnz == 0
    assert adata.uns["neighbors"] == uns["neighbors"]
    assert adata.uns["leiden"] == uns["leiden"]
    assert adata.uns["scalar"] == 42


def test_engine_uns_dataframe_roundtrip(base_h5ad):
    """pandas DataFrames in uns (e.g. rank_genes_groups) round-trip."""
    from core.utils import write_h5ad_incremental

    _, path = base_h5ad
    df = pd.DataFrame(
        {
            "names": ["GENE_A", "GENE_B"],
            "scores": [1.5, 0.5],
        }
    )
    write_h5ad_incremental(path, uns={"rank_genes_groups": df})
    adata = sc.read(path)
    assert adata.uns["rank_genes_groups"].equals(df)


def test_engine_overwrite_semantics(base_h5ad):
    """Re-running replaces stale keys (old X_umap / leiden_* must not survive)."""
    from core.utils import write_h5ad_incremental

    _, path = base_h5ad
    n = 300
    first = np.full((n, 2), 1.0)
    second = np.full((n, 2), 7.0)
    write_h5ad_incremental(path, obsm={"X_umap": first})
    write_h5ad_incremental(path, obsm={"X_umap": second})
    adata = sc.read(path)
    assert np.all(adata.obsm["X_umap"] == 7.0)

    df1 = pd.DataFrame({"leiden": pd.Categorical(["a", "b", "c"] * (n // 3))})
    df2 = pd.DataFrame({"leiden": np.arange(n)})
    write_h5ad_incremental(path, obs=df1)
    write_h5ad_incremental(path, obs=df2)
    adata = sc.read(path)
    assert adata.obs["leiden"].dtype == np.int64


def test_engine_X_untouched(base_h5ad):
    """X content and compression must be identical after incremental writes."""
    from core.utils import write_h5ad_incremental

    base, path = base_h5ad
    n = 300
    rng = np.random.default_rng(2)
    write_h5ad_incremental(
        path,
        obs=pd.DataFrame({"k": np.arange(n)}),
        obsm={"X_umap": rng.random((n, 2))},
        uns={"note": "x"},
    )
    adata = sc.read(path)
    assert np.array_equal(adata.X, base.X)
    with h5py.File(path, "r") as f:
        assert f["X"].compression == "gzip"
        # base file has an empty /layers group and no /raw — untouched by appends
        assert len(f["layers"].keys()) == 0
        assert "raw" not in f


def test_engine_rejects_missing_file(tmp_path):
    """Appending to a non-existent file must fail loudly."""
    from core.utils import write_h5ad_incremental

    missing = str(tmp_path / "nope.h5ad")
    with pytest.raises(FileNotFoundError):
        write_h5ad_incremental(missing, obs=pd.DataFrame({"a": [1]}))


def test_engine_obs_length_mismatch(base_h5ad):
    """Obs length must match the file's n_obs."""
    from core.utils import write_h5ad_incremental

    _, path = base_h5ad
    with pytest.raises(ValueError, match="n_obs"):
        write_h5ad_incremental(path, obs=pd.DataFrame({"a": [1, 2]}))


# ── 1.5: compression consistency ─────────────────────────────────────


def test_engine_gzip_file_uses_gzip_append(base_h5ad):
    """gzip h5ad → appended datasets are also gzip."""
    from core.utils import write_h5ad_incremental

    _, path = base_h5ad
    n = 300
    write_h5ad_incremental(path, obsm={"X_umap": np.random.default_rng(3).random((n, 2))})
    with h5py.File(path, "r") as f:
        assert f["obsm/X_umap"].compression == "gzip"


def test_engine_uncompressed_file_appends_uncompressed(base_h5ad_uncompressed):
    """Uncompressed (ATAC-style) h5ad → appended datasets stay uncompressed."""
    from core.utils import write_h5ad_incremental

    _, path = base_h5ad_uncompressed
    n = 300
    write_h5ad_incremental(path, obsm={"X_umap": np.random.default_rng(4).random((n, 2))})
    with h5py.File(path, "r") as f:
        assert f["X"].compression is None
        assert f["obsm/X_umap"].compression is None


# ── 1.2-B: in-place writeback with .bak backup/restore ───────────────


def test_inplace_writeback_success(base_h5ad):
    """Success unlinks the .bak; column is present."""
    from core.utils import write_obs_columns_inplace

    _, path = base_h5ad
    n = 300
    write_obs_columns_inplace(path, pd.DataFrame({"cell_subtype": np.arange(n)}))
    assert not os.path.exists(path + ".bak")
    adata = sc.read(path)
    assert "cell_subtype" in adata.obs.columns


def test_inplace_writeback_restores_on_failure(base_h5ad, monkeypatch):
    """Failed append restores the original file from .bak atomically."""
    import core.utils._io as io_mod
    from core.utils import write_obs_columns_inplace

    _, path = base_h5ad
    with open(path, "rb") as fh:
        original = fh.read()
    n = 300

    def _explode(h5ad_path, obs=None, obsm=None, obsp=None, uns=None, logger=None):
        with open(h5ad_path, "ab") as fh:
            fh.write(b"GARBAGE-CORRUPTED")
        raise RuntimeError("simulated append crash")

    monkeypatch.setattr(io_mod, "write_h5ad_incremental", _explode)
    with pytest.raises(RuntimeError, match="simulated append crash"):
        write_obs_columns_inplace(path, pd.DataFrame({"cell_subtype": np.arange(n)}))

    assert not os.path.exists(path + ".bak")
    with open(path, "rb") as fh:
        assert fh.read() == original


# ── 1.5: safe_write(delta_only=...) routing ───────────────────────────


def test_safe_write_delta_only_updates_in_place(base_h5ad):
    """delta_only=True appends to the existing file without touching X."""
    from core.utils import safe_write

    base, path = base_h5ad
    n = 300
    rng = np.random.default_rng(5)
    new = sc.AnnData(X=base.X.copy())
    new.obs["leiden"] = pd.Categorical([f"c{i % 5}" for i in range(n)])
    new.obsm["X_umap"] = rng.random((n, 2))

    safe_write(new, path, delta_only=True)
    adata = sc.read(path)
    assert "leiden" in adata.obs.columns
    assert np.allclose(adata.obsm["X_umap"], new.obsm["X_umap"])
    assert np.array_equal(adata.X, base.X)


def test_safe_write_delta_only_creates_when_missing(tmp_path):
    """delta_only=True with no target file falls back to the full write."""
    from core.utils import safe_write

    target = str(tmp_path / "new.h5ad")
    adata = sc.AnnData(X=np.random.default_rng(6).random((50, 5)))
    safe_write(adata, target, delta_only=True)
    assert os.path.exists(target)
    assert sc.read(target).n_obs == 50


def test_safe_write_default_full_path_unaffected(base_h5ad):
    """Existing full-write path (no delta_only) keeps working."""
    from core.utils import safe_write

    _, path = base_h5ad
    adata = sc.AnnData(X=np.random.default_rng(7).random((80, 5)))
    adata.obs["x"] = np.arange(80)
    safe_write(adata, path)
    assert sc.read(path).n_obs == 80
