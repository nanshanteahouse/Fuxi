"""Tests for core.utils._write_chunked — shared chunked CSR → h5ad writer.

Self-contained: all helpers are defined here (no shared fixtures), mirroring the
"random CSR 50k×2k @5% density + boolean cell_mask" acceptance in the plan.
Reference oracle: a naive scipy row selection (``X[mask]``) — the chunked writer
must produce byte-identical data/indices/indptr.
"""

import glob
import os

import h5py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

from core.utils._write_chunked import resolve_block_size, write_csr_chunked


# ── helpers ─────────────────────────────────────────────────────────────
def _make_csr(n_cells, n_genes, density, seed=0):
    """Random CSR with the requested density and no duplicate (row, col) coords."""
    rng = np.random.default_rng(seed)
    target = max(1, int(n_cells * n_genes * density))
    total = n_cells * n_genes
    frac = min(1.0, target / total)
    if frac >= 0.95:  # near-dense: build directly
        row, col = np.indices((n_cells, n_genes))
        flat = (row * n_genes + col).ravel()
    else:
        # Inverse coupon-collector: sample n_try WITH replacement so the
        # expected unique count (M·(1-e^(-k/M))) comfortably exceeds target.
        n_try = min(total, int(-np.log1p(-frac) * total * 1.05) + 64)
        row = rng.integers(0, n_cells, size=n_try, dtype=np.int64)
        col = rng.integers(0, n_genes, size=n_try, dtype=np.int64)
        flat = np.unique(row * n_genes + col)
        if len(flat) < target:
            raise RuntimeError("could not generate enough unique (row, col) entries")
        flat = np.sort(rng.choice(flat, size=target, replace=False))
    row, col = np.unravel_index(flat, (n_cells, n_genes))
    data = rng.uniform(0.5, 20.0, size=target).astype(np.float32)
    return sp.csr_matrix((data, (row, col)), shape=(n_cells, n_genes))


def _make_mask(n, keep_frac=0.5, seed=0, all_false=False, all_true=False):
    if all_false:
        return np.zeros(n, dtype=bool)
    if all_true:
        return np.ones(n, dtype=bool)
    rng = np.random.default_rng(seed)
    return rng.random(n) < keep_frac


def _assert_csr_equal(actual, expected):
    """Exact equality of shape/data/indices/indptr (NaN-aware on data)."""
    assert actual.shape == expected.shape
    assert np.array_equal(actual.data, expected.data, equal_nan=True)
    assert np.array_equal(actual.indices, expected.indices)
    assert np.array_equal(actual.indptr, expected.indptr)


def _dataset_filter_id(dset):
    """Return the first h5py filter id (1=gzip, 32015=zstd) or None."""
    plist = dset.id.get_create_plist()
    for i in range(plist.get_nfilters()):
        fid, _, _, _ = plist.get_filter(i)
        if fid:
            return fid
    return None


# ── block boundaries ────────────────────────────────────────────────────
def test_roundtrip_smaller_than_block(tmp_path):
    x = _make_csr(1_000, 300, density=0.2, seed=1)
    target = tmp_path / "small.h5ad"
    write_csr_chunked(str(target), x, block_size=50_000)
    ad = sc.read_h5ad(str(target))
    assert ad.shape == x.shape
    _assert_csr_equal(ad.X, x)


def test_roundtrip_larger_than_block(tmp_path):
    # 120_000 rows / 50_000 block → 3 blocks, last block partial.
    x = _make_csr(120_000, 500, density=0.05, seed=2)
    target = tmp_path / "big.h5ad"
    write_csr_chunked(str(target), x, block_size=50_000)
    ad = sc.read_h5ad(str(target))
    assert ad.shape == x.shape
    _assert_csr_equal(ad.X, x)


def test_roundtrip_auto_block_size(tmp_path):
    # block_size=None → auto-resolved from density; small matrix → one block.
    x = _make_csr(5_000, 1_000, density=0.1, seed=11)
    target = tmp_path / "auto.h5ad"
    write_csr_chunked(str(target), x)
    ad = sc.read_h5ad(str(target))
    _assert_csr_equal(ad.X, x)


# ── cell_mask variants ──────────────────────────────────────────────────
def test_roundtrip_no_mask(tmp_path):
    x = _make_csr(30_000, 500, density=0.1, seed=3)
    target = tmp_path / "nomask.h5ad"
    write_csr_chunked(str(target), x, cell_mask=None, block_size=10_000)
    ad = sc.read_h5ad(str(target))
    _assert_csr_equal(ad.X, x)


def test_roundtrip_partial_mask(tmp_path):
    # 50k×2k @5% (plan acceptance) with a 40%-keep mask, 3 blocks of 20k.
    x = _make_csr(50_000, 2_000, density=0.05, seed=4)
    mask = _make_mask(x.shape[0], keep_frac=0.4, seed=7)
    target = tmp_path / "masked.h5ad"
    write_csr_chunked(str(target), x, cell_mask=mask, block_size=20_000)
    ad = sc.read_h5ad(str(target))
    assert ad.shape == (int(mask.sum()), x.shape[1])
    _assert_csr_equal(ad.X, x[mask])  # naive scipy selection = reference


def test_roundtrip_all_true_mask(tmp_path):
    x = _make_csr(2_000, 300, density=0.2, seed=5)
    mask = _make_mask(x.shape[0], all_true=True)
    target = tmp_path / "alltrue.h5ad"
    write_csr_chunked(str(target), x, cell_mask=mask, block_size=50_000)
    ad = sc.read_h5ad(str(target))
    _assert_csr_equal(ad.X, x)


def test_roundtrip_all_false_mask(tmp_path):
    # empty-mask edge: 0-row matrix written, read succeeds.
    x = _make_csr(2_000, 300, density=0.2, seed=6)
    mask = _make_mask(x.shape[0], all_false=True)
    target = tmp_path / "allfalse.h5ad"
    write_csr_chunked(str(target), x, cell_mask=mask, block_size=50_000)
    ad = sc.read_h5ad(str(target))
    assert ad.shape == (0, x.shape[1])
    assert ad.X.nnz == 0


def test_all_zero_block_keeps_rows(tmp_path):
    # Entire first block (50k rows) is empty — its rows must still appear in
    # indptr so the output row count stays exact.
    n, m, block = 150_000, 100, 50_000
    x = _make_csr(n, m, density=0.01, seed=8)
    empty = sp.csr_matrix((block, m), dtype=np.float32)
    x = sp.vstack([empty, x[block:]]).tocsr()
    target = tmp_path / "zeroblock.h5ad"
    write_csr_chunked(str(target), x, block_size=block)
    ad = sc.read_h5ad(str(target))
    assert ad.shape == x.shape
    _assert_csr_equal(ad.X, x)


# ── metadata + NaN semantics ────────────────────────────────────────────
def test_obs_var_written(tmp_path):
    n, m = 500, 100
    x = _make_csr(n, m, 0.1, seed=9)
    obs = pd.DataFrame({"sample": ["a"] * n}, index=[f"cell{i}" for i in range(n)])
    var = pd.DataFrame(
        {"n_cells": np.arange(m, dtype=np.int64)}, index=[f"gene{i}" for i in range(m)]
    )
    target = tmp_path / "with_meta.h5ad"
    write_csr_chunked(str(target), x, obs=obs, var=var, block_size=200)
    ad = sc.read_h5ad(str(target))
    assert list(ad.obs_names) == list(obs.index)
    assert list(ad.var_names) == list(var.index)
    _assert_csr_equal(ad.X, x)


def test_roundtrip_preserves_nan(tmp_path):
    x = _make_csr(2_000, 400, density=0.1, seed=13)
    x.data[:100] = np.nan  # inject NaN; exact-0.0 stays (matches csr semantics)
    target = tmp_path / "nan.h5ad"
    write_csr_chunked(str(target), x, block_size=50_000)
    ad = sc.read_h5ad(str(target))
    _assert_csr_equal(ad.X, x)


def test_csr_encoding_attrs(tmp_path):
    x = _make_csr(1_000, 200, density=0.1, seed=15)
    target = tmp_path / "attrs.h5ad"
    write_csr_chunked(str(target), x, block_size=50_000)
    with h5py.File(str(target), "r") as f:
        assert f["X"].attrs["encoding-type"] == "csr_matrix"
        assert tuple(f["X"].attrs["shape"]) == x.shape
        assert f["X"]["data"].chunks == (65536,)
        assert f["X"]["indices"].chunks == (65536,)


# ── tmp cleanup + compression ───────────────────────────────────────────
def test_tmp_cleanup_on_failure(tmp_path):
    x = _make_csr(500, 100, density=0.1, seed=16)
    target = tmp_path / "out.h5ad"
    with pytest.raises(Exception):
        write_csr_chunked(str(target), x, compression="not-a-real-compressor")
    assert not os.path.exists(str(target))
    assert glob.glob(str(tmp_path / ".out.h5ad.tmp.*")) == []


def test_compression_gzip_passthrough(tmp_path):
    x = _make_csr(1_000, 200, density=0.1, seed=17)
    target = tmp_path / "gzip.h5ad"
    write_csr_chunked(str(target), x, compression="gzip", compression_opts=1, block_size=50_000)
    with h5py.File(str(target), "r") as f:
        for key in ("data", "indices", "indptr"):
            d = f["X"][key]
            assert d.compression == "gzip", key
            assert d.compression_opts == 1, key


def test_compression_zstd_passthrough(tmp_path):
    pytest.importorskip("hdf5plugin")
    x = _make_csr(1_000, 200, density=0.1, seed=18)
    target = tmp_path / "zstd.h5ad"
    write_csr_chunked(str(target), x, compression="zstd", block_size=50_000)
    with h5py.File(str(target), "r") as f:
        assert _dataset_filter_id(f["X"]["data"]) == 32015  # H5Z-ZSTD


def test_verify_false_writes_without_check(tmp_path):
    x = _make_csr(1_000, 100, density=0.1, seed=19)
    target = tmp_path / "noverify.h5ad"
    write_csr_chunked(str(target), x, verify=False, block_size=50_000)
    ad = sc.read_h5ad(str(target))
    _assert_csr_equal(ad.X, x)


# ── resolve_block_size (pure function) ──────────────────────────────────
def test_resolve_block_size_math_and_clamp():
    # Realistic density at 32 GiB and 98 GiB avail → both clamp to 500k
    # (mirrors the 02_qc anchor expectations for a large sparse matrix).
    for gb in (32, 98):
        bs = resolve_block_size((2_000_000, 20_000), 2_000_000 * 200, gb * 2**30, prefetch=1)
        assert bs == 500_000
    # Dense + huge avail → upper clamp.
    assert resolve_block_size((1_000, 1_000), 1_000_000, 100 * 2**30) == 500_000
    # Dense + tiny avail → lower clamp.
    assert resolve_block_size((1_000, 1_000), 1_000_000, 1 << 20) == 50_000
    # Zero-nnz matrix → max_block, no ZeroDivisionError.
    assert resolve_block_size((100, 100), 0, 32 * 2**30) == 500_000
    # Monotonic in avail_bytes.
    small = resolve_block_size((1_000_000, 10_000), 10_000_000, 8 * 2**30)
    large = resolve_block_size((1_000_000, 10_000), 10_000_000, 64 * 2**30)
    assert small <= large


def test_resolve_block_size_exact_math():
    # density = 1500 nnz/row, 32 GiB avail → raw formula value, no clamp hit.
    avail = 32 * 2**30
    bs = resolve_block_size((2_000_000, 2_000), 2_000_000 * 1_500, avail, prefetch=1)
    expected = int(avail * 0.4 / (1500 * 12 * 2))
    assert bs == expected
    assert 50_000 < expected < 500_000
