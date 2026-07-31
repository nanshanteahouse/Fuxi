"""Crash-recovery drills — plan h5ad-incremental-io task 11 ③.

Simulates interrupted appends (failure injection via monkeypatch / bad input —
never by killing real processes) and proves the two failure policies plus the
integrity gate:

a) **in-place mode** (``write_obs_columns_inplace``): a failed append must
   restore the shared checkpoint atomically from its ``.bak`` — byte-identical
   to the pre-write file. Two injections: the engine's own obs-length
   ``ValueError`` (rejected before any write) and a monkeypatched engine that
   leaks a partial dataset into the file before raising.
b) **copy+append mode** (``write_obs_columns_lightweight``): a failed append
   deletes the corrupt target; the intact source file survives untouched
   ("失败即删文件" — re-running the step regenerates the copy).
c) **integrity gate** (``verify_incremental_write``) intercepts hand-corrupted
   ("half-written") files: a missing dataset (truncated write) and a torn
   dataset (wrong values) both raise ``RuntimeError``.

File-level corruption detection (c) is what the function-level backup tests in
``test_io_incremental.py`` did NOT cover — this file closes that gap.
"""

from __future__ import annotations

import os
import shutil

import h5py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from core.utils import (
    write_h5ad_incremental,
    write_obs_columns_inplace,
    write_obs_columns_lightweight,
)
from core.utils._io_incremental import verify_incremental_write


@pytest.fixture
def base_h5ad(tmp_path):
    """A small gzip h5ad; returns (base_adata, path)."""
    rng = np.random.default_rng(3)
    adata = sc.AnnData(X=rng.random((300, 20)).astype(np.float32))
    adata.obs["batch"] = pd.Categorical(["a", "b", "c"] * 100)
    adata.obsm["X_pca"] = rng.random((300, 10))
    path = tmp_path / "base.h5ad"
    adata.write(str(path), compression="gzip")
    return adata, str(path)


# ── a) in-place mode: .bak atomic restore ──────────────────────────────


def test_inplace_restores_on_obs_length_mismatch(base_h5ad):
    """Engine's own length check (ValueError) → .bak restores original bytes.

    The length mismatch is rejected before any dataset write, so the restore is
    trivially byte-exact — but it still proves the writeback wrapper's atomic
    restore contract holds on the engine-native error path.
    """
    _, path = base_h5ad
    with open(path, "rb") as fh:
        original = fh.read()
    wrong_length = pd.DataFrame({"cell_subtype": np.arange(10)})  # file has 300 obs

    with pytest.raises(ValueError, match="does not match"):
        write_obs_columns_inplace(path, wrong_length)

    assert not os.path.exists(path + ".bak"), "backup must be consumed by the restore"
    with open(path, "rb") as fh:
        assert fh.read() == original, "file must be restored byte-identical"
    assert "cell_subtype" not in sc.read(path).obs, "no partial column may survive"


def test_inplace_restores_on_mid_append_corruption(base_h5ad, monkeypatch):
    """A crash mid-append (partial dataset leaked, then raise) → full restore.

    Simulates the HDF5 no-WAL worst case: a partially appended dataset is left
    inside the file. The wrapper must replace the whole file from .bak so the
    leaked dataset disappears with it.
    """
    import core.utils._io as io_mod

    _, path = base_h5ad
    with open(path, "rb") as fh:
        original = fh.read()

    def _explode(h5ad_path, obs=None, obsm=None, obsp=None, uns=None, logger=None):
        with h5py.File(h5ad_path, "a") as f:
            f["obs/evil_partial"] = np.arange(5)  # leaked partial dataset
        raise RuntimeError("simulated mid-append crash")

    monkeypatch.setattr(io_mod, "write_h5ad_incremental", _explode)
    with pytest.raises(RuntimeError, match="simulated mid-append crash"):
        write_obs_columns_inplace(path, pd.DataFrame({"cell_subtype": np.arange(300)}))

    assert not os.path.exists(path + ".bak")
    with open(path, "rb") as fh:
        assert fh.read() == original, "file must be restored byte-identical"
    with h5py.File(path, "r") as f:
        assert "evil_partial" not in f["obs"], "leaked dataset must be gone after restore"


# ── b) copy+append mode: delete corrupt copy, keep source ──────────────


def test_copy_append_failure_deletes_target_keeps_source(base_h5ad, monkeypatch):
    """Mode A: failed append removes the derived copy; source stays byte-identical."""
    import core.utils._io as io_mod

    _, src = base_h5ad
    target = os.path.join(os.path.dirname(src), "05_annotated.h5ad")
    shutil.copy2(src, target)  # the step's copy2-from-source step
    with open(src, "rb") as fh:
        original_src = fh.read()

    def _explode(h5ad_path, obs=None, obsm=None, obsp=None, uns=None, logger=None):
        raise RuntimeError("simulated append crash")

    monkeypatch.setattr(io_mod, "write_h5ad_incremental", _explode)
    with pytest.raises(RuntimeError, match="simulated append crash"):
        write_obs_columns_lightweight(target, pd.DataFrame({"ct": np.arange(300)}))

    assert not os.path.exists(target), "corrupt copy must be deleted (re-run regenerates)"
    with open(src, "rb") as fh:
        assert fh.read() == original_src, "source file must be untouched"


# ── c) integrity gate: verify_incremental_write intercepts corruption ──


def test_verify_detects_missing_dataset(base_h5ad):
    """Half-written file: the written dataset is gone → verify raises."""
    _, path = base_h5ad
    df = pd.DataFrame({"cell_type": pd.Categorical(["t1", "t2", "t3"] * 100)})
    write_h5ad_incremental(path, obs=df)

    with h5py.File(path, "a") as f:
        del f["obs/cell_type"]  # simulate a write that never completed the dataset

    with pytest.raises(RuntimeError, match="missing from /obs"):
        verify_incremental_write(path, obs=df)


def test_verify_detects_corrupted_values(base_h5ad):
    """Half-written file: dataset exists but holds torn values → verify raises."""
    _, path = base_h5ad
    df = pd.DataFrame({"cell_type": pd.Categorical([f"c{i % 4}" for i in range(300)])})
    write_h5ad_incremental(path, obs=df)

    with h5py.File(path, "a") as f:
        codes = f["obs/cell_type/codes"]
        codes[...] = np.zeros_like(codes[...])  # torn write: wrong codes on disk

    with pytest.raises(RuntimeError, match="does not round-trip"):
        verify_incremental_write(path, obs=df)


def test_verify_passes_on_intact_append(base_h5ad):
    """Positive control: an untouched append passes the integrity gate."""
    _, path = base_h5ad
    rng = np.random.default_rng(1)
    df = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["t1", "t2", "t3"] * 100),
            "lab": pd.array([f"x{i}" if i % 4 else None for i in range(300)], dtype="string"),
        }
    )
    payload = {"obs": df, "obsm": {"X_umap": rng.random((300, 2))}}
    write_h5ad_incremental(path, **payload)
    verify_incremental_write(path, **payload)  # must not raise
