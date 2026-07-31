"""Tests for copy_h5ad — opportunistic reflink (FICLONE) with copy2 fallback.

Covers the three branches of :func:`core.utils._io.copy_h5ad`:
  1. real-filesystem correctness — byte-identical + anndata-readable copy
     on whatever fs the tmpdir lives on (ext4/tmpfs usually lack FICLONE,
     so this exercises the copy2 fallback; on btrfs/xfs it exercises the
     reflink path);
  2. mocked ioctl success → the reflink branch runs and copy2 is never
     reached (branch logic, no byte verification required by spec);
  3. mocked ioctl OSError → silent fallback to copy2 with an identical file;
  4. no ``fcntl`` (non-POSIX) → plain copy2.

The ``.bak`` backup path of write_obs_columns_inplace is covered by the
existing regression tests in tests/test_core/test_io_incremental.py.
"""

from __future__ import annotations

import builtins
import errno
import fcntl
import os
import shutil

import numpy as np
import pytest
import scanpy as sc

from core.utils._io import FICLONE, copy_h5ad


@pytest.fixture
def src_h5ad(tmp_path):
    """A small gzip h5ad source file; returns (base_adata, path)."""
    rng = np.random.default_rng(0)
    adata = sc.AnnData(X=rng.random((120, 15)))
    adata.obs["batch"] = [f"b{i % 3}" for i in range(120)]
    adata.obsm["X_pca"] = rng.random((120, 5))
    path = tmp_path / "src.h5ad"
    adata.write(str(path), compression="gzip")
    return adata, path


def test_copy_h5ad_byte_identical_and_readable(src_h5ad, tmp_path):
    """Real-fs copy is byte-identical and opens with scanpy (reflink or copy2)."""
    _, src = src_h5ad
    dst = tmp_path / "dst.h5ad"
    copy_h5ad(str(src), str(dst))
    assert dst.read_bytes() == src.read_bytes()
    got = sc.read(str(dst))
    assert got.shape == (120, 15)


def test_copy_h5ad_takes_reflink_path_when_ioctl_succeeds(src_h5ad, tmp_path, monkeypatch):
    """Mocked successful ioctl → FICLONE branch runs; copy2 is never reached."""
    _, src = src_h5ad
    dst = tmp_path / "dst.h5ad"
    copy2_calls: list[tuple] = []
    monkeypatch.setattr(shutil, "copy2", lambda s, d: copy2_calls.append((s, d)))

    def _fake_ioctl(fd, request, arg):
        assert request == FICLONE
        assert isinstance(arg, int)  # src fd passed as the ioctl argument
        fd.write(src.read_bytes())  # emulate a real clone into the dst fd
        return 0

    monkeypatch.setattr(fcntl, "ioctl", _fake_ioctl)
    copy_h5ad(str(src), str(dst))
    assert copy2_calls == []
    assert dst.read_bytes() == src.read_bytes()


def test_copy_h5ad_falls_back_to_copy2_on_ioctl_failure(src_h5ad, tmp_path, monkeypatch):
    """ioctl OSError → silent copy2 fallback with identical bytes + mode."""
    _, src = src_h5ad
    dst = tmp_path / "dst.h5ad"

    def _fail_ioctl(*_args):
        raise OSError(errno.EINVAL, "FICLONE not supported on this filesystem")

    monkeypatch.setattr(fcntl, "ioctl", _fail_ioctl)
    copy_h5ad(str(src), str(dst))
    assert dst.read_bytes() == src.read_bytes()
    assert os.stat(dst).st_mode == os.stat(src).st_mode


def test_copy_h5ad_without_fcntl_falls_back(src_h5ad, tmp_path, monkeypatch):
    """No fcntl module (non-POSIX) → plain copy2, content identical."""
    _, src = src_h5ad  # build the source before patching __import__
    dst = tmp_path / "dst.h5ad"
    real_import = builtins.__import__

    def _no_fcntl(name, *args, **kwargs):
        if name == "fcntl":
            raise ImportError("no fcntl (non-POSIX platform)")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_fcntl)
    copy_h5ad(str(src), str(dst))
    assert dst.read_bytes() == src.read_bytes()
