"""Compression resolution tests for ``core.utils.safe_write``.

Regression for the per-step compression clobbering bug (metis G1): the value
resolved from ``cfg.per_step_h5ad_compression[step_alias]`` was overwritten by
``cfg.h5ad_compression`` when ``compression_override`` was ``None``.

Resolution order under test (highest first):
    ``compression_override`` > per-step value (when ``step_alias`` is given
    and present in ``per_step_h5ad_compression``) > ``cfg.h5ad_compression``
    > caller default ("gzip").

All helpers are self-contained; no conftest fixtures required.
"""

from types import SimpleNamespace

import h5py
import numpy as np
import scanpy as sc
import scipy.sparse as sp

from core.config.schema import Config
from core.utils import _io as io_mod
from core.utils import safe_write


def _make_cfg(h5ad_compression="gzip", per_step=None, opts=None, tmpdir=None):
    """Minimal stand-in for the Pydantic Config (only the attrs safe_write reads)."""
    per_step = {} if per_step is None else per_step
    return SimpleNamespace(
        h5ad_compression=h5ad_compression,
        per_step_h5ad_compression=dict(per_step),
        h5ad_compression_opts=opts,
        h5ad_tempdir=str(tmpdir),
        verify_write_integrity=False,
    )


def _write(path, cfg, *, step_alias=None, compression_override=None, compression="gzip"):
    """Write a tiny sparse AnnData via safe_write and return the target path."""
    adata = sc.AnnData(X=sp.csr_matrix(np.arange(12, dtype=np.float32).reshape(3, 4)))
    safe_write(
        adata,
        path,
        cfg=cfg,
        step_alias=step_alias,
        compression_override=compression_override,
        compression=compression,
    )
    return path


def _x_compression(path):
    """Compression family actually used for the stored /X/data dataset."""
    with h5py.File(path, "r") as f:
        dset = f["X/data"]
        if "32015" in dset._filters:  # zstd filter registered by hdf5plugin
            return "zstd"
        return dset.compression


def _assert_resolved(path, expected):
    """Assert the resolved compression for *expected*, env-adaptive for zstd.

    When zstd is chosen but hdf5plugin is absent, ``safe_write`` falls back to
    plain gzip at h5py's DEFAULT level (no level-1 forcing) — so "zstd" is
    accepted as the gzip fallback in that environment.
    """
    got = _x_compression(path)
    if expected == "zstd" and io_mod.hdf5plugin is None:
        assert got == "gzip", f"zstd fallback should be gzip, got {got!r}"
        return
    assert got == expected, f"expected {expected!r}, got {got!r}"


def _has_zstd():
    return io_mod.hdf5plugin is not None


# ── Failing-first regression: per-step zstd must beat global gzip ──────
def test_per_step_zstd_wins_over_global_gzip(tmp_path):
    """per_step_h5ad_compression={"raw": "zstd"} + step_alias="raw" must NOT
    be clobbered back to the global h5ad_compression="gzip"."""
    cfg = _make_cfg(h5ad_compression="gzip", per_step={"raw": "zstd"}, tmpdir=tmp_path)
    path = _write(str(tmp_path / "raw.h5ad"), cfg, step_alias="raw")
    _assert_resolved(path, "zstd")


# ── Resolution order: the three precedence levels ──────────────────────
def test_compression_override_beats_per_step(tmp_path):
    """compression_override (highest) beats a present per-step value."""
    cfg = _make_cfg(per_step={"raw": "zstd"}, tmpdir=tmp_path)
    path = _write(
        str(tmp_path / "override.h5ad"), cfg, step_alias="raw", compression_override="gzip"
    )
    _assert_resolved(path, "gzip")


def test_per_step_beats_global_cfg(tmp_path):
    """per-step value beats cfg.h5ad_compression."""
    cfg = _make_cfg(h5ad_compression="lzf", per_step={"raw": "zstd"}, tmpdir=tmp_path)
    path = _write(str(tmp_path / "perstep.h5ad"), cfg, step_alias="raw")
    _assert_resolved(path, "zstd")


def test_global_cfg_beats_caller_default(tmp_path):
    """cfg.h5ad_compression beats the caller's default compression="gzip"."""
    cfg = _make_cfg(h5ad_compression="lzf", tmpdir=tmp_path)
    path = _write(str(tmp_path / "global.h5ad"), cfg)
    _assert_resolved(path, "lzf")


def test_step_alias_absent_in_per_step_falls_through_to_global(tmp_path):
    """step_alias not present in per_step → cfg.h5ad_compression wins."""
    cfg = _make_cfg(h5ad_compression="lzf", per_step={"raw": "zstd"}, tmpdir=tmp_path)
    path = _write(str(tmp_path / "absent.h5ad"), cfg, step_alias="integrated")
    _assert_resolved(path, "lzf")


def test_cfg_absent_uses_caller_default(tmp_path):
    """No cfg → caller default ("gzip") is used."""
    path = _write(str(tmp_path / "default.h5ad"), cfg=None)
    _assert_resolved(path, "gzip")


# ── zstd fallback without hdf5plugin (metis G9) ─────────────────────────
def test_zstd_falls_back_to_gzip_without_hdf5plugin(tmp_path, monkeypatch):
    """zstd chosen + hdf5plugin absent → plain gzip (no level assertion)."""
    monkeypatch.setattr(io_mod, "hdf5plugin", None)
    cfg = _make_cfg(h5ad_compression="gzip", per_step={"raw": "zstd"}, tmpdir=tmp_path)
    path = _write(str(tmp_path / "fallback.h5ad"), cfg, step_alias="raw")
    _assert_resolved(path, "zstd")  # maps to the gzip-fallback branch above


# ── Integration: inspect the stored dataset's filter ───────────────────
def test_integration_stored_dataset_uses_zstd_filter(tmp_path):
    """safe_write(..., cfg=cfg, step_alias="raw") stores X/data with the zstd
    filter (id 32015) when hdf5plugin is importable; gzip fallback otherwise."""
    cfg = _make_cfg(h5ad_compression="gzip", per_step={"raw": "zstd"}, tmpdir=tmp_path)
    path = str(tmp_path / "integrated.h5ad")
    adata = sc.AnnData(X=sp.csr_matrix(np.arange(12, dtype=np.float32).reshape(3, 4)))
    safe_write(adata, path, cfg=cfg, step_alias="raw")

    with h5py.File(path, "r") as f:
        dset = f["X/data"]
        if _has_zstd():
            assert "32015" in dset._filters, f"expected zstd filter 32015, got {dset._filters}"
            assert dset.compression != "gzip"
        else:
            # env-adaptive: no hdf5plugin → plain gzip fallback, default level
            assert dset.compression == "gzip"
            assert "32015" not in dset._filters


# ── Schema default ──────────────────────────────────────────────────────
def test_schema_default_per_step_contains_raw():
    """per_step_h5ad_compression schema default gains "raw": "zstd"."""
    per_step = Config().per_step_h5ad_compression
    assert per_step.get("raw") == "zstd"
    assert per_step.get("integrated") == "gzip"
