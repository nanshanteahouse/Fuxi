"""DrvFs (/mnt) incremental-write stability tests — plan h5ad-incremental-io task 11 ①.

These tests run against the REAL WSL /mnt filesystem at ``$FUXI_DATA_ROOT``
(typically ``/mnt/e/neurobiology``, a DrvFs/9p mount). They are the dedicated
DrvFs QA requested by the plan's risk table ("WSL /mnt DrvFs 下 h5py 'a' mode
高风险 → 优先专项测试"): h5py ``'a'``-mode appends must behave identically to
local filesystem appends, and ``copy_h5ad`` must silently fall back to
``shutil.copy2`` for the ``.bak`` backup (DrvFs has no FICLONE/reflink).

Skip policy: skipped ONLY when ``FUXI_DATA_ROOT`` is unset or not writable —
deliberately NOT gated by ``SKIP_SLOW_TESTS`` (this is the DrvFs acceptance
suite).

Cleanup contract: every test creates ``$FUXI_DATA_ROOT/_h5ad_incr_qa_<pid>``
and removes it in ``finally``, asserting the directory is gone afterwards.
``test_no_leftover_qa_dirs`` fails loudly if any earlier crashed run leaked one.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from core.utils import write_h5ad_incremental, write_obs_columns_inplace
from core.utils._io import FICLONE

_N_QA_PREFIX = "_h5ad_incr_qa_"


def _drvfs_root() -> Path | None:
    """Return the DrvFs root, or None when the test must be skipped."""
    root = os.environ.get("FUXI_DATA_ROOT")
    if not root:
        return None
    p = Path(root)
    if not p.is_dir():
        return None
    probe = p / f"{_N_QA_PREFIX}writeprobe_{os.getpid()}"
    try:
        probe.write_bytes(b"probe")
        probe.unlink()
    except OSError:
        return None
    return p


@pytest.fixture
def drvfs_tmpdir():
    """A fresh temp dir under $FUXI_DATA_ROOT, always removed afterwards."""
    root = _drvfs_root()
    if root is None:
        pytest.skip("FUXI_DATA_ROOT unset / not writable — DrvFs QA requires /mnt")
    tmp = root / f"{_N_QA_PREFIX}{os.getpid()}"
    tmp.mkdir(parents=True, exist_ok=False)
    try:
        yield tmp
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        assert not tmp.exists(), f"DrvFs QA temp dir leaked: {tmp}"


def _make_base_h5ad(path: Path) -> sc.AnnData:
    """Write a small gzip-compressed h5ad (RNA layout) and return the AnnData."""
    rng = np.random.default_rng(7)
    adata = sc.AnnData(X=rng.random((200, 40)).astype(np.float32))
    adata.obs["batch"] = pd.Categorical([f"b{i % 4}" for i in range(200)])
    adata.obsm["X_pca"] = rng.random((200, 10))
    adata.write(str(path), compression="gzip")
    return adata


def _append_payload(rng) -> dict:
    """The identical incremental payload applied to local + DrvFs copies."""
    n = 200
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical([f"T{i % 3}" for i in range(n)]),
            "lab": pd.array([f"m{i}" if i % 5 else None for i in range(n)], dtype="string"),
            "score": np.linspace(0.0, 1.0, n),
        },
        index=range(n),
    )
    return {
        "obs": obs,
        "obsm": {"X_umap": rng.random((n, 2))},
        "uns": {"leiden_params": {"resolution": 0.8, "n_neighbors": 15}},
    }


def _assert_h5ad_equivalent(local: sc.AnnData, drvfs: sc.AnnData, payload: dict) -> None:
    """Local vs DrvFs read-back must be identical on every written key + X."""
    assert set(drvfs.obs.columns) == set(local.obs.columns), (
        f"obs columns differ on DrvFs: {set(drvfs.obs.columns) ^ set(local.obs.columns)}"
    )
    for col in payload["obs"].columns:
        lv, dv = local.obs[col], drvfs.obs[col]
        assert list(lv.astype(str)) == list(dv.astype(str)), f"obs '{col}' differs on DrvFs"
        assert str(lv.dtype) == str(dv.dtype), f"obs '{col}' dtype differs on DrvFs"
    assert np.array_equal(local.obsm["X_umap"], drvfs.obsm["X_umap"]), "X_umap differs"
    assert drvfs.uns["leiden_params"] == payload["uns"]["leiden_params"]
    assert np.array_equal(np.asarray(local.X), np.asarray(drvfs.X)), "X must be untouched"
    assert np.array_equal(local.obsm["X_pca"], drvfs.obsm["X_pca"]), "X_pca must be untouched"


def test_drvfs_incremental_append_matches_local(drvfs_tmpdir, tmp_path):
    """h5py 'a'-mode append on DrvFs round-trips identically to the local twin.

    Same base h5ad (gzip) → same incremental payload (obs categorical + nullable
    string + float, obsm, uns) → read-back equality on values, dtypes and
    untouched X/obsm. This is the risk-table item: RNA gzip h5ad appends must be
    stable on WSL /mnt.
    """
    local_base = tmp_path / "base.h5ad"
    drvfs_base = drvfs_tmpdir / "base.h5ad"
    _make_base_h5ad(local_base)
    _make_base_h5ad(drvfs_base)

    payload = _append_payload(np.random.default_rng(11))
    write_h5ad_incremental(str(local_base), logger=None, **payload)
    write_h5ad_incremental(str(drvfs_base), logger=None, **payload)

    local = sc.read(str(local_base))
    drvfs = sc.read(str(drvfs_base))
    _assert_h5ad_equivalent(local, drvfs, payload)

    # Appends on the gzip DrvFs file must keep gzip compression (uniform layout).
    with h5py.File(drvfs_base, "r") as f:
        assert f["X"].compression == "gzip"
        assert f["obs/cell_type/codes"].compression == "gzip"
        assert f["obsm/X_umap"].compression == "gzip"


def test_drvfs_inplace_writeback_bak_falls_back_to_copy2(drvfs_tmpdir, monkeypatch):
    """write_obs_columns_inplace on DrvFs: .bak copy silently falls back to copy2.

    DrvFs does not support FICLONE (verified by an ioctl probe against a scratch
    file in the same temp dir), so ``copy_h5ad`` must take the copy2 fallback for
    the ``.bak`` — proven by recording every ``shutil.copy2`` call whose dst is
    the backup. Success path then removes the backup and leaves X untouched.
    """
    base = drvfs_tmpdir / "05_annotated.h5ad"
    _make_base_h5ad(base)
    x_before = np.asarray(sc.read(str(base)).X)
    before = base.read_bytes()

    # ── Probe: FICLONE is unsupported on this DrvFs dir (the premise of copy2) ──
    probe_src, probe_dst = drvfs_tmpdir / "p_src", drvfs_tmpdir / "p_dst"
    probe_src.write_bytes(b"x" * 8192)
    import fcntl

    try:
        with open(probe_src, "rb") as s, open(probe_dst, "wb") as d:
            fcntl.ioctl(d, FICLONE, s.fileno())
        reflink_supported = True
    except OSError:
        reflink_supported = False
    finally:
        probe_src.unlink(missing_ok=True)
        probe_dst.unlink(missing_ok=True)

    # ── Record .bak copies made through the real copy2 (not reflink) ──────────
    bak_copies: list[tuple[str, str]] = []
    real_copy2 = shutil.copy2

    def _counting_copy2(src, dst, *args, **kwargs):
        if str(dst).endswith(".bak"):
            bak_copies.append((str(src), str(dst)))
        return real_copy2(src, dst, *args, **kwargs)

    monkeypatch.setattr(shutil, "copy2", _counting_copy2)

    n = 200
    write_obs_columns_inplace(
        str(base),
        pd.DataFrame({"cell_subtype": pd.Categorical([f"s{i % 6}" for i in range(n)])}),
    )

    # Reflink unsupported → the .bak was produced by copy2 exactly once.
    if reflink_supported:
        assert not bak_copies, "reflink path should have avoided copy2"
    else:
        assert len(bak_copies) == 1, f"expected exactly one copy2 .bak, got {bak_copies}"

    # Success contract: backup removed, column appended, X untouched.
    assert not (drvfs_tmpdir / "05_annotated.h5ad.bak").exists()
    after = sc.read(str(base))
    assert "cell_subtype" in after.obs.columns
    assert np.array_equal(np.asarray(after.X), x_before), "X must be untouched"
    assert base.read_bytes() != before  # the in-place append changed the file


def test_no_leftover_qa_dirs():
    """No ``_h5ad_incr_qa_*`` dirs may remain under $FUXI_DATA_ROOT.

    Guards the plan's cleanup acceptance item: a crashed test run that leaked
    its temp dir must fail loudly on the next suite execution.
    """
    root = _drvfs_root()
    if root is None:
        pytest.skip("FUXI_DATA_ROOT unset / not writable — DrvFs QA requires /mnt")
    leftovers = sorted(p for p in root.glob(f"{_N_QA_PREFIX}*") if p.is_dir())
    assert not leftovers, f"leaked DrvFs QA dirs from earlier runs: {leftovers}"
