"""Tests for the ambient-correction block in rna/steps/00_load.py.

Covers
------
- ``_run_ambient_correction`` with ``cfg.ambient.run=True`` writes
  ``ambient_removed.h5ad`` without raising TypeError for both
  ``cellbender`` and ``soupx`` methods (regression for the dropped
  invalid ``file_type="h5ad"`` kwarg on :func:`core.utils.safe_write`).
- missing dependency → logged skip, no crash, no file written (existing
  pass-through behaviour preserved).

The ambient ``safe_write`` call really executes: a tiny synthetic AnnData
is written to a real h5ad in ``tmp_path`` and read back.
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_ambient_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════


def _make_cfg(tmp_path: Path, method: str) -> MagicMock:
    """Minimal config with concrete values safe_write needs (no mocks inside)."""
    cfg = MagicMock()
    cfg.ambient.run = True
    cfg.ambient.method = method
    cfg.raw_h5ad = str(tmp_path / "raw.h5ad")
    # concrete safe_write inputs — MagicMock defaults would poison the write
    cfg.per_step_h5ad_compression = {}
    cfg.h5ad_compression = "gzip"
    cfg.h5ad_compression_opts = None
    cfg.h5ad_tempdir = str(tmp_path)
    cfg.verify_write_integrity = False
    return cfg


def _make_adata() -> anndata.AnnData:
    rng = np.random.RandomState(7)
    x = csr_matrix(rng.randint(0, 5, size=(3, 4)).astype(np.float32))
    obs = pd.DataFrame({"sample": ["a", "b", "c"]}, index=["c1", "c2", "c3"])
    var = pd.DataFrame(index=[f"GENE_{i}" for i in range(4)])
    return anndata.AnnData(X=x, obs=obs, var=var)


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_load_ambient")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


# ═════════════════════════════════════════════════════════════════════
#  ambient safe_write path
# ═════════════════════════════════════════════════════════════════════


class TestAmbientWrite:
    @pytest.mark.parametrize("method", ["cellbender", "soupx"])
    def test_writes_ambient_h5ad_without_typeerror(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, method: str
    ) -> None:
        # make require_* succeed (avoids the installed-dependency branch)
        monkeypatch.setattr(f"core.utils._optional.require_{method}", lambda: None)
        cfg = _make_cfg(tmp_path, method)
        adata = _make_adata()

        # must NOT raise TypeError (the ambient safe_write really executes)
        _mod._run_ambient_correction(adata, cfg, _make_logger())

        out = tmp_path / "ambient_removed.h5ad"
        assert out.exists(), f"ambient_removed.h5ad was not written for {method}"
        reread = anndata.read_h5ad(out)
        assert reread.shape == adata.shape
        assert list(reread.obs_names) == list(adata.obs_names)


class TestAmbientMissingDependency:
    def test_skips_without_crash(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        def _missing() -> None:
            raise ImportError("cellbender is required ...")

        monkeypatch.setattr("core.utils._optional.require_cellbender", _missing)
        cfg = _make_cfg(tmp_path, "cellbender")

        _mod._run_ambient_correction(_make_adata(), cfg, _make_logger())

        assert not (tmp_path / "ambient_removed.h5ad").exists()
