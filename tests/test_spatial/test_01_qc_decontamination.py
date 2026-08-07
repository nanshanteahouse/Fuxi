"""TDD red-phase tests for the optional ambient-RNA decontamination stage in
``spatial/steps/01_qc.py``.

Branch under test (plan ``.omo/plans/spatial-pipeline-rewrite-phase2.md``
todo 5): a pre-QC, opt-in decontamination hook keyed off
``cfg.spatial.decontamination`` (``none`` | ``cellbender`` | ``soupx``).  The
schema field lands in a later wave, so the implementation reads it via
``getattr(cfg.spatial, "decontamination", "none")`` and MUST leave the existing
QC flow byte-identical when the value is ``"none"`` (the default).

Contract pinned here
---------------------
- ``run_decontamination(adata, cfg, log)`` is a no-op for ``"none"`` — the
  CellBender / SoupX runners are never invoked and ``uns`` is untouched.
- ``_decontamination_method(cfg)`` is backward compatible: a config without the
  ``decontamination`` attribute (or without ``cfg.spatial``) resolves to
  ``"none"``.
- ``cellbender`` branch: ``require_cellbender()`` gates the real tool.  On
  ``ImportError`` → ``log.warning`` + graceful skip (no crash, counts
  untouched, ``uns['decontamination'] == {status: 'skipped', ...}``).  When the
  dependency is present the real ``_run_cellbender_on_counts`` is invoked and
  the returned decontaminated matrix replaces ``X`` (original counts preserved
  in ``layers['raw']``).
- ``soupx`` branch: identical gate semantics via ``require_soupx()``.
- Unknown methods degrade to a ``log.warning`` and continue (no crash).

The real CellBender/SoupX executions are NOT exercised (CellBender GPU training
is far too slow for unit tests); the real invocation helpers
``_run_cellbender_on_counts`` / ``_run_soupx_on_counts`` are monkeypatched to
return a known decontaminated matrix.  Missing-dependency gating is verified
against the actual ``core.utils._optional`` guards.
"""

from __future__ import annotations

import importlib
import importlib.util
import logging
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import anndata
import numpy as np
import pytest
from scipy.sparse import csr_matrix

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "spatial", "steps", "01_qc.py")
_spec = importlib.util.spec_from_file_location(
    "spatial.steps._01_qc_decontamination_test", _STEP_PATH
)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════


def _make_cfg(tmp_path: Path, method: str) -> MagicMock:
    """Minimal config: only ``spatial.decontamination`` is consulted."""
    cfg = MagicMock()
    cfg.spatial = MagicMock()
    cfg.spatial.decontamination = method
    return cfg


def _make_adata(seed: int = 7) -> anndata.AnnData:
    rng = np.random.RandomState(seed)
    x = csr_matrix(rng.randint(0, 5, size=(6, 8)).astype(np.float32))
    adata = anndata.AnnData(X=x)
    adata.obs_names = [f"s{i}" for i in range(6)]
    adata.var_names = [f"G{i}" for i in range(8)]
    return adata


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_01_qc_decontamination")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


def _missing_dep(pkg: str):
    """Factory for a ``require_*`` guard that raises like the real one."""

    def _raise() -> None:
        raise ImportError(
            f"{pkg} is required for ambient RNA removal. Install with: pip install fuxi[{pkg}]"
        )

    return _raise


# ═════════════════════════════════════════════════════════════════════
#  Gate: "none" (default) must keep the existing QC flow untouched
# ═════════════════════════════════════════════════════════════════════


class TestNoneGate:
    def test_none_is_noop_and_never_calls_runners(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        cfg = _make_cfg(tmp_path, "none")
        adata = _make_adata()
        cellbender = MagicMock()
        soupx = MagicMock()
        monkeypatch.setattr(_mod, "_decontaminate_cellbender", cellbender)
        monkeypatch.setattr(_mod, "_decontaminate_soupx", soupx)

        _mod.run_decontamination(adata, cfg, _make_logger())

        cellbender.assert_not_called()
        soupx.assert_not_called()
        assert "decontamination" not in adata.uns
        # X and layers untouched — existing QC flow sees the raw object
        assert "raw" not in adata.layers
        assert adata.n_obs == 6 and adata.n_vars == 8

    def test_method_resolves_none_when_attribute_missing(self) -> None:
        # schema field lands in a later wave — old configs must resolve "none"
        cfg = MagicMock()
        cfg.spatial = SimpleNamespace()  # no `decontamination` attribute yet
        assert _mod._decontamination_method(cfg) == "none"

    def test_method_resolves_none_when_spatial_missing(self) -> None:
        cfg = SimpleNamespace()  # no cfg.spatial at all (defensive)
        assert _mod._decontamination_method(cfg) == "none"

    def test_method_respects_explicit_value(self) -> None:
        cfg = MagicMock()
        cfg.spatial = SimpleNamespace(decontamination="cellbender")
        assert _mod._decontamination_method(cfg) == "cellbender"


# ═════════════════════════════════════════════════════════════════════
#  CellBender branch
# ═════════════════════════════════════════════════════════════════════


class TestCellbender:
    def test_missing_dependency_skips_with_warning(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("core.utils._optional.require_cellbender", _missing_dep("cellbender"))
        adata = _make_adata()
        x_before = adata.X.toarray().copy()
        log = MagicMock()
        cfg = _make_cfg(tmp_path, "cellbender")

        _mod.run_decontamination(adata, cfg, log)

        log.warning.assert_called_once()
        meta = adata.uns["decontamination"]
        assert meta["status"] == "skipped"
        assert meta["method"] == "cellbender"
        assert "cellbender is required" in meta["reason"]
        # counts untouched — graceful degradation, no crash
        assert np.array_equal(adata.X.toarray(), x_before)
        assert "raw" not in adata.layers

    def test_calls_cellbender_and_replaces_counts(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("core.utils._optional.require_cellbender", lambda: None)
        adata = _make_adata()
        clean = csr_matrix(adata.X.toarray() * 0.5).astype(np.float32)
        monkeypatch.setattr(_mod, "_run_cellbender_on_counts", lambda adata, cfg, log: clean)
        log = _make_logger()

        _mod.run_decontamination(adata, _make_cfg(tmp_path, "cellbender"), log)

        meta = adata.uns["decontamination"]
        assert meta == {"status": "completed", "method": "cellbender"}
        assert "raw" in adata.layers  # original counts preserved
        np.testing.assert_allclose(adata.layers["raw"].toarray(), _make_adata().X.toarray())
        np.testing.assert_allclose(adata.X.toarray(), clean.toarray())

    def test_failure_keeps_original_counts(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("core.utils._optional.require_cellbender", lambda: None)
        adata = _make_adata()
        x_before = adata.X.toarray().copy()

        def _boom(adata, cfg, log):
            raise RuntimeError("cellbender crashed")

        monkeypatch.setattr(_mod, "_run_cellbender_on_counts", _boom)

        _mod.run_decontamination(adata, _make_cfg(tmp_path, "cellbender"), _make_logger())

        assert adata.uns["decontamination"]["status"] == "failed"
        assert np.array_equal(adata.X.toarray(), x_before)
        assert "raw" not in adata.layers


# ═════════════════════════════════════════════════════════════════════
#  SoupX branch
# ═════════════════════════════════════════════════════════════════════


class TestSoupx:
    def test_missing_dependency_skips_with_warning(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("core.utils._optional.require_soupx", _missing_dep("soupx"))
        adata = _make_adata()
        x_before = adata.X.toarray().copy()
        log = MagicMock()

        _mod.run_decontamination(adata, _make_cfg(tmp_path, "soupx"), log)

        log.warning.assert_called_once()
        meta = adata.uns["decontamination"]
        assert meta["status"] == "skipped"
        assert meta["method"] == "soupx"
        assert "soupx is required" in meta["reason"]
        assert np.array_equal(adata.X.toarray(), x_before)

    def test_calls_soupx_and_replaces_counts(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("core.utils._optional.require_soupx", lambda: None)
        adata = _make_adata()
        clean = csr_matrix(adata.X.toarray() * 0.8).astype(np.float32)
        monkeypatch.setattr(_mod, "_run_soupx_on_counts", lambda adata, cfg, log: clean)
        log = _make_logger()

        _mod.run_decontamination(adata, _make_cfg(tmp_path, "soupx"), log)

        meta = adata.uns["decontamination"]
        assert meta == {"status": "completed", "method": "soupx"}
        assert "raw" in adata.layers
        np.testing.assert_allclose(adata.X.toarray(), clean.toarray())


# ═════════════════════════════════════════════════════════════════════
#  Unknown method + dispatcher robustness
# ═════════════════════════════════════════════════════════════════════


class TestDispatcher:
    def test_unknown_method_warns_and_continues(self, tmp_path: Path) -> None:
        adata = _make_adata()
        log = MagicMock()

        _mod.run_decontamination(adata, _make_cfg(tmp_path, "bogus"), log)

        log.warning.assert_called_once()
        assert "decontamination" not in adata.uns

    def test_cellbender_dispatched_through_run_decontamination(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("core.utils._optional.require_cellbender", lambda: None)
        adata = _make_adata()
        clean = csr_matrix(adata.X.toarray() * 0.25).astype(np.float32)
        monkeypatch.setattr(_mod, "_run_cellbender_on_counts", lambda adata, cfg, log: clean)

        _mod.run_decontamination(adata, _make_cfg(tmp_path, "cellbender"), _make_logger())

        assert adata.uns["decontamination"]["status"] == "completed"
        assert adata.uns["decontamination"]["method"] == "cellbender"
