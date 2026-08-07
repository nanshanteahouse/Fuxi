"""TDD RED-phase tests for full-gene SVG in ``spatial/steps/07_spatial_stats.py``.

Branch under test (plan ``.omo/plans/spatial-pipeline-rewrite-phase2.md``
todo 9): the Moran's I SVG screening MUST run on the **full gene set** drawn
from ``adata.raw`` (stored by 03_normalize BEFORE HVG subsetting), not on the
HVG-subset ``adata.X``.  Contract pinned here:

- ``run_spatial_autocorr`` feeds ``sq.gr.spatial_autocorr`` a matrix whose
  ``n_vars`` equals the full-gene count in ``adata.raw`` — NOT the HVG count
  in ``adata.X``.
- If ``adata.raw.X`` holds raw counts (max > 50), the matrix is
  ``normalize_total`` + ``log1p``-ed first so Moran's I runs on the same
  transform as the existing X path.
- Missing ``adata.raw`` → ``log.warning`` + fall back to ``adata.X``
  (graceful, no crash).
- A var present in ``adata.X`` but absent from ``adata.raw`` is dropped
  without crashing the SVG computation.
- ``svg_rankings.csv`` contains the full-gene ranking; ``uns['svg']`` records
  ``{'n_genes_tested', 'method': 'moran', 'n_top', 'moran_percentile'}``.
- Checkpoint-before-plot: ``06_svg.h5ad`` is written BEFORE the top-SVG
  spatial plot.
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from unittest.mock import MagicMock

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "spatial", "steps", "07_spatial_stats.py")
_spec = importlib.util.spec_from_file_location("spatial.steps._07_svg_full_gene_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

_N_FULL = 100
_N_HVG = 20


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════


def _make_raw_counts(n_spots: int, n_genes: int, seed: int, as_counts: bool) -> np.ndarray:
    rng = np.random.RandomState(seed)
    if as_counts:
        # guaranteed max > 50 so the counts-detection heuristic fires
        counts = rng.randint(0, 200, size=(n_spots, n_genes)).astype(np.float32) + 1
    else:
        counts = rng.poisson(lam=3.0, size=(n_spots, n_genes)).astype(np.float32)
    # guarantee every row has ≥1 count (normalize_total divides by row sums)
    counts = np.maximum(counts, 1.0)
    return counts


def _make_adata(
    n_spots: int = 30,
    n_full: int = _N_FULL,
    n_hvg: int = _N_HVG,
    seed: int = 7,
    raw_is_counts: bool = False,
) -> anndata.AnnData:
    """Synthetic AnnData: ``adata.raw`` = full gene set, ``adata.X`` = HVG subset."""
    full_counts = _make_raw_counts(n_spots, n_full, seed, raw_is_counts)

    raw = anndata.AnnData(X=sp.csr_matrix(full_counts))
    raw.var_names = [f"G{i}" for i in range(n_full)]
    raw.obs_names = [f"s{i}" for i in range(n_spots)]

    hvg_idx = np.arange(n_hvg)  # deterministic first genes → top-SVG overlap with X
    adata = anndata.AnnData(X=sp.csr_matrix(full_counts[:, hvg_idx].copy()))
    adata.var_names = [f"G{i}" for i in hvg_idx]
    adata.obs_names = [f"s{i}" for i in range(n_spots)]
    adata.raw = raw

    # spatial graph (keeps run_spatial_autocorr from rebuilding)
    conn = sp.eye(n_spots, format="csr")
    adata.obsp["spatial_connectivities"] = conn
    return adata


def _make_cfg(
    tmp_path,
    *,
    svg_n_top: int = 2000,
    moran_percentile: int = 90,
    run_autocorr: bool = True,
) -> MagicMock:
    cfg = MagicMock()
    cfg.spatial = MagicMock()
    cfg.spatial.run_autocorr = run_autocorr
    cfg.spatial.svg_n_top = svg_n_top
    cfg.spatial.moran_percentile = moran_percentile
    cfg.spatial.neighbors_n = 6
    cfg.spatial.neighbors_radius = 0
    cfg.execution = MagicMock()
    cfg.execution.n_jobs = 1
    cfg.normalization = MagicMock()
    cfg.normalization.normalize_target_sum = 10000
    cfg.table_dir = str(tmp_path)
    cfg.figure_dir = str(tmp_path)
    cfg.h5ad_dir = str(tmp_path)
    cfg.log_dir = str(tmp_path)
    cfg.plot = MagicMock()
    cfg.plot.umap_panel_size = (4, 4)
    cfg.plot.figure_dpi = 150
    cfg.plot.figure_format = "png"
    return cfg


def _fake_spatial_autocorr(captured: dict) -> callable:
    """Return a squidpy-like spatial_autocorr mock that records its input and
    writes ``uns['moranI']`` (index = whatever var_names it received)."""

    def _fake(adata, **kwargs):
        captured["n_vars"] = adata.n_vars
        captured["var_names"] = list(adata.var_names)
        captured["x_max"] = adata.X.max()
        rng = np.random.RandomState(0)
        n = adata.n_vars
        df = pd.DataFrame(
            {
                "I": rng.rand(n),
                "pval_norm": rng.rand(n),
                "pval_sim": rng.rand(n),
                "pval_bh": rng.rand(n),
                "pval_fdr_bh": rng.rand(n),
            },
            index=adata.var_names,
        )
        adata.uns["moranI"] = df

    return _fake


# ═════════════════════════════════════════════════════════════════════
#  Tests
# ═════════════════════════════════════════════════════════════════════


class TestFullGeneSVG:
    def test_spatial_autocorr_receives_full_gene_matrix(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """SVG input MUST be the full gene set (raw), NOT the HVG subset."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        cfg = _make_cfg(tmp_path)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))

        result = _mod.run_spatial_autocorr(adata, cfg, MagicMock())

        assert result is not None
        assert captured["n_vars"] == _N_FULL, (
            f"spatial_autocorr got {captured['n_vars']} genes; expected full gene set "
            f"{_N_FULL}, not the {_N_HVG}-gene HVG subset"
        )

    def test_raw_counts_log1p_normalized_before_moran(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """raw.X holding counts (max>50) is normalize_total+log1p'd before Moran's I."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG, raw_is_counts=True)
        cfg = _make_cfg(tmp_path)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))

        result = _mod.run_spatial_autocorr(adata, cfg, MagicMock())

        assert result is not None
        assert captured["n_vars"] == _N_FULL
        assert captured["x_max"] < 50, (
            f"Moran's I received max expression {captured['x_max']} — raw counts were "
            "NOT log1p-normalized before spatial autocorrelation"
        )

    def test_falls_back_to_x_with_warning_when_raw_missing(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """No adata.raw → warn + use adata.X (HVG) without crashing."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        del adata.raw
        cfg = _make_cfg(tmp_path)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))
        log = MagicMock()

        result = _mod.run_spatial_autocorr(adata, cfg, log)

        assert result is not None
        assert captured["n_vars"] == _N_HVG, "fallback must use adata.X (HVG subset)"
        log.warning.assert_called()

    def test_gene_missing_in_x_handled_gracefully(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A var present in X but absent from raw must not crash SVG computation."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        raw = adata.raw.to_adata()  # Raw → AnnData so rebuilt.raw is settable
        # inject an extra gene into X that raw does not contain
        ghost = sp.csr_matrix(np.ones((adata.n_obs, 1), dtype=np.float32))
        rebuilt = anndata.AnnData(X=sp.hstack([adata.X, ghost]))
        rebuilt.var_names = list(adata.var_names) + ["GHOST_GENE"]
        rebuilt.obs_names = adata.obs_names
        rebuilt.raw = raw
        rebuilt.obsp["spatial_connectivities"] = adata.obsp["spatial_connectivities"]
        cfg = _make_cfg(tmp_path)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))

        result = _mod.run_spatial_autocorr(rebuilt, cfg, MagicMock())

        assert result is not None
        # GHOST_GENE is dropped — SVG runs on the raw gene set only
        assert captured["n_vars"] == _N_FULL
        assert "GHOST_GENE" not in captured["var_names"]

    def test_svg_rankings_csv_contains_full_gene_ranking(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """svg_rankings.csv must contain the full-gene ranking (not just HVG)."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        cfg = _make_cfg(tmp_path, svg_n_top=2000)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))

        _mod.run_spatial_autocorr(adata, cfg, MagicMock())

        csv_path = os.path.join(str(tmp_path), "svg_rankings.csv")
        assert os.path.exists(csv_path), f"svg_rankings.csv not written at {csv_path}"
        df = pd.read_csv(csv_path, index_col=0)
        assert len(df) == _N_FULL, f"csv has {len(df)} rows; expected full-gene {_N_FULL}"
        assert len(df) > _N_HVG, "csv must NOT be limited to the HVG subset"

    def test_uns_svg_metadata(self, tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
        """uns['svg'] records gene count, method, and top-N."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        cfg = _make_cfg(tmp_path)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))

        _mod.run_spatial_autocorr(adata, cfg, MagicMock())

        meta = adata.uns["svg"]
        assert meta["n_genes_tested"] == _N_FULL
        assert meta["method"] == "moran"
        assert meta["n_top"] == min(cfg.spatial.svg_n_top, _N_FULL)

    def test_moran_percentile_default_preserved(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """moran_percentile (default 90) survives into uns['svg'] metadata."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        cfg = _make_cfg(tmp_path, moran_percentile=90)
        captured: dict = {}
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_spatial_autocorr(captured))

        _mod.run_spatial_autocorr(adata, cfg, MagicMock())

        assert adata.uns["svg"]["moran_percentile"] == 90


class TestMainFlow:
    def _patched_main(
        self,
        adata: anndata.AnnData,
        cfg: MagicMock,
        tmp_path,
        monkeypatch: pytest.MonkeyPatch,
        call_log: list,
    ) -> None:
        (tmp_path / "05_annotated.h5ad").touch()

        def _fake_autocorr(adata_obj, **kwargs):
            rng = np.random.RandomState(0)
            n = adata_obj.n_vars
            df = pd.DataFrame(
                {"I": rng.rand(n), "pval_sim": rng.rand(n), "pval_bh": rng.rand(n)},
                index=adata_obj.var_names,
            )
            adata_obj.uns["moranI"] = df

        monkeypatch.setattr(
            _mod.argparse.ArgumentParser,
            "parse_args",
            lambda self: argparse.Namespace(config="/tmp/test.yaml"),
        )
        monkeypatch.setattr(_mod, "resolve_config", lambda *a, **k: cfg)
        monkeypatch.setattr(_mod, "setup_logger", lambda *a, **k: MagicMock())
        monkeypatch.setattr(_mod.sc, "read", lambda *a, **k: adata)
        monkeypatch.setattr(_mod.sq.gr, "spatial_autocorr", _fake_autocorr)
        monkeypatch.setattr(_mod, "run_de_per_cluster", lambda *a, **k: None)
        monkeypatch.setattr(_mod, "run_nhood_enrichment", lambda *a, **k: None)
        monkeypatch.setattr(_mod, "run_co_occurrence", lambda *a, **k: None)
        monkeypatch.setattr(_mod, "run_interaction_matrix", lambda *a, **k: None)

        def _record_save(obj, *a, **kw):
            call_log.append("save")

        def _record_plot(*a, **kw):
            call_log.append("plot")
            return None

        monkeypatch.setattr(_mod, "safe_write", _record_save)
        monkeypatch.setattr(_mod, "save_figure", _record_plot)
        monkeypatch.setattr(_mod.sc.pl, "umap", _record_plot)
        monkeypatch.setattr(_mod.sq.pl, "spatial_scatter", _record_plot)
        _mod.main()

    def test_main_checkpoint_before_plot(self, tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
        """06_svg.h5ad checkpoint MUST be written before the top-SVG plot."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        adata.obsm["spatial"] = np.random.RandomState(1).randn(adata.n_obs, 2)
        cfg = _make_cfg(tmp_path)
        call_log: list = []

        self._patched_main(adata, cfg, tmp_path, monkeypatch, call_log)

        assert "save" in call_log, "06_svg.h5ad checkpoint was never written"
        assert "plot" in call_log, "expected the top-SVG spatial plot to be reached"
        assert call_log.index("save") < call_log.index("plot"), (
            "checkpoint must be written BEFORE the top-SVG spatial plot "
            "(checkpoint-before-plot hard convention)"
        )

    def test_main_marks_spatially_variable_from_full_gene_ranking(
        self, tmp_path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """var['spatially_variable'] derives from the full-gene SVG ranking."""
        adata = _make_adata(n_full=_N_FULL, n_hvg=_N_HVG)
        adata.obsm["spatial"] = np.random.RandomState(1).randn(adata.n_obs, 2)
        cfg = _make_cfg(tmp_path, svg_n_top=5)
        call_log: list = []

        self._patched_main(adata, cfg, tmp_path, monkeypatch, call_log)

        assert "spatially_variable" in adata.var
        assert adata.uns["svg"]["n_top"] == 5
        # full-gene ranking has 100 genes; only those present in the HVG X are marked
        assert int(adata.var["spatially_variable"].sum()) <= min(5, _N_HVG)
