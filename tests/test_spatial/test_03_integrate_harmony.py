"""Tests for spatial/steps/03_normalize.py — multi-slide Harmony integration.

TDD RED phase (plan ``.omo/plans/spatial-pipeline-rewrite-phase2.md`` todo 6).

The spatial Step 03 inserts an optional batch-integration phase between PCA
and the checkpoint save, mirroring the RNA 03_integrate conventions:

- Only runs when ``cfg.integration.method == "harmony"`` (read via
  ``getattr(cfg.integration, "method", "none")`` for backward compatibility —
  spatial default is "none", i.e. no integration).
- Only runs on multi-slide data (obs ``batch_key`` has >= 2 unique values);
  single-slide data skips with a log.info.
- GPU-first via ``core.utils.gpu_harmony``; CPU fallback via
  ``harmonypy.run_harmony``. On the CPU path the run report is scanned for
  "perfectly collinear" warnings → the corrected embedding is NOT applied and
  ``uns['harmony_skipped'] = {'reason': 'collinearity', ...}`` is recorded
  (cross-batch critical fix T3 — notes/engineering/2026-07-15_cross_batch_critical_fixes.md).
- n_pcs min() guard: downstream PCA/neighbor consumers must never slice
  beyond the actual ``X_integrated`` width
  (notes/engineering/2026-07-28_n_pcs_use_scvi_dimension_mismatch.md).
- checkpoint-before-plot hard convention: ``safe_write`` of ``03_processed.h5ad``
  MUST precede any comparison plotting
  (notes/engineering/2026-07-30_checkpoint_before_plot.md).
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

# ── Ensure repo root is on sys.path ──────────────────────────────────────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 03_normalize module via file path ───────────────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "spatial", "steps", "03_normalize.py")
_spec = importlib.util.spec_from_file_location(
    "spatial.steps._03_normalize_test",
    _STEP_PATH,
)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════


def _make_adata(
    n_cells: int = 30,
    n_genes: int = 100,
    n_unique_samples: int = 2,
    seed: int = 42,
) -> AnnData:
    """Create a minimal multi-slide AnnData with X_pca pre-computed.

    obs['sample'] has ``n_unique_samples`` unique values (2 = multi-slide,
    1 = single-slide) — the column added by spatial 00_load for merged slides.
    """
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    if n_unique_samples >= 2:
        adata.obs["sample"] = rng.choice(["S1", "S2"], n_cells)
    else:
        adata.obs["sample"] = ["S1"] * n_cells
    adata.obsm["X_pca"] = rng.randn(n_cells, 50)
    adata.var["highly_variable"] = [True] * n_genes
    adata.uns["pca"] = {"variance_ratio": np.zeros(n_genes)}
    return adata


def _make_cfg(
    method: str = "harmony",
    h5ad_dir: str | None = None,
) -> MagicMock:
    """Create a Config mock sufficient for the integration + save section."""
    cfg = MagicMock()

    # Integration settings (RNA-style enum; spatial default "none")
    cfg.integration = MagicMock()
    cfg.integration.method = method
    cfg.integration.batch_key = "sample"
    cfg.integration.max_iter = 20

    # PCA settings
    cfg.pca = MagicMock()
    cfg.pca.n_pcs_use = 50

    # Execution settings (device="cpu" forces the harmonypy CPU path in tests)
    cfg.execution = MagicMock()
    cfg.execution.device = "cpu"
    cfg.execution.random_seed = 42

    # HVG settings
    cfg.hvg = MagicMock()
    cfg.hvg.n_top_genes = 100
    cfg.hvg.flavor = "seurat_v3"
    cfg.hvg.batch_key = None
    cfg.hvg.forced_genes = []

    # Marker settings (force-include block reads these)
    cfg.marker = MagicMock()
    cfg.marker.marker_dict = {}

    # Normalization settings
    cfg.normalization = MagicMock()
    cfg.normalization.normalize_target_sum = 10000

    # Spatial neighbor settings
    cfg.spatial = MagicMock()
    cfg.spatial.neighbors_n = 6
    cfg.spatial.neighbors_radius = 0

    # I/O paths
    cfg.figure_dir = "/tmp"
    cfg.h5ad_dir = h5ad_dir or "/tmp"
    cfg.log_dir = "/tmp"

    # Plot config (comparison plot reads these)
    cfg.plot = MagicMock()
    cfg.plot.figure_dpi = 150
    cfg.plot.figure_format = "png"

    return cfg


def _make_report(
    warnings: list[str] | None = None,
) -> MagicMock:
    """Create a batch-diagnosis-style report mock with optional warnings."""
    report = MagicMock()
    report.batch_cols = ["sample"]
    report.biology_cols = []
    report.ambiguous_cols = []
    report.warnings = warnings or []
    return report


def _make_harmony_result(
    n_cells: int = 30,
    n_dims: int = 50,
    warnings: list[str] | None = None,
) -> MagicMock:
    """Create a harmonypy run result mock (Z_corr + optional collinear report)."""
    ho = MagicMock()
    ho.Z_corr = np.random.RandomState(42).randn(n_cells, n_dims)
    ho.report = _make_report(warnings=warnings)
    return ho


def _fake_spatial_neighbors(adata: AnnData, *args, **kwargs) -> None:
    """Side-effect for the patched sq.gr.spatial_neighbors: populate obsp.

    The step verifies ``spatial_connectivities`` after building the graph, so
    the mock must actually create it.
    """
    eye = sp.eye(adata.n_obs, format="csr")
    adata.obsp["spatial_connectivities"] = eye
    adata.obsp["spatial_distances"] = eye


def _ensure_input(tmp_path) -> None:
    """Create the Step 02/01 input file so input resolution succeeds.

    sc.read is patched, so the empty file is sufficient.
    """
    (tmp_path / "02_image.h5ad").touch()


def _capture_adata_on_save(captured: list) -> callable:
    """Return a safe_write side-effect that captures the final adata."""

    def _side_effect(adata, *args, **kwargs):
        captured.append(adata)

    return _side_effect


def _common_patches(adata: AnnData, cfg: MagicMock):
    """Return the shared patch chain (I/O + scanpy/squidpy ops)."""
    return [
        patch.object(
            _mod.argparse.ArgumentParser,
            "parse_args",
            return_value=argparse.Namespace(config="/tmp/test.yaml"),
        ),
        patch.object(_mod, "resolve_config", return_value=cfg),
        patch.object(_mod, "setup_logger", return_value=MagicMock()),
        patch.object(_mod.sc, "read", return_value=adata),
        patch.object(_mod.sc.pp, "highly_variable_genes"),
        patch.object(_mod.sc.pp, "normalize_total"),
        patch.object(_mod.sc.pp, "log1p"),
        patch.object(_mod.sc.pp, "pca"),
        patch.object(
            _mod.sq.gr,
            "spatial_neighbors",
            side_effect=_fake_spatial_neighbors,
        ),
    ]


def _run_main(
    adata: AnnData,
    cfg: MagicMock,
    extra_patches: list,
) -> None:
    """Enter all patch contexts via ExitStack, then run main().

    ``*``-unpacking inside a parenthesized ``with`` is not supported, so the
    shared patch chain plus the per-test extra patches are stacked manually.
    """
    from contextlib import ExitStack

    with ExitStack() as stack:
        for p in _common_patches(adata, cfg) + list(extra_patches):
            stack.enter_context(p)
        _mod.main()


# ═══════════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════════


def test_collinearity_guard_skips_harmony_and_sets_flag(tmp_path) -> None:
    """Perfectly collinear batch_key → integration skipped + uns flag set.

    Given:  multi-slide data, method="harmony", device="cpu".
    When:   run_harmony reports "perfectly collinear".
    Then:   X_integrated NOT created, uns['harmony_skipped'] set,
            GPU gpu_harmony NOT invoked.
    """
    adata = _make_adata()
    cfg = _make_cfg(method="harmony", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    collinear_ho = _make_harmony_result(
        warnings=[
            "Column 'sample' is perfectly collinear with 'genotype' (V=1.0) — redundant column.",
        ]
    )

    captured: list[AnnData] = []
    gpu_mock = MagicMock()
    _run_main(
        adata,
        cfg,
        [
            patch("harmonypy.run_harmony", return_value=collinear_ho),
            patch.object(_mod, "gpu_harmony", gpu_mock),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
        ],
    )

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]

    # Guard must have fired
    assert "harmony_skipped" in result.uns, "harmony_skipped should be set when guard fires"
    assert result.uns["harmony_skipped"]["reason"] == "collinearity"
    assert len(result.uns["harmony_skipped"]["warnings"]) > 0

    # Corrected embedding must NOT be applied when the guard fires
    assert "X_integrated" not in result.obsm, (
        "X_integrated should NOT be created when the collinearity guard aborts"
    )
    gpu_mock.assert_not_called()


def test_harmony_success_creates_x_integrated_with_batch_key(tmp_path) -> None:
    """Clean run → X_integrated created, harmonize called with the batch_key.

    Given:  multi-slide data, method="harmony", no collinearity warnings.
    When:   main() runs through the integration phase.
    Then:   obsm['X_integrated'] present; run_harmony called with
            vars_use=["sample"].
    """
    adata = _make_adata()
    cfg = _make_cfg(method="harmony", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    ok_ho = _make_harmony_result(n_dims=50)

    captured: list[AnnData] = []
    run_mock = MagicMock(return_value=ok_ho)
    _run_main(
        adata,
        cfg,
        [
            patch("harmonypy.run_harmony", run_mock),
            patch.object(_mod.sc.pl, "embedding", return_value=None),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
        ],
    )

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]

    # Harmony must have run with the batch_key
    run_mock.assert_called_once()
    assert run_mock.call_args.kwargs["vars_use"] == ["sample"], (
        "run_harmony must be called with vars_use=[batch_key]"
    )
    assert "X_integrated" in result.obsm, "X_integrated should be created on success"
    assert result.obsm["X_integrated"].shape == (adata.n_obs, 50)
    assert "harmony_skipped" not in result.uns, (
        "harmony_skipped should NOT be set when Harmony succeeds"
    )


def test_n_pcs_min_guard_caps_below_n_pcs_use(tmp_path) -> None:
    """X_integrated narrower than n_pcs_use → min() guard keeps downstream safe.

    Given:  multi-slide data, method="harmony", fake output has only 30 dims
            while cfg.pca.n_pcs_use=50.
    When:   main() runs.
    Then:   No crash; X_integrated keeps its 30 dims; _effective_n_pcs caps
            downstream PCA/neighbor n_pcs at min(50, 30) = 30.
    """
    adata = _make_adata()
    cfg = _make_cfg(method="harmony", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    narrow_ho = _make_harmony_result(n_dims=30)

    captured: list[AnnData] = []
    _run_main(
        adata,
        cfg,
        [
            patch("harmonypy.run_harmony", return_value=narrow_ho),
            patch.object(_mod.sc.pl, "embedding", return_value=None),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
        ],
    )

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]

    # The narrow embedding is kept as-is; the guard must cap downstream n_pcs
    assert result.obsm["X_integrated"].shape == (adata.n_obs, 30)
    assert _mod._effective_n_pcs(cfg.pca.n_pcs_use, result.obsm["X_integrated"].shape[1]) == 30
    # Guard is a no-op when the embedding is wide enough
    assert _mod._effective_n_pcs(cfg.pca.n_pcs_use, 50) == 50


def test_single_slide_skips_integration(tmp_path) -> None:
    """obs batch_key has 1 unique value → no integration, no X_integrated.

    Given:  single-slide data (obs['sample'] all "S1"), method="harmony".
    When:   main() runs.
    Then:   X_integrated NOT created; run_harmony NOT called; X_pca stays primary.
    """
    adata = _make_adata(n_unique_samples=1)
    cfg = _make_cfg(method="harmony", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    run_mock = MagicMock()
    _run_main(
        adata,
        cfg,
        [
            patch("harmonypy.run_harmony", run_mock),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
        ],
    )

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]
    assert "X_integrated" not in result.obsm, "single-slide data must skip Harmony integration"
    assert "X_pca" in result.obsm, "X_pca remains the primary representation"
    run_mock.assert_not_called()


def test_method_not_harmony_skips_integration(tmp_path) -> None:
    """integration.method == "none" (spatial default) → no harmonize call.

    Given:  multi-slide data, method="none".
    When:   main() runs.
    Then:   X_integrated NOT created; run_harmony NOT called; X_pca stays primary.
    """
    adata = _make_adata()
    cfg = _make_cfg(method="none", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    run_mock = MagicMock()
    _run_main(
        adata,
        cfg,
        [
            patch("harmonypy.run_harmony", run_mock),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
        ],
    )

    assert len(captured) == 1, "safe_write should have been called once"
    result = captured[0]
    assert "X_integrated" not in result.obsm, "method != 'harmony' must not create X_integrated"
    assert "X_pca" in result.obsm
    run_mock.assert_not_called()


def test_checkpoint_before_plot(tmp_path) -> None:
    """safe_write of 03_processed.h5ad MUST precede the comparison plot.

    Given:  multi-slide data, method="harmony", successful run.
    When:   main() runs.
    Then:   call order is save → plot → plot (checkpoint first).
    """
    adata = _make_adata()
    cfg = _make_cfg(method="harmony", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    ok_ho = _make_harmony_result(n_dims=50)

    call_log: list[str] = []

    def _record_save(adata, *args, **kwargs):
        call_log.append("save")

    def _record_plot(*args, **kwargs):
        call_log.append("plot")
        return None

    _run_main(
        adata,
        cfg,
        [
            patch("harmonypy.run_harmony", return_value=ok_ho),
            patch.object(_mod.sc.pl, "embedding", side_effect=_record_plot),
            patch.object(_mod, "safe_write", side_effect=_record_save),
        ],
    )

    assert call_log.count("save") == 1, "safe_write should be called exactly once"
    assert call_log.count("plot") == 2, "comparison plot draws before and after views"
    assert call_log.index("save") < call_log.index("plot"), (
        "checkpoint must be written BEFORE the comparison plot"
    )
