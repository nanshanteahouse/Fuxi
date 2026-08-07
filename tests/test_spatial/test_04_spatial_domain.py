"""Tests for spatial/steps/04_cluster.py — spatial-domain clustering.

TDD RED phase (plan ``.omo/plans/spatial-pipeline-rewrite-phase2.md`` todo 7).

Round-1 multi-metric grid clustering (n_neighbors × resolution + UMAP sweep)
is untouched; this wave ADDS an optional spatial-domain label driven by
``CFG.spatial.domain_method``:

- ``"none"`` (off) — existing behavior unchanged, no spatial-domain column.
- ``"leiden_spatial"`` (default, zero extra deps) — Leiden on the
  ``spatial_connectivities`` graph at ``domain_resolution``.
- ``"stagate"`` (optional GPU deep model) — gated import; on ImportError it
  logs a warning and falls back to ``leiden_spatial`` (never crashes).

Spatial domains are a *supplementary* label: ``obs['spatial_domain']`` plus
``uns['spatial_domain']`` metadata. The round-1 ``obs['leiden']`` cluster
column is preserved.

Hard conventions:
- checkpoint-before-plot: ``safe_write`` of the domain-bearing checkpoint MUST
  precede any spatial-domain plotting
  (notes/engineering/2026-07-30_checkpoint_before_plot.md).
- n_pcs min() guard: the STAGATE embedding may be narrower than
  ``cfg.pca.n_pcs_use``; downstream neighbors n_pcs is capped
  (notes/engineering/2026-07-28_n_pcs_use_scvi_dimension_mismatch.md).
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

# ── Load the 04_cluster module via file path ─────────────────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "spatial", "steps", "04_cluster.py")
_spec = importlib.util.spec_from_file_location(
    "spatial.steps._04_cluster_test",
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
    *,
    with_connectivities: bool = True,
    seed: int = 42,
) -> AnnData:
    """Minimal spot-level AnnData as Step 03 would hand over.

    Includes a pre-populated ``obs['leiden']`` column simulating the round-1
    multi-metric auto-lock — the tests assert it survives the domain phase.
    """
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obsm["X_pca"] = rng.randn(n_cells, 50)
    # Round-1 multi-metric result (must never be overwritten by the domain phase)
    adata.obs["leiden"] = [str(i % 4) for i in range(n_cells)]
    if with_connectivities:
        eye = sp.eye(n_cells, format="csr")
        adata.obsp["spatial_connectivities"] = eye
        adata.obsp["spatial_distances"] = eye
    return adata


def _make_cfg(
    domain_method: str = "leiden_spatial",
    domain_resolution: float = 1.0,
    h5ad_dir: str | None = None,
) -> MagicMock:
    """Config mock sufficient for the full main() flow (grid + domain phase)."""
    cfg = MagicMock()

    # Round-1 grid clustering settings
    cfg.clustering = MagicMock()
    cfg.clustering.param_grid_n_neighbors = [15]
    cfg.clustering.param_grid_resolutions = [0.5]
    cfg.clustering.multi_metric_adaptive_resolution = False
    cfg.clustering.umap_min_dist = 0.3
    cfg.clustering.umap_spread = 1.0
    cfg.clustering.leiden_flavor = "igraph"
    cfg.clustering.leiden_n_iterations = 2
    cfg.clustering.cluster_selection_method = "multi_metric"
    cfg.clustering.best_resolution = 1.0
    cfg.clustering.best_n_neighbors = 0
    cfg.clustering.umap_paga_init = False

    # Spatial-domain settings (schema fields arrive in Wave 3; getattr default)
    cfg.spatial = MagicMock()
    cfg.spatial.domain_method = domain_method
    cfg.spatial.domain_resolution = domain_resolution

    cfg.pca = MagicMock()
    cfg.pca.n_pcs_use = 50

    cfg.execution = MagicMock()
    cfg.execution.device = "cpu"
    cfg.execution.random_seed = 42
    cfg.execution.n_jobs = 4

    cfg.plot = MagicMock()
    cfg.plot.figure_dpi = 150
    cfg.plot.figure_format = "png"
    cfg.plot.umap_panel_size = (4, 4)

    base = h5ad_dir or "/tmp"
    cfg.figure_dir = base
    cfg.h5ad_dir = base
    cfg.log_dir = base
    cfg.table_dir = base
    return cfg


def _make_results_summary() -> list[dict]:
    """Single grid entry whose auto-lock keys are absent from the adata.

    Step 04 then skips its round-1 cluster checkpoint write (leiden_15_0.5 /
    umap_15_0.5 not materialized), so the only clustering-phase checkpoint is
    the spatial-domain one — keeping the save/plot call-order assertions clean.
    """
    return [
        {
            "n_neighbors": 15,
            "resolution": 0.5,
            "cluster_key": "leiden_15_0.5",
            "score": 0.42,
        }
    ]


def _fake_leiden(
    adata: AnnData,
    *args,
    resolution: float = 1.0,
    key_added: str = "leiden",
    **kwargs,
) -> None:
    """Side-effect for the patched sc.tl.leiden: materialize the label column."""
    adata.obs[key_added] = [str(i % 3) for i in range(adata.n_obs)]


def _capture_adata_on_save(captured: list) -> callable:
    """Return a safe_write side-effect that captures the written adata."""

    def _side_effect(adata, *args, **kwargs):
        captured.append(adata)

    return _side_effect


def _common_patches(adata: AnnData, cfg: MagicMock):
    """Round-1 flow is mocked out wholesale; only the domain phase is real."""
    return [
        patch.object(
            _mod.argparse.ArgumentParser,
            "parse_args",
            return_value=argparse.Namespace(config="/tmp/test.yaml"),
        ),
        patch.object(_mod, "resolve_config", return_value=cfg),
        patch.object(_mod, "setup_logger", return_value=MagicMock()),
        patch.object(_mod.sc, "read", return_value=adata),
        # Round-1 grid clustering + multi-metric selection (fully mocked)
        patch.object(_mod, "grid_search_clustering", return_value=_make_results_summary()),
        patch("core.cluster.evaluation._detect_granularity", return_value="tissue"),
        patch("core.cluster.evaluation.enrich_grid_results"),
        patch.object(_mod, "select_best_params", return_value=(15, 0.5, "multi_metric", "test")),
        patch.object(_mod, "select_best_umap_params", return_value=(0.3, 1.0, "convex_hull", [])),
        # Plotting / I/O plumbing
        patch.object(_mod.sc.tl, "umap"),
        patch.object(_mod.sc.pl, "umap"),
        patch.object(_mod, "safe_plot"),
        patch.object(_mod.plt, "subplots", return_value=(MagicMock(), [[MagicMock()]])),
        patch.object(_mod.plt, "close"),
        patch.object(_mod, "save_figure"),
    ]


def _ensure_input(tmp_path) -> None:
    """Create the Step 03 input file so the input-resolution check succeeds.

    sc.read is patched, so the empty file is sufficient."""
    (tmp_path / "03_processed.h5ad").touch()


def _run_main(
    adata: AnnData,
    cfg: MagicMock,
    extra_patches: list,
) -> None:
    """Enter all patch contexts via ExitStack, then run main().

    ``*``-unpacking inside a parenthesized ``with`` is not supported, so the
    shared patch chain plus the per-test extra patches are stacked manually.
    Later patches override earlier ones on the same target.
    """
    from contextlib import ExitStack

    with ExitStack() as stack:
        for p in _common_patches(adata, cfg) + list(extra_patches):
            stack.enter_context(p)
        _mod.main()


# ═══════════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════════


def test_domain_none_preserves_existing_behavior(tmp_path) -> None:
    """domain_method='none' → no spatial domain, no Leiden on the spatial graph.

    Round-1 multi-metric behavior stays untouched: obs['leiden'] intact,
    no 'spatial_domain' column, no uns metadata, no spatial-graph leiden call.
    """
    adata = _make_adata()
    cfg = _make_cfg(domain_method="none", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    leiden_mock = MagicMock(side_effect=_fake_leiden)

    _run_main(
        adata,
        cfg,
        [
            patch.object(_mod.sc.tl, "leiden", leiden_mock),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save([])),
        ],
    )

    assert adata.obs["leiden"].tolist() == [str(i % 4) for i in range(30)]
    assert "spatial_domain" not in adata.obs, "domain_method='none' must not add a domain column"
    assert "spatial_domain" not in adata.uns, "domain_method='none' must not set domain metadata"
    leiden_mock.assert_not_called()


def test_domain_leiden_spatial_adds_labels_on_spatial_graph(tmp_path) -> None:
    """domain_method='leiden_spatial' → Leiden on spatial_connectivities.

    obs['spatial_domain'] present with the configured resolution, uns metadata
    recorded, and the round-1 obs['leiden'] column preserved.
    """
    adata = _make_adata()
    cfg = _make_cfg(domain_method="leiden_spatial", domain_resolution=1.0, h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    leiden_mock = MagicMock(side_effect=_fake_leiden)
    captured: list[AnnData] = []

    _run_main(
        adata,
        cfg,
        [
            patch.object(_mod.sc.tl, "leiden", leiden_mock),
            patch.object(_mod, "safe_write", side_effect=_capture_adata_on_save(captured)),
        ],
    )

    # Round-1 cluster column untouched by the domain phase
    assert adata.obs["leiden"].tolist() == [str(i % 4) for i in range(30)]

    # Leiden ran on the spatial graph with the configured resolution
    assert leiden_mock.called, "leiden_spatial must run sc.tl.leiden"
    kwargs = leiden_mock.call_args.kwargs
    assert kwargs["adjacency"] is adata.obsp["spatial_connectivities"], (
        "leiden must run on the spatial_connectivities graph"
    )
    assert kwargs["resolution"] == 1.0
    assert kwargs["key_added"] == "spatial_domain"

    # Domain label + metadata
    assert "spatial_domain" in adata.obs
    assert len(adata.obs["spatial_domain"]) == adata.n_obs
    assert adata.uns["spatial_domain"]["method"] == "leiden_spatial"
    assert adata.uns["spatial_domain"]["resolution"] == 1.0

    # Domain-bearing checkpoint written (first save carries obs + uns)
    assert len(captured) >= 1, "safe_write must persist the spatial-domain checkpoint"
    saved = captured[0]
    assert "spatial_domain" in saved.obs
    assert saved.uns["spatial_domain"]["method"] == "leiden_spatial"


def test_domain_resolution_respected(tmp_path) -> None:
    """domain_resolution is passed through to the spatial leiden call."""
    adata = _make_adata()
    cfg = _make_cfg(domain_method="leiden_spatial", domain_resolution=2.0, h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    leiden_mock = MagicMock(side_effect=_fake_leiden)

    _run_main(
        adata,
        cfg,
        [patch.object(_mod.sc.tl, "leiden", leiden_mock)],
    )

    assert leiden_mock.called
    assert leiden_mock.call_args.kwargs["resolution"] == 2.0


def test_domain_stagate_unavailable_falls_back_to_leiden_spatial(tmp_path) -> None:
    """STAGATE import gate raises → log.warning + graceful leiden_spatial fallback.

    The run must NOT crash; uns['spatial_domain'] records the fallback reason.
    """
    adata = _make_adata()
    cfg = _make_cfg(domain_method="stagate", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    log = MagicMock()
    leiden_mock = MagicMock(side_effect=_fake_leiden)

    _run_main(
        adata,
        cfg,
        [
            patch.object(
                _mod, "_require_stagate", side_effect=ImportError("STAGATE not installed")
            ),
            patch.object(_mod, "setup_logger", return_value=log),
            patch.object(_mod.sc.tl, "leiden", leiden_mock),
        ],
    )

    assert log.warning.called, "unavailable STAGATE must be logged as a warning"
    assert leiden_mock.called, "fallback must run leiden_spatial on the spatial graph"
    assert leiden_mock.call_args.kwargs["adjacency"] is adata.obsp["spatial_connectivities"]
    assert adata.uns["spatial_domain"] == {
        "method": "leiden_spatial",
        "reason": "stagate unavailable",
    }
    assert "spatial_domain" in adata.obs


def test_domain_stagate_success_clusters_embedding(tmp_path) -> None:
    """STAGATE available (mocked) → spatial net + train → Leiden on embedding.

    The min() n_pcs guard caps the STAGATE neighbors n_pcs at the embedding
    width (20 < cfg.pca.n_pcs_use=50).
    """
    adata = _make_adata()
    adata.obsm["STAGATE"] = np.random.RandomState(7).randn(adata.n_obs, 20)
    cfg = _make_cfg(domain_method="stagate", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    stagate_mod = MagicMock()
    neighbors_mock = MagicMock()
    leiden_mock = MagicMock(side_effect=_fake_leiden)

    _run_main(
        adata,
        cfg,
        [
            patch.object(_mod, "_require_stagate"),  # import gate passes
            patch.dict("sys.modules", {"STAGATE": stagate_mod}),
            patch.object(_mod.sc.pp, "neighbors", neighbors_mock),
            patch.object(_mod.sc.tl, "leiden", leiden_mock),
        ],
    )

    stagate_mod.Cal_Spatial_Net.assert_called_once()
    stagate_mod.train_STAGATE.assert_called_once()

    # Embedding clustered with the n_pcs min() guard applied
    neighbors_mock.assert_called_once()
    assert neighbors_mock.call_args.kwargs["use_rep"] == "STAGATE"
    assert neighbors_mock.call_args.kwargs["n_pcs"] == 20, (
        "n_pcs must be capped to the STAGATE embedding width (min() guard)"
    )

    assert adata.uns["spatial_domain"]["method"] == "stagate"
    assert "spatial_domain" in adata.obs
    assert adata.obs["leiden"].tolist() == [str(i % 4) for i in range(30)]


def test_domain_missing_connectivities_graceful(tmp_path) -> None:
    """No spatial_connectivities → log.warning, no crash, domain skipped.

    uns['spatial_domain'] records the skip reason; obs is left untouched.
    """
    adata = _make_adata(with_connectivities=False)
    cfg = _make_cfg(domain_method="leiden_spatial", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
    log = MagicMock()
    leiden_mock = MagicMock()

    _run_main(
        adata,
        cfg,
        [
            patch.object(_mod, "setup_logger", return_value=log),
            patch.object(_mod.sc.tl, "leiden", leiden_mock),
        ],
    )

    assert log.warning.called, "missing spatial graph must be logged as a warning"
    leiden_mock.assert_not_called()
    assert adata.uns["spatial_domain"] == {
        "method": "none",
        "reason": "spatial_connectivities missing",
    }
    assert "spatial_domain" not in adata.obs


def test_checkpoint_before_domain_plot(tmp_path) -> None:
    """safe_write of the domain checkpoint MUST precede the domain plot.

    Hard convention: checkpoint-before-plot (no plotting over unsaved labels).
    """
    adata = _make_adata()
    cfg = _make_cfg(domain_method="leiden_spatial", h5ad_dir=str(tmp_path))
    _ensure_input(tmp_path)
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
            patch.object(_mod.sc.tl, "leiden", MagicMock(side_effect=_fake_leiden)),
            patch.object(_mod, "safe_write", side_effect=_record_save),
            patch.object(_mod, "safe_plot", side_effect=_record_plot),
        ],
    )

    assert call_log.count("save") >= 1, "domain checkpoint must be written"
    assert call_log.count("plot") == 1, "exactly one spatial-domain plot expected"
    assert call_log.index("save") < call_log.index("plot"), (
        "checkpoint must be written BEFORE the domain plot"
    )


def test_effective_n_pcs_min_guard() -> None:
    """_effective_n_pcs caps a requested n_pcs at the embedding width."""
    assert _mod._effective_n_pcs(50, 20) == 20
    assert _mod._effective_n_pcs(50, 60) == 50
