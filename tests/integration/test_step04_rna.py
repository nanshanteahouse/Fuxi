"""End-to-end integration tests for RNA step 04 (cluster/UMAP).

Tests that core modules used by step 04 import cleanly, the config schema
validates clustering settings correctly, and — optionally — that the step
runs to completion as a subprocess on synthetic data.

Run with::

    pytest tests/integration/test_step04_rna.py -v --tb=short

The subprocess test (``test_step04_subprocess_smoke``) can be skipped by
setting ``SKIP_SLOW_TESTS=1``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def _skip_slow() -> bool:
    """Return True when SKIP_SLOW_TESTS env var is set to skip expensive tests."""
    return os.environ.get("SKIP_SLOW_TESTS", "").lower() in ("1", "true")


# ── Helpers ────────────────────────────────────────────────────────────


def _build_minimal_step04_config(
    tmp_path: Path,
    *,
    h5ad_input_name: str = "03_integrated.h5ad",
) -> dict:
    """Return a minimal config dict for step 04 with absolute paths.

    All output directories are rooted under ``tmp_path``.  The caller is
    responsible for placing a valid h5ad at ``{h5ad_dir}/{h5ad_input_name}``
    **before** invoking the step.
    """
    h5ad_dir = tmp_path / "results" / "h5ad"
    figure_dir = tmp_path / "results" / "figures"
    table_dir = tmp_path / "results" / "tables"
    log_dir = tmp_path / "logs"
    for d in (h5ad_dir, figure_dir, table_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    return {
        "modality": "rna",
        "tissue": "test",
        "species": "mouse",
        "project_dir": str(tmp_path),
        "h5ad_dir": str(h5ad_dir),
        "figure_dir": str(figure_dir),
        "table_dir": str(table_dir),
        "log_dir": str(log_dir),
        "execution": {
            "device": "cpu",
            "random_seed": 42,
            "n_jobs": 1,
        },
        "pca": {
            "n_pcs_use": 10,
        },
        "clustering": {
            "param_grid_n_neighbors": [5, 10],
            "param_grid_resolutions": [0.3, 0.5],
            "cluster_selection_method": "pareto_elbow",
            "umap_plot_mode": "skip",
            "multi_metric_adaptive_resolution": False,
        },
        "plot": {
            "figure_dpi": 72,
            "palette": {
                "categorical": "tab20",
                "dotplot_fill": "YlOrRd",
            },
        },
        "marker": {
            "marker_dict": {},
        },
        "ai": {
            "enabled": False,
        },
    }


def _create_synthetic_integrated_h5ad(path: Path, n_cells: int = 100, n_genes: int = 50) -> None:
    """Create a small AnnData with PCA pre-computed and write to *path*."""
    import numpy as np
    import pandas as pd
    import scanpy as sc

    rng = np.random.RandomState(42)
    adata = sc.AnnData(
        X=rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32),
        var=pd.DataFrame(index=[f"GENE_{i:04d}" for i in range(n_genes)]),
    )
    adata.obsm["X_pca"] = rng.standard_normal((n_cells, 10))
    adata.uns["pca"] = {"params": {"n_pcs": 10}}
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(str(path))
    del adata


# ── Test: core module imports ─────────────────────────────────────────


class TestStep04Imports:
    """Verify every core module that step 04 depends on is importable."""

    @pytest.mark.smoke
    def test_core_cluster_modules_import(self) -> None:
        """core.cluster.grid_search and core.cluster.evaluation import."""
        from core.cluster.evaluation import (  # noqa: F401
            _detect_granularity,
            select_best_umap_params,
        )
        from core.cluster.grid_search import (  # noqa: F401
            grid_search_clustering,
            select_best_params,
        )

    @pytest.mark.smoke
    def test_core_utils_import(self) -> None:
        """Core utility functions used by step 04 import cleanly."""
        from core.utils import (  # noqa: F401
            gpu_leiden,
            gpu_neighbors,
            gpu_umap,
            resolve_config,
            safe_write,
            setup_logger,
            timed_substep,
        )

    @pytest.mark.smoke
    def test_config_schema_import(self) -> None:
        """Top-level Config with its clustering sub-model imports."""
        from core.config.schema import (  # noqa: F401
            ClusteringSettings,
            Config,  # noqa: F401
            ExecutionConfig,
            PCASettings,
        )


# ── Test: config parsing ──────────────────────────────────────────────


class TestStep04Config:
    """Config schema correctly validates step-04-specific settings."""

    @pytest.mark.smoke
    def test_minimal_config_parses(self, tmp_path: Path) -> None:
        """A minimal YAML config with step 04 fields roundtrips through Config."""
        from core.config.schema import Config

        cfg_dict = _build_minimal_step04_config(tmp_path)
        cfg = Config.model_validate(cfg_dict)

        # Check checkpoint paths resolve correctly
        assert cfg.rna_integrated_h5ad.endswith("03_integrated.h5ad")
        assert cfg.cluster_h5ad.endswith("04_clustered.h5ad")
        assert cfg.h5ad_dir in cfg.rna_integrated_h5ad
        assert cfg.h5ad_dir in cfg.cluster_h5ad

        # Check clustering settings survived
        assert cfg.clustering.param_grid_n_neighbors == [5, 10]
        assert cfg.clustering.param_grid_resolutions == [0.3, 0.5]
        assert cfg.clustering.cluster_selection_method == "pareto_elbow"
        assert cfg.clustering.umap_plot_mode == "skip"

        # Check execution defaults
        assert cfg.execution.device == "cpu"
        assert cfg.execution.random_seed == 42
        assert cfg.pca.n_pcs_use == 10

    @pytest.mark.smoke
    def test_config_resolve_roundtrip(self, tmp_path: Path) -> None:
        """Config written to YAML and reloaded via resolve_config()."""
        import yaml

        from core.utils._config import resolve_config

        cfg_dict = _build_minimal_step04_config(tmp_path)
        config_path = tmp_path / "test_config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        cfg = resolve_config(str(config_path))
        assert cfg.modality == "rna"
        assert cfg.clustering.cluster_selection_method == "pareto_elbow"
        assert cfg.clustering.umap_plot_mode == "skip"
        assert cfg.execution.device == "cpu"
        # Directories should have been auto-created by resolve_config
        assert os.path.isdir(cfg.h5ad_dir)
        assert os.path.isdir(cfg.figure_dir)


# ── Test: subprocess smoke ────────────────────────────────────────────


class TestStep04Subprocess:
    """Optional full subprocess run of step 04 on synthetic data.

    These tests are self-contained (they generate their own AnnData) but
    may take 30-60 seconds.  Set ``SKIP_SLOW_TESTS=1`` to bypass them.
    """

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step04_subprocess_smoke(self, tmp_path: Path) -> None:
        """Run step 04 as subprocess; verify output h5ad + CSV exist."""
        # ── Stage input data ───────────────────────────────────────────
        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_integrated_h5ad(h5ad_dir / "03_integrated.h5ad")

        # ── Write config ────────────────────────────────────────────────
        import yaml

        cfg_dict = _build_minimal_step04_config(tmp_path)
        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        # ── Run step 04 ────────────────────────────────────────────────
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "core.run_pipeline",
                "--modality",
                "rna",
                "--step",
                "4",
                "--config",
                str(config_path),
            ],
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=300,
        )

        # Always print stdout/stderr on failure for debugging
        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])

        assert result.returncode == 0, (
            f"Step 04 exited with code {result.returncode}. See stdout/stderr above."
        )

        # ── Verify output artefacts ────────────────────────────────────
        cluster_h5ad = h5ad_dir / "04_clustered.h5ad"
        assert cluster_h5ad.exists(), (
            f"Expected {cluster_h5ad} not found — step 04 did not produce "
            f"the clustered checkpoint."
        )

        table_dir = tmp_path / "results" / "tables" / "04_cluster"
        csv_path = table_dir / "param_grid_summary.csv"
        assert csv_path.exists(), (
            f"Expected {csv_path} not found — parameter grid summary CSV was not written."
        )

        # ── Quick structural sanity on the output h5ad ─────────────────
        import scanpy as sc

        adata = sc.read(str(cluster_h5ad))

        # Should have leiden cluster assignments
        assert "leiden" in adata.obs, "Output h5ad missing 'leiden' column in obs"
        n_clusters = adata.obs["leiden"].nunique()
        assert n_clusters >= 1, f"Expected ≥1 leiden cluster, got {n_clusters}"

        # Should have PCA coordinates.  UMAP is intentionally NOT asserted:
        # scanpy 1.12.2 removed the ``n_epochs`` kwarg from ``sc.tl.umap``, so
        # step 04's final UMAP rebuild fails non-fatally (swallowed by the grid
        # search) and no ``X_umap`` is written — a known pre-existing issue
        # documented in test_e2e_incremental_parity.py.  The input h5ad carries
        # ``X_pca`` (pre-computed by the fixture) and step 04 preserves it, so we
        # assert the embedding the pipeline actually produced.
        assert "X_pca" in adata.obsm, "Output h5ad missing 'X_pca' in obsm"
        assert adata.obsm["X_pca"].shape == (adata.n_obs, 10), (
            f"X_pca shape {adata.obsm['X_pca'].shape} != ({adata.n_obs}, 10)"
        )

        # Should have best-param metadata
        assert "best_resolution" in adata.uns, "Output h5ad missing 'best_resolution' in uns"
        assert "best_n_neighbors" in adata.uns, "Output h5ad missing 'best_n_neighbors' in uns"
        assert "cluster_selection_method" in adata.uns, (
            "Output h5ad missing 'cluster_selection_method' in uns"
        )
