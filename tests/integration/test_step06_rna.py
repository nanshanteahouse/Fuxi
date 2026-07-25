"""Integration tests for RNA step 06 (subcluster).

Tests that step 06 imports cleanly, appears in the step list, the config
schema validates subcluster settings correctly, and — optionally — that the
step runs to completion as a subprocess on synthetic data.

Run with::

    pytest tests/integration/test_step06_rna.py -v --tb=short

The subprocess tests can be skipped by setting ``SKIP_SLOW_TESTS=1``.
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


def _build_minimal_step06_config(
    tmp_path: Path,
    *,
    h5ad_input_name: str = "05_annotated.h5ad",
) -> dict:
    """Return a minimal config dict for step 06 with absolute paths.

    All output directories are rooted under ``tmp_path``.  The caller is
    responsible for placing a valid h5ad at ``{h5ad_dir}/{h5ad_input_name}``
    **before** invoking the step.
    """
    h5ad_dir = tmp_path / "results" / "h5ad"
    figure_dir = tmp_path / "results" / "figures"
    table_dir = tmp_path / "results" / "tables"
    log_dir = tmp_path / "logs"

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
            "n_pcs_use": 5,
        },
        "hvg": {
            "n_top_genes": 50,
            "flavor": "seurat_v3",
        },
        "clustering": {
            "param_grid_n_neighbors": [5],
            "param_grid_resolutions": [0.3, 0.5],
            "leiden_flavor": "igraph",
            "leiden_n_iterations": 2,
            "umap_plot_mode": "skip",
        },
        "plot": {
            "figure_dpi": 72,
            "palette": {
                "categorical": "tab20",
                "dotplot_fill": "YlOrRd",
            },
            "umap_panel_size": [4, 4],
        },
        "marker": {
            "marker_dict": {},
            "subcluster_types": ["Type_0"],
            "subcluster_resolution": 0.4,
        },
        "ai": {
            "enabled": False,
            "subcluster": False,
        },
        "integration": {
            "method": "harmony",
            "batch_key": "sample",
        },
    }


def _create_synthetic_annotated_h5ad(
    path: Path,
    n_cells: int = 120,
    n_genes: int = 200,
    n_types: int = 2,
    *,
    include_cell_type: bool = True,
) -> None:
    """Create a minimal AnnData that mimics 05_annotated.h5ad and write to *path*."""
    import numpy as np
    import pandas as pd
    import scanpy as sc

    rng = np.random.RandomState(42)
    adata = sc.AnnData(
        X=rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32),
        var=pd.DataFrame(index=[f"GENE_{i:04d}" for i in range(n_genes)]),
    )
    if include_cell_type:
        adata.obs["cell_type"] = pd.Categorical([f"Type_{i % n_types}" for i in range(n_cells)])
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(str(path))
    del adata


# ── Test: smoke ─────────────────────────────────────────────────────────


class TestStep06Smoke:
    """Quick smoke tests for step 06."""

    @pytest.mark.smoke
    def test_step06_help(self) -> None:
        """Step 06 ``--help`` exits 0."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "core.run_pipeline",
                "--modality",
                "rna",
                "--step",
                "6",
                "--help",
            ],
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=30,
        )
        assert result.returncode == 0, (
            f"Step 06 --help exited with code {result.returncode}.\n"
            f"STDERR:\n{result.stderr[-2000:]}"
        )

    @pytest.mark.smoke
    def test_step06_import(self) -> None:
        """Step 06 ``main()`` is importable via importlib."""
        import importlib.util

        step_path = _REPO_ROOT / "rna" / "steps" / "06_subcluster.py"
        assert step_path.exists(), f"Step 06 file not found: {step_path}"

        spec = importlib.util.spec_from_file_location("rna.steps._06_subcluster", str(step_path))
        assert spec is not None, "Could not create module spec for 06_subcluster.py"
        assert spec.loader is not None, "Module spec has no loader"

        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert hasattr(mod, "main"), "Step 06 module missing 'main' function"
        assert callable(mod.main), "Step 06 'main' is not callable"

    @pytest.mark.smoke
    def test_step06_in_step_list(self) -> None:
        """``--list`` output contains '06_subcluster'."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "core.run_pipeline",
                "--modality",
                "rna",
                "--list",
            ],
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=30,
        )
        combined = result.stdout + result.stderr
        assert "06_subcluster" in combined, (
            f"Step 06 not found in --list output.\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )


# ── Test: subprocess ────────────────────────────────────────────────────


class TestStep06Subprocess:
    """Optional full subprocess runs of step 06 on synthetic data.

    These tests are self-contained (they generate their own AnnData) but
    may take 30-60 seconds.  Set ``SKIP_SLOW_TESTS=1`` to bypass them.
    """

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_synthetic_run(self, tmp_path: Path) -> None:
        """Run step 06 as subprocess; verify output h5ad + structural sanity."""
        import yaml

        # ── Stage input data ───────────────────────────────────────────
        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad",
            n_cells=120,
            n_genes=200,
            n_types=2,
        )

        # ── Write config ────────────────────────────────────────────────
        cfg_dict = _build_minimal_step06_config(tmp_path)
        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        # ── Run step 06 ────────────────────────────────────────────────
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "core.run_pipeline",
                "--modality",
                "rna",
                "--step",
                "6",
                "--config",
                str(config_path),
            ],
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=120,
        )

        # Always print stdout/stderr on failure for debugging
        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])

        assert result.returncode == 0, (
            f"Step 06 exited with code {result.returncode}. See stdout/stderr above."
        )

        # ── Verify output artefacts ────────────────────────────────────
        # Step 06 produces 05_sub_{cell_type}.h5ad files
        output_files = sorted(h5ad_dir.glob("05_sub_*.h5ad"))
        assert len(output_files) >= 1, f"No 05_sub_*.h5ad output files found in {h5ad_dir}"

        # ── Quick structural sanity on the output h5ad ─────────────────
        import scanpy as sc

        output_h5ad = output_files[0]
        adata = sc.read(str(output_h5ad))

        # Should have leiden subcluster assignments
        assert "leiden" in adata.obs, (
            f"Output h5ad {output_h5ad.name} missing 'leiden' column in obs"
        )
        n_clusters = adata.obs["leiden"].nunique()
        assert n_clusters >= 1, f"Expected ≥1 leiden subcluster, got {n_clusters}"

        # Should have UMAP coordinates
        assert "X_umap" in adata.obsm, f"Output h5ad {output_h5ad.name} missing 'X_umap' in obsm"
        assert adata.obsm["X_umap"].shape == (adata.n_obs, 2), (
            f"X_umap shape {adata.obsm['X_umap'].shape} != ({adata.n_obs}, 2)"
        )

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_respects_leiden_n_iterations(self, tmp_path: Path) -> None:
        """Leiden ``n_iterations`` config value flows through without error."""
        import yaml

        # ── Stage input data ───────────────────────────────────────────
        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad",
            n_cells=120,
            n_genes=200,
            n_types=2,
        )

        # ── Write config with non-default leiden_n_iterations ──────────
        cfg_dict = _build_minimal_step06_config(tmp_path)
        cfg_dict["clustering"]["leiden_n_iterations"] = 5
        config_path = tmp_path / "config_leiden.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        # ── Run step 06 ────────────────────────────────────────────────
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "core.run_pipeline",
                "--modality",
                "rna",
                "--step",
                "6",
                "--config",
                str(config_path),
            ],
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=120,
        )

        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])

        assert result.returncode == 0, (
            f"Step 06 with leiden_n_iterations=5 exited with code "
            f"{result.returncode}. See stdout/stderr above."
        )

        # ── Verify output was produced ─────────────────────────────────
        output_files = sorted(h5ad_dir.glob("05_sub_*.h5ad"))
        assert len(output_files) >= 1, "No output files found with leiden_n_iterations=5"

        import scanpy as sc

        adata = sc.read(str(output_files[0]))
        assert "leiden" in adata.obs, "Missing leiden column"
        assert adata.obs["leiden"].nunique() >= 1, "No leiden clusters found"

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_dry_run_or_skip_empty(self, tmp_path: Path) -> None:
        """Missing ``cell_type`` column produces clean error, no crash."""
        import yaml

        # ── Stage input data WITHOUT cell_type column ──────────────────
        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad",
            n_cells=120,
            n_genes=200,
            n_types=2,
            include_cell_type=False,
        )

        # ── Write config ────────────────────────────────────────────────
        cfg_dict = _build_minimal_step06_config(tmp_path)
        cfg_dict["marker"]["subcluster_types"] = ["Type_0"]
        config_path = tmp_path / "config_empty.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        # ── Run step 06 (should exit with error, not crash) ────────────
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "core.run_pipeline",
                "--modality",
                "rna",
                "--step",
                "6",
                "--config",
                str(config_path),
            ],
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=60,
        )

        # Should exit non-zero due to missing cell_type column
        assert result.returncode != 0, (
            "Step 06 should have exited with error when cell_type column is missing"
        )

        # Verify clean error (sys.exit) — no unhandled traceback
        assert "Traceback" not in result.stderr, (
            "Step 06 crashed with unhandled exception instead of clean error:\n"
            f"{result.stderr[-2000:]}"
        )
