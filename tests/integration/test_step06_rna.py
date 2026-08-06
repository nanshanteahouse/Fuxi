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
from unittest.mock import Mock, patch

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
            "figure_format": "pdf",  # explicit project override (global.yaml defaults to png)
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
    cell_subtype: str | None = None,
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

    # Plan todo 8: parent X_umap so the pre-subcluster subtype UMAP has
    # coordinates to draw. Values are irrelevant — safe_plot only needs the
    # embedding to exist, and the step 06 reprocessing block recomputes
    # X_umap on the subset later anyway.
    adata.obsm["X_umap"] = rng.randn(n_cells, 2).astype(np.float32)

    if cell_subtype is not None:
        # Plan todo 7 fixture: a pre-existing cell_subtype column with MIXED
        # values so per-type writeback sees BOTH classes (resolved labels to
        # preserve + "unresolved" cells that may be overwritten). Pattern i % 4
        # yields ~50/50 resolved/unresolved within EVERY cell_type.
        if cell_subtype == "mixed2":
            # Plan todo 8 variant: TWO distinct resolved labels per cell type so
            # the guarded pre-subcluster subtype UMAP (≥2 distinct values outside
            # {"unresolved", "N/A", ""}) fires within BOTH Type_0 and Type_1.
            # Pattern i % 3 gives every cell type exactly {RGC_Foxp2, RGC_W3,
            # unresolved}.
            subtypes = [
                "unresolved" if i % 3 == 0 else "RGC_Foxp2" if i % 3 == 1 else "RGC_W3"
                for i in range(n_cells)
            ]
        else:
            subtypes = ["RGC_Foxp2" if i % 4 in (1, 2) else "unresolved" for i in range(n_cells)]
        adata.obs["cell_subtype"] = pd.Categorical(subtypes)
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

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_plots_cell_subtype_umap(self, tmp_path: Path) -> None:
        """Plan todo 8: pre-subcluster UMAP colored by cell_subtype → PDF produced."""
        import yaml

        # ── Stage input data WITH cell_subtype + X_umap ────────────────
        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad",
            n_cells=120,
            n_genes=200,
            n_types=2,
            cell_subtype="mixed2",
        )

        # ── Write config ────────────────────────────────────────────────
        cfg_dict = _build_minimal_step06_config(tmp_path)
        config_path = tmp_path / "config_subtype.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        # ── Run step 06 ────────────────────────────────────────────────
        result = _run_step06(config_path)
        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])
        assert result.returncode == 0, (
            f"Step 06 exited {result.returncode}. See stdout/stderr above."
        )

        # ── Verify the pre-subcluster subtype UMAP exists (format from config) ──
        fig_path = (
            tmp_path
            / "results"
            / "figures"
            / "06_subcluster"
            / f"sub_Type_0_cell_subtype_umap.{cfg_dict['plot']['figure_format']}"
        )
        assert fig_path.exists(), (
            f"expected cell_subtype UMAP at {fig_path} — not produced by step 06"
        )

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_skips_subtype_umap_when_column_absent(self, tmp_path: Path) -> None:
        """Plan todo 8: no cell_subtype column → clean skip, NO subtype UMAP PDF."""
        import yaml

        # ── Stage input data WITHOUT cell_subtype ───────────────────────
        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad",
            n_cells=120,
            n_genes=200,
            n_types=2,
        )

        # ── Write config ────────────────────────────────────────────────
        cfg_dict = _build_minimal_step06_config(tmp_path)
        config_path = tmp_path / "config_nosubtype.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        # ── Run step 06 ────────────────────────────────────────────────
        result = _run_step06(config_path)
        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])
        assert result.returncode == 0, (
            f"Step 06 exited {result.returncode}. See stdout/stderr above."
        )

        # ── Verify the skip path produced no subtype UMAP (any format) ──
        _fmt = cfg_dict["plot"]["figure_format"]
        no_subtype_fig = (
            tmp_path
            / "results"
            / "figures"
            / "06_subcluster"
            / f"sub_Type_0_cell_subtype_umap.{_fmt}"
        )
        assert not no_subtype_fig.exists(), (
            f"subtype UMAP must NOT exist without cell_subtype: {no_subtype_fig}"
        )


# ── Test: writeback in-place (function-level, Item 1.3) ─────────────────


def _load_step06_module():
    """Load rna/steps/06_subcluster.py as a module (fresh object per call)."""
    import importlib.util

    step_path = _REPO_ROOT / "rna" / "steps" / "06_subcluster.py"
    spec = importlib.util.spec_from_file_location("rna.steps._06_subcluster_t", str(step_path))
    assert spec is not None and spec.loader is not None, "cannot build spec for 06_subcluster.py"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _DuckCfg:
    """Minimal config duck for auto_writeback (only touches these fields)."""

    def __init__(self, incremental_io: bool = True) -> None:
        self.incremental_io = incremental_io
        self.verify_write_integrity = False


class TestStep06Writeback:
    """Function-level tests for auto_writeback's obs write path (Item 1.3)."""

    @staticmethod
    def _make_main(
        tmp_path: Path, n_cells: int = 60, n_genes: int = 50, cell_subtype: str | None = None
    ):
        path = tmp_path / "results" / "h5ad" / "05_annotated.h5ad"
        _create_synthetic_annotated_h5ad(
            path,
            n_cells=n_cells,
            n_genes=n_genes,
            n_types=2,
            cell_subtype=cell_subtype,
        )
        import scanpy as sc

        return path, sc.read(str(path))

    @staticmethod
    def _make_sub(main):
        sub = main[main.obs["cell_type"].astype(str) == "Type_0"].copy()
        import pandas as pd

        half = sub.n_obs // 2
        sub.obs["sub_ai_label"] = pd.Categorical(
            ["Subcluster_0"] * half + ["Subcluster_1"] * (sub.n_obs - half)
        )
        return sub

    def test_writeback_inplace_passes_only_cell_subtype(self, tmp_path: Path) -> None:
        """In-place path passes a DataFrame with exactly the cell_subtype column."""
        mod = _load_step06_module()
        main_path, main = self._make_main(tmp_path)
        sub = self._make_sub(main)

        with (
            patch("core.utils.write_obs_columns_inplace") as mock_inplace,
            patch.object(mod, "safe_write") as mock_safe,
        ):
            n = mod.auto_writeback(
                sub, "Type_0", str(main_path), log=None, cfg=_DuckCfg(incremental_io=True)
            )

        assert n == sub.n_obs, f"expected {sub.n_obs} Type_0 cells written back, got {n}"
        mock_inplace.assert_called_once()
        written = mock_inplace.call_args.args[1]
        assert list(written.columns) == ["cell_subtype"], "only cell_subtype may be written back"
        assert len(written) == main.n_obs, "obs_df must be full-length (aligned with file n_obs)"
        mock_safe.assert_not_called()

    def test_writeback_inplace_preserves_main_file(self, tmp_path: Path) -> None:
        """Real run: X and pre-existing obs columns survive; cell_subtype is the only change."""
        import numpy as np

        mod = _load_step06_module()
        main_path, before = self._make_main(tmp_path)
        sub = self._make_sub(before)
        x_before = np.asarray(before.X)

        n = mod.auto_writeback(
            sub, "Type_0", str(main_path), log=None, cfg=_DuckCfg(incremental_io=True)
        )

        assert n == sub.n_obs
        import scanpy as sc

        after = sc.read(str(main_path))
        assert set(after.obs.columns) == set(before.obs.columns) | {"cell_subtype"}
        assert np.array_equal(np.asarray(after.X), x_before), "X must not change"
        assert (
            after.obs["cell_type"].astype(str).to_numpy().tolist()
            == before.obs["cell_type"].astype(str).to_numpy().tolist()
        )
        got = (
            after.obs.loc[sub.obs_names.astype(str), "cell_subtype"]
            .astype(str)
            .to_numpy()
            .tolist()
        )
        assert got == sub.obs["sub_ai_label"].astype(str).to_numpy().tolist()

    def test_writeback_full_fallback_when_incremental_disabled(self, tmp_path: Path) -> None:
        """incremental_io=False falls back to the full safe_write path."""
        mod = _load_step06_module()
        main_path, main = self._make_main(tmp_path)
        sub = self._make_sub(main)

        with (
            patch("core.utils.write_obs_columns_inplace") as mock_inplace,
            patch.object(mod, "safe_write") as mock_safe,
        ):
            n = mod.auto_writeback(
                sub, "Type_0", str(main_path), log=None, cfg=_DuckCfg(incremental_io=False)
            )

        assert n == sub.n_obs
        mock_safe.assert_called_once()
        mock_inplace.assert_not_called()

    def test_writeback_noop_without_ai_label(self, tmp_path: Path) -> None:
        """No sub_ai_label → auto_writeback is a no-op (no write of any kind)."""
        mod = _load_step06_module()
        main_path, main = self._make_main(tmp_path)
        sub = main[main.obs["cell_type"].astype(str) == "Type_0"].copy()

        with (
            patch("core.utils.write_obs_columns_inplace") as mock_inplace,
            patch.object(mod, "safe_write") as mock_safe,
        ):
            n = mod.auto_writeback(
                sub, "Type_0", str(main_path), log=None, cfg=_DuckCfg(incremental_io=True)
            )

        assert n == 0
        mock_inplace.assert_not_called()
        mock_safe.assert_not_called()

    def test_writeback_preserves_resolved_subtype_labels(self, tmp_path: Path) -> None:
        """Plan todo 7 / decision D6: resolved cell_subtype survives; only unresolved cells overwritten."""
        mod = _load_step06_module()
        main_path, before = self._make_main(tmp_path, cell_subtype="mixed")
        sub = self._make_sub(before)

        # Classify Type_0 barcodes from the pre-writeback column.
        before_sub = (
            before.obs.loc[sub.obs_names.astype(str), "cell_subtype"].astype(str).to_dict()
        )
        resolved_bcs = [bc for bc, v in before_sub.items() if v not in ("unresolved", "N/A")]
        unresolved_bcs = [bc for bc, v in before_sub.items() if v == "unresolved"]
        assert resolved_bcs and unresolved_bcs, "mixed fixture must contain both classes"

        n = mod.auto_writeback(
            sub, "Type_0", str(main_path), log=None, cfg=_DuckCfg(incremental_io=True)
        )
        # n_updated counts ONLY overwritten (unresolved) cells.
        assert n == len(unresolved_bcs), f"expected {len(unresolved_bcs)} updated, got {n}"

        import scanpy as sc

        after = sc.read(str(main_path))
        after_sub = after.obs["cell_subtype"].astype(str).to_dict()
        # (a) resolved cells KEEP their value after writeback.
        for bc in resolved_bcs:
            assert after_sub[bc] == before_sub[bc], f"resolved label for {bc} was overwritten"
        # (b) unresolved cells get the new sub_ai_label.
        expected = sub.obs["sub_ai_label"].astype(str).to_dict()
        for bc in unresolved_bcs:
            assert after_sub[bc] == expected[bc], f"unresolved cell {bc} did not get sub_ai_label"

    def test_writeback_missing_column_bootstraps_and_overwrites(self, tmp_path: Path) -> None:
        """(c) Missing cell_subtype → bootstrap (copy of cell_type) → parent-copy cells ARE overwritten."""
        mod = _load_step06_module()
        main_path, before = self._make_main(tmp_path)  # no cell_subtype column
        assert "cell_subtype" not in before.obs
        sub = self._make_sub(before)

        n = mod.auto_writeback(
            sub, "Type_0", str(main_path), log=None, cfg=_DuckCfg(incremental_io=True)
        )
        assert n == sub.n_obs, (
            f"all bootstrapped (parent-copy) cells should be overwritten, got {n}"
        )
        import scanpy as sc

        after = sc.read(str(main_path))
        got = (
            after.obs.loc[sub.obs_names.astype(str), "cell_subtype"]
            .astype(str)
            .to_numpy()
            .tolist()
        )
        assert got == sub.obs["sub_ai_label"].astype(str).to_numpy().tolist()


# ── Test: KB-constrained AI annotation (plan todo 6) ────────────────────


def _make_ai_sub(n_cells: int = 8, n_leiden: int = 2):
    """Minimal AnnData subset with a 'leiden' column for AI annotation tests."""
    import numpy as np
    import scanpy as sc

    sub = sc.AnnData(X=np.zeros((n_cells, 5), dtype=np.float32))
    sub.obs["leiden"] = np.array([str(i % n_leiden) for i in range(n_cells)])
    return sub


class TestStep06KbConstrainedAnnotation:
    """Plan todo 6: KB-constrained AI prompt injection + label canonicalization."""

    _CANDIDATES = ["RGC_Foxp2", "RGC_W3"]

    @staticmethod
    def _cfg(tissue_kb: str = "retina", subcluster_kb_constrained: bool = True):
        from types import SimpleNamespace

        return SimpleNamespace(
            tissue_kb=tissue_kb,
            ai=SimpleNamespace(subcluster_kb_constrained=subcluster_kb_constrained),
        )

    @staticmethod
    def _args(cell_type: str = "RGC"):
        from types import SimpleNamespace

        return SimpleNamespace(cell_type=cell_type)

    def test_prompt_injects_kb_candidates_when_they_exist(self):
        """kb_candidates receives build_subtype_candidates output when non-empty."""
        mod = _load_step06_module()
        sub = _make_ai_sub()
        cfg = self._cfg()
        args = self._args()
        log = Mock()
        ai_response = (
            '{"0": {"cell_type": "RGC", "subtype": "Foxp2"}, '
            '"1": {"cell_type": "RGC", "subtype": "W3"}}'
        )
        with (
            patch("core.kb.load_kb", return_value={"_hierarchy": {}}),
            patch.object(mod, "build_subtype_candidates", return_value=self._CANDIDATES),
            patch(
                "core.ai.prompts.build_annotation_prompt", return_value=("sys", "user")
            ) as mock_bap,
            patch("core.ai.caller.ai_query", return_value=ai_response),
            patch.object(mod, "safe_plot"),
        ):
            mod._ai_subcluster_annotation(sub, cfg, args, log)

        mock_bap.assert_called_once()
        kwargs = mock_bap.call_args.kwargs
        assert kwargs["kb_candidates"] == self._CANDIDATES
        # subcluster_kb_constrained defaults True → unconstrained=False.
        assert kwargs["unconstrained"] is False

    def test_unconstrained_true_when_subcluster_kb_constrained_false(self):
        """Flag False → unconstrained=True is passed to build_annotation_prompt."""
        mod = _load_step06_module()
        sub = _make_ai_sub()
        cfg = self._cfg(subcluster_kb_constrained=False)
        args = self._args()
        log = Mock()
        with (
            patch("core.kb.load_kb", return_value={"_hierarchy": {}}),
            patch.object(mod, "build_subtype_candidates", return_value=self._CANDIDATES),
            patch(
                "core.ai.prompts.build_annotation_prompt", return_value=("sys", "user")
            ) as mock_bap,
            patch(
                "core.ai.caller.ai_query",
                return_value='{"0": {"cell_type": "RGC", "subtype": "Foxp2"}}',
            ),
            patch.object(mod, "safe_plot"),
        ):
            mod._ai_subcluster_annotation(sub, cfg, args, log)

        mock_bap.assert_called_once()
        assert mock_bap.call_args.kwargs["kb_candidates"] == self._CANDIDATES
        assert mock_bap.call_args.kwargs["unconstrained"] is True

    def test_no_kb_candidates_when_candidates_empty(self):
        """Empty candidates → no kb_candidates, current behavior unchanged."""
        mod = _load_step06_module()
        sub = _make_ai_sub()
        cfg = self._cfg()
        args = self._args()
        log = Mock()
        with (
            patch("core.kb.load_kb", return_value={"_hierarchy": {}}),
            patch.object(mod, "build_subtype_candidates", return_value=[]),
            patch(
                "core.ai.prompts.build_annotation_prompt", return_value=("sys", "user")
            ) as mock_bap,
            patch(
                "core.ai.caller.ai_query",
                return_value='{"0": {"cell_type": "RGC", "subtype": "Foxp2"}}',
            ),
            patch.object(mod, "safe_plot"),
        ):
            mod._ai_subcluster_annotation(sub, cfg, args, log)

        mock_bap.assert_called_once()
        assert "kb_candidates" not in mock_bap.call_args.kwargs
        assert "unconstrained" not in mock_bap.call_args.kwargs

    def test_load_kb_failure_falls_back_unconstrained(self):
        """load_kb raising ValueError → warn, no kb_candidates, AI call still made."""
        mod = _load_step06_module()
        sub = _make_ai_sub()
        cfg = self._cfg(tissue_kb="ghost")
        args = self._args()
        log = Mock()
        with (
            patch("core.kb.load_kb", side_effect=ValueError("unsupported tissue")),
            patch(
                "core.ai.prompts.build_annotation_prompt", return_value=("sys", "user")
            ) as mock_bap,
            patch(
                "core.ai.caller.ai_query",
                return_value='{"0": {"cell_type": "RGC", "subtype": "Foxp2"}}',
            ) as mock_aq,
            patch.object(mod, "safe_plot"),
        ):
            mod._ai_subcluster_annotation(sub, cfg, args, log)

        log.warning.assert_any_call("KB load failed — unconstrained subcluster annotation")
        mock_bap.assert_called_once()
        assert "kb_candidates" not in mock_bap.call_args.kwargs
        assert "unconstrained" not in mock_bap.call_args.kwargs
        mock_aq.assert_called_once()

    # ── Canonicalization of the already-sanitized label ──

    def test_canonicalize_casefold_label_to_candidate(self):
        mod = _load_step06_module()
        log = Mock()
        assert mod.canonicalize_subtype_label("rgc_foxp2", self._CANDIDATES, log) == "RGC_Foxp2"
        log.info.assert_not_called()

    def test_canonicalize_delimiter_to_underscore(self):
        mod = _load_step06_module()
        log = Mock()
        assert mod.canonicalize_subtype_label("RGC-Foxp2", self._CANDIDATES, log) == "RGC_Foxp2"

    def test_canonicalize_parentheses_stripped(self):
        mod = _load_step06_module()
        log = Mock()
        assert mod.canonicalize_subtype_label("RGC(Foxp2)", self._CANDIDATES, log) == "RGC_Foxp2"

    def test_canonicalize_word_order_reversal_kept_raw_and_logs_novel(self):
        mod = _load_step06_module()
        log = Mock()
        assert mod.canonicalize_subtype_label("Foxp2_RGC", self._CANDIDATES, log) == "Foxp2_RGC"
        log.info.assert_any_call("NOVEL: %s", "Foxp2_RGC")

    def test_ai_labels_canonicalized_through_parse_loop(self):
        """The parse loop applies canonicalization to final sanitized labels."""
        mod = _load_step06_module()
        sub = _make_ai_sub()
        cfg = self._cfg()
        args = self._args()
        log = Mock()
        ai_response = (
            '{"0": {"cell_type": "RGC", "subtype": "foxp2"}, '
            '"1": {"cell_type": "RGC", "subtype": "Foxp2"}}'
        )
        with (
            patch("core.kb.load_kb", return_value={"_hierarchy": {}}),
            patch.object(mod, "build_subtype_candidates", return_value=self._CANDIDATES),
            patch("core.ai.prompts.build_annotation_prompt", return_value=("sys", "user")),
            patch("core.ai.caller.ai_query", return_value=ai_response),
            patch.object(mod, "safe_plot"),
        ):
            mod._ai_subcluster_annotation(sub, cfg, args, log)

        labels = dict(zip(sub.obs["leiden"], sub.obs["sub_ai_label"]))
        # "RGC_foxp2" (casefold) and "RGC_Foxp2" (exact) both land on the candidate.
        assert labels["0"] == "RGC_Foxp2"
        assert labels["1"] == "RGC_Foxp2"


# ── Test: sentinel contract (Item 1.6) ─────────────────────────────────


def _ai_enabled_config(tmp_path: Path) -> dict:
    """Step-06 config that deterministically reaches a real writeback.

    build_annotation_prompt is called with precomputed_rank=True but the
    subset never ran rank_genes_groups → KeyError → the step falls back to
    numeric ``sub_ai_label`` → auto_writeback actually writes. The API endpoint
    is a closed local port so even a hypothetical live call fails fast.
    """
    cfg_dict = _build_minimal_step06_config(tmp_path)
    cfg_dict["ai"] = {
        "enabled": True,
        "subcluster": True,
        "api_base": "http://127.0.0.1:1",
        "api_key": "sk-test-invalid",
        "model": "test-model",
        "max_tokens": 256,
        "temperature": 0.1,
        "timeout": 5,
        "thinking_enabled": False,
    }
    return cfg_dict


def _run_step06(config_path: Path, timeout: int = 120) -> subprocess.CompletedProcess:
    """Run step 06 as a subprocess via the pipeline runner."""
    return subprocess.run(
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
        timeout=timeout,
    )


class TestStep06Sentinel:
    """Item 1.6 sentinel contract: written only after a real writeback."""

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_sentinel_written_after_writeback(self, tmp_path: Path) -> None:
        """Incremental (default) mode writes a non-empty sentinel after writeback."""
        import yaml

        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad", n_cells=120, n_genes=200, n_types=2
        )
        cfg_dict = _ai_enabled_config(tmp_path)
        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        result = _run_step06(config_path)
        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])
        assert result.returncode == 0, f"Step 06 exited {result.returncode}. See output above."

        sentinel = h5ad_dir / "05_annotated.h5ad.step06_done"
        assert sentinel.exists(), f"sentinel {sentinel} not created"
        assert sentinel.read_text().strip() == "done", "sentinel must be non-empty ('done')"

        import scanpy as sc

        main = sc.read(str(h5ad_dir / "05_annotated.h5ad"))
        assert "cell_subtype" in main.obs, "writeback must have created cell_subtype"
        changed = main.obs["cell_subtype"].astype(str) != main.obs["cell_type"].astype(str)
        assert changed.any(), "writeback must have changed cell_subtype on subclustered cells"
        # In-place writeback removes its .bak after success
        assert not (h5ad_dir / "05_annotated.h5ad.bak").exists()

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_sentinel_full_write_when_incremental_disabled(self, tmp_path: Path) -> None:
        """incremental_io=False: full safe_write still writes the sentinel, same result."""
        import yaml

        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad", n_cells=120, n_genes=200, n_types=2
        )
        cfg_dict = _ai_enabled_config(tmp_path)
        cfg_dict["incremental_io"] = False
        config_path = tmp_path / "config_full.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        result = _run_step06(config_path)
        if result.returncode != 0:
            print("=== STDOUT ===")
            print(result.stdout[-3000:])
            print("=== STDERR ===")
            print(result.stderr[-3000:])
        assert result.returncode == 0, f"Step 06 exited {result.returncode}. See output above."

        sentinel = h5ad_dir / "05_annotated.h5ad.step06_done"
        assert sentinel.exists(), f"sentinel {sentinel} not created"
        assert sentinel.read_text().strip() == "done", "sentinel must be non-empty ('done')"

        import scanpy as sc

        # Results consistent with in-place mode: cell_subtype mirrors sub_ai_label
        sub_files = sorted(h5ad_dir.glob("05_sub_*.h5ad"))
        assert sub_files, "no 05_sub_*.h5ad output"
        sub_adata = sc.read(str(sub_files[0]))
        labels = dict(
            zip(sub_adata.obs_names.astype(str), sub_adata.obs["sub_ai_label"].astype(str))
        )
        main = sc.read(str(h5ad_dir / "05_annotated.h5ad"))
        type0 = main.obs["cell_type"].astype(str) == "Type_0"
        got = main.obs.loc[type0, "cell_subtype"].astype(str)
        mismatches = [bc for bc, v in got.items() if labels[bc] != v]
        assert not mismatches, (
            f"{len(mismatches)} Type_0 cells disagree with sub_ai_label (e.g. {mismatches[:3]})"
        )
        # Full write leaves no .bak behind
        assert not (h5ad_dir / "05_annotated.h5ad.bak").exists()

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_step06_exit2_skip_writes_no_sentinel(self, tmp_path: Path) -> None:
        """No --cell-type and no subcluster_types → exit 2 skip, NO sentinel."""
        import yaml

        h5ad_dir = tmp_path / "results" / "h5ad"
        _create_synthetic_annotated_h5ad(
            h5ad_dir / "05_annotated.h5ad", n_cells=120, n_genes=200, n_types=2
        )
        cfg_dict = _build_minimal_step06_config(tmp_path)
        cfg_dict["marker"]["subcluster_types"] = []
        config_path = tmp_path / "config_skip.yaml"
        with open(config_path, "w") as f:
            yaml.dump(cfg_dict, f)

        result = _run_step06(config_path, timeout=60)
        # runner absorbs exit 2 as "skipped" and exits 0 itself
        assert result.returncode == 0, (
            f"Step 06 exit-2 skip should be absorbed by the runner, got {result.returncode}"
        )
        assert not (h5ad_dir / "05_annotated.h5ad.step06_done").exists(), (
            "sentinel must NOT be written when subcluster is skipped (exit 2)"
        )
