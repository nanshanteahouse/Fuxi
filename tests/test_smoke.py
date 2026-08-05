"""Smoke tests — fast, no-heavy-dependency checks that gate CI and pre-push.

All tests here MUST:
  - Complete in < 3 seconds each (total suite < 15s)
  - Not require network, LLM, or real data (FUXI_DATA_ROOT optional)
  - Not depend on other test files

Tag: @pytest.mark.smoke
Run:  python -m pytest tests/ -m smoke -x --tb=short
"""

import concurrent.futures
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.smoke

REPO_ROOT = Path(__file__).resolve().parent.parent


# ── Core imports ────────────────────────────────────────────────────────


def test_core_config_imports() -> None:
    """Config classes import cleanly."""
    from core.config import Config  # noqa: F401
    from core.config.schema import DataInputConfig, QCSettings  # noqa: F401


def test_core_run_pipeline_module() -> None:
    """Pipeline runner module is importable."""
    import core.run_pipeline  # noqa: F401


# ── Modality imports ─────────────────────────────────────────────────────


def test_rna_module_imports() -> None:
    """scRNA-seq module imports without crash."""
    import rna  # noqa: F401


def test_atac_module_imports() -> None:
    """scATAC-seq module imports without crash."""
    import atac  # noqa: F401


def test_spatial_module_imports() -> None:
    """Spatial transcriptomics module imports without crash."""
    import spatial  # noqa: F401


# ── Pipeline discovery ───────────────────────────────────────────────────
# NOTE: All 3 modality --list invocations run in parallel via ThreadPoolExecutor.
# Each subprocess pays ~18s Python startup + anndata import on WSL2; serial = 54s,
# parallel = ~20s. Total smoke suite target: <30s to avoid GitHub SSH idle
# timeout (~60-90s) during pre-push hook.


def test_pipeline_list_all_modalities_parallel() -> None:
    """All 3 modality --list commands exit cleanly and print step list.

    Replaces three separate test_pipeline_list_{rna,atac,spatial} tests that ran
    sequentially. Parallelization cuts total wall from ~54s to ~20s, fitting
    within GitHub's SSH idle timeout during pre-push hook execution.
    """
    modalities = ["rna", "atac", "spatial"]

    def _run_one(modality: str) -> tuple[str, subprocess.CompletedProcess]:
        result = subprocess.run(
            [sys.executable, "-m", "core.run_pipeline", "--modality", modality, "--list"],
            capture_output=True,
            text=True,
            cwd=str(REPO_ROOT),
            timeout=60,
        )
        return modality, result

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        results = dict(executor.map(_run_one, modalities))

    for modality, result in results.items():
        assert result.returncode == 0, f"{modality} --list failed: {result.stderr}"
        assert "00_load" in result.stdout, f"{modality}: missing '00_load' in output"


# ── Config template discovery ────────────────────────────────────────────


def test_scaffold_renders_starter_configs() -> None:
    """Every format spec renders to a loadable starter config (smoke)."""
    from core.config.scaffold import render_template_text
    from core.preprocess.config_specs import materialized_specs

    specs = materialized_specs()
    assert len(specs) >= 4, f"Expected ≥4 format specs, found {len(specs)}"
    for spec in specs:
        data = yaml.safe_load(render_template_text(spec))
        assert data["data_format"] == spec.data_format


# ── Pyproject integrity ──────────────────────────────────────────────────


def test_pyproject_loads() -> None:
    """pyproject.toml is valid TOML."""
    import tomllib  # Python 3.11+

    raw = (REPO_ROOT / "pyproject.toml").read_text()
    data = tomllib.loads(raw)
    assert "tool" in data
    assert "pytest" in data["tool"]
    assert "ruff" in data["tool"]


# ── TUI app can be instantiated ──────────────────────────────────────────


@pytest.mark.asyncio
async def test_tui_app_creates() -> None:
    """FuxiTUI app object creates without crash."""
    from core.tui.app import FuxiTUI

    app = FuxiTUI()
    assert app is not None
