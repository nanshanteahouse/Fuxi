"""Smoke tests — fast, no-heavy-dependency checks that gate CI and pre-push.

All tests here MUST:
  - Complete in < 3 seconds each (total suite < 15s)
  - Not require network, LLM, or real data (FUXI_DATA_ROOT optional)
  - Not depend on other test files

Tag: @pytest.mark.smoke
Run:  python -m pytest tests/ -m smoke -x --tb=short
"""

import subprocess
import sys
from pathlib import Path

import pytest

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


def test_pipeline_list_rna() -> None:
    """run_pipeline --modality rna --list exits cleanly."""
    result = subprocess.run(
        [sys.executable, "-m", "core.run_pipeline", "--modality", "rna", "--list"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        timeout=15,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert "00_load" in result.stdout


def test_pipeline_list_atac() -> None:
    """run_pipeline --modality atac --list exits cleanly."""
    result = subprocess.run(
        [sys.executable, "-m", "core.run_pipeline", "--modality", "atac", "--list"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        timeout=15,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"


def test_pipeline_list_spatial() -> None:
    """run_pipeline --modality spatial --list exits cleanly."""
    result = subprocess.run(
        [sys.executable, "-m", "core.run_pipeline", "--modality", "spatial", "--list"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        timeout=15,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"


# ── Config template discovery ────────────────────────────────────────────


def test_config_templates_exist() -> None:
    """All config templates are present on disk."""
    template_dir = REPO_ROOT / "templates" / "config_templates"
    templates = list(template_dir.glob("config_*.yaml"))
    assert len(templates) >= 4, f"Expected ≥4 templates, found {len(templates)}"


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
