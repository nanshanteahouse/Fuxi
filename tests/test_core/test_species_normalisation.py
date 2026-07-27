"""Tests for P0–P3 species normalisation fixes.

Three test layers:
  1. ``resolve_config`` — YAML-driven config loading (tests 1–3)
  2. ``Config`` model_validator — direct Pydantic construction (tests 4–5)
  3. ``load_cell_cycle_genes`` alias — ``macaque`` → ``macaca`` (test 6)
  4. ``_normalise_species`` — 6 parametrised alias mappings (test 7)

P0–P3 gaps documented by failing parametrised cases:
  - ``"macaca mulatta"`` (lowercase with space, not underscored)
  - ``"cynomolgus"`` (accepted pipeline key but no normalisation entry)
  - ``"cynomolgus_macaque"`` (no normalisation entry at all)

These are *expected* to fail until the corresponding entries are added
to ``_SPECIES_NORMALISE`` in ``core/preprocess/format_detector.py``.
"""

from pathlib import Path

import pytest
import yaml

from core.config.schema import Config
from core.kb.cell_cycle import load_cell_cycle_genes
from core.preprocess.format_detector import _normalise_species
from core.utils._config import resolve_config

# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def cfg_path(tmp_path: Path) -> Path:
    """Write a minimal config YAML with ``species: macaca_fascicularis``."""
    data = {
        "species": "macaca_fascicularis",
        "project_dir": str(tmp_path),
    }
    path = tmp_path / "config.yaml"
    with open(path, "w") as f:
        yaml.dump(data, f)
    return path


@pytest.fixture
def cfg_path_macaque(tmp_path: Path) -> Path:
    """Write a minimal config YAML with ``species: macaque``."""
    data = {
        "species": "macaque",
        "project_dir": str(tmp_path),
    }
    path = tmp_path / "config.yaml"
    with open(path, "w") as f:
        yaml.dump(data, f)
    return path


# ═══════════════════════════════════════════════════════════════════════
# 1) resolve_config canonicalisation
# ═══════════════════════════════════════════════════════════════════════


def test_resolve_config_canonicalises_macaca_fascicularis(cfg_path: Path):
    """``resolve_config`` normalises ``macaca_fascicularis`` → ``macaque``."""
    cfg = resolve_config(str(cfg_path))
    assert cfg.species == "macaque"


def test_resolve_config_syncs_grn_species(cfg_path_macaque: Path):
    """GRN species default ``"human"`` is overridden to ``cfg.species``."""
    cfg = resolve_config(str(cfg_path_macaque))
    # The sync logic in resolve_config writes cfg.species into grn.species
    # when grn.species is still the default "human" and cfg.species ≠ "human".
    assert cfg.grn.species == "macaque"


def test_resolve_config_syncs_enrichment_organism(cfg_path_macaque: Path):
    """Enrichment organism default ``"human"`` is overridden to ``cfg.species``."""
    cfg = resolve_config(str(cfg_path_macaque))
    assert cfg.enrichment.organism == "macaque"


# ═══════════════════════════════════════════════════════════════════════
# 2) Config model_validator
# ═══════════════════════════════════════════════════════════════════════


def test_model_validator_normalises_on_construct():
    """``Config(species="macaca_fascicularis")`` → ``cfg.species == "macaque"``."""
    cfg = Config(species="macaca_fascicularis")
    assert cfg.species == "macaque"


def test_model_validator_idempotent():
    """Already-canonical ``"macaque"`` stays unchanged."""
    cfg = Config(species="macaque")
    assert cfg.species == "macaque"


# ═══════════════════════════════════════════════════════════════════════
# 3) cell_cycle alias: macaque → macaca
# ═══════════════════════════════════════════════════════════════════════


def test_cell_cycle_alias_macaque_to_macaca():
    """``load_cell_cycle_genes("macaque")`` redirects to the ``macaca`` module."""
    s_genes, g2m_genes = load_cell_cycle_genes("macaque")
    assert isinstance(s_genes, list)
    assert isinstance(g2m_genes, list)
    assert len(s_genes) > 0, "S-phase gene list should be non-empty"
    assert len(g2m_genes) > 0, "G2/M-phase gene list should be non-empty"


# ═══════════════════════════════════════════════════════════════════════
# 4) Parametrised alias mappings → macaque
# ═══════════════════════════════════════════════════════════════════════
#
# P0–P3 gaps (require adding entries to _SPECIES_NORMALISE):
#   "macaca mulatta"     — lowercase space form (missing)
#   "cynomolgus"         — not yet in normalisation dict
#   "cynomolgus_macaque" — not yet in normalisation dict
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "alias",
    [
        pytest.param("macaca_fascicularis"),  # P0 — underscored form
        pytest.param("Macaca fascicularis"),  # P0 — Latin binomial
        pytest.param(
            "macaca mulatta",  # P1 — lowercase space
            marks=pytest.mark.xfail(
                strict=False,
                reason="'macaca mulatta' not yet in _SPECIES_NORMALISE",
            ),
        ),
        pytest.param("Macaca_mulatta"),  # P1 — titlecase underscore
        pytest.param(
            "cynomolgus",  # P2 — common alias
            marks=pytest.mark.xfail(
                strict=False,
                reason="'cynomolgus' not yet in _SPECIES_NORMALISE",
            ),
        ),
        pytest.param(
            "cynomolgus_macaque",  # P3 — compound alias
            marks=pytest.mark.xfail(
                strict=False,
                reason="'cynomolgus_macaque' not yet in _SPECIES_NORMALISE",
            ),
        ),
    ],
)
def test_parametrised_aliases(alias: str):
    """All non-canonical species strings normalise to ``"macaque"``."""
    result = _normalise_species(alias)
    assert result == "macaque", f"{alias!r} → {result!r}, expected 'macaque'"
