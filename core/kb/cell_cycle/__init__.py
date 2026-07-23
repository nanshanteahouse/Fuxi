"""Cell cycle gene lists per species for sc.tl.score_genes_cell_cycle.

Each species module (``human``, ``mouse``, ``macaca``, ``zebrafish``) exports
``S_GENES`` and ``G2M_GENES`` as Python lists of gene symbols.

Usage::

    from core.kb.cell_cycle import load_cell_cycle_genes
    s_genes, g2m_genes = load_cell_cycle_genes("human")
"""

import importlib
import logging
from typing import List, Tuple

import yaml

logger = logging.getLogger(__name__)

_AVAILABLE_SPECIES = {"human", "mouse", "macaca", "zebrafish"}


def load_cell_cycle_genes(species: str) -> Tuple[List[str], List[str]]:
    """Load S-phase and G2M-phase gene lists for a species.

    Parameters
    ----------
    species : str
        Species identifier (e.g. ``"human"``, ``"mouse"``).

    Returns
    -------
    tuple of (list[str], list[str])
        ``(s_genes, g2m_genes)`` — lists of gene symbols for S phase
        and G2/M phase respectively.

    Raises
    ------
    ValueError
        If *species* is not supported.
    """
    species = species.lower()
    if species not in _AVAILABLE_SPECIES:
        raise ValueError(
            f"Unsupported species '{species}' for cell-cycle scoring. "
            f"Supported: {', '.join(sorted(_AVAILABLE_SPECIES))}. "
            "Set cfg.normalization.score_cell_cycle=False or add species YAML."
        )
    mod = importlib.import_module(f".{species}", "core.kb.cell_cycle")
    return (mod.S_GENES, mod.G2M_GENES)


def _load_yaml_species(species: str) -> Tuple[List[str], List[str]]:
    """Load cell-cycle gene lists from a YAML file for a given species.

    This is used internally by each species module (``human.py``,
    ``mouse.py``, etc.) to load lists from the corresponding ``.yaml`` file.
    """
    import os

    _dir = os.path.dirname(os.path.abspath(__file__))
    yaml_path = os.path.join(_dir, f"{species}.yaml")
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)
    s_genes: List[str] = data.get("S_GENES", [])
    g2m_genes: List[str] = data.get("G2M_GENES", [])
    return (s_genes, g2m_genes)
