"""KB — Knowledge Base loaders for supported tissues.

Usage::

    from core.kb import load_kb
    kb = load_kb("retina")
"""

import importlib
import logging
import os

import pandas as pd
import yaml

_log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════
#  Tissue discovery
# ═══════════════════════════════════════════════════════════════════════


def _discover_tissues() -> dict[str, str]:
    """Scan subdirectories for tissue modules."""
    _dir = os.path.dirname(os.path.abspath(__file__))
    tissues: dict[str, str] = {}
    for entry in sorted(os.listdir(_dir)):
        path = os.path.join(_dir, entry)
        if not os.path.isdir(path):
            continue
        if entry.startswith("_") or entry.startswith("."):
            continue
        init = os.path.join(path, "__init__.py")
        if os.path.isfile(init):
            tissues[entry] = entry
    return tissues


_AVAILABLE_TISSUES = _discover_tissues()

# ═══════════════════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════════════════


def load_kb(tissue_name: str):
    """Load the Knowledge Base for a given tissue.

    Parameters
    ----------
    tissue_name : str
        Tissue identifier (e.g. ``"retina"``).

    Returns
    -------
    dict
        KB dict consumable by ``utils.marker_scoring`` and
        ``utils.evidence_fusion``.

    Raises
    ------
    ValueError
        If the tissue name is not supported.
    """
    if tissue_name not in _AVAILABLE_TISSUES:
        raise ValueError(
            f"Unsupported tissue KB: '{tissue_name}'. "
            f"Available: {', '.join(sorted(_AVAILABLE_TISSUES))}"
        )
    mod = importlib.import_module(f".{tissue_name}", __package__)
    kb_attr = f"{tissue_name}_expert_kb"
    return getattr(mod, kb_attr)


def load_all_kb_markers(tissue_name: str) -> set[str]:
    """Extract a flat set of all marker gene symbols from a tissue KB.

    Parameters
    ----------
    tissue_name : str
        Tissue identifier (e.g. ``"retina"``).

    Returns
    -------
    set[str]
        Flat set of all marker gene symbols (UPPERCASE).
        Returns empty set if tissue is unsupported or KB module is missing.
    """
    try:
        kb = load_kb(tissue_name)
    except (ValueError, ImportError):
        _log.warning("Could not load KB markers for %s", tissue_name)
        return set()
    markers: set[str] = set()
    for entry in kb.values():
        if isinstance(entry, dict) and "markers" in entry:
            for key in ("confirm", "add", "refine"):
                markers.update(entry["markers"].get(key, {}).keys())
    return {g.upper() for g in markers}


def load_adjacency(tissue_name: str) -> pd.DataFrame:
    """Load the anatomical adjacency matrix for a given tissue.

    Parameters
    ----------
    tissue_name : str
        Tissue identifier (e.g. ``"retina"``).

    Returns
    -------
    pd.DataFrame
        Adjacency table with columns ``["source", "target", "adjacency_type"]``.
        Returns an empty DataFrame (same columns) if the tissue has no
        adjacency module or is not supported.
    """
    try:
        adj_mod = importlib.import_module(f".{tissue_name}.adjacency", __package__)
        return pd.DataFrame(
            adj_mod.ADJACENCY,
            columns=["source", "target", "adjacency_type"],
        )
    except ImportError:
        pass
    return pd.DataFrame(columns=["source", "target", "adjacency_type"])


def load_pathway_relevance(tissue_name: str) -> dict:
    """Load tissue-specific pathway relevance metadata.

    Parameters
    ----------
    tissue_name : str
        Tissue identifier (e.g. ``"retina"``).

    Returns
    -------
    dict
        Dict with keys:
        - key_pathways (list[str])
        - generic_pathways (list[str])
        - kb_pathway_markers (dict[str, list[str]])
        Returns empty dict for unsupported tissues.
    """
    try:
        pr_mod = importlib.import_module(f".{tissue_name}.pathway_relevance", __package__)
        prefix = tissue_name.upper()
        return {
            "key_pathways": list(getattr(pr_mod, f"{prefix}_KEY_PATHWAYS", [])),
            "generic_pathways": list(getattr(pr_mod, f"{prefix}_GENERIC_PATHWAYS", [])),
            "kb_pathway_markers": dict(getattr(pr_mod, f"{prefix}_KB_PATHWAY_MARKERS", {})),
        }
    except ImportError:
        pass
    return {}


def load_synonyms(tissue_name: str) -> dict:
    """Load cell-type synonyms for a given tissue.

    Parameters
    ----------
    tissue_name : str
        Tissue identifier (e.g. ``"retina"``).

    Returns
    -------
    dict
        ``{canonical_key: {"display_name": str, "synonyms": list[str]}}``
        or empty dict if unavailable.
    """
    try:
        cfg_mod = importlib.import_module(f".{tissue_name}.config", __package__)
        syn_path = getattr(cfg_mod, "SYNONYMS_PATH", None)
        if syn_path and os.path.isfile(syn_path):
            with open(syn_path) as fh:
                return yaml.safe_load(fh)
    except ImportError:
        pass
    return {}
