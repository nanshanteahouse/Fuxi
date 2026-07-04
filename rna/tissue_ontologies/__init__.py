"""tissue_ontologies — Knowledge Base loaders for supported tissues.

Usage::

    from rna.tissue_ontologies import load_kb
    kb = load_kb("retina")
"""
import pandas as pd



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
    if tissue_name == "retina":
        from .retina import retina_expert_kb
        return retina_expert_kb

    raise ValueError(
        f"Unsupported tissue KB: '{tissue_name}'. "
        f"Available: retina"
    )

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
        if tissue_name == "retina":
            from .retina.adjacency import ADJACENCY
            return pd.DataFrame(
                ADJACENCY,
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
        - key_pathways (list[str]): 该组织关键通路白名单
        - generic_pathways (list[str]): 在该组织场景下通用的通路黑名单
        - kb_pathway_markers (dict[str, list[str]]): 通路→标记基因映射
    Returns empty dict (never crashes) for unsupported tissues.
    """
    try:
        if tissue_name == "retina":
            from .retina.pathway_relevance import (
                RETINA_KEY_PATHWAYS, RETINA_GENERIC_PATHWAYS,
                RETINA_KB_PATHWAY_MARKERS,
            )
            return {
                "key_pathways": list(RETINA_KEY_PATHWAYS),
                "generic_pathways": list(RETINA_GENERIC_PATHWAYS),
                "kb_pathway_markers": dict(RETINA_KB_PATHWAY_MARKERS),
            }
    except ImportError:
        pass
    return {}
