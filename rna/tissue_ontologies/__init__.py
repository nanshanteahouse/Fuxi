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
