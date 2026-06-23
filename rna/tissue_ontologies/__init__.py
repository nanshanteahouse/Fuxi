"""tissue_ontologies — Knowledge Base loaders for supported tissues.

Usage::

    from tissue_ontologies import load_kb
    kb = load_kb("retina")
"""


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
        from tissue_ontologies.retina import retina_expert_kb
        return retina_expert_kb

    raise ValueError(
        f"Unsupported tissue KB: '{tissue_name}'. "
        f"Available: retina"
    )
