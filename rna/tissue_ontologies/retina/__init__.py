"""Retina cell-type Knowledge Base.

Import the unified retina KB directly::

    from tissue_ontologies.retina import retina_expert_kb

The KB is built at import time by merging 7 source publications.
Validation::

    from tissue_ontologies.retina.validate import validate_kb
    is_ok, errors = validate_kb(retina_expert_kb)
"""

from .merge import retina_expert_kb

__all__ = ["retina_expert_kb"]
