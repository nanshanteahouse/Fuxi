"""Retina cell-type Knowledge Base.

Import the unified retina KB directly::

    from core.kb.retina import retina_expert_kb

The KB is built at import time by merging curated source publications.

Validation::

    from core.kb.validate import validate_kb
    is_ok, errors = validate_kb(retina_expert_kb)
"""

from ..merge import build_tissue_kb
from .config import TYPE_ALIASES, SOURCES_DIR, HIERARCHY_PATH

# ── Build KB once at import time ────────────────────────────────────
retina_expert_kb = build_tissue_kb(
    SOURCES_DIR,
    type_aliases=TYPE_ALIASES,
    hierarchy_yaml_path=HIERARCHY_PATH,
)

__all__ = ["retina_expert_kb"]
