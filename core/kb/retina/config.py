"""
retina/config.py — Retina-specific configuration for the tissue ontology engine.

Provides constants consumed by :mod:`core.kb.merge` to build
the unified retina expert knowledge base.
"""

import os

# ── Type-key synonym table for cross-source merging ──────────────────
# Maps source-internal type keys to canonical KB names.
TYPE_ALIASES: dict[str, str] = {
    "Retinal_Ganglion_Cell": "RGC",
}

# ── Paths (relative to this file's directory) ────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))

SOURCES_DIR: str = os.path.join(_HERE, "sources")
HIERARCHY_PATH: str = os.path.join(_HERE, "hierarchy.yaml")
SYNONYMS_PATH: str = os.path.join(_HERE, "synonyms.yaml")
ADJACENCY_PATH: str = os.path.join(_HERE, "adjacency.py")
PATHWAY_PATH: str = os.path.join(_HERE, "pathway_relevance.py")
