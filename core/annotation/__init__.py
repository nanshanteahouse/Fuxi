"""core.annotation — Shared cell-type annotation engine.

Provides unified annotation, name standardization, and marker-based scoring
across RNA, ATAC, and Spatial modalities.
"""

from core.annotation.engine import run_unified_annotation
from core.annotation.standardizer import StandardOntology, map_annotations
from core.annotation.scoring import _normalize_gene_name

__all__ = [
    "run_unified_annotation",
    "StandardOntology",
    "map_annotations",
    "_normalize_gene_name",
]
