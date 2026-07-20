"""core.annotation — Shared cell-type annotation engine.

Provides unified annotation, name standardization, and marker-based scoring
across RNA, ATAC, and Spatial modalities.
"""

from core.annotation.engine import run_unified_annotation
from core.annotation.scoring import _normalize_gene_name
from core.annotation.standardizer import StandardOntology

__all__ = [
    "run_unified_annotation",
    "StandardOntology",
    "_normalize_gene_name",
]
