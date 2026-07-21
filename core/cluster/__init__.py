"""core.cluster — Shared clustering infrastructure across modalities.

Provides modality-agnostic grid-search orchestration and parameter evaluation.
"""

from core.cluster.evaluation import (
    _compute_cluster_coherence,
    _compute_splitting_gain,
    _compute_stability,
    _detect_granularity,
    _select_de_gated,
    select_best_params,
    select_best_umap_params,
)

__all__ = [
    "select_best_params",
    "_compute_stability",
    "_compute_cluster_coherence",
    "_compute_splitting_gain",
    "_detect_granularity",
    "_select_de_gated",
    "select_best_umap_params",
]
