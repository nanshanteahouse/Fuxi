"""Cluster parameter selection methods for the Fuxi pipeline.

Provides objective, quantitative selection of the best (n_neighbors, resolution)
from a grid search summary, supporting four methods: 'pareto_elbow',
'silhouette', 'multi_metric' (default), and None (manual).

Exports:
    select_best_params(results_summary, method, best_resolution=None, best_n_neighbors=0, multi_metric_weights=None)
        -> (best_n, best_r, method_label, reason_str)
    _compute_stability(adata, resolution, ...) -> float
    _compute_cluster_coherence(adata, cluster_key, per_cell_scores, ...) -> float
"""

# Re-export from submodules — preserves the ``from core.cluster.evaluation import X`` pattern.
from .enrichment import DEFAULT_MULTI_METRIC_WEIGHTS, enrich_grid_results
from .granularity import _compute_pairwise_de_markers, _detect_granularity, _select_de_gated
from .selection import _select_multi_metric, select_best_params, select_best_umap_params
from .stability import _compute_cluster_coherence, _compute_splitting_gain, _compute_stability

__all__ = [
    "select_best_params",
    "select_best_umap_params",
    "enrich_grid_results",
    "_detect_granularity",
    "_select_de_gated",
    "_compute_stability",
    "_compute_cluster_coherence",
    "_compute_splitting_gain",
    "_compute_pairwise_de_markers",
    "_select_multi_metric",
    "DEFAULT_MULTI_METRIC_WEIGHTS",
]
