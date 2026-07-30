"""Grid-search result enrichment submodule.

Enriches each entry in a grid-search results summary with multi-metric scores
(stability, cluster coherence, kb annotatable rate, splitting gain).
"""

import logging

from .stability import _compute_cluster_coherence, _compute_splitting_gain, _compute_stability

logger = logging.getLogger(__name__)

DEFAULT_MULTI_METRIC_WEIGHTS = {
    "silhouette": 0.15,
    "stability": 0.20,
    "cluster_coherence": 0.45,
    "splitting_gain": 0.15,
    "kb_annotatable_rate": 0.05,
}


def enrich_grid_results(
    adata,
    results_summary,
    cfg,
    compute_per_cell_scores_fn=None,
    log=None,
    use_rep=None,
):
    """Enrich *results_summary* entries with multi-metric scores.

    For each (n_neighbors, resolution) entry, computes:
      - stability_score  (via cross-seed ARI)
      - cluster_coherence (via marker gene per-cell scores)
      - kb_annotatable_rate (via tissue knowledge base, if *cfg.tissue_kb* is set)
      - splitting_gain (per n_neighbors group)

    The list is modified **in place** and returned. The caller is
    responsible for persistence (e.g. via ``param_grid_summary.csv``).

    Parameters
    ----------
    adata : AnnData
    results_summary : list[dict]
        Each dict has keys: *n_neighbors*, *resolution*, *cluster_key*.
        Enriched with *stability_score*, *cluster_coherence*, etc. IN PLACE.
    cfg : Config
    compute_per_cell_scores_fn : callable or None
        ``(adata, cfg) -> dict[str, np.ndarray]``
        Pre-computed per-cell scores for each cell type.
        When ``None`` and ``cfg.marker.marker_dict`` has entries, the
        default marker-based logic (``sc.tl.score_genes``) is used.
    log : Logger or None
    use_rep : str or None
        Key in ``adata.obsm`` for the embedding (e.g. ``'X_pca'``).
        Required when the KNN graph needs rebuilding.

    Returns
    -------
    list[dict]
        The same *results_summary* (modified in place).
    """
    import logging as _logging

    import scanpy as sc

    if log is None:
        log = _logging.getLogger(__name__)

    marker_dict = getattr(cfg.marker, "marker_dict", None) or {}
    has_markers = bool(marker_dict)
    n_stab_seeds = getattr(cfg.clustering, "stability_n_seeds", 12)
    dominance_threshold = getattr(cfg.clustering, "multi_metric_coverage_ratio_threshold", 2.5)
    leiden_flavor = getattr(cfg.clustering, "leiden_flavor", "igraph")
    device = getattr(cfg.execution, "device", "cpu")

    # Group results by n_neighbors
    by_n = {}
    for r in results_summary:
        n_ = r.get("n_neighbors")
        by_n.setdefault(n_, []).append(r)

    for n_val, group in by_n.items():
        # --- Check if KNN rebuild is needed (optimisation) ---
        _need_rebuild_knn = True
        try:
            _nb_params = adata.uns.get("neighbors", {}).get("params", {}) or {}
            if (
                _nb_params.get("n_neighbors") == n_val
                and use_rep is not None
                and _nb_params.get("use_rep") == use_rep
                and "connectivities" in adata.obsp
            ):
                _need_rebuild_knn = False
        except Exception:
            pass

        if _need_rebuild_knn:
            if use_rep is None:
                log.warning(
                    "KNN rebuild needed but use_rep is None -- skipping n_neighbors=%d",
                    n_val,
                )
                continue
            try:
                # Ensure obsm data is CPU-resident (grid search may leave GPU arrays)
                _rep = adata.obsm[use_rep]
                if hasattr(_rep, "get"):  # cupy array → numpy
                    import numpy as np

                    adata.obsm[use_rep] = _rep.get()  # cupy requires explicit .get()
                sc.pp.neighbors(
                    adata,
                    n_neighbors=n_val,
                    n_pcs=min(cfg.pca.n_pcs_use, adata.obsm[use_rep].shape[1]),
                    use_rep=use_rep,
                    random_state=cfg.execution.random_seed,
                )
            except Exception as e:
                log.warning(
                    "KNN rebuild failed for n_neighbors=%d: %s -- skipping group",
                    n_val,
                    e,
                )
                continue

        # --- Pre-compute per-cell scores ---
        per_cell_scores = {}
        if compute_per_cell_scores_fn is not None:
            try:
                scores = compute_per_cell_scores_fn(adata, cfg)
                if scores is not None:
                    per_cell_scores = scores
            except Exception as e:
                log.warning(
                    "per_cell_scores computation failed: %s -- falling back to no markers",
                    e,
                )
        elif has_markers and adata.raw is not None:
            # Default marker_dict-based scoring (RNA and Spatial common case)
            from anndata import utils as anndata_utils

            adata.raw._var.index = anndata_utils.make_index_unique(adata.raw._var.index, join="-")
            try:
                for ct, genes in marker_dict.items():
                    valid_genes = [g for g in genes if g in adata.raw.var_names]
                    if valid_genes:
                        sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f"_score_{ct}")
                        per_cell_scores[ct] = adata.obs[f"_score_{ct}"].values.copy()
                # Clean up temporary score columns
                for col in list(adata.obs.columns):
                    if col.startswith("_score_") and col in adata.obs.columns:
                        adata.obs.drop(columns=[col], inplace=True)
            except Exception as e:
                log.warning(
                    "Marker score pre-computation failed: %s -- falling back to no markers",
                    e,
                )
                per_cell_scores = {}
        elif has_markers and adata.raw is None:
            log.warning(
                "adata.raw is None -- cannot compute marker coverage. Degrading to silhouette+stability only."
            )
            has_markers = False

        # --- Compute stability + marker coverage for each entry ---
        for entry in group:
            try:
                resolution = entry["resolution"]
                ck = entry["cluster_key"]

                entry["stability_score"] = _compute_stability(
                    adata,
                    resolution=resolution,
                    leiden_flavor=leiden_flavor,
                    n_seeds=n_stab_seeds,
                    device=device,
                    cfg=cfg,
                )

                if per_cell_scores:
                    entry["cluster_coherence"] = _compute_cluster_coherence(
                        adata,
                        ck,
                        per_cell_scores,
                        dominance_threshold=dominance_threshold,
                    )

                # --- KB annotatable rate ---
                if getattr(cfg, "tissue_kb", None) and per_cell_scores:
                    labels = adata.obs[ck].values
                    unique_clusters = np.unique(labels)
                    n_total = len(unique_clusters)
                    n_annotatable = 0
                    for cl in unique_clusters:
                        mask_ = labels == cl
                        best_score = 0.0
                        for ct_ in per_cell_scores:
                            scores = per_cell_scores[ct_]
                            if scores is not None and len(scores) == len(labels):
                                mean_val = float(np.mean(scores[mask_]))
                                if mean_val > best_score:
                                    best_score = mean_val
                        if best_score > 0.5:
                            n_annotatable += 1
                    rate = n_annotatable / n_total if n_total > 0 else 0.0
                    entry["kb_annotatable_rate"] = rate
                    log.info("KB annotatable rate: %.3f", rate)

            except Exception as e:
                log.warning(
                    "Enrichment failed for n_neighbors=%d, resolution=%.1f: %s",
                    entry.get("n_neighbors"),
                    entry.get("resolution"),
                    e,
                )
                entry["stability_score"] = None
                entry["cluster_coherence"] = None
                entry["kb_annotatable_rate"] = None

        # --- Compute splitting_gain for this n_neighbors group ---
        if len(group) >= 2:
            group_sorted = sorted(group, key=lambda e: e.get("resolution", 0.0))
            gains = _compute_splitting_gain(group_sorted)
            for entry in group:
                entry["splitting_gain"] = gains.get(entry["resolution"], 0.0)

    return results_summary
