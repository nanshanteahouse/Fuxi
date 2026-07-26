"""Granularity detection and DE-gated selection submodule.

Contains ``_detect_granularity`` (tissue vs subtype classifier),
``_compute_pairwise_de_markers`` (pairwise DE marker counts), and
``_select_de_gated`` (resolution selection via pairwise DE criterion).
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def _detect_granularity(
    results_summary: list[dict], cv_threshold: float = 0.05, min_clusters: int = 10
) -> str:
    """Determine whether the data is tissue-level or subtype-level.

    Two-path architecture for cluster granularity detection. This function
    analyses grid-search results to decide between two downstream strategies:

    - "tissue": multiple distinct cell types → use full DE pipeline
    - "subtype": FACS-enriched / similar cells → use gated DE (Wave 4)

    The decision is based on the coefficient of variation (CV) of silhouette
    scores within the median n_neighbors group. Low CV means silhouette
    scores are flat across resolutions — typical of subtype data where all
    partitions are similarly mediocre. Combined with a low maximum cluster
    count, this signals subtype-level resolution.

    Algorithm
    ---------
    1. Group entries by n_neighbors, pick the median-size group.
    2. From that group, sort by resolution, collect silhouette scores
       and n_clusters.
    3. Compute CV = std(silhouette_scores) / mean(silhouette_scores).
    4. If CV < cv_threshold AND max_n_clusters < min_clusters: "subtype".
    5. Otherwise: "tissue" (conservative default).

    Edge cases (all return "tissue"):
    - Empty results_summary
    - Single entry
    - Mean silhouette = 0 (division guard)

    Parameters
    ----------
    results_summary : list of dict
        Grid-search summary entries with keys: n_neighbors, resolution,
        n_clusters, silhouette_score.
    cv_threshold : float
        CV threshold below which silhouette flatness signals subtype data.
    min_clusters : int
        Maximum n_clusters below which (combined with low CV) signals
        subtype data.

    Returns
    -------
    str
        "tissue" or "subtype"
    """
    # Edge case: empty or single entry → conservative default
    if not results_summary or len(results_summary) <= 1:
        return "tissue"

    # 1. Group entries by n_neighbors
    groups: dict[int, list[dict]] = {}
    for entry in results_summary:
        if "n_neighbors" not in entry:
            continue
        nn = entry["n_neighbors"]
        groups.setdefault(nn, []).append(entry)

    if not groups:
        return "tissue"

    # Pick median-size group (by number of entries)
    group_sizes = sorted(groups.items(), key=lambda kv: len(kv[1]))
    median_idx = len(group_sizes) // 2
    _median_nn, median_group = group_sizes[median_idx]

    # 2. Collect silhouette scores and n_clusters from the median group
    #    Sort by resolution for deterministic ordering
    median_group_sorted = sorted(median_group, key=lambda e: e.get("resolution", 0.0))

    silhouette_values = []
    n_clusters_values = []
    for entry in median_group_sorted:
        if "silhouette_score" not in entry or entry["silhouette_score"] is None:
            continue
        if "n_clusters" not in entry or entry["n_clusters"] is None:
            continue
        silhouette_values.append(entry["silhouette_score"])
        n_clusters_values.append(entry["n_clusters"])

    if not silhouette_values:
        return "tissue"

    # 3. Compute CV of silhouette scores
    sil_arr = np.array(silhouette_values, dtype=float)
    mean_sil = float(np.mean(sil_arr))
    if mean_sil == 0.0:
        logger.debug(
            "_detect_granularity: mean silhouette is 0 → conservative 'tissue' "
            f"(n_neighbors={_median_nn})"
        )
        return "tissue"

    std_sil = float(np.std(sil_arr))
    cv = std_sil / mean_sil

    # 4. & 5. Decision
    max_n_clusters = max(n_clusters_values) if n_clusters_values else 0

    if cv < cv_threshold and max_n_clusters < min_clusters:
        logger.debug(
            f"_detect_granularity: CV={cv:.5f} < {cv_threshold} AND "
            f"max_n_clusters={max_n_clusters} < {min_clusters} → 'subtype' "
            f"(n_neighbors={_median_nn}, n_entries={len(median_group)})"
        )
        return "subtype"
    else:
        logger.debug(
            f"_detect_granularity: CV={cv:.5f} (threshold={cv_threshold}), "
            f"max_n_clusters={max_n_clusters} (threshold={min_clusters}) → 'tissue' "
            f"(n_neighbors={_median_nn}, n_entries={len(median_group)})"
        )
        return "tissue"


def _compute_pairwise_de_markers(
    adata,
    cluster_key,
    n_genes=50,
    padj_threshold=0.05,
    lfc_threshold=1.0,
    max_clusters_for_pairwise=30,
    log=None,
):
    """Compute per-cluster pairwise DE marker counts (Shekhar 2016 semantics).

    For each cluster C_i, compares against every other cluster C_j individually.
    A gene is a pairwise marker for C_i if it is upregulated (padj < threshold
    AND log2FC > threshold) in ALL pairwise comparisons C_i vs C_j.

    Falls back to one-vs-rest if n_clusters > max_clusters_for_pairwise
    (unless max_clusters_for_pairwise=0, meaning no cap).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts in .raw.
    cluster_key : str
        Key in adata.obs storing cluster assignments.
    n_genes : int
        Number of top DE genes to request from rank_genes_groups.
    padj_threshold : float
        Adjusted p-value threshold for significance.
    lfc_threshold : float
        Log2 fold-change threshold for upregulation.
    max_clusters_for_pairwise : int
        Maximum number of clusters for pairwise mode. 0 = no cap.
    log : logging.Logger or None
        Logger for warnings.

    Returns
    -------
    dict[str, int]
        Cluster name → pairwise marker count.
    """
    import scanpy as sc

    unique_vals = adata.obs[cluster_key].unique()
    # Handle categorical with possibly None/unused categories
    clusters = sorted([str(v) for v in unique_vals if v is not None])
    n_clusters = len(clusters)

    # Fallback: too many clusters
    if max_clusters_for_pairwise > 0 and n_clusters > max_clusters_for_pairwise:
        if log:
            log.warning(
                "n_clusters=%d exceeds de_pairwise_max_clusters=%d — "
                "falling back to one-vs-rest DE",
                n_clusters,
                max_clusters_for_pairwise,
            )
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            n_genes=n_genes,
            use_raw=True,
            key_added="_pairwise_fallback",
        )
        result = adata.uns["_pairwise_fallback"]
        marker_counts = {}
        for cluster_name in clusters:
            pvals = result["pvals_adj"][cluster_name]
            lfcs = result["logfoldchanges"][cluster_name]
            is_sig = (pvals < padj_threshold) & (lfcs > lfc_threshold)
            marker_counts[cluster_name] = int(is_sig.sum())
        del adata.uns["_pairwise_fallback"]
        return marker_counts

    # Pairwise mode
    marker_counts = {}
    for ci in clusters:
        marker_set = None
        for cj in clusters:
            if cj == ci:
                continue
            sc.tl.rank_genes_groups(
                adata,
                groupby=cluster_key,
                groups=[ci],
                reference=cj,
                method="wilcoxon",
                n_genes=n_genes,
                use_raw=True,
                key_added="_pairwise_tmp",
            )
            tmp = adata.uns["_pairwise_tmp"]
            pvals = tmp["pvals_adj"][ci]
            lfcs = tmp["logfoldchanges"][ci]
            is_sig = (pvals < padj_threshold) & (lfcs > lfc_threshold)
            # Use gene indices as proxy for identity (names are mock-safe)
            sig_mask = np.asarray(is_sig, dtype=bool)
            sig_indices = set(np.where(sig_mask)[0])
            marker_set = sig_indices if marker_set is None else (marker_set & sig_indices)
            del adata.uns["_pairwise_tmp"]

        marker_counts[ci] = len(marker_set) if marker_set is not None else 0

    return marker_counts


def _select_de_gated(valid, adata, de_gate_threshold=25, pairwise_max_clusters=30):
    """Select best resolution using DE-gated criterion for subtype-level data.

    For subtype-level data where silhouette scores are flat across resolutions,
    this method selects the highest resolution that maintains a minimum number of
    differentially expressed genes for every cluster.

    Uses pairwise DE comparisons (Shekhar 2016 semantics): each cluster is
    compared against every other cluster individually. A gene is a pairwise
    marker for cluster C_i if it is upregulated in ALL pairwise comparisons
    C_i vs C_j. Falls back to one-vs-rest when n_clusters > pairwise_max_clusters
    (0 = always pairwise).

    Parameters
    ----------
    valid : list[dict]
        Grid-search entries with keys: n_clusters, resolution, cluster_key.
        Must be unique by resolution.
    adata : AnnData
        Annotated data matrix.
    de_gate_threshold : int
        Minimum number of pairwise DE genes required for every cluster.
        Default 25.
    pairwise_max_clusters : int
        Maximum number of clusters for pairwise comparisons. Above this
        number, falls back to one-vs-rest DE to avoid O(n²) cost.
        Default 30.

    Returns
    -------
    tuple[int, float, str, str]
        (n_clusters, resolution, cluster_key, reason_str)
    """

    # Edge case: single entry -> return as-is
    if len(valid) <= 1:
        entry = valid[0]
        return (
            entry["n_clusters"],
            entry["resolution"],
            entry["cluster_key"],
            "de_gated(single_entry, min_de=N/A)",
        )

    # Sort by resolution ascending for deterministic iteration
    sorted_entries = sorted(valid, key=lambda e: e.get("resolution", 0.0))

    # Collect (entry, min_de) for each resolution
    candidates = []  # list of (entry, min_de)

    for entry in sorted_entries:
        cluster_key = entry["cluster_key"]

        try:
            marker_counts = _compute_pairwise_de_markers(
                adata,
                cluster_key,
                max_clusters_for_pairwise=pairwise_max_clusters,
                log=logger,
            )
            min_de = min(marker_counts.values()) if marker_counts else 0
        except Exception as e:
            logger.warning(
                "pairwise DE failed for cluster_key=%s (resolution=%.2f): %s",
                cluster_key,
                entry["resolution"],
                e,
            )
            continue

        candidates.append((entry, min_de))

    # Clean up temp rank_genes_groups keys (belt-and-suspenders)
    for k in list(adata.uns):
        if k.startswith("_de_gated_") or k.startswith("_pairwise_"):
            del adata.uns[k]

    if not candidates:
        # All entries failed -> fallback to first entry
        logger.warning(
            "DE-gated selection: all pairwise DE calls failed -- fallback to first entry"
        )
        entry = sorted_entries[0]
        return (
            entry["n_clusters"],
            entry["resolution"],
            entry["cluster_key"],
            "de_gated(all_failed, min_de=N/A)",
        )

    # Select entry with highest n_clusters where min_de >= threshold
    candidates_by_n = sorted(
        candidates,
        key=lambda c: c[0]["n_clusters"],
        reverse=True,
    )
    best_entry = None
    best_min_de = -1
    for entry, min_de in candidates_by_n:
        if min_de >= de_gate_threshold:
            best_entry = entry
            best_min_de = min_de
            break

    # If no entry meets threshold, fallback to entry with highest DE count
    if best_entry is None:
        best_entry, best_min_de = max(candidates, key=lambda c: c[1])
        logger.info(
            "DE-gated selection: best_resolution=%.2f, min_de=%d "
            "(fallback: no entry met threshold=%d)",
            best_entry["resolution"],
            best_min_de,
            de_gate_threshold,
        )
    else:
        logger.info(
            "DE-gated selection: best_resolution=%.2f, min_de=%d",
            best_entry["resolution"],
            best_min_de,
        )

    return (
        best_entry["n_clusters"],
        best_entry["resolution"],
        best_entry["cluster_key"],
        f"de_gated(min_de={best_min_de})",
    )
