"""Stability scoring submodule for cluster parameter evaluation.

Contains ``_compute_stability`` (cross-seed ARI),
``_compute_cluster_coherence`` (per-cluster marker coherence), and
``_compute_splitting_gain`` (new clusters per resolution step).
"""

import logging
from typing import Literal

import numpy as np

logger = logging.getLogger(__name__)


def _compute_stability(
    adata,
    resolution,
    leiden_flavor: Literal["leidenalg", "igraph"] = "igraph",
    n_seeds=5,
    base_seed=42,
    device: str = "cpu",
    cfg=None,
    n_iterations=-1,
):
    """Compute cross-seed clustering stability via pairwise ARI.

    Re-runs Leiden clustering with *n_seeds* different random seeds at the
    given *resolution*, then computes pairwise adjusted Rand index between
    all label sets.  Returns mean ARI as a stability score.

    Parameters
    ----------
    adata : AnnData
    resolution : float
        Resolution parameter for Leiden clustering.
    leiden_flavor : str
        Flavour of leiden algorithm (default 'igraph').
    n_seeds : int
        Number of different seeds to run.
    base_seed : int
        Starting seed value.

    Returns
    -------
    float
        Mean pairwise ARI.  1.0 when *n_seeds* <= 1 or all clusterings
        are identical.
    """
    # Resolve n_iterations from config when available
    if cfg is not None:
        n_iterations = getattr(cfg.clustering, "stability_leiden_n_iterations", n_iterations)

    if n_seeds <= 1:
        return float("nan")

    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score

    seeds = range(base_seed, base_seed + n_seeds)
    temp_keys = [f"_temp_stab_{i}" for i in range(n_seeds)]

    label_sets, consecutive_aris = [], []
    for i, seed in enumerate(seeds):
        key = temp_keys[i]
        try:
            if device != "cpu":
                # Lazy import to avoid circular dependency at module load.
                from core.utils._gpu import gpu_leiden

                gpu_leiden(
                    adata,
                    resolution=resolution,
                    key_added=key,
                    random_state=seed,
                    # cuGraph doesn't accept flavor/directed/n_iterations;
                    # gpu_leiden strips them on the GPU path automatically.
                    flavor=leiden_flavor,
                    n_iterations=n_iterations,
                    directed=False,
                )
            else:
                sc.tl.leiden(
                    adata,
                    resolution=resolution,
                    key_added=key,
                    random_state=seed,
                    flavor=leiden_flavor,
                    n_iterations=n_iterations,
                    directed=False,
                )
            label_sets.append(adata.obs[key].values)
        except Exception as e:
            logger.warning(f"Stability seed {seed} failed: {e}")
        finally:
            if key in adata.obs.columns:
                del adata.obs[key]
        # Early termination: check if last 2 consecutive ARIs are 1.0
        if len(label_sets) >= 2:
            ari = adjusted_rand_score(label_sets[-2], label_sets[-1])
            consecutive_aris.append(ari)
            if len(consecutive_aris) >= 2 and all(a == 1.0 for a in consecutive_aris[-2:]):
                logger.info(
                    f"Early termination at seed {seed + 1}/{n_seeds}: "
                    "last 2 consecutive ARIs are 1.0"
                )
                break

    # Final cleanup: remove any remaining temp columns
    for col in list(adata.obs.columns):
        if col.startswith("_temp_stab_"):
            del adata.obs[col]

    n_runs = len(label_sets)
    if n_runs <= 1:
        return float("nan")

    # Compute pairwise ARI
    aris = []
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            ari = adjusted_rand_score(label_sets[i], label_sets[j])
            aris.append(ari)

    return float(np.mean(aris))


def _compute_cluster_coherence(
    adata,
    cluster_key,
    per_cell_scores,
    dominance_threshold=2.5,
    min_expression=0.05,
    min_cluster_size=10,
    log=None,
):
    """Per-cluster marker coherence metric that peaks at intermediate resolution.

    For each cluster label in adata.obs[cluster_key]:
      - Compute mean of per_cell_scores per cell type
      - Find best_score (max) and second_best_score
      - A cluster is "coherent" if:
        (a) best_score > min_expression (not noise)
        (b) best_score / max(second_best, 1e-10) > dominance_threshold
      - This ensures the cluster has ONE clearly dominant cell type

    coherence = n_coherent / n_valid_clusters
    where n_valid_clusters excludes clusters smaller than min_cluster_size.

    Parameters
    ----------
    adata : AnnData
    cluster_key : str
        Key in ``adata.obs`` for cluster labels.
    per_cell_scores : dict[str, np.ndarray]
        Cell-type → per-cell score array (shape: n_cells,).  Pre-computed
        by the caller (e.g. via ``sc.tl.score_genes()``).
    dominance_threshold : float
        Minimum ratio of best score to second-best for coherence.
    min_expression : float or str
        Minimum mean expression for a cluster to be considered (not noise).
        If a string starting with 'p' (e.g. 'p25'), interpreted as a percentile
        of per-cell scores for the best-matching cell type.
    min_cluster_size : int
        Minimum number of cells a cluster must have to be evaluated.
        Clusters below this threshold are skipped.
    log : Logger or None

    Returns
    -------
    float
        Fraction of clusters that are coherent (0.0-1.0).
        NaN if *per_cell_scores* is empty or has no valid entries.
    """
    if not per_cell_scores:
        return float("nan")

    # Check for any valid (non-None) entries
    has_valid = any(v is not None for v in per_cell_scores.values())
    if not has_valid:
        return float("nan")

    cluster_labels = adata.obs[cluster_key].values
    unique_clusters = np.unique(cluster_labels)
    n_clusters = len(unique_clusters)

    if n_clusters <= 1:
        return 1.0

    cell_types = list(per_cell_scores.keys())
    n_cells = len(cluster_labels)
    n_coherent = 0
    n_valid_clusters = 0
    _log = log or logger

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        cluster_size = int(np.sum(mask))

        # Skip clusters smaller than min_cluster_size
        if cluster_size < min_cluster_size:
            _log.debug(
                "Skipping cluster %s in coherence: size %d < min_cluster_size %d",
                cluster,
                cluster_size,
                min_cluster_size,
            )
            continue

        n_valid_clusters += 1
        mean_scores = []
        for ct in cell_types:
            scores = per_cell_scores[ct]
            if scores is not None and len(scores) == n_cells:
                mean_scores.append((ct, float(np.mean(scores[mask]))))

        if not mean_scores:
            continue

        # Sort by mean score descending
        mean_scores.sort(key=lambda x: x[1], reverse=True)
        top_score = mean_scores[0][1]
        best_ct = mean_scores[0][0]

        # Resolve threshold (float or percentile string)
        threshold = min_expression
        if isinstance(min_expression, str) and min_expression.startswith("p"):
            pct = int(min_expression[1:])
            if best_ct in per_cell_scores and per_cell_scores[best_ct] is not None:
                threshold = float(np.percentile(per_cell_scores[best_ct], pct))

        if top_score <= threshold:
            continue

        if len(mean_scores) < 2:
            return float("nan")
        else:
            second_score = mean_scores[1][1]
            if top_score / max(second_score, 1e-10) > dominance_threshold:
                n_coherent += 1

    if n_valid_clusters == 0:
        return float("nan")

    return n_coherent / n_valid_clusters


def _compute_splitting_gain(valid_by_resolution: list[dict]) -> dict[float, float]:
    """Compute splitting gain per resolution.

    Splitting gain measures how many new clusters are created per unit
    resolution increase. For each resolution r_i (except the lowest):
        splitting_gain(r_i) = max(0, (n_clusters(r_i) - n_clusters(r_{i-1})) / (r_i - r_{i-1}))

    Parameters
    ----------
    valid_by_resolution : list of dict
        List of results_summary entries sorted by resolution ascending,
        all from the SAME n_neighbors group.

    Returns
    -------
    dict[float, float]
        {resolution_value: splitting_gain_value}
        Empty dict if fewer than 2 resolutions.
    """
    if len(valid_by_resolution) < 2:
        return {}

    # Ensure sorted by resolution
    sorted_entries = sorted(valid_by_resolution, key=lambda e: e.get("resolution", 0.0))

    gains = {}
    for i in range(1, len(sorted_entries)):
        prev = sorted_entries[i - 1]
        curr = sorted_entries[i]
        r_prev = prev.get("resolution")
        r_curr = curr.get("resolution")
        k_prev = prev.get("n_clusters")
        k_curr = curr.get("n_clusters")

        if k_prev is None or k_curr is None or r_prev is None or r_curr is None:
            continue

        delta_k = k_curr - k_prev
        delta_r = r_curr - r_prev

        if delta_r > 0:
            gain = max(0.0, delta_k / delta_r)
        else:
            gain = 0.0

        gains[r_curr] = gain

    return gains
