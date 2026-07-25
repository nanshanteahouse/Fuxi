# ruff: noqa: N806
#!/usr/bin/env python3
"""
Generate synthetic h5ad fixtures for unit testing.

Produces 4 files in tests/fixtures/:
    synthetic_rna.h5ad                5k cells, well-separated (10 blobs)
    synthetic_overlapping.h5ad        5k cells, moderate overlap
    synthetic_severely_overlapping.h5ad  5k cells, severe overlap
    synthetic_100k.h5ad              100k cells, well-separated (funnel tests)

Each h5ad contains:
  - X:          Poisson count matrix (n_cells × n_genes)
  - obsm[X_pca]:  PCA-like embedding (30-dim)
  - obs[batch]:   3 batch labels (categorical)
  - obs[true_cell_type]:  ground-truth cluster labels
  - var[highly_variable]:  True for top 500 genes
  - uns[marker_dict]:  dict[cell_type_name → list[marker_gene_names]]

Usage:
    python tests/fixtures/generate_synthetic.py

Dependencies:  scanpy, numpy, scikit-learn (already in project)
"""

from __future__ import annotations

import os

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.datasets import make_blobs
from sklearn.metrics.pairwise import euclidean_distances

# ── Paths ────────────────────────────────────────────────────────────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_FIXTURES_DIR = os.path.join(_REPO_ROOT, "tests", "fixtures")

# ── Constants ────────────────────────────────────────────────────────
RANDOM_STATE = 42
N_GENES = 2000
N_PCS = 30  # PCA dimensions
N_CLUSTERS = 10
N_BATCHES = 3
N_HVG = 500  # number of highly-variable genes
N_MARKERS_PER_TYPE = 5
N_HVG_POOL = 200  # select markers from top N HVG


def _ensure_min_distance(
    centroids: np.ndarray,
    target_min_dist: float,
    rng: np.random.RandomState,
) -> np.ndarray:
    """Scale centroids so that the minimum pairwise distance >= target_min_dist.

    Args:
        centroids:  (n_centers, n_features) array of cluster centers.
        target_min_dist:  Desired minimum Euclidean distance between any two
            centroids (in standard-deviation units — i.e. meaningfully
            comparable to ``cluster_std``).
        rng:  NumPy RandomState for reproducibility.

    Returns:
        Scaled centroids with the required minimum separation.
    """
    dists = euclidean_distances(centroids)
    np.fill_diagonal(dists, np.inf)
    current_min = dists.min()
    if current_min < target_min_dist * 1.005:
        scale = target_min_dist / current_min * 1.005  # small margin
        centroids = centroids * scale
    return centroids


def _random_orthogonal(d: int, rng: np.random.RandomState) -> np.ndarray:
    """Return a random d×d orthogonal matrix (rotation/reflection)."""
    Q, _ = np.linalg.qr(rng.randn(d, d))
    return Q


def generate_dataset(
    n_cells_per_blob: int = 500,
    cluster_std: float = 1.0,
    min_centroid_distance: float = 5.0,
    seed: int = RANDOM_STATE,
    label_prefix: str = "",
) -> ad.AnnData:
    """Generate a synthetic single-cell dataset as an AnnData object.

    Parameters
    ----------
    n_cells_per_blob:
        Number of cells per Gaussian blob.
    cluster_std:
        Standard deviation of each Gaussian blob.
    min_centroid_distance:
        Minimum Euclidean distance enforced between any two blob centroids.
        Together with ``cluster_std`` this controls the degree of cluster
        separation / overlap.
    seed:
        Random seed for reproducibility.
    label_prefix:
        Optional prefix for cell barcodes (e.g. ``"well_"``).

    Returns
    -------
    AnnData object with ``X``, ``obsm["X_pca"]``, ``obs["batch"]``,
    ``obs["true_cell_type"]``, ``var["highly_variable"]``,
    ``uns["marker_dict"]``.
    """
    rng = np.random.RandomState(seed)
    n_cells = N_CLUSTERS * n_cells_per_blob

    # ── 1. Generate cluster centroids via make_blobs ──────────────────
    # Use make_blobs to produce the true cluster centers (requirement).
    # We pass n_samples=N_CLUSTERS so we get one sample per cluster;
    # the true internal centers are returned via return_centers=True.
    _, _, centroids = make_blobs(
        n_samples=N_CLUSTERS,
        n_features=N_PCS,
        centers=N_CLUSTERS,
        cluster_std=1.0,
        center_box=(-10, 10),
        return_centers=True,
        random_state=rng.randint(0, 2**31),
    )
    centroids = _ensure_min_distance(centroids, min_centroid_distance, rng)

    # Cluster assignments (evenly split)
    labels = np.repeat(np.arange(N_CLUSTERS), n_cells_per_blob)

    # ── 2. Generate cell positions (PCA embedding) ───────────────────
    X_pca = np.zeros((n_cells, N_PCS), dtype=np.float32)
    for i in range(N_CLUSTERS):
        mask = labels == i
        n_i = mask.sum()
        X_pca[mask] = centroids[i] + rng.normal(0, cluster_std, (n_i, N_PCS))

    # ── 3. Assign batches with rotation effects ──────────────────────
    # Each batch gets a random 10×10 orthogonal transformation applied
    # to the first 10 principal components, simulating technical batch
    # variation while preserving within-batch geometry.
    batch_labels = rng.randint(0, N_BATCHES, size=n_cells)
    for b in range(N_BATCHES):
        mask = batch_labels == b
        Q = _random_orthogonal(10, rng)
        # `mask` is a boolean array — combined indexing with a slice works
        # as a single __setitem__ in NumPy, so this assignment is valid.
        rotated = X_pca[mask, :10].copy() @ Q.T
        X_pca[mask, :10] = rotated

    # ── 4. Generate expression data (Poisson) ────────────────────────
    # Per-cluster mean expression profiles:  (N_CLUSTERS, N_GENES)
    cluster_means = np.abs(rng.randn(N_CLUSTERS, N_GENES)) * 10.0 + 1.0

    X = np.zeros((n_cells, N_GENES), dtype=np.float32)
    for i in range(N_CLUSTERS):
        mask = labels == i
        n_i = mask.sum()
        # Cell-specific scaling (captures library-size variation)
        cell_scale = np.abs(rng.randn(n_i, 1)) * 0.2 + 1.0
        lam = cluster_means[i] * cell_scale
        X[mask] = rng.poisson(lam).astype(np.float32)

    # ── 5. Build AnnData ─────────────────────────────────────────────
    adata = ad.AnnData(X)
    adata.obs_names = [f"{label_prefix}cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(N_GENES)]

    # PCA embedding
    adata.obsm["X_pca"] = X_pca

    # Observations
    adata.obs["batch"] = pd.Categorical([f"batch_{b}" for b in batch_labels])
    adata.obs["true_cell_type"] = pd.Categorical([f"cell_type_{i}" for i in labels])

    # Highly variable genes (first N_HVG genes)
    hvg = np.zeros(N_GENES, dtype=bool)
    hvg[:N_HVG] = True
    adata.var["highly_variable"] = hvg

    # Marker dictionary — N_MARKERS_PER_TYPE markers per cell type
    # selected randomly from the top N_HVG_POOL gene indices.
    marker_dict: dict[str, list[str]] = {}
    gene_names = adata.var_names.tolist()
    pool_indices = list(range(min(N_HVG_POOL, N_GENES)))
    for i in range(N_CLUSTERS):
        rng.shuffle(pool_indices)
        markers = [gene_names[j] for j in pool_indices[:N_MARKERS_PER_TYPE]]
        marker_dict[f"cell_type_{i}"] = markers

    adata.uns["marker_dict"] = marker_dict

    return adata


def _verify(adata: ad.AnnData, expected_cells: int, label: str) -> None:
    """Sanity-check an AnnData fixture."""
    assert adata.shape[0] == expected_cells, (
        f"{label}: expected {expected_cells} cells, got {adata.shape[0]}"
    )
    assert adata.shape[1] == N_GENES, f"{label}: expected {N_GENES} genes, got {adata.shape[1]}"

    assert "X_pca" in adata.obsm, f"{label}: missing obsm['X_pca']"
    assert adata.obsm["X_pca"].shape == (expected_cells, N_PCS), (
        f"{label}: obsm['X_pca'] shape mismatch"
    )

    assert "batch" in adata.obs, f"{label}: missing obs['batch']"
    assert adata.obs["batch"].dtype.name == "category", f"{label}: obs['batch'] not categorical"
    assert adata.obs["batch"].nunique() == N_BATCHES, f"{label}: expected {N_BATCHES} batches"

    assert "true_cell_type" in adata.obs, f"{label}: missing obs['true_cell_type']"
    assert adata.obs["true_cell_type"].dtype.name == "category", (
        f"{label}: obs['true_cell_type'] not categorical"
    )
    assert adata.obs["true_cell_type"].nunique() == N_CLUSTERS, (
        f"{label}: expected {N_CLUSTERS} cell types"
    )

    assert "highly_variable" in adata.var, f"{label}: missing var['highly_variable']"
    assert adata.var["highly_variable"].sum() == N_HVG, (
        f"{label}: expected {N_HVG} HVGs, got {adata.var['highly_variable'].sum()}"
    )

    assert "marker_dict" in adata.uns, f"{label}: missing uns['marker_dict']"
    assert len(adata.uns["marker_dict"]) == N_CLUSTERS, (
        f"{label}: marker_dict has {len(adata.uns['marker_dict'])} keys, expected {N_CLUSTERS}"
    )
    for ct, markers in adata.uns["marker_dict"].items():
        assert len(markers) == N_MARKERS_PER_TYPE, (
            f"{label}: cell_type {ct} has {len(markers)} markers, expected {N_MARKERS_PER_TYPE}"
        )

    print(f"  ✓ {label}: {adata.shape[0]} cells × {adata.shape[1]} genes — all checks passed")


def main() -> None:
    os.makedirs(_FIXTURES_DIR, exist_ok=True)

    # ── synthetic_rna.h5ad (well-separated, 5k cells) ────────────────
    print("Generating synthetic_rna.h5ad  (5k cells, well-separated, 5σ) ...")
    adata = generate_dataset(
        n_cells_per_blob=500,
        cluster_std=1.0,
        min_centroid_distance=5.0,
        label_prefix="well_",
    )
    _verify(adata, expected_cells=5000, label="synthetic_rna")
    adata.write(os.path.join(_FIXTURES_DIR, "synthetic_rna.h5ad"), compression="gzip")
    print()

    # ── synthetic_overlapping.h5ad (moderate overlap, 5k cells) ──────
    print("Generating synthetic_overlapping.h5ad  (5k cells, 2.5σ) ...")
    adata = generate_dataset(
        n_cells_per_blob=500,
        cluster_std=1.0,
        min_centroid_distance=2.5,
        label_prefix="overlap_",
    )
    _verify(adata, expected_cells=5000, label="synthetic_overlapping")
    adata.write(
        os.path.join(_FIXTURES_DIR, "synthetic_overlapping.h5ad"),
        compression="gzip",
    )
    print()

    # ── synthetic_severely_overlapping.h5ad (severe overlap) ────────
    print("Generating synthetic_severely_overlapping.h5ad  (5k cells, 1.25σ) ...")
    adata = generate_dataset(
        n_cells_per_blob=500,
        cluster_std=1.0,
        min_centroid_distance=1.25,
        label_prefix="severe_",
    )
    _verify(adata, expected_cells=5000, label="synthetic_severely_overlapping")
    adata.write(
        os.path.join(_FIXTURES_DIR, "synthetic_severely_overlapping.h5ad"),
        compression="gzip",
    )
    print()

    # ── synthetic_100k.h5ad (large, well-separated) ──────────────────
    print("Generating synthetic_100k.h5ad  (100k cells, well-separated) ...")
    adata = generate_dataset(
        n_cells_per_blob=10000,
        cluster_std=1.0,
        min_centroid_distance=5.0,
        label_prefix="100k_",
    )
    _verify(adata, expected_cells=100_000, label="synthetic_100k")
    adata.write(os.path.join(_FIXTURES_DIR, "synthetic_100k.h5ad"), compression="gzip")
    print()

    # ── Summary ───────────────────────────────────────────────────────
    print("All synthetic datasets generated successfully.\n")
    print(f"{'File':42s} {'Size':>8s}  {'Cells':>7s}  {'Genes':>6s}")
    print("-" * 66)
    for fname in (
        "synthetic_rna.h5ad",
        "synthetic_overlapping.h5ad",
        "synthetic_severely_overlapping.h5ad",
        "synthetic_100k.h5ad",
    ):
        fpath = os.path.join(_FIXTURES_DIR, fname)
        sz = os.path.getsize(fpath)
        print(f"{fname:42s} {sz / 1024**2:7.1f} MB")


if __name__ == "__main__":
    main()
