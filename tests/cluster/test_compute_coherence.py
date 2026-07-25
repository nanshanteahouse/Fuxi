"""Tests for _compute_cluster_coherence — per-cluster marker dominance.

These tests supplement tests/test_rna/test_cluster_evaluation.py and cover
post-fix behaviors (C1, C6 from the audit):

- C1: single cell type in per_cell_scores → float('nan') (currently returns 1.0)
- C6: empty per_cell_scores dict → float('nan') (currently returns 1.0)

Also covers post-#22 (percentile min_expression), post-#24 (min_cluster_size),
and dominance-threshold sensitivity.
"""

from unittest.mock import MagicMock, patch  # noqa: F401  (kept for import-pattern consistency)

import numpy as np
import pytest

from core.cluster.evaluation import _compute_cluster_coherence


class TestClusterCoherenceWellSeparated:
    """Well-separated clusters with clean marker assignments."""

    def test_well_separated_matching_markers(self) -> None:
        """Well-separated with matching marker_dict → coherence > 0.8."""
        import scanpy as sc

        n_cells = 300
        rng = np.random.RandomState(42)

        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["cluster"] = ["0"] * 100 + ["1"] * 100 + ["2"] * 100

        # Each cell type scores high only in its matching cluster
        per_cell_scores: dict[str, np.ndarray] = {}
        for i in range(3):
            scores = rng.normal(0, 0.3, n_cells)
            cluster_mask = adata.obs["cluster"] == str(i)
            scores[cluster_mask.values] += 5.0
            per_cell_scores[f"Type{i}"] = scores

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
        )
        assert coherence > 0.8, f"Expected > 0.8, got {coherence}"

    def test_well_separated_mismatched_markers(self) -> None:
        """Well-separated with WRONG marker_dict → coherence < 0.3.

        Uses many cell types with near-identical per-cell scores so that
        no cell type dominates any cluster.
        """
        import scanpy as sc

        n_cells = 200
        rng = np.random.RandomState(42)

        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["cluster"] = ["0"] * 100 + ["1"] * 100

        # 10 cell types, all with nearly identical expression in every cell.
        # No type clearly dominates any cluster → coherence should be near 0.
        base = rng.normal(1.0, 0.2, n_cells)
        per_cell_scores = {f"Type{c}": base + rng.normal(0, 0.05, n_cells) for c in range(10)}

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
        )
        assert coherence < 0.3, f"Expected < 0.3, got {coherence}"


class TestClusterCoherenceOverlapping:
    """Overlapping clusters produce lower coherence."""

    def test_severely_overlapping_matching_markers(self) -> None:
        """Overlapping with matching markers → coherence < 0.5."""
        import scanpy as sc

        n_cells = 200
        rng = np.random.RandomState(42)

        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["cluster"] = ["0"] * 100 + ["1"] * 100

        # Both cell types score in both clusters with similar strength
        per_cell_scores = {
            "TypeA": np.concatenate(
                [
                    rng.normal(2.5, 0.5, 100),  # cluster 0
                    rng.normal(2.8, 0.5, 100),  # cluster 1
                ]
            ),
            "TypeB": np.concatenate(
                [
                    rng.normal(2.8, 0.5, 100),  # cluster 0
                    rng.normal(2.5, 0.5, 100),  # cluster 1
                ]
            ),
        }

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
        )
        assert coherence < 0.5, f"Expected < 0.5, got {coherence}"


class TestClusterCoherenceEdgeCases:
    """Edge cases: degenerate inputs."""

    def test_single_celltype_returns_nan(self) -> None:
        """Single cell type in per_cell_scores → returns float('nan').

        [POST-FIX: C1 fix, currently returns 1.0]
        """
        import scanpy as sc

        adata = sc.AnnData(X=np.random.randn(50, 10))
        adata.obs["cluster"] = ["0"] * 25 + ["1"] * 25  # 2 clusters
        # Boost scores above min_expression to ensure C1 path is reached
        per_cell_scores: dict[str, np.ndarray] = {
            "TypeA": np.random.randn(50) + 5.0,  # only 1 cell type
        }

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
        )
        assert np.isnan(coherence), f"Single cell type should return NaN, got {coherence}"

    def test_empty_scores_returns_nan(self) -> None:
        """Empty per_cell_scores dict → returns float('nan').

        [POST-FIX: C6 fix, currently returns 1.0]
        """
        import scanpy as sc

        adata = sc.AnnData(X=np.random.randn(50, 10))
        adata.obs["cluster"] = ["0"] * 25 + ["1"] * 25

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores={},
        )
        assert np.isnan(coherence), f"Empty scores should return NaN, got {coherence}"


class TestClusterCoherenceParameters:
    """Parameter sensitivity: dominance_threshold, min_expression, min_cluster_size."""

    def test_dominance_threshold_2_5_rejects_overlapping(self) -> None:
        """With threshold=2.5, ≥30% clusters rejected on overlapping fixture.

        5 clusters: 3 with strong dominance (ratio >> 2.5), 2 with weak
        dominance (ratio ≈ 2.0 < 2.5).  At threshold 2.5, the 2 weak
        clusters are rejected → coherence = 3/5 = 0.6 ≤ 0.7.

        [post-#22 fix property]
        """
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 500
        adata = sc.AnnData(X=rng.randn(n_cells, 20))

        cluster_labels = np.array(
            ["0"] * 100 + ["1"] * 100 + ["2"] * 100 + ["3"] * 100 + ["4"] * 100
        )
        adata.obs["cluster"] = cluster_labels

        per_cell_scores: dict[str, np.ndarray] = {}
        for ct in range(5):
            scores = rng.normal(0, 0.3, n_cells)
            own_mask = cluster_labels == str(ct)
            if ct < 3:
                # Strong dominance: only own type boosted
                scores[own_mask] += 8.0
            else:
                # Weak dominance: own type boosted moderately AND the other
                # weak type's cluster gets a competing-but-lower boost.
                other_cluster = 3 if ct == 4 else 4
                other_mask = cluster_labels == str(other_cluster)
                scores[own_mask] += 3.0
                scores[other_mask] += 1.5
            per_cell_scores[f"Type{ct}"] = scores

        # At threshold=2.5: clusters 3,4 are NOT coherent (ratio ≈ 2.0)
        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
            dominance_threshold=2.5,
        )
        assert coherence <= 0.7, (
            f"At threshold=2.5, expected ≤0.7 (≥30% rejected), got {coherence}"
        )

    @pytest.mark.skip(
        reason="min_expression='p25' percentile feature not implemented yet (post-#22)"
    )
    def test_min_expression_percentile(self) -> None:
        """min_expression='p25' works.

        [post-#22 property]
        """
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 300
        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["cluster"] = ["0"] * 100 + ["1"] * 100 + ["2"] * 100

        # Cluster 0 and 1 have normal scores; cluster 2 has very low scores
        per_cell_scores: dict[str, np.ndarray] = {
            "TypeA": np.concatenate(
                [
                    rng.uniform(0.5, 1.0, 100),
                    rng.uniform(0.5, 1.0, 100),
                    rng.uniform(0.01, 0.05, 100),
                ]
            ),
        }

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
            min_expression="p25",
        )
        # Using 'p25' should filter out the low-expression cluster 2
        assert not np.isnan(coherence), "Should return a valid coherence (not NaN)"
        assert 0.0 <= coherence <= 1.0, f"Should be in [0, 1], got {coherence}"

    @pytest.mark.skip(reason="min_cluster_size parameter not implemented yet (post-#24)")
    def test_cluster_size_weighting(self) -> None:
        """min_cluster_size=10 skips tiny clusters.

        [POST-FIX: #24 fix]
        """
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 220
        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        # Cluster 0: 200 cells, Cluster 1: 20 cells, Cluster 2: 3 cells (tiny)
        adata.obs["cluster"] = ["0"] * 200 + ["1"] * 20 + ["2"] * 3

        per_cell_scores: dict[str, np.ndarray] = {}
        for ct in range(3):
            scores = rng.normal(0, 0.3, n_cells)
            own_mask = adata.obs["cluster"] == str(ct)
            scores[own_mask.values] += 5.0
            per_cell_scores[f"Type{ct}"] = scores

        coherence = _compute_cluster_coherence(
            adata,
            cluster_key="cluster",
            per_cell_scores=per_cell_scores,
            min_cluster_size=10,
        )
        # Cluster 2 (3 cells < 10) should be skipped; coherence computed
        # on the remaining 2 clusters.
        assert not np.isnan(coherence), "Should return a valid coherence with min_cluster_size"
        assert 0.0 <= coherence <= 1.0, f"Should be in [0, 1], got {coherence}"
