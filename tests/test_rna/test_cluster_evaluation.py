"""Numerical tests for rna/utils/cluster_evaluation.py."""

import logging

import pytest

from rna.utils.cluster_evaluation import (
    select_best_params,
    select_best_umap_params,
)


class TestClusterEvaluationImport:
    """Verify that cluster-evaluation symbols are importable."""

    def test_import_cluster_evaluation(self) -> None:
        assert select_best_params is not None
        assert select_best_umap_params is not None

    def test_select_best_params_callable(self) -> None:
        """select_best_params should be a callable function."""
        assert callable(select_best_params)


class TestSelectBestParams:
    """Numerical assertions for select_best_params."""

    def test_pareto_elbow_with_known_points(self) -> None:
        """3 Pareto points in (k, ss) space; elbow picks the middle point.

        Points: (3, 0.4), (6, 0.55), (9, 0.6) — all Pareto-undominated.
        Normalized: k_norm = [0, 0.5, 1.0], s_norm = [0, 0.75, 1.0].
        Dist to ideal (k=0, s=1): [1.0, 0.559, 1.0] -> argmin=1.
        Expected best: n_neighbors=20, resolution=0.8.
        """
        results = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 6, "silhouette_score": 0.55,
            },
            {
                "n_neighbors": 30, "resolution": 1.0,
                "n_clusters": 9, "silhouette_score": 0.6,
            },
        ]
        best_n, best_r, method, reason = select_best_params(
            results, method="pareto_elbow",
        )

        assert best_n == 20
        assert best_r == pytest.approx(0.8)
        assert method == "pareto_elbow"
        assert "dist_to_ideal" in reason

    def test_pareto_single_point(self) -> None:
        """Single result -> single_pareto_point code path."""
        results = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
            },
        ]
        best_n, best_r, method, reason = select_best_params(
            results, method="pareto_elbow",
        )

        assert best_n == 10
        assert best_r == pytest.approx(0.5)
        assert "single_pareto_point" in reason

    def test_silhouette_method_selects_max(self) -> None:
        """method='silhouette' -> entry with highest silhouette_score."""
        results = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 5, "silhouette_score": 0.7,
            },
            {
                "n_neighbors": 30, "resolution": 1.0,
                "n_clusters": 8, "silhouette_score": 0.6,
            },
        ]
        best_n, best_r, method, _reason = select_best_params(
            results, method="silhouette",
        )

        assert best_n == 20
        assert best_r == pytest.approx(0.8)
        assert method == "silhouette"

    def test_no_valid_scores_raises_value_error(self) -> None:
        """Empty or all-NaN results_summary -> ValueError."""
        with pytest.raises(ValueError, match="No valid silhouette scores"):
            select_best_params([], method="silhouette")

        with pytest.raises(ValueError, match="No valid silhouette scores"):
            select_best_params(
                [
                    {
                        "n_neighbors": 10, "resolution": 0.5,
                        "n_clusters": 3, "silhouette_score": float("nan"),
                    },
                ],
                method="silhouette",
            )

    def test_manual_method_with_resolution(self) -> None:
        """method=None with best_resolution selects correct row."""
        results = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 6, "silhouette_score": 0.7,
            },
            {
                "n_neighbors": 30, "resolution": 1.2,
                "n_clusters": 9, "silhouette_score": 0.6,
            },
        ]
        best_n, best_r, method, _reason = select_best_params(
            results, method=None, best_resolution=0.8,
        )

        assert best_n == 20
        assert best_r == pytest.approx(0.8)
        assert method == "manual"


class MockCFG:
    """Minimal CFG mock for select_best_umap_params manual mode."""
    umap_min_dist = 0.5
    umap_spread = 1.5


class TestSelectBestUmapParams:
    """Numerical assertions for select_best_umap_params."""

    def test_manual_mode_returns_cfg_values(self) -> None:
        """method=None -> returns CFG.umap_min_dist and CFG.umap_spread directly."""
        logger = logging.getLogger(__name__)

        best_md, best_sp, method_label, sweep_results = (
            select_best_umap_params(
                adata=None,
                best_n=10,
                min_dist_grid=None,
                spread_grid=None,
                method=None,
                CFG=MockCFG(),
                use_rep="X_pca",
                log=logger,
            )
        )

        assert best_md == pytest.approx(0.5)
        assert best_sp == pytest.approx(1.5)
        assert method_label == "manual"
        assert sweep_results == []
