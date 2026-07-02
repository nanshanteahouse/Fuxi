#!/usr/bin/env python3
"""Test harness for ``core/clustering.py`` — shared clustering interface.

Uses mock AnnData objects and synthetic callables to verify that
grid_search_clustering, umap_sweep, and select_best_params behave correctly
without depending on real Scanpy / SnapATAC2 clustering runs.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

from core.clustering import (
    grid_search_clustering,
    select_best_params,
    umap_sweep,
)


# ---------------------------------------------------------------------------
#  Test fixtures — mock data & callables
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_adata() -> AnnData:
    """Small synthetic AnnData with 100 cells × 20 genes and a PCA embedding."""
    rng = np.random.RandomState(42)
    X = rng.normal(size=(100, 20))
    adata = AnnData(X)
    adata.obsm["X_pca"] = rng.normal(size=(100, 10))
    return adata


# ---------------------------------------------------------------------------
#  Mock callables (injectable via lambda or class)
# ---------------------------------------------------------------------------


def _mock_neighbor_fn(adata: AnnData, **kwargs: object) -> None:
    """No-op: neighbour computation stub."""
    adata.uns["_mock_neighbor_called"] = True


def _mock_umap_fn(adata: AnnData, **kwargs: object) -> np.ndarray:
    """Return random 2-D embedding and attach it to adata."""
    rng = np.random.RandomState(hash(str(kwargs)) % (2**31))
    coords = rng.normal(size=(adata.n_obs, 2))
    adata.obsm["X_umap"] = coords
    return coords


_counter: int = 0


def _mock_clusterer(adata: AnnData, resolution: float = 1.0, **kwargs: object) -> str:
    """Create 3 clusters deterministically based on resolution."""
    global _counter
    _counter += 1
    key = f"leiden_{resolution}_{_counter}"
    labels = np.array([i % 3 for i in range(adata.n_obs)], dtype=str)
    adata.obs[key] = labels
    return key


def _mock_eval_fn(adata: AnnData, cluster_key: str, **kwargs: object) -> float:
    """Return a fixed silhouette-like score based on cluster count."""
    labels = adata.obs[cluster_key]
    n = labels.nunique()
    return 1.0 - (n - 1) * 0.05  # fewer clusters = higher score


# ---------------------------------------------------------------------------
#  grid_search_clustering tests
# ---------------------------------------------------------------------------


class TestGridSearchClustering:
    """Unit tests for grid_search_clustering."""

    def test_empty_param_grid_raises(self, mock_adata: AnnData) -> None:
        """Empty param_grid must raise ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            grid_search_clustering(mock_adata, {}, _mock_clusterer)

    def test_flat_grid_no_grouping(self, mock_adata: AnnData) -> None:
        """Flat grid without group_key evaluates every combination."""
        param_grid: dict = {"resolution": [0.5, 1.0, 1.5]}
        results = grid_search_clustering(
            mock_adata,
            param_grid,
            clusterer=_mock_clusterer,
            evaluation_fn=_mock_eval_fn,
            random_seed=42,
        )
        assert len(results) == 3
        for r in results:
            assert "resolution" in r
            assert r["n_clusters"] == 3
            assert isinstance(r["score"], float)
            assert r["score"] > 0.8

    def test_grouped_grid(self, mock_adata: AnnData) -> None:
        """Grouped grid: n_neighbors groups, resolution is inner loop."""
        param_grid: dict = {
            "n_neighbors": [15, 20],
            "resolution": [0.5, 1.0],
        }
        neighbor_calls: list = []

        def tracking_neighbor(adata: AnnData, **kw: object) -> None:
            neighbor_calls.append(kw.get("n_neighbors"))
            _mock_neighbor_fn(adata, **kw)

        results = grid_search_clustering(
            mock_adata,
            param_grid,
            clusterer=_mock_clusterer,
            neighbor_fn=tracking_neighbor,
            umap_fn=_mock_umap_fn,
            evaluation_fn=_mock_eval_fn,
            group_key="n_neighbors",
            random_seed=42,
        )
        # 2 n_neighbors × 2 resolutions = 4 results
        assert len(results) == 4
        # Neighbour called exactly once per n_neighbors value
        assert neighbor_calls == [15, 20]

        # Every result has both param keys
        for r in results:
            assert "n_neighbors" in r
            assert "resolution" in r
            assert "score" in r

    def test_invalid_group_key_raises(self, mock_adata: AnnData) -> None:
        """group_key not in param_grid must raise ValueError."""
        param_grid: dict = {"resolution": [0.5, 1.0]}
        with pytest.raises(ValueError, match="group_key"):
            grid_search_clustering(
                mock_adata,
                param_grid,
                clusterer=_mock_clusterer,
                group_key="n_neighbors",
            )

    def test_clusterer_failure_is_recorded(self, mock_adata: AnnData) -> None:
        """When clusterer raises, the combo is recorded with error field."""

        def failing_clusterer(adata: AnnData, **kw: object) -> str:
            raise RuntimeError("boom")

        param_grid: dict = {"resolution": [0.5]}
        results = grid_search_clustering(
            mock_adata, param_grid, clusterer=failing_clusterer, random_seed=42
        )
        assert len(results) == 1
        assert results[0]["error"] == "clusterer failed"

    def test_neighbor_failure_skips_group(self, mock_adata: AnnData) -> None:
        """When neighbor_fn raises, the entire group is skipped."""

        def failing_neighbor(adata: AnnData, **kw: object) -> None:
            raise RuntimeError("knn failed")

        param_grid: dict = {"n_neighbors": [15], "resolution": [0.5, 1.0]}
        results = grid_search_clustering(
            mock_adata,
            param_grid,
            clusterer=_mock_clusterer,
            neighbor_fn=failing_neighbor,
            group_key="n_neighbors",
            random_seed=42,
        )
        assert len(results) == 0  # entire group skipped

    def test_umap_failure_non_fatal(self, mock_adata: AnnData) -> None:
        """UMAP failure should not prevent clustering."""

        def failing_umap(adata: AnnData, **kw: object) -> None:
            raise RuntimeError("umap failed")

        param_grid: dict = {"n_neighbors": [15], "resolution": [0.5]}
        results = grid_search_clustering(
            mock_adata,
            param_grid,
            clusterer=_mock_clusterer,
            umap_fn=failing_umap,
            group_key="n_neighbors",
            random_seed=42,
        )
        # Clustering still proceeds; UMAP failure is non-fatal
        assert len(results) == 1
        assert results[0]["n_clusters"] == 3

    def test_no_evaluation_fn(self, mock_adata: AnnData) -> None:
        """When evaluation_fn is None, score is omitted."""
        param_grid: dict = {"resolution": [0.5]}
        results = grid_search_clustering(
            mock_adata, param_grid, clusterer=_mock_clusterer, random_seed=42
        )
        assert len(results) == 1
        assert "score" not in results[0]
        assert results[0]["n_clusters"] == 3

    def test_results_shape_matches_grid(self, mock_adata: AnnData) -> None:
        """Number of results equals product of param value counts."""
        param_grid: dict = {
            "n_neighbors": [15, 20, 30],
            "resolution": [0.5, 1.0, 1.5],
        }
        results = grid_search_clustering(
            mock_adata,
            param_grid,
            clusterer=_mock_clusterer,
            neighbor_fn=_mock_neighbor_fn,
            umap_fn=_mock_umap_fn,
            group_key="n_neighbors",
            random_seed=42,
        )
        expected_count = 3 * 3  # 9 combinations
        assert len(results) == expected_count


# ---------------------------------------------------------------------------
#  select_best_params tests
# ---------------------------------------------------------------------------


class TestSelectBestParams:
    """Regression tests for the Pareto-elbow selection re-export."""

    def test_import_is_callable(self) -> None:
        """select_best_params is callable (re-export from rna.utils)."""
        assert callable(select_best_params)

    def test_pareto_elbow_selection(self) -> None:
        """Smoke test: known grid returned by grid_search_clustering."""
        summary = [
            {"n_neighbors": 15, "resolution": 0.5, "n_clusters": 5, "silhouette_score": 0.60},
            {"n_neighbors": 15, "resolution": 1.0, "n_clusters": 8, "silhouette_score": 0.55},
            {"n_neighbors": 20, "resolution": 0.5, "n_clusters": 6, "silhouette_score": 0.62},
            {"n_neighbors": 20, "resolution": 1.0, "n_clusters": 9, "silhouette_score": 0.50},
        ]
        best_n, best_r, method_label, reason = select_best_params(summary, method="pareto_elbow")
        assert isinstance(best_n, int)
        assert isinstance(best_r, float)
        assert method_label == "pareto_elbow"
        assert isinstance(reason, str)

    def test_silhouette_selection(self) -> None:
        """Silhouette method picks the highest score."""
        summary = [
            {"n_neighbors": 15, "resolution": 0.5, "n_clusters": 3, "silhouette_score": 0.60},
            {"n_neighbors": 20, "resolution": 0.5, "n_clusters": 4, "silhouette_score": 0.75},
            {"n_neighbors": 30, "resolution": 0.5, "n_clusters": 5, "silhouette_score": 0.55},
        ]
        best_n, best_r, method_label, _ = select_best_params(summary, method="silhouette")
        assert best_n == 20
        assert best_r == 0.5
        assert method_label == "silhouette"

    def test_no_valid_silhouette_raises(self) -> None:
        """All-None silhouette scores must raise ValueError."""
        summary = [
            {"n_neighbors": 15, "resolution": 0.5, "n_clusters": 3, "silhouette_score": None},
        ]
        with pytest.raises(ValueError, match="No valid silhouette"):
            select_best_params(summary)

    def test_unknown_method_raises(self) -> None:
        """Unknown method name must raise ValueError."""
        summary = [
            {"n_neighbors": 15, "resolution": 0.5, "n_clusters": 3, "silhouette_score": 0.60},
        ]
        with pytest.raises(ValueError, match="Unknown cluster_selection_method"):
            select_best_params(summary, method="magic")


# ---------------------------------------------------------------------------
#  umap_sweep tests
# ---------------------------------------------------------------------------


class TestUmapSweep:
    """Unit tests for umap_sweep."""

    def test_empty_param_sweep_raises(self, mock_adata: AnnData) -> None:
        """Empty param_sweep must raise ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            umap_sweep(mock_adata, {})

    def test_basic_sweep(self, mock_adata: AnnData) -> None:
        """Basic 2×1 sweep returns correct number of results."""
        param_sweep: dict = {"min_dist": [0.1, 0.3], "spread": [1.0]}

        def area_eval(adata: AnnData, coords: np.ndarray, **kw: object) -> float:
            return float(np.prod(coords.max(axis=0) - coords.min(axis=0)))

        results = umap_sweep(
            mock_adata,
            param_sweep,
            umap_fn=_mock_umap_fn,
            evaluation_fn=area_eval,
            random_seed=42,
        )
        assert len(results) == 2  # 2 min_dist × 1 spread
        for r in results:
            assert "min_dist" in r
            assert "spread" in r
            assert isinstance(r["score"], float)
            assert r["score"] > 0

    def test_umap_failure_recorded(self, mock_adata: AnnData) -> None:
        """When umap_fn raises, the entry records the error."""

        def failing_umap(adata: AnnData, **kw: object) -> np.ndarray:
            raise RuntimeError("no UMAP for you")

        param_sweep: dict = {"min_dist": [0.1]}
        results = umap_sweep(
            mock_adata, param_sweep, umap_fn=failing_umap, random_seed=42
        )
        assert len(results) == 1
        assert results[0]["score"] is None
        assert results[0]["error"] == "UMAP failed"

    def test_no_umap_fn_fallback(self, mock_adata: AnnData) -> None:
        """When umap_fn is None, existing X_umap in obsm is used."""
        mock_adata.obsm["X_umap"] = np.random.RandomState(42).normal(size=(100, 2))
        param_sweep: dict = {"min_dist": [0.1]}

        results = umap_sweep(mock_adata, param_sweep, random_seed=42)
        assert len(results) == 1
        assert "min_dist" in results[0]
        # No score because evaluation_fn is None
        assert "score" not in results[0]

    def test_noop_umap_fn_returns_existing(self, mock_adata: AnnData) -> None:
        """_noop_umap_fn returns whatever is in adata.obsm['X_umap']."""
        expected = np.random.RandomState(42).normal(size=(100, 2))
        mock_adata.obsm["X_umap"] = expected
        from core.clustering import _noop_umap_fn

        result = _noop_umap_fn(mock_adata)
        np.testing.assert_array_equal(result, expected)

    def test_sweep_shape_matches_product(self, mock_adata: AnnData) -> None:
        """Result count equals product of sweep dimensions."""
        param_sweep: dict = {"min_dist": [0.1, 0.3], "spread": [0.5, 1.0, 2.0]}
        results = umap_sweep(
            mock_adata, param_sweep, umap_fn=_mock_umap_fn, random_seed=42
        )
        assert len(results) == 2 * 3  # 6 combinations


# ---------------------------------------------------------------------------
#  Integration-style tests
# ---------------------------------------------------------------------------


class TestEndToEnd:
    """End-to-end: mock grid search → select_best_params → verify."""

    def test_full_workflow_scanpy_like(self, mock_adata: AnnData) -> None:
        """Simulate a Scanpy-like workflow: neighbor → UMAP → Leiden per resolution."""
        param_grid: dict = {
            "n_neighbors": [15, 20],
            "resolution": [0.5, 1.0, 2.0],
        }

        # Build summary in the shape select_best_params expects
        raw = grid_search_clustering(
            mock_adata,
            param_grid,
            clusterer=_mock_clusterer,
            neighbor_fn=_mock_neighbor_fn,
            umap_fn=_mock_umap_fn,
            evaluation_fn=_mock_eval_fn,
            group_key="n_neighbors",
            random_seed=42,
        )

        # Convert to the summary format expected by select_best_params
        summary = []
        for r in raw:
            summary.append({
                "n_neighbors": r["n_neighbors"],
                "resolution": r["resolution"],
                "n_clusters": r["n_clusters"],
                "silhouette_score": r["score"],
            })

        best_n, best_r, method_label, reason = select_best_params(summary, method="pareto_elbow")
        assert best_n in (15, 20)
        assert best_r in (0.5, 1.0, 2.0)
        assert method_label == "pareto_elbow"
        assert "dist_to_ideal" in reason or "single_pareto_point" in reason
