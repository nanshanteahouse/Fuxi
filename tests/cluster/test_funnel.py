"""Tests for funnel grid search — progressive subsample for large datasets.

Tests cover:
1. Large dataset triggers funnel (subsample occurs when n_obs > size)
2. Small dataset skips funnel (full data returned when n_obs <= size)
3. target_n_clusters overrides funnel (caller bypass)
4. Rare cluster preservation via subsample_stratified
5. Per-batch proportional subsampling
6. funnel_lineage written to adata.uns with required keys
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from core.cluster.funnel import run_funnel_grid_search, subsample_stratified


def _should_use_funnel(cfg: object) -> bool:
    """Return True if funnel grid search should be used by the caller.

    When target_n_clusters is set the step script bypasses funnel and uses
    target-based clustering instead.
    """
    target = getattr(cfg.clustering, "target_n_clusters", None)  # type: ignore[union-attr]
    return target is None


# ---------------------------------------------------------------------------
# Test 1 — large dataset triggers funnel
# ---------------------------------------------------------------------------


class TestFunnelTrigger:
    """Funnel triggers only when n_obs > subsample_size."""

    def test_large_dataset_triggers_funnel(self) -> None:
        """n_obs (300) > size (100) → subsample_stratified returns a subset."""
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 300
        adata = sc.AnnData(rng.randn(n_cells, 20))
        adata.obsm["X_pca"] = rng.randn(n_cells, 10)
        adata.obs["batch"] = ["batch_0"] * n_cells

        sub, idx = subsample_stratified(adata, size=100, use_rep="X_pca")
        assert sub.n_obs < n_cells, "Large dataset should be subsampled"
        assert len(idx) < n_cells
        assert len(idx) > 0
        assert sub.n_obs == len(idx), "Subset must match idx length"
        assert sub is not adata, "Subsample must be a copy"

    def test_small_dataset_skips_funnel(self) -> None:
        """n_obs (30) <= size (50) → subsample_stratified returns full data."""
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 30
        adata = sc.AnnData(rng.randn(n_cells, 20))
        adata.obsm["X_pca"] = rng.randn(n_cells, 10)

        sub, idx = subsample_stratified(adata, size=50, use_rep="X_pca")
        assert sub is adata, "Small dataset must return the same object"
        assert len(idx) == n_cells
        assert np.array_equal(idx, np.arange(n_cells))


# ---------------------------------------------------------------------------
# Test 3 — target_n_clusters overrides funnel
# ---------------------------------------------------------------------------


class TestTargetOverride:
    """When target_n_clusters is set, caller must bypass funnel."""

    @pytest.mark.parametrize(
        ("cfg", "expected"),
        [
            (
                SimpleNamespace(clustering=SimpleNamespace(target_n_clusters=10)),
                False,
            ),
            (
                SimpleNamespace(clustering=SimpleNamespace(target_n_clusters=0)),
                False,
            ),
            (
                SimpleNamespace(clustering=SimpleNamespace(funnel_subsample_size=50000)),
                True,
            ),
        ],
    )
    def test_target_overrides_funnel(self, cfg: object, expected: bool) -> None:
        """_should_use_funnel returns False when target_n_clusters is set."""
        assert _should_use_funnel(cfg) is expected


# ---------------------------------------------------------------------------
# Test 4 — rare cluster preservation
# ---------------------------------------------------------------------------


class TestRareClusterPreservation:
    """subsample_stratified preserves rare clusters (≥5 cells)."""

    def test_subsample_preserves_rare_clusters(self) -> None:
        """Rare cluster (15 cells of 300) keeps ≥5 cells after subsampling."""
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 300
        adata = sc.AnnData(rng.randn(n_cells, 20))
        adata.obsm["X_pca"] = rng.randn(n_cells, 10)
        # Single batch to keep things deterministic
        adata.obs["batch"] = ["batch_0"] * n_cells
        # 285 common cells + 15 rare cells (rare = small enough that KMeans
        # may undersample them, triggering the preservation logic)
        adata.obs["true_cell_type"] = ["common"] * 285 + ["rare"] * 15

        sub, idx = subsample_stratified(adata, size=150, use_rep="X_pca")
        rare_count = int((sub.obs["true_cell_type"] == "rare").sum())
        assert rare_count >= 5, f"Rare cluster should have ≥5 cells in subsample, got {rare_count}"
        # The rare cluster must not be lost entirely
        assert "rare" in sub.obs["true_cell_type"].values, "Rare cluster missing"
        # Common cluster proportions should still be reasonable
        total = sub.n_obs
        assert total == len(idx)


# ---------------------------------------------------------------------------
# Test 5 — per-batch proportional subsampling
# ---------------------------------------------------------------------------


class TestBatchProportional:
    """Batch proportions are maintained after subsampling."""

    def test_subsample_per_batch_proportional(self) -> None:
        """Batch proportions in subsample match original within ±5 %."""
        import scanpy as sc

        rng = np.random.RandomState(42)
        n_cells = 200
        adata = sc.AnnData(rng.randn(n_cells, 20))
        adata.obsm["X_pca"] = rng.randn(n_cells, 10)
        # Two batches with 75/25 split
        adata.obs["batch"] = ["batch_0"] * 150 + ["batch_1"] * 50
        # Add cell-type column so rare-cluster preservation runs (it must
        # not interfere with batch proportions for well-sampled clusters)
        adata.obs["true_cell_type"] = ["ct_0"] * n_cells

        sub, idx = subsample_stratified(adata, size=100, use_rep="X_pca")
        total = sub.n_obs
        assert total == len(idx)

        orig_pct_0 = 150 / 200 * 100  # 75 %
        sub_pct_0 = (sub.obs["batch"] == "batch_0").sum() / total * 100
        sub_pct_1 = (sub.obs["batch"] == "batch_1").sum() / total * 100

        assert abs(sub_pct_0 - orig_pct_0) <= 5, (
            f"batch_0 proportion {sub_pct_0:.1f}% deviates >5% from original {orig_pct_0:.1f}%"
        )
        # batch_1 is the complement: 100 % - sub_pct_0
        assert abs(sub_pct_1 - 25.0) <= 5, (
            f"batch_1 proportion {sub_pct_1:.1f}% deviates >5% from original 25%"
        )


# ---------------------------------------------------------------------------
# Test 6 — funnel_lineage written
# ---------------------------------------------------------------------------


class TestFunnelLineage:
    """run_funnel_grid_search writes funnel_lineage to adata.uns."""

    @patch("core.cluster.funnel._validate_on_full")
    def test_funnel_lineage_written(self, mock_validate: object, synthetic_adata: object) -> None:
        """Lineage dict contains all required keys after funnel run."""

        cfg = SimpleNamespace(
            clustering=SimpleNamespace(
                funnel_subsample_size=50000,
                funnel_top_k=3,
                multi_metric_weights=None,
            ),
            pca=SimpleNamespace(n_pcs_use=30),
            execution=SimpleNamespace(random_seed=42),
        )

        # full_grid_fn — returns enriched-style results (silhouette_score
        # is the only non-zero metric; _compute_composite_scores degrades
        # gracefully and produces a valid composite_score).
        def mock_full_grid_fn(adata: object, cfg: object) -> list[dict[str, object]]:
            return [
                {
                    "n_neighbors": 15,
                    "resolution": 0.5,
                    "cluster_key": "leiden_15_0.5",
                    "n_clusters": 5,
                    "silhouette_score": 0.65,
                },
                {
                    "n_neighbors": 20,
                    "resolution": 1.0,
                    "cluster_key": "leiden_20_1.0",
                    "n_clusters": 8,
                    "silhouette_score": 0.72,
                },
            ]

        # Mock the expensive full-data validation step
        mock_validate.return_value = [
            {
                "n_neighbors": 15,
                "resolution": 0.5,
                "cluster_key": "funnel_15_0.5",
                "n_clusters": 5,
                "silhouette_score": 0.65,
                "composite_score": 0.70,
            },
            {
                "n_neighbors": 20,
                "resolution": 1.0,
                "cluster_key": "funnel_20_1.0",
                "n_clusters": 8,
                "silhouette_score": 0.72,
                "composite_score": 0.85,
            },
        ]

        # synthetic_adata has 5k cells; subsample_size=50k ensures
        # n_obs <= size so _validate_on_full receives full data quickly.
        run_funnel_grid_search(
            synthetic_adata,
            cfg,
            mock_full_grid_fn,  # type: ignore[arg-type]
        )

        assert "funnel_lineage" in synthetic_adata.uns
        lineage = synthetic_adata.uns["funnel_lineage"]

        required_keys = (
            "subsample_n",
            "subsample_idx",
            "subsample_top_k_composites",
            "full_top_k_composites",
            "best_key",
            "composite_delta_full_minus_sub",
        )
        for key in required_keys:
            assert key in lineage, f"Missing required key: {key!r}"

        # Type/sanity checks
        assert isinstance(lineage["subsample_n"], int)
        assert isinstance(lineage["subsample_idx"], list)
        assert len(lineage["subsample_idx"]) > 0
        assert isinstance(lineage["subsample_top_k_composites"], list)
        assert isinstance(lineage["full_top_k_composites"], list)
        assert isinstance(lineage["best_key"], str)
        assert isinstance(lineage["composite_delta_full_minus_sub"], float)
        # best_key should reference a funnel_* cluster key from full validation
        assert lineage["best_key"].startswith("funnel_")
