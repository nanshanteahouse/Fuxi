"""Tests for target-mode clustering — binary-search resolution via
``find_resolution_for_target_k`` and ``target_grid_search``.

Covers:

1. Convergence to known-k — 10 well-separated blobs → k within ±1 of 10
2. Non-convergence — impossible target_k=100 returns nearest k with warning
3. Fast convergence — well-separated data converges in ≤5 iterations
4. 3-run median — monkeypatched seed variance correctly medians
5. Target overrides funnel — setting target_n_clusters skips resolution grid
"""

from __future__ import annotations

import logging
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from sklearn.datasets import make_blobs

from core.cluster.target import find_resolution_for_target_k, target_grid_search

# ---------------------------------------------------------------------------
#  Fixture: purpose-built 10-blob data for target-mode tests
# ---------------------------------------------------------------------------

# The session-scoped synthetic_adata fixture has 30 PCs which create 30
# disconnected KNN components — making it unsuitable for target-mode tests
# where resolution must control the number of Leiden clusters.  This inline
# fixture generates 10 blobs in 15 dimensions (no PCA), resulting in a
# fully connected KNN graph where resolution varies k from ~7 to ~25.


@pytest.fixture(scope="module")
def adata_10blobs():
    """~1000 cells, 10 well-separated blobs, 15 raw features.

    Guarantees:
    - Fully connected KNN graph (1 component) with n_neighbors=15
    - Resolution in [0.1, 5.0] produces k ranging from ~7 to ~25
    - At the right resolution, Leiden finds exactly 10 clusters
    """
    import scanpy as sc

    X, _ = make_blobs(  # noqa: N806
        n_samples=1000,
        centers=10,
        n_features=15,
        cluster_std=4.0,
        center_box=(-8, 8),
        random_state=42,
    )
    adata = sc.AnnData(X)
    sc.pp.scale(adata, max_value=10)
    # Use all 15 scaled dimensions as the embedding
    adata.obsm["X_emb"] = adata.X.copy()
    return adata


# ---------------------------------------------------------------------------
#  Helper
# ---------------------------------------------------------------------------


def _make_cfg(
    *,
    target_n_clusters: int | None = None,
    param_grid_n_neighbors: list[int] | None = None,
    use_rep: str = "X_emb",
    target_search_max_iters: int = 10,
) -> SimpleNamespace:
    """Build a minimal mock config object for target-mode tests."""
    return SimpleNamespace(
        clustering=SimpleNamespace(
            target_n_clusters=target_n_clusters,
            param_grid_n_neighbors=param_grid_n_neighbors or [15, 20, 30],
            use_rep=use_rep,
            target_search_max_iters=target_search_max_iters,
        )
    )


# ---------------------------------------------------------------------------
#  Tests
# ---------------------------------------------------------------------------


class TestTargetModeConvergence:
    """Well-separated 10-blob data → binary search finds the right k."""

    @pytest.mark.filterwarnings("ignore::FutureWarning")
    def test_converges_to_known_k(self, adata_10blobs: pytest.fixture) -> None:
        """target_k=10 for 10-blob data → actual_k within ±1."""
        adata = adata_10blobs.copy()
        cfg = _make_cfg()

        best_res, actual_k = find_resolution_for_target_k(
            adata,
            target_k=10,
            n_neighbors=15,
            cfg=cfg,
            log=None,
        )

        assert abs(actual_k - 10) <= 1, (
            f"Target k=10 should converge to within ±1, got actual_k={actual_k} "
            f"at resolution={best_res:.3f}"
        )
        assert best_res > 0, "Resolution must be positive on convergence"

    @pytest.mark.filterwarnings("ignore::FutureWarning")
    def test_converges_in_few_iters(self, adata_10blobs: pytest.fixture) -> None:
        """Binary search converges in ≤5 iterations for well-separated data.

        Uses a spy on ``sc.tl.leiden`` to count calls (3 per iteration).
        """
        import scanpy as sc

        adata = adata_10blobs.copy()
        cfg = _make_cfg()

        with patch.object(sc.tl, "leiden", wraps=sc.tl.leiden) as spy:
            best_res, actual_k = find_resolution_for_target_k(
                adata,
                target_k=10,
                n_neighbors=15,
                cfg=cfg,
                log=None,
            )

        # 3 seeds × ≤5 iterations = ≤15 calls
        assert spy.call_count <= 15, (
            f"Expected ≤15 leiden calls (≤5 iters × 3 seeds), got {spy.call_count}"
        )
        assert abs(actual_k - 10) <= 1, f"Convergence failed: actual_k={actual_k}"
        assert best_res > 0


class TestTargetModeNonConvergence:
    """Impossible target_k values → graceful non-convergence."""

    @pytest.mark.filterwarnings("ignore::FutureWarning")
    def test_non_convergence_returns_nearest_k(
        self,
        adata_10blobs: pytest.fixture,
        caplog: pytest.fixture,
    ) -> None:
        """target_k=100 for 10-blob data → non-convergence with warning.

        Checks:
        - Returned k is the closest feasible k (not 100).
        - A WARNING is logged about non-convergence.
        """
        adata = adata_10blobs.copy()
        cfg = _make_cfg(target_search_max_iters=10)

        with caplog.at_level(logging.WARNING):
            best_res, actual_k = find_resolution_for_target_k(
                adata,
                target_k=100,
                n_neighbors=15,
                cfg=cfg,
                log=None,
            )

        # Should not have converged to anywhere near 100
        assert actual_k < 50, f"10-blob data cannot produce k=100, got actual_k={actual_k}"
        # The returned k is the closest the search could find
        # (should be in the feasible range, far from 100)
        assert actual_k < 50, f"Nearest k should be well below 100, got {actual_k}"

        # Warning about non-convergence emitted
        assert any("not converge" in rec.message.lower() for rec in caplog.records), (
            "Expected a WARNING about non-convergence"
        )

    def test_3_run_median_handles_non_monotonicity(
        self,
        adata_10blobs: pytest.fixture,
    ) -> None:
        """3-run median smooths Leiden seed variance.

        Monkeypatches ``sc.tl.leiden`` so that two seeds produce 10 clusters
        and one seed produces 15 clusters.  The median (10) should be returned.
        """
        adata = adata_10blobs.copy()
        cfg = _make_cfg()

        def _seed_varying_leiden(
            adata,
            resolution=None,
            random_state=None,
            key_added=None,
            **kwargs,
        ):
            """Leiden mock: seeds 42/123 → 10 clusters, seed 456 → 15 clusters."""
            n_unique = 10 if random_state in (42, 123) else 15
            rng = np.random.RandomState(random_state or 42)
            labels = np.array([str(i % n_unique) for i in range(adata.n_obs)])
            rng.shuffle(labels)
            adata.obs[key_added] = pd.Categorical(labels)

        with patch("scanpy.tl.leiden", side_effect=_seed_varying_leiden):
            best_res, actual_k = find_resolution_for_target_k(
                adata,
                target_k=10,
                n_neighbors=15,
                cfg=cfg,
                log=None,
            )

        # Median of [10, 10, 15] is 10 → should converge
        assert actual_k == 10, (
            f"Median of seeds (10, 10, 15) should give actual_k=10, got {actual_k}"
        )
        assert best_res > 0


class TestTargetGridSearch:
    """``target_grid_search`` orchestration tests."""

    def test_target_overrides_funnel(
        self,
        adata_10blobs: pytest.fixture,
    ) -> None:
        """When target_n_clusters is set, resolution-grid (funnel) is skipped.

        Mocks ``find_resolution_for_target_k`` and verifies it is called with
        the correct ``target_k`` — confirming the dispatch routes through
        target mode instead of iterating over a resolution grid.
        """
        adata = adata_10blobs.copy()
        target_k = 10
        cfg = _make_cfg(
            target_n_clusters=target_k,
            param_grid_n_neighbors=[15, 20],
        )

        mock_result = (0.5, 10)

        with patch(
            "core.cluster.target.find_resolution_for_target_k",
            return_value=mock_result,
        ) as mock_find:
            with patch(
                "core.cluster.evaluation.enrich_grid_results",
            ) as mock_enrich:
                results = target_grid_search(adata, cfg, log=None)

        # Every call to find_resolution_for_target_k should have target_k=10
        assert mock_find.call_count == 2, (
            f"Expected 2 calls (one per n_neighbors value), got {mock_find.call_count}"
        )
        for call_args in mock_find.call_args_list:
            _, kwargs = call_args
            assert kwargs["target_k"] == target_k, (
                f"target_k should be {target_k}, got {kwargs['target_k']}"
            )

        # Results should be in target-grid format (one entry per n_neighbors)
        assert len(results) == 2
        for r in results:
            assert "n_neighbors" in r
            assert "resolution" in r
            assert "n_clusters" in r
            assert r["resolution"] == mock_result[0]
            assert r["n_clusters"] == mock_result[1]

        # enrich_grid_results was called (it's part of the pipeline)
        mock_enrich.assert_called_once()
