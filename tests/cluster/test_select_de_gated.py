"""Post-fix tests for _select_de_gated — D2/D4/D3 bug proofs.

These tests supplement (do NOT replace) ``tests/test_rna/test_cluster_evaluation.py``
and validate behaviors that the current code does not yet support:

- D2: ``adata.uns['rank_genes_groups']`` default slot pollution → post-fix uses
  ``key_added`` and cleans up.
- D4: fallback picks lowest resolution → post-fix picks entry with max DE count.
- D3: cache check is dead code → post-fix verifies ``rank_genes_groups`` is
  called for every entry regardless of cached state.
"""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from core.cluster.evaluation import _select_de_gated

# ── Helper ────────────────────────────────────────────────────────────────────────


def _make_mock_rank_results(
    n_clusters: int,
    n_genes: int = 50,
    n_de_genes: int = 30,
    rng: np.random.RandomState | None = None,
) -> dict:
    """Create mock rank_genes_groups results dict with structured arrays.

    ``padj < 0.05`` AND ``lfc > 1.0`` holds for exactly *n_de_genes* per group.
    """
    if rng is None:
        rng = np.random.RandomState(42)
    group_names = [str(i) for i in range(n_clusters)]
    dtype = np.dtype([(g, np.float64) for g in group_names])

    pvals_adj = np.zeros(n_genes, dtype=dtype)
    logfoldchanges = np.zeros(n_genes, dtype=dtype)

    for g in group_names:
        de_mask = np.arange(n_genes) < n_de_genes
        pvals_adj[g][de_mask] = rng.uniform(0.001, 0.04, n_de_genes)
        pvals_adj[g][~de_mask] = rng.uniform(0.05, 0.5, n_genes - n_de_genes)
        logfoldchanges[g][de_mask] = rng.uniform(1.1, 3.0, n_de_genes)
        logfoldchanges[g][~de_mask] = rng.uniform(-0.5, 0.8, n_genes - n_de_genes)

    return {
        "pvals_adj": pvals_adj,
        "logfoldchanges": logfoldchanges,
    }


# ── Tests ─────────────────────────────────────────────────────────────────────────


class TestSelectDeGatedPostFix:
    """Post-fix behavior tests for _select_de_gated."""

    # ── Test 1: well-separated picks highest resolution ───────────────────────────

    def test_well_separated_selects_highest_resolution(self) -> None:
        """Synthetic 3 entries with increasing DE counts → picks highest n_clusters meeting threshold."""
        valid = [
            {
                "n_clusters": 5,
                "resolution": 0.5,
                "cluster_key": "leiden_0.5",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 10,
                "resolution": 1.0,
                "cluster_key": "leiden_1.0",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 20,
                "resolution": 2.0,
                "cluster_key": "leiden_2.0",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()
        uns_store = {}
        adata.uns = uns_store

        rng = np.random.RandomState(42)
        de_by_key = {"leiden_0.5": 10, "leiden_1.0": 20, "leiden_2.0": 30}
        n_clusters_by_key = {"leiden_0.5": 5, "leiden_1.0": 10, "leiden_2.0": 20}

        def rank_side_effect(adata, groupby, **kwargs):
            nk = n_clusters_by_key.get(groupby, 5)
            nd = de_by_key.get(groupby, 0)
            mock_results = _make_mock_rank_results(nk, n_genes=50, n_de_genes=nd, rng=rng)
            uns_store[f"_de_gated_{groupby}"] = mock_results

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            result = _select_de_gated(valid, adata, de_gate_threshold=25)

        n_clusters, resolution, cluster_key, reason = result
        assert n_clusters == 20
        assert resolution == pytest.approx(2.0)
        assert cluster_key == "leiden_2.0"
        assert "de_gated" in reason

    # ── Test 2: fallback picks max DE (D4 post-fix) ───────────────────────────────

    def test_severely_overlapping_fallback_max_de(self) -> None:
        """All below threshold → fallback picks entry with max DE count, NOT lowest resolution.

        POST-FIX: D4 — current code returns ``candidates[0]`` (lowest resolution),
        post-fix returns entry with highest ``min_pairwise_de``.
        """
        valid = [
            {
                "n_clusters": 5,
                "resolution": 0.5,
                "cluster_key": "leiden_0.5",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 10,
                "resolution": 1.0,
                "cluster_key": "leiden_1.0",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 20,
                "resolution": 2.0,
                "cluster_key": "leiden_2.0",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()
        uns_store = {}
        adata.uns = uns_store

        rng = np.random.RandomState(42)
        # DE counts increase with resolution but ALL below threshold=25
        de_by_key = {"leiden_0.5": 3, "leiden_1.0": 10, "leiden_2.0": 20}
        n_clusters_by_key = {"leiden_0.5": 5, "leiden_1.0": 10, "leiden_2.0": 20}

        def rank_side_effect(adata, groupby, **kwargs):
            nk = n_clusters_by_key.get(groupby, 5)
            nd = de_by_key.get(groupby, 0)
            mock_results = _make_mock_rank_results(nk, n_genes=50, n_de_genes=nd, rng=rng)
            uns_store[f"_de_gated_{groupby}"] = mock_results

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            result = _select_de_gated(valid, adata, de_gate_threshold=25)

        n_clusters, resolution, cluster_key, reason = result
        # D4 post-fix: entry with MAX DE count → leiden_2.0 (DE=20)
        assert n_clusters == 20, (
            f"D4 post-fix: expected max-DE entry (n_clusters=20, DE=20), "
            f"got n_clusters={n_clusters}"
        )
        assert resolution == pytest.approx(2.0)
        assert cluster_key == "leiden_2.0"
        assert "de_gated" in reason

    # ── Test 3: uns not polluted after call (D2 post-fix) ─────────────────────────

    def test_uns_not_polluted_after_call(self) -> None:
        """``adata.uns['rank_genes_groups']`` remains empty after call; temporary keys cleaned up.

        POST-FIX: D2 — current code writes to the default ``'rank_genes_groups'``
        slot; post-fix uses ``key_added=f'_de_gated_{cluster_key}'`` and cleans up
        after iteration.
        """
        valid = [
            {
                "n_clusters": 5,
                "resolution": 0.5,
                "cluster_key": "leiden_0.5",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 10,
                "resolution": 1.0,
                "cluster_key": "leiden_1.0",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()
        uns_store: dict = {}
        adata.uns = uns_store

        rng = np.random.RandomState(42)

        def rank_side_effect(adata, groupby, **kwargs):
            mock_results = _make_mock_rank_results(
                n_clusters=10, n_genes=50, n_de_genes=30, rng=rng
            )
            # Post-fix writes via key_added parameter
            uns_store[f"_de_gated_{groupby}"] = mock_results

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            _select_de_gated(valid, adata, de_gate_threshold=25)

        # D2 post-fix: default 'rank_genes_groups' slot must NOT be polluted
        assert "rank_genes_groups" not in adata.uns, (
            "D2 post-fix: default 'rank_genes_groups' slot should not be written to"
        )

    # ── Test 4: cache check is no-op (D3 post-fix) ────────────────────────────────

    def test_cache_check_is_noop(self) -> None:
        """``rank_genes_groups`` is called for EVERY entry regardless of cache.

        POST-FIX: D3 — current code has a dead cache check that never matches in
        practice (each entry has a different ``cluster_key``). Post-fix removes this
        and always recomputes.
        """
        valid = [
            {
                "n_clusters": 5,
                "resolution": 0.5,
                "cluster_key": "leiden_0.5",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 10,
                "resolution": 1.0,
                "cluster_key": "leiden_1.0",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 20,
                "resolution": 2.0,
                "cluster_key": "leiden_2.0",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()

        # Pre-populate adata.uns with something that LOOKS like a cache hit
        # (matches the first entry's cluster_key). Post-fix should IGNORE it.
        adata.uns = {
            "rank_genes_groups": {
                "params": {"groupby": "leiden_0.5"},
                "pvals_adj": None,
                "logfoldchanges": None,
            }
        }

        rng = np.random.RandomState(42)

        def rank_side_effect(adata, groupby, **kwargs):
            nk = {"leiden_0.5": 5, "leiden_1.0": 10, "leiden_2.0": 20}[groupby]
            mock_results = _make_mock_rank_results(nk, n_genes=50, n_de_genes=30, rng=rng)
            adata.uns[f"_de_gated_{groupby}"] = mock_results

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            _select_de_gated(valid, adata, de_gate_threshold=25)

        # D3 post-fix: called for EVERY entry, cache irrelevant
        assert mock_rank.call_count == 3, (
            f"D3 post-fix: expected rank_genes_groups called 3 times "
            f"(once per entry), got {mock_rank.call_count}"
        )

    # ── Test 5: configurable threshold ────────────────────────────────────────────

    def test_threshold_25_configurable(self) -> None:
        """``de_gate_threshold=10`` works, picks correct entry at lower bar."""
        valid = [
            {
                "n_clusters": 5,
                "resolution": 0.5,
                "cluster_key": "leiden_0.5",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 10,
                "resolution": 1.0,
                "cluster_key": "leiden_1.0",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 20,
                "resolution": 2.0,
                "cluster_key": "leiden_2.0",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()
        uns_store = {}
        adata.uns = uns_store

        rng = np.random.RandomState(42)
        de_by_key = {"leiden_0.5": 5, "leiden_1.0": 12, "leiden_2.0": 30}
        n_clusters_by_key = {"leiden_0.5": 5, "leiden_1.0": 10, "leiden_2.0": 20}

        def rank_side_effect(adata, groupby, **kwargs):
            nk = n_clusters_by_key.get(groupby, 5)
            nd = de_by_key.get(groupby, 0)
            mock_results = _make_mock_rank_results(nk, n_genes=50, n_de_genes=nd, rng=rng)
            uns_store[f"_de_gated_{groupby}"] = mock_results

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            # de_gate_threshold=10 → leiden_1.0 (DE=12) and leiden_2.0 (DE=30) pass
            result = _select_de_gated(valid, adata, de_gate_threshold=10)

        n_clusters, resolution, cluster_key, reason = result
        # Highest n_clusters meeting threshold=10 → leiden_2.0
        assert n_clusters == 20
        assert resolution == pytest.approx(2.0)
        assert cluster_key == "leiden_2.0"
        assert "de_gated" in reason

    # ── Test 6: threshold met at k=10 ─────────────────────────────────────────────

    def test_threshold_met_for_k10(self) -> None:
        """When threshold met at k=10 entry, returns k=10 not k=8."""
        valid = [
            {
                "n_clusters": 8,
                "resolution": 0.8,
                "cluster_key": "leiden_0.8",
                "n_neighbors": 15,
            },
            {
                "n_clusters": 10,
                "resolution": 1.0,
                "cluster_key": "leiden_1.0",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()
        uns_store = {}
        adata.uns = uns_store

        rng = np.random.RandomState(42)
        de_by_key = {"leiden_0.8": 10, "leiden_1.0": 30}
        n_clusters_by_key = {"leiden_0.8": 8, "leiden_1.0": 10}

        def rank_side_effect(adata, groupby, **kwargs):
            nk = n_clusters_by_key.get(groupby, 5)
            nd = de_by_key.get(groupby, 0)
            mock_results = _make_mock_rank_results(nk, n_genes=50, n_de_genes=nd, rng=rng)
            uns_store[f"_de_gated_{groupby}"] = mock_results

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            result = _select_de_gated(valid, adata, de_gate_threshold=25)

        n_clusters, resolution, cluster_key, reason = result
        # k=10 meets threshold (DE=30) and is highest → selected
        assert n_clusters == 10
        assert resolution == pytest.approx(1.0)
        assert cluster_key == "leiden_1.0"
        assert "de_gated" in reason
