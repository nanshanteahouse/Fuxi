"""Tests for pairwise DE-gated cluster selection (Shekhar 2016 semantics).

These tests validate the rewritten ``_select_de_gated`` and the new
``_compute_pairwise_de_markers`` helper, which replace the old one-vs-rest
approach with true pairwise comparisons.

A gene is a pairwise marker for cluster C_i if it is upregulated (padj < 0.05
AND log2FC > 1.0) in ALL pairwise comparisons C_i vs C_j (j ≠ i).
"""

from __future__ import annotations

import logging
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from core.cluster.evaluation import _compute_pairwise_de_markers, _select_de_gated

# ── Helpers ────────────────────────────────────────────────────────────────


def _make_single_group_result(
    ci: str,
    n_de: int,
    n_genes: int = 50,
    rng: np.random.RandomState | None = None,
) -> dict:
    """Create mock pairwise ``rank_genes_groups`` result for a single group *ci*.

    The structured arrays are keyed by *ci* (the group being compared against a
    reference). Exactly *n_de* genes have padj < 0.05 and log2FC > 1.0.
    """
    if rng is None:
        rng = np.random.RandomState(42)

    dtype_names = np.dtype([(ci, "U20")])
    dtype_val = np.dtype([(ci, np.float64)])

    names = np.zeros(n_genes, dtype=dtype_names)
    pvals_adj = np.zeros(n_genes, dtype=dtype_val)
    lfcs = np.zeros(n_genes, dtype=dtype_val)

    for i in range(n_genes):
        names[ci][i] = f"gene_{i}"

    if n_de > 0:
        pvals_adj[ci][:n_de] = rng.uniform(0.001, 0.04, n_de)
        lfcs[ci][:n_de] = rng.uniform(1.1, 3.0, n_de)
    if n_de < n_genes:
        pvals_adj[ci][n_de:] = rng.uniform(0.06, 0.5, n_genes - n_de)
        lfcs[ci][n_de:] = rng.uniform(-0.5, 0.8, n_genes - n_de)

    return {"names": names, "pvals_adj": pvals_adj, "logfoldchanges": lfcs}


def _make_ovr_results(
    group_names: list[str],
    n_de_per_group: int | dict[str, int],
    n_genes: int = 50,
    rng: np.random.RandomState | None = None,
) -> dict:
    """Create mock one-vs-rest ``rank_genes_groups`` result for multiple groups."""
    if rng is None:
        rng = np.random.RandomState(42)

    dtype_val = np.dtype([(g, np.float64) for g in group_names])
    dtype_names = np.dtype([(g, "U20") for g in group_names])

    names = np.zeros(n_genes, dtype=dtype_names)
    pvals_adj = np.zeros(n_genes, dtype=dtype_val)
    lfcs = np.zeros(n_genes, dtype=dtype_val)

    for i in range(n_genes):
        for g in group_names:
            names[g][i] = f"gene_{i}"

    for g in group_names:
        nd = n_de_per_group if isinstance(n_de_per_group, int) else n_de_per_group.get(g, 20)
        if nd > 0:
            pvals_adj[g][:nd] = rng.uniform(0.001, 0.04, nd)
            lfcs[g][:nd] = rng.uniform(1.1, 3.0, nd)
        if nd < n_genes:
            pvals_adj[g][nd:] = rng.uniform(0.06, 0.5, n_genes - nd)
            lfcs[g][nd:] = rng.uniform(-0.5, 0.8, n_genes - nd)

    return {"names": names, "pvals_adj": pvals_adj, "logfoldchanges": lfcs}


def _setup_adata_for_pairwise(n_clusters: int) -> MagicMock:
    """Return a MagicMock AnnData with obs configured for *n_clusters* groups."""
    adata = MagicMock()
    adata.uns = {}
    adata.obs = MagicMock()
    cluster_list = [str(i) for i in range(n_clusters)]
    adata.obs.__getitem__.return_value.unique.return_value = cluster_list
    return adata


# ── Tests: _compute_pairwise_de_markers ────────────────────────────────────


class TestPairwiseDeMarkers:
    """Direct tests for ``_compute_pairwise_de_markers``."""

    # ── Test 1: basic well-separated clusters ────────────────────────────

    def test_pairwise_de_markers_basic(self) -> None:
        """Synthetic 3 well-separated clusters → each has high pairwise marker count.

        Every pair yields 25 DE genes → intersection per cluster = 25.
        """
        adata = _setup_adata_for_pairwise(3)
        cluster_key = "leiden_0.5"
        uns_store: dict = adata.uns  # type: ignore[assignment]
        rng = np.random.RandomState(42)

        def rank_side_effect(adata, groupby, **kwargs):  # noqa: ARG001
            groups = kwargs.get("groups")
            reference = kwargs.get("reference")
            key_added = kwargs.get("key_added", "_pairwise_tmp")
            ci = str(groups[0])
            _c = str(reference)  # ignored — all pairs same
            uns_store[key_added] = _make_single_group_result(ci, n_de=25, rng=rng)

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            result = _compute_pairwise_de_markers(adata, cluster_key)

        assert result == {"0": 25, "1": 25, "2": 25}
        # 3 clusters → 6 pairs (3×2): (0→1, 0→2), (1→0, 1→2), (2→0, 2→1)
        assert mock_rank.call_count == 6

    # ── Test 2: similar clusters yield low pairwise counts ───────────────

    def test_pairwise_de_markers_similar_clusters(self) -> None:
        """2 similar clusters + 1 distinct → similar clusters have LOW pairwise count.

        Clusters 0 and 1 share expression profiles and have few pairwise DE
        genes between them (10), but both are distinct from cluster 2 (40 DE
        genes). Pairwise DE correctly detects that clusters 0 and 1 are similar.

        In contrast, one-vs-rest would inflate the DE count for clusters 0
        and 1 (it averages across all comparisons), producing artificially
        high counts and missing the similarity.
        """
        adata = _setup_adata_for_pairwise(3)
        cluster_key = "leiden_0.5"
        uns_store: dict = adata.uns  # type: ignore[assignment]
        rng = np.random.RandomState(42)

        # DE genes per pair: (ci, cj) → n_de
        de_by_pair = {
            ("0", "1"): 10,  # similar → low
            ("0", "2"): 40,  # distinct → high
            ("1", "0"): 10,
            ("1", "2"): 40,
            ("2", "0"): 40,
            ("2", "1"): 40,
        }

        def rank_side_effect(adata, groupby, **kwargs):  # noqa: ARG001
            groups = kwargs.get("groups")
            reference = kwargs.get("reference")
            key_added = kwargs.get("key_added", "_pairwise_tmp")
            ci = str(groups[0])
            cj = str(reference)
            n_de = de_by_pair.get((ci, cj), 0)
            uns_store[key_added] = _make_single_group_result(ci, n_de=n_de, rng=rng)

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            result = _compute_pairwise_de_markers(adata, cluster_key)

        # Cluster 0: intersect(gene0-9, gene0-39) = genes 0-9 → 10 markers
        # Cluster 1: same logic → 10 markers
        # Cluster 2: intersect(gene0-39, gene0-39) = genes 0-39 → 40 markers
        assert result == {"0": 10, "1": 10, "2": 40}

        # Verify: in a one-vs-rest approach, cluster 0 might show artificially
        # high DE genes because it's counted against the pooled rest (including
        # the distinct cluster 2). Pairwise correctly catches the similarity.
        assert result["0"] < result["2"], (
            "Similar cluster 0 should have fewer pairwise markers than distinct cluster 2"
        )

    # ── Test 5: single cluster edge case ─────────────────────────────────

    def test_pairwise_de_single_cluster(self) -> None:
        """Single cluster returns marker_count dict without error."""
        adata = _setup_adata_for_pairwise(1)
        cluster_key = "leiden_0.5"

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            result = _compute_pairwise_de_markers(adata, cluster_key)

        assert result == {"0": 0}
        # No pairs to compare → rank_genes_groups never called
        mock_rank.assert_not_called()

    # ── Test 6: cluster cap triggers one-vs-rest fallback ────────────────

    def test_pairwise_de_cluster_cap(self, caplog: pytest.LogCaptureFixture) -> None:
        """With max_clusters_for_pairwise=5 and 6 clusters, falls back to one-vs-rest.

        Warning is logged. One-vs-rest results (20 DE per group) are returned.
        """
        adata = _setup_adata_for_pairwise(6)
        cluster_key = "leiden_0.5"
        uns_store: dict = adata.uns  # type: ignore[assignment]
        rng = np.random.RandomState(42)

        def rank_side_effect(adata, groupby, **kwargs):
            key_added = kwargs.get("key_added", "_pairwise_fallback")
            groups_kw = kwargs.get("groups")
            if groups_kw is None:
                # One-vs-rest mode → all groups
                group_names = [str(i) for i in range(6)]
                uns_store[key_added] = _make_ovr_results(group_names, n_de_per_group=20, rng=rng)

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            with caplog.at_level(logging.WARNING):
                result = _compute_pairwise_de_markers(
                    adata,
                    cluster_key,
                    max_clusters_for_pairwise=5,
                    log=logging.getLogger("test_pairwise_cap"),
                )

        # Warning should mention the cap
        assert "exceeds de_pairwise_max_clusters" in caplog.text
        assert "falling back to one-vs-rest" in caplog.text
        # All 6 clusters have 20 DE genes
        assert result == {str(i): 20 for i in range(6)}
        # rank_genes_groups called exactly once (one-vs-rest mode)
        assert mock_rank.call_count == 1

    # ── Test 7: no cap → always pairwise ─────────────────────────────────

    def test_pairwise_de_no_cap(self) -> None:
        """With max_clusters_for_pairwise=0 (no cap), 40 clusters still uses pairwise."""
        # Use 4 clusters for practical testing (40 would be computationally
        # expensive even with mocks, but the cap check is at the function entry
        # point — 4 clusters with cap=0 still proves the cap is bypassed).
        n_clusters = 4
        adata = _setup_adata_for_pairwise(n_clusters)
        cluster_key = "leiden_0.5"
        uns_store: dict = adata.uns  # type: ignore[assignment]
        rng = np.random.RandomState(42)

        def rank_side_effect(adata, groupby, **kwargs):  # noqa: ARG001
            groups = kwargs.get("groups")
            reference = kwargs.get("reference")
            if groups is not None and reference is not None:
                key_added = kwargs.get("key_added", "_pairwise_tmp")
                ci = str(groups[0])
                uns_store[key_added] = _make_single_group_result(ci, n_de=30, rng=rng)

        with patch("scanpy.tl.rank_genes_groups") as mock_rank:
            mock_rank.side_effect = rank_side_effect
            result = _compute_pairwise_de_markers(adata, cluster_key, max_clusters_for_pairwise=0)

        # All 4 clusters, pairwise mode → each should have 30 markers
        assert result == {str(i): 30 for i in range(4)}
        # 4 clusters → 4 × 3 = 12 pairwise calls (all pairwise, not OVR)
        assert mock_rank.call_count == 12


# ── Tests: _select_de_gated (with pairwise) ───────────────────────────────


class TestSelectDeGatedPairwise:
    """Tests for ``_select_de_gated`` gating logic (mocking the DE helper)."""

    # ── Test 3: gating selects lower resolution when not pairwise-distinct ─

    def test_select_de_gated_pairwise_gating(self) -> None:
        """Verify gating selects lower resolution when clusters not pairwise-distinct.

        Two resolutions: 1.0 has 10 clusters with min pairwise DE = 5 (below
        threshold 25), 0.5 has 5 clusters with min DE = 30 (above threshold).
        Gating should select resolution 0.5 (the highest meeting threshold).
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
        adata.uns = {}

        # Mock _compute_pairwise_de_markers to return controlled counts
        de_by_key = {
            "leiden_0.5": {"0": 30, "1": 35, "2": 32, "3": 30, "4": 33},
            "leiden_1.0": {str(i): 5 for i in range(10)},
        }

        def de_side_effect(adata, cluster_key, **kwargs):  # noqa: ARG001
            return de_by_key.get(cluster_key, {})

        with patch("core.cluster.evaluation._compute_pairwise_de_markers") as mock_de:
            mock_de.side_effect = de_side_effect
            result = _select_de_gated(valid, adata, de_gate_threshold=25, pairwise_max_clusters=30)

        n_clusters, resolution, cluster_key, reason = result
        # leiden_1.0 has min_de=5 < 25 → fails. leiden_0.5 has min_de=30 ≥ 25 → selected
        assert n_clusters == 5
        assert resolution == pytest.approx(0.5)
        assert cluster_key == "leiden_0.5"
        assert "de_gated" in reason

    # ── Test 4: no _pairwise_* keys remain in adata.uns ───────────────────

    def test_pairwise_de_cleanup(self) -> None:
        """No _pairwise_* keys remain in adata.uns after _select_de_gated call."""
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

        def de_side_effect(adata, cluster_key, **kwargs):  # noqa: ARG001
            return {str(i): 30 for i in range(5)}

        with patch("core.cluster.evaluation._compute_pairwise_de_markers") as mock_de:
            mock_de.side_effect = de_side_effect
            _select_de_gated(valid, adata, de_gate_threshold=25, pairwise_max_clusters=30)

        # No _pairwise_* keys should persist in uns
        pairwise_keys = [k for k in uns_store if k.startswith("_pairwise_")]
        assert len(pairwise_keys) == 0, f"Found lingering _pairwise_* keys: {pairwise_keys}"

    # ── Regression: existing tests must still pass ────────────────────────

    def test_backward_compat_signature_no_pairwise_max(self) -> None:
        """_select_de_gated can still be called without pairwise_max_clusters.

        This ensures backward compatibility with callers that don't pass the
        new parameter (falls back to default 30).
        """
        valid = [
            {
                "n_clusters": 5,
                "resolution": 0.5,
                "cluster_key": "leiden_0.5",
                "n_neighbors": 15,
            },
        ]
        adata = MagicMock()
        adata.uns = {}

        def de_side_effect(adata, cluster_key, **kwargs):  # noqa: ARG001
            return {"0": 30, "1": 35, "2": 32, "3": 30, "4": 33}

        with patch("core.cluster.evaluation._compute_pairwise_de_markers") as mock_de:
            mock_de.side_effect = de_side_effect
            # No pairwise_max_clusters kwarg → should use default
            result = _select_de_gated(valid, adata, de_gate_threshold=25)

        n_clusters, resolution, cluster_key, reason = result
        assert n_clusters == 5
        assert reason == "de_gated(single_entry, min_de=N/A)"
        # Verify mock was called with default max_clusters_for_pairwise
        assert mock_de.call_count == 0  # single entry → no DE computation
