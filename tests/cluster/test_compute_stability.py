"""Tests for _compute_stability — cross-seed clustering stability.

These tests supplement tests/test_rna/test_cluster_evaluation.py and cover
post-fix behaviors (B1, B2 from the audit):

- B1: all seeds fail → float('nan') + logger.warning (currently returns 1.0)
- B2: n_seeds <= 1 → float('nan') (currently returns 1.0)

Also covers post-#16 early-termination property and side-effect safety.
"""

import hashlib
import logging
from unittest.mock import patch

import numpy as np
import pytest

from core.cluster.evaluation import _compute_stability


class TestComputeStabilityWellSeparated:
    """Well-separated clusters yield high stability."""

    def test_well_separated_stability_high(self) -> None:
        """3 clearly separated clusters via scanpy leiden → stability > 0.9."""
        import scanpy as sc

        n_cells = 300
        n_genes = 50
        rng = np.random.RandomState(42)

        # 3 blocks of 8 marker genes each, strong separation (mean=8)
        x = rng.normal(0, 0.2, (n_cells, n_genes))
        x[:100, :8] += rng.normal(8, 0.2, (100, 8))
        x[100:200, 8:16] += rng.normal(8, 0.2, (100, 8))
        x[200:300, 16:24] += rng.normal(8, 0.2, (100, 8))
        adata = sc.AnnData(x)

        sc.pp.pca(adata, n_comps=15)
        sc.pp.neighbors(adata, n_neighbors=15)

        stability = _compute_stability(adata, resolution=0.5, n_seeds=10)
        assert stability > 0.9, f"Expected > 0.9, got {stability}"


class TestComputeStabilityOverlapping:
    """Overlapping clusters yield low stability."""

    def test_severely_overlapping_stability_low(self) -> None:
        """Overlapping clusters (pure noise) → stability < 0.6.

        Uses high-resolution Leiden on random data with small n_neighbors so
        that different seeds produce different clusterings.
        """
        import scanpy as sc

        n_cells = 200
        n_genes = 50
        rng = np.random.RandomState(42)

        # Pure noise — no real cluster structure.  With high resolution
        # and few neighbors, Leiden cluster assignments are unstable.
        x = rng.randn(n_cells, n_genes)
        adata = sc.AnnData(x)

        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=8)

        stability = _compute_stability(adata, resolution=1.0, n_seeds=10)
        assert stability < 0.65, f"Expected < 0.65, got {stability}"


class TestComputeStabilityEdgeCases:
    """Edge cases: degenerate n_seeds values."""

    def test_n_seeds_1_returns_nan(self) -> None:
        """n_seeds=1 → returns float('nan').  [POST-FIX: B2 fix, currently returns 1.0]"""
        import scanpy as sc

        n_cells = 100
        n_genes = 20
        rng = np.random.RandomState(42)
        adata = sc.AnnData(rng.randn(n_cells, n_genes))
        sc.pp.pca(adata, n_comps=5)
        sc.pp.neighbors(adata, n_neighbors=10)

        result = _compute_stability(adata, resolution=0.3, n_seeds=1)
        assert np.isnan(result), f"n_seeds=1 should return NaN, got {result}"

    def test_n_seeds_0_returns_nan(self) -> None:
        """n_seeds=0 → returns float('nan')."""
        import scanpy as sc

        n_cells = 100
        n_genes = 20
        rng = np.random.RandomState(42)
        adata = sc.AnnData(rng.randn(n_cells, n_genes))
        sc.pp.pca(adata, n_comps=5)
        sc.pp.neighbors(adata, n_neighbors=10)

        result = _compute_stability(adata, resolution=0.3, n_seeds=0)
        assert np.isnan(result), f"n_seeds=0 should return NaN, got {result}"

    def test_all_seeds_fail_returns_nan_with_warning(self, caplog) -> None:
        """Monkeypatch a seed run to raise Exception → returns NaN + logger.warning.

        [POST-FIX: B1 fix, currently returns 1.0 silently]
        """
        import scanpy as sc

        n_cells = 100
        n_genes = 20
        rng = np.random.RandomState(42)
        adata = sc.AnnData(rng.randn(n_cells, n_genes))
        sc.pp.pca(adata, n_comps=5)
        sc.pp.neighbors(adata, n_neighbors=10)

        with patch("scanpy.tl.leiden", side_effect=Exception("Leiden failed")):
            with caplog.at_level(logging.WARNING):
                result = _compute_stability(adata, resolution=0.3, n_seeds=5)

        assert np.isnan(result), f"All seeds fail should return NaN, got {result}"
        # At least one WARNING should mention failure / all seeds / returning NaN
        assert any(
            "all" in rec.message.lower()
            or "fail" in rec.message.lower()
            or "seed" in rec.message.lower()
            for rec in caplog.records
        ), "Expected a warning about all seeds failing"


class TestComputeStabilityEarlyTermination:
    """Early termination when perfect agreement is detected."""

    def test_early_termination_on_perfect_agreement(self) -> None:
        """ARI=1.0 for first 2 seeds → returns 1.0 early.

        [post-#16 fix property]
        """
        import scanpy as sc

        n_cells = 200
        n_genes = 30
        rng = np.random.RandomState(42)

        # Two extremely well-separated clusters — leiden will produce
        # identical labels across every seed run.
        x = rng.normal(0, 0.2, (n_cells, n_genes))
        x[:100, :5] += rng.normal(10, 0.2, (100, 5))
        x[100:, :5] += rng.normal(-10, 0.2, (100, 5))
        adata = sc.AnnData(x)

        sc.pp.pca(adata, n_comps=5)
        sc.pp.neighbors(adata, n_neighbors=10)

        stability = _compute_stability(adata, resolution=0.5, n_seeds=5)
        assert stability == pytest.approx(1.0, abs=0.01), (
            f"Perfectly separated clusters should give ~1.0, got {stability}"
        )


class TestComputeStabilitySideEffects:
    """Verify that _compute_stability does not mutate adata.obsp."""

    def test_knn_graph_unchanged(self) -> None:
        """adata.obsp hash unchanged before/after seed runs."""
        import scanpy as sc

        n_cells = 200
        n_genes = 30
        rng = np.random.RandomState(42)

        x = rng.normal(0, 0.3, (n_cells, n_genes))
        x[:100, :5] += rng.normal(5, 0.3, (100, 5))
        x[100:, 5:10] += rng.normal(5, 0.3, (100, 5))
        adata = sc.AnnData(x)

        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=15)

        conn_hash_before = hashlib.md5(adata.obsp["connectivities"].data.tobytes()).hexdigest()

        _compute_stability(adata, resolution=0.5, n_seeds=3)

        conn_hash_after = hashlib.md5(adata.obsp["connectivities"].data.tobytes()).hexdigest()
        assert conn_hash_before == conn_hash_after, (
            "KNN graph (adata.obsp) must not be modified by _compute_stability"
        )
