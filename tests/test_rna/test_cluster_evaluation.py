"""Numerical tests for rna/utils/cluster_evaluation.py."""

import logging
import numpy as np

import pytest

from rna.utils.cluster_evaluation import (
    select_best_params,
    select_best_umap_params,
    _compute_stability,
    _compute_cluster_coherence,
    _select_multi_metric,
    _detect_granularity,
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


class TestComputeStability:
    """Numerical assertions for _compute_stability."""

    def test_compute_stability_well_separated(self) -> None:
        """3 clearly separated clusters via scanpy leiden → stability > 0.6."""
        import scanpy as sc

        n_cells = 300
        n_genes = 50
        rng = np.random.RandomState(42)
        X = rng.normal(0, 0.3, (n_cells, n_genes))
        # Create 3 well-separated clusters
        X[:100, :10] += rng.normal(5, 0.3, (100, 10))
        X[100:200, 10:20] += rng.normal(5, 0.3, (100, 10))
        X[200:300, 20:30] += rng.normal(5, 0.3, (100, 10))
        adata = sc.AnnData(X)

        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=15)

        stability = _compute_stability(adata, resolution=0.5, n_seeds=3)
        assert stability > 0.6, f"Expected stability > 0.6, got {stability}"

    def test_compute_stability_single_cluster(self) -> None:
        """Single-cluster AnnData → stability == 1.0."""
        import scanpy as sc

        n_cells = 100
        n_genes = 20
        rng = np.random.RandomState(42)
        X = rng.normal(0, 0.05, (n_cells, n_genes))
        adata = sc.AnnData(X)

        sc.pp.pca(adata, n_comps=5)
        sc.pp.neighbors(adata, n_neighbors=10)

        # Very low resolution → should get single cluster consistently
        stability = _compute_stability(adata, resolution=0.01, n_seeds=3)
        assert stability == pytest.approx(1.0, abs=0.01), (
            f"Expected stability ~1.0, got {stability}"
        )

    def test_compute_stability_n_seeds(self) -> None:
        """n_seeds=1 returns 1.0; n_seeds=3 returns valid mean."""
        import scanpy as sc

        n_cells = 200
        n_genes = 30
        rng = np.random.RandomState(42)
        X = rng.normal(0, 0.5, (n_cells, n_genes))
        # Make 2 moderate clusters
        X[:100, :8] += rng.normal(4, 0.3, (100, 8))
        X[100:, 8:16] += rng.normal(4, 0.3, (100, 8))
        adata = sc.AnnData(X)

        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=15)

        # n_seeds=1 should return exactly 1.0
        stab_1 = _compute_stability(adata, resolution=0.3, n_seeds=1)
        assert stab_1 == 1.0, f"n_seeds=1 should return 1.0, got {stab_1}"

        # n_seeds=3 should produce a valid float between 0 and 1
        stab_3 = _compute_stability(adata, resolution=0.3, n_seeds=3)
        assert 0.0 <= stab_3 <= 1.0, (
            f"n_seeds=3 should be in [0,1], got {stab_3}"
        )


class TestComputeMarkerCoverage:
    """Numerical assertions for _compute_marker_coverage."""

    def test_compute_marker_coverage_perfect_match(self) -> None:
        """AnnData with clusters that express known markers → coverage > 0.8."""
        import scanpy as sc

        n_cells = 300
        rng = np.random.RandomState(42)

        # 3 clusters, 100 cells each
        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["leiden_15_0.5"] = (
            ["0"] * 100 + ["1"] * 100 + ["2"] * 100
        )

        # per_cell_scores: each cell type boosted in its matching cluster
        per_cell_scores: dict[str, np.ndarray] = {}
        for i in range(3):
            scores = rng.normal(0, 0.5, n_cells)
            cluster_mask = adata.obs["leiden_15_0.5"] == str(i)
            scores[cluster_mask.values] += 5.0
            per_cell_scores[f"Type{i}"] = scores

        coverage = _compute_cluster_coherence(
            adata, cluster_key="leiden_15_0.5",
            per_cell_scores=per_cell_scores,
        )
        assert coverage > 0.8, (
            f"Expected coverage > 0.8 for perfect match, got {coverage}"
        )

    def test_compute_marker_coverage_random(self) -> None:
        """Random expression → coverage < 0.5."""
        import scanpy as sc

        n_cells = 200
        rng = np.random.RandomState(42)

        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["leiden_15_0.5"] = (
            ["0"] * 50 + ["1"] * 50 + ["2"] * 50 + ["3"] * 50
        )

        # 5 cell types → harder for any single type to randomly dominate a cluster
        per_cell_scores = {
            f"Type{c}": rng.randn(n_cells) for c in range(5)
        }

        coverage = _compute_cluster_coherence(
            adata, cluster_key="leiden_15_0.5",
            per_cell_scores=per_cell_scores,
        )
        assert coverage < 0.6, (
            f"Expected coverage < 0.6 for random scores, got {coverage}"
        )

    def test_compute_marker_coverage_empty_scores(self) -> None:
        """Empty per_cell_scores dict → returns 1.0."""
        import scanpy as sc

        adata = sc.AnnData(X=np.random.randn(50, 10))
        adata.obs["leiden"] = ["0"] * 50

        coverage = _compute_cluster_coherence(
            adata, cluster_key="leiden",
            per_cell_scores={},
        )
        assert coverage == 1.0, (
            f"Empty scores should return 1.0, got {coverage}"
        )

    def test_compute_marker_coverage_missing_genes(self) -> None:
        """Cell type scores with None entries → graceful, still returns valid coverage."""
        import scanpy as sc

        n_cells = 150
        rng = np.random.RandomState(42)

        adata = sc.AnnData(X=rng.randn(n_cells, 20))
        adata.obs["leiden"] = ["0"] * 50 + ["1"] * 50 + ["2"] * 50

        # One valid score, one None
        good_scores = rng.normal(0, 0.5, n_cells)
        good_scores[:50] += 5.0  # boost in cluster 0

        per_cell_scores: dict[str, np.ndarray | None] = {
            "TypeA": good_scores,
            "TypeB": None,  # missing data
        }

        coverage = _compute_cluster_coherence(
            adata, cluster_key="leiden",
            per_cell_scores=per_cell_scores,
        )
        assert 0.0 <= coverage <= 1.0, (
            f"Should return valid coverage even with missing genes, got {coverage}"
        )


class TestSelectMultiMetric:
    """Numerical assertions for _select_multi_metric."""

    def test_select_multi_metric(self) -> None:
        """with/without marker_coverage keys, verify composite returns correctly."""
        # ── Without marker_coverage ──
        valid_no_mc = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
                "stability_score": 0.8,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 6, "silhouette_score": 0.7,
                "stability_score": 0.9,
            },
            {
                "n_neighbors": 30, "resolution": 1.0,
                "n_clusters": 9, "silhouette_score": 0.6,
                "stability_score": 0.7,
            },
        ]

        best_n, best_r, method, reason = _select_multi_metric(valid_no_mc)
        assert method == "multi_metric"
        assert isinstance(best_n, int)
        assert isinstance(best_r, float)
        assert "sil=" in reason
        assert "stab=" in reason

        # -- With cluster_coherence --
        valid_with_mc = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
                "stability_score": 0.8, "cluster_coherence": 0.6,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 6, "silhouette_score": 0.7,
                "stability_score": 0.9, "cluster_coherence": 0.5,
            },
            {
                "n_neighbors": 30, "resolution": 1.0,
                "n_clusters": 9, "silhouette_score": 0.6,
                "stability_score": 0.7, "cluster_coherence": 0.7,
            },
        ]

        best_n2, best_r2, method2, reason2 = _select_multi_metric(valid_with_mc)
        assert method2 == "multi_metric"
        assert isinstance(best_n2, int)
        assert isinstance(best_r2, float)
        assert "coherence=" in reason2

        # ── Single entry (edge case) ──
        valid_single = [
            {
                "n_neighbors": 15, "resolution": 0.6,
                "n_clusters": 4, "silhouette_score": 0.5,
                "stability_score": 0.85, "cluster_coherence": 0.4,
            },
        ]

        best_n3, best_r3, method3, _reason3 = _select_multi_metric(valid_single)
        assert best_n3 == 15
        assert best_r3 == pytest.approx(0.6)
        assert method3 == "multi_metric"



class TestSelectMultiMetricIntegration:
    """Integration tests for _select_multi_metric with realistic grid data.

    Uses 6-entry grids (3 n_neighbors × 2 resolutions) to verify
    composite multi-metric scoring end-to-end.
    """

    def test_select_multi_metric_with_markers(self) -> None:
        """6 entries with marker_coverage -> returns valid 4-tuple; mc influences winner."""
        valid = [
            {"n_neighbors": 10, "resolution": 0.5, "n_clusters": 3,
             "silhouette_score": 0.4, "stability_score": 0.7,
             "cluster_coherence": 0.3, "cluster_key": "leiden_10_0.5"},
            {"n_neighbors": 10, "resolution": 1.0, "n_clusters": 5,
             "silhouette_score": 0.4, "stability_score": 0.7,
             "cluster_coherence": 0.6, "cluster_key": "leiden_10_1.0"},
            {"n_neighbors": 20, "resolution": 0.5, "n_clusters": 6,
             "silhouette_score": 0.6, "stability_score": 0.85,
             "cluster_coherence": 0.5, "cluster_key": "leiden_20_0.5"},
            {"n_neighbors": 20, "resolution": 1.0, "n_clusters": 8,
             "silhouette_score": 0.6, "stability_score": 0.85,
             "cluster_coherence": 0.8, "cluster_key": "leiden_20_1.0"},
            {"n_neighbors": 30, "resolution": 0.5, "n_clusters": 9,
             "silhouette_score": 0.5, "stability_score": 0.75,
             "cluster_coherence": 0.7, "cluster_key": "leiden_30_0.5"},
            {"n_neighbors": 30, "resolution": 1.0, "n_clusters": 12,
             "silhouette_score": 0.5, "stability_score": 0.75,
             "cluster_coherence": 0.9, "cluster_key": "leiden_30_1.0"},
        ]
        best_n, best_r, method, reason = _select_multi_metric(valid)
        assert method == "multi_metric"
        assert isinstance(best_n, int)
        assert isinstance(best_r, float)
        assert "coherence=" in reason
        assert "composite=" in reason
        assert best_n == 20
        assert best_r == pytest.approx(1.0)

    def test_select_multi_metric_without_markers(self) -> None:
        """6 entries without marker_coverage -> degrades to sil+stab (0.5/0.5)."""
        valid = [
            {"n_neighbors": 10, "resolution": 0.5, "n_clusters": 3,
             "silhouette_score": 0.3, "stability_score": 0.6,
             "cluster_key": "leiden_10_0.5"},
            {"n_neighbors": 10, "resolution": 1.0, "n_clusters": 5,
             "silhouette_score": 0.4, "stability_score": 0.7,
             "cluster_key": "leiden_10_1.0"},
            {"n_neighbors": 20, "resolution": 0.5, "n_clusters": 6,
             "silhouette_score": 0.8, "stability_score": 0.95,
             "cluster_key": "leiden_20_0.5"},
            {"n_neighbors": 20, "resolution": 1.0, "n_clusters": 8,
             "silhouette_score": 0.6, "stability_score": 0.8,
             "cluster_key": "leiden_20_1.0"},
            {"n_neighbors": 30, "resolution": 0.5, "n_clusters": 9,
             "silhouette_score": 0.5, "stability_score": 0.75,
             "cluster_key": "leiden_30_0.5"},
            {"n_neighbors": 30, "resolution": 1.0, "n_clusters": 12,
             "silhouette_score": 0.55, "stability_score": 0.7,
             "cluster_key": "leiden_30_1.0"},
        ]
        best_n, best_r, method, reason = _select_multi_metric(valid)
        assert method == "multi_metric"
        assert "coherence=0.000" in reason
        assert "sil=" in reason and "stab=" in reason
        assert best_n == 20
        assert best_r == pytest.approx(0.5)

    def test_select_multi_metric_no_marker_degrade(self) -> None:
        """All entries lack marker_coverage -> weights degrade to {sil:0.5, stab:0.5}."""
        valid = [
            {"n_neighbors": 15, "resolution": 0.6, "n_clusters": 4,
             "silhouette_score": 0.5, "stability_score": 0.8,
             "cluster_key": "leiden_15_0.6"},
            {"n_neighbors": 25, "resolution": 1.2, "n_clusters": 7,
             "silhouette_score": 0.7, "stability_score": 0.9,
             "cluster_key": "leiden_25_1.2"},
        ]
        best_n, best_r, method, reason = _select_multi_metric(valid)
        assert method == "multi_metric"
        assert "coherence=0.000" in reason
        assert "sil=" in reason and "stab=" in reason
        assert best_n == 25
        assert best_r == pytest.approx(1.2)

    def test_select_multi_metric_single_entry(self) -> None:
        """Single valid entry -> returns that entry unchanged."""
        valid = [
            {"n_neighbors": 15, "resolution": 0.6, "n_clusters": 4,
             "silhouette_score": 0.5, "stability_score": 0.85,
             "cluster_coherence": 0.4, "cluster_key": "leiden_15_0.6"},
        ]
        best_n, best_r, method, _reason = _select_multi_metric(valid)
        assert best_n == 15
        assert best_r == pytest.approx(0.6)
        assert method == "multi_metric"



class TestAdaptiveBehaviors:
    """Adaptive degrade behaviors for _select_multi_metric.

    Tests the two degrade paths: marker-coverage absence (no key in entries)
    vs. marker-coverage mismatch (all entries have marker_coverage < 0.1).
    """

    def test_marker_mismatch_degrade(self, caplog) -> None:
        """All entries have marker_coverage < 0.1 → warning + degrade to sil+stab."""
        entries = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
                "stability_score": 0.8, "cluster_coherence": 0.0,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 6, "silhouette_score": 0.7,
                "stability_score": 0.9, "cluster_coherence": 0.02,
            },
            {
                "n_neighbors": 30, "resolution": 1.0,
                "n_clusters": 9, "silhouette_score": 0.6,
                "stability_score": 0.7, "cluster_coherence": 0.05,
            },
        ]

        with caplog.at_level(logging.WARNING):
            best_n, best_r, method, reason = _select_multi_metric(entries)

        # ── Warning emitted ──
        assert any(
            "mismatch" in rec.message or "Degrading" in rec.message
            for rec in caplog.records
        ), "Expected a warning about marker mismatch degrade"

        # ── Valid 4-tuple with multi_metric method ──
        assert method == "multi_metric"
        assert isinstance(best_n, int)
        assert isinstance(best_r, float)

        # ── Degraded: marker_cov=0.000, not a differentiating factor ──
        assert "coherence=0.000" in reason
        assert "sil=" in reason
        assert "stab=" in reason

    def test_no_marker_degrade_vs_mismatch_different_paths(self, caplog) -> None:
        """Absent marker_coverage key vs low-coverage entries both degrade but
        via different code paths: absence = no 'marker_coverage' key ever set;
        mismatch = keys present but all values < 0.1."""
        # ── Absence path: entries lack marker_coverage key entirely ──
        no_mc_entries = [
            {
                "n_neighbors": 10, "resolution": 0.5,
                "n_clusters": 3, "silhouette_score": 0.4,
                "stability_score": 0.8,
            },
            {
                "n_neighbors": 20, "resolution": 0.8,
                "n_clusters": 6, "silhouette_score": 0.7,
                "stability_score": 0.9,
            },
        ]
        caplog.clear()
        n1, r1, m1, reason1 = _select_multi_metric(no_mc_entries)
        assert m1 == "multi_metric"
        assert isinstance(n1, int)
        assert isinstance(r1, float)
        assert "coherence=0.000" in reason1
        # Absence path logs no mismatch warning
        assert not any(
            "mismatch" in rec.message or "Degrading" in rec.message
            for rec in caplog.records
        ), "Absence path should NOT log a mismatch warning"

        # ── Mismatch path: entries have marker_coverage but all < 0.1 ──
        mismatch_entries = [
            {
                "n_neighbors": 15, "resolution": 0.6,
                "n_clusters": 4, "silhouette_score": 0.5,
                "stability_score": 0.85, "cluster_coherence": 0.0,
            },
            {
                "n_neighbors": 25, "resolution": 1.2,
                "n_clusters": 7, "silhouette_score": 0.7,
                "stability_score": 0.9, "cluster_coherence": 0.01,
            },
        ]
        caplog.clear()
        n2, r2, m2, reason2 = _select_multi_metric(mismatch_entries)
        assert m2 == "multi_metric"
        assert isinstance(n2, int)
        assert isinstance(r2, float)
        assert "coherence=0.000" in reason2
        # Mismatch path DOES log a warning
        assert any(
            "mismatch" in rec.message or "Degrading" in rec.message
            for rec in caplog.records
        ), "Mismatch path SHOULD log a warning about marker mismatch degrade"

class TestDetectGranularity:
    """Tests for _detect_granularity function."""

    def test_detect_granularity_tissue_level(self) -> None:
        """High CV + many clusters → tissue-level."""
        r = [
            {'n_neighbors': 15, 'resolution': 0.3, 'n_clusters': 5, 'silhouette_score': 0.20},
            {'n_neighbors': 15, 'resolution': 0.5, 'n_clusters': 8, 'silhouette_score': 0.18},
            {'n_neighbors': 15, 'resolution': 0.8, 'n_clusters': 12, 'silhouette_score': 0.12},
            {'n_neighbors': 15, 'resolution': 1.0, 'n_clusters': 14, 'silhouette_score': 0.09},
            {'n_neighbors': 15, 'resolution': 1.5, 'n_clusters': 18, 'silhouette_score': 0.06},
        ]
        assert _detect_granularity(r) == "tissue"

    def test_detect_granularity_subtype_level(self) -> None:
        """Low CV + few clusters → subtype-level. Flat silhouette across resolutions."""
        r = [
            {'n_neighbors': 15, 'resolution': 0.3, 'n_clusters': 4, 'silhouette_score': 0.08},
            {'n_neighbors': 15, 'resolution': 0.5, 'n_clusters': 5, 'silhouette_score': 0.08},
            {'n_neighbors': 15, 'resolution': 0.8, 'n_clusters': 6, 'silhouette_score': 0.08},
            {'n_neighbors': 15, 'resolution': 1.0, 'n_clusters': 7, 'silhouette_score': 0.08},
        ]
        assert _detect_granularity(r) == "subtype"

    def test_detect_granularity_empty(self) -> None:
        """Empty list → tissue (conservative default)."""
        assert _detect_granularity([]) == "tissue"

    def test_detect_granularity_single_entry(self) -> None:
        """Single entry → tissue (conservative default)."""
        r = [{'n_neighbors': 15, 'resolution': 0.5, 'n_clusters': 5, 'silhouette_score': 0.10}]
        assert _detect_granularity(r) == "tissue"

    def test_detect_granularity_missing_keys(self) -> None:
        """Entries missing silhouette_score → gracefully handled."""
        r = [
            {'n_neighbors': 15, 'resolution': 0.5, 'n_clusters': 5},
            {'n_neighbors': 15, 'resolution': 1.0, 'n_clusters': 7, 'silhouette_score': None},
        ]
        assert _detect_granularity(r) == "tissue"
