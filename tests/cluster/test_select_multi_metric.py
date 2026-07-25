"""Post-fix tests for _select_multi_metric — M4 + weight degrade path proofs.

These tests supplement (do NOT replace) ``tests/test_rna/test_cluster_evaluation.py``
and validate behaviors that the current code does not yet fully support:

- Weight-explicit tests: verify that composite weights in the reason string
  sum to 1.0 across all degrade paths (5-metric, 4-metric, 3-metric).
- Coherence mismatch degrade: verify the degrade shows in the reason string
  with correct weight sum.
- M4: 3-tier resolution recommendation — currently only logged; post-fix
  returns it embedded in the result (reason dict).
"""

import logging
import re

from core.cluster.evaluation import _select_multi_metric

# ── Helpers ────────────────────────────────────────────────────────────────────────


def _parse_weights(reason: str) -> dict[str, float] | None:
    """Parse weight sub-string from the reason string.

    Post-fix format expected: ``weights=sil:0.20,stab:0.20,coh:0.30,split:0.20,kb:0.10``
    """
    m = re.search(r"weights=([\w.:,]+)", reason)
    if not m:
        return None
    parts = m.group(1).split(",")
    weights: dict[str, float] = {}
    for part in parts:
        kv = part.split(":")
        if len(kv) == 2:
            weights[kv[0].strip()] = float(kv[1].strip())
    return weights


def _sum_weights(reason: str) -> float:
    """Parse weights from reason and return their sum (0.0 if not found)."""
    w = _parse_weights(reason)
    if w is None:
        return 0.0
    return sum(w.values())


# ── Tests ─────────────────────────────────────────────────────────────────────────


class TestSelectMultiMetricPostFixWeights:
    """Post-fix weight-explicit tests for _select_multi_metric."""

    # ── Test 1: 5-metric weights sum to 1.0 ───────────────────────────────────────

    def test_5_metric_weights_sum_to_one(self) -> None:
        """Entries with all 5 metrics → composite weights in reason sum to 1.0."""
        valid = [
            {
                "n_neighbors": 10,
                "resolution": 0.5,
                "n_clusters": 3,
                "silhouette_score": 0.3,
                "stability_score": 0.6,
                "cluster_coherence": 0.4,
                "splitting_gain": 1.0,
                "kb_annotatable_rate": 0.2,
            },
            {
                "n_neighbors": 20,
                "resolution": 0.8,
                "n_clusters": 6,
                "silhouette_score": 0.7,
                "stability_score": 0.9,
                "cluster_coherence": 0.8,
                "splitting_gain": 3.0,
                "kb_annotatable_rate": 0.9,
            },
            {
                "n_neighbors": 30,
                "resolution": 1.0,
                "n_clusters": 9,
                "silhouette_score": 0.5,
                "stability_score": 0.7,
                "cluster_coherence": 0.6,
                "splitting_gain": 2.0,
                "kb_annotatable_rate": 0.5,
            },
        ]

        _n, _r, method, reason = _select_multi_metric(valid)

        assert method == "multi_metric"
        weight_sum = _sum_weights(reason)
        assert abs(weight_sum - 1.0) < 0.01, (
            f"POST-FIX: 5-metric weights should sum to 1.0, got {weight_sum:.4f} "
            f"(raw reason: {reason[:200]})"
        )

    # ── Test 2: 4-metric weights (no kb_annotatable_rate) sum to 1.0 ──────────────

    def test_4_metric_weights_sum_to_one(self) -> None:
        """Entries without ``kb_annotatable_rate`` → 4-metric degrade weights sum to 1.0."""
        valid = [
            {
                "n_neighbors": 10,
                "resolution": 0.5,
                "n_clusters": 3,
                "silhouette_score": 0.3,
                "stability_score": 0.6,
                "cluster_coherence": 0.4,
                "splitting_gain": 1.0,
            },
            {
                "n_neighbors": 20,
                "resolution": 0.8,
                "n_clusters": 6,
                "silhouette_score": 0.7,
                "stability_score": 0.9,
                "cluster_coherence": 0.8,
                "splitting_gain": 3.0,
            },
            {
                "n_neighbors": 30,
                "resolution": 1.0,
                "n_clusters": 9,
                "silhouette_score": 0.5,
                "stability_score": 0.7,
                "cluster_coherence": 0.6,
                "splitting_gain": 2.0,
            },
        ]

        _n, _r, method, reason = _select_multi_metric(valid)

        assert method == "multi_metric"
        weight_sum = _sum_weights(reason)
        assert abs(weight_sum - 1.0) < 0.01, (
            f"POST-FIX: 4-metric weights (no kb_rate) should sum to 1.0, "
            f"got {weight_sum:.4f} (raw reason: {reason[:200]})"
        )

    # ── Test 3: 3-metric weights (no coherence) sum to 1.0 ────────────────────────

    def test_3_metric_weights(self) -> None:
        """Entries without ``cluster_coherence`` → degrade to sil+stab weights sum to 1.0.

        POST-FIX: P3-A degrade path — when no entry has ``cluster_coherence``,
        weights degrade to ``{silhouette: 0.5, stability: 0.5}``.
        """
        valid = [
            {
                "n_neighbors": 10,
                "resolution": 0.5,
                "n_clusters": 3,
                "silhouette_score": 0.3,
                "stability_score": 0.6,
            },
            {
                "n_neighbors": 20,
                "resolution": 0.8,
                "n_clusters": 6,
                "silhouette_score": 0.7,
                "stability_score": 0.9,
            },
            {
                "n_neighbors": 30,
                "resolution": 1.0,
                "n_clusters": 9,
                "silhouette_score": 0.5,
                "stability_score": 0.7,
            },
        ]

        _n, _r, method, reason = _select_multi_metric(valid)

        assert method == "multi_metric"
        weight_sum = _sum_weights(reason)
        assert abs(weight_sum - 1.0) < 0.01, (
            f"POST-FIX: 3-metric degrade weights (sil+stab only) should sum to 1.0, "
            f"got {weight_sum:.4f} (raw reason: {reason[:200]})"
        )

    # ── Test 4: coherence mismatch degrade shows in reason ────────────────────────

    def test_coherence_mismatch_degrades(self, caplog) -> None:
        """All entries have ``cluster_coherence < 0.1`` → mismatch degrade triggered.

        POST-FIX: degrade shows in reason string with correct weight sum.
        """
        valid = [
            {
                "n_neighbors": 10,
                "resolution": 0.5,
                "n_clusters": 3,
                "silhouette_score": 0.4,
                "stability_score": 0.8,
                "cluster_coherence": 0.0,
            },
            {
                "n_neighbors": 20,
                "resolution": 0.8,
                "n_clusters": 6,
                "silhouette_score": 0.7,
                "stability_score": 0.9,
                "cluster_coherence": 0.02,
            },
            {
                "n_neighbors": 30,
                "resolution": 1.0,
                "n_clusters": 9,
                "silhouette_score": 0.6,
                "stability_score": 0.7,
                "cluster_coherence": 0.05,
            },
        ]

        with caplog.at_level(logging.WARNING):
            _n, _r, method, reason = _select_multi_metric(valid)

        assert method == "multi_metric"

        # Warning emitted about mismatch degrade
        assert any(
            "mismatch" in rec.message or "Degrading" in rec.message for rec in caplog.records
        ), "Expected warning about coherence mismatch degrade"

        # POST-FIX: reason contains weight info with correct sum
        weight_sum = _sum_weights(reason)
        assert abs(weight_sum - 1.0) < 0.01, (
            f"POST-FIX: after coherence mismatch degrade, weights should sum to 1.0, "
            f"got {weight_sum:.4f} (raw reason: {reason[:200]})"
        )

    # ── Test 5: 3-tier recommendation in result (M4 post-fix) ─────────────────────

    def test_tier_recommendations_in_result(self) -> None:
        """3-tier resolution recommendation returned in the result.

        POST-FIX: M4 — current code only logs the 3-tier info; post-fix returns
        it embedded in the reason dict / string so callers can inspect it.
        """
        valid = [
            {
                "n_neighbors": 15,
                "resolution": 0.1,
                "n_clusters": 3,
                "silhouette_score": 0.05,
                "stability_score": 0.70,
            },
            {
                "n_neighbors": 15,
                "resolution": 0.3,
                "n_clusters": 5,
                "silhouette_score": 0.30,
                "stability_score": 0.80,
            },
            {
                "n_neighbors": 15,
                "resolution": 0.5,
                "n_clusters": 8,
                "silhouette_score": 0.70,
                "stability_score": 0.85,
            },
            {
                "n_neighbors": 15,
                "resolution": 0.8,
                "n_clusters": 10,
                "silhouette_score": 0.60,
                "stability_score": 0.90,
            },
            {
                "n_neighbors": 15,
                "resolution": 1.0,
                "n_clusters": 12,
                "silhouette_score": 0.20,
                "stability_score": 0.88,
            },
        ]

        result = _select_multi_metric(valid)

        # Return signature still a 4-tuple
        assert isinstance(result, tuple)
        assert len(result) == 4

        n_neighbors, resolution, method, reason = result
        assert method == "multi_metric"
        assert isinstance(n_neighbors, int)
        assert isinstance(resolution, float)
        assert isinstance(reason, str)

        # POST-FIX M4: reason contains coarse / balanced / fine tier info
        # Check for structural patterns that identify the 3 tiers
        has_coarse = "coarse" in reason.lower()
        has_balanced = "balanced" in reason.lower()
        has_fine = "fine" in reason.lower()

        assert has_coarse and has_balanced and has_fine, (
            f"POST-FIX M4: reason should contain 3-tier recommendation "
            f"(coarse, balanced, fine). Got reason: {reason[:300]}"
        )

        # Verify all three have resolution values (not nan/0)
        tier_pattern = re.compile(
            r"coarse:\s*r=([\d.]+|nan)\s*\(k=(\d+)\)\s*"
            r"/\s*balanced:\s*r=([\d.]+|nan)\s*\(k=(\d+)\)\s*"
            r"/\s*fine:\s*r=([\d.]+|nan)\s*\(k=(\d+)\)",
            re.IGNORECASE,
        )
        m = tier_pattern.search(reason)
        if m:
            # resolutions should be positive
            for r_val in [m.group(1), m.group(3), m.group(5)]:
                if r_val != "nan":
                    assert float(r_val) > 0
            # k should be positive
            for k_val in [m.group(2), m.group(4), m.group(6)]:
                assert int(k_val) > 0, f"Expected k > 0, got {k_val}"
