"""Numerical tests for rna/utils/evidence_fusion.py."""

import pytest

from rna.utils.evidence_fusion import (
    fuse_all_clusters,
    FusionDecision,
    DiagnosticInfo,
)


class TestEvidenceFusionImport:
    """Verify that evidence-fusion symbols are importable."""

    def test_import_evidence_fusion(self) -> None:
        assert fuse_all_clusters is not None
        assert FusionDecision is not None
        assert DiagnosticInfo is not None

    def test_fusion_decision_namedtuple_fields(self) -> None:
        """FusionDecision namedtuple should have the expected fields."""
        fields = set(FusionDecision._fields)
        expected = {
            "cell_type",
            "confidence",
            "score",
            "method",
            "n_markers_found",
            "ai_agreed",
            "ai_suggested",
            "explanation",
            "alternative_rules",
            "diagnostic",
        }
        assert fields == expected, f"Missing fields: {expected - fields}"


class TestFuseAllClusters:
    """Numerical assertions for fuse_all_clusters.

    Uses bare floats as marker scores (fuse_evidence._resolve_score
    handles both Score objects and floats).
    """

    def test_high_confidence_from_marker_scores(self) -> None:
        """Best score >= 0.7 -> 'high' confidence."""
        all_scores = {
            "0": {"T_cell": 0.85, "B_cell": 0.30},
        }
        all_rules = {"0": None}
        decisions = fuse_all_clusters(all_scores, all_rules)

        assert len(decisions) == 1
        d = decisions[0]
        assert d.cell_type == "T_cell"
        assert d.confidence == "high"
        assert d.score == pytest.approx(0.85)
        assert d.method == "marker_scoring_high"

    def test_medium_and_low_confidence(self) -> None:
        """Scores in [0.5, 0.7) -> 'medium'; [0.25, 0.5) -> 'low'."""
        all_scores = {
            "0": {"T_cell": 0.60, "B_cell": 0.40},
            "1": {"B_cell": 0.35},
        }
        all_rules = {"0": None, "1": None}
        decisions = fuse_all_clusters(all_scores, all_rules)

        assert decisions[0].confidence == "medium"
        assert decisions[0].cell_type == "T_cell"
        assert decisions[0].score == pytest.approx(0.60)

        assert decisions[1].confidence == "low"
        assert decisions[1].cell_type == "B_cell"
        assert decisions[1].score == pytest.approx(0.35)

    def test_empty_scores_returns_unknown(self) -> None:
        """Empty marker_score dict -> confidence 'unknown', cell_type 'Unknown'."""
        all_scores = {"0": {}}
        all_rules = {"0": None}
        decisions = fuse_all_clusters(all_scores, all_rules)

        assert len(decisions) == 1
        d = decisions[0]
        assert d.cell_type == "Unknown"
        assert d.confidence == "unknown"
        assert d.score == pytest.approx(0.0)
        assert d.method == "unknown"

    def test_return_quality_metadata(self) -> None:
        """When return_quality=True, second element is a quality dict."""
        all_scores = {
            "0": {"T_cell": 0.85},
            "1": {"B_cell": 0.60},
        }
        all_rules = {"0": None, "1": None}
        _, quality = fuse_all_clusters(
            all_scores, all_rules, return_quality=True,
        )
        assert quality["total"] == 2
        assert quality["unknown"] == 0
        assert quality["annotated_by_scoring"] == 2
        assert quality["annotated_by_rule"] == 0
