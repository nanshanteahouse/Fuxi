"""Numerical tests for rna/utils/marker_scoring.py."""

import pandas as pd
import pytest

from rna.utils.marker_scoring import (
    score_cluster_against_kb,
    Score,
    resolve_expert_rule_params,
    apply_expert_rules,
    detect_low_quality_cluster,
    _negative_marker_penalty,
)


class TestMarkerScoringImport:
    """Verify that marker-scoring symbols are importable."""

    def test_import_marker_scoring(self) -> None:
        assert score_cluster_against_kb is not None
        assert Score is not None
        assert resolve_expert_rule_params is not None
        assert apply_expert_rules is not None
        assert detect_low_quality_cluster is not None
        assert _negative_marker_penalty is not None

    def test_score_namedtuple_fields(self) -> None:
        """Score namedtuple should have the expected fields."""
        fields = list(Score._fields)
        assert "score" in fields
        assert "p_value" in fields
        assert "method" in fields
        assert "n_markers_found" in fields
        assert "negative_penalty" in fields


class TestScoreClusterAgainstKB:
    """Numerical assertions for score_cluster_against_kb.

    Uses a pre-built KB lookup dict to avoid the _build_kb_lookup
    code path — tests focus on the Fisher exact scoring, cosine
    similarity, and confidence multiplier logic.
    """

    # Two cell types with disjoint positive markers.
    KB_LOOKUP = {
        "T_cell": {
            "positive": ["T1", "T2", "T3", "T4", "T5"],
            "negative": [],
            "species": ["human"],
            "synonyms": [],
            "parent": "Lymphocyte",
            "marker_weights": {},
            "consensus_levels": {},
        },
        "B_cell": {
            "positive": ["B1", "B2", "B3"],
            "negative": [],
            "species": ["human"],
            "synonyms": [],
            "parent": "Lymphocyte",
            "marker_weights": {},
            "consensus_levels": {},
        },
    }

    def test_fisher_exact_known_enrichment(self) -> None:
        """Fisher exact test with known overlap -> expected score > 0.7.

        T_cell has 5/5 markers in top-20 (a=5, b=0, c=0, d=3).
        P = C(5,5)*C(3,0)/C(8,5) = 1/56 approx 0.01786; 1 - p = 0.9821.
        conf_mult = 0.8 (5 type markers, not > 5).
        final approx 0.9821 * 0.8 = 0.7857.
        """
        cluster_markers = pd.DataFrame({
            "names": ["T1", "T2", "T3", "T4", "T5"]
                     + [f"X{i}" for i in range(15)],
            "logfoldchanges": [3.0 - i * 0.1 for i in range(20)],
            "pvals_adj": [1e-10] * 20,
        })
        result = score_cluster_against_kb(self.KB_LOOKUP, cluster_markers)

        t = result["T_cell"]
        assert t.score == pytest.approx(0.7857, rel=1e-3), f"T_cell score {t.score}"
        assert t.p_value == pytest.approx(1 / 56, rel=1e-3)
        assert t.method == "hypergeometric"
        assert t.n_markers_found == 5
        assert not t.negative_penalty

        # B_cell has zero overlap -> score 0, p=1
        b = result["B_cell"]
        assert b.score == pytest.approx(0.0, abs=1e-6)
        assert b.p_value == pytest.approx(1.0, abs=1e-6)
        assert b.n_markers_found == 0

    def test_no_kb_overlap_returns_zero(self) -> None:
        """Cluster markers not overlapping any KB type -> all scores 0."""
        cluster_markers = pd.DataFrame({
            "names": [f"GENE{i}" for i in range(20)],
            "logfoldchanges": [1.0] * 20,
            "pvals_adj": [0.05] * 20,
        })
        result = score_cluster_against_kb(self.KB_LOOKUP, cluster_markers)
        for cell_type, sc in result.items():
            assert sc.score == pytest.approx(0.0, abs=1e-6), (
                f"{cell_type} score={sc.score}"
            )
            assert sc.p_value == pytest.approx(1.0, abs=1e-6)

    def test_no_positive_markers_score_zero(self) -> None:
        """Cell type with empty positive marker list -> Score(0, 1, 'none', 0, False)."""
        kb = {
            "Null_type": {
                "positive": [],
                "negative": [],
                "species": ["human"],
                "synonyms": [],
                "parent": "",
                "marker_weights": {},
                "consensus_levels": {},
            },
        }
        cluster_markers = pd.DataFrame({
            "names": ["T1", "T2", "T3"] + [f"X{i}" for i in range(17)],
            "logfoldchanges": [2.0] * 20,
            "pvals_adj": [1e-5] * 20,
        })
        result = score_cluster_against_kb(kb, cluster_markers)
        sc = result["Null_type"]
        assert sc.score == pytest.approx(0.0)
        assert sc.p_value == pytest.approx(1.0)
        assert sc.method == "none"
        assert sc.n_markers_found == 0
        assert not sc.negative_penalty
