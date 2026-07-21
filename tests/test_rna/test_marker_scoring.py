"""Numerical tests for rna/utils/marker_scoring.py."""

import pandas as pd
import pytest

from core.annotation.scoring import (
    Score,
    _negative_marker_penalty,
    detect_low_quality_cluster,
    score_cluster_against_kb,
)
from rna.utils.marker_expert_rules import (
    apply_expert_rules,
    resolve_expert_rule_params,
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
        cluster_markers = pd.DataFrame(
            {
                "names": ["T1", "T2", "T3", "T4", "T5"] + [f"X{i}" for i in range(15)],
                "logfoldchanges": [3.0 - i * 0.1 for i in range(20)],
                "pvals_adj": [1e-10] * 20,
            }
        )
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
        cluster_markers = pd.DataFrame(
            {
                "names": [f"GENE{i}" for i in range(20)],
                "logfoldchanges": [1.0] * 20,
                "pvals_adj": [0.05] * 20,
            }
        )
        result = score_cluster_against_kb(self.KB_LOOKUP, cluster_markers)
        for cell_type, sc in result.items():
            assert sc.score == pytest.approx(0.0, abs=1e-6), f"{cell_type} score={sc.score}"
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
        cluster_markers = pd.DataFrame(
            {
                "names": ["T1", "T2", "T3"] + [f"X{i}" for i in range(17)],
                "logfoldchanges": [2.0] * 20,
                "pvals_adj": [1e-5] * 20,
            }
        )
        result = score_cluster_against_kb(kb, cluster_markers)
        sc = result["Null_type"]
        assert sc.score == pytest.approx(0.0)
        assert sc.p_value == pytest.approx(1.0)
        assert sc.method == "none"
        assert sc.n_markers_found == 0
        assert not sc.negative_penalty


class TestPhylogeneticWeighting:
    """Phylogenetic weight integration via score_cluster_against_kb.

    Verifies that ``target_class`` / ``target_order`` parameters produce
    the correct multiplicative weight on the final score.

    KB entries include taxonomic fields (class/order/classes) so the
    phylogenetic weighting code path is exercised.
    """

    # Each cell type has **unique** markers (no cross-overlap),
    # which makes the Fisher p-values easy to reason about.
    # Extra markers outside the top-20 ensure a non-degenerate
    # background (d > 0).
    KB_WITH_TAXONOMY = {
        "Mammal_Cell": {
            "positive": ["M1", "M2", "M3", "M4", "M5"] + [f"Mm_extra_{i}" for i in range(10)],
            "negative": [],
            "species": ["human"],
            "synonyms": [],
            "parent": "Lymphocyte",
            "marker_weights": {},
            "consensus_levels": {},
            "class": "Mammalia",
            "order": "Primates",
            "classes": ["Mammalia"],
            "orders": ["Primates"],
        },
        "Bird_Cell": {
            "positive": ["B1", "B2", "B3"] + [f"Bb_extra_{i}" for i in range(10)],
            "negative": [],
            "species": ["chicken"],
            "synonyms": [],
            "parent": "Lymphocyte",
            "marker_weights": {},
            "consensus_levels": {},
            "class": "Aves",
            "order": "Galliformes",
            "classes": ["Aves"],
            "orders": ["Galliformes"],
        },
    }

    MARKER_DF = pd.DataFrame(
        {
            "names": ["M1", "M2", "M3", "M4", "M5", "B1", "B2", "B3"]
            + [f"X{i}" for i in range(12)],
            "logfoldchanges": [3.0 - i * 0.1 for i in range(20)],
            "pvals_adj": [1e-10] * 20,
        }
    )

    def _unweighted_score(self, scores: dict) -> float:
        """Return the raw score from a result dict without phylogenetic weight."""
        raw = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
        )
        return raw[sorted(raw.keys())[0]].score

    def test_same_class_weight(self) -> None:
        """Same class (Mammalia -> Mammalia) → phylogenetic weight = 1.0.

        The weighted score should equal the unweighted score.
        """
        unweighted = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
        )["Mammal_Cell"].score
        weighted = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
            target_class="Mammalia",
        )["Mammal_Cell"].score
        assert weighted == pytest.approx(unweighted, rel=1e-6), (
            f"Same-class weight {weighted} != unweighted {unweighted}"
        )

    def test_same_class_same_order_weight(self) -> None:
        """Same class + same order (Primates -> Primates) → weight = 1.0.

        Matching Mammalia + Primates gets the full unweighted score.
        """
        unweighted = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
        )["Mammal_Cell"].score
        weighted = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
            target_class="Mammalia",
            target_order="Primates",
        )["Mammal_Cell"].score
        assert weighted == pytest.approx(unweighted, rel=1e-6), (
            f"Same-class+order weight {weighted} != unweighted {unweighted}"
        )

    def test_cross_class_weight_bird_to_mammal(self) -> None:
        """Cross-class (Aves -> Mammalia) → weight < 1.0.

        Bird_Cell (class Aves) scored against Mammalia target should
        receive a phylogenetic penalty (weight = 0.6 for single-class
        source in a different class).
        """
        unweighted = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
        )["Bird_Cell"].score
        weighted = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
            target_class="Mammalia",
        )["Bird_Cell"].score
        assert weighted < unweighted, (
            f"Cross-class weight {weighted} should be < unweighted {unweighted}"
        )
        # With single-class Aves -> Mammalia, weight = 0.6,
        # so weighted should be exactly 0.6x unweighted.
        assert weighted == pytest.approx(unweighted * 0.6, rel=1e-6), (
            f"Cross-class weight {weighted} != {unweighted} * 0.6 = {unweighted * 0.6}"
        )

    def test_cross_class_weight_decreases_with_distance(self) -> None:
        """Cross-class weight < same-class weight for the same input data.

        Bird_Cell (Aves) against Mammalia target gets a lower score
        than Mammal_Cell (Mammalia) against the same Mammalia target.
        """
        result = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
            target_class="Mammalia",
        )
        mammal_score = result["Mammal_Cell"].score
        bird_score = result["Bird_Cell"].score
        assert bird_score < mammal_score, (
            f"Bird_Cell ({bird_score:.4f}) should score lower than "
            f"Mammal_Cell ({mammal_score:.4f}) under Mammalia target"
        )

    def test_no_target_class_no_phylogenetic_effect(self) -> None:
        """With no target_class, phylogenetic weighting is skipped entirely.

        Even with taxonomic fields in the KB, not setting target_class
        means all scores are unweighted.
        """
        result = score_cluster_against_kb(
            self.KB_WITH_TAXONOMY,
            self.MARKER_DF,
        )
        bird_score = result["Bird_Cell"].score
        # Bird_Cell markers are in the top-20, so score should be > 0
        assert bird_score > 0, f"Bird_Cell score should be > 0 without filtering, got {bird_score}"
