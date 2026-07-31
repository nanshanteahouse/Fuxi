"""Tests for rna/utils/tiered_annotation.py - hierarchical label resolution."""

from core.annotation.scoring import Score
from rna.utils.tiered_annotation import format_subtype_candidates, resolve_tiered_label


def _s(score: float) -> Score:
    """Build a minimal Score for testing."""
    return Score(
        score=score,
        p_value=0.05,
        method="test",
        n_markers_found=1,
        negative_penalty=False,
    )


HIERARCHY = {
    "categories": {
        "Neuron": {
            "members": ["RGC", "Amacrine_Cell", "RGC_Foxp2", "RGC_Alpha"],
            "subtypes": {
                "RGC": {"members": ["RGC_Foxp2", "RGC_Alpha"]},
            },
        }
    },
    "incompatible_transitions": [],
}


def _kb_lookup(foxp2_consensus: str = "medium", foxp2_private: list | None = None) -> dict:
    """Build a KB lookup with RGC (no subtypes in KB usage here) and RGC_Foxp2."""
    if foxp2_private is None:
        foxp2_private = ["FOXP2"]
    return {
        "Broad_Neuron": {
            "parent": "",
            "_private_markers": [],
            "consensus_levels": {},
            "marker_weights": {},
        },
        "RGC": {
            "parent": "Broad_Neuron",
            "_private_markers": ["RBPMS", "POU4F1"],
            "consensus_levels": {"RBPMS": "high", "POU4F1": "high"},
            "marker_weights": {"RBPMS": 10},
        },
        "Amacrine_Cell": {
            "parent": "Broad_Neuron",
            "_private_markers": ["GAD1"],
            "consensus_levels": {"GAD1": "high"},
            "marker_weights": {"GAD1": 5},
        },
        "RGC_Foxp2": {
            "parent": "RGC",
            "_private_markers": foxp2_private,
            "consensus_levels": {"FOXP2": foxp2_consensus},
            "marker_weights": {"FOXP2": 2},
        },
        "RGC_Alpha": {
            "parent": "RGC",
            "_private_markers": [],
            "consensus_levels": {},
            "marker_weights": {},
        },
    }


class TestResolveTieredLabel:
    def test_l2_no_subtypes_na(self) -> None:
        """L2 with no subtypes in KB → tier=L2, resolution=na."""
        scores = {
            "Broad_Neuron": _s(0.7),
            "Amacrine_Cell": _s(0.6),
            "RGC": _s(0.5),
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["GAD1"])
        assert label == "Amacrine_Cell"
        assert ev["tier"] == "L2"
        assert ev["subtype_resolution"] == "na"

    def test_l3_resolved_all_gates_pass(self) -> None:
        """Subtype passes all 3 gates → tier=L3, resolution=resolved."""
        scores = {
            "Broad_Neuron": _s(0.7),
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.58),  # within delta of RGC (0.60-0.08=0.52)
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2", "RBPMS"])
        assert label == "RGC_Foxp2"
        assert ev["tier"] == "L3"
        assert ev["subtype_resolution"] == "resolved"
        assert "RGC_Foxp2" in ev["available_subtypes"]

    def test_unresolved_score_gate_fail(self) -> None:
        """Subtype score too far below L2 → gate A fails → unresolved."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.40),  # 0.40 < 0.52 → gate A fails
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2"])
        assert label == "RGC"
        assert ev["tier"] == "L2"
        assert ev["subtype_resolution"] == "unresolved"

    def test_unresolved_private_marker_gate_fail(self) -> None:
        """Score OK but no private marker hit → gate B fails → unresolved."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.58),
        }
        # FOXP2 NOT in cluster genes → gate B fails
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["RBPMS"])
        assert label == "RGC"
        assert ev["subtype_resolution"] == "unresolved"

    def test_unresolved_consensus_gate_fail(self) -> None:
        """Private marker hit but consensus=low < floor=medium → gate C fails."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.58),
        }
        kb = _kb_lookup(foxp2_consensus="low")  # FOXP2 only low consensus
        label, ev = resolve_tiered_label(scores, HIERARCHY, kb, ["FOXP2"])
        assert label == "RGC"
        assert ev["subtype_resolution"] == "unresolved"

    def test_l3_below_min_score(self) -> None:
        """Subtype within delta but below min_score → gate A fails."""
        scores = {
            "RGC": _s(0.30),
            "RGC_Foxp2": _s(0.24),  # < 0.25 min_score
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2"])
        assert label == "RGC"
        assert ev["subtype_resolution"] == "unresolved"

    def test_no_l2_candidate_returns_empty(self) -> None:
        """Only Broad_* and subtypes scored, no L2 → empty label."""
        scores = {
            "Broad_Neuron": _s(0.7),
            "RGC_Foxp2": _s(0.58),
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2"])
        assert label == ""
        assert ev["subtype_resolution"] == "na"

    def test_l2_unresolved_when_subtype_not_scored(self) -> None:
        """L2 with subtypes defined but no subtype scored → unresolved."""
        scores = {
            "RGC": _s(0.60),
            # RGC_Foxp2 / RGC_Alpha not in scores
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["RBPMS"])
        assert label == "RGC"
        assert ev["tier"] == "L2"
        assert ev["subtype_resolution"] == "unresolved"

    def test_evidence_carries_consensus_and_n_sources(self) -> None:
        """Evidence dict should carry consensus + n_sources for the chosen type."""
        scores = {"RGC": _s(0.6), "RGC_Foxp2": _s(0.58)}
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2", "RBPMS"])
        assert label == "RGC_Foxp2"
        # RGC_Foxp2 hit FOXP2 (medium consensus)
        assert ev["consensus"] == "medium"
        assert ev["n_sources"] == 2  # FOXP2 marker_weight=2

    def test_subtype_with_empty_private_markers_fails_gate_b(self) -> None:
        """Subtype with _private_markers=[] always fails gate B."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Alpha": _s(0.59),  # within delta
        }
        kb = _kb_lookup()
        label, ev = resolve_tiered_label(scores, HIERARCHY, kb, ["RBPMS", "FOXP2"])
        assert label == "RGC"
        assert ev["subtype_resolution"] == "unresolved"

    def test_custom_thresholds_override(self) -> None:
        """Caller can override gates (e.g. relaxed delta)."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.48),  # fails default delta(0.08) but passes 0.15
        }
        label, ev = resolve_tiered_label(
            scores, HIERARCHY, _kb_lookup(), ["FOXP2"], subtype_delta=0.15
        )
        assert label == "RGC_Foxp2"
        assert ev["subtype_resolution"] == "resolved"

    def test_subtype_candidates_present_on_resolved_winner(self) -> None:
        """Winner (looked up by type, not [0]) has failed_gates=[]; sorted desc."""
        scores = {
            "Broad_Neuron": _s(0.7),
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.58),  # passes all gates
            "RGC_Alpha": _s(0.62),  # higher score but _private_markers=[] → gate B fails
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2", "RBPMS"])
        assert label == "RGC_Foxp2"
        cands = ev["subtype_candidates"]
        # All scored L3s of the best L2 present, sorted by score descending.
        assert [c["type"] for c in cands] == ["RGC_Alpha", "RGC_Foxp2"]
        assert [c["score"] for c in cands] == [0.62, 0.58]
        # Winner entry: failed_gates == [] and private-marker hits recorded.
        winner = next(c for c in cands if c["type"] == "RGC_Foxp2")
        assert winner["failed_gates"] == []
        assert winner["private_markers_hit"] == 1
        # Higher-scoring gate-failing subtype occupies index 0 and carries its failure.
        alpha = cands[0]
        assert alpha["type"] == "RGC_Alpha"
        assert "B" in alpha["failed_gates"]
        assert alpha["private_markers_hit"] == 0

    def test_subtype_candidates_gate_a_fail(self) -> None:
        """Score too far below L2 → gate A recorded as failed."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.40),  # 0.40 < 0.52 → gate A fails
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2"])
        assert label == "RGC"
        cands = ev["subtype_candidates"]
        assert len(cands) == 1
        assert cands[0]["type"] == "RGC_Foxp2"
        assert cands[0]["failed_gates"] == ["A"]

    def test_subtype_candidates_gate_b_fail(self) -> None:
        """Private marker absent from cluster genes → gate B recorded as failed."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.58),
        }
        # FOXP2 NOT in cluster genes → gate B fails (and gate C too, no hit markers).
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["RBPMS"])
        assert label == "RGC"
        cands = ev["subtype_candidates"]
        assert cands[0]["type"] == "RGC_Foxp2"
        assert "B" in cands[0]["failed_gates"]

    def test_subtype_candidates_gate_c_fail(self) -> None:
        """Consensus below floor → gate C recorded as failed."""
        scores = {
            "RGC": _s(0.60),
            "RGC_Foxp2": _s(0.58),
        }
        kb = _kb_lookup(foxp2_consensus="low")  # FOXP2 low < medium floor
        label, ev = resolve_tiered_label(scores, HIERARCHY, kb, ["FOXP2"])
        assert label == "RGC"
        cands = ev["subtype_candidates"]
        assert cands[0]["type"] == "RGC_Foxp2"
        assert cands[0]["failed_gates"] == ["C"]

    def test_subtype_candidates_empty_when_l2_has_no_subtypes(self) -> None:
        """L2 without KB subtypes → subtype_candidates == [] (na path)."""
        scores = {
            "Broad_Neuron": _s(0.7),
            "Amacrine_Cell": _s(0.6),
            "RGC": _s(0.5),
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["GAD1"])
        assert label == "Amacrine_Cell"
        assert ev["subtype_candidates"] == []

    def test_subtype_candidates_empty_when_no_l2_candidate(self) -> None:
        """No L2 candidate → inline dict path carries empty subtype_candidates."""
        scores = {
            "Broad_Neuron": _s(0.7),
            "RGC_Foxp2": _s(0.58),
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["FOXP2"])
        assert label == ""
        assert ev["subtype_candidates"] == []

    def test_subtype_candidates_empty_when_none_scored(self) -> None:
        """L2 subtypes exist in KB but none scored → candidates [] while available non-empty."""
        scores = {
            "RGC": _s(0.60),
            # RGC_Foxp2 / RGC_Alpha not scored
        }
        label, ev = resolve_tiered_label(scores, HIERARCHY, _kb_lookup(), ["RBPMS"])
        assert label == "RGC"
        assert ev["subtype_resolution"] == "unresolved"
        assert ev["subtype_candidates"] == []
        assert len(ev["available_subtypes"]) == 2  # RGC_Foxp2, RGC_Alpha


class TestFormatSubtypeCandidates:
    def test_joins_entries_with_gates(self) -> None:
        cands = [
            {"type": "RGC_Foxp2", "score": 0.58, "failed_gates": ["B", "C"]},
            {"type": "RGC_Alpha", "score": 0.55, "failed_gates": ["C"]},
        ]
        assert format_subtype_candidates(cands) == "RGC_Foxp2:0.58[B,C]; RGC_Alpha:0.55[C]"

    def test_winner_omits_brackets(self) -> None:
        cands = [{"type": "RGC_Foxp2", "score": 0.58, "failed_gates": []}]
        assert format_subtype_candidates(cands) == "RGC_Foxp2:0.58"

    def test_empty_input(self) -> None:
        assert format_subtype_candidates([]) == ""

    def test_score_rounded_to_2dp(self) -> None:
        cands = [{"type": "RGC_Foxp2", "score": 0.577, "failed_gates": []}]
        assert format_subtype_candidates(cands) == "RGC_Foxp2:0.58"
