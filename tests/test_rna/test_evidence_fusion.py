"""Tests for rna/utils/evidence_fusion.py — private helper
_is_transition_state and public API fuse_all_clusters/FusionDecision/DiagnosticInfo."""

import pytest

from rna.utils.evidence_fusion import (
    DiagnosticInfo,
    FusionDecision,
    _is_transition_state,
    fuse_all_clusters,
)

# ── Public API ────────────────────────────────────────────────────────


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
            "cell_category",
            "confidence",
            "score",
            "method",
            "n_markers_found",
            "ai_agreed",
            "ai_suggested",
            "explanation",
            "alternative_rules",
            "diagnostic",
            "tier",
            "consensus",
            "n_sources",
            "subtype_resolution",
        }
        missing = expected - fields
        extra = fields - expected
        assert fields == expected, f"missing={missing}, extra={extra}"


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
            all_scores,
            all_rules,
            return_quality=True,
        )
        assert quality["total"] == 2
        assert quality["unknown"] == 0
        assert quality["annotated_by_scoring"] == 2
        assert quality["annotated_by_rule"] == 0

    # ── Multi-peak downgrade (P0.6) ───────────────────────────────────

    def test_multi_peak_three_tied_types_ambiguous(self) -> None:
        """>= 3 types with score >= 0.9 on a confirmed transition candidate
        (same parent + shared marker) → Unknown/ambiguous with full
        top_competitors (every tied type, not just top-3).
        """
        all_scores = {
            "0": {"RGC": 1.0, "Amacrine": 1.0, "Bipolar": 1.0, "HC": 0.8},
        }
        all_rules = {"0": None}
        # Top-2 (RGC/Amacrine) share Broad_Neuron + SNCG → confirmed
        # transition candidate, so the multi-peak check applies.
        kb = {
            "RGC": {
                "parent": "Broad_Neuron",
                "markers": {"confirm": {"RBPMS": ["src"], "SNCG": ["src"]}, "add": {}},
            },
            "Amacrine": {
                "parent": "Broad_Neuron",
                "markers": {"confirm": {"TFAP2A": ["src"], "SNCG": ["src"]}, "add": {}},
            },
            "Bipolar": {
                "parent": "Broad_Neuron",
                "markers": {"confirm": {"VSX2": ["src"], "SNCG": ["src"]}, "add": {}},
            },
        }
        decisions = fuse_all_clusters(
            all_scores,
            all_rules,
            kb=kb,
            allows_transitions=True,
        )

        assert len(decisions) == 1
        d = decisions[0]
        assert d.cell_type == "Unknown"
        assert d.confidence == "unknown"
        assert d.method == "ambiguous"
        assert d.score == pytest.approx(1.0)
        assert d.diagnostic is not None
        assert d.diagnostic.category == "ambiguous"
        competitors = {c["cell_type"] for c in d.diagnostic.top_competitors}
        assert competitors == {"RGC", "Amacrine", "Bipolar"}
        assert all(c["score"] >= 0.9 for c in d.diagnostic.top_competitors)

    def test_multi_peak_two_tied_types_still_transition(self) -> None:
        """Only 2 types tied (below multi_peak_min_types=3) → the transition
        path still runs. Exercises fuse_all_clusters → fuse_evidence
        pass-through of the multi-peak parameters.
        """
        all_scores = {
            "0": {"RGC": 0.95, "Amacrine": 0.93},
        }
        all_rules = {"0": None}
        kb = {
            "RGC": {
                "parent": "Broad_Neuron",
                "markers": {
                    "confirm": {"RBPMS": ["src"], "SNCG": ["src"]},
                    "add": {},
                },
            },
            "Amacrine": {
                "parent": "Broad_Neuron",
                "markers": {
                    "confirm": {"TFAP2A": ["src"], "SNCG": ["src"]},
                    "add": {},
                },
            },
        }
        decisions = fuse_all_clusters(
            all_scores,
            all_rules,
            kb=kb,
            allows_transitions=True,
        )

        d = decisions[0]
        assert d.method == "transition_state"
        assert d.cell_type == "transitional: RGC/Amacrine"

    def test_multi_peak_high_min_types_no_downgrade(self) -> None:
        """multi_peak_min_types=6 with only 3 tied types → no downgrade;
        normal logic applies (transition_state here).
        """
        all_scores = {
            "0": {"RGC": 1.0, "Amacrine": 1.0, "Bipolar": 1.0, "HC": 0.8},
        }
        all_rules = {"0": None}
        kb = {
            "RGC": {
                "parent": "Broad_Neuron",
                "markers": {"confirm": {"RBPMS": ["src"], "SNCG": ["src"]}, "add": {}},
            },
            "Amacrine": {
                "parent": "Broad_Neuron",
                "markers": {"confirm": {"TFAP2A": ["src"], "SNCG": ["src"]}, "add": {}},
            },
            "Bipolar": {
                "parent": "Broad_Neuron",
                "markers": {"confirm": {"VSX2": ["src"], "SNCG": ["src"]}, "add": {}},
            },
        }
        decisions = fuse_all_clusters(
            all_scores,
            all_rules,
            kb=kb,
            allows_transitions=True,
            multi_peak_min_types=6,
        )

        d = decisions[0]
        assert d.method == "transition_state"
        assert d.cell_type == "transitional: RGC/Amacrine"

    def test_multi_peak_non_candidate_not_degraded(self) -> None:
        """Multi-peak saturation on a NON-transition candidate (top-2 have
        different parents) must NOT downgrade — the cluster keeps its
        marker-scoring label. Pins the scope-narrowing semantics.
        """
        all_scores = {
            "0": {
                "Muller_Glia": 1.0,
                "Proliferating_RPC": 1.0,
                "NRPC": 1.0,
                "RPC": 0.95,
            },
        }
        all_rules = {"0": None}
        kb = {
            "Muller_Glia": {"parent": "Broad_Glia"},
            "Proliferating_RPC": {"parent": "Broad_Progenitor"},
            "NRPC": {"parent": "Broad_Progenitor"},
            "RPC": {"parent": "Broad_Progenitor"},
        }
        decisions = fuse_all_clusters(
            all_scores,
            all_rules,
            kb=kb,
            allows_transitions=True,
        )

        d = decisions[0]
        assert d.method != "ambiguous"
        assert d.method == "marker_scoring_high"
        assert d.cell_type == "Muller_Glia"
        assert d.confidence == "high"


# ── Private helper: _is_transition_state ──────────────────────────────


class TestIsTransitionState:
    """Edge-case coverage for the transition-state detection function."""

    # ── Happy path ───────────────────────────────────────────────────

    def test_transition_normal(self) -> None:
        """delta < threshold, same parent, top score >= 0.25 → (top1, top2).

        Both types carry shared type-specific markers (SNCG) so the P0.4
        private-marker gate passes.
        """
        scores = {"RGC": 0.45, "Amacrine_Cell": 0.35}
        kb = {
            "RGC": {
                "parent": "Broad_Neuron",
                "markers": {
                    "confirm": {"RBPMS": ["src"], "SNCG": ["src"]},
                    "add": {},
                },
            },
            "Amacrine_Cell": {
                "parent": "Broad_Neuron",
                "markers": {
                    "confirm": {"TFAP2A": ["src"], "SNCG": ["src"]},
                    "add": {},
                },
            },
        }
        result = _is_transition_state(scores, kb)
        assert result is not None
        # Sorted descending: RGC (0.45) first, Amacrine_Cell (0.35) second
        assert result[0] == "RGC"
        assert result[1] == "Amacrine_Cell"

    # ── P0.4 private-marker gate ──────────────────────────────────────

    def test_transition_markers_wrapped_by_parent_no_private(self) -> None:
        """P0.4: both types' markers fully covered by parent confirm and no
        shared ``_private_markers`` → None (the gate must not be skipped).
        """
        scores = {"NRPC_Cone_BC_fate": 0.45, "NRPC_RGC_fate": 0.35}
        kb = {
            "Broad_Progenitor": {
                "markers": {
                    "confirm": {"VSX2": ["src"], "SOX2": ["src"], "HES1": ["src"]},
                    "add": {},
                },
            },
            "NRPC_Cone_BC_fate": {
                "parent": "Broad_Progenitor",
                "markers": {
                    "confirm": {"VSX2": ["src"], "SOX2": ["src"]},
                    "add": {},
                },
            },
            "NRPC_RGC_fate": {
                "parent": "Broad_Progenitor",
                "markers": {
                    "confirm": {"HES1": ["src"]},
                    "add": {},
                },
            },
        }
        result = _is_transition_state(scores, kb)
        assert result is None

    def test_transition_markers_wrapped_by_parent_shared_private(self) -> None:
        """P0.4: markers covered by parent but both types share a KB
        ``_private_markers`` entry → transition accepted (tuple returned).
        """
        scores = {"NRPC_Cone_BC_fate": 0.45, "NRPC_RGC_fate": 0.35}
        kb = {
            "Broad_Progenitor": {
                "markers": {
                    "confirm": {"VSX2": ["src"], "SOX2": ["src"], "HES1": ["src"]},
                    "add": {},
                },
            },
            "NRPC_Cone_BC_fate": {
                "parent": "Broad_Progenitor",
                "markers": {
                    "confirm": {"VSX2": ["src"], "SOX2": ["src"]},
                    "add": {},
                },
                "_private_markers": ["PRDM1"],
            },
            "NRPC_RGC_fate": {
                "parent": "Broad_Progenitor",
                "markers": {
                    "confirm": {"HES1": ["src"]},
                    "add": {},
                },
                "_private_markers": ["PRDM1"],
            },
        }
        result = _is_transition_state(scores, kb)
        assert result is not None
        assert result[0] == "NRPC_Cone_BC_fate"
        assert result[1] == "NRPC_RGC_fate"

    # ── Delta threshold ──────────────────────────────────────────────

    def test_transition_delta_too_large(self) -> None:
        """delta >= 0.15 → None (gap too large to be a transition)."""
        scores = {"RGC": 0.50, "Amacrine_Cell": 0.30}
        kb = {
            "RGC": {"parent": "Broad_Neuron"},
            "Amacrine_Cell": {"parent": "Broad_Neuron"},
        }
        result = _is_transition_state(scores, kb)
        assert result is None

    # ── Parent mismatch ──────────────────────────────────────────────

    def test_transition_mismatched_parents(self) -> None:
        """Same score delta but different parents → None."""
        scores = {"RGC": 0.45, "Amacrine_Cell": 0.35}
        kb = {
            "RGC": {"parent": "Broad_Neuron"},
            "Amacrine_Cell": {"parent": "Glia"},  # different parent
        }
        result = _is_transition_state(scores, kb)
        assert result is None

    # ── Score floor ──────────────────────────────────────────────────

    def test_transition_below_score_floor(self) -> None:
        """Top score < 0.25 → None (not confident enough to label)."""
        scores = {"RGC": 0.20, "Amacrine_Cell": 0.15}
        kb = {
            "RGC": {"parent": "Broad_Neuron"},
            "Amacrine_Cell": {"parent": "Broad_Neuron"},
        }
        result = _is_transition_state(scores, kb)
        assert result is None

    # ── Insufficient entries ─────────────────────────────────────────

    def test_transition_single_entry(self) -> None:
        """Only 1 type in marker_scores → None (< 2 entries)."""
        scores = {"RGC": 0.45}
        result = _is_transition_state(scores, {})
        assert result is None

    def test_transition_empty_dict(self) -> None:
        """Empty marker_scores dict → None."""
        result = _is_transition_state({}, {})
        assert result is None

    # ── KB edge cases ────────────────────────────────────────────────

    def test_transition_kb_none(self) -> None:
        """kb is None → AttributeError (function does not guard against None).

        ``_is_transition_state`` calls ``kb.get(...)`` directly without a
        None check, so passing ``kb=None`` raises ``AttributeError``.
        This test documents that contract — callers must ensure kb is a dict.
        """
        scores = {"RGC": 0.45, "Amacrine_Cell": 0.35}
        with pytest.raises(AttributeError):
            _is_transition_state(scores, None)

    def test_transition_missing_parent_field(self) -> None:
        """One type has no 'parent' in KB → parent is '' → None."""
        scores = {"RGC": 0.45, "Amacrine_Cell": 0.35}
        kb = {
            "RGC": {"parent": "Broad_Neuron"},
            "Amacrine_Cell": {},  # no parent field
        }
        result = _is_transition_state(scores, kb)
        # parent2 defaults to '' → parent1 != parent2 → None
        assert result is None
