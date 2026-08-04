"""Tests for rna/utils/evidence_fusion.py — private helper
_is_transition_state and public API fuse_all_clusters/FusionDecision/DiagnosticInfo."""

import pytest

from core.annotation.potency import KADPConfig
from core.annotation.scoring import Score
from rna.utils.evidence_fusion import (
    DiagnosticInfo,
    FusionDecision,
    _is_transition_state,
    fuse_all_clusters,
    fuse_evidence,
    harmonize_label,
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
            "review_reason",
            "potency",
            "source_votes",
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
        d = decisions[0]
        assert d.method != "ambiguous"
        assert d.method == "marker_scoring_high"
        assert d.cell_type == "Muller_Glia"
        assert d.confidence == "high"


# ── D5 weak-rule arbitration ─────────────────────────────────────────


class TestD5WeakRuleArbitration:
    """D5 arbitration in fuse_all_clusters → fuse_evidence.

    Covers the three Tier-0 branches: corroborated rule hits keep the legacy
    full-confidence early return; uncorroborated single-marker rule hits
    either yield to a strong multi-marker competitor or return at low
    confidence flagged for engine-side review.
    """

    @staticmethod
    def _make_rule(action, markers, corroborators, corroborated):
        """Build a KB expert-rule dict carrying apply_expert_rules metadata."""
        return {
            "priority": 10,
            "action": action,
            "condition": {
                "markers_present": markers,
                "corroborators": corroborators,
            },
            "corroborated": corroborated,
            "corroborators_hit": [],
        }

    @staticmethod
    def _score(value, n_markers, evidence_type, consensus):
        return Score(
            value,
            0.001,
            "hypergeometric",
            n_markers,
            False,
            evidence_type=evidence_type,
            consensus=consensus,
        )

    def test_corroborated_rule_keeps_full_confidence(self) -> None:
        """Given a corroborated rule hit, When fused, Then the cluster keeps
        the legacy full-confidence label (confidence='rule')."""
        all_scores = {
            "0": {"RGC_Alpha": self._score(0.80, 2, "multi_marker", "gold")},
        }
        all_rules = {
            "0": (
                "RGC_Alpha",
                [self._make_rule("RGC_Alpha", {"SPP1": 1.0}, ["RBPMS"], True)],
            ),
        }
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.cell_type == "RGC_Alpha"
        assert d.confidence == "rule"
        assert d.method == "expert_rule"
        assert d.review_reason == ""

    def test_uncorroborated_single_yields_to_strong_competitor(self) -> None:
        """Given an uncorroborated single-marker rule and a top-2 multi_marker
        competitor with consensus >= high, When fused, Then marker scoring wins
        and the rule label is recorded in the explanation as 'weak rule
        overridden' (the cluster is NOT labelled by the rule)."""
        all_scores = {
            "0": {
                "RGC_Alpha": self._score(0.30, 1, "single_marker", "medium"),
                "RGC": self._score(0.80, 3, "multi_marker", "gold"),
            },
        }
        all_rules = {
            "0": (
                "RGC_Alpha",
                [self._make_rule("RGC_Alpha", {"SPP1": 1.0}, ["RBPMS", "NEFH"], False)],
            ),
        }
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.cell_type == "RGC"
        assert d.method == "marker_scoring_high"
        assert d.confidence == "high"
        assert "weak rule overridden" in d.explanation

    def test_uncorroborated_single_without_competitor_low_and_review(self) -> None:
        """Given an uncorroborated single-marker rule with no strong
        multi-marker competitor, When fused, Then the rule label is kept at low
        confidence and flagged review_reason='single_marker_rule' for the
        engine-side review queue."""
        all_scores = {
            "0": {
                "RGC_Alpha": self._score(0.30, 1, "single_marker", "medium"),
                "RGC": self._score(0.45, 1, "single_marker", "low"),
            },
        }
        all_rules = {
            "0": (
                "RGC_Alpha",
                [self._make_rule("RGC_Alpha", {"SPP1": 1.0}, ["RBPMS", "NEFH"], False)],
            ),
        }
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.cell_type == "RGC_Alpha"
        assert d.confidence == "low"
        assert d.method == "expert_rule"
        assert d.review_reason == "single_marker_rule"

    def test_bare_rule_string_without_metadata_stays_legacy(self) -> None:
        """Callers passing a bare rule action (no matched-rule metadata) keep
        the legacy full-confidence behavior: corroboration is undeterminable,
        so the hit is treated as corroborated."""
        all_scores = {
            "0": {"RGC_Alpha": self._score(0.80, 2, "multi_marker", "gold")},
        }
        all_rules = {"0": "RGC_Alpha"}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.cell_type == "RGC_Alpha"
        assert d.confidence == "rule"
        assert d.method == "expert_rule"

    def test_uncorroborated_multi_marker_rule_not_arbitrated(self) -> None:
        """D5 scopes arbitration to SINGLE-marker rules: an uncorroborated
        rule requiring two independent markers keeps the legacy
        full-confidence behavior."""
        all_scores = {
            "0": {"RGC": self._score(0.80, 2, "multi_marker", "gold")},
        }
        all_rules = {
            "0": (
                "RGC",
                [self._make_rule("RGC", {"RBPMS": 1.0, "SNCG": 0.5}, ["NEFL"], False)],
            ),
        }
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.cell_type == "RGC"
        assert d.confidence == "rule"
        assert d.method == "expert_rule"


# ── D2 weak-evidence cap ─────────────────────────────────────────────


class TestD2WeakEvidenceCap:
    """D2 cap in fuse_all_clusters → fuse_evidence.

    Confidence is derived from *evidence strength*, not from the tier
    alone: a marker-scoring winner whose evidence_type is in WEAK_EVIDENCE
    (single-marker / window-padding / weak-multi / zero-evidence) or is
    ai_only is capped at confidence='low' regardless of the tier mapping,
    and review_reason=evidence_type flags the cluster for the engine-side
    review queue.
    """

    @staticmethod
    def _score(value, n_markers, evidence_type, consensus):
        return Score(
            value,
            0.001,
            "hypergeometric",
            n_markers,
            False,
            evidence_type=evidence_type,
            consensus=consensus,
        )

    def test_single_marker_high_tier_capped_to_low(self) -> None:
        """Given a score >= 0.7 (marker_scoring_high) whose winner carries
        evidence_type='single_marker', When fused, Then confidence is forced
        to 'low' and review_reason='single_marker' flags the review queue."""
        all_scores = {
            "0": {"RGC_Alpha": self._score(0.85, 1, "single_marker", "gold")},
        }
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.cell_type == "RGC_Alpha"
        assert d.method == "marker_scoring_high"
        assert d.confidence == "low"
        assert d.review_reason == "single_marker"

    def test_window_padding_high_tier_capped_to_low(self) -> None:
        """Given a high-tier winner whose evidence exists only because the
        window was padded (> 20), When fused, Then confidence is 'low' with
        review_reason='window_padding'."""
        all_scores = {
            "0": {"RGC": self._score(0.80, 2, "window_padding", "medium")},
        }
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.method == "marker_scoring_high"
        assert d.confidence == "low"
        assert d.review_reason == "window_padding"

    def test_strong_multi_marker_not_capped(self) -> None:
        """Given a strong winner (multi_marker + gold consensus) at high
        tier, When fused, Then confidence keeps the tier mapping ('high') and
        no review_reason is set."""
        all_scores = {
            "0": {"RGC": self._score(0.85, 3, "multi_marker", "gold")},
        }
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.method == "marker_scoring_high"
        assert d.confidence == "high"
        assert d.review_reason == ""

    def test_weak_multi_medium_tier_capped_to_low(self) -> None:
        """Given a medium-tier winner with evidence_type='weak_multi', When
        fused, Then confidence is capped to 'low' and review_reason set."""
        all_scores = {
            "0": {"Bipolar": self._score(0.60, 2, "weak_multi", "medium")},
        }
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.method == "marker_scoring_medium"
        assert d.confidence == "low"
        assert d.review_reason == "weak_multi"

    def test_zero_evidence_capped_to_low(self) -> None:
        """zero_evidence is part of WEAK_EVIDENCE: capped at low + review."""
        all_scores = {
            "0": {"Amacrine": self._score(0.75, 0, "zero_evidence", "")},
        }
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.confidence == "low"
        assert d.review_reason == "zero_evidence"

    def test_ai_only_capped_at_tier_level(self) -> None:
        """ai_only joins the capping union: a tier-level Score carrying it
        is capped to low (normally set engine-side by task 9; the union check
        is validated here)."""
        all_scores = {
            "0": {"RGC": self._score(0.80, 2, "ai_only", "gold")},
        }
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.confidence == "low"
        assert d.review_reason == "ai_only"

    def test_bare_float_winner_not_capped(self) -> None:
        """Bare-float score entries (no Score metadata) have empty
        evidence_type → not capped; tier confidence mapping applies."""
        all_scores = {"0": {"T_cell": 0.85}}
        all_rules = {"0": None}
        d = fuse_all_clusters(all_scores, all_rules)[0]

        assert d.confidence == "high"
        assert d.review_reason == ""


# ── Private helper: _is_transition_state ──────────────────────────────


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


# ── KADP developmental-potency naming branch (plan todo 4) ───────────
# A multi-peak ambiguous candidate is named via the KADP potency axis when
# the KB hierarchy (Progenitor pole) dominates the marker scores.  The exit
# semantics are written in stone: a KADP *miss* returns the candidate
# decision byte-identical to the baseline (never falling into the tier loop).
# Each test KB carries both ``_hierarchy`` (for pole derivation) and the type
# entries (for _is_transition_state confirmation).


def _hcat(label, members):
    return {"label": label, "members": list(members), "markers": {}, "subtypes": {}}


_KADP_KB = {
    "_hierarchy": {
        "categories": {
            "Progenitor": _hcat(
                "Progenitor", ["NRPC", "Proliferating_RPC", "Photoreceptor_Precursor"]
            ),
            "Neuron": _hcat("Neuron", ["Rod_Photoreceptor", "Cone_Photoreceptor", "Bipolar_Cell"]),
            "Glia": _hcat("Glia", ["Muller_Glia"]),
            "Non-neural": _hcat("Non-neural", ["RPE"]),
        }
    },
    "Progenitor": {"parent": "", "markers": {"confirm": {"PAX6": ["s1"]}}},
    "NRPC": {
        "parent": "Progenitor",
        "markers": {"confirm": {"VSX2": ["s1"], "PAX6": ["s1"]}},
        "_private_markers": ["VSX2"],
    },
    "Proliferating_RPC": {
        "parent": "Progenitor",
        "markers": {"confirm": {"VSX2": ["s1"], "PCNA": ["s1"]}},
        "_private_markers": ["VSX2"],
    },
    "Photoreceptor": {"parent": "", "markers": {"confirm": {"RHO": ["s1"]}}},
    "Photoreceptor_Precursor": {
        "parent": "Photoreceptor",
        "markers": {"confirm": {"RHO": ["s1"], "NRL": ["s1"]}},
        "_private_markers": ["RHO"],
    },
    "Rod_Photoreceptor": {
        "parent": "Photoreceptor",
        "markers": {"confirm": {"RHO": ["s1"]}},
        "_private_markers": ["RHO"],
    },
    "Cone_Photoreceptor": {
        "parent": "Photoreceptor",
        "markers": {"confirm": {"RHO": ["s1"]}},
        "_private_markers": ["RHO"],
    },
    "Bipolar": {"parent": "", "markers": {"confirm": {"VSX2": ["s1"]}}},
    "Bipolar_Cell": {
        "parent": "Bipolar",
        "markers": {"confirm": {"VSX2": ["s1"], "OTX2": ["s1"]}},
        "_private_markers": ["VSX2"],
    },
    "Muller_Glia": {"parent": "Glia", "markers": {"confirm": {"RLBP1": ["s1"]}}},
    "RPE": {"parent": "", "markers": {"confirm": {"RPE65": ["s1"]}}},
}

# Progenitor-dominant multi-peak: 3 tied types >= 0.9, terminal pole low.
#   -> ambiguous candidate, then KADP names NRPC (ratio ~3.17 >= 2.0).
_KADP_HIT_SCORES = {
    "NRPC": 0.95,
    "Proliferating_RPC": 0.93,
    "Photoreceptor_Precursor": 0.91,
    "Rod_Photoreceptor": 0.30,
}

# Terminal-dominant multi-peak: 3 tied types >= 0.9, progenitor pole weak.
#   -> ambiguous candidate, potency fails (max_prog=0.30 < max_term=0.95),
#      KADP misses and the candidate is returned byte-identical.
_KADP_MISS_SCORES = {
    "Rod_Photoreceptor": 0.95,
    "Cone_Photoreceptor": 0.93,
    "Bipolar_Cell": 0.91,
    "NRPC": 0.30,
}

_KADP_CFG = KADPConfig(enabled=True)


class TestKADPNameBranch:
    """KADP developmental-potency naming in fuse_evidence (plan todo 4).

    Covers the full TDD matrix: hit, miss (exit semantics — candidate
    returned byte-identical to the disabled baseline), disabled, adult
    context, dev_mode trigger, _replace field inheritance and the
    three-value potency dict write.
    """

    # ── KADP hit ─────────────────────────────────────────────────────

    def test_kadp_hit_names_progenitor(self) -> None:
        """High potency ambiguous candidate -> named precursor with
        method=developmental_potency, confidence=medium, review_reason=kadp_precursor.
        """
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=_KADP_CFG,
        )
        assert d.method == "developmental_potency"
        assert d.cell_type == "NRPC"
        assert d.confidence == "medium"
        assert d.review_reason == "kadp_precursor"
        assert d.diagnostic is not None
        assert d.diagnostic.category == "developmental_potency"

    def test_kadp_hit_potency_three_value_dict(self) -> None:
        """potency tail field carries the full {'ratio','abs','gap'} dict."""
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=_KADP_CFG,
        )
        assert d.potency is not None
        assert set(d.potency) == {"ratio", "abs", "gap"}
        assert d.potency["abs"] == pytest.approx(0.95)
        assert d.potency["gap"] == pytest.approx(0.95 - 0.30)
        assert d.potency["ratio"] == pytest.approx(0.95 / 0.30)
        assert d.source_votes is None

    def test_kadp_hit_replace_inherits_candidate_fields(self) -> None:
        """The named decision is built via _replace: score / n_markers_found /
        ai_* fields are inherited from the ambiguous candidate, not recreated.
        """
        base = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            ai_suggestion="Rod Cell",
        )
        assert base.method == "ambiguous"
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            ai_suggestion="Rod Cell",
            kadp_cfg=_KADP_CFG,
        )
        assert d.method == "developmental_potency"
        # Inherited via _replace — identical to the candidate's values.
        assert d.score == base.score == pytest.approx(0.95)
        assert d.n_markers_found == base.n_markers_found == 0
        assert d.ai_agreed == base.ai_agreed is False
        assert d.ai_suggested == base.ai_suggested == "Rod Cell"
        assert d.alternative_rules == base.alternative_rules == []
        # KADP-specific overrides.
        assert d.cell_type == "NRPC"
        assert d.review_reason == "kadp_precursor"

    def test_kadp_hit_top_competitors_dict_form(self) -> None:
        """diagnostic.top_competitors stays in {'cell_type','score'} dict form
        (the same shape the corpus ambiguous diagnostic uses).
        """
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=_KADP_CFG,
        )
        comps = d.diagnostic.top_competitors
        assert len(comps) >= 3
        assert all(isinstance(c, dict) for c in comps)
        assert all({"cell_type", "score"} <= set(c) for c in comps)
        assert {c["cell_type"] for c in comps} == {
            "NRPC",
            "Proliferating_RPC",
            "Photoreceptor_Precursor",
        }

    # ── KADP miss — exit semantics (Momus r2 MAJOR-1) ─────────────────

    def test_kadp_miss_returns_candidate_byte_identical(self) -> None:
        """kadp_cfg.enabled=True with a low-potency candidate must return
        the SAME ambiguous decision as the disabled baseline — full-field
        equality including the diagnostic (no drift, no tier-loop fallthrough).
        """
        disabled = fuse_evidence(
            marker_scores=dict(_KADP_MISS_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
        )
        enabled = fuse_evidence(
            marker_scores=dict(_KADP_MISS_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=_KADP_CFG,
        )
        assert enabled == disabled
        assert enabled.method == "ambiguous"
        assert enabled.cell_type == "Unknown"
        assert enabled.confidence == "unknown"
        assert enabled.diagnostic is not None
        assert enabled.diagnostic.category == "ambiguous"
        assert enabled.potency is None
        assert enabled.source_votes is None

    def test_kadp_miss_never_reaches_tier_loop(self) -> None:
        """A KADP miss must not fall into the tier loop: the ambiguous
        candidate keeps its ambiguous method (no marker_scoring label).
        """
        d = fuse_evidence(
            marker_scores=dict(_KADP_MISS_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=_KADP_CFG,
        )
        assert d.method == "ambiguous"
        assert not d.method.startswith("marker_scoring")

    # ── Disabled / adult / dev_mode ──────────────────────────────────

    def test_kadp_disabled_default_params(self) -> None:
        """kadp_cfg=None (default) keeps baseline behavior: the same input
        stays ambiguous without any developmental naming.
        """
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
        )
        assert d.method == "ambiguous"
        assert d.cell_type == "Unknown"
        assert d.review_reason == ""
        assert d.potency is None

    def test_kadp_disabled_explicit_enabled_false(self) -> None:
        """KADPConfig(enabled=False) behaves exactly like kadp_cfg=None."""
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=KADPConfig(enabled=False),
        )
        assert d.method == "ambiguous"
        assert d.potency is None

    def test_adult_context_no_kadp_trigger(self) -> None:
        """allows_transitions=False (adult tissue) disables the whole
        transition block — KADP never fires; normal tier logic applies.
        """
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=False,
            kadp_cfg=_KADP_CFG,
        )
        assert d.method == "marker_scoring_high"
        assert d.cell_type == "NRPC"
        assert d.review_reason == ""
        assert d.potency is None

    def test_dev_mode_transitions_trigger_kadp(self) -> None:
        """engine.py sets allows_transitions = bool(transition_context) or
        _dev_mode (marker.developmental_mode).  Through that single flag
        dev_mode activates the KADP branch: naming fires.
        """
        d = fuse_evidence(
            marker_scores=dict(_KADP_HIT_SCORES),
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,  # == dev-mode activation path
            kadp_cfg=_KADP_CFG,
        )
        assert d.method == "developmental_potency"
        assert d.cell_type == "NRPC"

    # ── transition_state candidate is NOT eligible for KADP ──────────

    def test_transition_candidate_not_kadp_eligible(self) -> None:
        """KADP names only ambiguous (multi-peak) candidates.  A genuine
        transition_state candidate stays untouched even with kadp enabled.
        """
        scores = {"NRPC": 0.60, "Proliferating_RPC": 0.55, "Rod_Photoreceptor": 0.10}
        d = fuse_evidence(
            marker_scores=scores,
            expert_rule_result=None,
            kb=_KADP_KB,
            allows_transitions=True,
            kadp_cfg=_KADP_CFG,
        )
        assert d.method == "transition_state"
        assert d.cell_type == "transitional: NRPC/Proliferating_RPC"
        assert d.potency is None
        assert d.review_reason == ""


# ═══════════════════════════════════════════════════════════════════════
#  Todo 8 — shared label harmonization (parallel A/B evaluation)
#  harmonize_label + fuse_all_clusters celltypist forwarding + rate
# ═══════════════════════════════════════════════════════════════════════


_HARMONIZE_KB = {
    "RGC": {
        "markers": {"confirm": {"RBPMS": ["s"], "SNCG": ["s"]}, "add": {}},
        "negative_markers": [],
    },
    "Muller_Glia": {
        "markers": {"confirm": {"RLBP1": ["s"]}, "add": {}},
        "negative_markers": [],
    },
    "Rod_Photoreceptor": {
        "markers": {"confirm": {"RHO": ["s"]}, "add": {}},
        "negative_markers": [],
    },
    "Amacrine_Cell": {
        "markers": {"confirm": {"TFAP2A": ["s"]}, "add": {}},
        "negative_markers": [],
    },
    "Amacrine_Precursor": {
        "markers": {"confirm": {"BARHL2": ["s"]}, "add": {}},
        "negative_markers": [],
    },
    "expert_rules": {},
    "_hierarchy": {},
}


_HARMONIZE_SYNONYMS = {
    "RGC": {
        "display_name": "Retinal Ganglion Cell",
        "synonyms": ["Retinal Ganglion Cell", "RGC", "rgc", "ganglion cell"],
    },
    "Muller_Glia": {
        "display_name": "Müller Glia",
        "synonyms": ["Muller Glia", "Müller Glia"],
    },
    "Rod_Photoreceptor": {
        "display_name": "Rod Photoreceptor",
        "synonyms": ["Rod Photoreceptor", "rod photoreceptor"],
    },
    "Amacrine_Cell": {
        "display_name": "Amacrine Cell",
        "synonyms": ["Amacrine Cell", "Amacrine", "AC"],
    },
    "Amacrine_Precursor": {
        "display_name": "Amacrine Precursor",
        "synonyms": ["AC", "amacrine precursor"],
    },
}


class TestHarmonizeLabel:
    """todo 8: shared harmonize_label chain — parallel A/B evaluation.

    Path A = KB type-key exact match (normalized); Path B = synonyms reverse
    lookup.  Both are evaluated independently; a conflict between the two
    (different candidate sets) abstains with ``None`` (Oracle r3 MAJOR 3 —
    sequential evaluation would let 'RPC' short-circuit on the KB key and
    never reach the synonyms path).
    """

    # ── Single-path hits ──────────────────────────────────────────

    def test_synonym_hit_resolves_canonical(self) -> None:
        """Path B only: 'Retinal Ganglion Cell' is an RGC synonym."""
        assert (
            harmonize_label("Retinal Ganglion Cell", _HARMONIZE_KB, _HARMONIZE_SYNONYMS) == "RGC"
        )

    def test_normalization_hit_resolves_kb_key(self) -> None:
        """Path A: 'Müller Glia' normalises (NFKD) to the KB key and both
        paths agree on 'Muller_Glia'."""
        assert harmonize_label("Müller Glia", _HARMONIZE_KB, _HARMONIZE_SYNONYMS) == "Muller_Glia"
        # punctuation/case-insensitive path-A match, paths agree
        assert (
            harmonize_label("Rod-Photoreceptor", _HARMONIZE_KB, _HARMONIZE_SYNONYMS)
            == "Rod_Photoreceptor"
        )

    def test_unresolvable_label_discarded(self) -> None:
        """No KB key and no synonym list contains the label → None."""
        assert harmonize_label("Foobar Cell", _HARMONIZE_KB, _HARMONIZE_SYNONYMS) is None

    # ── Abstain paths ─────────────────────────────────────────────

    def test_ambiguity_synonym_tie_abstains(self) -> None:
        """Path B alone: 'AC' is a synonym of both Amacrine_Cell and
        Amacrine_Precursor, and their KB marker counts tie (1 each) → None.
        """
        assert harmonize_label("AC", _HARMONIZE_KB, _HARMONIZE_SYNONYMS) is None

    def test_a_b_conflict_abstains(self) -> None:
        """'RPC' is BOTH a KB type key and a Broad_Progenitor synonym.

        Sequential evaluation would hit the KB key first and never reach the
        synonyms path — parallel evaluation sees A={RPC} vs B={RPC,
        Broad_Progenitor} and abstains (Oracle r3 MAJOR 3).
        """
        kb = {
            "RPC": {"markers": {"confirm": {"A": ["s"], "B": ["s"]}, "add": {}}},
            "NRPC": {"markers": {"confirm": {"C": ["s"]}, "add": {}}},
            "expert_rules": {},
        }
        synonyms = {
            "RPC": {
                "display_name": "Retinal Progenitor Cell",
                "synonyms": ["RPC", "rpc", "retinal progenitor cell"],
            },
            "Broad_Progenitor": {
                "display_name": "Progenitor",
                "synonyms": ["Progenitor", "RPC"],
            },
        }
        assert harmonize_label("RPC", kb, synonyms) is None

    def test_empty_or_none_label_returns_none(self) -> None:
        """Empty/None labels never resolve."""
        assert harmonize_label("", _HARMONIZE_KB, _HARMONIZE_SYNONYMS) is None
        assert harmonize_label(None, _HARMONIZE_KB, _HARMONIZE_SYNONYMS) is None

    def test_missing_kb_or_synonyms_tolerated(self) -> None:
        """None kb / None synonyms must not crash — the synonym path still
        resolves without a kb and path A still resolves without synonyms."""
        # 'AC' needs the KB for the marker-count tie-break → None, no crash
        assert harmonize_label("AC", None, _HARMONIZE_SYNONYMS) is None
        assert harmonize_label("Rod_Photoreceptor", _HARMONIZE_KB, None) == "Rod_Photoreceptor"

    # ── Real retina data (the mandated ambiguity pair) ────────────

    def test_real_retina_rpc_abstains_and_aliases_resolve(self) -> None:
        """Real retina KB + synonyms: 'RPC' (KB key + Broad_Progenitor
        synonym, synonyms.yaml L608-618) abstains; well-known aliases resolve.
        """
        from core.kb import load_kb, load_synonyms

        kb = load_kb("retina")
        synonyms = load_synonyms("retina")
        assert harmonize_label("RPC", kb, synonyms) is None
        assert harmonize_label("Retinal Ganglion Cell", kb, synonyms) == "RGC"
        assert harmonize_label("Rod Photoreceptor", kb, synonyms) == "Rod_Photoreceptor"
        assert harmonize_label("Müller Glia", kb, synonyms) == "Muller_Glia"
        assert harmonize_label("definitely-not-a-cell", kb, synonyms) is None


class TestFuseAllClustersHarmonization:
    """todo 8: celltypist labels are harmonized per cluster, forwarded as
    ``celltypist_suggestion``, and contribute ``harmonization_rate`` to the
    quality dict (None when no raw labels exist — zero-division guard).
    """

    _ALL_SCORES = {
        "0": {"RGC": 0.85},
        "1": {"Amacrine_Cell": 0.60},
        "2": {"Muller_Glia": 0.90},
    }
    _ALL_RULES = {"0": None, "1": None, "2": None}

    def test_celltypist_harmonized_and_forwarded(self) -> None:
        """Resolvable CellTypist labels are forwarded as the canonical
        celltypist_suggestion; unresolvable ones abstain (None)."""
        from unittest.mock import patch

        celltypist_results = {
            "0": "Retinal Ganglion Cell",  # → RGC (synonym hit)
            "1": "Amacrine",  # → Amacrine_Cell (synonym hit)
            "2": "Complete nonsense",  # → None (abstain)
        }
        captured: list = []
        real = fuse_evidence

        def _spy(*args, **kwargs):
            captured.append(kwargs.get("celltypist_suggestion"))
            return real(*args, **kwargs)

        with patch("rna.utils.evidence_fusion.fuse_evidence", side_effect=_spy):
            decisions, quality = fuse_all_clusters(
                dict(self._ALL_SCORES),
                dict(self._ALL_RULES),
                kb=_HARMONIZE_KB,
                celltypist_results=dict(celltypist_results),
                synonyms=_HARMONIZE_SYNONYMS,
                return_quality=True,
            )
        assert len(decisions) == 3
        assert captured == ["RGC", "Amacrine_Cell", None], captured
        assert quality["celltypist"] is True
        assert quality["harmonization_rate"] == pytest.approx(2 / 3)

    def test_harmonization_rate_none_when_no_raw_labels(self) -> None:
        """Empty/None celltypist_results → harmonization_rate None (never
        a division by zero), and quality['celltypist'] stays False."""
        _, q_empty = fuse_all_clusters(
            dict(self._ALL_SCORES),
            dict(self._ALL_RULES),
            kb=_HARMONIZE_KB,
            celltypist_results={},
            synonyms=_HARMONIZE_SYNONYMS,
            return_quality=True,
        )
        assert q_empty["harmonization_rate"] is None
        assert q_empty["celltypist"] is False

        _, q_none = fuse_all_clusters(
            dict(self._ALL_SCORES),
            dict(self._ALL_RULES),
            kb=_HARMONIZE_KB,
            celltypist_results=None,
            synonyms=_HARMONIZE_SYNONYMS,
            return_quality=True,
        )
        assert q_none["harmonization_rate"] is None

    def test_celltypist_all_unresolvable_rate_zero(self) -> None:
        """Raw labels exist but none align to the KB → rate 0.0 (not None),
        and every cluster's celltypist_suggestion abstains."""
        from unittest.mock import patch

        captured: list = []
        real = fuse_evidence

        def _spy(*args, **kwargs):
            captured.append(kwargs.get("celltypist_suggestion"))
            return real(*args, **kwargs)

        with patch("rna.utils.evidence_fusion.fuse_evidence", side_effect=_spy):
            _, quality = fuse_all_clusters(
                dict(self._ALL_SCORES),
                dict(self._ALL_RULES),
                kb=_HARMONIZE_KB,
                celltypist_results={"0": "zzz", "1": "yyy"},
                synonyms=_HARMONIZE_SYNONYMS,
                return_quality=True,
            )
        assert quality["harmonization_rate"] == 0.0
        # one fuse_evidence call per cluster (3) — every cluster abstains
        assert captured == [None, None, None]
