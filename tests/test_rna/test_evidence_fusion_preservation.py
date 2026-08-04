"""Behavior-preservation corpus for rna/utils/evidence_fusion.py (plan todo 1, F2).

Proves that `fuse_evidence` with ALL default parameters (no kadp/metc config,
no celltypist, no AI) produces decisions byte-identical to the frozen
baseline — i.e. the append-only FusionDecision extension (potency /
source_votes tail fields) and the new METCConfig dataclass introduce zero
behavior drift.  Expected values are hard-coded from a baseline capture of
the unmodified implementation.

Coverage required by the plan (todo 1b):
  * multi-peak tie        -> ambiguous
  * transition candidate  -> transition_state
  * low-score             -> marker_scoring_low
  * expert-rule           -> rule (corroborated + uncorroborated single)
  * empty scores          -> unknown
  * no kb                 -> transition detection disabled
  * plus high/medium, weak-evidence cap, low-quality, AI agree/disagree.
"""

import dataclasses

import pytest

from core.annotation.scoring import Score
from rna.utils.evidence_fusion import (
    DiagnosticInfo,
    FusionDecision,
    METCConfig,
    fuse_evidence,
)

# ── Representative retina-ish KB for transition detection ─────────────────
KB = {
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
    "Rod_Photoreceptor": {
        "parent": "Photoreceptor",
        "markers": {"confirm": {"RHO": ["s1"], "NRL": ["s1"]}},
        "_private_markers": ["RHO"],
    },
    "Ganglion_Cell": {"parent": "", "markers": {"confirm": {"RBPMS": ["s1"]}}},
    "RGC": {
        "parent": "Ganglion_Cell",
        "markers": {"confirm": {"RBPMS": ["s1"], "POU4F1": ["s1"]}},
        "_private_markers": ["RBPMS"],
    },
    "Bipolar": {"parent": "", "markers": {"confirm": {"VSX2": ["s1"]}}},
    "Bipolar_Cell": {
        "parent": "Bipolar",
        "markers": {"confirm": {"VSX2": ["s1"], "OTX2": ["s1"]}},
        "_private_markers": ["VSX2"],
    },
}

_SINGLE_MARKER_RULE = [{"action": "RGC", "condition": {"markers_present": {"RBPMS": ["RBPMS"]}}}]
_UNCORROBORATED_SINGLE_RULE = [
    {
        "action": "RGC",
        "corroborated": False,
        "condition": {"markers_present": {"RBPMS": ["RBPMS"]}},
    }
]


# ═══════════════════════════════════════════════════════════════════════
#  F2 — behavior-preservation corpus (fuse_evidence, all defaults)
# ═══════════════════════════════════════════════════════════════════════
# Each corpus entry: (name, fuse_evidence kwargs, hard-coded expected decision).
# Expected objects are written with the SAME default tail fields
# (potency=None, source_votes=None); equality is asserted field-by-field.
CORPUS = [
    pytest.param(
        dict(
            marker_scores={"NRPC": 0.95, "Proliferating_RPC": 0.93, "Rod_Photoreceptor": 0.91},
            expert_rule_result=None,
            kb=KB,
            allows_transitions=True,
        ),
        FusionDecision(
            cell_type="Unknown",
            confidence="unknown",
            score=0.95,
            method="ambiguous",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="Marker scoring selected Unknown with score 0.950",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category="ambiguous",
                top_competitors=[
                    {"cell_type": "NRPC", "score": 0.95},
                    {"cell_type": "Proliferating_RPC", "score": 0.93},
                    {"cell_type": "Rod_Photoreceptor", "score": 0.91},
                ],
                detail=(
                    "Multi-peak tie: 3 types >= 0.9 "
                    "(top: NRPC=0.950, Proliferating_RPC=0.930, Rod_Photoreceptor=0.910)"
                ),
            ),
        ),
        id="multi_peak_tie_ambiguous",
    ),
    pytest.param(
        dict(
            marker_scores={"NRPC": 0.60, "Proliferating_RPC": 0.55, "Rod_Photoreceptor": 0.10},
            expert_rule_result=None,
            kb=KB,
            allows_transitions=True,
        ),
        FusionDecision(
            cell_type="transitional: NRPC/Proliferating_RPC",
            confidence="transition",
            score=0.60,
            method="transition_state",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation=(
                "Transitional state detected: transitional: NRPC/Proliferating_RPC\n"
                "  Confidence: transition  Score: 0.6000\n"
                "  Top-2 scores within 0.050, shared lineage Progenitor\n"
                "  NRPC: 0.600\n"
                "  Proliferating_RPC: 0.550"
            ),
            alternative_rules=[],
            diagnostic=None,
        ),
        id="transition_candidate",
    ),
    pytest.param(
        dict(
            marker_scores={"RGC": 0.85, "Bipolar_Cell": 0.30},
            expert_rule_result=None,
            kb=KB,
            allows_transitions=True,
        ),
        FusionDecision(
            cell_type="RGC",
            confidence="high",
            score=0.85,
            method="marker_scoring_high",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="Marker scoring selected RGC with score 0.850",
            alternative_rules=[],
            diagnostic=None,
        ),
        id="high_confidence_no_transition",
    ),
    pytest.param(
        dict(marker_scores={"RGC": 0.60, "Bipolar_Cell": 0.40}, expert_rule_result=None),
        FusionDecision(
            cell_type="RGC",
            confidence="medium",
            score=0.60,
            method="marker_scoring_medium",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="Marker scoring selected RGC with score 0.600",
            alternative_rules=[],
            diagnostic=None,
        ),
        id="medium_confidence",
    ),
    pytest.param(
        dict(marker_scores={"Bipolar_Cell": 0.35}, expert_rule_result=None),
        FusionDecision(
            cell_type="Bipolar_Cell",
            confidence="low",
            score=0.35,
            method="marker_scoring_low",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="Marker scoring selected Bipolar_Cell with score 0.350",
            alternative_rules=[],
            diagnostic=None,
        ),
        id="low_score",
    ),
    pytest.param(
        dict(
            marker_scores={"RGC": 0.80, "Bipolar_Cell": 0.30},
            expert_rule_result="RGC",
            alternative_rules=_SINGLE_MARKER_RULE,
        ),
        FusionDecision(
            cell_type="RGC",
            confidence="rule",
            score=0.80,
            method="expert_rule",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="Expert rule matched: RGC (marker score: 0.800)",
            alternative_rules=_SINGLE_MARKER_RULE,
            diagnostic=None,
        ),
        id="expert_rule_corroborated",
    ),
    pytest.param(
        dict(
            marker_scores={"RGC": 0.30, "Bipolar_Cell": 0.20},
            expert_rule_result="RGC",
            alternative_rules=_UNCORROBORATED_SINGLE_RULE,
        ),
        FusionDecision(
            cell_type="RGC",
            confidence="low",
            score=0.30,
            method="expert_rule",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation=(
                "Expert rule matched 'RGC' but its corroborators were absent from the "
                "DE top-N subset. Uncorroborated single-marker rule — labelled at low "
                "confidence. (marker score: 0.300)"
            ),
            alternative_rules=_UNCORROBORATED_SINGLE_RULE,
            diagnostic=None,
            review_reason="single_marker_rule",
        ),
        id="expert_rule_uncorroborated_single",
    ),
    pytest.param(
        dict(marker_scores={}, expert_rule_result=None),
        FusionDecision(
            cell_type="Unknown",
            confidence="unknown",
            score=0.0,
            method="unknown",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="No marker scores available for this cluster.",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category="true_unknown",
                top_competitors=[],
                detail="No marker scores calculated — empty or missing data.",
            ),
        ),
        id="empty_scores",
    ),
    pytest.param(
        dict(marker_scores={"Bipolar_Cell": 0.15}, expert_rule_result=None),
        FusionDecision(
            cell_type="Unknown",
            confidence="unknown",
            score=0.15,
            method="unknown",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="No cell type could be confidently assigned (best match: Bipolar_Cell, score: 0.150)",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category="weak_signal",
                top_competitors=[("Bipolar_Cell", 0.15)],
                detail="Best score 0.1500 below 0.25 threshold (best type: Bipolar_Cell)",
            ),
        ),
        id="weak_signal_below_threshold",
    ),
    pytest.param(
        dict(
            marker_scores={"NRPC": 0.60, "Proliferating_RPC": 0.55},
            expert_rule_result=None,
            kb=None,
            allows_transitions=True,
        ),
        FusionDecision(
            cell_type="NRPC",
            confidence="medium",
            score=0.60,
            method="marker_scoring_medium",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="Marker scoring selected NRPC with score 0.600",
            alternative_rules=[],
            diagnostic=None,
        ),
        id="no_kb_disables_transition",
    ),
    pytest.param(
        dict(
            marker_scores={
                "RGC": Score(
                    score=0.85,
                    p_value=0.001,
                    method="hypergeometric",
                    n_markers_found=1,
                    negative_penalty=False,
                    evidence_type="single_marker",
                )
            },
            expert_rule_result=None,
        ),
        FusionDecision(
            cell_type="RGC",
            confidence="low",
            score=0.85,
            method="marker_scoring_high",
            n_markers_found=1,
            ai_agreed=False,
            ai_suggested="",
            explanation="Marker scoring selected RGC with score 0.850 (1 KB markers found in cluster top-20)",
            alternative_rules=[],
            diagnostic=None,
            review_reason="single_marker",
        ),
        id="weak_evidence_cap_single_marker",
    ),
    pytest.param(
        dict(
            marker_scores={"RGC": 0.15, "Bipolar_Cell": 0.10},
            expert_rule_result=None,
            low_quality_reason="low_marker_overlap",
        ),
        FusionDecision(
            cell_type="Unknown",
            confidence="unknown",
            score=0.15,
            method="unknown",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="No cell type could be confidently assigned (best match: RGC, score: 0.150)",
            alternative_rules=[],
            diagnostic=DiagnosticInfo(
                category="low_quality_data",
                top_competitors=[("RGC", 0.15), ("Bipolar_Cell", 0.10)],
                detail="Cluster flagged as low-quality: low_marker_overlap",
            ),
        ),
        id="low_quality_cluster",
    ),
    pytest.param(
        dict(
            marker_scores={"RGC": 0.85, "Bipolar_Cell": 0.30},
            expert_rule_result=None,
            ai_suggestion="Bipolar Cell",
        ),
        FusionDecision(
            cell_type="RGC",
            confidence="high",
            score=0.85,
            method="marker_scoring_high",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="Bipolar Cell",
            explanation=(
                "Marker scoring selected RGC with score 0.850 — "
                "AI suggested 'Bipolar Cell' (different from marker-based result)"
            ),
            alternative_rules=[],
            diagnostic=None,
        ),
        id="ai_disagreement",
    ),
    pytest.param(
        dict(
            marker_scores={"RGC": 0.85, "Bipolar_Cell": 0.30},
            expert_rule_result=None,
            ai_suggestion="RGC",
        ),
        FusionDecision(
            cell_type="RGC",
            confidence="high",
            score=0.85,
            method="marker_scoring_high",
            n_markers_found=0,
            ai_agreed=True,
            ai_suggested="RGC",
            explanation="Marker scoring selected RGC with score 0.850 — AI agreed with this assignment",
            alternative_rules=[],
            diagnostic=None,
        ),
        id="ai_agreement",
    ),
]


class TestCorpusPreservation:
    """Full-field equality between default-parameter decisions and the frozen baseline."""

    @pytest.mark.parametrize(("kwargs", "expected"), CORPUS)
    def test_decision_matches_baseline(self, kwargs: dict, expected: FusionDecision) -> None:
        decision = fuse_evidence(**kwargs)

        # All 11 baseline fields (incl. nested diagnostic) equal, byte-for-byte.
        assert decision == expected

        # The new append-only tail fields must stay None under default config
        # (kadp_cfg=None / metc_cfg=None / no celltypist / no ai).
        assert decision.potency is None
        assert decision.source_votes is None

    def test_corpus_has_required_scenarios(self) -> None:
        """The corpus covers every scenario the plan mandates (>= 10 inputs)."""
        ids = {p.id for p in CORPUS}
        assert len(CORPUS) >= 10
        for required in (
            "multi_peak_tie_ambiguous",
            "transition_candidate",
            "low_score",
            "expert_rule_corroborated",
            "empty_scores",
            "no_kb_disables_transition",
        ):
            assert required in ids, f"missing required corpus scenario: {required}"


# ═══════════════════════════════════════════════════════════════════════
#  F8 — FusionDecision append-only tail fields (field-order regression)
# ═══════════════════════════════════════════════════════════════════════
class TestFusionDecisionTailFields:
    def test_tail_fields_appended_in_order(self) -> None:
        """potency/source_votes must be the LAST two fields (no mid insertion)."""
        assert list(FusionDecision._fields)[-2:] == ["potency", "source_votes"]

    def test_9_positional_fallback_defaults(self) -> None:
        """The 9-positional fallback still constructs and the tail fields default None.

        Mirrors the exact fallback literal at evidence_fusion.py:896-898.
        """
        d = FusionDecision(
            "Unknown", "unknown", 0.0, "unknown", 0, False, "", "Fallback: no tier matched.", []
        )
        assert d.cell_type == "Unknown"
        assert d.confidence == "unknown"
        assert d.method == "unknown"
        assert d.potency is None
        assert d.source_votes is None

    def test_explicit_tail_field_values(self) -> None:
        """Callers may set potency / source_votes via keyword (used by KADP/METC later)."""
        d = FusionDecision(
            cell_type="Unknown",
            confidence="unknown",
            score=0.5,
            method="ambiguous",
            n_markers_found=0,
            ai_agreed=False,
            ai_suggested="",
            explanation="x",
            alternative_rules=[],
            potency={"ratio": 3.0, "abs": 0.9, "gap": 0.6},
            source_votes={"marker": "NRPC", "expert": None, "ai": None, "celltypist": None},
        )
        assert d.potency == {"ratio": 3.0, "abs": 0.9, "gap": 0.6}
        assert d.source_votes == {"marker": "NRPC", "expert": None, "ai": None, "celltypist": None}


# ═══════════════════════════════════════════════════════════════════════
#  METCConfig dataclass (todo 1c) — defaults locked for the plan
# ═══════════════════════════════════════════════════════════════════════
class TestMETCConfig:
    def test_is_dataclass(self) -> None:
        assert dataclasses.is_dataclass(METCConfig)

    def test_defaults(self) -> None:
        cfg = METCConfig()
        assert cfg.enabled is False
        assert cfg.min_sources == 3
        assert cfg.min_distinct_transition == 3

    def test_explicit_construction(self) -> None:
        cfg = METCConfig(enabled=True, min_sources=4, min_distinct_transition=2)
        assert cfg.enabled is True
        assert cfg.min_sources == 4
        assert cfg.min_distinct_transition == 2

    def test_partial_override_keeps_defaults(self) -> None:
        cfg = METCConfig(min_sources=2)
        assert cfg.enabled is False
        assert cfg.min_sources == 2
        assert cfg.min_distinct_transition == 3
