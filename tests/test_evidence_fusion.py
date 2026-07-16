"""Unit tests for _is_transition_state in rna/utils/evidence_fusion.py."""

import pytest

from rna.utils.evidence_fusion import _is_transition_state


class TestIsTransitionState:
    """Edge-case coverage for the transition-state detection function."""

    # ── Happy path ───────────────────────────────────────────────────

    def test_transition_normal(self) -> None:
        """delta < threshold, same parent, top score >= 0.25 → (top1, top2)."""
        scores = {"RGC": 0.45, "Amacrine_Cell": 0.35}
        kb = {
            "RGC": {"parent": "Broad_Neuron"},
            "Amacrine_Cell": {"parent": "Broad_Neuron"},
        }
        result = _is_transition_state(scores, kb)
        assert result is not None
        # Sorted descending: RGC (0.45) first, Amacrine_Cell (0.35) second
        assert result[0] == "RGC"
        assert result[1] == "Amacrine_Cell"

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
