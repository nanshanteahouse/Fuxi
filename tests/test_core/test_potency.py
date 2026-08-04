"""TDD tests for ``core/annotation/potency.py`` — KADP developmental potency.

Three potency variants (ratio / abs-with-guard / gap) are computed over the
progenitor vs terminal poles derived from the KB hierarchy.  All functions
are pure — no IO, no rna-layer imports.
"""

import math

import pytest
import yaml

from core.annotation.potency import (
    KADPConfig,
    PotencyResult,
    compute_potency,
    derive_developmental_poles,
    evaluate_passes,
    filter_pole_scores,
    to_potency_dict,
)
from core.annotation.scoring import Score

# ── Helpers ────────────────────────────────────────────────────────────


def _cat(label, members):
    return {"label": label, "members": list(members), "markers": {}, "subtypes": {}}


def _mini_kb(
    prog=("RPC", "NRPC", "Photoreceptor_Precursor"),
    neuron=("Rod_Photoreceptor", "Bipolar_Cell"),
    glia=("Muller_Glia",),
    non_neural=("RPE",),
):
    return {
        "_hierarchy": {
            "categories": {
                "Progenitor": _cat("Progenitor", prog),
                "Neuron": _cat("Neuron", neuron),
                "Glia": _cat("Glia", glia),
                "Non-neural": _cat("Non-neural", non_neural),
            }
        }
    }


def _s(score):
    return Score(score, 1.0, "none", 0, False)


# ── KADPConfig defaults ────────────────────────────────────────────────


def test_kadp_config_defaults():
    cfg = KADPConfig()
    assert cfg.enabled is False
    assert cfg.ratio_threshold == 2.0
    assert cfg.abs_threshold == 0.6
    assert cfg.gap_threshold == 0.1
    assert cfg.use_gap_criterion is False
    assert cfg.epsilon == 1e-9


def test_kadp_config_override_and_immutable():
    cfg = KADPConfig(enabled=True, ratio_threshold=3.0, use_gap_criterion=True)
    assert cfg.enabled is True
    assert cfg.ratio_threshold == 3.0
    assert cfg.use_gap_criterion is True
    with pytest.raises(Exception):
        cfg.ratio_threshold = 5.0  # frozen dataclass


# ── derive_developmental_poles ─────────────────────────────────────────


def test_derive_poles_from_categories():
    prog, term = derive_developmental_poles(_mini_kb())
    assert prog == {"RPC", "NRPC", "Photoreceptor_Precursor"}
    assert term == {"Rod_Photoreceptor", "Bipolar_Cell", "Muller_Glia", "RPE"}


def test_derive_poles_missing_hierarchy_returns_empty():
    for kb in ({}, None, {"markers": {}}):
        prog, term = derive_developmental_poles(kb)
        assert prog == set()
        assert term == set()


def test_derive_poles_partial_categories():
    # Glia / Non-neural absent → terminal pole = Neuron only; no crash.
    kb = {
        "_hierarchy": {
            "categories": {
                "Progenitor": _cat("Progenitor", ["RPC"]),
                "Neuron": _cat("Neuron", ["RGC"]),
            }
        }
    }
    prog, term = derive_developmental_poles(kb)
    assert prog == {"RPC"}
    assert term == {"RGC"}


def test_derive_poles_retina_hierarchy_contract():
    """Pin the real retina hierarchy: 15 progenitor members incl. the MG
    reprogramming entries, and terminal = Neuron ∪ Glia ∪ Non-neural."""
    with open("core/kb/retina/hierarchy.yaml", encoding="utf-8") as fh:
        hierarchy = yaml.safe_load(fh)
    prog, term = derive_developmental_poles(
        {"_hierarchy": {"categories": hierarchy["categories"]}}
    )

    expected_prog = {
        "RPC",
        "Proliferating_RPC",
        "Neonatal_RPCs",
        "NRPC",
        "NRPC_RGC_fate",
        "NRPC_AC_HC_fate",
        "NRPC_Cone_BC_fate",
        "NRPC_Rod_fate",
        "Photoreceptor_Precursor",
        "Amacrine_Precursor",
        "Bipolar_Precursor",
        "Developing_AC_HC_Precursors",
        "Developing_BC_Photo_Precursors",
        "ASCL1_Reprogrammed_MG",
        "Proliferating_MG",
    }
    assert prog == expected_prog
    assert "Muller_Glia" in term
    assert "Rod_Photoreceptor" in term
    assert "RPE" in term
    assert not (prog & term)  # poles disjoint


# ── filter_pole_scores ─────────────────────────────────────────────────


def test_filter_pole_scores_drops_ghost_members():
    pole = {"RPC", "NRPC"}  # NRPC is a ghost — absent from marker_scores
    filtered = filter_pole_scores(pole, {"RPC": 0.5})
    assert filtered == {"RPC": 0.5}


def test_filter_pole_scores_positive_only():
    pole = {"RPC", "NRPC", "Photoreceptor_Precursor"}
    marker_scores = {
        "RPC": 0.5,
        "NRPC": 0.0,  # score == 0 → dropped
        "Photoreceptor_Precursor": _s(0.0),  # zero-scored Score object → dropped
    }
    filtered = filter_pole_scores(pole, marker_scores)
    assert filtered == {"RPC": 0.5}


def test_filter_pole_scores_accepts_score_objects():
    pole = {"RPC", "NRPC"}
    marker_scores = {"RPC": _s(0.7), "NRPC": 0.3}  # mix Score objects + raw floats
    filtered = filter_pole_scores(pole, marker_scores)
    assert filtered == {"RPC": 0.7, "NRPC": 0.3}


def test_filter_pole_scores_empty_inputs():
    assert filter_pole_scores(set(), {"RPC": 0.5}) == {}
    assert filter_pole_scores({"RPC"}, {}) == {}
    assert filter_pole_scores(set(), {}) == {}


# ── compute_potency — three variants ───────────────────────────────────


def test_compute_potency_happy_progenitor_dominant():
    cfg = KADPConfig()
    res = compute_potency(
        {"RPC": 0.9, "Rod_Photoreceptor": 0.3},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert math.isclose(res.max_prog, 0.9)
    assert math.isclose(res.max_term, 0.3)
    assert math.isclose(res.ratio, 3.0)
    assert math.isclose(res.gap, 0.6)
    assert res.best_progenitor_type == "RPC"
    assert res.passes_ratio is True
    assert res.passes_gap is True
    assert res.passes_abs is True
    assert evaluate_passes(res, cfg) is True


def test_compute_potency_saturated_abs_guard():
    """Saturated 1.0/1.0: ratio=1.0 fails, gap=0 fails, and the abs guard
    (max_prog > max_term) keeps passes_abs False."""
    cfg = KADPConfig()
    res = compute_potency(
        {"RPC": 1.0, "Rod_Photoreceptor": 1.0},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert math.isclose(res.ratio, 1.0)
    assert math.isclose(res.gap, 0.0)
    assert res.passes_ratio is False
    assert res.passes_gap is False
    assert res.passes_abs is False  # abs guard fires even though 1.0 >= 0.6
    assert evaluate_passes(res, cfg) is False


def test_compute_potency_div_by_zero_uses_epsilon():
    """max_term == 0 (terminal pole empty/zero) → ratio = max_prog / epsilon."""
    cfg = KADPConfig()
    res = compute_potency(
        {"RPC": 0.5},
        ({"RPC"}, set()),
        cfg,
    )
    assert res.max_term == 0.0
    assert math.isclose(res.ratio, 0.5 / cfg.epsilon)
    assert res.passes_ratio is True


def test_compute_potency_double_zero():
    cfg = KADPConfig()
    res = compute_potency({}, (set(), set()), cfg)
    assert res.max_prog == 0.0
    assert res.max_term == 0.0
    assert res.ratio == 0.0
    assert res.gap == 0.0
    assert res.best_progenitor_type is None
    assert res.passes_ratio is False
    assert res.passes_gap is False
    assert res.passes_abs is False
    assert evaluate_passes(res, cfg) is False


def test_compute_potency_empty_marker_scores():
    cfg = KADPConfig()
    res = compute_potency({}, ({"RPC", "NRPC"}, {"Muller_Glia"}), cfg)
    assert res.max_prog == 0.0
    assert res.max_term == 0.0
    assert res.best_progenitor_type is None
    assert evaluate_passes(res, cfg) is False


def test_compute_potency_best_progenitor_argmax():
    cfg = KADPConfig()
    res = compute_potency(
        {"RPC": 0.4, "NRPC": 0.6, "Photoreceptor_Precursor": 0.2, "Rod_Photoreceptor": 0.1},
        ({"RPC", "NRPC", "Photoreceptor_Precursor"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert res.best_progenitor_type == "NRPC"
    # Tie-break is deterministic: same score → lexicographically larger name
    # ("NRPC" < "RPC").
    res_tie = compute_potency(
        {"RPC": 0.5, "NRPC": 0.5, "Rod_Photoreceptor": 0.1},
        ({"RPC", "NRPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert res_tie.best_progenitor_type == "RPC"


def test_compute_potency_ratio_boundary():
    cfg = KADPConfig()  # ratio_threshold = 2.0
    at = compute_potency(
        {"RPC": 0.4, "Rod_Photoreceptor": 0.2},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert at.passes_ratio is True  # exactly 2.0
    assert at.passes_abs is False  # 0.4 < 0.6 → isolated
    below = compute_potency(
        {"RPC": 0.39, "Rod_Photoreceptor": 0.2},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert below.passes_ratio is False  # 1.95 < 2.0


def test_compute_potency_gap_boundary():
    # 0.25 is exactly representable in binary → exact-boundary check is stable.
    cfg = KADPConfig(gap_threshold=0.25)
    at = compute_potency(
        {"RPC": 0.5, "Rod_Photoreceptor": 0.25},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert at.gap == 0.25
    assert at.passes_gap is True  # gap == threshold → passes
    assert at.passes_abs is False  # 0.5 < 0.6 → isolated
    below = compute_potency(
        {"RPC": 0.45, "Rod_Photoreceptor": 0.25},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert below.gap == 0.2
    assert below.passes_gap is False  # 0.2 < 0.25


def test_compute_potency_abs_boundary():
    # abs threshold isolated by raising the other two thresholds.
    cfg = KADPConfig(ratio_threshold=10.0, gap_threshold=10.0)
    at = compute_potency(
        {"RPC": 0.6, "Rod_Photoreceptor": 0.2},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert at.passes_abs is True  # 0.6 >= 0.6 and 0.6 > 0.2
    assert at.passes_ratio is False
    assert at.passes_gap is False
    below = compute_potency(
        {"RPC": 0.59, "Rod_Photoreceptor": 0.2},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert below.passes_abs is False


def test_compute_potency_no_naming_when_zero():
    """Failure scenario: all progenitor scores zero → no naming, all variants
    fail regardless of terminal scores."""
    cfg = KADPConfig()
    res = compute_potency(
        {"RPC": 0.0, "Rod_Photoreceptor": 0.5},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        cfg,
    )
    assert res.max_prog == 0.0
    assert res.best_progenitor_type is None
    assert res.passes_ratio is False
    assert res.passes_gap is False
    assert res.passes_abs is False
    assert evaluate_passes(res, cfg) is False


# ── evaluate_passes (hard-coded combination) ───────────────────────────


def test_evaluate_passes_combination():
    # use_gap_criterion=False (default): gap alone is insufficient.
    gap_only = PotencyResult(0.3, 0.2, 1.5, 0.1, "RPC", False, False, True)
    assert evaluate_passes(gap_only, KADPConfig()) is False
    # use_gap_criterion=True: gap-only passes.
    assert evaluate_passes(gap_only, KADPConfig(use_gap_criterion=True)) is True
    # ratio-only passes regardless of use_gap_criterion.
    ratio_only = PotencyResult(0.8, 0.2, 4.0, 0.6, "RPC", True, False, False)
    assert evaluate_passes(ratio_only, KADPConfig()) is True
    assert evaluate_passes(ratio_only, KADPConfig(use_gap_criterion=True)) is True
    # abs-only passes regardless of use_gap_criterion.
    abs_only = PotencyResult(0.7, 0.3, 7 / 3, 0.4, "RPC", False, True, False)
    assert evaluate_passes(abs_only, KADPConfig()) is True
    # all-false → False.
    all_false = PotencyResult(0.2, 0.8, 0.25, -0.6, None, False, False, False)
    assert evaluate_passes(all_false, KADPConfig(use_gap_criterion=True)) is False


# ── to_potency_dict ────────────────────────────────────────────────────


def test_to_potency_dict_roundtrip():
    res = compute_potency(
        {"RPC": 0.9, "Rod_Photoreceptor": 0.3},
        ({"RPC"}, {"Rod_Photoreceptor"}),
        KADPConfig(),
    )
    d = res.to_potency_dict()
    assert d == {"ratio": pytest.approx(3.0), "abs": pytest.approx(0.9), "gap": pytest.approx(0.6)}
    # Three values are preserved — never a single float.
    assert set(d) == {"ratio", "abs", "gap"}
    # Module-level convenience mirrors the method.
    assert to_potency_dict(res) == d


def test_to_potency_dict_zero_case():
    res = compute_potency({}, (set(), set()), KADPConfig())
    assert res.to_potency_dict() == {"ratio": 0.0, "abs": 0.0, "gap": 0.0}
