"""Tests for rna/annotation_engine.py — use_raw fix in rank_genes_groups.

T1 (P0-CRITICAL) from cross-batch-critical-fixes plan:
  rank_genes_groups must pass use_raw=True when adata.raw exists,
  with a null-guard for adata without .raw.
"""

import copy as _copy_mod
import json
import logging
import os
import sys
import tempfile
import types as _types_mod
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

import rna.utils.evidence_fusion as _ef
from core.annotation.engine import (
    _apply_canonical_expression_fallback,
    _cluster_marker_pcts,
    _flag_ai_only_decisions,
    _map_cell_state,
    _top_consensus_markers,
    _write_quality_report,
    run_unified_annotation,
)
from core.annotation.scoring import Score
from rna.utils.evidence_fusion import DiagnosticInfo, FusionDecision

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


def test_T1_use_raw_uses_raw_layer() -> None:
    """Happy path: use_raw=True surfaces genes only present in .raw.

    Simulates the scenario rank_genes_groups runs in ``run_unified_annotation``
    after the fix: .X has 100 HVG-subset genes, .raw has 1000 full genes.
    KB marker "FOO" exists only in .raw.  With ``use_raw=True``,
    rank_genes_groups must find FOO in DE results.
    """
    n_cells = 30
    n_hvg = 100
    n_full = 1000
    rng = np.random.RandomState(42)

    # .X: 100 HVG-subset genes (no "FOO")
    x = rng.poisson(lam=1.0, size=(n_cells, n_hvg)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"HVG_{i}" for i in range(n_hvg)]

    # .raw: 1000 full genes, last one is "FOO"
    raw_x = rng.poisson(lam=1.0, size=(n_cells, n_full)).astype(np.float32)
    raw = AnnData(raw_x)
    raw_var_names = [f"GENE_{i}" for i in range(n_full)]
    raw_var_names[-1] = "FOO"
    raw.var_names = raw_var_names
    adata.raw = raw

    # Add leiden clusters
    leiden = rng.choice(["0", "1"], n_cells)
    adata.obs["leiden"] = leiden

    # Make FOO strongly differentially expressed in cluster 0 vs 1
    cluster_0 = leiden == "0"
    cluster_1 = leiden == "1"
    adata.raw.X[cluster_0, -1] = 100.0  # FOO is the last column
    adata.raw.X[cluster_1, -1] = 0.0

    # ── The fix: use_raw=True when .raw exists ─────────────────────────
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        use_raw=True,
    )

    # Assert FOO (only in .raw) appears in DE results for cluster 0
    df = sc.get.rank_genes_groups_df(adata, group="0")
    top_genes = df["names"].tolist()
    assert "FOO" in top_genes, (
        f"FOO (only in .raw) must appear in DE results with use_raw=True. "
        f"Top-20 genes: {top_genes[:20]}"
    )


def test_T1_use_raw_null_guard() -> None:
    """Failure path: null guard prevents crash when adata.raw is None.

    Without the null guard (bare ``use_raw=True``), scanpy raises a KeyError
    on ``adata.raw``.  With the conditional ``use_raw=True if adata.raw
    is not None else None``, it silently falls back to .X.
    """
    n_cells = 30
    n_genes = 100
    rng = np.random.RandomState(42)

    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs["leiden"] = rng.choice(["0", "1"], n_cells)

    # .raw is NOT set — null-guard scenario
    assert adata.raw is None

    # This must not raise: the conditional skips use_raw when .raw is None
    use_raw = True if adata.raw is not None else None
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        use_raw=use_raw,
    )

    # Verify results exist (used .X, no crash)
    df = sc.get.rank_genes_groups_df(adata, group="0")
    assert len(df) > 0, "DE results should be present when falling back to .X"


# ═══════════════════════════════════════════════════════════════════════
#  T4 — Case-insensitive gene-name matching and zero-score warning
# ═══════════════════════════════════════════════════════════════════════
def _make_zero_scores(n_clusters: int = 5) -> dict:
    """Build all_scores dict where every cluster has zero KB hits."""
    return {
        str(i): {
            "CT": Score(
                score=0.0, p_value=1.0, method="none", n_markers_found=0, negative_penalty=False
            )
        }
        for i in range(n_clusters)
    }


def _make_kb_lowercase() -> dict:
    """Build a minimal KB with lowercase marker keys (case-mismatch scenario)."""
    return {
        "CT": {
            "markers": {
                "confirm": {"rho": ["PMID1"]},
                "add": {"gnat1": ["PMID2"]},
            },
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
    }


def _make_marker_df(n_clusters: int = 5) -> pd.DataFrame:
    """Build marker_df with all-uppercase DE gene names."""
    rows = []
    for cl in range(n_clusters):
        genes = ["RHO", "GNAT1", "GENE3", "GENE4", "GENE5"]
        for i, g in enumerate(genes):
            rows.append(
                {
                    "names": g,
                    "logfoldchanges": 5.0 - i * 0.5,
                    "pvals_adj": 1e-50,
                    "cluster": str(cl),
                }
            )
    return pd.DataFrame(rows)


def _make_logger() -> logging.Logger:
    """Create a logger for _check_zero_scores_and_retry."""
    log = logging.getLogger("test_T4")
    log.setLevel(logging.DEBUG)
    log.addHandler(logging.NullHandler())
    return log


def _check_zero_scores_and_retry_wrapper(
    kb,
    all_scores,
    marker_df,
    clusters,
    species,
    target_class,
    target_order,
    tissue_kb,
    logger,
    cfg=None,
):
    """Lazy-import and call _check_zero_scores_and_retry.

    Avoids circular-import issues at module level by importing only when called.
    """
    from core.annotation.engine import _check_zero_scores_and_retry

    if cfg is None:
        # Minimal stub: production reads CFG.marker.candidate_pool_expand_steps
        # (engine.py:132) to size the case-insensitive retry pass.
        cfg = SimpleNamespace(marker=SimpleNamespace(candidate_pool_expand_steps=[1]))

    return _check_zero_scores_and_retry(
        kb,
        all_scores,
        marker_df,
        clusters,
        species,
        target_class,
        target_order,
        tissue_kb,
        cfg,
        logger,
    )


def test_T4_case_insensitive_retry_succeeds() -> None:
    """Lowercase DE genes + uppercase KB -> retry fixes zero scores.

    all_scores starts with zero KB hits; after _check_zero_scores_and_retry
    uppercases KB keys and re-runs scoring, the retried scores have hits.
    """
    all_scores = _make_zero_scores(5)
    kb = _make_kb_lowercase()
    marker_df = _make_marker_df(5)
    clusters = [str(i) for i in range(5)]
    logger = _make_logger()

    # Mock score_cluster_against_kb to return hits on retry
    with patch(
        "core.annotation.scoring.score_cluster_against_kb",
        return_value={"CT": Score(0.85, 0.001, "hypergeometric", 2, False)},
    ):
        result_scores, total_hits, n_clusters = _check_zero_scores_and_retry_wrapper(
            kb,
            all_scores,
            marker_df,
            clusters,
            species="human",
            target_class="",
            target_order="",
            tissue_kb="test_kb",
            logger=logger,
        )

    assert total_hits > 0, f"Expected retry to improve hits, got total_hits={total_hits}"
    assert n_clusters == 5
    # The original all_scores must be unchanged
    for v in all_scores.values():
        assert list(v.values())[0].n_markers_found == 0, "Original all_scores mutated"


def test_T4_case_insensitive_skip_when_already_matching() -> None:
    """Already matching scores — retry should not trigger.

    When total_hits > 0, _check_zero_scores_and_retry must skip the
    retry block entirely and return the original all_scores unchanged.
    """
    all_scores = {
        str(i): {
            "CT": Score(
                score=0.5,
                p_value=0.01,
                method="hypergeometric",
                n_markers_found=1,
                negative_penalty=False,
            )
        }
        for i in range(5)
    }
    kb = _make_kb_lowercase()
    marker_df = _make_marker_df(5)
    clusters = [str(i) for i in range(5)]
    logger = _make_logger()

    # Patch score_cluster_against_kb to track if it gets called
    with patch("core.annotation.scoring.score_cluster_against_kb") as mock_sc:
        result_scores, total_hits, n_clusters = _check_zero_scores_and_retry_wrapper(
            kb,
            all_scores,
            marker_df,
            clusters,
            species="human",
            target_class="",
            target_order="",
            tissue_kb="test_kb",
            logger=logger,
        )

    mock_sc.assert_not_called()
    assert total_hits == 5, f"Expected total_hits=5 (1 per cluster), got {total_hits}"
    assert n_clusters == 5


def test_T4_zero_score_warning_fires() -> None:
    """Species mismatch — zero-score ERROR diagnostic fires.

    When retry does not improve total_hits, the function must log.error
    with a diagnostic hint about species coverage.
    """
    all_scores = _make_zero_scores(5)
    kb = _make_kb_lowercase()
    marker_df = _make_marker_df(5)
    clusters = [str(i) for i in range(5)]

    # Mock score_cluster_against_kb to return ZERO hits (retry doesn't help)
    with patch(
        "core.annotation.scoring.score_cluster_against_kb",
        return_value={"CT": Score(0.0, 1.0, "none", 0, False)},
    ):
        logger = logging.getLogger("test_T4_warning")
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.NullHandler())
        with patch.object(logger, "error") as mock_err:
            _check_zero_scores_and_retry_wrapper(
                kb,
                all_scores,
                marker_df,
                clusters,
                species="danio_rerio",
                target_class="",
                target_order="",
                tissue_kb="test_kb",
                logger=logger,
            )

    mock_err.assert_called_once()
    msg = mock_err.call_args[0][0]
    assert "species mismatch" in msg or "missing cell types" in msg, (
        f"Diagnostic hint missing in ERROR message: {msg}"
    )


def test_T4_kb_dict_not_mutated() -> None:
    """Deep copy prevents side effects — original KB unchanged after retry.

    Verifies that _check_zero_scores_and_retry does not mutate the original
    KB dict (lowercase marker keys should remain lowercase).
    """
    all_scores = _make_zero_scores(5)
    kb = _make_kb_lowercase()
    kb_original = _copy_mod.deepcopy(kb)
    marker_df = _make_marker_df(5)
    clusters = [str(i) for i in range(5)]
    logger = _make_logger()

    with patch(
        "core.annotation.scoring.score_cluster_against_kb",
        return_value={"CT": Score(0.85, 0.001, "hypergeometric", 2, False)},
    ):
        _check_zero_scores_and_retry_wrapper(
            kb,
            all_scores,
            marker_df,
            clusters,
            species="human",
            target_class="",
            target_order="",
            tissue_kb="test_kb",
            logger=logger,
        )

    # Original KB must be identical
    assert kb == kb_original, "Original KB dict was mutated"
    # Specifically check lowercase keys survived
    assert "rho" in kb["CT"]["markers"]["confirm"], (
        "KB marker key 'rho' became uppercase — original KB mutated"
    )


# ═══════════════════════════════════════════════════════════════════════

#  T9 — transition-annotation-p0: cell_state map (D5) + quality report

#        new fields (D6 review_queue / D7 transition_clusters / D8 kb_coverage)

# ═══════════════════════════════════════════════════════════════════════


def _make_decision(
    method="marker_scoring", cell_type="Rod Photoreceptor", confidence="high", diagnostic=None
):
    """Minimal decision stub.



    ``_map_cell_state`` only reads method/cell_type/confidence;

    ``_write_quality_report`` additionally reads ``.diagnostic`` for the

    ambiguous review_queue row.  SimpleNamespace suffices — the functions

    never touch the other FusionDecision fields.

    """

    return SimpleNamespace(
        method=method, cell_type=cell_type, confidence=confidence, diagnostic=diagnostic
    )


def test_T9_map_cell_state_transition_state() -> None:
    """D5 row 1: method == 'transition_state' → 'transient_transitional'."""

    d = _make_decision(method="transition_state", cell_type="transitional: RGC/Amacrine")

    assert _map_cell_state(d, "Broad_Neuron") == "transient_transitional"


def test_T9_map_cell_state_ambiguous() -> None:
    """D5 row 2: method == 'ambiguous' → 'N/A' (downgraded, awaiting review)."""

    d = _make_decision(method="ambiguous", confidence="unknown")

    assert _map_cell_state(d, "Broad_Neuron") == "N/A"


def test_T9_map_cell_state_cycling() -> None:
    """D5 row 3: cell_type contains 'Proliferating' → 'cycling'."""

    d = _make_decision(cell_type="Proliferating_RPC")

    assert _map_cell_state(d, "Broad_Progenitor") == "cycling"


def test_T9_map_cell_state_committed_precursor() -> None:
    """D5 row 4: cell_category == 'Broad_Progenitor' (non-proliferating) → 'committed_precursor'."""

    d = _make_decision(cell_type="RPC")

    assert _map_cell_state(d, "Broad_Progenitor") == "committed_precursor"


def test_T9_map_cell_state_terminal_high_neuron() -> None:
    """D5 row 5a: confidence high + Broad_Neuron → 'terminal'."""

    d = _make_decision(confidence="high", cell_type="RGC")

    assert _map_cell_state(d, "Broad_Neuron") == "terminal"


def test_T9_map_cell_state_terminal_medium_glia() -> None:
    """D5 row 5b: confidence medium + Broad_Glia → 'terminal'."""

    d = _make_decision(confidence="medium", cell_type="Müller Glia")

    assert _map_cell_state(d, "Broad_Glia") == "terminal"


def test_T9_map_cell_state_terminal_high_non_neural() -> None:
    """D5 row 5c: confidence high + Broad_Non-neural → 'terminal'."""

    d = _make_decision(confidence="high", cell_type="Endothelial")

    assert _map_cell_state(d, "Broad_Non-neural") == "terminal"


def test_T9_map_cell_state_low_na() -> None:
    """D5 row 6: confidence 'low' → 'N/A' even under a terminal broad category."""

    d = _make_decision(confidence="low", cell_type="RGC")

    assert _map_cell_state(d, "Broad_Neuron") == "N/A"


def test_T9_map_cell_state_unknown_na() -> None:
    """D5 fallback row: unknown/other confidence → 'N/A'."""

    d = _make_decision(confidence="unknown", cell_type="RGC")

    assert _map_cell_state(d, "Broad_Neuron") == "N/A"


def test_T9_map_cell_state_proliferating_precedes_progenitor() -> None:
    """D5 order sensitivity: Proliferating row wins before Broad_Progenitor row."""

    d = _make_decision(cell_type="Proliferating_MG")

    assert _map_cell_state(d, "Broad_Progenitor") == "cycling"


def test_T9_map_cell_state_developmental_potency() -> None:
    """KADP: method == 'developmental_potency' → 'differentiating' (plan todo 5)."""

    d = _make_decision(method="developmental_potency", cell_type="NRPC", confidence="medium")

    assert _map_cell_state(d, "Broad_Progenitor") == "differentiating"


def test_T9_map_cell_state_kadp_precedes_proliferating_row() -> None:
    """Row order: the differentiating row sits BEFORE the 'Proliferating' row.

    A KADP-named type that happens to contain 'Proliferating' (e.g.
    Proliferating_RPC) must still map to 'differentiating', not 'cycling'."""

    d = _make_decision(
        method="developmental_potency", cell_type="Proliferating_RPC", confidence="medium"
    )

    assert _map_cell_state(d, "Broad_Progenitor") == "differentiating"


def _make_kb() -> dict:
    """Minimal retina-like KB: 4 fine types + hierarchy/expert/Broad_* containers."""

    return {
        "RGC": {
            "markers": {"confirm": {"RBPMS": ["PMID1"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Amacrine": {
            "markers": {"confirm": {"TFAP2A": ["PMID2"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Rod Photoreceptor": {
            "markers": {"confirm": {"RHO": ["PMID3"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Müller Glia": {
            "markers": {"confirm": {"RLBP1": ["PMID4"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "_hierarchy": {"categories": ["Progenitor", "Neuron", "Glia", "Non-neural"]},
        "expert_rules": {},
        "Broad_Neuron": {},
    }


def test_T9_quality_report_new_fields() -> None:
    """D6/D7/D8: review_queue + transition_clusters detail + kb_coverage in the JSON.



    Given: decision_map with one transition_state decision

           (cell_type='transitional: RGC/Amacrine') and one ambiguous decision

           (DiagnosticInfo category='ambiguous' with full top_competitors).

    When:  _write_quality_report writes 05_annotation_quality.json.

    Then:  transition_clusters lists {cluster, pair}; review_queue carries the

           ambiguous cluster with n_tied_types/top_types; kb_coverage carries

           annotated_types/kb_types_unannotated/ghost_endpoints.

    """

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Production calls _write_quality_report BEFORE marker_validation is set

        # (05_annotate_major.py:627), so the column is absent → pass_cells is a

        # pure Python int (0) and the JSON stays numpy-scalar-free.

        adata = SimpleNamespace(
            obs={},  # no marker_validation column — mirrors production timing
            n_obs=10,
        )

        ann_records = [
            {"cluster": "0", "reasoning": "best match", "ai_agreed": True},
            {"cluster": "1", "reasoning": "also matched rules: X", "ai_agreed": False},
            {"cluster": "2", "reasoning": "best match", "ai_agreed": True},
        ]

        fusion_quality = {"annotated_by_rule": 0, "annotated_by_scoring": 2, "unknown": 1}

        cell_category_map = {"0": "Broad_Neuron", "1": "", "2": "Broad_Neuron"}

        decision_map = {
            "0": _make_decision(
                method="transition_state",
                cell_type="transitional: RGC/Amacrine",
                confidence="transition",
            ),
            "1": _make_decision(
                method="ambiguous",
                cell_type="Unknown",
                confidence="unknown",
                diagnostic=DiagnosticInfo(
                    category="ambiguous",
                    top_competitors=[
                        {"cell_type": "RGC", "score": 1.0},
                        {"cell_type": "Amacrine", "score": 0.98},
                        {"cell_type": "Cone", "score": 0.95},
                    ],
                    detail="Multi-peak tie: 3 types >= 0.9 (top: RGC=1.000, Amacrine=0.980, Cone=0.950)",
                ),
            ),
            "2": _make_decision(method="marker_scoring", cell_type="RGC", confidence="high"),
        }

        cfg = SimpleNamespace(table_dir=tmp_dir)

        logger = _make_logger()

        _write_quality_report(
            adata,
            ann_records,
            fusion_quality,
            cell_category_map,
            decision_map,
            cfg,
            logger,
            kb=_make_kb(),
        )

        quality_path = os.path.join(tmp_dir, "05_annotation_quality.json")

        assert os.path.exists(quality_path)

        with open(quality_path) as f:
            quality = json.load(f)

        # D7 — transition_clusters is a {cluster, pair} detail list

        assert quality["transition_clusters"] == [{"cluster": "0", "pair": "RGC/Amacrine"}]

        # D6 — review_queue carries the ambiguous cluster with tie detail

        assert [q["cluster"] for q in quality["review_queue"]] == ["1"]

        entry = quality["review_queue"][0]

        # Task 10 (D6): ambiguous (multi-peak) entries carry reason "ambiguous"
        assert entry["reason"] == "ambiguous"

        assert entry["n_tied_types"] == 3

        assert entry["top_types"] == ["RGC", "Amacrine", "Cone"]

        # D8 — kb_coverage has the three keys

        cov = quality["kb_coverage"]

        assert set(cov) == {"annotated_types", "kb_types_unannotated", "ghost_endpoints"}

        assert cov["annotated_types"] == sorted({"transitional: RGC/Amacrine", "Unknown", "RGC"})

        # RGC is annotated as cluster '2', so only Amacrine is a ghost endpoint

        assert cov["ghost_endpoints"] == ["Amacrine"]

        assert cov["kb_types_unannotated"] == sorted(
            {"Amacrine", "Rod Photoreceptor", "Müller Glia"}
        )


def test_T9_quality_report_kb_none_backward_compat() -> None:
    """D8: kb=None must not raise and yields empty kb_types_unannotated."""

    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = SimpleNamespace(
            obs={},  # no marker_validation column — mirrors production timing
            n_obs=10,
        )

        decision_map = {
            "0": _make_decision(
                method="transition_state",
                cell_type="transitional: RGC/Amacrine",
                confidence="transition",
            )
        }

        cfg = SimpleNamespace(table_dir=tmp_dir)

        logger = _make_logger()

        # must not raise when kb is omitted

        _write_quality_report(adata, [], {}, {}, decision_map, cfg, logger)

        with open(os.path.join(tmp_dir, "05_annotation_quality.json")) as f:
            quality = json.load(f)

        assert quality["kb_coverage"]["kb_types_unannotated"] == []

        assert quality["kb_coverage"]["ghost_endpoints"] == ["Amacrine", "RGC"]


# ═══════════════════════════════════════════════════════════════════════
#  T10 — D3 canonical-expression fallback + ai_only audit
# ═══════════════════════════════════════════════════════════════════════


def _make_fusion_decision(**overrides):
    """Build a FusionDecision stub for the D3 audit functions."""
    base = dict(
        cell_type="RGC",
        confidence="high",
        score=0.85,
        method="marker_scoring_high",
        n_markers_found=3,
        ai_agreed=False,
        ai_suggested="",
        explanation="test",
        alternative_rules=[],
        diagnostic=None,
    )
    base.update(overrides)
    return FusionDecision(**base)


def _make_consensus_kb() -> dict:
    """Minimal KB: one multi-source type (RGC) + one single-source type (RGC_Alpha)."""
    return {
        "RGC": {
            "markers": {
                "confirm": {
                    "RBPMS": ["s1", "s2", "s3"],
                    "NEFL": ["s1", "s2"],
                    "POU4F1": ["s1", "s2", "s3", "s4"],
                }
            },
            "consensus_levels": {"RBPMS": "gold", "NEFL": "high", "POU4F1": "gold"},
        },
        "RGC_Alpha": {
            "markers": {"confirm": {"SPP1": ["tran2019"]}},
            "consensus_levels": {"SPP1": "low"},
        },
    }


def _make_expression_adata(cl_pcts: dict) -> AnnData:
    """AnnData where cluster '0' (10 cells) expresses genes at given pcts.

    ``cl_pcts: {gene: fraction}`` of the 10 cluster-0 cells carry raw count 1
    for that gene; cluster '1' (5 cells) never expresses anything.
    """
    genes = list(cl_pcts)
    n_cluster, n_other = 10, 5
    x = np.zeros((n_cluster + n_other, len(genes)), dtype=np.float32)
    for j, g in enumerate(genes):
        n_expr = int(round(cl_pcts[g] * n_cluster))
        x[:n_expr, j] = 1.0
    adata = AnnData(x)
    adata.var_names = genes
    adata.obs["leiden"] = ["0"] * n_cluster + ["1"] * n_other
    raw = AnnData(x.copy())
    raw.var_names = genes
    adata.raw = raw
    return adata


def test_D3_top_consensus_markers_orders_and_filters() -> None:
    """gold/high markers qualify, ordered by rank; single-source type → empty."""
    kb = _make_consensus_kb()
    assert _top_consensus_markers(kb, "RGC") == ["RBPMS", "POU4F1", "NEFL"]
    assert _top_consensus_markers(kb, "RGC_Alpha") == []


def test_D3_canonical_fallback_downgrades_confident_label() -> None:
    """All top-consensus markers ~0 pct → confidence forced low + review reason."""
    kb = _make_consensus_kb()
    adata = _make_expression_adata({"RBPMS": 0.0, "NEFL": 0.0, "POU4F1": 0.0})
    dm = {"0": _make_fusion_decision()}
    cfg = SimpleNamespace(annotation=SimpleNamespace(canonical_pct_floor=0.05))
    reasons = _apply_canonical_expression_fallback(adata, kb, dm, cfg, _make_logger())
    assert reasons == {"0": "no_canonical_expression"}
    assert dm["0"].confidence == "low"
    assert dm["0"].review_reason == "no_canonical_expression"
    assert "no_canonical_expression" in dm["0"].explanation


def test_D3_canonical_fallback_keeps_when_markers_expressed() -> None:
    """Markers expressed above the floor → decision untouched."""
    kb = _make_consensus_kb()
    adata = _make_expression_adata({"RBPMS": 0.9, "NEFL": 0.8, "POU4F1": 0.7})
    dm = {"0": _make_fusion_decision()}
    cfg = SimpleNamespace(annotation=SimpleNamespace(canonical_pct_floor=0.05))
    reasons = _apply_canonical_expression_fallback(adata, kb, dm, cfg, _make_logger())
    assert reasons == {}
    assert dm["0"].confidence == "high"
    assert dm["0"].review_reason == ""


def test_D3_canonical_fallback_skips_when_raw_missing() -> None:
    """adata.raw is None → silent skip (mirrors _ribo_fallback_pct_scores)."""
    kb = _make_consensus_kb()
    adata = _make_expression_adata({"RBPMS": 0.0, "NEFL": 0.0, "POU4F1": 0.0})
    adata.raw = None
    dm = {"0": _make_fusion_decision()}
    cfg = SimpleNamespace(annotation=SimpleNamespace(canonical_pct_floor=0.05))
    reasons = _apply_canonical_expression_fallback(adata, kb, dm, cfg, _make_logger())
    assert reasons == {}
    assert dm["0"].confidence == "high"


def test_D3_canonical_fallback_skips_single_source_type() -> None:
    """Winning type with no consensus>=2 markers → empty set, naturally skipped."""
    kb = _make_consensus_kb()
    adata = _make_expression_adata({"SPP1": 0.0})
    dm = {"0": _make_fusion_decision(cell_type="RGC_Alpha")}
    cfg = SimpleNamespace(annotation=SimpleNamespace(canonical_pct_floor=0.05))
    reasons = _apply_canonical_expression_fallback(adata, kb, dm, cfg, _make_logger())
    assert reasons == {}
    assert dm["0"].confidence == "high"


def test_D3_canonical_fallback_default_floor_when_cfg_absent() -> None:
    """Missing CFG.annotation.canonical_pct_floor → 0.05 getattr fallback."""
    kb = _make_consensus_kb()
    adata = _make_expression_adata({"RBPMS": 0.0, "NEFL": 0.0, "POU4F1": 0.0})
    dm = {"0": _make_fusion_decision()}
    reasons = _apply_canonical_expression_fallback(
        adata, kb, dm, SimpleNamespace(), _make_logger()
    )
    assert reasons == {"0": "no_canonical_expression"}
    assert dm["0"].confidence == "low"


def test_D3_cluster_marker_pcts_csr_pattern() -> None:
    """Pct fractions computed per gene via raw counts > 0 (CSR row-slice pattern)."""
    adata = _make_expression_adata({"RBPMS": 0.5, "NEFL": 0.0, "POU4F1": 1.0})
    pcts = _cluster_marker_pcts(adata, "0", ["RBPMS", "NEFL", "POU4F1"])
    assert pcts is not None
    np.testing.assert_allclose(pcts, [0.5, 0.0, 1.0], atol=1e-6)
    # raw missing → None (silent skip contract)
    adata.raw = None
    assert _cluster_marker_pcts(adata, "0", ["RBPMS"]) is None


def test_D3_ai_only_forces_low_confidence() -> None:
    """ai_unconstrained (medium, no KB support) → forced low + ai_only reason."""
    dm = {
        "0": _make_fusion_decision(
            cell_type="Foo",
            method="ai_unconstrained",
            confidence="medium",
            score=0.0,
            n_markers_found=0,
            ai_agreed=True,
            ai_suggested="Foo",
        )
    }
    reasons = _flag_ai_only_decisions(dm, _make_logger())
    assert reasons == {"0": "ai_only"}
    assert dm["0"].confidence == "low"
    assert dm["0"].review_reason == "ai_only"


def test_D3_ai_agreed_with_kb_support_not_ai_only() -> None:
    """AI agreeing on a KB-backed decision (score/markers present) is not ai_only."""
    dm = {"0": _make_fusion_decision(ai_agreed=True, ai_suggested="RGC")}
    reasons = _flag_ai_only_decisions(dm, _make_logger())
    assert reasons == {}
    assert dm["0"].confidence == "high"
    assert dm["0"].review_reason == ""


def test_T9_quality_report_review_queue_includes_review_reason() -> None:
    """D3-downgraded decision (review_reason set) enters the quality-report review_queue."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = SimpleNamespace(obs={}, n_obs=10)
        ann_records = [
            {"cluster": "0", "reasoning": "best match", "ai_agreed": True},
        ]
        fusion_quality = {"annotated_by_rule": 0, "annotated_by_scoring": 1, "unknown": 0}
        cell_category_map = {"0": "Broad_Neuron"}
        decision_map = {
            "0": _make_fusion_decision(confidence="low", review_reason="no_canonical_expression"),
        }
        cfg = SimpleNamespace(table_dir=tmp_dir)
        _write_quality_report(
            adata,
            ann_records,
            fusion_quality,
            cell_category_map,
            decision_map,
            cfg,
            _make_logger(),
            kb=_make_consensus_kb(),
        )
        with open(os.path.join(tmp_dir, "05_annotation_quality.json")) as f:
            quality = json.load(f)
        assert [q["cluster"] for q in quality["review_queue"]] == ["0"]
        entry = quality["review_queue"][0]
        # Task 10 (D6): downgraded entries carry the decision's review_reason;
        # non-multi-peak entries have empty tie detail.
        assert entry["reason"] == "no_canonical_expression"
        assert entry["n_tied_types"] == 0
        assert entry["top_types"] == []


def test_T10_review_queue_review_reason_beats_ambiguous() -> None:
    """D6 merge rule: a decision with both review_reason and ambiguous keeps review_reason."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        adata = SimpleNamespace(obs={}, n_obs=10)
        ann_records = [{"cluster": "0", "reasoning": "best match", "ai_agreed": True}]
        fusion_quality = {"annotated_by_rule": 0, "annotated_by_scoring": 1, "unknown": 0}
        cell_category_map = {"0": "Broad_Neuron"}
        decision_map = {
            "0": _make_fusion_decision(
                method="ambiguous",
                confidence="unknown",
                review_reason="zero_evidence",
                diagnostic=DiagnosticInfo(
                    category="ambiguous",
                    top_competitors=[{"cell_type": "RGC", "score": 1.0}],
                    detail="tie",
                ),
            ),
        }
        cfg = SimpleNamespace(table_dir=tmp_dir)
        _write_quality_report(
            adata,
            ann_records,
            fusion_quality,
            cell_category_map,
            decision_map,
            cfg,
            _make_logger(),
            kb=_make_consensus_kb(),
        )
        with open(os.path.join(tmp_dir, "05_annotation_quality.json")) as f:
            quality = json.load(f)
        entry = quality["review_queue"][0]
        assert entry["reason"] == "zero_evidence"


def test_T10_review_queue_reason_passthrough() -> None:
    """D6: downgraded decisions surface their review_reason verbatim in review_queue.

    Covers every non-ambiguous reason in the D6 enum, plus the old-format
    entry (no reason) staying loadable.
    """
    reasons = [
        "single_marker_rule",
        "window_padding",
        "weak_multi",
        "zero_evidence",
        "ai_only",
        "no_canonical_expression",
    ]
    for reason in reasons:
        with tempfile.TemporaryDirectory() as tmp_dir:
            adata = SimpleNamespace(obs={}, n_obs=10)
            ann_records = [{"cluster": "0", "reasoning": "best match", "ai_agreed": True}]
            fusion_quality = {"annotated_by_rule": 0, "annotated_by_scoring": 1, "unknown": 0}
            cell_category_map = {"0": "Broad_Neuron"}
            decision_map = {
                "0": _make_fusion_decision(confidence="low", review_reason=reason),
            }
            cfg = SimpleNamespace(table_dir=tmp_dir)
            _write_quality_report(
                adata,
                ann_records,
                fusion_quality,
                cell_category_map,
                decision_map,
                cfg,
                _make_logger(),
                kb=_make_consensus_kb(),
            )
            with open(os.path.join(tmp_dir, "05_annotation_quality.json")) as f:
                quality = json.load(f)
            entry = quality["review_queue"][0]
            assert entry["cluster"] == "0"
            assert entry["reason"] == reason
            # non-multi-peak: no tie detail
            assert entry["n_tied_types"] == 0
            assert entry["top_types"] == []


# ═══════════════════════════════════════════════════════════════════════
#  F1 — AI-fallback gate: ambiguous + transition_state inclusion
#  (annotation-kadp-metc todo 2).  The two-segment ``low_conf_clusters``
#  selection must pull ambiguous / transition_state candidates (incl.
#  confidence="transition") into the AI fallback only when kadp/metc are
#  enabled; with both off it stays byte-identical to the baseline
#  ``confidence in ("low","unknown") and method != "ambiguous"`` rule.
# ═══════════════════════════════════════════════════════════════════════


def _make_ai_gate_adata() -> AnnData:
    """3 leiden clusters with cluster-distinct expression (real rank_genes_groups)."""
    rng = np.random.RandomState(0)
    var_names = ["RBPMS", "POU4F1", "ONECUT1", "GENE_01", "GENE_02", "GENE_03"]
    n_cells = 90
    x = rng.poisson(lam=1.0, size=(n_cells, len(var_names))).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = var_names
    adata.obs["leiden"] = np.repeat(["0", "1", "2"], 30)
    c0 = adata.obs["leiden"] == "0"
    c1 = adata.obs["leiden"] == "1"
    adata.X[c0, 0] = adata.X[c0, 0] + rng.poisson(lam=20.0, size=int(c0.sum())).astype(np.float32)
    adata.X[c0, 1] = adata.X[c0, 1] + rng.poisson(lam=20.0, size=int(c0.sum())).astype(np.float32)
    adata.X[c1, 2] = adata.X[c1, 2] + rng.poisson(lam=20.0, size=int(c1.sum())).astype(np.float32)
    adata.obsm["X_umap"] = rng.randn(n_cells, 2).astype(np.float32)
    return adata


def _make_ai_gate_cfg(
    tmp_path, kadp_enabled: bool = False, metc_enabled: bool = False
) -> SimpleNamespace:
    """Plain-value config: AI on, celltypist off, optional kadp/metc flags."""
    table_dir = tmp_path / "results" / "tables"
    figure_dir = tmp_path / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    return SimpleNamespace(
        tissue_kb="retina",
        tissue="retina",
        species="human",
        target_class="",
        target_order="",
        tissue_maturity="",
        table_dir=str(table_dir),
        figure_dir=str(figure_dir),
        interactive=False,
        marker=SimpleNamespace(
            candidate_pool_expand_steps=[5],
            expert_rule_top_n=0,
            expert_rule_pval_cutoff=0.0,
            expert_rule_strictness="default",
            developmental_mode=False,
        ),
        annotation=SimpleNamespace(
            method="kb_unified",
            celltypist=SimpleNamespace(enabled=False, model="", majority_voting=False),
            kadp_enabled=kadp_enabled,
            metc_enabled=metc_enabled,
        ),
        ai=SimpleNamespace(enabled=True, ai_annotation=True, unconstrained_annotation=False),
        plot=SimpleNamespace(figure_dpi=150, figure_format="pdf", figure_transparent=True),
        execution=SimpleNamespace(random_seed=42),
    )


def _make_ai_gate_kb() -> dict:
    """Minimal retina-like KB without ``_hierarchy`` (skips the tiered block)."""
    return {
        "RGC": {
            "markers": {"confirm": {"RBPMS": ["PMID1"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Amacrine": {
            "markers": {"confirm": {"TFAP2A": ["PMID2"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "expert_rules": {},
        "Broad_Neuron": {},
    }


def _ai_gate_canned_scores() -> dict:
    """Non-zero canned scores: retry path never fires, deterministic fusion input."""
    return {
        "RGC": Score(
            score=0.60,
            p_value=0.05,
            method="test",
            n_markers_found=2,
            negative_penalty=False,
            tier="L2",
            private_markers_hit=0,
            consensus="high",
            n_sources=1,
        ),
        "Amacrine": Score(
            score=0.55,
            p_value=0.05,
            method="test",
            n_markers_found=1,
            negative_penalty=False,
            tier="L2",
            private_markers_hit=0,
            consensus="",
            n_sources=1,
        ),
    }


def _ai_gate_decisions() -> list:
    """Canned fusion output for clusters 0/1/2: ambiguous, transition_state, high."""
    ambiguous = FusionDecision(
        cell_type="Unknown",
        confidence="low",
        score=0.0,
        method="ambiguous",
        n_markers_found=0,
        ai_agreed=False,
        ai_suggested="",
        explanation="Multi-peak tie (RGC/Amacrine 0.9).",
        alternative_rules=[],
        diagnostic=DiagnosticInfo(
            category="ambiguous",
            top_competitors=[
                {"cell_type": "RGC", "score": 0.9},
                {"cell_type": "Amacrine", "score": 0.9},
            ],
            detail="multi-peak",
        ),
    )
    transition = FusionDecision(
        cell_type="transitional: RGC/Amacrine",
        confidence="transition",
        score=0.5,
        method="transition_state",
        n_markers_found=2,
        ai_agreed=False,
        ai_suggested="",
        explanation="Top-2 same lineage, delta below threshold.",
        alternative_rules=[],
    )
    high = FusionDecision(
        cell_type="RGC",
        confidence="high",
        score=0.9,
        method="marker_scoring_high",
        n_markers_found=3,
        ai_agreed=False,
        ai_suggested="",
        explanation="Strong RBPMS/POU4F1.",
        alternative_rules=[],
        cell_category="Broad_Neuron",
        tier="L2",
    )
    return [ambiguous, transition, high]


def _run_ai_gate_engine(adata, cfg, ai_response: str):
    """Drive the real ``run_unified_annotation`` with canned fusion + AI."""
    decisions = _ai_gate_decisions()
    captured_ai: dict = {}
    fuse_calls: list = []

    def _fake_fuse_all(*args, **kwargs):
        fuse_calls.append(kwargs)
        if "ai_results" in kwargs:
            captured_ai.update(kwargs["ai_results"])
        if len(fuse_calls) == 1:
            return list(decisions), {}
        return list(decisions), {}

    logger = MagicMock()
    with (
        patch("core.kb.load_kb", return_value=_make_ai_gate_kb()),
        patch(
            "core.annotation.scoring.score_cluster_against_kb",
            side_effect=lambda *a, **k: dict(_ai_gate_canned_scores()),
        ),
        patch("rna.utils.evidence_fusion.fuse_all_clusters", side_effect=_fake_fuse_all),
        patch("core.ai.caller.ai_query", return_value=ai_response) as ai_mock,
        patch("core.annotation.engine.safe_plot"),
    ):
        result = run_unified_annotation(adata, cfg, logger)
    return result, captured_ai, fuse_calls, ai_mock


def test_ai_fallback_gate_includes_ambiguous_decision(tmp_path) -> None:
    """F1: ``method='ambiguous'`` enters low_conf_clusters when kadp is on.

    Enabled → AI fallback fires and the ambiguous cluster receives an
    ``ai_results`` entry; disabled → baseline parity (excluded, no AI call).
    """
    ai_response = json.dumps(
        {
            "0": {"cell_type": "Amacrine Cell"},
            "1": {"cell_type": "RGC"},
        }
    )

    # ── enabled (kadp_enabled=True) ──
    cfg_on = _make_ai_gate_cfg(tmp_path / "on", kadp_enabled=True)
    result, captured_ai, fuse_calls, ai_mock = _run_ai_gate_engine(
        _make_ai_gate_adata(), cfg_on, ai_response
    )
    assert result is not None
    assert ai_mock.called, "ambiguous candidate must be gated into the AI fallback"
    assert len(fuse_calls) == 2, "second-pass fusion with ai_results must run"
    assert captured_ai.get("0") == "Amacrine_Cell", captured_ai
    assert captured_ai.get("1") == "RGC", captured_ai
    assert captured_ai.get("1") == "RGC", captured_ai

    # ── disabled (baseline parity) ──
    cfg_off = _make_ai_gate_cfg(tmp_path / "off")
    result_off, captured_ai_off, fuse_calls_off, ai_mock_off = _run_ai_gate_engine(
        _make_ai_gate_adata(), cfg_off, ai_response
    )
    assert result_off is not None
    assert not ai_mock_off.called, "baseline: ambiguous must stay out of low_conf_clusters"
    assert captured_ai_off == {}
    assert len(fuse_calls_off) == 1, "no second-pass fusion when fallback is skipped"


def test_ai_fallback_gate_includes_transition_state_decision(tmp_path) -> None:
    """F1: ``method='transition_state'`` (confidence='transition') enters the
    AI fallback when metc is on — the confidence value the baseline filter
    never matched.  Disabled → baseline parity (excluded, no AI call).
    """
    ai_response = json.dumps(
        {
            "0": {"cell_type": "Amacrine Cell"},
            "1": {"cell_type": "RGC"},
        }
    )

    # ── enabled (metc_enabled=True) ──
    cfg_on = _make_ai_gate_cfg(tmp_path / "on", metc_enabled=True)
    result, captured_ai, fuse_calls, ai_mock = _run_ai_gate_engine(
        _make_ai_gate_adata(), cfg_on, ai_response
    )
    assert result is not None
    assert ai_mock.called, "transition_state candidate must be gated into the AI fallback"
    assert len(fuse_calls) == 2
    assert captured_ai.get("1") == "RGC", captured_ai
    assert captured_ai.get("0") == "Amacrine_Cell", captured_ai

    # ── disabled (baseline parity) ──
    cfg_off = _make_ai_gate_cfg(tmp_path / "off")
    result_off, captured_ai_off, fuse_calls_off, ai_mock_off = _run_ai_gate_engine(
        _make_ai_gate_adata(), cfg_off, ai_response
    )
    assert result_off is not None
    assert not ai_mock_off.called, "baseline: transition_state must stay out of low_conf_clusters"
    assert captured_ai_off == {}
    assert len(fuse_calls_off) == 1


# ═══════════════════════════════════════════════════════════════════════
#  Todo 8 — AI labels share the harmonization chain before the 2nd fuse.
# ═══════════════════════════════════════════════════════════════════════


def test_ai_labels_harmonized_before_second_fuse(tmp_path) -> None:
    """AI labels pass through harmonize_label (same chain as CellTypist).

    Resolvable labels are replaced by their canonical KB name; unresolvable
    ones are dropped so that cluster's ai_suggestion abstains (Oracle r2
    MAJOR 1 — an unparseable AI vote must not inflate METC distinct)."""
    ai_response = json.dumps(
        {
            "0": {"cell_type": "Retinal Ganglion Cell"},  # → RGC (synonym hit)
            "1": {"cell_type": "NotARealCellType"},  # → dropped (abstain)
        }
    )
    cfg = _make_ai_gate_cfg(tmp_path / "ai_harm", metc_enabled=True)
    result, captured_ai, fuse_calls, ai_mock = _run_ai_gate_engine(
        _make_ai_gate_adata(), cfg, ai_response
    )
    assert result is not None
    assert ai_mock.called
    assert len(fuse_calls) == 2, "second-pass fusion must receive harmonized ai_results"
    assert captured_ai.get("0") == "RGC", captured_ai
    assert "1" not in captured_ai, "unresolvable AI label must be dropped (source abstains)"


# ═══════════════════════════════════════════════════════════════════════
#  Todo 5 — engine/schema KADP wiring: dual-call mirroring + obs-level
#  four-value assertions + tiered-block exemption + AI second-pass keep.
#  The engine drives the REAL run_unified_annotation with canned fusion
#  output; the KADPConfig is captured from BOTH fuse_all_clusters calls.
# ═══════════════════════════════════════════════════════════════════════


def _make_kadp_wiring_cfg(
    tmp_path,
    *,
    kadp_enabled: bool = True,
    ratio: float = 2.0,
    abs_: float = 0.6,
    gap: float = 0.1,
    use_gap: bool = False,
) -> SimpleNamespace:
    """Developing-tissue config: AI on, celltypist off, KADP on with
    explicit thresholds so the captured KADPConfig can be value-asserted."""
    table_dir = tmp_path / "results" / "tables"
    figure_dir = tmp_path / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    return SimpleNamespace(
        tissue_kb="retina",
        tissue="retina",
        species="human",
        target_class="",
        target_order="",
        tissue_maturity="developing",
        table_dir=str(table_dir),
        figure_dir=str(figure_dir),
        interactive=False,
        marker=SimpleNamespace(
            candidate_pool_expand_steps=[5],
            expert_rule_top_n=0,
            expert_rule_pval_cutoff=0.0,
            expert_rule_strictness="default",
            developmental_mode=False,
        ),
        annotation=SimpleNamespace(
            method="kb_unified",
            celltypist=SimpleNamespace(enabled=False, model="", majority_voting=False),
            kadp_enabled=kadp_enabled,
            kadp_ratio_threshold=ratio,
            kadp_abs_threshold=abs_,
            kadp_gap_threshold=gap,
            use_gap_criterion=use_gap,
            metc_enabled=False,
        ),
        ai=SimpleNamespace(enabled=True, ai_annotation=True, unconstrained_annotation=False),
        plot=SimpleNamespace(figure_dpi=150, figure_format="pdf", figure_transparent=True),
        execution=SimpleNamespace(random_seed=42),
    )


def _make_kadp_tiered_kb() -> dict:
    """Retina-like KB WITH ``_hierarchy`` so the engine's tiered block runs.

    The Progenitor/Neuron/Glia/Non-neural category members drive both the
    category guard and the tiered block; the KADP cluster (NRPC) must be
    exempted from the tiered block's _replace entirely."""
    return {
        "RGC": {
            "markers": {"confirm": {"RBPMS": ["PMID1"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Amacrine": {
            "markers": {"confirm": {"TFAP2A": ["PMID2"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Rod Photoreceptor": {
            "markers": {"confirm": {"RHO": ["PMID3"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "Müller Glia": {
            "markers": {"confirm": {"RLBP1": ["PMID4"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "NRPC": {
            "markers": {"confirm": {"NRPCMARK": ["PMID5"]}},
            "negative_markers": [],
            "species": ["human"],
            "synonyms": [],
        },
        "_hierarchy": {
            "categories": {
                "Progenitor": {"members": ["NRPC"], "subtypes": {}},
                "Neuron": {
                    "members": ["RGC", "Amacrine", "Rod Photoreceptor"],
                    "subtypes": {},
                },
                "Glia": {"members": ["Müller Glia"], "subtypes": {}},
                "Non-neural": {"members": [], "subtypes": {}},
            }
        },
        "expert_rules": {},
        "Broad_Progenitor": {},
        "Broad_Neuron": {},
        "Broad_Glia": {},
        "Broad_Non-neural": {},
    }


_KADP_POTENCY = {"ratio": 3.17, "abs": 0.95, "gap": 0.65}


def _kadp_wiring_decisions() -> list:
    """Canned fusion output: KADP (cluster 0), low (cluster 1), high (2).

    The KADP decision carries ``ai_agreed=True`` + ``ai_suggested='NRPC'`` to
    simulate the AI second pass agreeing with the KADP name — this makes the
    ``_clean_method`` ordering test meaningful (developmental_potency must
    win over the marker_scoring+ai suffix)."""
    kadp = FusionDecision(
        cell_type="NRPC",
        confidence="medium",
        score=0.95,
        method="developmental_potency",
        n_markers_found=3,
        ai_agreed=True,
        ai_suggested="NRPC",
        explanation="Developmental potency named 'NRPC' as differentiating precursor.",
        alternative_rules=[],
        diagnostic=DiagnosticInfo(
            category="developmental_potency",
            top_competitors=[{"cell_type": "NRPC", "score": 0.95}],
            detail="Developmental potency naming -- max progenitor 0.950 vs max terminal 0.300.",
        ),
        review_reason="kadp_precursor",
        potency=dict(_KADP_POTENCY),
        source_votes=None,
    )
    low = FusionDecision(
        cell_type="RGC",
        confidence="low",
        score=0.6,
        method="marker_scoring",
        n_markers_found=2,
        ai_agreed=False,
        ai_suggested="",
        explanation="Weak RBPMS signal.",
        alternative_rules=[],
        tier="L2",
    )
    high = FusionDecision(
        cell_type="Amacrine",
        confidence="high",
        score=0.9,
        method="marker_scoring_high",
        n_markers_found=3,
        ai_agreed=False,
        ai_suggested="",
        explanation="Strong TFAP2A.",
        alternative_rules=[],
        cell_category="Broad_Neuron",
        tier="L2",
    )
    return [kadp, low, high]


def _run_kadp_wiring_engine(adata, cfg, ai_response: str):
    """Drive the real ``run_unified_annotation``; capture both fuse calls."""
    decisions = _kadp_wiring_decisions()
    fuse_calls: list = []

    def _fake_fuse_all(*args, **kwargs):
        fuse_calls.append(kwargs)
        if len(fuse_calls) == 1:
            return list(decisions), {}
        return list(decisions), {}

    logger = MagicMock()
    with (
        patch("core.kb.load_kb", return_value=_make_kadp_tiered_kb()),
        patch(
            "core.annotation.scoring.score_cluster_against_kb",
            side_effect=lambda *a, **k: dict(_ai_gate_canned_scores()),
        ),
        patch("rna.utils.evidence_fusion.fuse_all_clusters", side_effect=_fake_fuse_all),
        patch("core.ai.caller.ai_query", return_value=ai_response),
        patch("core.annotation.engine.safe_plot"),
    ):
        result = run_unified_annotation(adata, cfg, logger)
    return result, adata, fuse_calls


def test_kadp_wiring_dual_call_mirror_and_obs(tmp_path) -> None:
    """todo 5: BOTH fuse_all_clusters calls receive the same KADPConfig, and
    the KADP decision survives the AI second pass with cell_state=differentiating,
    annot_method=developmental_potency, tiered metadata untouched, potency in
    annot_evidence / ann_records / cell_metadata, and a kadp_precursor review.
    """
    cfg = _make_kadp_wiring_cfg(tmp_path / "kadp")
    ai_response = json.dumps({"0": {"cell_type": "NRPC"}, "1": {"cell_type": "RGC"}})
    result, adata, fuse_calls = _run_kadp_wiring_engine(_make_ai_gate_adata(), cfg, ai_response)

    # ── AI fallback fired → second pass ran ──
    assert len(fuse_calls) == 2, "a low-confidence cluster must trigger the AI second pass"

    # ── dual-call mirroring (Oracle r1 BLOCKER 1) ──
    kadp_cfg_0 = fuse_calls[0].get("kadp_cfg")
    kadp_cfg_1 = fuse_calls[1].get("kadp_cfg")
    assert kadp_cfg_0 is not None, "first fuse_all_clusters call must receive kadp_cfg"
    assert kadp_cfg_0 is kadp_cfg_1, "both calls must share the SAME KADPConfig instance"
    assert kadp_cfg_0.enabled is True
    assert kadp_cfg_0.ratio_threshold == 2.0
    assert kadp_cfg_0.abs_threshold == 0.6
    assert kadp_cfg_0.gap_threshold == 0.1
    assert kadp_cfg_0.use_gap_criterion is False

    # ── KADP naming preserved after the AI second pass ──
    assert result is not None
    assert result["0"].method == "developmental_potency"
    assert result["0"].cell_type == "NRPC"

    # ── tiered-block exemption (Oracle r3 MAJOR 2): tiered metadata untouched ──
    d0 = result["0"]
    assert d0.tier == "", f"tiered block must not overwrite KADP tier: {d0.tier!r}"
    assert d0.consensus == "", f"KADP consensus must stay default: {d0.consensus!r}"
    assert d0.n_sources == 0, f"KADP n_sources must stay default: {d0.n_sources!r}"
    assert d0.subtype_resolution == "", (
        f"KADP subtype_resolution must stay default: {d0.subtype_resolution!r}"
    )

    # ── obs-level four-value assertions ──
    mask0 = adata.obs["leiden"].astype(str) == "0"
    assert adata.obs.loc[mask0, "cell_state"].unique().tolist() == ["differentiating"]
    assert adata.obs.loc[mask0, "annot_method"].unique().tolist() == ["developmental_potency"], (
        "_clean_method must keep KADP annot_method (not marker_scoring+ai)"
    )
    ev0 = json.loads(adata.obs.loc[mask0, "annot_evidence"].iloc[0])
    assert ev0["potency"] == _KADP_POTENCY, ev0

    # ── ann_records / cell_type_annotations.csv contain potency ──
    ann_df = pd.read_csv(os.path.join(cfg.table_dir, "cell_type_annotations.csv"))
    row0 = ann_df[ann_df["cluster"].astype(str) == "0"].iloc[0]
    assert json.loads(row0["potency"]) == _KADP_POTENCY
    row1 = ann_df[ann_df["cluster"].astype(str) == "1"].iloc[0]
    assert pd.isna(row1["potency"]) or row1["potency"] == "", (
        "non-KADP ann_records potency must be empty"
    )
    # ── cell_metadata.csv potency column (single JSON-string column) ──
    meta_df = pd.read_csv(os.path.join(cfg.table_dir, "cell_metadata.csv"))
    assert meta_df.loc[mask0.values, "potency"].iloc[0] == json.dumps(_KADP_POTENCY)
    non_kadp = meta_df.loc[~mask0.values, "potency"]
    assert non_kadp.isna().all(), "non-KADP cell_metadata potency must be empty"
    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    with open(quality_path, encoding="utf-8") as f:
        quality = json.load(f)
    q0 = next(q for q in quality["review_queue"] if q["cluster"] == "0")
    assert q0["reason"] == "kadp_precursor", q0
    assert q0["n_tied_types"] == 0, q0
    assert q0["top_types"] == [], q0


# ═══════════════════════════════════════════════════════════════════════
#  Todo 10 — engine/schema METC wiring: dual-call mirroring + second-pass
#  quality capture (Oracle r3 MINOR 5) + review_queue metc reasons +
#  tie-gate extension (F13) + source_votes in annot_evidence +
#  harmonization_rate plumbing (None guard).
# ═══════════════════════════════════════════════════════════════════════


def _make_metc_wiring_cfg(
    tmp_path,
    *,
    metc_enabled: bool = True,
    min_sources: int = 3,
    min_distinct: int = 3,
) -> SimpleNamespace:
    """Developing-tissue config: AI on, celltypist off, METC on with
    explicit thresholds so the captured METCConfig can be value-asserted."""
    table_dir = tmp_path / "results" / "tables"
    figure_dir = tmp_path / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    return SimpleNamespace(
        tissue_kb="retina",
        tissue="retina",
        species="human",
        target_class="",
        target_order="",
        tissue_maturity="developing",
        table_dir=str(table_dir),
        figure_dir=str(figure_dir),
        interactive=False,
        marker=SimpleNamespace(
            candidate_pool_expand_steps=[5],
            expert_rule_top_n=0,
            expert_rule_pval_cutoff=0.0,
            expert_rule_strictness="default",
            developmental_mode=False,
        ),
        annotation=SimpleNamespace(
            method="kb_unified",
            celltypist=SimpleNamespace(enabled=False, model="", majority_voting=False),
            kadp_enabled=False,
            metc_enabled=metc_enabled,
            metc_min_sources=min_sources,
            metc_min_distinct_transition=min_distinct,
        ),
        ai=SimpleNamespace(enabled=True, ai_annotation=True, unconstrained_annotation=False),
        plot=SimpleNamespace(figure_dpi=150, figure_format="pdf", figure_transparent=True),
        execution=SimpleNamespace(random_seed=42),
    )


_METC_SOURCE_VOTES = {
    "divergent": {"marker": "RGC", "expert": None, "ai": "Amacrine", "celltypist": "NRPC"},
    "2way": {"marker": "RGC", "expert": None, "ai": "Amacrine", "celltypist": None},
    "consensus": {"marker": "RGC", "expert": None, "ai": "RGC", "celltypist": "RGC"},
}

_METC_QUALITY_FIRST = {
    "annotated_by_rule": 0,
    "annotated_by_scoring": 1,
    "unknown": 0,
    "ambiguity": 0,
    "ai_agreed": 0,
    "total": 3,
    "diagnostic_summary": {},
    "celltypist": True,
    "harmonization_rate": 0.0,
}

_METC_QUALITY_SECOND = {
    "annotated_by_rule": 0,
    "annotated_by_scoring": 1,
    "unknown": 0,
    "ambiguity": 0,
    "ai_agreed": 0,
    "total": 3,
    "diagnostic_summary": {},
    "celltypist": True,
    "harmonization_rate": 0.6667,
}


def _metc_wiring_decisions() -> list:
    """Canned fusion output covering the three METC review reasons.

    Cluster 0 → metc_divergent (transition_state), 1 → metc_2way (ambiguous),
    2 → metc_consensus (rescued marker_scoring).  Each carries a fresh
    ``source_votes`` dict so annot_evidence can round-trip it."""
    divergent = FusionDecision(
        cell_type="transitional: RGC/Amacrine",
        confidence="transition",
        score=0.0,
        method="transition_state",
        n_markers_found=0,
        ai_agreed=False,
        ai_suggested="",
        explanation="METC divergent: 3 sources split across 3 labels.",
        alternative_rules=[],
        diagnostic=DiagnosticInfo(
            category="metc_divergent",
            top_competitors=[
                {"cell_type": "RGC", "score": 1},
                {"cell_type": "Amacrine", "score": 1},
                {"cell_type": "NRPC", "score": 1},
            ],
            detail="METC divergent votes.",
        ),
        review_reason="metc_divergent",
        source_votes=dict(_METC_SOURCE_VOTES["divergent"]),
    )
    two_way = FusionDecision(
        cell_type="Unknown",
        confidence="low",
        score=0.0,
        method="ambiguous",
        n_markers_found=0,
        ai_agreed=False,
        ai_suggested="",
        explanation="METC 2-way split: 'RGC' (1 vote) vs 'Amacrine' (1 vote).",
        alternative_rules=[],
        diagnostic=DiagnosticInfo(
            category="metc_2way",
            top_competitors=[
                {"cell_type": "RGC", "score": 1},
                {"cell_type": "Amacrine", "score": 1},
            ],
            detail="METC 2-way votes.",
        ),
        review_reason="metc_2way",
        source_votes=dict(_METC_SOURCE_VOTES["2way"]),
    )
    consensus = FusionDecision(
        cell_type="RGC",
        confidence="medium",
        score=0.0,
        method="marker_scoring",
        n_markers_found=0,
        ai_agreed=False,
        ai_suggested="",
        explanation="METC consensus: all 3 sources agree on 'RGC'.",
        alternative_rules=[],
        diagnostic=DiagnosticInfo(
            category="metc_consensus",
            top_competitors=[],
            detail="METC consensus votes.",
        ),
        review_reason="metc_consensus",
        source_votes=dict(_METC_SOURCE_VOTES["consensus"]),
    )
    return [divergent, two_way, consensus]


def _run_metc_wiring_engine(adata, cfg, ai_response: str, second_quality=None):
    """Drive the real ``run_unified_annotation``; capture BOTH fuse calls and
    keep the first/second-pass quality dicts distinct so the second-pass
    capture can be asserted on the written 05_annotation_quality.json."""
    decisions = _metc_wiring_decisions()
    fuse_calls: list = []
    q2 = second_quality if second_quality is not None else dict(_METC_QUALITY_SECOND)

    def _fake_fuse_all(*args, **kwargs):
        fuse_calls.append(kwargs)
        if len(fuse_calls) == 1:
            return list(decisions), dict(_METC_QUALITY_FIRST)
        return list(decisions), dict(q2)

    logger = MagicMock()
    with (
        patch("core.kb.load_kb", return_value=_make_ai_gate_kb()),
        patch(
            "core.annotation.scoring.score_cluster_against_kb",
            side_effect=lambda *a, **k: dict(_ai_gate_canned_scores()),
        ),
        patch("rna.utils.evidence_fusion.fuse_all_clusters", side_effect=_fake_fuse_all),
        patch("core.ai.caller.ai_query", return_value=ai_response),
        patch("core.annotation.engine.safe_plot"),
    ):
        result = run_unified_annotation(adata, cfg, logger)
    return result, adata, fuse_calls


def test_metc_wiring_dual_call_mirror_and_second_pass_quality(tmp_path) -> None:
    """todo 10: BOTH fuse_all_clusters calls receive the SAME METCConfig
    instance and ``return_quality=True``; the quality written to
    05_annotation_quality.json is the SECOND (AI-enhanced) pass's quality
    (Oracle r3 MINOR 5) carrying the harmonization_rate key; the review_queue
    surfaces metc_divergent/metc_2way/metc_consensus with tie detail (F13);
    annot_evidence carries source_votes."""
    cfg = _make_metc_wiring_cfg(tmp_path / "metc")
    ai_response = json.dumps(
        {
            "0": {"cell_type": "RGC"},
            "1": {"cell_type": "Amacrine"},
            "2": {"cell_type": "RGC"},
        }
    )
    result, adata, fuse_calls = _run_metc_wiring_engine(_make_ai_gate_adata(), cfg, ai_response)

    # ── AI fallback fired → second pass ran ──
    assert len(fuse_calls) == 2, "low-conf/transition candidates must trigger the AI second pass"

    # ── dual-call mirroring (Oracle r1 BLOCKER 1, extended to METC) ──
    metc_cfg_0 = fuse_calls[0].get("metc_cfg")
    metc_cfg_1 = fuse_calls[1].get("metc_cfg")
    assert metc_cfg_0 is not None, "first fuse_all_clusters call must receive metc_cfg"
    assert metc_cfg_0 is metc_cfg_1, "both calls must share the SAME METCConfig instance"
    assert metc_cfg_0.enabled is True
    assert metc_cfg_0.min_sources == 3
    assert metc_cfg_0.min_distinct_transition == 3

    # ── both calls request quality; the SECOND one is written (Oracle r3 MINOR 5) ──
    assert fuse_calls[0].get("return_quality") is True
    assert fuse_calls[1].get("return_quality") is True
    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    with open(quality_path, encoding="utf-8") as f:
        quality = json.load(f)
    assert quality["harmonization_rate"] == _METC_QUALITY_SECOND["harmonization_rate"], (
        "written quality must come from the SECOND (AI-enhanced) fuse pass"
    )

    # ── review_queue: metc_divergent reason + tie detail populated (F13) ──
    q0 = next(q for q in quality["review_queue"] if q["cluster"] == "0")
    assert q0["reason"] == "metc_divergent", q0
    assert q0["n_tied_types"] == 3, q0
    assert q0["top_types"] == ["RGC", "Amacrine", "NRPC"], q0

    # ── review_queue: metc_2way / metc_consensus reasons ──
    q1 = next(q for q in quality["review_queue"] if q["cluster"] == "1")
    assert q1["reason"] == "metc_2way", q1
    assert q1["n_tied_types"] == 2, q1
    q2 = next(q for q in quality["review_queue"] if q["cluster"] == "2")
    assert q2["reason"] == "metc_consensus", q2
    assert q2["n_tied_types"] == 0, q2

    # ── obs-level: metc method + annot_evidence source_votes ──
    assert result is not None
    assert result["0"].method == "transition_state"
    assert result["0"].review_reason == "metc_divergent"
    mask0 = adata.obs["leiden"].astype(str) == "0"
    assert adata.obs.loc[mask0, "annot_method"].unique().tolist() == ["transition_state"]
    ev0 = json.loads(adata.obs.loc[mask0, "annot_evidence"].iloc[0])
    assert ev0["source_votes"] == _METC_SOURCE_VOTES["divergent"], ev0
    mask1 = adata.obs["leiden"].astype(str) == "1"
    ev1 = json.loads(adata.obs.loc[mask1, "annot_evidence"].iloc[0])
    assert ev1["source_votes"] == _METC_SOURCE_VOTES["2way"], ev1


def test_metc_quality_harmonization_rate_none_guard(tmp_path) -> None:
    """No celltypist labels → second-pass quality carries harmonization_rate
    None and the written JSON keeps it null (Oracle r3 MINOR 6) — the engine
    plumbing passes the None through without crashing."""
    cfg = _make_metc_wiring_cfg(tmp_path / "metc_none")
    ai_response = json.dumps({"0": {"cell_type": "RGC"}, "1": {"cell_type": "Amacrine"}})
    q2_none = dict(_METC_QUALITY_SECOND)
    q2_none["harmonization_rate"] = None
    _run_metc_wiring_engine(_make_ai_gate_adata(), cfg, ai_response, second_quality=q2_none)
    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    with open(quality_path, encoding="utf-8") as f:
        quality = json.load(f)
    assert quality["harmonization_rate"] is None, quality


# ═══════════════════════════════════════════════════════════════════════
#  Todo 11 — CellTypist AnnotationResult label capture (engine fix).
#
#  celltypist >= 1.6 returns an AnnotationResult WITHOUT mutating adata.obs.
#  The engine must capture ``_res = celltypist.annotate(...)`` and read per-
#  cell labels from ``_res.predicted_labels`` (a DataFrame whose column is
#  "majority_voting" when majority_voting=True else "predicted_labels", index
#  aligned with adata.obs_names).  Pre-fix the result was discarded and the
#  ``_label_col in adata.obs`` guard was ALWAYS False -> celltypist_results
#  stayed {} -> the CellTypist METC source abstained (n_spoke=2<3).
#
#  The engine imports celltypist LAZILY inside run_unified_annotation, and
#  importing the real package under pytest's filterwarnings=error raises on
#  scanpy's ``__version__`` FutureWarning.  These tests therefore inject a
#  FAKE ``celltypist`` module tree into ``sys.modules`` so the lazy import
#  resolves to the mock without ever touching the installed package.
# ═══════════════════════════════════════════════════════════════════════


def _make_celltypist_cfg(tmp_path, *, majority_voting: bool = True) -> SimpleNamespace:
    """AI-on, celltypist-ON config (defaults kadp/metc off)."""
    cfg = _make_ai_gate_cfg(tmp_path)
    cfg.annotation.celltypist = SimpleNamespace(
        enabled=True,
        model="Fetal_Human_Retina.pkl",
        majority_voting=majority_voting,
    )
    return cfg


def _celltypist_canned_scores() -> dict:
    """Low non-zero scores: real fusion yields confidence 'low' so the AI
    fallback fires and the SECOND-pass quality (with harmonization_rate) is
    the written one."""
    return {
        "RGC": Score(
            score=0.30,
            p_value=0.05,
            method="test",
            n_markers_found=1,
            negative_penalty=False,
            tier="L3",
            private_markers_hit=0,
            consensus="low",
            n_sources=1,
        ),
        "Amacrine": Score(
            score=0.25,
            p_value=0.05,
            method="test",
            n_markers_found=1,
            negative_penalty=False,
            tier="L3",
            private_markers_hit=0,
            consensus="low",
            n_sources=1,
        ),
    }


class _FakeCelltypistResult:
    """Minimal stand-in for celltypist.AnnotationResult (>= 1.6)."""

    def __init__(self, predicted_labels):
        self.predicted_labels = predicted_labels


def _celltypist_labels_df(adata, column: str, per_cluster: dict) -> pd.DataFrame:
    """Build a predicted_labels DataFrame aligned with adata.obs_names.

    ``per_cluster`` maps cluster string -> label; the label is repeated for
    every cell of that cluster so the engine's per-cluster mode() returns
    exactly that label."""
    labels = []
    for cl in adata.obs["leiden"].astype(str):
        labels.append(per_cluster[cl])
    return pd.DataFrame({column: labels}, index=adata.obs_names)


def _importlib_util_spec(name: str):
    """A minimal ModuleSpec so importlib.util.find_spec does not raise."""
    import importlib.machinery as _machinery

    return _machinery.ModuleSpec(name, loader=None)


def _inject_fake_celltypist(annotate_mock):
    """Register a fake ``celltypist`` package tree in sys.modules.

    The engine lazily does ``import celltypist`` then
    ``celltypist.models.Model.load(model=...)`` and
    ``celltypist.annotate(...)``.  Provide a ModuleType-based fake so nothing
    from the real package is ever imported (avoiding the scanpy FutureWarning
    that pytest's filterwarnings=error turns into an exception).  A
    ``__spec__`` is set on the fake modules so ``importlib.util.find_spec``
    (called by ``require_celltypist``) does not raise.
    """
    fake_root = _types_mod.ModuleType("celltypist")
    fake_models = _types_mod.ModuleType("celltypist.models")
    fake_root.__spec__ = _importlib_util_spec("celltypist")
    fake_models.__spec__ = _importlib_util_spec("celltypist.models")
    fake_model_cls = MagicMock(name="Model")
    fake_model_cls.load = MagicMock(return_value=MagicMock(name="loaded_model"))
    fake_models.Model = fake_model_cls
    fake_root.models = fake_models
    fake_root.annotate = annotate_mock
    return {
        "celltypist": fake_root,
        "celltypist.models": fake_models,
    }


def _run_celltypist_capture_engine(
    adata,
    cfg,
    fake_result=None,
    ai_response: str = "{}",
    annotate_mock=None,
):
    """Drive the REAL run_unified_annotation with celltypist enabled.

    ``celltypist.annotate`` is mocked via an injected fake package tree (see
    ``_inject_fake_celltypist``) — pass either ``fake_result`` (returned by
    the mock) or a pre-built ``annotate_mock`` (e.g. raising).
    fuse_all_clusters is a spy that captures the ``celltypist_results`` kwarg
    AND delegates to the real function so the quality dict (celltypist /
    harmonization_rate) is computed by the frozen t8/t9 code."""
    fuse_calls = []
    qualities = []
    real_fuse = _ef.fuse_all_clusters

    def _spy_fuse(*args, **kwargs):
        fuse_calls.append(kwargs)
        out = real_fuse(*args, **kwargs)
        if isinstance(out, tuple) and len(out) == 2:
            qualities.append(out[1])
        return out

    if annotate_mock is None:
        annotate_mock = MagicMock(return_value=fake_result)
    fake_modules = _inject_fake_celltypist(annotate_mock)

    logger = MagicMock()
    with (
        patch.dict(sys.modules, fake_modules),
        patch("core.kb.load_kb", return_value=_make_ai_gate_kb()),
        patch("core.kb.load_synonyms", return_value={}),
        patch(
            "core.annotation.scoring.score_cluster_against_kb",
            side_effect=lambda *a, **k: dict(_celltypist_canned_scores()),
        ),
        patch("rna.utils.evidence_fusion.fuse_all_clusters", side_effect=_spy_fuse),
        patch("core.ai.caller.ai_query", return_value=ai_response),
        patch("core.annotation.engine.safe_plot"),
    ):
        result = run_unified_annotation(adata, cfg, logger)
    return result, fuse_calls, qualities


def test_celltypist_capture_majority_voting_column(tmp_path) -> None:
    """todo 11 fix: with majority_voting=True the engine must read
    ``_res.predicted_labels['majority_voting']`` and fill celltypist_results
    with the per-cluster mode label.  Pre-fix this stayed {} because the
    result was discarded and `_label_col in adata.obs` was always False."""
    cfg = _make_celltypist_cfg(tmp_path / "ct_mv", majority_voting=True)
    adata = _make_ai_gate_adata()
    fake = _FakeCelltypistResult(
        _celltypist_labels_df(adata, "majority_voting", {"0": "RGC", "1": "Amacrine", "2": "RGC"})
    )
    ai_response = json.dumps(
        {
            "0": {"cell_type": "RGC"},
            "1": {"cell_type": "Amacrine"},
            "2": {"cell_type": "RGC"},
        }
    )
    result, fuse_calls, qualities = _run_celltypist_capture_engine(adata, cfg, fake, ai_response)

    assert result is not None
    expected_ct = {"0": "RGC", "1": "Amacrine", "2": "RGC"}
    assert fuse_calls, "fuse_all_clusters must be called"
    assert fuse_calls[0].get("celltypist_results") == expected_ct, fuse_calls[0]
    # both passes mirror the same captured dict (Oracle r1 BLOCKER 1)
    assert fuse_calls[-1].get("celltypist_results") == expected_ct
    # real fusion quality: celltypist True + numeric harmonization_rate
    assert qualities, "return_quality must be requested"
    assert qualities[-1]["celltypist"] is True
    assert isinstance(qualities[-1]["harmonization_rate"], float)
    assert qualities[-1]["harmonization_rate"] > 0
    # written JSON carries the numeric harmonization_rate
    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    with open(quality_path, encoding="utf-8") as f:
        written = json.load(f)
    assert written["harmonization_rate"] == qualities[-1]["harmonization_rate"]


def test_celltypist_capture_predicted_labels_column(tmp_path) -> None:
    """todo 11 fix: majority_voting=False path uses the "predicted_labels"
    column of the AnnotationResult."""
    cfg = _make_celltypist_cfg(tmp_path / "ct_pl", majority_voting=False)
    adata = _make_ai_gate_adata()
    fake = _FakeCelltypistResult(
        _celltypist_labels_df(adata, "predicted_labels", {"0": "RGC", "1": "RGC", "2": "Amacrine"})
    )
    ai_response = json.dumps(
        {
            "0": {"cell_type": "RGC"},
            "1": {"cell_type": "RGC"},
            "2": {"cell_type": "Amacrine"},
        }
    )
    result, fuse_calls, qualities = _run_celltypist_capture_engine(adata, cfg, fake, ai_response)

    assert result is not None
    assert fuse_calls[0].get("celltypist_results") == {"0": "RGC", "1": "RGC", "2": "Amacrine"}
    assert qualities[-1]["celltypist"] is True
    assert qualities[-1]["harmonization_rate"] > 0


def test_celltypist_capture_degrades_when_labels_missing(tmp_path) -> None:
    """todo 11: empty/absent predicted_labels -> celltypist_results stays {},
    no crash, quality celltypist False (the existing degrade path)."""
    cfg = _make_celltypist_cfg(tmp_path / "ct_deg", majority_voting=True)
    adata = _make_ai_gate_adata()
    ai_response = json.dumps({})

    # absent predicted_labels column: DataFrame without the expected column
    empty_df = pd.DataFrame(index=adata.obs_names)
    fake = _FakeCelltypistResult(empty_df)
    result, fuse_calls, qualities = _run_celltypist_capture_engine(adata, cfg, fake, ai_response)
    assert result is not None
    assert fuse_calls[0].get("celltypist_results") == {}, "no labels -> source abstains"
    assert qualities[-1]["celltypist"] is False
    assert qualities[-1]["harmonization_rate"] is None

    # annotate raising -> caught, no crash, empty results
    def _boom(*a, **k):
        raise RuntimeError("model download failed")

    cfg2 = _make_celltypist_cfg(tmp_path / "ct_boom", majority_voting=True)
    result2, fuse_calls2, _ = _run_celltypist_capture_engine(
        adata, cfg2, annotate_mock=MagicMock(side_effect=_boom)
    )
    assert result2 is not None, "annotate failure must degrade, not crash"
    assert fuse_calls2[0].get("celltypist_results") == {}, "caught -> source abstains"
