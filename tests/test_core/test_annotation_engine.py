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
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from core.annotation.engine import _map_cell_state, _write_quality_report
from core.annotation.scoring import Score
from rna.utils.evidence_fusion import DiagnosticInfo

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
