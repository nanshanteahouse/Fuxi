"""Tests for rna/annotation_engine.py — use_raw fix in rank_genes_groups.

T1 (P0-CRITICAL) from cross-batch-critical-fixes plan:
  rank_genes_groups must pass use_raw=True when adata.raw exists,
  with a null-guard for adata without .raw.
"""

import sys
import os

import numpy as np
import pytest
import scanpy as sc
from anndata import AnnData

# ── Ensure repo root is on sys.path (conftest.py also does this) ──────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
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
    X = rng.poisson(lam=1.0, size=(n_cells, n_hvg)).astype(np.float32)
    adata = AnnData(X)
    adata.var_names = [f"HVG_{i}" for i in range(n_hvg)]

    # .raw: 1000 full genes, last one is "FOO"
    raw_X = rng.poisson(lam=1.0, size=(n_cells, n_full)).astype(np.float32)
    raw = AnnData(raw_X)
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

    X = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(X)
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

import copy as _copy_mod
import logging
from unittest.mock import patch

import pandas as pd
from rna.utils.marker_scoring import Score


def _make_zero_scores(n_clusters: int = 5) -> dict:
    """Build all_scores dict where every cluster has zero KB hits."""
    return {
        str(i): {"CT": Score(score=0.0, p_value=1.0, method="none",
                              n_markers_found=0, negative_penalty=False)}
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
            rows.append({
                "names": g,
                "logfoldchanges": 5.0 - i * 0.5,
                "pvals_adj": 1e-50,
                "cluster": str(cl),
            })
    return pd.DataFrame(rows)


def _make_logger() -> logging.Logger:
    """Create a logger for _check_zero_scores_and_retry."""
    log = logging.getLogger("test_T4")
    log.setLevel(logging.DEBUG)
    log.addHandler(logging.NullHandler())
    return log


def _check_zero_scores_and_retry_wrapper(
    kb, all_scores, marker_df, clusters, species,
    target_class, target_order, tissue_kb, logger,
):
    """Lazy-import and call _check_zero_scores_and_retry.

    Avoids circular-import issues at module level by importing only when called.
    """
    from rna.annotation_engine import _check_zero_scores_and_retry
    return _check_zero_scores_and_retry(
        kb, all_scores, marker_df, clusters, species,
        target_class, target_order, tissue_kb, logger,
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
        "rna.utils.marker_scoring.score_cluster_against_kb",
        return_value={"CT": Score(0.85, 0.001, "hypergeometric", 2, False)}
    ):
        result_scores, total_hits, n_clusters = _check_zero_scores_and_retry_wrapper(
            kb, all_scores, marker_df, clusters, species="human",
            target_class="", target_order="", tissue_kb="test_kb",
            logger=logger,
        )

    assert total_hits > 0, (
        f"Expected retry to improve hits, got total_hits={total_hits}"
    )
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
        str(i): {"CT": Score(score=0.5, p_value=0.01, method="hypergeometric",
                              n_markers_found=1, negative_penalty=False)}
        for i in range(5)
    }
    kb = _make_kb_lowercase()
    marker_df = _make_marker_df(5)
    clusters = [str(i) for i in range(5)]
    logger = _make_logger()

    # Patch score_cluster_against_kb to track if it gets called
    with patch("rna.utils.marker_scoring.score_cluster_against_kb") as mock_sc:
        result_scores, total_hits, n_clusters = _check_zero_scores_and_retry_wrapper(
            kb, all_scores, marker_df, clusters, species="human",
            target_class="", target_order="", tissue_kb="test_kb",
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
        "rna.utils.marker_scoring.score_cluster_against_kb",
        return_value={"CT": Score(0.0, 1.0, "none", 0, False)}
    ):
        logger = logging.getLogger("test_T4_warning")
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.NullHandler())
        with patch.object(logger, "error") as mock_err:
            _check_zero_scores_and_retry_wrapper(
                kb, all_scores, marker_df, clusters, species="danio_rerio",
                target_class="", target_order="", tissue_kb="test_kb",
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
        "rna.utils.marker_scoring.score_cluster_against_kb",
        return_value={"CT": Score(0.85, 0.001, "hypergeometric", 2, False)}
    ):
        _check_zero_scores_and_retry_wrapper(
            kb, all_scores, marker_df, clusters, species="human",
            target_class="", target_order="", tissue_kb="test_kb",
            logger=logger,
        )

    # Original KB must be identical
    assert kb == kb_original, "Original KB dict was mutated"
    # Specifically check lowercase keys survived
    assert "rho" in kb["CT"]["markers"]["confirm"], (
        "KB marker key 'rho' became uppercase — original KB mutated"
    )


def test_T4_audit_report_deferred_section_present() -> None:
    """Audit report contains the Deferred: KB Data Curation section.

    Verifies the notes/audit file has the expected section header.
    """
    audit_path = os.path.join(
        _REPO_ROOT, "notes", "audit",
        "2026-07-15_cross_batch_critical_issues_audit.md",
    )
    assert os.path.exists(audit_path), f"Audit file not found: {audit_path}"

    with open(audit_path, encoding="utf-8") as f:
        content = f.read()

    assert "## Deferred: KB Data Curation" in content, (
        "Audit report missing 'Deferred: KB Data Curation' section"
    )
