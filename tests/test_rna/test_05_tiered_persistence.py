"""Tests for tiered near-miss subtype-candidate persistence in the annotation engine.

Todo 2 (tiered-subtype-reuse): ``run_unified_annotation`` must persist the
``subtype_candidates`` list computed by ``resolve_tiered_label`` into each
cell's ``annot_evidence`` obs JSON (new key) and into
``cell_type_annotations.csv`` (new column). Both additions are guarded on
``if tiered_candidates:`` so that a KB without ``_hierarchy`` behaves
byte-for-byte unchanged (no key in JSON, no column in CSV).
"""

import json
import os
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
from anndata import AnnData

from core.annotation.engine import run_unified_annotation
from core.annotation.scoring import Score


def _minimal_kb(with_hierarchy: bool = True) -> dict:
    """Minimal retina-like KB with ``_hierarchy`` + subtype private markers.

    RGC_Foxp2's private marker ``FOXP2`` is deliberately NOT present in the
    synthetic gene matrix, so gate B (private-marker hit) and gate C
    (consensus among hit private markers) fail *by construction* for every
    cluster — gates B/C are re-derived from the real cluster top genes, never
    from the canned ``Score`` fields.
    """
    kb = {
        "Broad_Neuron": {
            "parent": "",
            "markers": {},
            "_private_markers": [],
            "consensus_levels": {},
            "marker_weights": {},
        },
        "RGC": {
            "parent": "Broad_Neuron",
            "markers": {},
            "_private_markers": ["RBPMS", "POU4F1"],
            "consensus_levels": {"RBPMS": "high", "POU4F1": "high"},
            "marker_weights": {"RBPMS": 10},
        },
        "RGC_Foxp2": {
            "parent": "RGC",
            "markers": {},
            "_private_markers": ["FOXP2"],
            "consensus_levels": {"FOXP2": "high"},
            "marker_weights": {"FOXP2": 2},
        },
        "RGC_Alpha": {
            "parent": "RGC",
            "markers": {},
            "_private_markers": [],
            "consensus_levels": {},
            "marker_weights": {},
        },
    }
    if with_hierarchy:
        kb["_hierarchy"] = {
            "categories": {
                "Neuron": {
                    "members": ["RGC", "RGC_Foxp2", "RGC_Alpha"],
                    "subtypes": {
                        "RGC": {"members": ["RGC_Foxp2", "RGC_Alpha"]},
                    },
                }
            },
            "incompatible_transitions": [],
        }
    return kb


def _make_cfg(tmp_path) -> SimpleNamespace:
    """Plain-value config (no MagicMock truthiness surprises) for the engine."""
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
        ),
        ai=SimpleNamespace(enabled=False, ai_annotation=False, unconstrained_annotation=False),
        plot=SimpleNamespace(figure_dpi=150, figure_format="pdf", figure_transparent=True),
        execution=SimpleNamespace(random_seed=42),
    )


def _make_adata() -> AnnData:
    """3 leiden clusters × 30 cells with cluster-distinct expression.

    ``var_names`` intentionally excludes ``FOXP2`` (the RGC_Foxp2 private
    marker) so gate B fails deterministically for every cluster.
    """
    rng = np.random.RandomState(0)
    var_names = [
        "RBPMS",
        "POU4F1",
        "ONECUT1",
        "GENE_01",
        "GENE_02",
        "GENE_03",
        "GENE_04",
        "GENE_05",
    ]
    n_cells = 90
    x = rng.poisson(lam=1.0, size=(n_cells, len(var_names))).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = var_names
    adata.obs["leiden"] = np.repeat(["0", "1", "2"], 30)

    c0 = adata.obs["leiden"] == "0"
    c1 = adata.obs["leiden"] == "1"
    # Cluster 0: RBPMS/POU4F1 high; cluster 1: ONECUT1 high; cluster 2: filler only.
    adata.X[c0, 0] = adata.X[c0, 0] + rng.poisson(lam=20.0, size=int(c0.sum())).astype(np.float32)
    adata.X[c0, 1] = adata.X[c0, 1] + rng.poisson(lam=20.0, size=int(c0.sum())).astype(np.float32)
    adata.X[c1, 2] = adata.X[c1, 2] + rng.poisson(lam=20.0, size=int(c1.sum())).astype(np.float32)

    adata.obsm["X_umap"] = rng.randn(n_cells, 2).astype(np.float32)
    return adata


def _canned_scores() -> dict[str, Score]:
    """Same canned Score dict for EVERY cluster (non-zero → no retry, deterministic)."""
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
        "RGC_Foxp2": Score(
            score=0.58,
            p_value=0.05,
            method="test",
            n_markers_found=1,
            negative_penalty=False,
            tier="L3",
            private_markers_hit=0,
            consensus="",
            n_sources=1,
        ),
        "RGC_Alpha": Score(
            score=0.55,
            p_value=0.05,
            method="test",
            n_markers_found=1,
            negative_penalty=False,
            tier="L3",
            private_markers_hit=0,
            consensus="",
            n_sources=1,
        ),
    }


def _run_engine(adata: AnnData, cfg, kb: dict):
    logger = MagicMock()
    with (
        patch("core.kb.load_kb", return_value=kb),
        patch(
            "core.annotation.scoring.score_cluster_against_kb",
            side_effect=lambda *a, **k: dict(_canned_scores()),
        ),
        patch("core.annotation.engine.safe_plot"),
    ):
        return run_unified_annotation(adata, cfg, logger)


def test_subtype_candidates_persist_to_annot_evidence_and_csv(tmp_path) -> None:
    """Every cell's annot_evidence JSON + the CSV carry subtype_candidates.

    The unresolved near-miss candidate renders as ``RGC_Foxp2:0.58[B,C]``
    (gate B = FOXP2 absent from the matrix; gate C = no hit private markers).
    """
    adata = _make_adata()
    cfg = _make_cfg(tmp_path)

    result = _run_engine(adata, cfg, _minimal_kb(with_hierarchy=True))

    assert result is not None

    # (a) every cell's annot_evidence JSON parses and carries the key
    for ev_json in adata.obs["annot_evidence"]:
        ev = json.loads(ev_json)
        assert "subtype_candidates" in ev, f"missing key in {ev!r}"
        cands = ev["subtype_candidates"]
        assert cands, "expected non-empty candidates for every cluster"
        assert {c["type"] for c in cands} == {"RGC_Foxp2", "RGC_Alpha"}, cands

    # (b) cell_type_annotations.csv has the subtype_candidates column
    csv_path = os.path.join(cfg.table_dir, "cell_type_annotations.csv")
    assert os.path.exists(csv_path), csv_path
    df = pd.read_csv(csv_path)
    assert "subtype_candidates" in df.columns, list(df.columns)

    # (c) at least one row renders the near-miss candidate with failed gates
    rendered = df["subtype_candidates"].astype(str).tolist()
    assert any("RGC_Foxp2:0.58[B,C]" in v for v in rendered), rendered


def test_no_hierarchy_leaves_annot_evidence_and_csv_unchanged(tmp_path) -> None:
    """Without ``_hierarchy`` the new key/column must not appear at all."""
    adata = _make_adata()
    cfg = _make_cfg(tmp_path)

    result = _run_engine(adata, cfg, _minimal_kb(with_hierarchy=False))

    assert result is not None

    for ev_json in adata.obs["annot_evidence"]:
        ev = json.loads(ev_json)
        assert "subtype_candidates" not in ev, f"unexpected key in {ev!r}"

    csv_path = os.path.join(cfg.table_dir, "cell_type_annotations.csv")
    assert os.path.exists(csv_path), csv_path
    df = pd.read_csv(csv_path)
    assert "subtype_candidates" not in df.columns, list(df.columns)
