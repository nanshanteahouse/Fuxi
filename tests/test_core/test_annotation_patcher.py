"""Tests for ``rna/utils/annotation_patcher.py`` CSV rewrite column preservation.

Todo 3 of the tiered-subtype-reuse plan: ``_rewrite_annotation_csv`` rebuilds
``cell_type_annotations.csv`` from ``adata.obs`` using only the core columns
(cluster/cell_type/confidence/method/reasoning), silently dropping any other
column. Since todo 2 the annotation engine writes a ``subtype_candidates``
column into this CSV — the rewrite must carry those values through, tolerating
an old CSV that is missing entirely or lacks the column (pre-change output).
"""

import os
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from rna.utils.annotation_patcher import apply_annotation_patches

CORE_COLUMNS = ["cluster", "cell_type", "confidence", "method", "reasoning"]


def _make_adata():
    """Small 3-cluster synthetic AnnData matching the patcher's obs contract."""
    adata = ad.AnnData(X=np.zeros((6, 4)))
    adata.obs["leiden"] = ["0", "0", "1", "1", "2", "2"]
    adata.obs["cell_type"] = ["RGC", "RGC", "Rod", "Rod", "Muller", "Muller"]
    return adata


def _make_cfg(table_dir):
    return SimpleNamespace(table_dir=str(table_dir))


def _seed_old_csv(table_dir, columns, rows):
    table_dir.mkdir(parents=True, exist_ok=True)
    path = os.path.join(str(table_dir), "cell_type_annotations.csv")
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)
    return path


def _run_rewrite(tmp_path, seed=None):
    """Run the patcher CSV rewrite path (via ``apply_annotation_patches``)."""
    table_dir = tmp_path / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    if seed is not None:
        _seed_old_csv(table_dir, seed[0], seed[1])
    cfg = _make_cfg(table_dir)
    apply_annotation_patches(_make_adata(), {"0": "RGC (patched)"}, cfg=cfg)
    return pd.read_csv(os.path.join(str(table_dir), "cell_type_annotations.csv"))


# ═══════════════════════════════════════════════════════════════════════
#  Baseline characterization — passes on unchanged code (green first)
# ═══════════════════════════════════════════════════════════════════════


def test_baseline_rewrite_roundtrips_core_columns(tmp_path):
    """Characterization: the 5 core columns survive a rewrite round-trip.

    Pins the CURRENT rewrite contract. (Today the output schema is exactly
    these 5 columns — extra columns such as ``subtype_candidates`` are
    silently dropped; that dropping is the bug pinned by the tests below.)
    """
    out = _run_rewrite(tmp_path)
    assert all(c in out.columns for c in CORE_COLUMNS)
    assert out["cluster"].astype(str).tolist() == ["0", "1", "2"]
    assert out["cell_type"].tolist() == ["RGC (patched)", "Rod", "Muller"]


# ═══════════════════════════════════════════════════════════════════════
#  Red-first: subtype_candidates must survive the rewrite
# ═══════════════════════════════════════════════════════════════════════


def test_rewrite_preserves_subtype_candidates_from_old_csv(tmp_path):
    """(a) Seed CSV WITH ``subtype_candidates`` → rewrite → column + values survive."""
    seed_columns = CORE_COLUMNS + ["subtype_candidates"]
    seed_rows = [
        ["0", "RGC", "high", "kb_unified", "top", "RGC_Foxp2:0.58[B,C]"],
        ["1", "Rod", "high", "kb_unified", "top", ""],
        ["2", "Muller", "high", "kb_unified", "top", "Muller_Glia:0.42"],
    ]
    out = _run_rewrite(tmp_path, seed=(seed_columns, seed_rows))
    assert "subtype_candidates" in out.columns
    # pandas reads empty strings back as NaN — normalize for comparison.
    got = dict(zip(out["cluster"].astype(str), out["subtype_candidates"].fillna("")))
    assert got["0"] == "RGC_Foxp2:0.58[B,C]"
    assert got["1"] == ""
    assert got["2"] == "Muller_Glia:0.42"


def test_rewrite_with_missing_old_csv_adds_empty_column(tmp_path):
    """(b) Failure case 1: no old CSV → rewrite succeeds, new column empty."""
    out = _run_rewrite(tmp_path)
    assert "subtype_candidates" in out.columns
    assert out["subtype_candidates"].fillna("").eq("").all()


def test_rewrite_with_old_csv_missing_column_is_tolerated(tmp_path):
    """(c) Failure case 2: old CSV WITHOUT the column (pre-change) → succeeds,
    empty column, no exception."""
    seed_columns = CORE_COLUMNS
    seed_rows = [
        ["0", "RGC", "high", "kb_unified", "top"],
        ["1", "Rod", "high", "kb_unified", "top"],
        ["2", "Muller", "high", "kb_unified", "top"],
    ]
    out = _run_rewrite(tmp_path, seed=(seed_columns, seed_rows))
    assert "subtype_candidates" in out.columns
    assert out["subtype_candidates"].fillna("").eq("").all()
