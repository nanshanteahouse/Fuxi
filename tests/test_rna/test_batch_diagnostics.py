#!/usr/bin/env python3
"""Tests for batch_diagnostics subsampling & simplified categorization (S2)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from rna.utils.batch_diagnostics import (
    _compute_cramer_v,
    diagnose_batch_candidates,
)

sc.settings.verbosity = 0


def _make_adata(n_cells: int = 20_000, rng_seed: int = 42):
    rng = np.random.RandomState(rng_seed)
    n_genes = 30
    x = rng.randn(n_cells, n_genes).astype(np.float32)
    tissue = rng.choice(["TissueA", "TissueB"], n_cells)
    x[tissue == "TissueA", :6] += 2.0
    x[tissue == "TissueB", :6] -= 1.0
    sample = rng.choice([f"S{i}" for i in range(12)], n_cells)
    adata = sc.AnnData(
        X=sp.csr_matrix(x),
        obs=pd.DataFrame(
            {"sample": pd.Categorical(sample), "tissue": pd.Categorical(tissue)},
            index=[f"c{i}" for i in range(n_cells)],
        ),
    )
    sc.pp.pca(adata, n_comps=10, random_state=rng_seed)
    return adata


def test_subsample_preserves_classification() -> None:
    """Full vs max_cells=5000 must yield the same batch/biology judgment."""
    adata = _make_adata()
    full = diagnose_batch_candidates(adata)
    sub = diagnose_batch_candidates(adata, max_cells=5000)
    assert set(full.batch_cols) == set(sub.batch_cols), (
        f"batch mismatch: {full.batch_cols} vs {sub.batch_cols}"
    )
    assert set(full.biology_cols) == set(sub.biology_cols), (
        f"biology mismatch: {full.biology_cols} vs {sub.biology_cols}"
    )
    assert "sample" in sub.batch_cols
    assert "tissue" in sub.biology_cols


def test_subsample_reproducible() -> None:
    """Same random_state -> identical judgments."""
    adata = _make_adata()
    r1 = diagnose_batch_candidates(adata, max_cells=5000, random_state=7)
    r2 = diagnose_batch_candidates(adata, max_cells=5000, random_state=7)
    assert [(d.column, d.judgment) for d in r1.column_diagnoses] == [
        (d.column, d.judgment) for d in r2.column_diagnoses
    ]


def test_cramer_v_perfect_collinear() -> None:
    """Identical columns -> V=1.0 (regression: return line was truncated)."""
    a = pd.Series(pd.Categorical(["X", "Y", "X", "Y"]))
    b = pd.Series(pd.Categorical(["X", "Y", "X", "Y"]))
    assert _compute_cramer_v(a, b) == 1.0


def test_cramer_v_none_collinear() -> None:
    """Independent columns -> V well below 1.0."""
    a = pd.Series(pd.Categorical(["X", "X", "Y", "Y"]))
    b = pd.Series(pd.Categorical(["A", "B", "A", "B"]))
    v = _compute_cramer_v(a, b)
    assert v is not None and v < 1.0
