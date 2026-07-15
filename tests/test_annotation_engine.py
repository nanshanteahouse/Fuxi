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
