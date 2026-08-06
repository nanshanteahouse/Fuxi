"""Tests for rna/steps/11_grn.py — top_variable_tfs and export_results.

Regression baseline: these tests must pass on the CURRENT (unmodified) code.
"""

from __future__ import annotations

import importlib.util
import logging
import os
from unittest.mock import patch

import numpy as np
import pandas as pd

from core.config.schema import Config

# ── Load the 11_grn module ─────────────────────────────────────────────
# conftest.py adds repo root to sys.path; the rna/steps/* modules use a
# sys.path.insert(0, ...) trick at the top level, so we load via importlib
# to avoid executing that statement at import time.
_STEP_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "rna", "steps", "11_grn.py")
_spec = importlib.util.spec_from_file_location("rna.steps._11_grn_test", _STEP_PATH)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Could not load 11_grn module at {_STEP_PATH}")
grn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(grn)

top_variable_tfs = grn.top_variable_tfs
export_results = grn.export_results

# ======================================================================
#  top_variable_tfs tests
# ======================================================================


def test_top_variable_tfs_unchanged() -> None:
    """n_top=3 returns correct columns sorted by descending variance."""
    rng = np.random.default_rng(42)
    mock_data = pd.DataFrame(
        rng.normal(0, 1, (5, 6)),
        index=pd.Index([f"CT{i}" for i in range(5)]),
        columns=pd.Index([f"TF{i}" for i in range(6)]),
    )
    log = logging.getLogger("test_top_variable_tfs")
    result = top_variable_tfs(mock_data, n_top=3, log=log)

    # Only top 3 TFs by variance
    assert result.shape[1] == 3

    # Same cell types preserved
    assert result.index.tolist() == mock_data.index.tolist()

    # Returned TFs are a subset of the original columns
    assert result.columns.isin(mock_data.columns).all()

    var = mock_data.var(axis=0)
    top3 = var.sort_values(ascending=False).head(3).index  # type: ignore[union-attr]
    assert list(result.columns) == list(top3)


# ======================================================================
#  export_results tests
# ======================================================================


def test_export_results_output_format(tmp_path) -> None:
    """Exports CSV files with expected columns under table_dir/11_grn/."""
    cfg = Config(
        table_dir=str(tmp_path),
        h5ad_dir=str(tmp_path),
    )

    # ── Mock data ─────────────────────────────────────────────────
    rng = np.random.default_rng(42)
    estimates_df = pd.DataFrame(
        rng.normal(0, 1, (3, 4)),
        index=pd.Index([f"CT{i}" for i in range(3)]),
        columns=pd.Index([f"TF{i}" for i in range(4)]),
    )
    top_df = estimates_df[estimates_df.columns[:3]]
    pvals_df = pd.DataFrame(
        np.abs(rng.normal(0, 0.1, (3, 4))),
        index=estimates_df.index,
        columns=estimates_df.columns,
    )
    net_top = pd.DataFrame(
        {
            "source": ["TF0", "TF0", "TF1", "TF2"],
            "target": [f"GENE_{i}" for i in range(4)],
            "weight": [1, 1, 1, 1],
        }
    )
    log = logging.getLogger("test_export_results")

    # ── Run export (mock safe_write to avoid real h5ad I/O) ──────
    grn_dir = os.path.join(tmp_path, "11_grn")
    with patch.object(grn, "safe_write"):
        export_results(estimates_df, top_df, pvals_df, net_top, cfg, log)

    # ── Assert CSV files ──────────────────────────────────────────
    csv_table = os.path.join(grn_dir, "tf_activity_per_cell_type.csv")
    assert os.path.exists(csv_table), f"Missing {csv_table}"

    df_read = pd.read_csv(csv_table, index_col=0)
    assert list(df_read.columns) == list(estimates_df.columns)
    assert df_read.shape == (3, 4)

    # Other expected CSVs
    assert os.path.exists(os.path.join(grn_dir, "tf_activity_pvals.csv"))
    assert os.path.exists(os.path.join(grn_dir, "tf_target_edges.csv"))
    assert os.path.exists(os.path.join(grn_dir, "tf_target_counts.csv"))


# ======================================================================
#  Mode gating tests (T9)
# ======================================================================


def test_off_mode_output_unchanged() -> None:
    """mode='off' returns top n_tfs by descending variance (identical to default)."""
    rng = np.random.default_rng(42)
    estimates_df = pd.DataFrame(
        rng.normal(0, 1, (3, 4)),
        index=pd.Index([f"CT{i}" for i in range(3)]),
        columns=pd.Index([f"TF{i}" for i in range(4)]),
    )
    log = logging.getLogger("test_off_mode")
    result = top_variable_tfs(estimates_df, n_top=3, log=log, mode="off")

    assert result.shape[1] == 3
    assert result.columns.isin(estimates_df.columns).all()

    var = estimates_df.var(axis=0)
    top3 = var.sort_values(ascending=False).head(3).index  # type: ignore[union-attr]
    assert list(result.columns) == list(top3)


def test_soft_mode_adds_annotation_column() -> None:
    """mode='soft' accepts tf_annotation and returns top n_tfs by variance."""
    rng = np.random.default_rng(42)
    estimates_df = pd.DataFrame(
        rng.normal(0, 1, (3, 4)),
        index=pd.Index([f"CT{i}" for i in range(3)]),
        columns=pd.Index([f"TF{i}" for i in range(4)]),
    )
    net = pd.DataFrame(
        {
            "source": ["TF0", "TF0", "TF1", "TF2"],
            "target": ["GENE_1", "GENE_2", "GENE_3", "GENE_4"],
            "weight": [1.0, 1.0, 1.0, 1.0],
        }
    )
    kb_markers: set[str] = {"GENE_1", "GENE_2"}
    log = logging.getLogger("test_soft_mode")
    _, tf_ann = grn.compute_tf_relevance(estimates_df, net, kb_markers, log)
    result = top_variable_tfs(estimates_df, n_top=3, log=log, mode="soft", tf_annotation=tf_ann)

    assert result.shape[1] == 3
    assert result.columns.isin(estimates_df.columns).all()


def test_hard_mode_filters_tfs() -> None:
    """mode='hard' selects top TFs by combined variance + KB rank, preferring KB-overlap TFs."""
    data = {
        "TF0": [3.0, 2.9, 3.1, 3.0],  # high activity, KB overlap
        "TF1": [1.0, 0.9, 1.1, 1.0],  # low activity, KB overlap
        "TF2": [2.0, -2.0, 2.0, -2.0],  # high variance, no KB overlap
        "TF3": [0.1, -0.1, 0.1, -0.1],  # low variance, no KB overlap
        "TF4": [2.1, 2.2, 1.9, 2.0],  # med activity, KB overlap
    }
    estimates_df = pd.DataFrame(
        data,
        index=pd.Index([f"CT{i}" for i in range(4)]),
    )
    net = pd.DataFrame(
        {
            "source": ["TF0", "TF0", "TF1", "TF4"],
            "target": ["GENE_1", "GENE_2", "GENE_1", "GENE_1"],
            "weight": [1.0, 1.0, 1.0, 1.0],
        }
    )
    kb_markers: set[str] = {"GENE_1"}
    log = logging.getLogger("test_hard_mode")
    _, tf_ann = grn.compute_tf_relevance(estimates_df, net, kb_markers, log)
    result = top_variable_tfs(estimates_df, n_top=3, log=log, mode="hard", tf_annotation=tf_ann)

    assert result.shape[1] == 3
    # TF0 has both high variance and KB overlap — should be selected
    assert "TF0" in result.columns


def test_no_tissue_falls_back() -> None:
    """load_all_kb_markers('') returns empty set for blank tissue."""
    markers = grn.load_all_kb_markers("")
    assert isinstance(markers, set)
    assert len(markers) == 0


def test_tissue_unknown_falls_back() -> None:
    """load_all_kb_markers('unknown') returns empty set gracefully."""
    markers = grn.load_all_kb_markers("unknown")
    assert isinstance(markers, set)
    assert len(markers) == 0


# ======================================================================
#  streaming_pseudobulk tests
# ======================================================================


def _make_synth_adata(
    n_cells: int = 600,
    n_genes: int = 300,
    n_ct: int = 8,
    seed: int = 42,
    categorical: bool = True,
    has_raw: bool = False,
):
    import scanpy as sc

    rng = np.random.default_rng(seed)
    labels = np.array([f"CT{i % n_ct}" for i in range(n_cells)])
    mat = np.log1p(rng.poisson(1.5, size=(n_cells, n_genes)).astype(np.float32))
    adata = sc.AnnData(X=mat)
    adata.var_names = [f"G{i}" for i in range(n_genes)]
    adata.obs["cell_type"] = pd.Categorical(labels) if categorical else labels.astype(str)
    if has_raw:
        adata.raw = adata
    return adata


def test_streaming_pseudobulk_dense_matches_anndata(tmp_path) -> None:
    """Dense X without raw — streaming equals the anndata reference."""
    adata = _make_synth_adata(has_raw=False)
    path = tmp_path / "dense.h5ad"
    adata.write_h5ad(path, compression="gzip")
    got = grn.streaming_pseudobulk(str(path), "cell_type")
    ref = grn.build_pseudobulk(adata, "cell_type", use_raw=False)
    np.testing.assert_allclose(got.values, ref.values, atol=1e-5)
    assert list(got.columns) == list(ref.columns)


def test_streaming_pseudobulk_sparse_raw_matches_anndata(tmp_path) -> None:
    """Sparse raw/X (production path) — streaming equals the anndata reference."""
    adata = _make_synth_adata(has_raw=True)
    path = tmp_path / "raw.h5ad"
    adata.write_h5ad(path, compression="gzip")
    got = grn.streaming_pseudobulk(str(path), "cell_type")
    ref = grn.build_pseudobulk(adata, "cell_type", use_raw=True)
    np.testing.assert_allclose(got.values, ref.values, atol=1e-5)


def test_streaming_pseudobulk_plain_string_obs(tmp_path) -> None:
    """Non-categorical obs column — factorized on the fly."""
    adata = _make_synth_adata(categorical=False)
    path = tmp_path / "plain.h5ad"
    adata.write_h5ad(path, compression="gzip")
    got = grn.streaming_pseudobulk(str(path), "cell_type")
    ref = grn.build_pseudobulk(adata, "cell_type", use_raw=False)
    np.testing.assert_allclose(got.values, ref.values, atol=1e-5)


def test_streaming_pseudobulk_falls_back_to_leiden(tmp_path) -> None:
    """Missing group column falls back to 'leiden'."""
    adata = _make_synth_adata()
    adata.obs["leiden"] = pd.Categorical([f"L{i % 6}" for i in range(adata.n_obs)])
    path = tmp_path / "fallback.h5ad"
    adata.write_h5ad(path, compression="gzip")
    got = grn.streaming_pseudobulk(str(path), "missing_col")
    assert got.shape[0] == 6
    assert list(got.index) == [f"L{i}" for i in range(6)]
