"""Tests for rna/steps/05_annotate_major.py — use_raw fix in ai_annotate.

T1 (P0-CRITICAL) from cross-batch-critical-fixes plan:
  ``ai_annotate()`` must call ``rank_genes_groups`` with ``use_raw=True``
  when ``adata.raw`` exists.
"""

import importlib.util
import json
import os
import shutil
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

# ── Ensure repo root is on sys.path (conftest.py also does this) ──────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 05_annotate_major module via file path ───────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "05_annotate_major.py")
_spec = importlib.util.spec_from_file_location("rna.steps._05_annotate_major_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
ai_annotate = _mod.ai_annotate
_warn_if_low_coverage = _mod._warn_if_low_coverage


# ── Shared helpers ────────────────────────────────────────────────────


def _make_test_adata(with_raw: bool = True) -> AnnData:
    """Create a synthetic AnnData with leiden clusters and optional .raw."""
    rng = np.random.RandomState(42)
    n_cells = 30
    n_genes = 100

    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs["leiden"] = rng.choice(["0", "1", "2"], n_cells)
    adata.obsm["X_umap"] = rng.randn(n_cells, 2)

    if with_raw:
        n_raw = 200
        raw_x = rng.poisson(lam=1.0, size=(n_cells, n_raw)).astype(np.float32)
        raw = AnnData(raw_x)
        raw.var_names = [f"RAW_{i}" for i in range(n_raw)]
        adata.raw = raw

    return adata


def _make_config() -> MagicMock:
    """Create a minimal Config mock sufficient for ai_annotate."""
    cfg = MagicMock()
    cfg.ai = MagicMock()
    cfg.ai.max_tokens = 4096
    cfg.ai.model = "test-model"
    cfg.tissue = "test_tissue"
    cfg.species = "test_species"
    cfg.table_dir = "/tmp"
    cfg.figure_dir = "/tmp"
    return cfg


def _setup_rgg_result(adata: AnnData) -> None:
    """Populate adata.uns['rank_genes_groups'] in scanpy's recarray format.

    Scanpy stores each metric as ``np.recarray`` with named object fields,
    one per group.  Shape is (n_top_genes,).
    """
    groups = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
    len(groups)
    n_top = 5

    # dtype: one object field per group name
    dtype = [(str(g), "O") for g in groups]

    # Build recarray rows: each row is a tuple with one value per group
    names_rows = []
    scores_rows = []
    pvals_rows = []
    pvals_adj_rows = []
    lfc_rows = []

    rng = np.random.RandomState(99)
    for i in range(n_top):
        names_rows.append(tuple(f"MARKER_{g}_{i}" for g in groups))
        scores_rows.append(tuple(float(rng.uniform(0.5, 5.0)) for _ in groups))
        pvals_rows.append(tuple(float(rng.uniform(1e-10, 0.05)) for _ in groups))
        pvals_adj_rows.append(tuple(float(rng.uniform(1e-8, 0.1)) for _ in groups))
        lfc_rows.append(tuple(float(rng.randn()) for _ in groups))

    adata.uns["rank_genes_groups"] = {
        "names": np.rec.array(names_rows, dtype=dtype),
        "scores": np.rec.array(scores_rows, dtype=dtype),
        "pvals": np.rec.array(pvals_rows, dtype=dtype),
        "pvals_adj": np.rec.array(pvals_adj_rows, dtype=dtype),
        "logfoldchanges": np.rec.array(lfc_rows, dtype=dtype),
        "params": {
            "groupby": "leiden",
            "method": "wilcoxon",
            "use_raw": True,
        },
    }


# ── Tests ─────────────────────────────────────────────────────────────


def test_T1_ai_annotate_use_raw() -> None:
    """Happy path: ai_annotate calls rank_genes_groups with use_raw=True.

    Given: synthetic adata with .leiden clusters and .raw (200 genes).
    When:  ai_annotate() runs.
    Then:  sc.tl.rank_genes_groups is called with use_raw=True.
    """
    adata = _make_test_adata(with_raw=True)
    cfg = _make_config()
    logger = MagicMock()

    # Track rank_genes_groups call kwargs and set up uns for downstream
    call_kwargs: dict = {}

    def _tracking_rgg(adata_, groupby="leiden", method="wilcoxon", **kwargs):
        call_kwargs.update(kwargs)
        _setup_rgg_result(adata_)

    ai_return = json.dumps(
        {
            "0": {
                "cell_type": "T cell",
                "state": "active",
                "subtype": "CD8+",
                "confidence": "high",
                "reasoning": "markers match",
            },
            "1": {
                "cell_type": "B cell",
                "state": "resting",
                "subtype": "naive",
                "confidence": "high",
                "reasoning": "markers match",
            },
            "2": {
                "cell_type": "NK cell",
                "state": "active",
                "subtype": "CD56+",
                "confidence": "high",
                "reasoning": "markers match",
            },
        }
    )

    with (
        patch.object(_mod.sc.tl, "rank_genes_groups", side_effect=_tracking_rgg),
        patch("core.ai_caller.ai_query", return_value=ai_return),
        patch.object(_mod, "safe_plot"),
    ):
        result = ai_annotate(adata, cfg, logger)

    assert result is not None, "ai_annotate should return a dict on success"
    assert call_kwargs.get("use_raw") is True, (
        f"rank_genes_groups should be called with use_raw=True, "
        f"got use_raw={call_kwargs.get('use_raw')!r}"
    )


# ═══════════════════════════════════════════════════════════════════════
# T2: Annotation quality gate (P0-CRITICAL)
# ═══════════════════════════════════════════════════════════════════════


def test_T2_pass_rate_gate_passes_when_threshold_met() -> None:
    """Gate does not abort when pass_rate >= threshold.

    Given: 80 cells, 64 with cell_type='T cell', 16 with 'Unknown' (80% pass).
    When:  _warn_if_low_coverage is called.
    Then:  No SystemExit is raised.
    """
    n = 80
    adata = AnnData(np.zeros((n, 10)))
    adata.obs["leiden"] = ["0"] * n  # original col (from 04_clustered.h5ad)
    adata.obs["cell_type"] = ["T cell"] * 64 + ["Unknown"] * 16
    original_cols = {"leiden"}

    cfg = MagicMock()
    cfg.marker = MagicMock()
    cfg.marker.quality_gate_min_pass_rate = 0.10
    log = MagicMock()

    # Should not raise — pass_rate=0.8 >= 0.10
    _warn_if_low_coverage(adata, cfg, log, original_cols)


def test_T2_pass_rate_gate_fires_at_zero_pass_rate(tmp_path) -> None:
    """Gate aborts when 100% of cells are Unknown and persists ONLY the new annotation cols.

    Given: adata with an original col (leiden) plus step-05 new cols
           (cell_type, cell_state), 0% pass rate.
    When:  _warn_if_low_coverage is called.
    Then:  SystemExit is raised and the lightweight writer receives a DataFrame
           with exactly the step-05 new columns — never the original ones.
    """
    n = 30
    adata = AnnData(np.zeros((n, 10)))
    adata.obs["leiden"] = ["0"] * n
    adata.obs["cell_type"] = ["Unknown"] * n
    adata.obs["cell_state"] = ["na"] * n
    original_cols = {"leiden"}

    annotated = tmp_path / "05_annotated.h5ad"
    annotated.touch()  # annotated exists → copy2 skipped, only append is exercised

    cfg = MagicMock()
    cfg.marker = MagicMock()
    cfg.marker.quality_gate_min_pass_rate = 0.10
    cfg.annotated_h5ad = str(annotated)
    cfg.cluster_h5ad = str(tmp_path / "04_clustered.h5ad")
    log = MagicMock()

    with (
        patch("core.utils.write_obs_columns_lightweight") as mock_write,
        pytest.raises(SystemExit),
    ):
        _warn_if_low_coverage(adata, cfg, log, original_cols)

    obs_df = mock_write.call_args.args[1]
    assert list(obs_df.columns) == ["cell_type", "cell_state"], (
        f"abort path must write only step-05 new columns, got {list(obs_df.columns)}"
    )
    assert mock_write.call_args.args[0] == str(annotated)


def test_T2_score_genes_path_gate_coverage_when_std_none(tmp_path) -> None:
    """Gate fires in the score_genes path even when std is None.

    Given: score_genes path (std=None), all cells annotated as Unknown.
    When:  _warn_if_low_coverage is called.
    Then:  SystemExit is raised regardless of standardizer availability,
           and only the new annotation column is passed to the writer.
    """
    adata = AnnData(np.zeros((30, 10)))
    adata.obs["leiden"] = ["0"] * 30
    adata.obs["cell_type"] = ["Unknown"] * 30
    original_cols = {"leiden"}

    annotated = tmp_path / "05_annotated.h5ad"
    annotated.touch()

    cfg = MagicMock()
    cfg.marker = MagicMock()
    cfg.marker.quality_gate_min_pass_rate = 0.10
    cfg.annotated_h5ad = str(annotated)
    cfg.cluster_h5ad = str(tmp_path / "04_clustered.h5ad")
    log = MagicMock()

    with (
        patch("core.utils.write_obs_columns_lightweight") as mock_write,
        pytest.raises(SystemExit),
    ):
        _warn_if_low_coverage(adata, cfg, log, original_cols)

    obs_df = mock_write.call_args.args[1]
    assert list(obs_df.columns) == ["cell_type"]


def test_append_path_writes_only_new_columns(tmp_path) -> None:
    """_write_lightweight appends ONLY step-05 new obs cols to the copied 04_clustered.

    Given: real 04_clustered.h5ad with original obs cols (leiden categorical,
           n_genes int, pct_mito float); 05_annotated already copied from it
           (copy+append mode); adata loaded from cluster + annotation cols added.
    When:  _write_lightweight(adata, cfg, log, original_cols) runs.
    Then:  read-back file has original cols untouched (same values/dtypes) and
           exactly the new annotation cols appended — no extras, no missing.
    """
    rng = np.random.RandomState(7)
    n = 20
    src = AnnData(rng.randn(n, 5).astype(np.float32))
    src.obs["leiden"] = pd.Categorical(["0", "1"] * (n // 2))
    src.obs["n_genes"] = rng.randint(100, 500, n)
    src.obs["pct_mito"] = rng.rand(n)

    cluster_path = tmp_path / "04_clustered.h5ad"
    annotated_path = tmp_path / "05_annotated.h5ad"
    src.write_h5ad(cluster_path)
    shutil.copy2(cluster_path, annotated_path)  # 05 starts as copy of 04

    adata = sc.read(cluster_path)
    original_cols = set(adata.obs.columns)
    adata.obs["cell_type"] = pd.Categorical(["T cell", "B cell"] * (n // 2))
    adata.obs["cell_state"] = ["active", "resting"] * (n // 2)
    adata.obs["annot_confidence"] = rng.rand(n)

    cfg = MagicMock()
    cfg.annotated_h5ad = str(annotated_path)
    cfg.cluster_h5ad = str(cluster_path)
    log = MagicMock()
    _mod._write_lightweight(adata, cfg, log, original_cols)

    out = sc.read(annotated_path)

    # Exactly the new annotation columns were appended — no extras, no missing
    expected_cols = original_cols | {"cell_type", "cell_state", "annot_confidence"}
    assert set(out.obs.columns) == expected_cols, (
        f"obs columns after append = {list(out.obs.columns)}; "
        f"expected exactly {sorted(expected_cols)}"
    )
    # Original columns preserved — same values, same dtypes (not rewritten)
    for col in ("leiden", "n_genes", "pct_mito"):
        assert out.obs[col].equals(src.obs[col]), f"original column {col!r} was altered"
        assert out.obs[col].dtype == src.obs[col].dtype, (
            f"original column {col!r} dtype changed: {out.obs[col].dtype} != {src.obs[col].dtype}"
        )
    # New columns appended with correct values and dtypes
    assert list(out.obs["cell_type"]) == list(adata.obs["cell_type"])
    assert isinstance(out.obs["cell_type"].dtype, pd.CategoricalDtype)
    assert list(out.obs["cell_state"]) == list(adata.obs["cell_state"])
    np.testing.assert_allclose(
        out.obs["annot_confidence"].values, adata.obs["annot_confidence"].values
    )
