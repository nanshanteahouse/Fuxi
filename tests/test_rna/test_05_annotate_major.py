"""Tests for rna/steps/05_annotate_major.py — use_raw fix in ai_annotate.

T1 (P0-CRITICAL) from cross-batch-critical-fixes plan:
  ``ai_annotate()`` must call ``rank_genes_groups`` with ``use_raw=True``
  when ``adata.raw`` exists.
"""

import importlib.util
import json
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
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
    adata.obs["cell_type"] = ["T cell"] * 64 + ["Unknown"] * 16

    cfg = MagicMock()
    cfg.marker = MagicMock()
    cfg.marker.quality_gate_min_pass_rate = 0.10
    log = MagicMock()

    # Should not raise — pass_rate=0.8 >= 0.10
    _warn_if_low_coverage(adata, cfg, log)


def test_T2_pass_rate_gate_fires_at_zero_pass_rate() -> None:
    """Gate aborts when 100% of cells are Unknown.

    Given: 30 cells, all with cell_type='Unknown' (0% pass).
    When:  _warn_if_low_coverage is called.
    Then:  SystemExit is raised and safe_write is called before exit.
    """
    n = 30
    adata = AnnData(np.zeros((n, 10)))
    adata.obs["cell_type"] = ["Unknown"] * n

    cfg = MagicMock()
    cfg.marker = MagicMock()
    cfg.marker.quality_gate_min_pass_rate = 0.10
    cfg.annotated_h5ad = "/tmp/test.h5ad"
    log = MagicMock()

    with patch.object(_mod, "safe_write") as mock_write, pytest.raises(SystemExit):
        _warn_if_low_coverage(adata, cfg, log)
    mock_write.assert_called_once_with(adata, "/tmp/test.h5ad", cfg=cfg)


def test_T2_score_genes_path_gate_coverage_when_std_none() -> None:
    """Gate fires in the score_genes path even when std is None.

    Given: score_genes path (std=None), all cells annotated as Unknown.
    When:  _warn_if_low_coverage is called.
    Then:  SystemExit is raised regardless of standardizer availability.
    """
    adata = AnnData(np.zeros((30, 10)))
    adata.obs["cell_type"] = ["Unknown"] * 30

    cfg = MagicMock()
    cfg.marker = MagicMock()
    cfg.marker.quality_gate_min_pass_rate = 0.10
    cfg.annotated_h5ad = "/tmp/test.h5ad"
    log = MagicMock()

    with patch.object(_mod, "safe_write"), pytest.raises(SystemExit):
        _warn_if_low_coverage(adata, cfg, log)
