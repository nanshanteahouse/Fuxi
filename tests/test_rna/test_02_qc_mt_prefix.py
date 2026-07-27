"""Tests for MT prefix auto-detection in rna/steps/02_qc.py.

P2-2: Extend MT prefix auto-detection from binary MT-↔mt- to
candidate list scan (MT-, mt-, Mt-).

Tests
-----
- test_auto_detect_uppercase     — MT-CO1 in var → detects MT-
- test_auto_detect_lowercase     — mt-Co1 in var  → detects mt-
- test_auto_detect_mixed_case    — Mt-Co1 in var  → detects Mt-
- test_no_mt_genes_warning       — no MT genes   → log.warning called
- test_mt_gene_list_takes_priority — no prefix match but mt_gene_list set
"""

from __future__ import annotations

import importlib
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
from anndata import AnnData

# ── Ensure repo root is on sys.path ──────────────────────────────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 02_qc module via file path ──────────────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "02_qc.py")
_spec = importlib.util.spec_from_file_location(
    "rna.steps._02_qc_test",
    _STEP_PATH,
)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═══════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════


def _make_cfg(
    mt_gene_pattern: str = "MT-",
    mt_gene_list: list[str] | None = None,
) -> MagicMock:
    """Create a Config mock with QC settings."""
    cfg = MagicMock()
    cfg.qc = MagicMock()
    cfg.qc.mt_gene_pattern = mt_gene_pattern
    cfg.qc.mt_gene_list = mt_gene_list
    return cfg


def _make_adata(
    gene_names: list[str] | None = None,
    n_cells: int = 5,
    n_genes: int = 10,
    seed: int = 42,
) -> AnnData:
    """Create a minimal AnnData."""
    if gene_names is not None:
        n_genes = len(gene_names)
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    if gene_names is not None:
        adata.var_names = gene_names
    else:
        adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    return adata


def _mock_calc_qc_metrics(adata, **kwargs):
    """Side effect for sc.pp.calculate_qc_metrics mock: set required obs columns."""
    n = adata.n_obs
    adata.obs["n_genes_by_counts"] = np.full(n, adata.n_vars)
    adata.obs["total_counts"] = np.full(n, 100.0)
    adata.obs["pct_counts_mt"] = np.zeros(n)
    adata.obs["pct_counts_ribo"] = np.zeros(n)
    adata.obs["log_genes_per_umi"] = np.full(n, 0.9)


# ═══════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════


class TestMtPrefixAutoDetect:
    """Auto-detect switches mt_gene_pattern when the default does not match."""

    def test_auto_detect_uppercase(self):
        """var_names containing MT-CO1 → detects MT- prefix."""
        adata = _make_adata(gene_names=["MT-CO1", "GENE_1", "GENE_2"])
        cfg = _make_cfg(mt_gene_pattern="MT-")
        log = MagicMock()

        # For initial check, "MT-" is the pattern and it does match MT-CO1
        # → no switch needed; verify pattern stays "MT-"
        with patch.object(_mod.sc.pp, "calculate_qc_metrics", side_effect=_mock_calc_qc_metrics):
            _mod.compute_qc_metrics(adata, cfg, log)

        assert cfg.qc.mt_gene_pattern == "MT-"
        assert adata.var["mt"].tolist() == [True, False, False]
        # Verify no auto-switch message was logged
        auto_switch_calls = [c for c in log.info.call_args_list if "Auto-switched" in str(c)]
        assert len(auto_switch_calls) == 0

    def test_auto_detect_lowercase(self):
        """var_names containing mt-Co1 → detects mt- prefix."""
        adata = _make_adata(gene_names=["mt-Co1", "GENE_1", "GENE_2"])
        cfg = _make_cfg(mt_gene_pattern="MT-")
        log = MagicMock()

        with patch.object(_mod.sc.pp, "calculate_qc_metrics", side_effect=_mock_calc_qc_metrics):
            _mod.compute_qc_metrics(adata, cfg, log)

        assert cfg.qc.mt_gene_pattern == "mt-"
        assert adata.var["mt"].tolist() == [True, False, False]
        log.info.assert_any_call("Auto-switched mt_gene_pattern: '%s' -> '%s'", "MT-", "mt-")

    def test_auto_detect_mixed_case(self):
        """var_names containing Mt-Co1 → detects Mt- prefix."""
        adata = _make_adata(gene_names=["Mt-Co1", "GENE_1", "GENE_2"])
        cfg = _make_cfg(mt_gene_pattern="MT-")
        log = MagicMock()

        with patch.object(_mod.sc.pp, "calculate_qc_metrics", side_effect=_mock_calc_qc_metrics):
            _mod.compute_qc_metrics(adata, cfg, log)

        assert cfg.qc.mt_gene_pattern == "Mt-"
        assert adata.var["mt"].tolist() == [True, False, False]
        log.info.assert_any_call("Auto-switched mt_gene_pattern: '%s' -> '%s'", "MT-", "Mt-")

    def test_no_mt_genes_warning(self):
        """No MT genes in var_names → log.warning called."""
        adata = _make_adata(gene_names=["GENE_A", "GENE_B", "GENE_C"])
        cfg = _make_cfg(mt_gene_pattern="MT-", mt_gene_list=[])
        log = MagicMock()

        with patch.object(_mod.sc.pp, "calculate_qc_metrics", side_effect=_mock_calc_qc_metrics):
            _mod.compute_qc_metrics(adata, cfg, log)

        # Pattern stays default since no candidate matched
        assert cfg.qc.mt_gene_pattern == "MT-"
        assert adata.var["mt"].tolist() == [False, False, False]
        # Should have called log.warning for no MT genes detected
        log.warning.assert_called_once()

    def test_mt_gene_list_takes_priority(self):
        """No prefix match but mt_gene_list set → no warning, mask via isin."""
        # var_names have no MT- / mt- / Mt- prefix
        adata = _make_adata(gene_names=["CO1", "GENE_1", "GENE_2"])
        # mt_gene_list explicitly lists the MT gene by exact name
        cfg = _make_cfg(mt_gene_pattern="MT-", mt_gene_list=["CO1"])
        log = MagicMock()

        with patch.object(_mod.sc.pp, "calculate_qc_metrics", side_effect=_mock_calc_qc_metrics):
            _mod.compute_qc_metrics(adata, cfg, log)

        # Pattern unchanged (no candidate matched)
        assert cfg.qc.mt_gene_pattern == "MT-"
        # mt_mask includes CO1 via isin (not via prefix)
        assert adata.var["mt"].tolist() == [True, False, False]
        # mt_gene_list non-empty → no warning
        log.warning.assert_not_called()
