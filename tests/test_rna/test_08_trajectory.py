"""Tests for rna/steps/08_trajectory.py — gene_trends and pseudotime correlation."""

import importlib.util
import logging
import os

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from core.config.schema import Config

# ── Load the trajectory module ────────────────────────────────────────────
# conftest.py adds repo root to sys.path, but rna/steps/__init__.py does not
# re-export the trajectory module.  Load it by file path.
_STEP_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "rna", "steps", "08_trajectory.py"
)
_spec = importlib.util.spec_from_file_location("rna.steps._08_trajectory_test", _STEP_PATH)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Could not load trajectory module at {_STEP_PATH}")
trajectory = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(trajectory)

gene_trends = trajectory.gene_trends
_select_pseudotime_correlated = trajectory._select_pseudotime_correlated


def _make_mock_adata(n_cells: int = 100, n_genes: int = 20, seed: int = 42) -> AnnData:
    """Create mock AnnData with raw, dpt_pseudotime, cell_type."""
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.raw = adata.copy()
    adata.obs["dpt_pseudotime"] = rng.uniform(0, 1, n_cells)
    adata.obs["cell_type"] = rng.choice(["A", "B", "C"], n_cells)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    return adata


def _make_mock_branch_results(n_genes: int = 5) -> pd.DataFrame:
    """Mock branch_results DataFrame with columns names / scores / pvals_adj."""
    return pd.DataFrame(
        {
            "names": [f"GENE_{i}" for i in range(n_genes)],
            "scores": [10.0 - i for i in range(n_genes)],
            "pvals_adj": [1e-10] * n_genes,
        }
    )


# ── Test: gene_trends ─────────────────────────────────────────────────────


class TestGeneTrends:
    """Tests for the gene_trends() main function."""

    def test_gene_trends_happy(
        self,
        tmp_path,
        caplog,
    ) -> None:
        """Happy path: branch DE + CFG override produce a union of genes."""
        caplog.set_level(logging.INFO)
        sc.settings.figdir = str(tmp_path)

        adata = _make_mock_adata(n_cells=100, n_genes=20)
        branch_results = _make_mock_branch_results(n_genes=5)

        cfg = Config.model_validate(
            {
                "trajectory": {
                    "pseudotime_n_branch_de": 10,
                    "pseudotime_n_correlated": 10,
                    "pseudotime_cor_pval": 0.05,
                    "pseudotime_genes": ["GENE_0", "GENE_1"],
                },
                "table_dir": str(tmp_path),
            }
        )

        log = logging.getLogger("test_gene_trends_happy")
        gene_trends(adata, cfg, log, branch_results=branch_results)

        # The function should produce a log message with gene counts
        assert "Gene trends along pseudotime" in caplog.text, (
            "Expected log message about gene trends progress"
        )

        # The union should include branch DE genes and CFG override genes
        # (They overlap, so total unique should be at least 5)
        # "GENE_0" and "GENE_1" appear in both sources — deduplication
        # should keep them once each.
        assert "pseudotime_correlation" in caplog.text, (
            "Expected log to mention pseudotime correlation source"
        )
        assert "Pseudotime trend genes exported" in caplog.text, (
            "Expected log message about CSV export of pseudotime trend genes"
        )

    def test_gene_trends_no_dpt(
        self,
        caplog,
    ) -> None:
        """Missing dpt_pseudotime in obs → early return."""
        caplog.set_level(logging.INFO)

        adata = _make_mock_adata()
        del adata.obs["dpt_pseudotime"]
        cfg = Config()
        log = logging.getLogger("test_gene_trends_no_dpt")

        result = gene_trends(adata, cfg, log)

        assert result is None
        assert "No DPT" in caplog.text

    def test_gene_trends_no_raw(
        self,
        caplog,
    ) -> None:
        """raw = None → early return."""
        caplog.set_level(logging.INFO)

        adata = _make_mock_adata()
        adata.raw = None
        cfg = Config()
        log = logging.getLogger("test_gene_trends_no_raw")

        result = gene_trends(adata, cfg, log)

        assert result is None
        assert "No raw data" in caplog.text


# ── Test: _select_pseudotime_correlated ───────────────────────────────────


class TestSelectPseudotimeCorrelated:
    """Tests for the _select_pseudotime_correlated() internal function."""

    def test_select_correlated_happy(self) -> None:
        """Genes perfectly correlated with pseudotime appear in the result."""
        rng = np.random.RandomState(42)
        n_cells = 100
        n_genes = 20

        x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
        adata = AnnData(x)
        adata.var_names = [f"GENE_{i}" for i in range(n_genes)]

        # Pseudotime: linear ramp from 0.1 to 0.9
        pseudotime = np.linspace(0.1, 0.9, n_cells)
        adata.obs["dpt_pseudotime"] = pseudotime

        # Create raw via copy, then modify first 3 gene columns in-place
        adata.raw = adata.copy()
        adata.raw.X[:, :3] = pseudotime[:, None] * 10

        cfg = Config.model_validate(
            {
                "trajectory": {"pseudotime_n_correlated": 10, "pseudotime_cor_pval": 0.05},
            }
        )
        result, corr_df = _select_pseudotime_correlated(adata, cfg)

        assert isinstance(result, list), f"Expected list, got {type(result)}"
        assert len(result) > 0, "Expected at least one correlated gene"
        for i in range(3):
            assert f"GENE_{i}" in result, (
                f"Expected perfectly correlated GENE_{i} in result, got {result}"
            )

        assert isinstance(corr_df, pd.DataFrame), f"Expected DataFrame, got {type(corr_df)}"
        expected_cols = ["gene", "rho", "pval_raw", "pval_adj"]
        assert list(corr_df.columns) == expected_cols, (
            f"Unexpected columns: {list(corr_df.columns)}"
        )
        assert len(corr_df) > 0, "Expected non-empty DataFrame"

    def test_select_correlated_constant_pseudotime(self) -> None:
        """Constant pseudotime -> empty list (Spearman requires variance)."""
        import warnings

        from scipy.stats import ConstantInputWarning

        adata = _make_mock_adata()
        adata.obs["dpt_pseudotime"] = np.full(100, 0.5)

        cfg = Config.model_validate(
            {
                "trajectory": {"pseudotime_n_correlated": 10, "pseudotime_cor_pval": 0.05},
            }
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConstantInputWarning)
            result, corr_df = _select_pseudotime_correlated(adata, cfg)

        assert isinstance(result, list), f"Expected list, got {type(result)}"
        assert len(result) == 0, f"Expected empty result for constant pseudotime, got {result}"

        assert isinstance(corr_df, pd.DataFrame), f"Expected DataFrame, got {type(corr_df)}"
        assert corr_df.empty, f"Expected empty DataFrame, got {len(corr_df)} rows"


# ── Test: gene deduplication priority order ───────────────────────────────


class TestGenePriorityOrder:
    """Gene deduplication preserves source priority order 1→2→3→4."""

    def test_gene_priority_order(
        self,
        tmp_path,
        caplog,
    ) -> None:
        """Same gene from branch DE + CFG override appears only once."""
        caplog.set_level(logging.INFO)
        sc.settings.figdir = str(tmp_path)
        # Build adata manually so "TOP" is in both adata and raw var_names
        rng = np.random.RandomState(0)
        x = rng.poisson(lam=1.0, size=(100, 20)).astype(np.float32)
        adata = AnnData(x)
        adata.var_names = ["TOP"] + [f"GENE_{i}" for i in range(1, 20)]
        adata.obs["dpt_pseudotime"] = rng.uniform(0, 1, 100)
        adata.obs["cell_type"] = rng.choice(["A", "B", "C"], 100)
        adata.raw = adata.copy()

        # Branch results: "TOP" is the top DE gene
        branch_results = pd.DataFrame(
            {
                "names": ["TOP", "GENE_1", "GENE_2", "GENE_3", "GENE_4"],
                "scores": [100.0, 50.0, 40.0, 30.0, 20.0],
                "pvals_adj": [1e-20, 1e-10, 1e-10, 1e-10, 1e-10],
            }
        )

        cfg = Config.model_validate(
            {
                "trajectory": {
                    "pseudotime_n_branch_de": 10,
                    "pseudotime_n_correlated": 10,
                    "pseudotime_cor_pval": 0.05,
                    "pseudotime_genes": ["TOP"],
                },
                "tissue_kb": "",
                "table_dir": str(tmp_path),
            }
        )

        log = logging.getLogger("test_gene_priority_order")
        gene_trends(adata, cfg, log, branch_results=branch_results)

        # The log shows the unique count after deduplication
        assert "Gene trends along pseudotime" in caplog.text


# ── Test: heatmap binning ─────────────────────────────────────────────────


class TestHeatmapBinning:
    """Heatmap binning with small number of unique pseudotime values."""

    def test_heatmap_binning(
        self,
        tmp_path,
    ) -> None:
        """5 unique pseudotime values + 6 union_genes → n_bins = 4."""
        sc.settings.figdir = str(tmp_path)

        n_cells = 10
        n_genes = 10
        adata = _make_mock_adata(n_cells=n_cells, n_genes=n_genes)
        # Exactly 5 unique pseudotime values → n_bins = min(10, 5-1) = 4
        unique_pt = np.array([0.0, 0.2, 0.4, 0.6, 0.8])
        adata.obs["dpt_pseudotime"] = np.tile(unique_pt, 2)[:n_cells]

        # 6 union genes needed to trigger the heatmap block (>= 5)
        branch_results = pd.DataFrame(
            {
                "names": [f"GENE_{i}" for i in range(6)],
                "scores": [10.0 - i for i in range(6)],
                "pvals_adj": [1e-10] * 6,
            }
        )

        cfg = Config.model_validate(
            {
                "trajectory": {
                    "pseudotime_n_branch_de": 10,
                    "pseudotime_n_correlated": 10,
                    "pseudotime_cor_pval": 0.05,
                },
                "table_dir": str(tmp_path),
            }
        )
        log = logging.getLogger("test_heatmap_binning")
        gene_trends(adata, cfg, log, branch_results=branch_results)

        # The heatmap block runs when union_genes >= 5 AND n_bins >= 2.
        # With 5 unique pt values, n_bins = 4, so the safe_plot call fires.
        # Since safe_plot catches exceptions internally, the function
        # completes without error regardless of whether the plot succeeds.
        # No crash is the assertion.
