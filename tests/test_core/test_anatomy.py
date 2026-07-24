"""Tests for core/anatomy.py — anatomical adjacency loading and CCI filtering.

Functions under test
--------------------
- ``core.anatomy.load_adjacency``     — primary entry point (custom file → tissue KB → empty)
- ``core.kb.load_adjacency`` — tissue knowledge-base lookup
- ``filter_cci_by_adjacency``         — annotate / filter CCI results against adjacency
"""

from __future__ import annotations

import pandas as pd
import pytest

from core.kb import load_adjacency as load_tissue_adj
from core.pipeline.anatomy import filter_cci_by_adjacency
from core.pipeline.anatomy import load_adjacency as load_adj

# ======================================================================
#  load_adjacency tests
# ======================================================================


class TestLoadAdjacency:
    """``load_adjacency`` — tissue KB and custom-file loading."""

    # ------------------------------------------------------------------
    #  Tissue KB: retina (known adjacency)
    # ------------------------------------------------------------------

    def test_load_adjacency_retina(self) -> None:
        """Retina adjacency returns ≥20 rows with correct columns & known pair."""
        df = load_tissue_adj("retina")
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["source", "target", "adjacency_type"]
        assert len(df) >= 20

        # Verify a specific known pair
        known = (df["source"] == "RPE") & (df["target"] == "Rod_Photoreceptor")
        assert known.any(), "Expected (RPE, Rod_Photoreceptor) in retina adjacency"
        assert (df.loc[known, "adjacency_type"] == "physical").all()

    # ------------------------------------------------------------------
    #  Tissue KB: unknown tissue → empty
    # ------------------------------------------------------------------

    def test_load_adjacency_unknown(self) -> None:
        """Unknown tissue returns an empty DataFrame — no exception."""
        df = load_tissue_adj("brain")
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["source", "target", "adjacency_type"]
        assert len(df) == 0

    # ------------------------------------------------------------------
    #  Custom CSV file
    # ------------------------------------------------------------------

    def test_load_adjacency_custom_file(self, tmp_path) -> None:
        """Custom CSV file is loaded and returned with correct shape."""
        csv_path = tmp_path / "adjacency.csv"
        csv_path.write_text(
            "source,target,adjacency_type\nA,X,synaptic\nB,Y,gap_junction\nC,Z,physical\n"
        )
        df = load_adj(tissue="", custom_file=str(csv_path))
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["source", "target", "adjacency_type"]
        assert len(df) == 3


# ======================================================================
#  filter_cci_by_adjacency tests
# ======================================================================


class TestFilterCCI:
    """``filter_cci_by_adjacency`` — three modes, bidirectionality, type filtering."""

    # ------------------------------------------------------------------
    #  Shared fixtures
    # ------------------------------------------------------------------

    @pytest.fixture
    def lr_res(self) -> pd.DataFrame:
        """Minimal LIANA result with 4 interactions."""
        return pd.DataFrame(
            {
                "source": ["A", "B", "A", "C"],
                "target": ["X", "Y", "Y", "Z"],
                "ligand": ["L1", "L2", "L3", "L4"],
            }
        )

    @pytest.fixture
    def adjacency(self) -> pd.DataFrame:
        """Adjacency table where A↔X and B↔Y are known contacts."""
        return pd.DataFrame(
            {
                "source": ["A", "B"],
                "target": ["X", "Y"],
                "adjacency_type": ["synaptic", "gap_junction"],
            }
        )

    # ------------------------------------------------------------------
    #  Hard mode
    # ------------------------------------------------------------------

    def test_filter_cci_hard(self, lr_res: pd.DataFrame, adjacency: pd.DataFrame) -> None:
        """Hard mode keeps only rows whose (source, target) exists in adjacency."""
        result = filter_cci_by_adjacency(lr_res, adjacency, mode="hard")
        assert len(result) == 2
        result = result.reset_index(drop=True)
        assert result.iloc[0]["source"] == "A"
        assert result.iloc[0]["target"] == "X"
        assert result.iloc[1]["source"] == "B"
        assert result.iloc[1]["target"] == "Y"

    # ------------------------------------------------------------------
    #  Soft mode
    # ------------------------------------------------------------------

    def test_filter_cci_soft(self, lr_res: pd.DataFrame, adjacency: pd.DataFrame) -> None:
        """Soft mode retains all rows and annotates with adjacent / adjacency_type."""
        result = filter_cci_by_adjacency(lr_res, adjacency, mode="soft")
        assert len(result) == 4
        assert "adjacent" in result.columns
        assert "adjacency_type" in result.columns
        # Rows 0–1 match; rows 2–3 do not
        assert list(result["adjacent"]) == [True, True, False, False]
        assert list(result["adjacency_type"]) == [
            "synaptic",
            "gap_junction",
            "",
            "",
        ]

    # ------------------------------------------------------------------
    #  Off mode
    # ------------------------------------------------------------------

    def test_filter_cci_off(self, lr_res: pd.DataFrame, adjacency: pd.DataFrame) -> None:
        """Off mode returns the identical object (no copy)."""
        result = filter_cci_by_adjacency(lr_res, adjacency, mode="off")
        assert result is lr_res

    # ------------------------------------------------------------------
    #  Bidirectional matching
    # ------------------------------------------------------------------

    def test_filter_cci_bidirectional(self) -> None:
        """A→X adjacency also matches X→A in CCI data."""
        adj = pd.DataFrame(
            {
                "source": ["A"],
                "target": ["X"],
                "adjacency_type": ["synaptic"],
            }
        )
        lr = pd.DataFrame(
            {
                "source": ["X"],
                "target": ["A"],
                "ligand": ["L1"],
            }
        )
        result = filter_cci_by_adjacency(lr, adj, mode="hard")
        assert len(result) == 1
        assert result.iloc[0]["source"] == "X"
        assert result.iloc[0]["target"] == "A"

    # ------------------------------------------------------------------
    #  adjacency_types filter
    # ------------------------------------------------------------------

    def test_filter_cci_adjacency_types(self) -> None:
        """adjacency_types restricts which adjacency edges are considered."""
        adj = pd.DataFrame(
            {
                "source": ["A", "B"],
                "target": ["X", "Y"],
                "adjacency_type": ["synaptic", "gap_junction"],
            }
        )
        lr = pd.DataFrame(
            {
                "source": ["A", "B"],
                "target": ["X", "Y"],
                "ligand": ["L1", "L2"],
            }
        )

        # Hard mode — only synaptic
        hard_result = filter_cci_by_adjacency(lr, adj, mode="hard", adjacency_types=["synaptic"])
        assert len(hard_result) == 1
        assert hard_result.iloc[0]["source"] == "A"
        assert hard_result.iloc[0]["target"] == "X"

        # Soft mode — only synaptic is adjacent
        soft_result = filter_cci_by_adjacency(lr, adj, mode="soft", adjacency_types=["synaptic"])
        assert len(soft_result) == 2
        assert list(soft_result["adjacent"]) == [True, False]
        assert list(soft_result["adjacency_type"]) == ["synaptic", ""]
