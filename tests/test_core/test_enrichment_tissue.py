"""Tests for core/enrichment_tissue.py — tissue-aware enrichment post-processor.

Functions under test
--------------------
- ``core.enrichment_tissue.compute_pathway_relevance``  — KB-overlap statistics
- ``core.enrichment_tissue.cluster_redundant_pathways``  — Jaccard-based greedy clustering
- ``core.enrichment_tissue.filter_enrichment_by_tissue`` — whitelist/blacklist/KB-score filtering
"""

from __future__ import annotations

import pandas as pd
import pytest

from core.pipeline.enrichment import (
    cluster_redundant_pathways,
    compute_pathway_relevance,
    filter_enrichment_by_tissue,
)

# ======================================================================
#  compute_pathway_relevance tests
# ======================================================================


class TestComputePathwayRelevance:
    """``compute_pathway_relevance`` — KB-overlap annotation and scoring."""

    # ------------------------------------------------------------------
    #  Normal case
    # ------------------------------------------------------------------

    def test_compute_pathway_relevance(self) -> None:
        """KB markers correctly annotate overlap count, ratio, genes, and score."""
        df = pd.DataFrame(
            {
                "Term": [
                    "Visual Perception",
                    "Ribosome",
                    "Phototransduction",
                ],
                "Adjusted P-value": [0.001, 0.01, 0.0001],
                "Overlap": [
                    "RHO,GNAT1,ACTB",
                    "RPL5,RPL10,RPS3",
                    "RHO,GNAT1,PDE6B,OPN1SW",
                ],
            }
        )
        kb_markers = {"RHO", "GNAT1", "PDE6B", "OPN1SW", "NRL", "VSX2"}
        result = compute_pathway_relevance(df, kb_markers)

        # -- Column presence --
        for col in (
            "kb_overlap_genes",
            "kb_overlap_count",
            "kb_overlap_ratio",
            "kb_relevance_score",
        ):
            assert col in result.columns, f"Missing column: {col}"

        # -- kb_overlap_count --
        assert list(result["kb_overlap_count"]) == [2, 0, 4]

        # -- kb_overlap_ratio --
        ratios = list(result["kb_overlap_ratio"])
        assert ratios[0] == pytest.approx(2 / 3)
        assert ratios[1] == pytest.approx(0.0)
        assert ratios[2] == pytest.approx(1.0)

        # -- kb_overlap_genes --
        assert result["kb_overlap_genes"].iloc[0] == ["RHO", "GNAT1"]
        assert result["kb_overlap_genes"].iloc[1] == []
        assert result["kb_overlap_genes"].iloc[2] == ["RHO", "GNAT1", "PDE6B", "OPN1SW"]

        # -- kb_relevance_score --
        # Visual Perception: -log10(0.001)*(1+2/3) = 3.0 * 1.666... = 5.0
        assert result["kb_relevance_score"].iloc[0] == pytest.approx(5.0, rel=1e-5)
        # Ribosome: -log10(0.01)*(1+0) = 2.0
        assert result["kb_relevance_score"].iloc[1] == pytest.approx(2.0, rel=1e-5)
        # Phototransduction: -log10(0.0001)*(1+1) = 4.0 * 2.0 = 8.0
        assert result["kb_relevance_score"].iloc[2] == pytest.approx(8.0, rel=1e-5)

    # ------------------------------------------------------------------
    #  Empty / None KB
    # ------------------------------------------------------------------

    def test_compute_pathway_relevance_empty_kb(self) -> None:
        """With None or empty KB, scores degrade to -log10(p-value) and overlap columns are 0."""
        df = pd.DataFrame(
            {
                "Term": [
                    "Visual Perception",
                    "Ribosome",
                    "Phototransduction",
                ],
                "Adjusted P-value": [0.001, 0.01, 0.0001],
                "Overlap": [
                    "RHO,GNAT1,ACTB",
                    "RPL5,RPL10,RPS3",
                    "RHO,GNAT1,PDE6B,OPN1SW",
                ],
            }
        )

        # Test with None
        result_none = compute_pathway_relevance(df, kb_markers=None)
        assert list(result_none["kb_overlap_count"]) == [0, 0, 0]
        assert list(result_none["kb_overlap_ratio"]) == [0.0, 0.0, 0.0]
        assert result_none["kb_relevance_score"].iloc[0] == pytest.approx(3.0, rel=1e-5)
        assert result_none["kb_relevance_score"].iloc[1] == pytest.approx(2.0, rel=1e-5)
        assert result_none["kb_relevance_score"].iloc[2] == pytest.approx(4.0, rel=1e-5)
        for genes in result_none["kb_overlap_genes"]:
            assert genes == []

        # Test with empty set (should behave identically)
        result_empty = compute_pathway_relevance(df, kb_markers=set())
        assert list(result_empty["kb_overlap_count"]) == [0, 0, 0]
        for genes in result_empty["kb_overlap_genes"]:
            assert genes == []

    # ------------------------------------------------------------------
    #  Missing columns
    # ------------------------------------------------------------------

    def test_compute_pathway_relevance_missing_overlap(self) -> None:
        """Missing 'Overlap' column: all KB overlap columns are 0 / empty."""
        df = pd.DataFrame(
            {
                "Term": ["PathwayA", "PathwayB"],
                "Adjusted P-value": [0.01, 0.001],
            }
        )
        result = compute_pathway_relevance(df, kb_markers={"RHO", "NRL"})

        assert list(result["kb_overlap_count"]) == [0, 0]
        assert list(result["kb_overlap_ratio"]) == [0.0, 0.0]
        for genes in result["kb_overlap_genes"]:
            assert genes == []
        # Score degrades to -log10(p-value)
        assert result["kb_relevance_score"].iloc[0] == pytest.approx(2.0, rel=1e-5)
        assert result["kb_relevance_score"].iloc[1] == pytest.approx(3.0, rel=1e-5)

    def test_compute_pathway_relevance_fallback_pvalue(self) -> None:
        """When 'Adjusted P-value' is absent, falls back to 'P-value'."""
        df = pd.DataFrame(
            {
                "Term": ["PathwayA", "PathwayB"],
                "P-value": [0.01, 0.001],
                "Overlap": ["GENE1,GENE2", "GENE3"],
            }
        )
        result = compute_pathway_relevance(df, kb_markers=None)
        # Should use P-value as fallback
        assert result["kb_relevance_score"].iloc[0] == pytest.approx(2.0, rel=1e-5)
        assert result["kb_relevance_score"].iloc[1] == pytest.approx(3.0, rel=1e-5)


# ======================================================================
#  cluster_redundant_pathways tests
# ======================================================================


class TestClusterRedundantPathways:
    """``cluster_redundant_pathways`` — greedy Jaccard-based clustering."""

    # ------------------------------------------------------------------
    #  Redundant pathways
    # ------------------------------------------------------------------

    def test_cluster_redundant_pathways(self) -> None:
        """Redundant pathways (Jaccard >= 0.6) are clustered; representatives selected."""
        df = pd.DataFrame(
            {
                "Term": [
                    "Visual Perception",
                    "Retina Development",
                    "Phototransduction Cascade",
                    "Ribosome Assembly",
                ],
                "Adjusted P-value": [0.01, 0.02, 0.03, 0.04],
                "Overlap": [
                    "A,B,C,D",
                    "A,B,E,F",
                    "A,B,C,D",
                    "X,Y,Z",
                ],
            }
        )
        result = cluster_redundant_pathways(df, similarity_threshold=0.6)

        # -- Column presence --
        assert "redundant_cluster_id" in result.columns
        assert "redundant_representative" in result.columns

        # Visual Perception (most significant, cluster 0, rep)
        # Phototransduction Cascade (Jaccard=1.0 with Visual → cluster 0, not rep)
        # Retina Development (Jaccard=0.33 with Visual → cluster 1, rep)
        # Ribosome Assembly (Jaccard=0 with all → cluster 2, rep)
        assert list(result["redundant_cluster_id"]) == [0, 1, 0, 2]
        assert list(result["redundant_representative"]) == [True, True, False, True]

    # ------------------------------------------------------------------
    #  No overlap
    # ------------------------------------------------------------------

    def test_cluster_redundant_pathways_no_overlap(self) -> None:
        """Disjoint gene sets yield separate clusters, all representatives."""
        df = pd.DataFrame(
            {
                "Term": ["PathA", "PathB", "PathC"],
                "Adjusted P-value": [0.01, 0.02, 0.03],
                "Overlap": ["A,B", "C,D", "E,F"],
            }
        )
        result = cluster_redundant_pathways(df, similarity_threshold=0.6)

        assert list(result["redundant_cluster_id"]) == [0, 1, 2]
        assert list(result["redundant_representative"]) == [True, True, True]

    # ------------------------------------------------------------------
    #  Single row — edge case
    # ------------------------------------------------------------------

    def test_cluster_redundant_pathways_single_row(self) -> None:
        """Single-row DataFrame: cluster 0, representative True."""
        df = pd.DataFrame(
            {
                "Term": ["Only Pathway"],
                "Adjusted P-value": [0.01],
                "Overlap": ["A,B,C"],
            }
        )
        result = cluster_redundant_pathways(df)
        assert list(result["redundant_cluster_id"]) == [0]
        assert list(result["redundant_representative"]) == [True]

    # ------------------------------------------------------------------
    #  Empty DataFrame — edge case
    # ------------------------------------------------------------------

    def test_cluster_redundant_pathways_empty(self) -> None:
        """Empty DataFrame returns an empty copy (no error)."""
        df = pd.DataFrame(columns=["Term", "Adjusted P-value", "Overlap"])
        result = cluster_redundant_pathways(df)
        assert len(result) == 0


# ======================================================================
#  filter_enrichment_by_tissue tests
# ======================================================================


class TestFilterEnrichmentByTissue:
    """``filter_enrichment_by_tissue`` — soft/hard/off modes, WL/BL priority."""

    # ------------------------------------------------------------------
    #  Shared fixture
    # ------------------------------------------------------------------

    @pytest.fixture
    def kb_df(self) -> pd.DataFrame:
        """DataFrame with kb_relevance_score for soft/hard mode tests."""
        return pd.DataFrame(
            {
                "Term": [
                    "Phototransduction",
                    "Ribosome",
                    "Neuron Projection",
                    "Cell Cycle",
                ],
                "kb_relevance_score": [0.5, 2.0, 3.0, 0.1],
            }
        )

    # ------------------------------------------------------------------
    #  Off mode
    # ------------------------------------------------------------------

    def test_filter_enrichment_by_tissue_off(self, kb_df: pd.DataFrame) -> None:
        """'off' mode returns the original DataFrame object unchanged."""
        result = filter_enrichment_by_tissue(kb_df, mode="off")
        assert result is kb_df  # same object, no copy
        assert list(result.columns) == ["Term", "kb_relevance_score"]

    # ------------------------------------------------------------------
    #  Soft mode
    # ------------------------------------------------------------------

    def test_filter_enrichment_by_tissue_soft(self, kb_df: pd.DataFrame) -> None:
        """'soft' mode annotates all rows; whitelist/blacklist/KB-score decisions correct."""
        result = filter_enrichment_by_tissue(
            kb_df,
            mode="soft",
            pathway_whitelist=["Phototransduction"],
            pathway_blacklist=["Ribosome"],
        )

        # All rows kept
        assert len(result) == 4
        assert "tissue_relevant" in result.columns
        assert "tissue_relevance_score" in result.columns

        # Row 0: whitelist → True
        assert result["tissue_relevant"].iloc[0]
        assert result["tissue_relevance_score"].iloc[0] == 1.0

        # Row 1: blacklist → False
        assert not result["tissue_relevant"].iloc[1]
        assert result["tissue_relevance_score"].iloc[1] == -1.0

        # Unbiased rows: [3.0, 0.1], median = 1.55
        # Row 2 (3.0 > 1.55): True, normalised = (3.0-0.1)/(3.0-0.1) = 1.0
        assert result["tissue_relevant"].iloc[2]
        assert result["tissue_relevance_score"].iloc[2] == pytest.approx(1.0, rel=1e-5)

        # Row 3 (0.1 ≯ 1.55): False, score=0.0
        assert not result["tissue_relevant"].iloc[3]
        assert result["tissue_relevance_score"].iloc[3] == 0.0

        # Verify the copy semantic (not the same object)
        assert result is not kb_df

    # ------------------------------------------------------------------
    #  Hard mode
    # ------------------------------------------------------------------

    def test_filter_enrichment_by_tissue_hard(self, kb_df: pd.DataFrame) -> None:
        """'hard' mode filters to only tissue-relevant rows."""
        result = filter_enrichment_by_tissue(
            kb_df,
            mode="hard",
            pathway_whitelist=["Phototransduction"],
            pathway_blacklist=["Ribosome"],
        )

        # Only rows 0 (whitelist) and 2 (KB score > median) survive
        assert len(result) == 2
        assert list(result["Term"]) == ["Phototransduction", "Neuron Projection"]
        assert list(result["tissue_relevant"]) == [True, True]

    # ------------------------------------------------------------------
    #  Whitelist overrides blacklist
    # ------------------------------------------------------------------

    def test_filter_enrichment_by_tissue_whitelist_overrides(self) -> None:
        """When a term matches both whitelist and blacklist, whitelist wins (kept)."""
        df = pd.DataFrame(
            {
                "Term": ["Phototransduction", "Ribosome"],
                "kb_relevance_score": [0.5, 0.5],
            }
        )
        # "Phototransduction" appears in both lists
        result = filter_enrichment_by_tissue(
            df,
            mode="hard",
            pathway_whitelist=["Phototransduction"],
            pathway_blacklist=["Phototransduction", "Ribosome"],
        )
        # Whitelist priority: Phototransduction kept; Ribosome filtered
        assert len(result) == 1
        assert result["Term"].iloc[0] == "Phototransduction"
        assert result["tissue_relevant"].iloc[0]

    # ------------------------------------------------------------------
    #  Invalid mode raises
    # ------------------------------------------------------------------

    def test_filter_enrichment_by_tissue_invalid_mode(self) -> None:
        """An unknown mode raises ValueError."""
        df = pd.DataFrame({"Term": ["PathwayA"]})
        with pytest.raises(ValueError, match="Unknown mode"):
            filter_enrichment_by_tissue(df, mode="unknown")


# ======================================================================
#  Integration: retina pathway relevance
# ======================================================================


class TestIntegrationRetina:
    """``load_pathway_relevance(\"retina\")`` — verifies actual data structure."""

    def test_integration_retina_pathway_relevance(self) -> None:
        """Retina pathway relevance dict has expected keys and sizes."""
        from core.kb import load_pathway_relevance

        pr = load_pathway_relevance("retina")

        # Three expected keys
        assert "key_pathways" in pr
        assert "generic_pathways" in pr
        assert "kb_pathway_markers" in pr

        # Key pathways: photoreception, development, synapses, etc.
        assert len(pr["key_pathways"]) >= 15
        # Generic pathways: translation, metabolism, DNA maintenance
        assert len(pr["generic_pathways"]) >= 10
        # KB markers: pathway → gene list mappings
        assert len(pr["kb_pathway_markers"]) >= 5

        # Verify specific known entries
        assert "Phototransduction" in pr["key_pathways"]
        assert "Ribosome" in pr["generic_pathways"]
        assert "phototransduction" in pr["kb_pathway_markers"]
        assert "RHO" in pr["kb_pathway_markers"]["phototransduction"]
