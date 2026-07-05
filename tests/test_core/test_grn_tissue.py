"""Tests for core/grn_tissue.py — TF activity annotation via KB marker overlap.

Function under test
--------------------
- ``core.grn_tissue.compute_tf_relevance``  — KB-overlap statistics and weighted activity
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from core.grn_tissue import compute_tf_relevance


# ======================================================================
#  compute_tf_relevance tests
# ======================================================================


class TestComputeTFRelevance:
    """``compute_tf_relevance`` — per-TF KB-overlap statistics and weighted activity."""

    # ------------------------------------------------------------------
    #  Normal overlap
    # ------------------------------------------------------------------

    def test_normal_overlap(self) -> None:
        """KB markers correctly annotated; overlap count/ratio per TF is correct."""
        rng = np.random.default_rng(42)
        activity_df = pd.DataFrame(
            rng.random((3, 3)),
            index=pd.Index(["CT0", "CT1", "CT2"]),
            columns=pd.Index(["TF_A", "TF_B", "TF_C"]),
        )
        net = pd.DataFrame({
            "source": ["TF_A", "TF_A", "TF_B", "TF_C"],
            "target": ["GENE_1", "GENE_2", "GENE_3", "GENE_4"],
        })
        kb_markers = {"GENE_1", "GENE_2", "GENE_4"}

        annotated_activity_df, tf_annotation = compute_tf_relevance(
            activity_df, net, kb_markers,
        )

        # -- Annotation table shape --
        assert len(tf_annotation) == 3

        # -- TF_A: 2/2 overlap --
        row_a = tf_annotation.loc[tf_annotation["tf"] == "TF_A"].iloc[0]
        assert row_a["kb_overlap_count"] == 2
        assert row_a["kb_overlap_ratio"] == pytest.approx(1.0)

        # -- TF_B: 0/1 overlap --
        row_b = tf_annotation.loc[tf_annotation["tf"] == "TF_B"].iloc[0]
        assert row_b["kb_overlap_count"] == 0
        assert row_b["kb_overlap_ratio"] == pytest.approx(0.0)

        # -- TF_C: 1/1 overlap --
        row_c = tf_annotation.loc[tf_annotation["tf"] == "TF_C"].iloc[0]
        assert row_c["kb_overlap_count"] == 1
        assert row_c["kb_overlap_ratio"] == pytest.approx(1.0)

        # -- Weighted activity preserves shape --
        assert annotated_activity_df.shape == (3, 3)

    # ------------------------------------------------------------------
    #  Empty / None KB fallback
    # ------------------------------------------------------------------

    def test_empty_kb_fallback(self) -> None:
        """With None KB, all overlap statistics are zero; weighted = abs(activity)."""
        rng = np.random.default_rng(42)
        activity_df = pd.DataFrame(
            rng.random((3, 3)),
            index=pd.Index(["CT0", "CT1", "CT2"]),
            columns=pd.Index(["TF_A", "TF_B", "TF_C"]),
        )
        net = pd.DataFrame({
            "source": ["TF_A", "TF_A", "TF_B", "TF_C"],
            "target": ["GENE_1", "GENE_2", "GENE_3", "GENE_4"],
        })

        annotated_activity_df, tf_annotation = compute_tf_relevance(
            activity_df, net, kb_markers=None,
        )

        assert (tf_annotation["kb_overlap_count"] == 0).all()
        assert (tf_annotation["kb_overlap_ratio"] == 0.0).all()
        # When ratio is 0, weighted = abs(activity) * 1.0
        pd.testing.assert_frame_equal(
            annotated_activity_df, activity_df.abs(),
        )

    # ------------------------------------------------------------------
    #  No overlap with KB
    # ------------------------------------------------------------------

    def test_no_overlap(self) -> None:
        """When kb_markers contains no matching genes, all overlap counts are 0."""
        rng = np.random.default_rng(42)
        activity_df = pd.DataFrame(
            rng.random((3, 3)),
            index=pd.Index(["CT0", "CT1", "CT2"]),
            columns=pd.Index(["TF_A", "TF_B", "TF_C"]),
        )
        net = pd.DataFrame({
            "source": ["TF_A", "TF_A", "TF_B", "TF_C"],
            "target": ["GENE_1", "GENE_2", "GENE_3", "GENE_4"],
        })
        kb_markers = {"NONEXISTENT_GENE"}

        _, tf_annotation = compute_tf_relevance(activity_df, net, kb_markers)

        assert (tf_annotation["kb_overlap_count"] == 0).all()

    # ------------------------------------------------------------------
    #  Empty regulon net
    # ------------------------------------------------------------------

    def test_empty_regulon_net(self) -> None:
        """An empty net yields zero overlap for all TFs."""
        rng = np.random.default_rng(42)
        activity_df = pd.DataFrame(
            rng.random((3, 2)),
            index=pd.Index(["CT0", "CT1", "CT2"]),
            columns=pd.Index(["TF_A", "TF_B"]),
        )
        net = pd.DataFrame(columns=pd.Index(["source", "target"]))
        kb_markers = {"GENE_1"}

        annotated_activity_df, tf_annotation = compute_tf_relevance(
            activity_df, net, kb_markers,
        )

        assert (tf_annotation["kb_overlap_count"] == 0).all()
        pd.testing.assert_frame_equal(
            annotated_activity_df, activity_df.abs(),
        )

    # ------------------------------------------------------------------
    #  Case normalisation
    # ------------------------------------------------------------------

    def test_case_normalization(self) -> None:
        """Mixed-case net targets are uppercased before KB lookup."""
        activity_df = pd.DataFrame(
            {"TF_A": [1.0], "TF_B": [1.0], "TF_C": [1.0]},
            index=pd.Index(["CT0"]),
        )
        net = pd.DataFrame({
            "source": ["TF_A", "TF_B", "TF_C"],
            "target": ["gene_1", "GENE_2", "Gene_3"],
        })
        kb_markers = {"GENE_1", "GENE_2", "GENE_3"}

        _, tf_annotation = compute_tf_relevance(activity_df, net, kb_markers)

        # All three TFs have a matching KB gene after case normalisation
        assert (tf_annotation["kb_overlap_count"] == 1).all()

    # ------------------------------------------------------------------
    #  Annotation table format
    # ------------------------------------------------------------------

    def test_tf_annotation_table_format(self) -> None:
        """Return tuple contains (activity, annotation); annotation has expected columns."""
        activity_df = pd.DataFrame(
            {"TF_A": [1.0], "TF_B": [1.0]},
            index=pd.Index(["CT0"]),
        )
        net = pd.DataFrame({
            "source": ["TF_A", "TF_B"],
            "target": ["GENE_1", "GENE_2"],
        })
        kb_markers = {"GENE_1", "GENE_2"}

        result = compute_tf_relevance(activity_df, net, kb_markers)

        # Two-element tuple
        assert len(result) == 2
        annotated_activity_df, tf_annotation = result

        # Annotation columns
        expected_cols = ["tf", "n_targets", "kb_overlap_count", "kb_overlap_ratio"]
        assert list(tf_annotation.columns) == expected_cols

        # One row per TF
        assert len(tf_annotation) == 2

    # ------------------------------------------------------------------
    #  No matching net rows
    # ------------------------------------------------------------------

    def test_no_net_matching_rows(self) -> None:
        """TFs with no rows in net have n_targets=0 and overlap=0."""
        activity_df = pd.DataFrame(
            {"TF_X": [1.0], "TF_Y": [1.0]},
            index=pd.Index(["CT0"]),
        )
        net = pd.DataFrame({
            "source": ["TF_Z"],
            "target": ["GENE_1"],
        })
        kb_markers = {"GENE_1"}

        _, tf_annotation = compute_tf_relevance(activity_df, net, kb_markers)

        assert len(tf_annotation) == 2
        for _, row in tf_annotation.iterrows():
            assert row["n_targets"] == 0
            assert row["kb_overlap_count"] == 0

    # ------------------------------------------------------------------
    #  Invalid activity_df
    # ------------------------------------------------------------------

    def test_invalid_activity_df(self) -> None:
        """Empty activity_df (no columns) returns empty results without crashing."""
        activity_df = pd.DataFrame()
        net = pd.DataFrame({
            "source": ["TF_A"],
            "target": ["GENE_1"],
        })
        kb_markers = {"GENE_1"}

        annotated_activity_df, tf_annotation = compute_tf_relevance(
            activity_df, net, kb_markers,
        )

        assert annotated_activity_df.empty
        assert tf_annotation.empty
        assert list(tf_annotation.columns) == [
            "tf", "n_targets", "kb_overlap_count", "kb_overlap_ratio",
        ]
