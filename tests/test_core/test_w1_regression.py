"""Regression tests for W1-01 through W1-03 fixes.

W1-01: safe_write WSL — tmp+mv strategy on /mnt/ paths
W1-02: format_cci_results descending sort — ascending=False parameter
W1-03: guess_genome pipeline keys — 11 known species map correctly
"""

import os
import shutil
from unittest.mock import MagicMock, patch

import pandas as pd

from core.utils import safe_write
from core.preprocess.format_detector import guess_genome
from core.interaction.cell_interaction import format_cci_results


# ── W1-01: safe_write WSL ────────────────────────────────────────────────


def test_safe_write_wsl():
    """W1-01: tmp+mv strategy on /mnt/ paths (WSL fix).

    Verifies that when target starts with /mnt/, the function:
      1. Creates tmpdir
      2. Writes to tmp path
      3. Moves tmp file to target
    """
    mock_adata = MagicMock()

    with patch("core.utils.os.makedirs") as mock_makedirs:
        with patch("core.utils.shutil.move") as mock_move:
            with patch("core.utils.os.path.getsize", return_value=1_000_000):
                safe_write(mock_adata, "/mnt/data/test.h5ad")

                # tmpdir created
                mock_makedirs.assert_called_once_with("/tmp/Fuxi", exist_ok=True)

                # adata.write called with tmp path (not target)
                mock_adata.write.assert_called_once_with(
                    "/tmp/Fuxi/test.h5ad", compression="gzip"
                )

                # shutil.move called with tmp -> target
                mock_move.assert_called_once_with(
                    "/tmp/Fuxi/test.h5ad", "/mnt/data/test.h5ad"
                )


# ── W1-02: format_cci_results descending sort ────────────────────────────


def test_format_cci_results_descending():
    """W1-02: top row has highest morans value when ascending=False.

    Creates a DataFrame with morans [0.1, 0.9, 0.5]; the 0.9 row must
    appear first after descending sort.
    """
    df = pd.DataFrame({
        "ligand":     ["A", "B", "C"],
        "receptor":   ["X", "Y", "Z"],
        "source":     ["S1", "S2", "S3"],
        "target":     ["T1", "T2", "T3"],
        "morans":     [0.1, 0.9, 0.5],
    })
    result = format_cci_results(df, ascending=False, pval_col="morans")
    assert result.iloc[0]["morans"] == 0.9, (
        f"Expected top morans=0.9, got {result.iloc[0]['morans']}"
    )


# ── W1-03: guess_genome pipeline keys ────────────────────────────────────


def test_guess_genome_pipeline_keys():
    """W1-03: all 11 known species return correct genome strings."""
    expected = {
        "human":      "hg38",
        "mouse":      "mm10",
        "rat":        "rn6",
        "zebrafish":  "danRer11",
        "cow":        "bosTau9",
        "pig":        "susScr11",
        "macaque":    "rheMac10",
        "chicken":    "galGal6",
        "drosophila": "dm6",
        "c_elegans":  "ce11",
        "frog":       "xenTro10",
    }
    for species, expected_genome in expected.items():
        result = guess_genome(species)
        assert result == expected_genome, (
            f"{species}: expected {expected_genome!r}, got {result!r}"
        )


def test_guess_genome_unknown_species():
    """W1-03: unknown species returns empty string (no crash)."""
    assert guess_genome("gorilla") == ""
    assert guess_genome("") == ""
    assert guess_genome("homo_sapiens") == ""  # pipeline key is "human"
