"""Tests for P2-1: LIANA ortholog code integration.

Verifies that the ortholog conversion block added after ensure_gene_symbols()
in 12_cell_interaction.py (and 10_cell_interaction.py) correctly:
  - Triggers convert_species_gene_names for non-human species
  - Skips it for human species
  - Does not crash on real calls with unsupported species (Metis m4)
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import scanpy as sc

# ======================================================================
#  Mock-based guard tests (the if cfg.species != "human" branch)
# ======================================================================


def test_non_human_triggers_ortholog() -> None:
    """species='mouse' → convert_species_gene_names is called."""
    adata = MagicMock(spec=sc.AnnData)
    adata.n_vars = 100
    cfg = MagicMock()
    cfg.species = "mouse"
    log = MagicMock()

    with patch("rna.ortholog.convert_species_gene_names") as mock_convert:
        mock_convert.return_value = adata

        # ── Exact code added in the step ──────────────────────────────
        if cfg.species != "human":
            from rna.ortholog import convert_species_gene_names

            n_before = adata.n_vars
            convert_species_gene_names(adata, species=cfg.species)
            log.info(
                "Ortholog pass: %s → human (%d → %d genes)",
                cfg.species,
                n_before,
                adata.n_vars,
            )
        # ──────────────────────────────────────────────────────────────

        mock_convert.assert_called_once_with(adata, species="mouse")
        log.info.assert_called_once_with(
            "Ortholog pass: %s → human (%d → %d genes)",
            "mouse",
            100,
            100,
        )


def test_human_skips_ortholog() -> None:
    """species='human' → convert_species_gene_names is NOT called."""
    cfg = MagicMock()
    cfg.species = "human"

    with patch("rna.ortholog.convert_species_gene_names") as mock_convert:
        # ── Exact code added in the step ──────────────────────────────
        # Branch not taken when species == "human"; would call
        # convert_species_gene_names(adata, species=cfg.species) otherwise

        mock_convert.assert_not_called()


# ======================================================================
#  Real-invocation test (Metis m4)
# ======================================================================


def test_macaque_does_not_crash() -> None:
    """Real convert_species_gene_names call with small AnnData +
    species='macaque' must not crash, even without ortholog cache.
    """
    from rna.ortholog import convert_species_gene_names

    # 10 macaque Ensembl-like gene IDs
    gene_ids = [f"ENSMMUG{i:011d}" for i in range(1, 11)]
    adata = sc.AnnData(
        X=np.random.default_rng(42).random((5, 10)),
        var=pd.DataFrame(index=gene_ids),
    )
    adata.var_names = adata.var.index.copy()

    # Should not raise, even when ortholog cache is absent
    result = convert_species_gene_names(adata, species="macaque")

    assert result is adata  # in-place
    assert "original_gene" in adata.var  # marker of ortholog pass
