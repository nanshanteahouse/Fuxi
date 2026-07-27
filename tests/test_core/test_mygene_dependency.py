"""Test mygene import fallback in cell_interaction.ensure_gene_symbols.

Verifies:
1. When mygene is unavailable and Ensembl IDs are present, an ImportError is raised.
2. When all var_names are gene symbols, no import is attempted (early return).
"""

from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from core.interaction.cell_interaction import ensure_gene_symbols


def _make_adata(var_names):
    """Create a lightweight mock AnnData with given var_names."""
    adata = MagicMock()
    adata.var_names = pd.Index(var_names)
    adata.n_vars = len(var_names)
    return adata


class TestMygeneDependency:
    """Verify mygene import guard and fallback behavior."""

    def test_import_error_fallback(self):
        """Ensembl IDs present + mygene unavailable -> ImportError."""
        adata = _make_adata(["ENSG00000123456", "ENSG00000789012"])
        with patch.dict("sys.modules", {"mygene": None}):
            with pytest.raises(ImportError) as exc_info:
                ensure_gene_symbols(adata)
        assert "mygene" in str(exc_info.value).lower()

    def test_normal_path_not_affected(self):
        """All gene symbols -> early return, adata unchanged."""
        adata = _make_adata(["RHO", "GNAT1", "RLBP1"])
        result = ensure_gene_symbols(adata)
        assert result is adata
