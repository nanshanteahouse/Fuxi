"""core.interaction — Shared cell-cell interaction (CCI) utilities.

Provides LIANA+ wrappers for RNA (permutation) and Spatial (bivariate) CCI analysis.
"""

from core.interaction.cell_interaction import (
    ensure_gene_symbols,
    format_cci_results,
    load_lr_database,
    run_cci_permutation,
    run_cci_spatial,
)

__all__ = [
    "ensure_gene_symbols",
    "load_lr_database",
    "run_cci_permutation",
    "run_cci_spatial",
    "format_cci_results",
]
