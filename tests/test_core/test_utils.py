"""Skeleton tests for core/utils.py."""
import os

from core.utils import (
    safe_write,
    safe_plot,
    setup_logger,
    resolve_config,
    validate_adata,
    repo_root,
    find_rna_h5ad,
    find_rna_marker_csv,
    load_scRNA_markers,
)


class TestCoreUtilsImport:
    """Verify that core utility symbols are importable."""

    def test_import_core_utils(self) -> None:
        """All core/utils.py public symbols should import without error."""
        # names — just ensure they resolved above
        assert safe_write is not None
        assert safe_plot is not None
        assert setup_logger is not None
        assert resolve_config is not None
        assert validate_adata is not None
        assert repo_root is not None
        assert find_rna_h5ad is not None
        assert find_rna_marker_csv is not None
        assert load_scRNA_markers is not None

    def test_repo_root_is_absolute(self) -> None:
        """repo_root() should return an absolute path that ends with Fuxi."""
        root = repo_root()
        assert os.path.isabs(root)
        assert root.endswith("Fuxi")
