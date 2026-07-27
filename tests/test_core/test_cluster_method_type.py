"""Tests: cluster_selection_method accepts None, str, or default."""

from core.config.global_config import GlobalClusteringConfig


class TestClusterSelectionMethodType:
    """Verify cluster_selection_method is Optional[str] — accepts None, str, or default."""

    def test_global_config_accepts_none(self):
        """GlobalClusteringConfig(cluster_selection_method=None) must NOT raise."""
        cfg = GlobalClusteringConfig(cluster_selection_method=None)
        assert cfg.cluster_selection_method is None

    def test_global_config_accepts_string(self):
        """GlobalClusteringConfig(cluster_selection_method="multi_metric") must work."""
        cfg = GlobalClusteringConfig(cluster_selection_method="multi_metric")
        assert cfg.cluster_selection_method == "multi_metric"

    def test_global_config_default(self):
        """GlobalClusteringConfig() default must be "multi_metric"."""
        cfg = GlobalClusteringConfig()
        assert cfg.cluster_selection_method == "multi_metric"
