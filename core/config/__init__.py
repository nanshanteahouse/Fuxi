"""core.config — Unified configuration for all modalities.

Provides Pydantic v2 Config schema, dataset.yaml model, and
field-to-topic mapping for config migration.
"""

from core.config.schema import (
    SILHOUETTE_SAMPLE_THRESHOLD,
    ClusteringSettings,
    Config,
    DataInputConfig,
    HVGSettings,
    QCSettings,
)

__all__ = [
    "Config",
    "DataInputConfig",
    "QCSettings",
    "HVGSettings",
    "ClusteringSettings",
    "SILHOUETTE_SAMPLE_THRESHOLD",
]
