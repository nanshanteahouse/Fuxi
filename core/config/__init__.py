"""core.config — Unified configuration for all modalities.

Provides Pydantic v2 Config schema, dataset.yaml model, and
field-to-topic mapping for config migration.
"""

from core.config.schema import (
    Config,
    DataInputConfig,
    QCConfig,
    HvgConfig,
    ClusteringConfig,
    SILHOUETTE_SAMPLE_THRESHOLD,
)

__all__ = [
    "Config",
    "DataInputConfig",
    "QCConfig",
    "HvgConfig",
    "ClusteringConfig",
    "SILHOUETTE_SAMPLE_THRESHOLD",
]
