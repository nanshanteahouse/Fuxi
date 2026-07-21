#!/usr/bin/env python3
"""
utils — Fuxi pipeline shared utilities
========================================

Sub-modules (private, import through this package):
  _path           — WSL detection, data_root, repo_root
  _io             — safe_write, safe_plot
  _logging        — setup_logger
  _config         — resolve_config, species validation, dataset.yaml helpers
  _cross_modality — find_rna_h5ad, find_rna_marker_csv, load_scRNA_markers
  _perf           — PerformanceReport, PerformanceSummary, monitor_performance
  _validation     — validate_adata, validate_pipeline_state
"""

import os
import platform

# ── WSL h5py file locking auto-detection ────────────────────────────
if "microsoft" in platform.release().lower():
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# ── Auto-load .env ──────────────────────────────────────────────────
try:
    from dotenv import load_dotenv

    _dotenv_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        ".env",
    )
    if os.path.exists(_dotenv_path):
        load_dotenv(_dotenv_path, override=True)
except ImportError:
    pass

# ── Re-exports ──────────────────────────────────────────────────────
from core.utils._config import (  # noqa: F401
    _KNOWN_SPECIES_KEYS,
    _find_dataset_yaml,
    _validate_species,
    resolve_config,
)
from core.utils._cross_modality import (  # noqa: F401
    find_rna_h5ad,
    find_rna_marker_csv,
    load_scRNA_markers,
)
from core.utils._io import safe_plot, safe_write  # noqa: F401
from core.utils._logging import setup_logger  # noqa: F401
from core.utils._path import (  # noqa: F401
    _DATA_ROOT_CACHE,
    _REPO_ROOT_CACHE,
    data_root,
    is_wsl,
    repo_root,
    wsl_to_win,
)
from core.utils._perf import (  # noqa: F401
    PerformanceReport,
    PerformanceSummary,
    monitor_performance,
)
from core.utils._validation import (  # noqa: F401
    _STEP_REQUIREMENTS,
    validate_adata,
    validate_pipeline_state,
)
