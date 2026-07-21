#!/usr/bin/env python3
"""
Fuxi Pre-Processing Pipeline
==============================
Automates the gap between "files downloaded from GEO" and
"pipeline-ready config + dataset manifest."

Entry point:
    from core.preprocess import run_preprocess
    run_preprocess('GSE12345')

CLI:
    python core/preprocess/preprocessor.py --gse GSE12345
"""

from .matrix_loader import TEMPLATE_MAP, _detect_primary_format, generate_config
from .metadata_parser import _resolve_input_dir, generate_dataset_yaml
from .preprocessor import main, run_preprocess

__all__ = [
    "run_preprocess",
    "main",
    "generate_dataset_yaml",
    "generate_config",
    "TEMPLATE_MAP",
    "_detect_primary_format",
    "_resolve_input_dir",
]
