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

from .preprocessor import run_preprocess, main
from .metadata_parser import generate_dataset_yaml, _resolve_input_dir
from .matrix_loader import generate_config, _detect_primary_format, TEMPLATE_MAP

__all__ = ['run_preprocess', 'main', 'generate_dataset_yaml', 'generate_config']
