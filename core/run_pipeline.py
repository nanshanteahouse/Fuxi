#!/usr/bin/env python3
"""Thin CLI wrapper — delegates to ``core.pipeline.runner`` for backward compatibility.

Usage (unchanged):
    python core/run_pipeline.py --modality rna --list
    python core/run_pipeline.py --modality atac --step 0 --config ...
"""

import sys
import os

# Ensure repo root is on sys.path
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from core.pipeline.runner import main

if __name__ == "__main__":
    main()
