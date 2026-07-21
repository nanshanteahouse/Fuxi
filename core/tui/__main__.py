#!/usr/bin/env python3
"""TUI entry point — launches ``FuxiTUI``.

Usage:
    python -m core.tui
"""

import os
import sys

# Ensure repo root is on sys.path
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from core.tui.app import FuxiTUI  # noqa: E402

if __name__ == "__main__":
    FuxiTUI().run()
