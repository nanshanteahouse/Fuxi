"""Path utilities — WSL detection, data root resolution, repo root."""

import os
import platform
from typing import Optional


# ── Module-level caches ──────────────────────────────────────────
_DATA_ROOT_CACHE: Optional[str] = None
_REPO_ROOT_CACHE: Optional[str] = None


def is_wsl() -> bool:
    """True when running inside Windows Subsystem for Linux."""
    return ('microsoft' in platform.uname().release.lower()
            and os.path.exists('/mnt/c'))


def data_root() -> str:
    """Absolute root of the raw data tree.

    Resolved in order of precedence:
      1. FUXI_DATA_ROOT  env var (canonical name)
      2. SCRNA_DATA_ROOT env var (legacy name, backward compat)

    If neither is set the function raises RuntimeError with setup
    instructions — this is intentional: every machine has its own
    data layout, so the path must be configured explicitly.

    Cached after first resolution.
    """
    global _DATA_ROOT_CACHE
    if _DATA_ROOT_CACHE is None:
        _DATA_ROOT_CACHE = (
            os.environ.get('FUXI_DATA_ROOT')
            or os.environ.get('SCRNA_DATA_ROOT')
        )
        if not _DATA_ROOT_CACHE:
            raise RuntimeError(
                "Data root not configured.\n"
                "  Set the FUXI_DATA_ROOT environment variable to the\n"
                "  directory containing your GEO dataset folders, e.g.:\n"
                '    export FUXI_DATA_ROOT=/mnt/e/data   # WSL\n'
                '    set FUXI_DATA_ROOT=E:/data          # Windows'
            )
    return _DATA_ROOT_CACHE


def repo_root() -> str:
    """Absolute path to this repository's root, located from __file__.

    Override (rare) via SCRNA_REPO_ROOT. Cached at import.
    """
    global _REPO_ROOT_CACHE
    if _REPO_ROOT_CACHE is None:
        _REPO_ROOT_CACHE = os.environ.get('SCRNA_REPO_ROOT') or os.path.abspath(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        )
    return _REPO_ROOT_CACHE


def wsl_to_win(path: str) -> str:
    """Translate /mnt/X/... -> X:/...; pass-through if not /mnt-prefixed."""
    if path.startswith('/mnt/') and len(path) > 5:
        return f"{path[5]}:/{path[6:]}"
    return path
