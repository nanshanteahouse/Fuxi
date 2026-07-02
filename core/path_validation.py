#!/usr/bin/env python3
"""
path_validation.py — Safe path traversal guards.

Provides ``validate_safe_path()`` to prevent directory-traversal attacks
when resolving user-supplied paths against a trusted base directory.
"""

import os


def validate_safe_path(path: str, base_dir: str) -> str:
    """Resolve *path* and assert it is within *base_dir*.

    Both paths are resolved to their real (canonical) absolute form via
    ``os.path.realpath()`` before comparison so that ``..`` components,
    symlinks, and relative references are neutralised.

    Parameters
    ----------
    path : str
        User-supplied path (may be absolute or relative).
    base_dir : str
        Trusted base directory that *path* must reside within.

    Returns
    -------
    str
        The resolved absolute path (for convenience).

    Raises
    ------
    ValueError
        If the resolved path does not start with the resolved *base_dir*,
        indicating a traversal attempt outside the trusted directory.
    FileNotFoundError
        If the resolved path does not exist on disk.
    """
    resolved_path = os.path.realpath(path)
    resolved_base = os.path.realpath(base_dir)

    if not os.path.exists(resolved_path):
        raise FileNotFoundError(
            f"Path does not exist after resolution: {resolved_path}"
        )

    if not resolved_path.startswith(resolved_base):
        raise ValueError(
            f"Path traversal detected: {resolved_path} is outside "
            f"the allowed base directory {resolved_base}"
        )

    return resolved_path
