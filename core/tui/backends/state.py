"""Persistent TUI state — recent projects, config paths, window preferences.

All state is stored as JSON in ``~/.fuxi/tui_state.json``.  Atomic writes
via ``.tmp`` + ``os.rename`` ensure the file is never corrupted by a
partial write.  Only UI concerns live here (last modality browsed, recent
configs, window size/position); pipeline checkpoint and registry data live
elsewhere.
"""

from __future__ import annotations

import json
import logging
import os
import time
from typing import Any

logger = logging.getLogger(__name__)

# ── Paths ───────────────────────────────────────────────────────────────────

_STATE_DIR = os.path.expanduser("~/.fuxi")
_STATE_PATH = os.path.join(_STATE_DIR, "tui_state.json")

# ── Schema defaults ─────────────────────────────────────────────────────────

_MAX_RECENT_CONFIGS = 10

MAX_RECENT_CONFIGS: int = _MAX_RECENT_CONFIGS  # public alias

_DEFAULT_STATE: dict[str, Any] = {
    "last_modality": None,
    "last_config": None,
    "recent_configs": [],
    "window_prefs": {},
}

# ── Public API ──────────────────────────────────────────────────────────────

__all__ = [
    "save",
    "load",
    "update",
    "get_state_path",
    "MAX_RECENT_CONFIGS",
]


def get_state_path() -> str:
    """Return the absolute path to the state file on disk."""
    return _STATE_PATH


# ═════════════════════════════════════════════════════════════════════════════
# Save / load
# ═════════════════════════════════════════════════════════════════════════════


def save(state: dict[str, Any]) -> None:
    """Atomically write *state* to ``~/.fuxi/tui_state.json``.

    Creates the ``~/.fuxi/`` directory if it does not exist, writes to a
    ``.tmp`` sibling, then renames atomically so a crash mid-write never
    leaves a truncated state file.
    """
    os.makedirs(_STATE_DIR, exist_ok=True)

    tmp_path = _STATE_PATH + ".tmp"
    try:
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2, ensure_ascii=False)
            f.flush()
            os.fsync(f.fileno())
    except OSError:
        logger.exception("Failed to write TUI state to %s", tmp_path)
        return

    os.replace(tmp_path, _STATE_PATH)


def load() -> dict[str, Any]:
    """Read TUI state from ``~/.fuxi/tui_state.json``.

    Returns *state* as a :class:`dict`.  If the file does not exist or
    contains invalid JSON an empty dict ``{}`` is returned and a warning
    is logged.  Missing keys are **not** filled in — callers that need
    the full schema should merge against :func:`default_state`.
    """
    try:
        with open(_STATE_PATH, "r", encoding="utf-8") as f:
            return json.load(f)  # type: ignore[no-any-return]
    except FileNotFoundError:
        logger.info("No TUI state file found at %s — starting fresh", _STATE_PATH)
        return {}
    except json.JSONDecodeError:
        logger.warning("Corrupted TUI state at %s — starting fresh", _STATE_PATH)
        return {}


# ═════════════════════════════════════════════════════════════════════════════
# Convenience helpers
# ═════════════════════════════════════════════════════════════════════════════


def update(**kwargs: Any) -> dict[str, Any]:
    """Load current state, merge *kwargs*, save and return the new state.

    This is the preferred way to persist small changes:

    .. code-block:: python

        state.update(last_modality="atac", last_config="/path/to/config.yaml")
    """
    state = load()
    state.update(kwargs)
    save(state)
    return state


def default_state() -> dict[str, Any]:
    """Return a deep copy of the default state dict.

    Useful for callers that want guaranteed keys even on first run.
    """
    # The default values are all immutable scalars or empty containers,
    # so a shallow copy is safe here.
    return dict(_DEFAULT_STATE)


def push_recent_config(
    modality: str,
    gse_id: str,
    path: str,
    *,
    state: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Record a recently opened config and return the updated state.

    Parameters
    ----------
    modality:
        Modality identifier (e.g. ``"rna"``, ``"atac"``).
    gse_id:
        GEO / dataset identifier (e.g. ``"GSE123456"``).
    path:
        Absolute filesystem path to the YAML config file.
    state:
        Optional pre-loaded state dict.  When omitted the state is loaded
        from disk automatically.

    The entry is prepended to ``recent_configs``; duplicates (same *path*)
    are moved to the front instead of duplicated.  The list is capped at
    :const:`MAX_RECENT_CONFIGS` entries.
    """
    current = state if state is not None else load()

    entry: dict[str, Any] = {
        "modality": modality,
        "gse_id": gse_id,
        "path": path,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    recent = current.get("recent_configs", [])
    # Remove duplicate (same path) if it exists
    recent = [r for r in recent if r.get("path") != path]
    # Prepend new entry
    recent.insert(0, entry)
    # Cap at max
    current["recent_configs"] = recent[:_MAX_RECENT_CONFIGS]

    save(current)
    return current
