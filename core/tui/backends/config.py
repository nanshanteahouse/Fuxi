"""Pydantic-driven config introspection for generating form widgets.

The field-walker utilities (``get_config_fields`` and friends) live in
``core.config.introspect`` — the shared schema-introspection layer used by
both the TUI form generator and the preprocess config scaffold.

This module keeps the widget-type mapping and YAML I/O helpers that are
TUI-specific, re-exporting the introspection functions for backward
compatibility."""

from typing import Any, Literal, Union, get_args, get_origin

import yaml
from pydantic.fields import FieldInfo

from core.config.introspect import (  # noqa: F401  re-exported for callers
    _maybe_resolve_forward_ref,
    _resolve_base_model,
    get_config_fields,
)
from core.config.schema import Config
from core.utils._config import resolve_config

# ═══════════════════════════════════════════════════════════════════
# Widget type mapping
# ═══════════════════════════════════════════════════════════════════


def field_to_widget_type(field_info: FieldInfo, field_name: str = "") -> str:
    """Map a Pydantic ``FieldInfo`` to a TUI widget type string.

    *field_name* is an optional dotted key (e.g. ``"ai.api_key"``) that allows
    the caller to annotate special widget types (e.g. password masking).
    When not provided, the mapping is purely annotation-based.

    Returns one of:
      ``"integer"``, ``"float"``, ``"text"``, ``"switch"``, ``"select"``,
      ``"password"``, ``"textarea"``

    Special cases
    -------------
    * ``api_key`` fields → ``"password"`` (masked input).
    * ``Optional[X]`` → widget type of ``X``.
    """
    # Sensitive-field detection by dotted name
    if field_name.endswith(".api_key") or field_name == "api_key":
        return "password"

    annotation = field_info.annotation
    if annotation is None:
        return "text"

    return _annotation_to_widget(annotation)


def _annotation_to_widget(annotation: Any) -> str:
    """Recursively resolve a type annotation to a widget name string."""
    origin = get_origin(annotation)

    # ── Optional / Union ──────────────────────────────────────────
    if origin is Union:
        args = get_args(annotation)
        non_none = [a for a in args if a is not type(None)]
        if len(non_none) == 1:
            # Optional[X] — unwrap to the inner type's widget
            return _annotation_to_widget(non_none[0])
        # Multi-branch Union is not supported in auto-forms; fall through
        return "text"

    # ── Literal (enum-like choices) ───────────────────────────────
    if origin is Literal:
        return "select"

    # ── Generic collections ───────────────────────────────────────
    if origin in (list,):
        return "text"  # comma-separated input

    if origin in (dict,):
        args = get_args(annotation)
        if len(args) == 2 and get_origin(args[1]) in (list,):
            return "textarea"  # key: val1,val2 format
        return "textarea"

    # ── Primitives ────────────────────────────────────────────────
    if annotation is int:
        return "integer"
    if annotation is float:
        return "float"
    if annotation is bool:
        return "switch"
    if annotation is str:
        return "text"

    # ── Fallback ──────────────────────────────────────────────────
    return "text"


# ═══════════════════════════════════════════════════════════════════
# YAML I/O helpers
# ═══════════════════════════════════════════════════════════════════


def load_yaml_config(path: str) -> Config:
    """Load a YAML config file and return a resolved ``Config`` instance.

    Thin wrapper around :func:`~core.utils._config.resolve_config`.
    """
    return resolve_config(path)


def save_yaml_config(config: Config, path: str) -> None:
    """Serialize *config* to YAML and write to *path*.

    Uses ``model_dump()`` so all declared fields (including nested sub-models)
    are included. Computed properties and private attributes are excluded.
    """
    data = config.model_dump()
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(
            data,
            f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )
