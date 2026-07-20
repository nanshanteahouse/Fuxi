"""Pydantic-driven config introspection for generating form widgets.

Provides utilities to introspect the Config Pydantic model tree
and produce field metadata that can drive automatic form generation
in TUI screens.
"""

from typing import Any, Dict, get_args, get_origin, Literal, Union

import yaml
from pydantic import BaseModel
from pydantic.fields import FieldInfo

from core.config.schema import Config
from core.utils._config import resolve_config


# ═══════════════════════════════════════════════════════════════════
# Field introspection
# ═══════════════════════════════════════════════════════════════════


def get_config_fields(
    model_cls: type[BaseModel] = Config,
    prefix: str = "",
) -> Dict[str, FieldInfo]:
    """Recursively introspect *model_cls* fields, returning a flat dict of
    ``{"dotted.key": FieldInfo}`` pairs.

    Sub-models that are ``BaseModel`` subclasses are recursed into with
    their field names prefixed by the parent field name (e.g. ``qc.min_genes``).
    ``Optional[BaseModel]`` is also unwrapped and recursed.

    Because every model in the Config tree uses ``extra="forbid"``, only
    declared fields are returned — no extra / unexpected keys leak through.
    """
    fields: Dict[str, FieldInfo] = {}

    for field_name, field_info in model_cls.model_fields.items():
        key = f"{prefix}{field_name}" if prefix else field_name
        annotation = field_info.annotation

        inner_model = _resolve_base_model(annotation)

        if inner_model is not None:
            # Recurse into the sub-model, extending the dotted prefix
            sub = get_config_fields(inner_model, prefix=f"{key}.")
            fields.update(sub)
        else:
            # Leaf field — keep as-is
            fields[key] = field_info

    return fields


def _resolve_base_model(annotation: Any) -> type[BaseModel] | None:
    """If *annotation* is (or ``Optional`` of) a ``BaseModel`` subclass, return
    it.  Otherwise return ``None``.

    Handles ``ForwardRef`` annotations (e.g. forward-declared nested models
    in the same module) by calling ``.evaluate()``.
    """
    # Resolve forward references (e.g. PseudobulkDESettings inside DESettings)
    resolved = _maybe_resolve_forward_ref(annotation)

    origin = get_origin(resolved)

    # Unwrap Optional[X] -> Union[X, None]
    if origin is Union:
        args = get_args(resolved)
        non_none = [a for a in args if a is not type(None)]
        if len(non_none) == 1:
            candidate = _maybe_resolve_forward_ref(non_none[0])
            if isinstance(candidate, type) and issubclass(candidate, BaseModel):
                return candidate
        return None

    # Direct BaseModel subclass
    if isinstance(resolved, type) and issubclass(resolved, BaseModel):
        return resolved

    return None


def _maybe_resolve_forward_ref(annotation: Any) -> Any:
    """If *annotation* is a ``ForwardRef``, evaluate and return the resolved
    type.  Otherwise return *annotation* unchanged.
    """
    from annotationlib import ForwardRef  # stdlib since Python 3.12
    if isinstance(annotation, ForwardRef):
        resolved = annotation.evaluate()
        if resolved is not None:
            return resolved
    return annotation


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
