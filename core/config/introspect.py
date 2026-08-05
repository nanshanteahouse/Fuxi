"""Pydantic-driven config introspection — shared field walker for the Config tree.

Provides utilities to introspect the ``Config`` Pydantic model tree and
produce a flat ``{"dotted.key": FieldInfo}`` map.  This is the single
source of field discovery for schema-driven tooling: the TUI config
editor form and the preprocess config scaffold both consume it.
"""

from typing import Any, Dict, Union, get_args, get_origin

from pydantic import BaseModel
from pydantic.fields import FieldInfo

from core.config.schema import Config

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
        # Pydantic ForwardRefs originate from core/config/schema.py.
        # Pass the module's globals so forward refs like "PseudobulkDESettings"
        # resolve correctly when evaluated.
        import sys

        mod = sys.modules.get("core.config.schema")
        globals_dict = vars(mod) if mod else None
        resolved = annotation.evaluate(globals=globals_dict)
        if resolved is not None:
            return resolved
    return annotation
