"""Config editor screen — form-based YAML config editor driven by Pydantic introspection.

Dynamically generates form widgets from the real ``Config`` schema
via :func:`~core.tui.backends.config.get_config_fields`.  Supports
load / save / open-in-editor workflows.
"""

from __future__ import annotations

import ast
import logging
import os
import subprocess
from typing import Any, get_args, get_origin, Literal

from pydantic.fields import FieldInfo

from textual.app import ComposeResult
from textual.containers import Horizontal, ScrollableContainer, Vertical
from textual.screen import Screen
from textual.widgets import (
    Button,
    Collapsible,
    Footer,
    Header,
    Input,
    Label,
    Select,
    Static,
    Switch,
    TextArea,
)

from core.tui.backends.config import (
    field_to_widget_type,
    get_config_fields,
    load_yaml_config,
    save_yaml_config,
)
from core.config.schema import Config

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

# Canonical display labels for sub-model sections
SECTION_LABELS: dict[str, str] = {
    "data_input": "Data Input",
    "sample_meta": "Sample Metadata",
    "qc": "QC",
    "scrublet": "Scrublet",
    "normalization": "Normalization",
    "hvg": "HVG",
    "pca": "PCA",
    "harmony": "Harmony",
    "clustering": "Clustering",
    "marker": "Marker",
    "de": "Differential Expression",
    "trajectory": "Trajectory",
    "enrichment": "Enrichment",
    "grn": "GRN",
    "cci": "Cell-Cell Interaction",
    "downsample": "Downsample",
    "spatial": "Spatial",
    "atac": "ATAC",
    "execution": "Execution",
    "ai": "AI",
}

# Canonical section ordering (all sections not listed are pushed to the end)
_CANONICAL_ORDER = list(SECTION_LABELS.keys())


def _section_name(key: str) -> str:
    """Extract the section name from a dotted key.

    ``"qc.min_genes"`` → ``"qc"``,  ``"modality"`` → ``"__root__"``.
    """
    return key.split(".")[0] if "." in key else "__root__"


def _display_name(key: str) -> str:
    """Human-readable label from a dotted field key.

    ``"qc.min_genes"`` → ``"Min Genes"``.
    """
    name = key.rsplit(".", 1)[-1]
    return name.replace("_", " ").title()


def _unflatten_config(flat: dict[str, Any]) -> dict[str, Any]:
    """Convert a flat ``{"a.b": v, "c": w}`` dict to nested ``{"a": {"b": v}, "c": w}``."""
    result: dict[str, Any] = {}
    for key, value in flat.items():
        parts = key.split(".")
        target = result
        for part in parts[:-1]:
            target = target.setdefault(part, {})
        target[parts[-1]] = value
    return result


# ═══════════════════════════════════════════════════════════════════════════
# Screen
# ═══════════════════════════════════════════════════════════════════════════


class ConfigEditorScreen(Screen):
    """Form-based YAML config editor driven by Pydantic introspection.

    Builds dynamic form widgets from the ``Config`` Pydantic model tree
    via :func:`~core.tui.backends.config.get_config_fields`.  Fields are
    grouped into collapsible sections by sub-model (QC, HVG, Clustering,
    etc.).  Three toolbar buttons drive the workflow:

    * **Load YAML** — read a path from the path input, load via
      ``load_yaml_config()``, populate all widgets.
    * **Save YAML** — collect widget values into a flat dict, unflatten to
      nested dict, construct ``Config(**nested)``, write via
      ``save_yaml_config()``.
    * **Open in $EDITOR** — launch the system ``$EDITOR`` (or ``vim``) on
      the loaded config path, then re-load the file on exit.
    """

    id = "config-editor"

    DEFAULT_CSS = """
    ConfigEditorScreen {
        height: 100%;
    }

    #editor-header {
        padding: 1 2;
        background: $bg-medium;
        border-bottom: solid $accent;
        text-style: bold;
        color: $accent;
        height: 3;
    }

    #toolbar {
        padding: 0 1;
        background: $bg-medium;
        border-bottom: solid $border;
        height: auto;
    }

    #toolbar Button {
        margin: 0 1;
    }

    #toolbar > #config-path-input {
        width: 1fr;
        margin: 0 0 0 2;
    }

    #form-container {
        height: 1fr;
        overflow-y: auto;
        padding: 1 2;
    }

    #status-bar {
        background: $bg-dark;
        color: $text-muted;
        padding: 0 2;
        height: 1;
        border-top: solid $border;
        text-style: italic;
    }

    ConfigEditorScreen Collapsible {
        margin: 1 0;
        border: solid $border;
    }

    ConfigEditorScreen Collapsible > .collapsible-content {
        padding: 1 2;
    }

    .field-row {
        height: auto;
        margin: 0 0 1 0;
    }

    .field-label {
        width: 30;
        padding: 0 1;
        text-style: bold;
        color: $text-secondary;
    }

    .field-input {
        width: 1fr;
    }

    .field-switch {
        margin: 0 1;
    }

    .field-select {
        width: 1fr;
    }

    .field-textarea {
        width: 1fr;
        height: 6;
    }

    #empty-message {
        color: $text-muted;
        text-style: italic;
        padding: 2 4;
        content-align: center middle;
        height: 100%;
    }

    /* ── Numeric validation feedback ──────────────────────────────── */
    .field-invalid {
        background: $error 10%;
        border: solid $error;
    }
    """

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._config: Config | None = None
        self._config_path: str | None = None
        # Introspect once — pure, no I/O
        self._last_status = ""
        self._fields: dict[str, FieldInfo] = get_config_fields()
        # Group by section for compose ordering
        self._section_keys: dict[str, list[str]] = {}
        for key in self._fields:
            section = _section_name(key)
            self._section_keys.setdefault(section, []).append(key)
        # Filled during compose
        self._field_widgets: dict[str, Any] = {}

    # ── compose ─────────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static(" Config Editor", id="editor-header")

        with Horizontal(id="toolbar"):
            yield Button("Load YAML", id="btn-load", variant="primary")
            yield Button("Save YAML", id="btn-save", variant="success")
            yield Button("Open in $EDITOR", id="btn-editor")
            yield Input(
                placeholder="Config file path …",
                id="config-path-input",
            )

        with ScrollableContainer(id="form-container"):
            yield from self._compose_sections()

        yield Static("Ready. Load a YAML file to edit.", id="status-bar")
        yield Footer()

    # ── section rendering ──────────────────────────────────────────────────

    def _compose_sections(self):
        """Yield collapsible section blocks for all field groups."""
        # ── Root (top-level) fields ────────────────────────────────────
        root_keys = self._section_keys.get("__root__", [])
        if root_keys:
            with Collapsible(title="General Settings", collapsed=False):
                with Vertical(classes="collapsible-content"):
                    for key in root_keys:
                        yield self._make_field_row(key)

        # ── Sub-model sections (sorted canonically) ────────────────────
        ordered = sorted(
            (s for s in self._section_keys if s != "__root__"),
            key=lambda s: _CANONICAL_ORDER.index(s) if s in _CANONICAL_ORDER else 999,
        )
        for section in ordered:
            label = SECTION_LABELS.get(section, section.replace("_", " ").title())
            keys = self._section_keys[section]
            with Collapsible(title=label, collapsed=True):
                with Vertical(classes="collapsible-content"):
                    for key in keys:
                        yield self._make_field_row(key)

    # ── field widget factory ───────────────────────────────────────────────

    def _make_field_row(self, key: str) -> Horizontal:
        """Build a labelled field row for the dotted *key*."""
        field_info = self._fields[key]
        widget_type = field_to_widget_type(field_info, field_name=key)
        label_text = _display_name(key)

        # Extract numeric constraints from Pydantic field metadata
        ge = le = None
        for item in field_info.metadata:
            if hasattr(item, "ge"):
                ge = item.ge
            if hasattr(item, "le"):
                le = item.le

        default_val = field_info.default
        widget = self._create_widget(key, widget_type, field_info, default_val, ge, le)

        self._field_widgets[key] = widget

        label = Label(label_text, classes="field-label")
        # Switch: widget on left, label on right
        if widget_type == "switch":
            return Horizontal(widget, label, classes="field-row")
        return Horizontal(label, widget, classes="field-row")

    def _create_widget(
        self,
        key: str,
        widget_type: str,
        field_info: FieldInfo,
        default_val: Any,
        ge: Any | None,
        le: Any | None,
    ) -> Any:
        """Instantiate the correct Textual widget for *widget_type*."""
        wid = f"field-{key.replace('.', '_')}"

        if widget_type == "integer":
            raw = self._safe_str(default_val)
            w = Input(value=raw, type="integer", id=wid, classes="field-input")
            # Store constraints for validation
            w._field_ge = ge
            w._field_le = le
            return w

        if widget_type == "float":
            raw = self._safe_str(default_val)
            w = Input(value=raw, type="number", id=wid, classes="field-input")
            w._field_ge = ge
            w._field_le = le
            return w

        if widget_type == "text":
            return Input(value=self._safe_str(default_val), id=wid, classes="field-input")

        if widget_type == "switch":
            return Switch(
                value=bool(default_val) if default_val is not None else False,
                id=wid,
                classes="field-switch",
            )

        if widget_type == "password":
            return Input(
                value=self._safe_str(default_val),
                password=True,
                id=wid,
                classes="field-input",
            )

        if widget_type == "select":
            choices = self._extract_literal_choices(field_info)
            return Select(
                options=choices,
                value=str(default_val) if default_val is not None else None,
                id=wid,
                classes="field-select",
                prompt="Select\u2026",
            )

        if widget_type == "textarea":
            return TextArea(
                text=self._serialize_complex(default_val),
                id=wid,
                classes="field-textarea",
            )

        # Fallback
        return Input(value=self._safe_str(default_val), id=wid, classes="field-input")

    # ── value helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _safe_str(val: Any) -> str:
        """Return a string representation, empty for ``None``."""
        if val is None:
            return ""
        if isinstance(val, bool):
            return ""
        return str(val)

    @staticmethod
    def _extract_literal_choices(field_info: FieldInfo) -> list[tuple[str, str]]:
        """Extract ``(display, value)`` options from a ``Literal`` annotation."""
        annotation = field_info.annotation
        origin = get_origin(annotation)
        if origin is Literal:
            args = get_args(annotation)
            return [(str(v), str(v)) for v in args]
        # Multi-type Union fallback: heuristically extract from field name
        # (rare in practice — mostly Literal-driven)
        return []

    @staticmethod
    def _serialize_complex(val: Any) -> str:
        """Serialise a list or dict for display in a ``TextArea``."""
        if isinstance(val, list):
            return ", ".join(str(v) for v in val)
        if isinstance(val, dict):
            return str(val)
        return str(val) if val is not None else ""

    @staticmethod
    def _parse_complex(text: str, field_info: FieldInfo) -> Any:
        """Parse a ``TextArea`` string back into a list or dict."""
        annotation = field_info.annotation
        origin = get_origin(annotation)

        if origin is list:
            return [item.strip() for item in text.split(",") if item.strip()]

        if origin is dict:
            text = text.strip()
            if not text:
                return {}
            try:
                return ast.literal_eval(text)
            except (ValueError, SyntaxError):
                return {}

        # Plain string fallback
        return text

    # ── button handlers ─────────────────────────────────────────────────────

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Route toolbar button presses."""
        bid = event.button.id
        if bid == "btn-load":
            self._on_load()
        elif bid == "btn-save":
            self._on_save()
        elif bid == "btn-editor":
            self._on_open_editor()

    # ── Load YAML ──────────────────────────────────────────────────────────

    def _on_load(self) -> None:
        """Read the path input and load the config file."""
        path_input = self.query_one("#config-path-input", Input)
        path = path_input.value.strip()

        if not path:
            self._update_status("[yellow]Enter a config file path first.[/]")
            return
        if not os.path.isfile(path):
            self._update_status(f"[red]File not found:[/] {path}")
            return

        try:
            config = load_yaml_config(path)
        except Exception as exc:
            logger.exception("Failed to load config")
            self._update_status(f"[red]Load failed:[/] {exc}")
            return

        self._config = config
        self._config_path = path
        self._populate_from_config(config)
        self._update_status(f"Loaded: {path}")

    def _populate_from_config(self, config: Config) -> None:
        """Walk the ``Config`` tree via dotted keys and set widget values."""
        for key, field_info in self._fields.items():
            parts = key.split(".")
            value: Any = config
            try:
                for part in parts:
                    value = getattr(value, part)
            except AttributeError:
                continue

            widget = self._field_widgets.get(key)
            if widget is None:
                continue
            try:
                self._set_widget_value(widget, key, value)
            except Exception as exc:
                logger.debug("Failed to set %s = %r: %s", key, value, exc)

    @staticmethod
    def _set_widget_value(widget: Any, key: str, value: Any) -> None:
        """Apply *value* to a form widget."""
        if isinstance(widget, Input):
            if value is None:
                widget.value = ""
            elif isinstance(value, bool):
                widget.value = "true" if value else "false"
            else:
                widget.value = str(value)
        elif isinstance(widget, Switch):
            widget.value = bool(value)
        elif isinstance(widget, Select):
            widget.value = str(value) if value is not None else None
        elif isinstance(widget, TextArea):
            if isinstance(value, (list, tuple)):
                widget.text = ", ".join(str(v) for v in value)
            elif isinstance(value, dict):
                widget.text = str(value)
            elif value is None:
                widget.text = ""
            else:
                widget.text = str(value)

    # ── Save YAML ──────────────────────────────────────────────────────────

    def _on_save(self) -> None:
        """Collect form values and write the YAML file."""
        path_input = self.query_one("#config-path-input", Input)
        path = path_input.value.strip()

        if not path:
            self._update_status("[yellow]Enter a save path first.[/]")
            return

        self._do_save(path)

    def _do_save(self, path: str) -> None:
        """Collect widget values, build a fresh ``Config``, and write."""
        try:
            flat: dict[str, Any] = {}
            for key, widget in self._field_widgets.items():
                flat[key] = self._get_widget_value(widget, key)

            nested = _unflatten_config(flat)
            config = Config(**nested)
            save_yaml_config(config, path)

            self._config = config
            self._config_path = path
            self._update_status(f"Saved: {path}")
        except Exception as exc:
            logger.exception("Failed to save config")
            self._update_status(f"[red]Save failed:[/] {exc}")

    def _get_widget_value(self, widget: Any, key: str) -> Any:
        """Extract the typed value from a form widget."""
        field_info = self._fields.get(key)
        if field_info is None:
            return None

        if isinstance(widget, Input):
            raw = widget.value.strip()
            widget_type = field_to_widget_type(field_info, field_name=key)

            if widget_type in ("integer",):
                return int(raw) if raw else field_info.default
            if widget_type in ("float",):
                return float(raw) if raw else field_info.default
            # text / password — handle bool & list round-trips
            if field_info.annotation is bool or raw.lower() in ("true", "false"):
                return raw.lower() == "true"
            if get_origin(field_info.annotation) is list:
                if raw:
                    stripped = raw.strip()
                    if stripped.startswith('[') and stripped.endswith(']'):
                        try:
                            return ast.literal_eval(stripped)
                        except (ValueError, SyntaxError):
                            pass
                    items = [item.strip().strip("\"'") for item in stripped.split(",") if item.strip()]
                    return items
                return field_info.default or []
            return raw

        if isinstance(widget, Switch):
            return widget.value

        if isinstance(widget, Select):
            return widget.value

        if isinstance(widget, TextArea):
            raw = widget.text.strip()
            if not raw:
                return field_info.default
            return self._parse_complex(raw, field_info)

        return None

    # ── Open in $EDITOR ─────────────────────────────────────────────────────

    def _on_open_editor(self) -> None:
        """Launch ``$EDITOR`` (or ``vim``) on the loaded config path."""
        if not self._config_path or not os.path.isfile(self._config_path):
            self._update_status("[yellow]No config loaded. Load or save a config first.[/]")
            return

        editor = os.environ.get("EDITOR", "vim")
        try:
            subprocess.call([editor, self._config_path])
        except Exception as exc:
            logger.exception("Failed to launch editor")
            self._update_status(f"[red]Editor failed:[/] {exc}")
            return

        # Re-load from disk after the editor exits
        # Re-load from disk after the editor exits
        self._on_load()

    # ── input validation ───────────────────────────────────────────────────

    def on_input_changed(self, event: Input.Changed) -> None:
        """Validate numeric constraints on the fly."""

        widget = event.input
        ge = getattr(widget, "_field_ge", None)
        le = getattr(widget, "_field_le", None)
        if ge is None and le is None:
            return  # no numeric constraints on this field
        raw = widget.value.strip()
        if not raw:
            return

        try:
            val = int(raw) if widget.type == "integer" else float(raw)
        except (ValueError, TypeError):
            return

        errors: list[str] = []
        if ge is not None and val < ge:
            errors.append(f"minimum {ge}")
        if le is not None and val > le:
            errors.append(f"maximum {le}")

        if errors:
            widget.classes = widget.classes + " field-invalid"
            self._update_status(
                f"[red]Validation:[/] {widget.id or 'field'} — {', '.join(errors)}"
            )
        else:
            # Remove invalid class if previously set
            classes = widget.classes.split()
            if "field-invalid" in classes:
                classes.remove("field-invalid")
                widget.classes = " ".join(classes)
            # Clear transient validation messages from status bar
            # (only if the current status starts with "Validation")
            status_bar = self.query_one("#status-bar", Static)
            if self._last_status.startswith("[red]Validation"):
                self._update_status(
                    f"Loaded: {self._config_path}" if self._config_path else "Ready."
                )

    # ── status helper ──────────────────────────────────────────────────────

    def _update_status(self, message: str) -> None:
        """Update the status bar text (best-effort)."""
        try:
            bar = self.query_one("#status-bar", Static)
            bar.update(message)
            self._last_status = message
        except Exception:
            logger.debug("Failed to update status bar", exc_info=True)
