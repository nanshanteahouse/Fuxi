"""Template scaffold — render committed starter-config YAML files from specs.

The committed files under ``templates/config_templates/`` are *generated
artifacts*: their only consumers are humans (manual config starting point)
and regression tests.  Runtime config generation never reads them — it
assembles from the specs directly (``core/preprocess/config_specs.py``).

Usage::

    python -m core.config scaffold            # regenerate templates/
    python -m core.config scaffold --check    # exit 1 if templates are stale
"""

from __future__ import annotations

import os
import re
import sys
from typing import Any, Tuple, Union

from core.preprocess.config_specs import FormatSpec, SpecComment, SpecField, materialized_specs

# ═══════════════════════════════════════════════════════════════════
# Scalar rendering
# ═══════════════════════════════════════════════════════════════════

_BARE_RE = re.compile(r"^[A-Za-z0-9_./+%@()-]+$")
_RESERVED = {"true", "false", "null", "~", "yes", "no", "on", "off"}


def _yaml_scalar(v: Any) -> str:
    """Render a scalar value as a YAML token (quoting when needed)."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if v is None:
        return ""
    if isinstance(v, list):
        return "[" + ", ".join(_yaml_scalar(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{}"
    s = str(v)
    if s == "":
        return '""'
    if _BARE_RE.match(s) and s.lower() not in _RESERVED:
        return s
    return f'"{s}"'


# ═══════════════════════════════════════════════════════════════════
# Template rendering
# ═══════════════════════════════════════════════════════════════════

GENERATED_HEADER = (
    "# ══════════════════════════════════════════════════════════════\n"
    "# AUTO-GENERATED from core/preprocess/config_specs.py\n"
    "# Do NOT edit by hand — run: python -m core.config scaffold\n"
    "# ══════════════════════════════════════════════════════════════\n"
)


def _render_items(items: Tuple[Union[SpecComment, SpecField], ...]) -> str:
    """Render spec items to YAML text with comments.

    ``SpecField`` indentation derives from the dotted path depth; parent
    keys are emitted when a field enters a new nested block.
    """
    out: list[str] = []
    open_stack: list[str] = []

    for item in items:
        if isinstance(item, SpecComment):
            indent = "  " * item.indent
            for line in item.lines:
                if line == "":
                    out.append("")
                elif line.startswith("#"):
                    out.append(indent + line)
                else:
                    out.append(indent + "# " + line)
            continue

        parts = item.path.split(".")
        parents = parts[:-1]
        leaf = parts[-1]

        if parents != open_stack:
            common = 0
            while (
                common < len(open_stack)
                and common < len(parents)
                and open_stack[common] == parents[common]
            ):
                common += 1
            for depth in range(common, len(parents)):
                out.append("  " * depth + f"{parents[depth]}:")
            open_stack = list(parents)

        indent = "  " * len(open_stack)
        if item.placeholder:
            value = '"{{' + item.placeholder + '}}"'
        else:
            value = _yaml_scalar(item.value)
        line = f"{indent}{leaf}: {value}"
        if item.comment:
            line += "   # " + item.comment
        out.append(line)

    return "\n".join(out) + "\n"


def render_template_text(spec: FormatSpec) -> str:
    """Render the full committed template text for *spec*."""
    return GENERATED_HEADER + _render_items(spec.items)


# ═══════════════════════════════════════════════════════════════════
# Scaffold I/O
# ═══════════════════════════════════════════════════════════════════


def _template_dir() -> str:
    """Absolute path to ``templates/config_templates/``."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(
        os.path.dirname(os.path.dirname(this_dir)), "templates", "config_templates"
    )


def write_templates(directory: str | None = None) -> list[str]:
    """Regenerate all committed template files; return written paths."""
    directory = directory or _template_dir()
    os.makedirs(directory, exist_ok=True)
    written: list[str] = []
    for spec in materialized_specs():
        path = os.path.join(directory, spec.template_name)
        with open(path, "w", encoding="utf-8") as f:
            f.write(render_template_text(spec))
        written.append(path)
    return written


def check_templates(directory: str | None = None) -> list[str]:
    """Return template paths whose committed content differs from the
    freshly rendered text (empty = everything is up to date)."""
    directory = directory or _template_dir()
    stale: list[str] = []
    for spec in materialized_specs():
        path = os.path.join(directory, spec.template_name)
        if not os.path.isfile(path):
            stale.append(f"{path} (missing)")
            continue
        with open(path, encoding="utf-8") as f:
            if f.read() != render_template_text(spec):
                stale.append(path)
    return stale


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    # Sanity guard: schema fields exist before rendering anything
    from core.preprocess.config_specs import validate_specs

    errors = validate_specs()
    if errors:
        print("config_specs validation failed:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 2

    check = "--check" in argv
    if check:
        stale = check_templates()
        if stale:
            for p in stale:
                print(f"[STALE] {p}")
            print("Run: python -m core.config scaffold")
            return 1
        print("templates up to date")
        return 0

    for p in write_templates():
        print(f"[WRITTEN] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
