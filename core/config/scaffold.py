"""Config scaffold — render starter-config YAML from schema specs on demand.

There are no committed template files: the specs in
``core/preprocess/config_specs.py`` are the single source, and this module
renders them whenever a human needs a starting point (manual dataset setup).

Usage::

    python -m core.config scaffold --list                    # available formats
    python -m core.config scaffold --format 10X_h5           # render to stdout
    python -m core.config scaffold --format 10X_mtx --out config_GSE.yaml
"""

from __future__ import annotations

import argparse
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
    """Render the full starter-config text for *spec*."""
    return GENERATED_HEADER + _render_items(spec.items)


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m core.config scaffold",
        description="Render starter config YAML from core/preprocess/config_specs.py specs.",
    )
    parser.add_argument("--list", action="store_true", help="list available format keys and exit")
    parser.add_argument(
        "--format", metavar="KEY", help="render a single format spec (default: all)"
    )
    parser.add_argument(
        "--out",
        metavar="PATH",
        help="write rendered YAML to PATH instead of stdout (requires --format)",
    )
    args = parser.parse_args(argv)

    # Sanity guard: schema fields exist before rendering anything
    from core.preprocess.config_specs import validate_specs

    errors = validate_specs()
    if errors:
        print("config_specs validation failed:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 2

    specs = materialized_specs()
    if args.list:
        for spec in specs:
            print(f"{spec.key:<16} {spec.modality:<9} {spec.data_format}")
        return 0

    if args.format is not None:
        specs = [s for s in specs if s.key == args.format]
        if not specs:
            print(f"unknown format key: {args.format}", file=sys.stderr)
            print(
                "available keys: " + ", ".join(s.key for s in materialized_specs()),
                file=sys.stderr,
            )
            return 1

    if args.out is not None:
        if len(specs) != 1:
            print("--out requires --format (one spec)", file=sys.stderr)
            return 1
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(render_template_text(specs[0]))
        print(f"[WRITTEN] {args.out}")
        return 0

    for spec in specs:
        sys.stdout.write(render_template_text(spec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
