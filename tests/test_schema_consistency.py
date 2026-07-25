"""Regression test for ClusteringSettings schema ↔ production code consistency.

Prevents the NF2-NF5 class of bugs where schema field names/defaults
don't match what production code actually reads via getattr / direct
attribute access.  Three checks:

1. Every ClusteringSettings field appears in ≥1 production .py file
   (dead-field detection — fields only read by YAML go in the allowlist).

2. Every ``getattr(cfg.clustering, "<field>", …)`` name exists in
   ``ClusteringSettings.model_fields`` (ghost-field detection).

3. Known ``getattr(…, fallback)`` fallback values match the schema
   defaults (drift detection).
"""

from __future__ import annotations

import ast
import pathlib
import re
import textwrap
from typing import Any

from core.config.schema import ClusteringSettings

# ── Paths ──────────────────────────────────────────────────────────────

_REPO = pathlib.Path(__file__).resolve().parent.parent
_SCAN_ROOTS = [
    _REPO / "core",
    _REPO / "rna",
    _REPO / "spatial",
    _REPO / "atac",
    _REPO / "bulk",
]
_EXCLUDED_SUBDIRS = {"tests", "__pycache__", ".git", ".mypy_cache", ".pytest_cache"}
# Exclude the schema file itself + conftest + __init__ (only re-exports).
_EXCLUDED_FILES = {"schema.py", "conftest.py", "__init__.py"}

# ── Test 1 allowlist ───────────────────────────────────────────────────
# Fields that are legitimately *not* referenced by name in production
# Python source (e.g. used only through YAML validation or Pydantic
# model_post_init).  Do NOT flag these as dead.

_TEST1_ALLOWLIST: set[str] = {
    "param_grid_resolutions",
    "param_grid_min_dist",
    "param_grid_spread",
    "best_resolution",
    "best_n_neighbors",
    "umap_color_by_batch",
    "umap_plot_mode",
    "batch_key_override",
    "umap_plot_max_cells",
    "kb_pass_rate",
    "kb_max_markers",
    "kb_mode",
    "kb_annotatable_rate",
    "n_pcs",
    "pca_random_state",
    "umap_random_state",
    "leiden_flavor",
    "umap_spread",
    # Fields that truly do *not* appear in any production .py file
    # as of 2026-07:
    "multi_metric_granularity_cv_threshold",
    "multi_metric_granularity_min_clusters",
    "param_grid_n_neighbors_adaptive",
}

# ── Helpers ────────────────────────────────────────────────────────────


def _collect_production_py_files() -> list[pathlib.Path]:
    """Walk scan roots, returning non-excluded .py files."""
    files: list[pathlib.Path] = []
    for root in _SCAN_ROOTS:
        if not root.is_dir():
            continue
        for py_file in root.rglob("*.py"):
            # Skip excluded subdirs anywhere in the path.
            parts = set(py_file.relative_to(_REPO).parts)
            if parts & _EXCLUDED_SUBDIRS:
                continue
            if py_file.name in _EXCLUDED_FILES:
                continue
            files.append(py_file)
    return files


def _safe_read(file_path: pathlib.Path) -> str:
    """Read file content, returning '' on any error."""
    try:
        return file_path.read_text(encoding="utf-8")
    except Exception:
        return ""


# ── Test 2: getattr extraction ─────────────────────────────────────────

# Matches: getattr(cfg.clustering, "field_name"  or  getattr(cfg.clustering, 'field_name'
_GETATTR_RE = re.compile(r'getattr\(\s*cfg\.clustering\s*,\s*["\']([a-zA-Z_][a-zA-Z0-9_]*)["\']')


def _extract_getattr_fields(
    files: list[pathlib.Path],
) -> dict[str, list[tuple[pathlib.Path, int]]]:
    """Return {field_name: [(file, line_no), ...]} for all getattr(cfg.clustering, ...) calls."""
    found: dict[str, list[tuple[pathlib.Path, int]]] = {}
    for fpath in files:
        content = _safe_read(fpath)
        for lineno, line in enumerate(content.splitlines(), start=1):
            match = _GETATTR_RE.search(line)
            if match:
                name = match.group(1)
                found.setdefault(name, []).append((fpath, lineno))
    return found


# ── Test 3: function parameter default resolution ──────────────────────


def _get_enclosing_function_default(
    file_path: pathlib.Path, target_line: int, param_name: str
) -> Any | None:
    """Parse *file_path* with AST; find the function that contains
    *target_line*; return the default value of *param_name*, or
    ``_UNRESOLVED`` if the default cannot be statically determined."""
    content = _safe_read(file_path)
    if not content:
        return _UNRESOLVED
    try:
        tree = ast.parse(content, filename=str(file_path))
    except SyntaxError:
        return _UNRESOLVED

    # Walk to find the innermost function containing target_line.
    target: ast.FunctionDef | ast.AsyncFunctionDef | None = None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.lineno <= target_line <= node.end_lineno:  # type: ignore[attr-defined]
                # Keep the innermost (tightest enclosing) function.
                if target is None or node.lineno >= target.lineno:
                    target = node

    if target is None:
        return _UNRESOLVED

    # The function's positional args form the first N elements of
    # target.args.args; defaults align to the *last* N args.
    args = target.args
    defaults = args.defaults or []
    # Number of positional args without defaults = total args - len(defaults)
    offset = len(args.args) - len(defaults)
    for i, arg in enumerate(args.args):
        if arg.arg == param_name and i >= offset:
            default_node = defaults[i - offset]
            return _eval_ast_literal(default_node)
    return _UNRESOLVED


class _Unresolved:
    """Sentinel for values that cannot be statically extracted from source."""

    def __repr__(self) -> str:
        return "<unresolved>"


_UNRESOLVED = _Unresolved()


def _eval_ast_literal(node: ast.expr) -> Any:
    """Best-effort evaluation of AST constant/literal nodes."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        if isinstance(node.operand, ast.Constant):
            return -node.operand.value
    if isinstance(node, ast.List):
        return [_eval_ast_literal(elt) for elt in node.elts]
    if isinstance(node, ast.Tuple):
        return tuple(_eval_ast_literal(elt) for elt in node.elts)
    if isinstance(node, ast.Dict):
        return {_eval_ast_literal(k): _eval_ast_literal(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.Name):
        # Common sentinels.
        if node.id == "None":
            return None
        if node.id == "True":
            return True
        if node.id == "False":
            return False
        # Variable reference — return name string for downstream resolution.
        return node.id


def _extract_getattr_fallback(file_path: pathlib.Path, line_no: int, field_name: str) -> Any:
    """Return the third argument (fallback) of the getattr call on *line_no*.

    Returns ``_UNRESOLVED`` on any parsing failure."""
    content = _safe_read(file_path)
    if not content:
        return _UNRESOLVED
    lines = content.splitlines()
    if line_no > len(lines):
        return _UNRESOLVED
    line = lines[line_no - 1]

    # pattern: getattr(cfg.clustering, "field", <FALLBACK>)
    pat = re.compile(
        rf'getattr\(\s*cfg\.clustering\s*,\s*["\']{re.escape(field_name)}["\']\s*,\s*(.+?)\s*\)'
    )
    m = pat.search(line)
    if not m:
        return _UNRESOLVED
    fallback_expr = m.group(1)

    # Try parsing as a Python expression.
    try:
        node = ast.parse(fallback_expr.strip(), mode="eval")
        return _eval_ast_literal(node.body)  # type: ignore[attr-defined]
    except SyntaxError:
        return fallback_expr.strip()  # raw string for diagnostics


# ═══════════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════════


class TestClusteringSchemaConsistency:
    """Verify ClusteringSettings fields match production code usage."""

    # ── Test 1: every field is read somewhere ──────────────────────────

    def test_all_clustering_fields_are_read(self) -> None:
        """Every ClusteringSettings field name appears in ≥1 production .py file."""
        schema_fields = set(ClusteringSettings.model_fields.keys())
        py_files = _collect_production_py_files()

        # Build a set of fields found in any .py file.
        found: set[str] = set()
        for fpath in py_files:
            content = _safe_read(fpath)
            for field_name in schema_fields - found:
                if field_name in content:
                    found.add(field_name)

        dead = schema_fields - found - _TEST1_ALLOWLIST
        assert not dead, textwrap.dedent(f"""\
            Dead ClusteringSettings field(s) — not referenced in any production .py file:

            {_fmt_field_list(dead)}

            These fields exist in the schema but aren't read by any production code.
            Options:
              • Remove from ClusteringSettings if truly unused.
              • Add to _TEST1_ALLOWLIST if the field is consumed through
                YAML config resolution, Pydantic model_post_init, or other
                indirect paths.
        """)

    # ── Test 2: no ghost fields ───────────────────────────────────────

    def test_no_ghost_clustering_fields(self) -> None:
        """No getattr(cfg.clustering, ...) references a non-existent field."""
        schema_fields = set(ClusteringSettings.model_fields.keys())
        py_files = _collect_production_py_files()
        getattr_usage = _extract_getattr_fields(py_files)

        ghosts: dict[str, list[tuple[pathlib.Path, int]]] = {}
        for name, locations in getattr_usage.items():
            if name not in schema_fields:
                ghosts[name] = locations

        if not ghosts:
            return  # all clean

        msg_lines = ["Ghost ClusteringSettings field(s) accessed via getattr():", ""]
        for name, locations in sorted(ghosts.items()):
            msg_lines.append(f"  • {name!r}  (used in {len(locations)} location(s)):")
            for fpath, lineno in locations:
                rel = fpath.relative_to(_REPO)
                msg_lines.append(f"      {rel}:{lineno}")
        msg_lines.append("")
        msg_lines.append(
            "Either add these fields to ClusteringSettings or remove the ghost getattr calls."
        )
        assert not ghosts, "\n".join(msg_lines)

    # ── Test 3: defaults match code fallbacks ──────────────────────────

    # Pairs: (field_name, file_relpath, line_no)
    _KNOWN_FALLBACK_PAIRS: list[tuple[str, str, int]] = [
        ("multi_metric_coverage_ratio_threshold", "core/cluster/evaluation.py", 1189),
        ("stability_n_seeds", "core/cluster/evaluation.py", 1188),
        ("stability_leiden_n_iterations", "core/cluster/evaluation.py", 301),
        ("leiden_flavor", "core/cluster/evaluation.py", 1190),
        ("umap_selection_metric", "rna/steps/04_cluster_umap.py", 567),
    ]

    def test_clustering_defaults_match_code_fallbacks(self) -> None:
        """Schema defaults must equal getattr(..., fallback) values."""
        schema_fields = ClusteringSettings.model_fields
        mismatches: list[str] = []

        for field_name, rel_path, line_no in self._KNOWN_FALLBACK_PAIRS:
            file_path = _REPO / rel_path
            schema_default = schema_fields[field_name].default

            code_fallback = _extract_getattr_fallback(file_path, line_no, field_name)

            # If the fallback is a variable name (e.g. ``n_iterations``),
            # resolve it to the enclosing function's parameter default.
            if isinstance(code_fallback, str) and code_fallback not in ("None", ""):
                resolved = _get_enclosing_function_default(file_path, line_no, code_fallback)
                if resolved is not _UNRESOLVED:
                    code_fallback = resolved

            if code_fallback is _UNRESOLVED:
                mismatches.append(
                    f"  {field_name}: schema={schema_default!r}, "
                    f"code_fallback=<unresolved> — {rel_path}:{line_no}"
                )
            elif schema_default != code_fallback:
                mismatches.append(
                    f"  {field_name}: schema={schema_default!r}, "
                    f"code_fallback={code_fallback!r} — {rel_path}:{line_no}"
                )

        assert not mismatches, (
            "Schema default / code fallback mismatch(es):\n\n"
            + "\n".join(mismatches)
            + "\n\nUpdate the schema default, the code fallback, or both "
            "so they are identical."
        )


# ── Formatting ─────────────────────────────────────────────────────────


def _fmt_field_list(fields: set[str]) -> str:
    """Return a sorted, indented list of field names."""
    return "\n".join(f"    • {f}" for f in sorted(fields))
