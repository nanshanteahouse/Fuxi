#!/usr/bin/env python3
"""
migrate_configs.py — .py → .yaml config migration script

Converts existing Fuxi project configs from Python dataclass assignments
(CFG.field = value) to YAML format consumable by the new Pydantic Config
model.  Uses AST-based parsing only — never exec() or import the config.

Usage:
    python core/migrate_configs.py --dry-run              # Preview all
    python core/migrate_configs.py --all                  # Migrate all
    python core/migrate_configs.py --gse GSE137846        # Single
    python core/migrate_configs.py --gse GSE137846 --dry-run
"""

# allow: SIZE_OK — standalone migration script (AST parser + YAML builder + CLI)
# at ~540 pure LOC; splitting would scatter tightly coupled logic across files.

import argparse
import ast
import os
import re
import shutil
import sys
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.config_field_map import FIELD_MAP, TOPIC_NAMES


# ── Topic section labels for YAML comment headers ──
TOPIC_LABELS: Dict[str, str] = {
    "data_input":     "Data input",
    "sample_meta":    "Sample / stage metadata",
    "qc":             "Quality control",
    "scrublet":       "Doublet detection (Scrublet)",
    "normalization":  "Normalization & cell-cycle",
    "hvg":            "Highly variable genes",
    "pca":            "Principal component analysis",
    "harmony":        "Harmony batch correction",
    "clustering":     "Clustering & UMAP",
    "marker":         "Cell-type markers & annotation",
    "de":             "Differential expression",
    "trajectory":     "Pseudotime / trajectory",
    "enrichment":     "Gene-set enrichment",
    "grn":            "Gene regulatory network",
    "cci":            "Cell-cell interaction",
    "downsample":     "Downsampling / subset",
    "spatial":        "Spatial transcriptomics",
    "atac":           "ATAC-specific",
    "execution":      "Execution environment",
    "ai":             "AI / LLM configuration",
}

# Ordered list of root-level fields (non-topic) for YAML output.
ROOT_FIELDS: List[str] = [
    "modality", "tissue", "species", "tissue_maturity",
    "expression_type", "data_format",
    "data_dir", "results_dir", "h5ad_dir", "figure_dir",
    "table_dir", "log_dir", "project_dir",
    "h5ad_compression", "h5ad_tempdir", "cleanup_intermediates",
    "perf_monitoring",
    "rna_h5ad", "rna_ref",
    "rna_marker_top_n", "rna_marker_pval_threshold", "rna_marker_logfc_min",
    "tissue_kb", "tissue_ontology",
    "target_class", "target_order",
]

TOPIC_ORDER: List[str] = list(TOPIC_NAMES.keys())


# ═══════════════════════════════════════════════════════════════════════
#  File discovery
# ═══════════════════════════════════════════════════════════════════════


def find_config_files(gse_id: Optional[str] = None) -> List[Tuple[str, str]]:
    """Find config_GSE*.py files under projects/.

    Skips historical backup directories (v1/, v2/).
    Returns list of (absolute_path, gse_id) tuples sorted by path.
    """
    results: List[Tuple[str, str]] = []
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    projects_dir = os.path.join(repo_root, "projects")

    if not os.path.isdir(projects_dir):
        print(f"Warning: projects dir not found at {projects_dir}")
        return results

    for root, _dirs, files in os.walk(projects_dir):
        rel = os.path.relpath(root, projects_dir)
        # Skip historical backup dirs
        if rel.startswith("rna/v1") or rel.startswith("rna/v2"):
            continue

        for fname in files:
            if not fname.startswith("config_GSE") or not fname.endswith(".py"):
                continue
            m = re.search(r"GSE(\d+)", fname)
            if not m:
                continue
            found = f"GSE{m.group(1)}"
            if gse_id and found != gse_id:
                continue
            results.append((os.path.join(root, fname), found))

    return sorted(results, key=lambda x: x[0])


# ═══════════════════════════════════════════════════════════════════════
#  AST‑based config parser
# ═══════════════════════════════════════════════════════════════════════


class _ConfigVisitor(ast.NodeVisitor):
    """AST visitor that extracts CFG.<field> = <value> assignments.

    Also tracks intermediate variable definitions so patterns like::

        DATA_DIR = os.path.join(data_root(), 'GSE138002')
        CFG.matrix_file = os.path.join(DATA_DIR, 'file.mtx.gz')

    can be resolved by substituting the variable definition.
    """

    def __init__(self, gse_id: str) -> None:
        super().__init__()
        self.gse_id = gse_id

        # flat field name → Python value
        self.fields: Dict[str, Any] = {}
        # flat field name → warning string or None
        self.warnings: Dict[str, Optional[str]] = {}
        # variable name → AST node (for intermediate var tracking)
        self._variables: Dict[str, ast.AST] = {}
        # set of seen flat names (for duplicate detection)
        self._seen: set = set()

    # ── dispatch ────────────────────────────────────────────────────

    def visit_Assign(self, node: ast.Assign) -> None:
        if len(node.targets) != 1:
            return

        target = node.targets[0]

        # Track intermediate variables (all-caps convention in configs)
        if isinstance(target, ast.Name) and not target.id.startswith("_"):
            self._variables[target.id] = node.value

        # CFG.x.y.z = value
        chain = self._get_cfg_chain(target)
        if chain is None:
            return

        flat_name = chain[-1]  # leaf attribute name
        value, warning = self._resolve_expr(node.value)

        # Last assignment wins
        self.fields[flat_name] = value
        self.warnings[flat_name] = warning

        if flat_name in self._seen:
            pass  # duplicate — last wins silently
        self._seen.add(flat_name)

    # ── CFG attribute chain extraction ──────────────────────────────

    @staticmethod
    def _get_cfg_chain(node: ast.AST) -> Optional[List[str]]:
        """Walk CFG.a.b.c → ['a', 'b', 'c']; return None if not CFG."""
        parts: List[str] = []
        cur = node
        while isinstance(cur, ast.Attribute):
            parts.append(cur.attr)
            cur = cur.value
        if isinstance(cur, ast.Name) and cur.id == "CFG":
            parts.reverse()
            return parts
        return None

    # ── expression resolution ───────────────────────────────────────

    def _resolve_expr(self, node: ast.AST) -> Tuple[Any, Optional[str]]:
        """Evaluate an AST expression -> (python_value, warning)."""
        # 1. Try ast.literal_eval (covers literals, lists, dicts)
        try:
            return ast.literal_eval(node), None
        except (ValueError, TypeError):
            pass

        # 2. os.path.join(data_root(), ...) patterns
        result = self._resolve_data_root_join(node)
        if result is not None:
            return result

        # 3. os.environ.get(...) patterns
        result = self._resolve_env_get(node)
        if result is not None:
            return result

        # 4. os.path.dirname(os.path.abspath(__file__)) -> project_dir
        if self._is_dirname_abspath_file(node):
            return "# AUTO (project_dir)", "# TODO: verify path"

        # 5. Variable reference (DATA_DIR etc.)
        if isinstance(node, ast.Name) and node.id in self._variables:
            return self._resolve_expr(self._variables[node.id])

        # 6. Nested os.path.join inside a variable
        if isinstance(node, ast.Call) and self._is_join_call(node):
            return self._resolve_join_with_vars(node)

        # 7. Fallback: unresolvable
        return f"# UNRESOLVED: {ast.dump(node)}", "# TODO: verify path"

    # ── pattern matchers ────────────────────────────────────────────

    @staticmethod
    def _is_data_root_call(node: ast.AST) -> bool:
        return (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "data_root"
            and not node.args
        )

    @staticmethod
    def _is_join_call(node: ast.AST) -> bool:
        return (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "join"
        )

    @staticmethod
    def _is_dirname_abspath_file(node: ast.AST) -> bool:
        """Match os.path.dirname(os.path.abspath(__file__))."""
        if not (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "dirname"
        ):
            return False
        inner = node.args[0] if node.args else None
        if not (
            isinstance(inner, ast.Call)
            and isinstance(inner.func, ast.Attribute)
            and inner.func.attr == "abspath"
        ):
            return False
        src = inner.args[0] if inner.args else None
        if isinstance(src, ast.Name) and src.id == "__file__":
            return True
        if isinstance(src, ast.Attribute) and src.attr == "__file__":
            return True
        return False

    @staticmethod
    def _call_args(node: ast.AST) -> list:
        """Safely get args from a Call node."""
        if isinstance(node, ast.Call):
            return list(node.args)
        return []

    # ── data_root( ) join resolver ──────────────────────────────────

    def _check_is_data_root_variable(self, node: ast.AST) -> bool:
        """Check if an AST node ultimately resolves to a data_root() call.

        Handles:
        - Direct: ``data_root()``
        - Variable: ``DATA_DIR`` where it was defined as ``data_root()``
        - Chain: ``DATA_DIR`` where ``DATA_DIR = os.path.join(data_root(), 'GSE_ID')``
        """
        if self._is_data_root_call(node):
            return True

        if isinstance(node, ast.Name) and node.id in self._variables:
            var_node = self._variables[node.id]
            if self._is_data_root_call(var_node):
                return True
            if isinstance(var_node, ast.Call) and self._is_join_call(var_node) and var_node.args:
                return self._check_is_data_root_variable(var_node.args[0])

        return False

    def _resolve_data_root_join(
        self, node: ast.AST
    ) -> Optional[Tuple[Any, Optional[str]]]:
        """Resolve os.path.join(data_root(), 'GSE_ID', 'rest/of/path').

        Returns (relative_tail, None) or None if the pattern does not match.
        """
        if not self._is_join_call(node):
            return None
        args = self._call_args(node)
        if len(args) < 2:
            return None

        # First arg must be data_root() or a variable resolving to it
        first = args[0]
        if not self._check_is_data_root_variable(first):
            return None

        # Extract string literal parts from remaining args
        path_parts: List[str] = []
        for arg in args[1:]:
            try:
                part = ast.literal_eval(arg)
            except (ValueError, TypeError):
                # arg might be a variable reference
                if isinstance(arg, ast.Name) and arg.id in self._variables:
                    sub_val, _ = self._resolve_expr(self._variables[arg.id])
                    if isinstance(sub_val, str):
                        part = sub_val
                    else:
                        return None
                else:
                    return None
            if not isinstance(part, str):
                return None
            path_parts.append(part)

        if not path_parts:
            return None

        # The first part after data_root should be the GSE ID — strip it
        # Use exact match to avoid stripping filenames like GSE138002_Final.mtx.gz
        if path_parts and re.match(r'^GSE\d+$', path_parts[0]):
            path_parts = path_parts[1:]

        tail = os.path.join(*path_parts) if path_parts else ""
        return tail, None

    # ── os.path.join with variable resolution ───────────────────────

    def _resolve_join_with_vars(
        self, node: ast.AST
    ) -> Tuple[Any, Optional[str]]:
        """Resolve a join whose args may contain variable references."""
        args = self._call_args(node)
        resolved_args: List[str] = []

        for arg in args:
            try:
                val = ast.literal_eval(arg)
                resolved_args.append(str(val))
            except (ValueError, TypeError):
                if isinstance(arg, ast.Name) and arg.id in self._variables:
                    sub_val, _ = self._resolve_expr(self._variables[arg.id])
                    if isinstance(sub_val, str):
                        resolved_args.append(sub_val)
                    else:
                        return (
                            f"# UNRESOLVED: {ast.dump(node)}",
                            "# TODO: verify path",
                        )
                else:
                    return (
                        f"# UNRESOLVED: {ast.dump(node)}",
                        "# TODO: verify path",
                    )

        # Strip leading data_root-like prefix
        if resolved_args and resolved_args[0].startswith("FUXI_DATA_ROOT"):
            resolved_args = resolved_args[1:]
        if resolved_args and re.match(r'^GSE\d+$', resolved_args[0]):
            resolved_args = resolved_args[1:]

        if resolved_args:
            return os.path.join(*resolved_args), None
        return "", None

    # ── os.environ.get resolver ─────────────────────────────────────

    def _resolve_env_get(
        self, node: ast.AST
    ) -> Optional[Tuple[Any, Optional[str]]]:
        """Resolve os.environ.get('VAR', 'default') -> (default, warning)."""
        if not (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "get"
        ):
            return None
        val = node.func.value
        if not (
            isinstance(val, ast.Attribute)
            and val.attr == "environ"
            and isinstance(val.value, ast.Name)
            and val.value.id == "os"
        ):
            return None

        if len(node.args) < 2:
            return None

        try:
            default_val = ast.literal_eval(node.args[1])
        except (ValueError, TypeError):
            return None

        try:
            var_name = ast.literal_eval(node.args[0])
        except (ValueError, TypeError):
            var_name = "?"

        return default_val, f"# from env {var_name}"


# ═══════════════════════════════════════════════════════════════════════
#  YAML structure builder
# ═══════════════════════════════════════════════════════════════════════


def build_yaml_structure(
    fields: Dict[str, Any],
    warnings: Dict[str, Optional[str]],
) -> Tuple[Dict[str, Any], Dict[str, Optional[str]], List[str]]:
    """Map flat field names to nested topic-based YAML structure.

    Returns (nested_dict, nested_warnings, unknown_fields).

    *nested_dict* uses topic keys for sections and plain keys for root-level
    fields.  *nested_warnings* uses dotted path keys (e.g. ``qc.min_genes``).
    """
    nested: Dict[str, Any] = {}
    nested_warnings: Dict[str, Optional[str]] = {}
    unknown: List[str] = []

    for flat_name, value in fields.items():
        if flat_name not in FIELD_MAP:
            unknown.append(flat_name)
            continue

        target_path = FIELD_MAP[flat_name]
        parts = target_path.split(".")

        if len(parts) == 1:
            # Root-level field
            nested[parts[0]] = value
            if warnings.get(flat_name):
                nested_warnings[parts[0]] = warnings[flat_name]
        else:
            topic = parts[0]
            leaf = ".".join(parts[1:])
            if topic not in nested:
                nested[topic] = {}
            elif not isinstance(nested[topic], dict):
                nested[topic] = {}
            nested[topic][leaf] = value
            full_path = f"{topic}.{leaf}"
            if warnings.get(flat_name):
                nested_warnings[full_path] = warnings[flat_name]

    return nested, nested_warnings, unknown


# ═══════════════════════════════════════════════════════════════════════
#  YAML serialisation (manual — full comment control)
# ═══════════════════════════════════════════════════════════════════════


def _yaml_scalar(value: Any) -> str:
    """Render a single Python value as a YAML scalar string."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str):
        if not value:
            return "''"
        if value in {
            "true", "false", "null", "yes", "no", "on", "off",
            "True", "False", "None", "Yes", "No", "On", "Off",
        }:
            return repr(value)
        needs_quoting = (
            value.startswith("#")
            or value.startswith("*")
            or value.startswith("&")
            or value.startswith("!")
            or value.startswith("{")
            or value.startswith("[")
            or value.startswith("'")
            or value.startswith('"')
            or value.startswith(" ")
            or value.endswith(" ")
            or ":" in value
            or "#" in value
        )
        return repr(value) if needs_quoting else value
    if isinstance(value, list):
        items = [_yaml_scalar(i) for i in value]
        return f"[{', '.join(items)}]"
    if isinstance(value, dict):
        items = [f"{k}: {_yaml_scalar(v)}" for k, v in value.items()]
        return f"{{{', '.join(items)}}}"
    return repr(value)


def _emit_line(
    lines: List[str],
    key: str,
    value: Any,
    warning: Optional[str],
    indent: int = 0,
) -> None:
    """Emit one ``key: value`` line, with an optional warning comment."""
    prefix = " " * indent
    scalar = _yaml_scalar(value)
    if warning:
        lines.append(f"{prefix}# {key}: {scalar}  {warning}")
        lines.append(f"{prefix}# {key}: {scalar}")
    else:
        lines.append(f"{prefix}{key}: {scalar}")


def to_yaml_text(
    data: Dict[str, Any],
    warnings: Dict[str, Optional[str]],
    unknown_fields: List[str],
) -> str:
    """Produce the complete YAML text with section comment headers."""
    lines: List[str] = []

    # Separate root vs topic sections
    root_vals: Dict[str, Any] = {}
    topic_vals: Dict[str, Dict[str, Any]] = {}
    for key, value in data.items():
        if isinstance(value, dict) and key in TOPIC_NAMES:
            topic_vals[key] = dict(value)
        else:
            root_vals[key] = value

    # ── Root-level fields in canonical order ──
    for field in ROOT_FIELDS:
        if field in root_vals:
            _emit_line(lines, field, root_vals.pop(field), warnings.get(field))
    for key, value in root_vals.items():
        _emit_line(lines, key, value, warnings.get(key))

    # ── Topic sections ──
    for topic in TOPIC_ORDER:
        if topic not in topic_vals:
            continue
        label = TOPIC_LABELS.get(topic, topic.replace("_", " ").title())
        lines.append("")
        lines.append(f"# ── {label} ──")
        lines.append(f"{topic}:")
        for key, value in topic_vals[topic].items():
            full = f"{topic}.{key}"
            _emit_line(lines, key, value, warnings.get(full), indent=2)

    # ── Migration notes footer ──
    notes: List[str] = []
    if unknown_fields:
        notes.append(
            f"Unknown fields (not in FIELD_MAP): {', '.join(unknown_fields)}"
        )
    if notes:
        lines.append("")
        lines.append("# ── Migration notes ──")
        for note in notes:
            lines.append(f"# {note}")

    lines.append("")
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
#  Migration orchestration
# ═══════════════════════════════════════════════════════════════════════


def migrate_one_config(filepath: str, dry_run: bool = False) -> bool:
    """Parse one .py config, write .yaml, rename .py to .py.bak.

    Returns True on success.
    """
    m = re.search(r"GSE(\d+)", filepath)
    if not m:
        print(f"  \u2717 Cannot extract GSE ID from {filepath}")
        return False
    gse_id = f"GSE{m.group(1)}"

    print(f"\n{'─' * 60}")
    print(f"  Config: {filepath}")
    print(f"  GSE ID: {gse_id}")

    # Read
    try:
        with open(filepath, encoding="utf-8") as fh:
            source = fh.read()
    except OSError as exc:
        print(f"  \u2717 Read error: {exc}")
        return False

    # Parse AST
    try:
        tree = ast.parse(source, filename=filepath)
    except SyntaxError as exc:
        print(f"  \u2717 Syntax error: {exc}")
        return False

    # Extract fields
    visitor = _ConfigVisitor(gse_id)
    visitor.visit(tree)

    if not visitor.fields:
        print(f"  \u26a0 No CFG assignments found")
        return False

    # Build YAML structure
    nested, nested_warnings, unknown = build_yaml_structure(
        visitor.fields, visitor.warnings
    )

    yaml_text = to_yaml_text(nested, nested_warnings, unknown)

    output_path = filepath.replace(".py", ".yaml")

    if dry_run:
        print(f"  \u2192 Would write: {output_path}")
        print(f"  \u2192 Would backup: {filepath} \u2192 {filepath}.bak")
        if unknown:
            print(f"  \u26a0 Unknown fields: {unknown}")
        print(f"\n  {'─' * 40}")
        print(f"  Generated YAML preview:")
        print(f"  {'─' * 40}")
        lines = yaml_text.splitlines()
        shown = 0
        for line in lines:
            if shown >= 60:
                remaining = len(lines) - 60
                print(f"  \u2026 ({remaining} more line{'s' if remaining != 1 else ''})")
                break
            print(f"  {line}")
            shown += 1
        print(f"  {'─' * 40}")
    else:
        try:
            with open(output_path, "w", encoding="utf-8") as fh:
                fh.write(yaml_text)
            print(f"  \u2713 Wrote {output_path}")
        except OSError as exc:
            print(f"  \u2717 Write error: {exc}")
            return False

        try:
            bak = filepath + ".bak"
            shutil.move(filepath, bak)
            print(f"  \u2713 Backed up \u2192 {bak}")
        except OSError as exc:
            print(f"  \u26a0 Backup failed: {exc}")

    return True


# ═══════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Migrate Fuxi project configs from .py to .yaml format"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--all",
        action="store_true",
        help="Migrate all config_GSE*.py files under projects/",
    )
    group.add_argument(
        "--gse",
        type=str,
        help="Migrate a single dataset (e.g. GSE137846)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview migration without writing anything",
    )

    args = parser.parse_args()

    if args.all:
        configs = find_config_files()
    else:
        gse = args.gse if args.gse.startswith("GSE") else f"GSE{args.gse}"
        configs = find_config_files(gse)

    if not configs:
        print("No config files found matching the given criteria.")
        sys.exit(0 if args.dry_run else 1)

    print(f"Found {len(configs)} config file(s):")
    for fp, _gse in configs:
        rel = os.path.relpath(
            fp, os.path.dirname(os.path.dirname(__file__))
        )
        print(f"  {_gse}: {rel}")

    success = 0
    failure = 0
    for fp, _gse in configs:
        ok = migrate_one_config(fp, dry_run=args.dry_run)
        if ok:
            success += 1
        else:
            failure += 1

    print(f"\n{'═' * 60}")
    if args.dry_run:
        print(f"DRY RUN: {success} config(s) would be migrated")
    else:
        print(f"Migration complete: {success} succeeded, {failure} failed")

    sys.exit(0 if failure == 0 else 1)


if __name__ == "__main__":
    main()
