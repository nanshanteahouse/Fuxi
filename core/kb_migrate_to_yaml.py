"""
kb_migrate_to_yaml.py — Convert retina KB Python source files to YAML format.

Reads all 10 active Python KB source files from
``rna/tissue_ontologies/retina/sources/``, extracts ``source_meta`` and
``markers``, and writes each to a YAML file following the schema in
``schema.yaml``.

Usage::

    python core/kb_migrate_to_yaml.py

Output
------
    10 files written to ``rna/tissue_ontologies/retina/sources/{name}.yaml``
    plus a round-trip verification summary printed to stdout.
"""

import importlib.util
import os
import sys
import yaml

# ── Repo root ─────────────────────────────────────────────────────────────────
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

# ── Paths ─────────────────────────────────────────────────────────────────────
SOURCES_DIR = os.path.join(
    REPO_ROOT, "rna", "tissue_ontologies", "retina", "sources"
)
SCHEMA_PATH = os.path.join(SOURCES_DIR, "schema.yaml")

# ── Default audit section (NEW* in schema) ────────────────────────────────────
EMPTY_AUDIT = {
    "expression_validated": [],
    "supplement_verified": [],
    "cross_species_validated": False,
    "last_audited": None,
    "flagged": False,
    "notes": "",
}


def load_source_module(filepath: str):
    """Dynamically import a Python source file.

    Returns (module, module_name) or (None, name) on failure.
    """
    module_name = os.path.basename(filepath)[:-3]
    spec = importlib.util.spec_from_file_location(module_name, filepath)
    if spec is None or spec.loader is None:
        print(f"  ⚠ Cannot load spec: {os.path.basename(filepath)}", file=sys.stderr)
        return None, module_name
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as exc:
        print(f"  ⚠ Error executing: {os.path.basename(filepath)} — {exc}", file=sys.stderr)
        return None, module_name
    return mod, module_name


def _normalize_refine(refine_data):
    """Normalise refine entries so they are always dicts.

    A source may store ``refine`` as a plain list of gene names
    (e.g. ``["GENE1", "GENE2"]``) or as a dict ``{"GENE": {...}}``.
    The YAML schema documents both forms; we preserve whatever the
    source defines.
    """
    if isinstance(refine_data, list):
        # List-of-strings form — simple marker list
        return refine_data
    if isinstance(refine_data, dict):
        # Dict form — may contain nested dicts for per-gene annotation
        # or flat string values.  Return as-is.
        return refine_data
    # Fallback
    return {}


def _sort_dict_genes(d):
    """Return a copy of *d* with keys sorted alphabetically.

    Used for ``confirm`` and ``add`` dicts so YAML output is reproducible.
    """
    return dict(sorted(d.items()))


def build_yaml_entry(source_mod):
    """Convert a loaded source module into a YAML-serialisable dict.

    The output dict mirrors the schema structure:
        source_meta, markers (per-type with confirm/add/refine/
        negative_markers/species/synonyms/audit), novel_types,
        expert_rules, conflicts.
    """
    source_meta = getattr(source_mod, "source_meta", {})
    markers_raw = getattr(source_mod, "markers", {})

    # ── Build per-type markers with the full set of fields ────────────
    markers_out = {}
    for cell_type, data in markers_raw.items():
        confirm = _sort_dict_genes(data.get("confirm", {}))
        add = _sort_dict_genes(data.get("add", {}))
        refine_raw = data.get("refine", {})
        refine = _normalize_refine(refine_raw)

        entry = {
            "confirm": confirm,
            "add": add,
            "refine": refine,
            "negative_markers": [],
            "species": [],
            "synonyms": [],
            "audit": dict(EMPTY_AUDIT),  # copy
        }
        markers_out[cell_type] = entry

    # ── Assemble top-level entry ──────────────────────────────────────
    entry = {
        "source_meta": source_meta,
        "markers": markers_out,
        "novel_types": getattr(source_mod, "novel_types", []),
        "expert_rules": getattr(source_mod, "expert_rules", []),
        "conflicts": getattr(source_mod, "conflicts", []),
    }
    return entry


def write_yaml(entry, yaml_path: str):
    """Serialize *entry* to a YAML file (sorted keys, readable style)."""
    with open(yaml_path, "w", encoding="utf-8") as fh:
        yaml.dump(
            entry,
            fh,
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=True,
        )


def discover_source_files():
    """Return sorted list of active ``.py`` source file paths.

    Skips files starting with ``_`` (``_TEMPLATE.py``), ``__init__.py``,
    and files whose first 20 lines contain a ``DISABLED`` marker, matching
    the logic in ``merge.py``.
    """
    sources = []
    for entry in sorted(os.listdir(SOURCES_DIR)):
        if not entry.endswith(".py"):
            continue
        if entry.startswith("_"):
            continue
        if entry == "__init__.py":
            continue

        filepath = os.path.join(SOURCES_DIR, entry)

        # Check for DISABLED marker
        with open(filepath, "r", encoding="utf-8") as fh:
            head = "".join(fh.readline() for _ in range(20))
        if "# DISABLED" in head or ("# NOTE" in head and "DISABLED" in head):
            print(f"  ⏭ Skipping disabled source: {entry}")
            continue

        sources.append(filepath)
    return sources


def run():
    """Main migration routine."""
    sources = discover_source_files()
    if not sources:
        print("No source files found.")
        sys.exit(1)

    print(f"Found {len(sources)} source file(s).\n")

    results = []
    for filepath in sources:
        fname = os.path.basename(filepath)
        name = fname[:-3]  # strip .py
        yaml_path = os.path.join(SOURCES_DIR, f"{name}.yaml")

        print(f"── {fname} ──")

        mod, _ = load_source_module(filepath)
        if mod is None:
            print(f"  ✗ Skipping (import failed)\n")
            continue

        source_meta = getattr(mod, "source_meta", {})
        markers = getattr(mod, "markers", {})
        if not source_meta:
            print(f"  ⚠ Empty source_meta — skipping\n")
            continue

        entry = build_yaml_entry(mod)
        write_yaml(entry, yaml_path)

        n_types = len(markers)
        n_confirm = sum(len(v.get("confirm", {})) for v in markers.values())
        n_add = sum(len(v.get("add", {})) for v in markers.values())

        print(f"  ✓ {name}.yaml — {n_types} cell types, "
              f"{n_confirm} confirm markers, {n_add} add markers")

        results.append((name, filepath, yaml_path, entry))

    # ── Round-trip verification ───────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Round-trip verification (Python dict → YAML → re-loaded)")
    print(f"{'=' * 60}\n")

    try:
        from deepdiff import DeepDiff
    except ImportError:
        print("⚠ deepdiff not installed — skipping round-trip verification.")
        print("  Install: pip install deepdiff")
        results.clear()

    all_ok = True
    for name, _, yaml_path, original_entry in results:
        with open(yaml_path, "r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh)

        if loaded is None:
            print(f"  ✗ {name}.yaml — empty/corrupt")
            all_ok = False
            continue

        diff = DeepDiff(original_entry, loaded, ignore_order=True)
        if diff:
            print(f"  ✗ {name}.yaml — ROUND-TRIP FAILED")
            print(f"    Diffs: {diff}")
            all_ok = False
        else:
            n_types = len(loaded.get("markers", {}))
            print(f"  ✓ {name}.yaml — round-trip OK ({n_types} cell types)")

    # ── Summary ───────────────────────────────────────────────────────
    if all_ok and results:
        total_types = sum(
            len(loaded.get("markers", {}))
            for _, _, yp, _ in results
            if (loaded := yaml.safe_load(open(yp, "r", encoding="utf-8")))
        )
        print(f"\n✅ All {len(results)} YAML files pass round-trip "
              f"({total_types} total cell types).")
    elif not all_ok:
        print(f"\n❌ Some files FAILED round-trip.")
        sys.exit(1)


if __name__ == "__main__":
    run()
