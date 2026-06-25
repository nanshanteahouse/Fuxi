"""rna.steps — cross-modality re-exports for spatial/ATAC pipelines.

The RNA step modules have numeric-starting filenames (e.g. 05_annotate_major.py)
that cannot be imported by name via ``importlib.import_module()``.  We load them
by file path via ``importlib.util`` and re-export the symbols needed by other
modalities.

Usage::

    from rna.steps import unified_annotate, run_ora, run_prerank
"""

import os, sys, importlib.util

# ═════════════════════════════════════════════════════════════════════════════
# Ensure ``rna/`` is on sys.path so that cross-modality callers of
# unified_annotate() can resolve ``from tissue_ontologies import load_kb``
# (tissue_ontologies lives under rna/, not at repo root).
# ═════════════════════════════════════════════════════════════════════════════
_rna_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _rna_dir not in sys.path:
    sys.path.insert(0, _rna_dir)


def _load_step_module(filename: str, mod_name: str):
    """Load a step .py file by path (num-starting filename workaround)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    spec = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── Annotation (Step 05) ────────────────────────────────────────────────────
_annotate = _load_step_module("05_annotate_major.py", "rna.steps._05_annotate_major")

# Name expected by spatial/steps/05_annotate.py (backward compat)
_run_unified_annotation = _annotate.unified_annotate
# Preferred name for new code
unified_annotate = _annotate.unified_annotate


# ── Enrichment (Step 09) ────────────────────────────────────────────────────
_enrichment = _load_step_module("09_enrichment.py", "rna.steps._09_enrichment")

run_ora = _enrichment.run_ora
run_prerank = _enrichment.run_prerank
