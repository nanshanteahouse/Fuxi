"""rna.steps — cross-modality re-exports for spatial pipeline.

The RNA step modules have numeric-starting filenames (e.g. 09_enrichment.py)
that cannot be imported by name via ``importlib.import_module()``.  We load them
by file path via ``importlib.util`` and re-export the symbols needed by other
modalities.

Note: ``unified_annotate`` has been moved to ``core.annotation.engine``.
Use ``from core.annotation.engine import run_unified_annotation`` directly.

Usage::

    from rna.steps import run_ora, run_prerank
"""

import importlib.util
import os


def _load_step_module(filename: str, mod_name: str):
    """Load a step .py file by path (num-starting filename workaround)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    spec = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── Enrichment (Step 09) ────────────────────────────────────────────────────
_enrichment = _load_step_module("09_enrichment.py", "rna.steps._09_enrichment")

run_ora = _enrichment.run_ora
run_prerank = _enrichment.run_prerank
