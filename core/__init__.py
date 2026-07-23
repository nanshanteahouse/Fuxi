#!/usr/bin/env python3
"""
core/ — Fuxi (伏羲) shared infrastructure for all modalities.

Sub-packages:
  - ai/           LLM caller + prompt templates
  - annotation/   Cell-type annotation engine + standardizer + marker scoring
  - cluster/      Grid-search clustering + parameter evaluation
  - config/       Unified Pydantic Config + dataset schema
  - interaction/  Cell-cell interaction (CCI) utilities
  - kb/           Tissue knowledge base (markers, adjacency, pathways)
  - paper/        Paper insights, registry, converter, cross-paper analysis
  - pipeline/     Pipeline runner, anatomy, enrichment, GRN, reproducibility
  - preprocess/   Format detection → config generation
  - utils/        I/O, logging, path resolution, validation, performance
"""

# True lazy re-exports — use __getattr__ to defer import until first access.
# This avoids triggering the runpy warning:
#   "'core.paper.registry' found in sys.modules after import of package
#   'core.paper', but prior to execution of 'core.paper.registry'"
# which was caused by eager `import core.pipeline.reproduce` during `import core`.

_LAZY_MODULES: dict[str, str] = {
    "ai_caller": "core.ai.caller",
    "run_reproduce": "core.pipeline.reproduce",
}


def __getattr__(name: str):
    if name in _LAZY_MODULES:
        import importlib

        mod = importlib.import_module(_LAZY_MODULES[name])
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'core' has no attribute '{name}'")
