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

# Lazy re-exports for unittest.mock.patch() compatibility.
# These allow patch("core.run_reproduce.subprocess.run") etc.
# to resolve via pkgutil.resolve_name without eager imports.
import core.ai.caller as ai_caller  # noqa: F401
import core.pipeline.reproduce as run_reproduce  # noqa: F401
