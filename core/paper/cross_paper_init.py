#!/usr/bin/env python3
"""
Cross-Paper Comparison Framework
=================================
Generic module for cross-paper single-cell analysis.
Supports arbitrary gene sets and pathway signatures, not restricted to RA.

Usage:
    from cross_paper.analyzer import CrossPaperAnalyzer, DatasetEntry

    # From YAML configs
    analyzer = CrossPaperAnalyzer.from_yaml(
        pathway_yaml="cross_paper/pathway_config.yaml",
        registry_yaml="cross_paper/dataset_registry.yaml",
    )
    analyzer.run()
    analyzer.compare_conditions()
"""

from __future__ import annotations

import os
import sys
from typing import Any

import yaml

# Add repo root so core/ is importable
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)


def load_yaml(path: str) -> dict[str, Any]:
    """Load a YAML file from a project-relative or absolute path."""
    if not os.path.isabs(path):
        path = os.path.join(_repo_root, path)
    with open(path, "r") as f:
        return yaml.safe_load(f)


def get_default_pathway_config() -> dict:
    """Return the default pathway gene set definitions."""
    return load_yaml(os.path.join(os.path.dirname(__file__), "pathway_config.yaml"))


def get_default_registry() -> dict:
    """Return the default dataset registry."""
    return load_yaml(os.path.join(os.path.dirname(__file__), "dataset_registry.yaml"))


__all__ = [
    "load_yaml",
    "get_default_pathway_config",
    "get_default_registry",
]
