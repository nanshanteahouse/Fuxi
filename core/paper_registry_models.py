#!/usr/bin/env python3
"""
core/paper_registry_models.py — Data models & YAML I/O for the PaperRegistry.

Re-exports:
  DatasetStatus, DatasetEntry, PaperEntry
  load_registry, save_registry
  _dataset_to_dict, _dict_to_dataset, _paper_to_dict, _dict_to_paper
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional

import yaml

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# Enums & Data classes
# ──────────────────────────────────────────────


class DatasetStatus(str, Enum):
    """Status of a GSE dataset relative to pipeline setup."""

    NOT_CONFIGURED = "not_configured"
    """Paper + data exist, but no pipeline config file."""
    CONFIG_EXISTS = "config_exists"
    """Paper + data + config all present."""
    PIPELINE_COMPLETE = "pipeline_complete"
    """Pipeline has been run."""
    DATA_NOT_DOWNLOADED = "data_not_downloaded"
    """Paper mentions GSE but data not in FUXI_DATA_ROOT / projects."""
    DATA_ONLY = "data_only"
    """Data / config exists but no paper links to this GSE."""
    UNKNOWN = "unknown"
    """Cannot determine state."""


@dataclass
class DatasetEntry:
    """A single GSE dataset linked to (or orphaned from) a paper."""

    gse_id: str
    config_path: str = ""
    status: DatasetStatus = DatasetStatus.NOT_CONFIGURED
    modality: str = "rna"  # rna | atac | spatial | multiome
    notes: str = ""
    experiments: Optional[list[ExperimentGroup]] = None


@dataclass
class PaperEntry:
    """A paper with its linked datasets and insight status."""

    pmid: str
    paper_dir: str  # relative path from projects/papers/
    title: str = ""
    journal: str = ""
    year: str = ""
    first_author: str = ""
    doi: str = ""
    datasets: list[DatasetEntry] = field(default_factory=list)
    insights_status: str = "generated"  # generated | pending | failed | no_geo


@dataclass
class ExperimentGroup:
    """An experimental sub-grouping within a dataset."""

    group_name: str
    sample_ids: list[str]
    subset_suffix: str
    modality: str
    status: DatasetStatus
    config_path: Optional[str] = None
    figures: list[str] = field(default_factory=list)


# ──────────────────────────────────────────────
# Serialization helpers
# ──────────────────────────────────────────────


def _dataset_to_dict(ds: DatasetEntry) -> dict[str, Any]:
    data: dict[str, Any] = {
        "gse_id": ds.gse_id,
        "config_path": ds.config_path,
        "status": ds.status.value,
        "modality": ds.modality,
        "notes": ds.notes,
    }
    if ds.experiments:
        data["experiments"] = [_exp_group_to_dict(g) for g in ds.experiments]
    return data


def _dict_to_dataset(data: dict[str, Any]) -> DatasetEntry:
    return DatasetEntry(
        gse_id=data["gse_id"],
        config_path=data.get("config_path", ""),
        status=DatasetStatus(data.get("status", "not_configured")),
        modality=data.get("modality", "rna"),
        notes=data.get("notes", ""),
        experiments=[_dict_to_exp_group(e) for e in data.get("experiments", [])]
        if "experiments" in data
        else None,
    )


def _exp_group_to_dict(group: ExperimentGroup) -> dict[str, Any]:
    data: dict[str, Any] = {
        "group_name": group.group_name,
        "sample_ids": group.sample_ids,
        "subset_suffix": group.subset_suffix,
        "modality": group.modality,
        "status": group.status.value,
    }
    if group.config_path is not None:
        data["config_path"] = group.config_path
    if group.figures:
        data["figures"] = group.figures
    return data


def _dict_to_exp_group(data: dict[str, Any]) -> ExperimentGroup:
    return ExperimentGroup(
        group_name=data["group_name"],
        sample_ids=data["sample_ids"],
        subset_suffix=data["subset_suffix"],
        modality=data["modality"],
        status=DatasetStatus(data["status"]),
        config_path=data.get("config_path"),
        figures=data.get("figures", []),
    )


def _paper_to_dict(entry: PaperEntry) -> dict[str, Any]:
    """Serialize a PaperEntry (incl. datasets) to a plain dict for YAML."""
    d: dict[str, Any] = {
        "pmid": entry.pmid,
        "paper_dir": entry.paper_dir,
        "title": entry.title,
        "journal": entry.journal,
        "year": entry.year,
        "first_author": entry.first_author,
        "doi": entry.doi,
        "insights_status": entry.insights_status,
    }
    d["datasets"] = [_dataset_to_dict(ds) for ds in entry.datasets]
    return d


def _dict_to_paper(data: dict[str, Any]) -> PaperEntry:
    """Reconstruct a PaperEntry from a plain dict (from YAML)."""
    datasets = [_dict_to_dataset(ds) for ds in data.get("datasets", [])]
    return PaperEntry(
        pmid=str(data.get("pmid", "")),
        paper_dir=data.get("paper_dir", ""),
        title=data.get("title", ""),
        journal=data.get("journal", ""),
        year=data.get("year", ""),
        first_author=data.get("first_author", ""),
        doi=data.get("doi", ""),
        datasets=datasets,
        insights_status=data.get("insights_status", "generated"),
    )


# ──────────────────────────────────────────────
# YAML I/O
# ──────────────────────────────────────────────


def load_registry(path: str = "projects/papers/registry.yaml") -> dict[str, Any]:
    """Load the registry YAML file.

    Returns a dict with at least a ``"papers"`` key.
    If the file does not exist or is empty, returns ``{"papers": []}``.
    """
    if not os.path.exists(path):
        return {"papers": []}
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if data is None:
            return {"papers": []}
        return data
    except Exception:
        logger.exception("Failed to load registry from %s", path)
        return {"papers": []}


def save_registry(
    registry: dict[str, Any],
    path: str = "projects/papers/registry.yaml",
) -> None:
    """Save the registry dict to a YAML file.

    Uses block style (``default_flow_style=False``) for human readability,
    matching the convention in ``core/dataset_schema.py``.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(registry, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
