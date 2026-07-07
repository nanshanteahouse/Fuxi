#!/usr/bin/env python3
"""
core/paper_registry.py — YAML-based registry linking papers to GSE datasets.

Provides:
  - Data models: PaperEntry, DatasetEntry, DatasetStatus (re-exported)
  - YAML I/O: load_registry, save_registry (re-exported)
  - Scanner: build_registry() builds a cross-reference from projects/papers/
    and projects/{rna,atac,spatial}/ into a single registry.yaml

Usage:
    from core.paper_registry import build_registry, save_registry

    registry = build_registry()
    save_registry(registry)
"""
from __future__ import annotations

import os
import sys
import logging
import re
from pathlib import Path
from typing import Any

# Ensure repo root is on sys.path for direct CLI invocation
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from core.paper_registry_models import (
    DatasetEntry,
    DatasetStatus,
    PaperEntry,
    _dataset_to_dict,
    _dict_to_dataset,
    _dict_to_paper,
    _paper_to_dict,
    load_registry,
    save_registry,
)

__all__ = [
    "DatasetStatus",
    "DatasetEntry",
    "PaperEntry",
    "load_registry",
    "save_registry",
    "build_registry",
    "detect_modality",
    "_scan_insights_yamls",
    "_scan_project_dirs",
    "_find_data_only_entries",
]

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────
# Modality detection
# ──────────────────────────────────────────────


def detect_modality(config_path: str) -> str:
    """Detect modality from a pipeline config file.

    Uses a regex scan for ``CFG.modality = "..."`` rather than importing
    the file (avoids side-effects). Returns the modality string or
    ``"unknown"`` if detection fails.
    """
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            content = f.read()
    except Exception:
        return "unknown"

    match = re.search(r"""CFG\.modality\s*=\s*["'](\w+)["']""", content)
    return match.group(1) if match else "unknown"


# ──────────────────────────────────────────────
# Scanners
# ──────────────────────────────────────────────

VALID_MODALITIES = ("rna", "atac", "spatial")


def _scan_insights_yamls(papers_dir: str) -> list[dict[str, Any]]:
    """Scan ``papers_dir/*/insights.yaml`` and extract paper metadata.

    Returns a list of plain dicts with keys:
      pmid, paper_dir, title, journal, year, first_author, doi, geo_ids.

    Papers without a ``pmid`` field or with corrupt YAML are skipped
    (a warning is logged).
    """
    papers: list[dict[str, Any]] = []
    root = Path(papers_dir)
    if not root.exists():
        return papers

    for subdir in sorted(root.iterdir()):
        if not subdir.is_dir():
            continue
        yaml_path = subdir / "insights.yaml"
        if not yaml_path.exists():
            continue

        try:
            data = _safe_load_yaml(yaml_path)
        except Exception as exc:
            logger.warning("Skipping %s: corrupt YAML (%s)", subdir.name, exc)
            continue

        if data is None:
            logger.warning("Skipping %s: empty insights.yaml", subdir.name)
            continue

        paper_meta = data.get("paper_meta", {})
        pmid = paper_meta.get("pmid")
        if not pmid:
            logger.warning("Skipping %s: no pmid in paper_meta", subdir.name)
            continue

        # Extract & deduplicate geo_ids
        data_access = data.get("data_access", {})
        raw_geo: list[str] = data_access.get("geo_ids") or []
        unique_geo = _dedup_list(raw_geo)

        papers.append({
            "pmid": str(pmid),
            "paper_dir": subdir.name,
            "title": paper_meta.get("title", ""),
            "journal": paper_meta.get("journal", ""),
            "year": paper_meta.get("year", ""),
            "first_author": paper_meta.get("first_author", ""),
            "doi": paper_meta.get("doi", ""),
            "geo_ids": unique_geo,
        })

    return papers


def _scan_project_dirs(projects_dir: str) -> list[dict[str, Any]]:
    """Scan ``projects/{rna,atac,spatial}/`` for GSE subdirectories.

    For each GSE directory, looks for a ``config_<GSE_ID>.py`` file
    (not variant configs like ``config_hard.py`` or ``config_mad.py``).

    Returns a list of dicts with keys:
      gse_id, config_path, modality.
    """
    entries: list[dict[str, Any]] = []
    root = Path(projects_dir)
    if not root.exists():
        return entries

    for modality in VALID_MODALITIES:
        mod_path = root / modality
        if not mod_path.exists():
            continue

        for child in sorted(mod_path.iterdir()):
            if not child.is_dir():
                continue

            gse_id = child.name

            # Prefer exact match: config_<GSE_ID>.py
            exact = child / f"config_{gse_id}.py"
            if exact.exists():
                config_path = str(exact)
            else:
                matches = sorted(child.glob("config_GSE*.py"))
                config_path = str(matches[0]) if matches else ""

            entries.append({
                "gse_id": gse_id,
                "config_path": config_path,
                "modality": modality,
            })

    return entries


def _find_data_only_entries(
    scanned_gses: list[dict[str, Any]],
    all_paper_gse_ids: set[str],
) -> list[DatasetEntry]:
    """Find GSE configs that have no paper linking to them.

    Only GSE entries with a non-empty ``config_path`` are considered
    (bare directories without a config file are not marked data-only).
    """
    data_only: list[DatasetEntry] = []
    for gse in scanned_gses:
        if gse["gse_id"] not in all_paper_gse_ids and gse["config_path"]:
            data_only.append(DatasetEntry(
                gse_id=gse["gse_id"],
                config_path=gse["config_path"],
                status=DatasetStatus.DATA_ONLY,
                modality=gse["modality"],
                notes="GSE data has config but no paper linked",
            ))
    return data_only


# ──────────────────────────────────────────────
# Small helpers
# ──────────────────────────────────────────────


def _safe_load_yaml(path: Path) -> Any:
    """Load a YAML file, raising on failure."""
    import yaml  # already a dependency
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _dedup_list(items: list[str]) -> list[str]:
    """Deduplicate a list while preserving order."""
    seen: set[str] = set()
    result: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result


# ──────────────────────────────────────────────
# Registry builder
# ──────────────────────────────────────────────


def build_registry(
    papers_dir: str = "projects/papers",
    projects_dir: str = "projects",
) -> dict[str, Any]:
    """Build the full registry by cross-referencing papers and GSE datasets.

    Steps:
      1. Scan ``papers_dir/*/insights.yaml`` for paper metadata & geo_ids.
      2. Scan ``projects_dir/{rna,atac,spatial}/`` for GSE directories & configs.
      3. Cross-reference: for each paper's geo_ids, look up matching GSE dirs.
      4. Determine per-dataset status (``config_exists`` / ``not_configured`` /
         ``data_not_downloaded``).
      5. Detect GSE configs not linked to any paper (``data_only``).

    Returns a dict ready for ``save_registry()``:

    .. code-block:: yaml

        papers:
          - pmid: "31653841"
            paper_dir: "2019_Menon_Nature_Com_..."
            ...
            datasets:
              - gse_id: "GSE107618"
                ...
        data_only_datasets:
          - gse_id: "GSE999999"
            ...
    """
    # Phase 1 — scan papers
    paper_infos = _scan_insights_yamls(papers_dir)

    # Phase 2 — scan project GSE directories
    gse_entries = _scan_project_dirs(projects_dir)

    # Index GSE entries by GSE ID (one ID may appear in multiple modalities)
    gse_by_id: dict[str, list[dict[str, Any]]] = {}
    for gse in gse_entries:
        gse_by_id.setdefault(gse["gse_id"], []).append(gse)

    # Collect all GSE IDs referenced by papers
    all_paper_gse_ids: set[str] = set()

    # Phase 3 — build PaperEntry list
    paper_entries: list[PaperEntry] = []
    for pinfo in paper_infos:
        datasets: list[DatasetEntry] = []

        for geo_id in pinfo["geo_ids"]:
            all_paper_gse_ids.add(geo_id)

            if geo_id in gse_by_id:
                for gse in gse_by_id[geo_id]:
                    config_path = gse["config_path"]
                    status = DatasetStatus.CONFIG_EXISTS if config_path else DatasetStatus.NOT_CONFIGURED

                    datasets.append(DatasetEntry(
                        gse_id=geo_id,
                        config_path=config_path,
                        status=status,
                        modality=gse["modality"],
                        notes="",
                    ))
            else:
                datasets.append(DatasetEntry(
                    gse_id=geo_id,
                    status=DatasetStatus.DATA_NOT_DOWNLOADED,
                    modality="rna",
                    notes="GSE dataset directory not found in projects/",
                ))

        insights_status = "no_geo" if not pinfo["geo_ids"] else "generated"

        paper_entries.append(PaperEntry(
            pmid=pinfo["pmid"],
            paper_dir=pinfo["paper_dir"],
            title=pinfo["title"],
            journal=pinfo["journal"],
            year=pinfo["year"],
            first_author=pinfo["first_author"],
            doi=pinfo["doi"],
            datasets=datasets,
            insights_status=insights_status,
        ))

    # Phase 4 — detect data-only GSE configs
    data_only = _find_data_only_entries(gse_entries, all_paper_gse_ids)

    # Phase 5 — assemble final registry
    registry: dict[str, Any] = {
        "papers": [_paper_to_dict(pe) for pe in paper_entries],
    }
    if data_only:
        registry["data_only_datasets"] = [_dataset_to_dict(ds) for ds in data_only]

    return registry


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────


def main() -> None:
    """CLI entry point for building/verifying PaperRegistry."""
    import argparse
    parser = argparse.ArgumentParser(description="PaperRegistry — paper ↔ GSE ↔ config linkage")
    parser.add_argument("--build", action="store_true", help="Build registry.yaml from projects/")
    parser.add_argument("--verify", action="store_true", help="Verify registry.yaml consistency")
    parser.add_argument("--dry-run", action="store_true", help="Preview without writing")
    parser.add_argument("--papers-dir", default="projects/papers", help="Papers directory")
    parser.add_argument("--projects-dir", default="projects", help="Projects root directory")
    parser.add_argument("--output", default="projects/papers/registry.yaml", help="Output registry path")
    args = parser.parse_args()

    if args.build:
        registry = build_registry(papers_dir=args.papers_dir, projects_dir=args.projects_dir)
        n_papers = len(registry.get("papers", []))
        n_datasets = sum(len(p.get("datasets", [])) for p in registry.get("papers", []))
        n_data_only = len(registry.get("data_only_datasets", []))
        print(f"Registry built: {n_papers} papers, {n_datasets} dataset links, {n_data_only} data-only")
        if not args.dry_run:
            save_registry(registry, args.output)
            print(f"Written: {args.output}")
    elif args.verify:
        try:
            registry = load_registry(args.output)
            n_papers = len(registry.get("papers", []))
            issues = 0
            for p in registry.get("papers", []):
                datasets = p.get("datasets", [])
                for d in datasets:
                    if d["status"] == "config_exists" and d.get("config_path"):
                        import os
                        if not os.path.exists(d["config_path"]):
                            print(f"  WARNING: {d['gse_id']} config not found: {d['config_path']}")
                            issues += 1
            if issues:
                print(f"Verify: {n_papers} papers, {issues} issue(s) found")
            else:
                print(f"Verify: {n_papers} papers, all consistent")
        except FileNotFoundError:
            print(f"Registry not found: {args.output}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
