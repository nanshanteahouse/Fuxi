#!/usr/bin/env python3
"""
core/run_reproduce.py — Paper result reproduction manager (P3).

Orchestrates the full pipeline execution for a paper's datasets:
  1. Load paper metadata from insights.yaml
  2. Cross-reference with the paper registry for config status
  3. For configured GSEs: run the pipeline via subprocess
  4. For unconfigured GSEs: report status as ``not_configured``
  5. Collect & return reproduction status per GSE

CLI Usage::

    python core/run_reproduce.py projects/papers/<paper_dir>
    python core/run_reproduce.py --all
    python core/run_reproduce.py projects/papers/<paper_dir> --dry-run
    python core/run_reproduce.py --gse GSE107618
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import os
import time
import argparse
from pathlib import Path
from typing import Any, Optional

import re

# Ensure repo root is on sys.path for core package imports
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

import yaml

from core.paper_registry import load_registry, detect_modality

REPRODUCE_TIMEOUT = 1800  # 30 minutes per GSE pipeline run

GSE_PATTERN = re.compile(r"GSE\d+")

# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────


def _detect_modality(config_path: str) -> str:
    """Read modality from an existing pipeline config file.

    Reuses ``core.paper_registry.detect_modality`` which does a regex scan
    for ``CFG.modality = "..."`` — safe (no import side-effects) and fast.
    Returns the modality string or ``"unknown"`` on failure.
    """
    return detect_modality(config_path)


def _extract_geo_ids(insights: dict, raw_text: str = "") -> list[str]:
    """Extract GEO accession IDs from the insights dict + optional raw text.

    Priority:
      1. ``data_access.geo_ids`` from the parsed YAML dict.
      2. If empty, regex scan the raw YAML text for ``GSE\\d+`` patterns.

    Returns a deduplicated list preserving first-occurrence order.
    """
    geo_ids: list[str] = insights.get("data_access", {}).get("geo_ids") or []

    if not geo_ids and raw_text:
        # Regex fallback — scan the raw file text for GSE patterns
        matches = GSE_PATTERN.findall(raw_text)
        # Deduplicate while preserving order
        seen: set[str] = set()
        geo_ids = [m for m in matches if not (m in seen or seen.add(m))]  # type: ignore[func-returns-value]

    return geo_ids


# ──────────────────────────────────────────────
# Pipeline execution
# ──────────────────────────────────────────────


def _run_pipeline_for_gse(
    _gse_id: str, config_path: str,
) -> dict[str, Any]:
    """Run the full pipeline for one GSE dataset as a subprocess.

    Args:
        gse_id:       GEO accession ID (e.g. ``GSE107618``).
        config_path:  Absolute path to the pipeline config file.

    Returns:
        A dict with keys:

        * ``status`` — ``"success"``, ``"failed"``, or ``"timeout"``.
        * ``modality`` — detected modality from the config file.
        * ``config_path`` — the config path (echoed back).
        * ``output`` — captured stdout.
        * ``error`` — captured stderr or error message.
        * ``duration_s`` — wall-clock seconds.
    """
    modality = _detect_modality(config_path)
    if modality == "unknown":
        return {
            "status": "failed",
            "config_path": config_path,
            "modality": "unknown",
            "output": "",
            "error": f"Cannot detect modality from config: {config_path}",
            "duration_s": 0.0,
        }

    cmd = [
        sys.executable,
        "core/run_pipeline.py",
        "--config",
        config_path,
        "--modality",
        modality,
    ]
    t0 = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=REPRODUCE_TIMEOUT,
        )
        elapsed = time.time() - t0
        return {
            "status": "success" if result.returncode == 0 else "failed",
            "config_path": config_path,
            "modality": modality,
            "output": result.stdout,
            "error": result.stderr,
            "duration_s": round(elapsed, 1),
        }
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        return {
            "status": "timeout",
            "config_path": config_path,
            "modality": modality,
            "output": "",
            "error": f"Timed out after {REPRODUCE_TIMEOUT}s",
            "duration_s": round(elapsed, 1),
        }
    except Exception as exc:
        elapsed = time.time() - t0
        return {
            "status": "failed",
            "config_path": config_path,
            "modality": modality,
            "output": "",
            "error": str(exc),
            "duration_s": round(elapsed, 1),
        }


# ──────────────────────────────────────────────
# Orchestration
# ──────────────────────────────────────────────


def run_reproduce(
    paper_dir: str,
    registry: Optional[dict] = None,
    gse_filter: Optional[str] = None,
    dry_run: bool = False,
) -> dict[str, dict[str, Any]]:
    """Reproduce (re-run) all datasets for a single paper.

    Workflow
    --------
    1. Load ``insights.yaml`` from *paper_dir*.
    2. Extract GEO accession IDs (dict + regex fallback).
    3. Look up the paper in the registry (matched by PMID).
    4. For each registered dataset entry:

       * ``config_exists`` with a valid path → run the pipeline.
       * ``not_configured`` or ``config_exists`` without path → skip with reason.
       * ``data_not_downloaded`` → skip with reason.
       * Any other status → skip with reason.

    5. If *dry_run* is ``True``, report the plan without executing anything.

    Args:
        paper_dir:   Path to the paper directory containing ``insights.yaml``.
        registry:    Pre-loaded registry dict (or ``None`` to auto-load).
        gse_filter:  Optional GSE ID to restrict processing to one dataset.
        dry_run:     When ``True``, preview the reproduction plan; do **not**
                     run any subprocess or preprocess step.

    Returns:
        A dict mapping ``gse_id → result_dict`` where each result has:

        * ``status`` — ``"success"`` | ``"failed"`` | ``"timeout"`` |
          ``"skipped"`` | ``"not_configured"`` | ``"dry_run"``.
        * ``config_path`` — path to the config file (if any).
        * ``modality`` — detected / registry modality.
        * ``reason`` — human-readable explanation (absent for ``dry_run``
          and pipeline-run results).

    Raises:
        ValueError: If ``insights.yaml`` does not exist or is empty/corrupt.
    """
    paper_path = Path(paper_dir)
    insights_path = paper_path / "insights.yaml"
    if not insights_path.exists():
        raise ValueError(f"No insights.yaml in {paper_dir}")

    with open(insights_path, encoding="utf-8") as f:
        raw_text = f.read()

    insights = yaml.safe_load(raw_text)
    if insights is None:
        raise ValueError(f"Empty or invalid insights.yaml in {paper_dir}")

    _geo_ids = _extract_geo_ids(insights, raw_text)  # extracted for potential future cross-validation
    paper_meta = insights.get("paper_meta", {})

    if registry is None:
        registry = load_registry()

    # Look up this paper in the registry by PMID
    paper_entry = None
    for p in registry.get("papers", []):
        if p.get("pmid") == paper_meta.get("pmid"):
            paper_entry = p
            break

    datasets: list[dict[str, Any]] = paper_entry.get("datasets", []) if paper_entry else []

    # Build paper context (species default + marker features from first figure)
    paper_context: dict[str, Any] = {
        "species": paper_meta.get("species", "human"),
        "marker_dict": {},
        "features": [],
    }
    for fig in insights.get("figures", []):
        if fig.get("features"):
            paper_context["marker_dict"] = {"extracted": fig["features"]}
            paper_context["features"] = fig["features"]
            break

    results: dict[str, dict[str, Any]] = {}
    for ds in datasets:
        gse_id = ds["gse_id"]
        if gse_filter and gse_id != gse_filter:
            continue

        status_enum = ds.get("status")
        config_path = ds.get("config_path", "")
        modality = ds.get("modality", "rna")

        if dry_run:
            results[gse_id] = {
                "status": "dry_run",
                "config_path": config_path,
                "modality": modality,
            }
            continue

        if status_enum == "config_exists" and config_path:
            result = _run_pipeline_for_gse(gse_id, config_path)
            # Ensure documented fields are always present
            result.setdefault("config_path", config_path)
            result.setdefault("modality", modality)
            results[gse_id] = result

        elif (
            status_enum == "not_configured"
            or (status_enum == "config_exists" and not config_path)
        ):
            results[gse_id] = {
                "status": "not_configured",
                "config_path": config_path,
                "modality": modality,
                "reason": "GSE needs config generation (use P2 first)",
            }

        elif status_enum == "data_not_downloaded":
            results[gse_id] = {
                "status": "skipped",
                "config_path": config_path,
                "modality": modality,
                "reason": "GSE data not downloaded",
            }

        else:
            results[gse_id] = {
                "status": "skipped",
                "config_path": config_path,
                "modality": modality,
                "reason": f"status={status_enum}",
            }

    return results


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fuxi (伏羲) — Paper reproduction manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python core/run_reproduce.py projects/papers/<paper_dir>\n"
            "  python core/run_reproduce.py --all\n"
            "  python core/run_reproduce.py --all --dry-run\n"
            "  python core/run_reproduce.py --gse GSE107618\n"
        ),
    )
    parser.add_argument(
        "paper_dir",
        nargs="?",
        help="Path to the paper directory (containing insights.yaml).",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Reproduce all papers in the registry.",
    )
    parser.add_argument(
        "--gse",
        type=str,
        default=None,
        help="Target a specific GSE ID (ignores others). "
        "Mutually exclusive with --all.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the reproduction plan without executing anything.",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    if args.all and args.paper_dir:
        parser.error("--all and paper_dir are mutually exclusive")
    if args.all and args.gse:
        parser.error("--all and --gse are mutually exclusive")

    registry = load_registry()
    papers_dir = Path("projects/papers")

    if args.all:
        papers: list[dict[str, Any]] = registry.get("papers", [])
    elif args.paper_dir:
        # paper_dir can be an absolute path or relative; treat as-is
        papers = [{"paper_dir": args.paper_dir}]
    else:
        parser.print_help()
        return

    if not papers:
        print("No papers found.")
        return

    all_results: dict[str, dict[str, dict[str, Any]]] = {}
    for p in papers:
        pd = p.get("paper_dir", "")
        paper_path = papers_dir / pd if pd and not os.path.isabs(pd) else Path(pd)
        try:
            results = run_reproduce(
                str(paper_path),
                registry=registry,
                gse_filter=args.gse,
                dry_run=args.dry_run,
            )
        except ValueError as exc:
            print(f"[SKIP] {pd}: {exc}")
            continue

        all_results[pd] = results

    # ── Print results table ───────────────────────────────────────────────
    header = f"{'Paper':35s} {'GSE':15s} {'Status':20s} {'Modality':10s} {'Config':45s}"
    sep = "-" * 125
    print(f"\n{header}")
    print(sep)
    for paper_name, results in all_results.items():
        if not results:
            print(f"{paper_name[:33]:35s} {'—':15s} {'no datasets':20s} {'':10s} {'':45s}")
            continue
        for gse_id, r in results.items():
            print(
                f"{paper_name[:33]:35s} "
                f"{gse_id:15s} "
                f"{r.get('status', '?'):20s} "
                f"{r.get('modality', '?'):10s} "
                f"{r.get('config_path', '')[:43]:45s}"
            )
    print()


if __name__ == "__main__":
    main()
