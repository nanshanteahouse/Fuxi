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
import time
from pathlib import Path
from typing import Any, Optional


# Ensure repo root is on sys.path for core package imports
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

import yaml

import logging
import tempfile

import shutil
from core.paper_registry import load_registry, detect_modality, save_registry, DatasetStatus
from core.dataset_schema import update_pipeline_status

from core.paper_registry_models import ExperimentGroup, _dict_to_exp_group

REPRODUCE_TIMEOUT = 1800  # 30 minutes per GSE pipeline run

GSE_PATTERN = re.compile(r"GSE\d+")

# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────

def _build_experiment_config_path(base_config_path: str, subset_suffix: str) -> str:
    """Build a config path for an experiment group by inserting *subset_suffix*.

    E.g. ``/path/config_GSE123.py`` + ``_T_cell`` -> ``/path/config_GSE123_T_cell.py``
    """
    if base_config_path.endswith(".py"):
        return base_config_path[:-3] + subset_suffix + ".py"
    return base_config_path + subset_suffix



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
    modality: str | None = None,
    experiment_group: ExperimentGroup | None = None,
) -> dict[str, Any]:
    """Run the full pipeline for one GSE dataset as a subprocess.

    Args:
        gse_id:           GEO accession ID (e.g. ``GSE107618``).
        config_path:      Absolute path to the pipeline config file.
        modality:         Explicit modality override. When ``None`` (default),
                          auto-detect from the config file for backward compat.
        experiment_group: Optional experiment group (used by W2.4 config
                          generation path). Stored but does not affect the
                          subprocess call in the basic path.

    Returns:
        A dict with keys:

        * ``status`` — ``"success"``, ``"failed"``, or ``"timeout"``.
        * ``modality`` — detected modality from the config file.
        * ``config_path`` — the config path (echoed back).
        * ``output`` — captured stdout.
        * ``error`` — captured stderr or error message.
        * ``duration_s`` — wall-clock seconds.
    """
    # Determine modality: explicit param overrides config-file detection
    if modality is None:
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
        _result = {
            "status": "success" if result.returncode == 0 else "failed",
            "config_path": config_path,
            "modality": modality,
            "output": result.stdout,
            "error": result.stderr,
            "duration_s": round(elapsed, 1),
        }
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        _result = {
            "status": "timeout",
            "config_path": config_path,
            "modality": modality,
            "output": "",
            "error": f"Timed out after {REPRODUCE_TIMEOUT}s",
            "duration_s": round(elapsed, 1),
        }
    except Exception as exc:
        elapsed = time.time() - t0
        _result = {
            "status": "failed",
            "config_path": config_path,
            "modality": modality,
            "output": "",
            "error": str(exc),
            "duration_s": round(elapsed, 1),
        }

    # Write pipeline status to dataset.yaml on success
    if _result["status"] == "success":
        try:
            config_dir = os.path.dirname(config_path)
            dataset_yaml = os.path.join(config_dir, "dataset.yaml")
            if not os.path.exists(dataset_yaml) and _gse_id:
                data_root = os.environ.get("FUXI_DATA_ROOT", "")
                if data_root:
                    alt_path = os.path.join(data_root, _gse_id, "dataset.yaml")
                    if os.path.exists(alt_path):
                        dataset_yaml = alt_path
            if os.path.exists(dataset_yaml):
                update_pipeline_status(dataset_yaml, modality, "completed")
        except Exception as e:
            _log = logging.getLogger(__name__)
            _log.warning(
                "Failed to update pipeline status in dataset.yaml: %s", e,
            )

    return _result


def _write_pipeline_status(
    registry_path: str,
    gse_id: str,
    result: dict[str, Any],
    paper_dir: str,
) -> None:
    """Incrementally write pipeline_complete status to registry.yaml.

    Uses atomic write (temp file + rename) to prevent corruption from
    partial writes.

    Args:
        registry_path: Path to the registry YAML file.
        gse_id: GEO accession ID.
        result: Result dict from ``_run_pipeline_for_gse()``; may contain
            ``group_name`` for experiment group results.
        paper_dir: Path to the paper directory (not used directly but
            kept for future extensibility).
    """
    logger = logging.getLogger(__name__)

    try:
        registry = load_registry(registry_path)
    except Exception:
        logger.warning("Cannot load registry from %s", registry_path)
        return

    # Find the DatasetEntry matching this GSE ID
    found = False
    for paper in registry.get("papers", []):
        for ds in paper.get("datasets", []):
            if ds.get("gse_id") != gse_id:
                continue
            group_name = result.get("group_name")
            if group_name and ds.get("experiments"):
                for exp in ds["experiments"]:
                    if exp.get("group_name") == group_name:
                        exp["status"] = DatasetStatus.PIPELINE_COMPLETE.value
                        found = True
                        break
            else:
                ds["status"] = DatasetStatus.PIPELINE_COMPLETE.value
                found = True
            break
        if found:
            break

    if not found:
        logger.warning(
            "GSE %s not found in registry at %s; "
            "cannot write pipeline_complete status",
            gse_id, registry_path,
        )
        return

    # Atomic write: write to temp file, then rename
    dir_path = os.path.dirname(registry_path) or "."
    os.makedirs(dir_path, exist_ok=True)

    tmp_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=dir_path, delete=False, mode="w", encoding="utf-8", suffix=".tmp",
        ) as f:
            tmp_path = f.name
            yaml.dump(
                registry, f,
                default_flow_style=False, sort_keys=False, allow_unicode=True,
            )
        os.rename(tmp_path, registry_path)
    except Exception:
        if tmp_path is not None and os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


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

    paper_meta = insights.get("paper_meta", {})

    if registry is None:
        registry = load_registry()

    # Compute registry path for incremental status writes
    registry_path = str(paper_path.parent / "registry.yaml")

    # Look up this paper in the registry by PMID
    paper_entry = None
    for p in registry.get("papers", []):
        if p.get("pmid") == paper_meta.get("pmid"):
            paper_entry = p
            break

    datasets: list[dict[str, Any]] = paper_entry.get("datasets", []) if paper_entry else []

    results: dict[str, dict[str, Any]] = {}
    for ds in datasets:
        gse_id = ds["gse_id"]
        if gse_filter and gse_id != gse_filter:
            continue

        status_enum = ds.get("status")
        config_path = ds.get("config_path", "")
        modality = ds.get("modality", "rna")

        # -- 3-layer experiment group dispatch --
        experiments_data = ds.get("experiments")

        if dry_run:
            if experiments_data:
                for exp_dict in experiments_data:
                    group = _dict_to_exp_group(exp_dict)
                    eg_modalities = (
                        ["rna", "atac"]
                        if group.modality == "multiome"
                        else [group.modality]
                    )
                    for eg_mod in eg_modalities:
                        exp_config_path = (
                            group.config_path
                            if group.config_path
                            else _build_experiment_config_path(config_path, group.subset_suffix)
                        )
                        results[f"{gse_id}_{group.group_name}_{eg_mod}"] = {
                            "status": "dry_run",
                            "config_path": exp_config_path,
                            "modality": eg_mod,
                        }
            else:
                results[gse_id] = {
                    "status": "dry_run",
                    "config_path": config_path,
                    "modality": modality,
                }
            continue

        if status_enum == "config_exists" and config_path:
            if experiments_data:
                for exp_dict in experiments_data:
                    group = _dict_to_exp_group(exp_dict)
                    eg_modalities = (
                        ["rna", "atac"]
                        if group.modality == "multiome"
                        else [group.modality]
                    )
                    for eg_mod in eg_modalities:
                        exp_config_path = (
                            group.config_path
                            if group.config_path
                            else _build_experiment_config_path(config_path, group.subset_suffix)
                        )
                        result = _run_pipeline_for_gse(
                            gse_id, exp_config_path,
                            modality=eg_mod, experiment_group=group,
                        )
                        result.setdefault("config_path", exp_config_path)
                        result.setdefault("modality", eg_mod)
                        result["group_name"] = group.group_name
                        results[f"{gse_id}_{group.group_name}_{eg_mod}"] = result
                        if result["status"] == "success" and not dry_run:
                            _write_pipeline_status(registry_path, gse_id, result, paper_dir)
            else:
                # Original single-config path (backward compatible)
                result = _run_pipeline_for_gse(gse_id, config_path, modality=modality)
                result.setdefault("config_path", config_path)
                result.setdefault("modality", modality)
                results[gse_id] = result
                if result["status"] == "success" and not dry_run:
                    _write_pipeline_status(registry_path, gse_id, result, paper_dir)
                results[gse_id] = result

        elif (
            status_enum == "not_configured"
            or (status_enum == "config_exists" and not config_path)
        ):
            if experiments_data:
                # Config generation for experiment groups
                from core.preprocess.preprocessor import run_preprocess
                from core.preprocess.matrix_loader import (_post_process_config, _resolve_project_dir)

                try:
                    retcode = run_preprocess(
                        gse_id=gse_id,
                        paper_context=insights,
                        force=False,
                        quiet=True,
                    )
                except Exception:
                    retcode = 1

                if retcode == 0:
                    # Base config was generated
                    proj_dir = _resolve_project_dir(modality, gse_id)
                    base_config_path = os.path.join(proj_dir, f'config_{gse_id}.py')

                    for exp_dict in experiments_data:
                        group = _dict_to_exp_group(exp_dict)
                        if group.config_path is None:
                            # Build suffixed config path and generate
                            copied_path = _build_experiment_config_path(
                                base_config_path, group.subset_suffix,
                            )
                            if os.path.exists(base_config_path):
                                shutil.copy2(base_config_path, copied_path)
                                _post_process_config(
                                    copied_path,
                                    paper_context={},
                                    inject={
                                        "sample_keep": group.sample_ids,
                                        "subset_suffix": group.subset_suffix,
                                    }
                                )
                            group.config_path = copied_path

                        eg_modalities = (
                            ["rna", "atac"]
                            if group.modality == "multiome"
                            else [group.modality]
                        )
                        for eg_mod in eg_modalities:
                            results[f"{gse_id}_{group.group_name}_{eg_mod}"] = {
                                "status": "configured",
                                "config_path": group.config_path,
                                "modality": eg_mod,
                                "reason": "Config generated for experiment group",
                            }
                else:
                    # run_preprocess failed — fall through to original behaviour
                    results[gse_id] = {
                        "status": "not_configured",
                        "config_path": config_path,
                        "modality": modality,
                        "reason": "GSE needs config generation (use P2 first)",
                    }
            else:
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
