"""MCP tools for pipeline execution (download, preprocess, run, insights).

These are the "active" tools — they trigger real work (downloading data,
generating configs, running pipeline steps). Long-running operations return
structured results so the AI agent can decide on next steps.

Uses a hybrid approach:
- download_dataset: calls Python API directly (download_gse returns structured dict)
- preprocess/run/insights: subprocess (CLI already handles orchestration)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from typing import Any

from mcp.server import MCPServer


def _repo_root() -> str:
    """Return the absolute path to the repo root."""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _run_subprocess(
    cmd: list[str],
    timeout: int = 600,
    description: str = "",
) -> dict[str, Any]:
    """Run a subprocess command and return structured result.

    Args:
        cmd: Command + args as list (e.g., [sys.executable, "script.py", "--gse", "X"])
        timeout: Max seconds before killing the process.
        description: Human-readable label for error messages.

    Returns:
        Dict with keys: ok (bool), stdout (str), stderr (str), exit_code (int),
        elapsed_s (float), timed_out (bool).
    """
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=_repo_root(),
        )
        elapsed = time.monotonic() - t0
        return {
            "ok": proc.returncode == 0,
            "stdout": proc.stdout[-5000:],  # Keep last 5KB to avoid context blowup
            "stderr": proc.stderr[-2000:] if proc.stderr else "",
            "exit_code": proc.returncode,
            "elapsed_s": round(elapsed, 1),
            "timed_out": False,
            "command": " ".join(cmd),
        }
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        return {
            "ok": False,
            "stdout": "",
            "stderr": f"Timeout after {timeout}s",
            "exit_code": -1,
            "elapsed_s": round(elapsed, 1),
            "timed_out": True,
            "command": " ".join(cmd),
        }
    except Exception as exc:
        elapsed = time.monotonic() - t0
        return {
            "ok": False,
            "stdout": "",
            "stderr": str(exc),
            "exit_code": -1,
            "elapsed_s": round(elapsed, 1),
            "timed_out": False,
            "command": " ".join(cmd),
        }


def register_execution_tools(server: MCPServer) -> None:
    """Register all execution MCP tools on the given server."""

    # ── download_dataset ────────────────────────────────────────────────

    @server.tool()
    async def download_dataset(gse: str, data_root: str = "", dry_run: bool = False) -> str:
        """Download a GEO dataset to FUXI_DATA_ROOT using wget/curl.

        Downloads supplementary files from NCBI GEO. Can take minutes to hours
        depending on dataset size. Call registry_status() first to check if
        data is already downloaded.

        Args:
            gse: GEO accession ID (e.g., "GSE123456")
            data_root: Override FUXI_DATA_ROOT (uses env var if empty)
            dry_run: If true, report what would happen without downloading

        Returns:
            JSON with download status, file counts, and next step suggestion.
        """
        gse = gse.upper()

        # Resolve data root
        root = data_root or os.environ.get("FUXI_DATA_ROOT", "")
        if not root:
            return json.dumps(
                {
                    "error": "FUXI_DATA_ROOT not set",
                    "action": "Set FUXI_DATA_ROOT environment variable or pass data_root parameter",
                },
                indent=2,
            )

        dest_dir = os.path.join(root, gse)

        try:
            from core.geo_downloader import download_gse, update_registry_after_download

            result = download_gse(
                gse_id=gse,
                dest_dir=dest_dir,
                dry_run=dry_run,
                skip_existing=True,
                fetch_meta=True,
                quiet=False,
            )

            # Update registry on success
            if not dry_run and result["failed"] == 0:
                try:
                    update_registry_after_download(gse, dry_run=False)
                    result["registry_updated"] = True
                except Exception:
                    result["registry_updated"] = False

            # Add next step
            if result["failed"] == 0 and not dry_run:
                result["next_step"] = "preprocess"
                result["next_action"] = (
                    f"Call preprocess_dataset(gse='{gse}') to detect format and generate config"
                )
            elif dry_run:
                result["next_step"] = "download"
                result["next_action"] = (
                    f"Call download_dataset(gse='{gse}') without dry_run to actually download"
                )
            else:
                result["next_step"] = "retry"
                result["next_action"] = (
                    "Some files failed. Check the 'files' list for errors and retry."
                )

            return json.dumps(result, indent=2, ensure_ascii=False)

        except ImportError as e:
            return json.dumps(
                {
                    "error": f"Cannot import geo_downloader: {e}",
                    "action": "Ensure wget or curl is installed, and run from the project root",
                },
                indent=2,
            )
        except Exception as e:
            return json.dumps(
                {
                    "error": str(e),
                    "gse": gse,
                    "action": "Check that FUXI_DATA_ROOT is set and writable",
                },
                indent=2,
            )

    # ── preprocess_dataset ──────────────────────────────────────────────

    @server.tool()
    async def preprocess_dataset(
        gse: str,
        modality: str = "",
        data_root: str = "",
        dry_run: bool = False,
        force: bool = False,
        download: bool = False,
    ) -> str:
        """Preprocess a downloaded GEO dataset: detect format and generate config YAML.

        Runs Phase 0-6 of the preprocessing pipeline: validates input,
        extracts archives, detects format, infers modality, generates
        dataset.yaml and config_GSE_ID.yaml.

        If data is not yet downloaded, pass download=true to auto-download first.

        Args:
            gse: GEO accession ID (e.g., "GSE12345")
            modality: Force modality ("rna"/"atac"/"spatial"/"multiome"), auto-detect if empty
            data_root: Override FUXI_DATA_ROOT
            dry_run: Detect and report only, don't write files
            force: Overwrite existing config files
            download: Auto-download from GEO before preprocessing

        Returns:
            JSON with preprocessing results, generated config path, and next step.
        """
        gse = gse.upper()

        cmd = [
            sys.executable,
            os.path.join(_repo_root(), "core", "preprocess", "preprocessor.py"),
            "--gse",
            gse,
        ]

        if modality:
            cmd.extend(["--modality", modality])
        if data_root:
            cmd.extend(["--data-root", data_root])
        if dry_run:
            cmd.append("--dry-run")
        if force:
            cmd.append("--force")
        if download:
            cmd.append("--download")

        result = _run_subprocess(cmd, timeout=300, description=f"Preprocess {gse}")

        # Try to find the generated config path
        config_path = ""
        if result["ok"]:
            for mod in ["rna", "atac", "spatial", "bulk"]:
                candidate = os.path.join(_repo_root(), "projects", mod, gse, f"config_{gse}.yaml")
                if os.path.exists(candidate):
                    config_path = f"projects/{mod}/{gse}/config_{gse}.yaml"
                    break

        output = {
            "ok": result["ok"],
            "gse": gse,
            "elapsed_s": result["elapsed_s"],
            "dry_run": dry_run,
        }

        if result["ok"]:
            if config_path:
                output["config_path"] = config_path
                output["next_step"] = "run_pipeline"
                output["next_action"] = (
                    f"Call list_steps() to see available steps, then "
                    f"run_step(modality='{modality or 'rna'}', step=0, config_path='{config_path}')"
                )
            else:
                output["next_step"] = "review"
                output["next_action"] = (
                    "Preprocessing completed but no config found. Review stdout for details."
                )
        else:
            output["error"] = result["stderr"] or result["stdout"] or "Unknown error"
            output["next_step"] = "troubleshoot"
            output["next_action"] = "Preprocessing failed. Review the error message above."

        # Include truncated stdout for agent context
        if result["stdout"]:
            # Extract key lines: successes, failures, generated paths
            key_lines = [
                line
                for line in result["stdout"].split("\n")
                if any(
                    kw in line.lower()
                    for kw in [
                        "generated",
                        "created",
                        "wrote",
                        "config",
                        "modality",
                        "error",
                        "warning",
                        "format",
                    ]
                )
            ]
            if key_lines:
                output["key_output"] = "\n".join(key_lines[-20:])

        return json.dumps(output, indent=2, ensure_ascii=False)

    # ── run_step ────────────────────────────────────────────────────────

    @server.tool()
    async def run_step(
        modality: str,
        step: int,
        config_path: str,
        cell_type: str = "",
        timeout: int = 1200,
    ) -> str:
        """Run a single pipeline step and return structured results.

        Use list_steps() to see available steps, then call this for each step
        in sequence. Long-running steps may take minutes. The result includes
        output files, metrics, and the next step number.

        Args:
            modality: One of "rna", "atac", "spatial", "bulk"
            step: Step number (0-indexed). Use list_steps() to see available steps.
            config_path: Path to YAML config (e.g., "projects/rna/GSE12345/config_GSE12345.yaml")
            cell_type: (RNA only) Cell type name for subclustering (step 6)
            timeout: Max seconds before timing out (default 1200, i.e. 20 min)

        Returns:
            JSON with step status, elapsed time, output paths, and next step.
        """
        cmd = [
            sys.executable,
            os.path.join(_repo_root(), "core", "run_pipeline.py"),
            "--modality",
            modality,
            "--step",
            str(step),
            "--config",
            os.path.join(_repo_root(), config_path),
        ]

        if cell_type and modality == "rna" and step == 6:
            cmd.extend(["--cell-type", cell_type])

        result = _run_subprocess(cmd, timeout=timeout, description=f"Run {modality} step {step}")

        output: dict[str, Any] = {
            "modality": modality,
            "step": step,
            "ok": result["ok"],
            "elapsed_s": result["elapsed_s"],
            "config_path": config_path,
        }

        if result["timed_out"]:
            output["status"] = "timeout"
            output["action"] = (
                f"Step {step} timed out. Try with a longer timeout, "
                f"or use resume_pipeline(modality='{modality}', config_path='{config_path}')"
            )
        elif result["ok"]:
            output["status"] = "completed"
            # Try to determine next step
            try:
                from core.ai.mcp_tools.pipeline import _STEPS

                steps = _STEPS.get(modality, [])
                max_step = len(steps) - 1
                output["next_step"] = step + 1 if step < max_step else None
                if output["next_step"] is not None:
                    output["next_action"] = (
                        f"Call run_step(modality='{modality}', step={output['next_step']}, "
                        f"config_path='{config_path}')"
                    )
                else:
                    output["next_action"] = (
                        "All steps completed. Check pipeline_status() for final results."
                    )
            except Exception:
                output["next_step"] = None

            # Extract key output lines
            key_lines = [
                line
                for line in result["stdout"].split("\n")
                if any(
                    kw in line.lower()
                    for kw in [
                        "completed",
                        "saved",
                        "wrote",
                        "output",
                        "n_cells",
                        "n_genes",
                        "metrics",
                    ]
                )
            ]
            if key_lines:
                output["key_output"] = "\n".join(key_lines[-10:])
        else:
            output["status"] = "failed"
            output["error"] = result["stderr"] or result["stdout"] or "Unknown error"
            output["action"] = (
                f"Step {step} failed. Review the error, fix the config, "
                f"and retry with run_step(modality='{modality}', step={step}, config_path='{config_path}')"
            )

        return json.dumps(output, indent=2, ensure_ascii=False)

    # ── run_pipeline ────────────────────────────────────────────────────

    @server.tool()
    async def run_pipeline(
        modality: str,
        config_path: str,
        resume: bool = True,
        timeout: int = 3600,
    ) -> str:
        """Run the full pipeline (or resume from checkpoint).

        Use --resume by default to pick up from the first incomplete step.
        For a fresh start, set resume=false. This may run for minutes to hours.

        Args:
            modality: One of "rna", "atac", "spatial", "bulk"
            config_path: Path to YAML config
            resume: If true (default), resume from checkpoint
            timeout: Max seconds (default 3600, i.e. 1 hour)

        Returns:
            JSON with completion status and summary.
        """
        cmd = [
            sys.executable,
            os.path.join(_repo_root(), "core", "run_pipeline.py"),
            "--modality",
            modality,
            "--config",
            os.path.join(_repo_root(), config_path),
        ]

        if resume:
            cmd.append("--resume")

        result = _run_subprocess(cmd, timeout=timeout, description=f"Run {modality} pipeline")

        output: dict[str, Any] = {
            "modality": modality,
            "ok": result["ok"],
            "elapsed_s": result["elapsed_s"],
            "config_path": config_path,
            "resumed": resume,
        }

        if result["timed_out"]:
            output["status"] = "partial"
            output["action"] = (
                "Pipeline timed out. Call pipeline_status() to check progress, "
                "then run_pipeline() again to resume from checkpoint."
            )
        elif result["ok"]:
            output["status"] = "completed"
            output["action"] = "Pipeline completed. Call pipeline_status() to review results."
        else:
            output["status"] = "failed"
            output["error"] = result["stderr"] or result["stdout"] or "Unknown error"
            output["action"] = "Pipeline failed. Review errors and fix before retrying."

        # Extract summary lines
        summary_lines = [
            line
            for line in result["stdout"].split("\n")
            if any(
                kw in line.lower()
                for kw in ["completed", "step", "all steps", "pipeline", "resume"]
            )
        ]
        if summary_lines:
            output["summary"] = "\n".join(summary_lines[-15:])

        return json.dumps(output, indent=2, ensure_ascii=False)

    # ── paper_insights ──────────────────────────────────────────────────

    @server.tool()
    async def paper_insights(pmid: str, methodology: bool = True) -> str:
        """Extract AI-powered structured insights from a paper by PMID.

        Downloads the paper, runs AI extraction to generate insights.yaml
        with metadata, figures, methods, and (optionally) methodology patterns.

        Use this to register a new paper in the knowledge base before
        analyzing its datasets.

        Args:
            pmid: PubMed ID (e.g., "31493975")
            methodology: If true (default), also extract 5-dim methodology patterns

        Returns:
            JSON with extraction status, output paths, and key findings.
        """
        cmd = [
            sys.executable,
            os.path.join(_repo_root(), "core", "paper", "insights.py"),
            "--pmid",
            pmid,
        ]

        if methodology:
            cmd.append("--methodology")

        result = _run_subprocess(cmd, timeout=600, description=f"Paper insights for PMID {pmid}")

        output: dict[str, Any] = {
            "pmid": pmid,
            "ok": result["ok"],
            "elapsed_s": result["elapsed_s"],
            "methodology_extracted": methodology,
        }

        if result["ok"]:
            # Look for the output insights.yaml
            insights_dir = os.path.join(_repo_root(), "projects", "papers", pmid)
            insights_yaml = os.path.join(insights_dir, "insights.yaml")
            if os.path.exists(insights_yaml):
                output["insights_path"] = f"projects/papers/{pmid}/insights.yaml"
                output["next_action"] = (
                    f"Call registry_status(pmid='{pmid}') to check linked datasets, "
                    "then register via: python -m core.paper.registry register --pmid {pmid}"
                )
            else:
                output["next_action"] = (
                    "Insights extracted but YAML path not confirmed. Check stdout."
                )

            # Extract key findings
            key_lines = [
                line
                for line in result["stdout"].split("\n")
                if any(
                    kw in line.lower()
                    for kw in [
                        "wrote",
                        "extracted",
                        "key_finding",
                        "archetype",
                        "toolbox",
                        "methodology",
                    ]
                )
            ]
            if key_lines:
                output["key_output"] = "\n".join(key_lines[-10:])
        else:
            output["error"] = result["stderr"] or result["stdout"] or "Unknown error"
            output["action"] = (
                f"Paper insights failed. Check that PMID {pmid} is valid and accessible. "
                "Try running manually: python core/paper/insights.py --pmid {pmid}"
            )

        return json.dumps(output, indent=2, ensure_ascii=False)
