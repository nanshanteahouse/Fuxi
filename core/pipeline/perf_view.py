"""perf_view — friendly viewer for perf_report.json.

Reads ``<results_dir>/perf_report.json`` produced by the pipeline runner and
prints a clean tabular summary of every recorded step (including killed /
failed / partial runs).

Usage::

    python -m core.pipeline.perf_view <results_dir>
    python -m core.pipeline.perf_view <path/to/perf_report.json>
    python -m core.pipeline.perf_view            # auto-discover via FUXI_DATA_ROOT

Exit codes:
    0 — perf_report.json found and parsed
    1 — file missing or unreadable
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Optional


def _format_duration(sec: float) -> str:
    """Human-friendly duration: 12.3s / 5m 23s / 1h 12m."""
    if sec < 60:
        return f"{sec:6.1f}s"
    if sec < 3600:
        return f"{int(sec // 60):2d}m {int(sec % 60):02d}s"
    h = int(sec // 3600)
    m = int((sec % 3600) // 60)
    return f"{h}h {m:02d}m"


def _format_mem(mib: float) -> str:
    """MiB → human: 512 MiB / 7.1 GiB."""
    if mib < 0:
        return "    —"
    if mib < 1024:
        return f"{mib:5.0f} MiB"
    return f"{mib / 1024:5.1f} GiB"


_STATUS_LABEL = {
    "completed": "OK",
    "killed": "KILLED",
    "failed": "FAILED",
    "skipped": "SKIP",
}


def _status_cell(exit_status: str) -> str:
    """Pad status to fixed width."""
    return _STATUS_LABEL.get(exit_status, exit_status or "?").ljust(8)


def render_report(data: dict, stream=sys.stdout) -> None:
    """Print a perf_report.json dict as a friendly table."""
    info = data.get("pipeline", {}) or {}
    steps = data.get("steps", []) or []
    summary = data.get("summary", {}) or {}

    # Header
    print("=" * 78, file=stream)
    modality = info.get("modality", "?")
    print(f"Fuxi {modality.upper()} Performance Report", file=stream)
    print("=" * 78, file=stream)

    # Pipeline metadata
    print(f"  config:            {info.get('config_path', '?')}", file=stream)
    print(f"  n_jobs:            {info.get('n_jobs', '?')}", file=stream)
    first_ts = info.get("first_run_timestamp", "—")
    last_ts = info.get("last_run_timestamp", "—")
    print(f"  first_run:         {first_ts}", file=stream)
    print(f"  last_run:          {last_ts}", file=stream)
    partial = info.get("partial", False)
    partial_tag = " [PARTIAL]" if partial else ""
    total_wall = info.get("total_wall_sec") or summary.get("total_wall_sec") or 0.0
    print(f"  total_wall:        {_format_duration(total_wall)}{partial_tag}", file=stream)
    print(file=stream)

    if not steps:
        print("  (no step records)", file=stream)
        return

    # Steps table
    header = f"  {'#':>3}  {'Step':<32} {'Wall':>10}  {'Peak RSS':>11}  {'GPU':>9}  {'Cells':>10}  {'Status':<8}"
    print(header, file=stream)
    print(
        f"  {'-' * 3}  {'-' * 32} {'-' * 10}  {'-' * 11}  {'-' * 9}  {'-' * 10}  {'-' * 8}",
        file=stream,
    )

    for s in steps:
        step_id = (s.get("step", "") or "").split(" ", 1)[0]
        step_desc = (s.get("step", "") or "").split(" ", 1)[1] if " " in s.get("step", "") else ""
        wall = s.get("wall_sec", 0.0) or 0.0
        rss = s.get("peak_rss_mib", -1.0)
        gpu = s.get("gpu_mem_mb", -1.0)
        cells = s.get("n_cells", 0) or 0
        status = s.get("exit_status", "completed")
        gpu_str = f"{gpu:5.0f} MiB" if gpu >= 0 else "       —"
        cells_str = f"{cells:>10,}" if cells > 0 else f"{'—':>10}"
        print(
            f"  {step_id:>3}  {step_desc[:32]:<32} {_format_duration(wall):>10}  "
            f"{_format_mem(rss):>11}  {gpu_str:>9}  {cells_str}  {_status_cell(status)}",
            file=stream,
        )

    print(
        f"  {'-' * 3}  {'-' * 32} {'-' * 10}  {'-' * 11}  {'-' * 9}  {'-' * 10}  {'-' * 8}",
        file=stream,
    )
    n_steps = len(steps)
    n_completed = sum(1 for s in steps if s.get("exit_status") == "completed")
    print(
        f"  Total: {n_steps} step record(s), {n_completed} completed. "
        f"Max peak RSS = {_format_mem(summary.get('max_peak_rss_mib', 0.0))} "
        f"({summary.get('max_peak_rss_step', '')})",
        file=stream,
    )

    # Footnote
    if partial:
        print(
            "\n  NOTE: this is a partial report (single-step run, --resume, or "
            "interrupted execution).\n        Some pipeline steps may be missing.",
            file=stream,
        )


def _resolve_path(arg: Optional[str]) -> Optional[str]:
    """Resolve CLI argument to perf_report.json path."""
    if arg is None:
        # Auto-discover via FUXI_DATA_ROOT? Not reliable; require explicit path.
        return None
    if os.path.isfile(arg):
        return arg
    if os.path.isdir(arg):
        candidate = os.path.join(arg, "perf_report.json")
        if os.path.isfile(candidate):
            return candidate
    return None


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m core.pipeline.perf_view",
        description="Friendly viewer for Fuxi perf_report.json",
    )
    parser.add_argument(
        "path",
        nargs="?",
        help="Path to perf_report.json or to a results dir containing it.",
    )
    args = parser.parse_args(argv)

    path = _resolve_path(args.path)
    if path is None:
        if args.path:
            print(f"Error: not a file or directory: {args.path}", file=sys.stderr)
        else:
            print("Usage: python -m core.pipeline.perf_view <path>", file=sys.stderr)
        return 1

    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"Error: failed to read {path}: {e}", file=sys.stderr)
        return 1

    render_report(data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
