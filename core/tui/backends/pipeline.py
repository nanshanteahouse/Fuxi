"""Async subprocess manager for running pipeline steps.

Provides:
- run_step: Run a single step as an async subprocess, yielding output lines.
- get_checkpoint_status: Detect which steps have completed checkpoints.
- get_step_dependency: Return the checkpoint file a step depends on.

All step execution delegates to ``core/run_pipeline.py`` so the TUI, CLI and
MCP share a single orchestration engine (``core.pipeline.runner``).
"""

from __future__ import annotations

import asyncio
import os
import sys
from typing import AsyncIterator

from core.pipeline.runner import (
    MODALITY_MAP,
    _get_step_dependency,
    find_first_incomplete,
)

# ═══════════════════════════════════════════════════════════════════
#  Repo root resolution
# ═══════════════════════════════════════════════════════════════════


def _resolve_repo_root() -> str:
    """Return absolute path to the repository root.

    This module lives at ``core/tui/backends/pipeline.py``, so the repo
    root is 4 ``dirname`` calls up:

        core/tui/backends/pipeline.py  →  core/tui/backends/
        core/tui/backends/              →  core/tui/
        core/tui/                       →  core/
        core/                           →  <repo-root>/
    """
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )


# ═══════════════════════════════════════════════════════════════════
#  Streaming subprocess runner
# ═══════════════════════════════════════════════════════════════════


def build_run_command(
    step_index: int,
    modality: str,
    config_path: str,
    extra_args: list[str] | None = None,
    *,
    cell_type: str | None = None,
) -> list[str]:
    """Build the ``core/run_pipeline.py`` CLI command for one step.

    Delegating to the CLI entry point makes TUI execution a thin wrapper
    over the same orchestration engine as CLI/MCP: config validation,
    BLAS thread limits, modality auto-discovery, perf_report recording
    and exit-code semantics all come from ``core/pipeline/runner.py``.
    """
    repo_root = _resolve_repo_root()
    cmd = [
        sys.executable,
        os.path.join(repo_root, "core", "run_pipeline.py"),
        "--modality",
        modality,
        "--step",
        str(step_index),
        "--config",
        os.path.abspath(config_path),
    ]
    if modality == "rna" and step_index == 6 and cell_type:
        cmd.extend(["--cell-type", cell_type])
    if extra_args:
        cmd.extend(extra_args)
    return cmd


async def run_step(
    step_index: int,
    modality: str,
    config_path: str,
    extra_args: list[str] | None = None,
    *,
    cell_type: str | None = None,
) -> AsyncIterator[dict]:
    """Run a single pipeline step as an async subprocess.

    Thin wrapper over ``core/run_pipeline.py --step N`` (see
    :func:`build_run_command`) so TUI behavior matches CLI/MCP.
    Output streams interleave and are yielded as ``{"type":
    "stdout"|"stderr", "data": line}`` dicts; the final item is
    ``{"type": "exit", "data": exit_code}`` (0 = completed, 2 =
    skipped, other = failed). If the consuming task is cancelled, the
    child process is terminated and no exit event is yielded.

    Parameters
    ----------
    step_index:
        0-based index into the modality's step list.
    modality:
        One of ``"rna"``, ``"atac"``, ``"spatial"``, ``"bulk"``.
    config_path:
        Absolute or relative path to the YAML config file.
    extra_args:
        Additional CLI arguments forwarded to ``run_pipeline.py``.
    cell_type:
        (RNA only) Cell type to subcluster — passed as ``--cell-type``
        for step 6 (step 06_subcluster.py).
    """
    cmd = build_run_command(
        step_index,
        modality,
        config_path,
        extra_args=extra_args,
        cell_type=cell_type,
    )
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=_resolve_repo_root(),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    try:
        stdout_gen = _tagged_lines("stdout", proc.stdout)
        stderr_gen = _tagged_lines("stderr", proc.stderr)
        async for item in _interleave(stdout_gen, stderr_gen):
            yield item
        exit_code = await proc.wait()
        yield {"type": "exit", "data": exit_code}
    finally:
        # Never leak the child: on CancelledError / GeneratorExit the
        # process is terminated (mirrors runner interrupt semantics).
        if proc.returncode is None:
            proc.terminate()
            try:
                await asyncio.wait_for(proc.wait(), timeout=5)
            except (asyncio.TimeoutError, ProcessLookupError):
                proc.kill()


async def _tagged_lines(
    tag: str,
    stream: asyncio.StreamReader | None,
) -> AsyncIterator[dict]:
    """Read lines from *stream*, yielding ``{"type": tag, "data": line}``."""
    if stream is None:
        return
    while True:
        line = await stream.readline()
        if not line:
            break
        yield {"type": tag, "data": line.decode("utf-8", errors="replace").rstrip("\r\n")}


async def _interleave(
    *generators: AsyncIterator[dict],
) -> AsyncIterator[dict]:
    """Interleave items from multiple async generators as they arrive.

    Each *generator* is consumed line-by-line; as soon as any generator
    produces a value it is yielded immediately.  Exhausted generators
    are removed from the scheduling loop.
    """
    # Map Task → generator for bookkeeping
    pending: dict[asyncio.Task[dict], AsyncIterator[dict]] = {}

    for gen in generators:
        task = asyncio.ensure_future(gen.__anext__())
        pending[task] = gen

    while pending:
        done, _ = await asyncio.wait(
            pending.keys(),
            return_when=asyncio.FIRST_COMPLETED,
        )
        for task in done:
            gen = pending.pop(task)
            try:
                yield task.result()
            except StopAsyncIteration:
                continue  # generator exhausted, no re-schedule

            # Schedule the next line from this generator
            new_task = asyncio.ensure_future(gen.__anext__())
            pending[new_task] = gen


# ═══════════════════════════════════════════════════════════════════
#  Checkpoint introspection
# ═══════════════════════════════════════════════════════════════════


def get_checkpoint_status(modality: str, h5ad_dir: str) -> list[bool]:
    """Return a boolean list, one per step, indicating checkpoint completeness.

    Steps before the first incomplete checkpoint are marked ``True``;
    the first incomplete step and all subsequent steps are ``False``.

    Delegates to :func:`runner.find_first_incomplete`.
    """
    mod = MODALITY_MAP[modality]
    first_incomplete = find_first_incomplete(
        h5ad_dir,
        mod["steps"],
        mod["checkpoints"],
        mod["write_checkpoints"],
        sentinels=mod["sentinels"],
    )
    n = len(mod["steps"])
    return [i < first_incomplete for i in range(n)]


def get_step_dependency(step_index: int, modality: str) -> str:
    """Return the checkpoint filename that *step_index* reads from.

    Delegates to :func:`runner._get_step_dependency`.
    """
    mod = MODALITY_MAP[modality]
    return _get_step_dependency(
        step_index,
        mod["steps"],
        mod["checkpoints"],
        modality,
    )
