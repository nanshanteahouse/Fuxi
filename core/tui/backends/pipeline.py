"""Async subprocess manager for running pipeline steps.

Provides:
- run_step: Run a single step as an async subprocess, yielding output lines.
- get_checkpoint_status: Detect which steps have completed checkpoints.
- get_step_dependency: Return the checkpoint file a step depends on.
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
from core.utils import _set_blas_env

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


async def run_step(
    step_index: int,
    modality: str,
    config_path: str,
    extra_args: list[str] | None = None,
    *,
    cell_type: str | None = None,
    annotate_method: str | None = None,
) -> AsyncIterator[dict]:
    """Run a single pipeline step as an async subprocess.

    Parameters
    ----------
    step_index:
        0-based index into the modality's step list.
    modality:
        One of ``"rna"``, ``"atac"``, ``"spatial"``.
    config_path:
        Absolute or relative path to the YAML config file.
    extra_args:
        Additional CLI arguments forwarded to the step script.
    cell_type:
        (RNA only) Cell type to subcluster — passed as ``--cell-type``
        to step 06_subcluster.py.
    annotate_method:
        (RNA only) Annotation method — passed as ``--annotate-method``
        to step 05_annotate_major.py.

    Yields
    ------
    dict
        With keys ``type`` (``"stdout"`` | ``"stderr"`` | ``"exit"``)
        and ``data`` (a decoded text line, or the exit code for
        ``"exit"`` events).
    """
    repo_root = _resolve_repo_root()
    mod = MODALITY_MAP[modality]
    _num, script, _desc = mod["steps"][step_index]

    script_path = os.path.join(repo_root, mod["dir"], "steps", script)

    # ── Build argument list ────────────────────────────────────────
    args_list = [f"--config={config_path}"]

    # RNA step-specific flags
    if modality == "rna":
        if step_index == 6 and cell_type:
            args_list.extend(["--cell-type", cell_type])
        if step_index == 5 and annotate_method:
            args_list.extend(["--annotate-method", annotate_method])

    # User-supplied extra arguments (e.g. --resume, --steps)
    if extra_args:
        args_list.extend(extra_args)

    # ── Environment: enforce BLAS thread limits ────────────────────
    _set_blas_env(4)
    proc_env = dict(os.environ)

    # ── Launch subprocess ──────────────────────────────────────────
    proc = await asyncio.create_subprocess_exec(
        sys.executable,
        script_path,
        *args_list,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        env=proc_env,
    )

    # ── Stream stdout / stderr concurrently ───────────────────────
    # Interleave lines from both streams so the TUI can show progress
    # in near-real-time regardless of which pipe produces output.
    stdout_gen = _tagged_lines("stdout", proc.stdout)
    stderr_gen = _tagged_lines("stderr", proc.stderr)

    async for item in _interleave(stdout_gen, stderr_gen):
        yield item

    exit_code = await proc.wait()
    yield {"type": "exit", "data": exit_code}


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
