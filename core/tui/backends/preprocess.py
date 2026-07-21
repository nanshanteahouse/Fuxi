#!/usr/bin/env python3
"""
preprocess.py — Async wrapper for data pre-processing.

Exposes the synchronous ``core/preprocess/preprocessor.py`` CLI and
``core.dataset_detector.scan_directory()`` as async coroutines suitable for
use from the TUI event loop.
"""

from __future__ import annotations

import asyncio
import os
import sys
from collections.abc import AsyncIterator
from typing import Optional

from core.dataset_detector import scan_directory
from core.utils import repo_root

# ── Preprocessor subprocess ──────────────────────────────────────────────


def _preprocessor_path() -> str:
    """Absolute path to the pre-processor entry-point script."""
    return os.path.join(repo_root(), "core", "preprocess", "preprocessor.py")


async def preprocess_async(
    gse_id: str,
    data_root: str,
    modality: Optional[str] = None,
    skip_extract: bool = False,
    query_ncbi: bool = False,
) -> AsyncIterator[str]:
    """Run the pre-processor as a subprocess and yield output lines live.

    Always passes ``--force`` (non-interactive mode).  Optional flags are
    appended only when their corresponding parameter is truth-y or non-None.

    Yields
    ------
    str
        One line of subprocess output at a time (stripped of trailing
        newline).

    Raises
    ------
    asyncio.CancelledError
        If the enclosing task is cancelled the subprocess is killed
        immediately.
    """
    cmd = [
        sys.executable,
        _preprocessor_path(),
        "--gse",
        gse_id,
        "--data-root",
        data_root,
        "--force",
    ]

    if modality is not None:
        cmd.extend(["--modality", modality])
    if skip_extract:
        cmd.append("--no-extract")
    if query_ncbi:
        cmd.append("--query-ncbi")

    process = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
    )

    try:
        assert process.stdout is not None
        async for raw_line in process.stdout:
            line = raw_line.decode("utf-8", errors="replace").rstrip("\n\r")
            yield line

        await process.wait()
    except asyncio.CancelledError:
        process.kill()
        await process.wait()
        raise

    if process.returncode != 0:
        raise RuntimeError(f"Pre-processor exited with code {process.returncode} (GSE: {gse_id})")


# ── Format detection (thread-pool wrapper) ────────────────────────────────


async def detect_formats_async(input_dir: str) -> dict:
    """Detect file formats in *input_dir* via the synchronous detector.

    Wraps ``core.dataset_detector.scan_directory()`` in
    ``asyncio.to_thread()`` so it does not block the event loop.

    Returns
    -------
    dict
        Detection result with keys ``modalities``, ``samples``, and
        ``unmatched_files``.
    """
    return await asyncio.to_thread(scan_directory, input_dir)
