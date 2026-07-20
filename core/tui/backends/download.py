"""Async wrappers for GEO data download operations.

Uses ``asyncio.create_subprocess_exec`` to run ``core/geo_downloader.py``
as a subprocess for the main download (streaming progress lines without
blocking the event loop), and wraps blocking library calls
(``list_suppl_files``, ``fetch_soft_metadata``) via
``asyncio.to_thread``.
"""

from __future__ import annotations

import asyncio
import os
import sys
from collections.abc import AsyncIterator

# Ensure repo root is on sys.path (same pattern as core/tui/__main__.py)
_repo_root = os.path.dirname(
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from core.geo_downloader import fetch_soft_metadata, list_suppl_files


async def download_gse_async(
    gse_id: str,
    dest_dir: str,
) -> AsyncIterator[str]:
    """Run ``geo_downloader.py`` as a subprocess, yielding stdout lines.

    Args:
        gse_id: GEO accession ID (e.g. ``GSE123456``).
        dest_dir: Data root directory passed as ``--data-root`` to the
            subprocess (usually ``$FUXI_DATA_ROOT``).

    Yields:
        Each line of stdout (stripped) as the subprocess produces it.

    Raises:
        asyncio.CancelledError: If the calling task is cancelled —
            the subprocess is killed before re-raising.
        RuntimeError: If the subprocess exits with a non-zero code.
    """
    script = os.path.join(_repo_root, "core", "geo_downloader.py")

    process = await asyncio.create_subprocess_exec(
        sys.executable,
        script,
        "--gse",
        gse_id,
        "--data-root",
        dest_dir,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
    )

    assert process.stdout is not None

    try:
        while True:
            line = await process.stdout.readline()
            if not line:
                break
            yield line.decode("utf-8", errors="replace").rstrip("\n\r")
    except asyncio.CancelledError:
        process.kill()
        await process.wait()
        raise

    await process.wait()
    if process.returncode and process.returncode != 0:
        raise RuntimeError(
            f"geo_downloader.py exited with code {process.returncode}"
        )


async def list_suppl_async(gse_id: str) -> list[dict]:
    """Fetch supplementary file listing in a thread.

    Wraps :func:`core.geo_downloader.list_suppl_files` in
    :func:`asyncio.to_thread`.

    Args:
        gse_id: GEO accession ID.

    Returns:
        List of file descriptors with download URLs attached.
    """
    return await asyncio.to_thread(list_suppl_files, gse_id)


async def fetch_meta_async(gse_id: str) -> dict:
    """Fetch SOFT-format metadata in a thread.

    Wraps :func:`core.geo_downloader.fetch_soft_metadata` in
    :func:`asyncio.to_thread`.

    Args:
        gse_id: GEO accession ID.

    Returns:
        Structured metadata dict.
    """
    return await asyncio.to_thread(fetch_soft_metadata, gse_id)
