"""Async wrappers for the paper/dataset registry API.

All functions are ``async def`` and wrap blocking I/O in ``asyncio.to_thread()``
to prevent blocking the Textual event loop.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

logger = logging.getLogger(__name__)

# ── Re-export enums for color-coded badges in the TUI ──
try:
    from core.paper.registry import (
        DatasetStatus,
        InsightStatus,
        ModalityStatus,
        _cmd_add_paper,
        _cmd_register_gse,
        load_master_registry,
        save_master_registry,
    )

    _HAS_REGISTRY = True
except ImportError:
    logger.warning("core.paper.registry not available — registry functions return empty data")
    InsightStatus = None  # type: ignore[misc]
    DatasetStatus = None  # type: ignore[misc]
    ModalityStatus = None  # type: ignore[misc]
    _HAS_REGISTRY = False

__all__ = [
    "get_all_papers_async",
    "get_all_datasets_async",
    "search_registry_async",
    "get_status_async",
    "register_paper_async",
    "register_dataset_async",
    "InsightStatus",
    "DatasetStatus",
    "ModalityStatus",
]


# ── helpers ──────────────────────────────────────────────────────────────


def _warn_no_registry() -> None:
    logger.warning("core.paper.registry is unavailable — returning empty result")


# ── public API ───────────────────────────────────────────────────────────


async def get_all_papers_async() -> list[Any]:
    """Return every registered paper entry as a list.

    Wraps ``load_master_registry()`` + ``registry.papers``.
    """
    if not _HAS_REGISTRY:
        _warn_no_registry()
        return []

    reg = await asyncio.to_thread(load_master_registry)
    return list(reg.papers)


async def get_all_datasets_async() -> dict[str, Any]:
    """Return every registered dataset keyed by GSE ID.

    Wraps ``load_master_registry()`` + ``registry.datasets``.
    """
    if not _HAS_REGISTRY:
        _warn_no_registry()
        return {}

    reg = await asyncio.to_thread(load_master_registry)
    return dict(reg.datasets)


async def search_registry_async(
    query: str,
) -> tuple[list[Any], dict[str, Any]]:
    """Search papers and datasets by keyword.

    Matches against GSE ID, PMID, title, journal, first author, and slug.
    Returns ``(matching_papers, matching_datasets)``.
    """
    if not _HAS_REGISTRY:
        _warn_no_registry()
        return [], {}

    reg = await asyncio.to_thread(load_master_registry)
    q = query.lower()

    matching_papers: list[Any] = []
    for p in reg.papers:
        fields = [
            p.paper_id or "",
            p.pmid or "",
            p.title or "",
            p.journal or "",
            p.first_author or "",
            p.slug or "",
            p.year or "",
        ]
        if any(q in f.lower() for f in fields):
            matching_papers.append(p)

    matching_datasets: dict[str, Any] = {}
    for ds_id, ds in reg.datasets.items():
        if q in ds_id.lower():
            matching_datasets[ds_id] = ds

    return matching_papers, matching_datasets


async def get_status_async(gse_id: str) -> dict[str, Any]:
    """Get comprehensive status dict for a GSE dataset.

    Constructs the result from the ``MasterRegistry`` API — **does NOT**
    call ``_cmd_status()`` (which is print-only).

    Returns a dict with keys:
    ``gse_id``, ``registered``, ``status``, ``type``, ``non_pipeline``,
    ``papers`` (list of linked paper summaries), ``modalities``, ``data_root``.
    """
    if not _HAS_REGISTRY:
        _warn_no_registry()
        return {"gse_id": gse_id.upper(), "registered": False, "status": "unknown"}

    gse_id = gse_id.upper()
    reg = await asyncio.to_thread(load_master_registry)
    ds = reg.get_dataset(gse_id)

    if ds is None:
        return {
            "gse_id": gse_id,
            "registered": False,
            "status": "not_registered",
            "papers": [],
            "modalities": [],
            "data_root": "",
        }

    # Gather linked papers via MasterRegistry API
    paper_links = reg.get_paper_links(gse_id)
    papers: list[dict[str, Any]] = []
    for p_id, role in paper_links:
        paper = reg.get_paper(p_id)
        papers.append(
            {
                "paper_id": p_id,
                "slug": paper.slug if paper else None,
                "title": paper.title if paper else None,
                "role": role.name,
            }
        )

    return {
        "gse_id": gse_id,
        "registered": True,
        "status": ds.status,
        "type": ds.type,
        "non_pipeline": ds.non_pipeline,
        "papers": papers,
        "modalities": list(ds.modalities.keys()) if ds.modalities else [],
        "data_root": ds.data_root or "",
    }


async def register_paper_async(pmid: str) -> Any:
    """Register a paper by PMID (wraps ``_cmd_add_paper``).

    Does HTTP to NCBI — runs in a thread so the event loop stays responsive.
    Returns the updated ``MasterRegistry`` or ``None`` on failure.
    """
    if not _HAS_REGISTRY:
        _warn_no_registry()
        return None

    reg = await asyncio.to_thread(load_master_registry)
    reg = await asyncio.to_thread(_cmd_add_paper, reg, pmid=pmid)
    await asyncio.to_thread(save_master_registry, reg)
    return reg


async def register_dataset_async(gse_id: str) -> Any:
    """Register a GSE dataset (wraps ``_cmd_register_gse``).

    Fetches SOFT metadata from NCBI and links to known papers.
    Returns the updated ``MasterRegistry`` or ``None`` on failure.
    """
    if not _HAS_REGISTRY:
        _warn_no_registry()
        return None

    reg = await asyncio.to_thread(load_master_registry)
    reg = await asyncio.to_thread(_cmd_register_gse, reg, gse_id)
    await asyncio.to_thread(save_master_registry, reg)
    return reg
