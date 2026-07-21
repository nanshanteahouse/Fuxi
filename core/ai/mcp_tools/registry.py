"""MCP tools for the paper/dataset registry (read-only queries).

Exposes: registry_status, list_papers, find_orphans
Each function wraps existing core.paper.registry APIs — no logic duplication.
"""

from __future__ import annotations

import json
import os
from typing import Any

from mcp.server import MCPServer


def _load_registry():
    """Lazy-load the master registry. Returns MasterRegistry."""
    from core.paper.registry import load_master_registry

    return load_master_registry()


def _fmt_status(gse_id: str | None, pmid: str | None) -> str:
    """Return a structured JSON status report mimicking the CLI `status` command."""
    reg = _load_registry()

    result: dict[str, Any] = {}

    if pmid:
        paper = reg.get_paper(pmid) or reg.get_paper_by_pmid(pmid)
        if paper is None:
            return json.dumps(
                {
                    "error": f"PMID {pmid} not found in registry",
                    "action": f"Run: python core/paper/insights.py --pmid {pmid} then register",
                },
                indent=2,
            )
        result["paper"] = {
            "paper_id": paper.paper_id,
            "pmid": paper.pmid,
            "slug": paper.slug,
            "title": paper.title,
            "journal": paper.journal,
            "year": paper.year,
            "first_author": paper.first_author,
            "doi": paper.doi,
        }
        linked = reg.get_dataset_links(paper.paper_id)
        datasets_out = []
        for ds_id, role in linked:
            ds_info = _ds_status_dict(reg, ds_id, role.name)
            datasets_out.append(ds_info)
        result["datasets"] = datasets_out
        if not datasets_out:
            result["note"] = "No linked datasets for this paper"

    elif gse_id:
        gse_id = gse_id.upper()
        ds_info = _ds_status_dict(reg, gse_id, None)
        result["dataset"] = ds_info
        # Also show linked papers
        paper_links = reg.get_paper_links(gse_id)
        if paper_links:
            papers_out = []
            for p_id, role in paper_links:
                p = reg.get_paper(p_id)
                papers_out.append(
                    {
                        "pmid": p_id,
                        "slug": p.slug if p else "?",
                        "title": (p.title or "?")[:80] if p else "?",
                        "role": role.name,
                    }
                )
            result["linked_papers"] = papers_out
    else:
        return json.dumps({"error": "Must specify gse or pmid"}, indent=2)

    return json.dumps(result, indent=2, ensure_ascii=False)


def _ds_status_dict(reg, gse_id: str, role: str | None) -> dict[str, Any]:
    """Build a per-dataset status dict."""
    ds = reg.get_dataset(gse_id)
    info: dict[str, Any] = {"gse_id": gse_id}

    if ds:
        info["registered"] = True
        info["registry_status"] = str(ds.status)
    else:
        info["registered"] = False
        info["action"] = f"python -m core.paper.registry register --gse {gse_id}"
        return info

    if role:
        info["role"] = role

    # Data status
    data_root = os.environ.get("FUXI_DATA_ROOT", "")
    data_dir = os.path.join(data_root, gse_id) if data_root else ""
    info["data_downloaded"] = bool(data_dir and os.path.isdir(data_dir))
    info["data_dir"] = data_dir if info["data_downloaded"] else None

    # Config status
    configs: dict[str, list[str]] = {}
    for mod in ["rna", "atac", "spatial", "bulk"]:
        config_dir = os.path.join("projects", mod, gse_id)
        if os.path.isdir(config_dir):
            cfgs = [
                f
                for f in os.listdir(config_dir)
                if f.endswith(".yaml") and f.startswith("config_")
            ]
            if cfgs:
                configs[mod] = cfgs
    info["configs"] = configs or None

    # Registry modality info
    if ds and ds.modalities:
        mod_info = {}
        for mod_key, mod_val in ds.modalities.items():
            cfgs_out = []
            for cfg in mod_val.configs:
                cfgs_out.append(
                    {
                        "path": cfg.path,
                        "pipeline_status": cfg.pipeline_status,
                    }
                )
            mod_info[mod_key] = cfgs_out
        info["modalities"] = mod_info

    # Next step suggestion
    if not info["data_downloaded"]:
        info["next_step"] = "download"
        info["next_action"] = f"Download data from GEO: {gse_id}"
    elif not configs:
        info["next_step"] = "generate_config"
        info["next_action"] = (
            f"Call generate_config(gse='{gse_id}', modality='rna') "
            f"or create projects/rna/{gse_id}/config_{gse_id}.yaml from template"
        )
    else:
        info["next_step"] = "run_pipeline"
        first_mod = next(iter(configs.keys()), "rna")
        first_cfg = configs[first_mod][0]
        info["next_action"] = (
            f"python core/run_pipeline.py --modality {first_mod} "
            f"--config projects/{first_mod}/{gse_id}/{first_cfg} --step 0"
        )

    return info


def _search_papers(query: str) -> str:
    """Fuzzy-search registered papers by keyword in title, author, journal, year, PMID, slug."""
    reg = _load_registry()
    q = query.lower().strip()

    def matches(p) -> bool:
        if not q:
            return True
        fields = [
            p.title or "",
            p.first_author or "",
            p.journal or "",
            str(p.year or ""),
            p.pmid or "",
            p.slug or "",
        ]
        return any(q in f.lower() for f in fields)

    results = [p for p in reg.papers if matches(p)]

    # If query looks like "author:X", filter by author
    if ":" in q:
        prefix, value = q.split(":", 1)
        value = value.strip()
        if prefix == "author":
            results = [p for p in results if value in (p.first_author or "").lower()]

    output = []
    for p in results:
        linked = reg.get_dataset_links(p.paper_id)
        output.append(
            {
                "paper_id": p.paper_id,
                "pmid": p.pmid,
                "slug": p.slug,
                "title": p.title,
                "journal": p.journal,
                "year": p.year,
                "first_author": p.first_author,
                "n_datasets": len(linked),
                "datasets": [ds_id for ds_id, _ in linked],
                "doi": p.doi,
            }
        )

    return json.dumps(
        {"query": query, "n_results": len(output), "results": output},
        indent=2,
        ensure_ascii=False,
    )


def register_registry_tools(server: MCPServer) -> None:
    """Register all registry-related MCP tools on the given server."""

    @server.tool()
    async def registry_status(gse: str = "", pmid: str = "") -> str:
        """Check registration, data, config, and pipeline status for a dataset or paper.

        Use this to discover what datasets exist, whether data is downloaded,
        whether configs are generated, and what the next pipeline step should be.

        Args:
            gse: GSE accession ID (e.g., "GSE123456")
            pmid: PubMed ID (e.g., "31493975")

        Returns:
            JSON with registration status, data availability, config paths,
            pipeline progress, and an "action" / "next_step" suggestion.
        """
        return _fmt_status(gse_id=gse or None, pmid=pmid or None)

    @server.tool()
    async def list_papers(query: str = "") -> str:
        """Search registered papers by keyword, author, journal, year, PMID, or slug.

        Use 'author:Name' to filter by first author.
        Leave query empty to list all registered papers.

        Args:
            query: Search term (empty = all papers). Supports "author:X" filter.

        Returns:
            JSON with matching papers, each showing title, PMID, linked datasets.
        """
        return _search_papers(query)

    @server.tool()
    async def find_orphans() -> str:
        """List orphan datasets — datasets with data/config but no linked paper.

        Returns:
            JSON with orphan dataset IDs, their modalities, and status.
        """
        reg = _load_registry()
        orphans = reg.find_orphans()

        if not orphans:
            return json.dumps({"n_orphans": 0, "message": "No orphan datasets found"}, indent=2)

        results = []
        for ds_id, ds in orphans:
            modalities = list(ds.modalities.keys()) if ds.modalities else []
            results.append(
                {
                    "gse_id": ds_id,
                    "status": str(ds.status),
                    "modalities": modalities,
                    "species": ds.species,
                    "tissue": ds.tissue,
                }
            )

        return json.dumps(
            {"n_orphans": len(results), "orphans": results},
            indent=2,
            ensure_ascii=False,
        )
