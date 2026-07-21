#!/usr/bin/env python3
"""Fuxi MCP server — exposes pipeline, registry, and analysis tools to AI agents.

Usage:
    python -m core.ai.mcp_server              # stdio mode (default, for MCP hosts)
    python -m core.ai.mcp_server --http 8080  # Streamable HTTP mode
    python -m core.ai.mcp_server --sse 8080   # SSE mode (legacy)

Phase 1: registry query + pipeline status (5 read-only tools)
Phase 2: download, preprocess, run, insights (5 execution tools)
"""

from __future__ import annotations

import argparse
import os
import sys

# Ensure repo root is on sys.path
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)


def _create_server():
    """Create and configure the Fuxi MCP server with all registered tools."""
    from mcp.server import MCPServer

    from core.ai.mcp_tools import (
        register_execution_tools,
        register_pipeline_tools,
        register_registry_tools,
    )

    server = MCPServer(
        "fuxi",
        instructions=(
            "Fuxi is a unified single-cell multi-omics pipeline for scRNA-seq (Scanpy), "
            "scATAC-seq (Snapatac2), and Spatial (Squidpy) analysis. "
            "Use this server to: "
            "1) Query the paper/dataset registry to discover available data, "
            "2) Check pipeline step definitions and progress, "
            "3) Download datasets, preprocess, and run pipeline steps. "
            "Always call registry_status() or list_papers() first to discover what datasets exist."
        ),
    )

    register_registry_tools(server)
    register_pipeline_tools(server)
    register_execution_tools(server)

    return server


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fuxi MCP server — AI agent interface for single-cell pipeline",
    )
    parser.add_argument(
        "--http",
        type=int,
        nargs="?",
        const=8000,
        default=None,
        metavar="PORT",
        help="Run in Streamable HTTP mode on the given port (default: 8000)",
    )
    parser.add_argument(
        "--sse",
        type=int,
        nargs="?",
        const=8000,
        default=None,
        metavar="PORT",
        help="Run in SSE mode on the given port (default: 8000, legacy)",
    )
    args = parser.parse_args()

    server = _create_server()

    import asyncio
    import logging

    logging.basicConfig(level=logging.WARNING)

    if args.http is not None:
        import uvicorn

        print(f"[fuxi-mcp] Starting HTTP server on port {args.http}...", file=sys.stderr)
        app = server.streamable_http_app()
        uvicorn.run(app, host="127.0.0.1", port=args.http, log_level="warning")
    elif args.sse is not None:
        import uvicorn

        print(f"[fuxi-mcp] Starting SSE server on port {args.sse}...", file=sys.stderr)
        app = server.sse_app()
        uvicorn.run(app, host="127.0.0.1", port=args.sse, log_level="warning")
    else:
        # Default: stdio mode
        asyncio.run(server.run_stdio_async())


if __name__ == "__main__":
    main()
