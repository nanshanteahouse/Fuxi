"""MCP tool modules — each module exposes async functions decorated as @mcp.tool()."""

from core.ai.mcp_tools.execution import register_execution_tools
from core.ai.mcp_tools.pipeline import register_pipeline_tools
from core.ai.mcp_tools.registry import register_registry_tools

__all__ = ["register_registry_tools", "register_pipeline_tools", "register_execution_tools"]
