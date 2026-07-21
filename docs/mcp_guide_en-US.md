# Fuxi MCP Server — Setup & Usage Guide

> For: **AI Agent users** | Connect Fuxi to Claude Desktop, VS Code Copilot, Cursor, or any MCP-compatible AI tool.

---

## Table of Contents

1. [What is MCP?](#1-what-is-mcp)
2. [What Can the Fuxi MCP Server Do?](#2-what-can-the-fuxi-mcp-server-do)
3. [Prerequisites](#3-prerequisites)
4. [Starting the MCP Server](#4-starting-the-mcp-server)
5. [Claude Desktop Setup](#5-claude-desktop-setup)
6. [VS Code / Cursor Setup](#6-vs-code--cursor-setup)
7. [Custom AI Agent Setup](#7-custom-ai-agent-setup)
8. [Tool Reference](#8-tool-reference)
9. [FAQ](#9-faq)

---

## 1. What is MCP?

**MCP** (Model Context Protocol) is an open standard that lets AI applications (Claude, ChatGPT, VS Code) securely connect to external tools and data sources.

Think of it as **USB-C for AI**. Instead of copy-pasting data for the AI to read, the AI can directly invoke your analysis tools, query results, and trigger computations.

The Fuxi MCP server acts as a translator — turning Fuxi's registry, pipeline status, data download, preprocessing, and step execution into "tools" that AI can understand and call.

## 2. What Can the Fuxi MCP Server Do?

The Fuxi MCP server provides **10 tools** in three categories:

### Query Tools (read-only, no side effects)

| Tool | Function |
|------|----------|
| `registry_status` | Check registration, data, config, and pipeline status for a GSE or PMID |
| `list_papers` | Search registered papers by keyword/author/year |
| `find_orphans` | List orphan datasets (data exists but not linked to a paper) |
| `list_steps` | List all pipeline steps for a given modality |
| `pipeline_status` | Check pipeline progress — which steps are complete, what's next |

### Execution Tools (trigger real operations)

| Tool | Function |
|------|----------|
| `download_dataset` | Download a GEO dataset (supports dry-run preview) |
| `preprocess_dataset` | Auto-detect format and generate config YAML |
| `run_step` | Run a single pipeline step |
| `run_pipeline` | Run the full pipeline (supports checkpoint resume) |
| `paper_insights` | Extract AI-powered structured insights from a paper |

### Typical Workflow

```
You (to AI): "Analyze GSE123456 for me"

AI automatically executes:
  registry_status(gse="GSE123456")
  → "Registered, data downloaded, no config → need to generate config"

  preprocess_dataset(gse="GSE123456")
  → "Generated config_GSE123456.yaml"

  list_steps(modality="rna")
  → "13 steps available"

  run_pipeline(modality="rna", config="projects/rna/GSE123456/config_GSE123456.yaml")
  → "All steps complete!"
```

> **Note**: `GSE123456` in the examples above is a placeholder — replace with your own GEO accession ID.

## 3. Prerequisites

1. **Python 3.14+** (Fuxi virtual environment activated)
2. **FUXI_DATA_ROOT** environment variable set (data root directory)
3. **LLM_API_KEY** environment variable set (for AI annotation, optional)
4. **mcp package installed**: `pip install "mcp==2.0.0b2"` (inside Fuxi venv)

Verify prerequisites:
```bash
source .venv/bin/activate
echo $FUXI_DATA_ROOT        # Should print your data root path
python -c "import mcp; print(mcp.__version__)"  # Should print version
```

## 4. Starting the MCP Server

### Option 1: stdio mode (recommended, for local AI tools)

```bash
source .venv/bin/activate
python -m core.ai.mcp_server
```

The server runs in the background, communicating with the AI host via standard input/output. Press `Ctrl+C` to exit.

### Option 2: HTTP mode (for remote clients)

```bash
source .venv/bin/activate
python -m core.ai.mcp_server --http 8080
```

The server listens on `http://127.0.0.1:8080`, accepting requests from remote MCP clients.

### Option 3: Dev debugging (using MCP Inspector)

```bash
source .venv/bin/activate
uv run mcp dev core/ai/mcp_server.py
```

Open the URL printed by the Inspector to test tools interactively in your browser.

## 5. Claude Desktop Setup

### 5.1 Locate the Config File

Claude Desktop's MCP configuration is stored at:
- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux**: `~/.config/Claude/claude_desktop_config.json`

### 5.2 Add the Fuxi Server

Edit `claude_desktop_config.json` and add under `mcpServers`:

```json
{
  "mcpServers": {
    "fuxi": {
      "command": "/mnt/d/Projects/Fuxi/.venv/bin/python",
      "args": [
        "-m",
        "core.ai.mcp_server"
      ],
      "cwd": "/mnt/d/Projects/Fuxi",
      "env": {
        "FUXI_DATA_ROOT": "/mnt/e/data",
        "LLM_API_KEY": "sk-your-key-here"
      }
    }
  }
}
```

> **Note**: Replace paths with your own. `command` points to `.venv/bin/python`, `cwd` points to the project root.

### 5.3 Verify

Restart Claude Desktop. In a new conversation, type:

> Help me look up retina-related papers in the Fuxi registry

Claude will automatically call `list_papers(query="retina")` and return results.

## 6. VS Code / Cursor Setup

VS Code and Cursor support MCP servers through extensions. Recommended: **Cline** or **Continue**.

### Using Cline

Add the MCP server in Cline settings:

```json
{
  "mcpServers": {
    "fuxi": {
      "command": "/mnt/d/Projects/Fuxi/.venv/bin/python",
      "args": ["-m", "core.ai.mcp_server"],
      "cwd": "/mnt/d/Projects/Fuxi",
      "env": {
        "FUXI_DATA_ROOT": "/mnt/e/data"
      }
    }
  }
}
```

### Using Continue

In Continue config (`~/.continue/config.json`):

```json
{
  "experimental": {
    "mcpServers": [
      {
        "name": "fuxi",
        "command": "/mnt/d/Projects/Fuxi/.venv/bin/python",
        "args": ["-m", "core.ai.mcp_server"],
        "cwd": "/mnt/d/Projects/Fuxi",
        "env": {
          "FUXI_DATA_ROOT": "/mnt/e/data"
        }
      }
    ]
  }
}
```

## 7. Custom AI Agent Setup

Any AI agent that supports the MCP protocol can connect to Fuxi.

### Via stdio (subprocess mode)

```python
import asyncio
from mcp import Client
from mcp.client.stdio import stdio_client

async def main():
    async with stdio_client(
        command="/mnt/d/Projects/Fuxi/.venv/bin/python",
        args=["-m", "core.ai.mcp_server"],
        cwd="/mnt/d/Projects/Fuxi",
    ) as (read, write):
        async with Client(read, write) as client:
            # List all available tools
            tools = await client.list_tools()
            for t in tools.tools:
                print(f"  {t.name}: {t.description}")

            # Call a tool
            result = await client.call_tool("registry_status", {"gse": "GSE123456"})
            print(result.content[0].text)

asyncio.run(main())
```

### Via HTTP (remote mode)

```python
async with Client("http://127.0.0.1:8080/mcp") as client:
    result = await client.call_tool("list_steps", {"modality": "rna"})
```

### Using in-memory client (test mode)

```python
from mcp import Client
from core.ai.mcp_server import _create_server

async def test():
    server = _create_server()
    async with Client(server) as client:
        result = await client.call_tool("find_orphans", {})
        print(result.content[0].text)
```

## 8. Tool Reference

### Input/Output Format

All tools return **JSON strings**. Every output includes `next_step` and `next_action` fields, enabling the AI agent to determine the next operation automatically.

### Complete Tool List

| # | Tool Name | Parameters | Returns |
|---|-----------|------------|---------|
| 1 | `registry_status` | `gse: str`, `pmid: str` (one required) | Registration status + data availability + config path + next step suggestion |
| 2 | `list_papers` | `query: str` (empty=all, supports `author:X` filter) | Matching papers with title/journal/year/linked datasets |
| 3 | `find_orphans` | none | Orphan dataset list with modalities and status |
| 4 | `list_steps` | `modality: str` (rna/atac/spatial/bulk) | Step number, script name, description |
| 5 | `pipeline_status` | `config_path: str` or `modality + gse` | Per-step completion status, next step number |
| 6 | `download_dataset` | `gse: str`, `dry_run: bool` | Download status, file counts, next step suggestion |
| 7 | `preprocess_dataset` | `gse: str`, `modality: str` (optional) | Preprocessing result, generated config path |
| 8 | `run_step` | `modality: str`, `step: int`, `config_path: str` | Step status, elapsed time, output paths, next step number |
| 9 | `run_pipeline` | `modality: str`, `config_path: str`, `resume: bool` | Completion status, summary |
| 10 | `paper_insights` | `pmid: str`, `methodology: bool` | Extraction status, insights.yaml path, key findings |

### Agent Auto-Decision Chain

Tool outputs are designed to **drive agent decision-making**:

```json
{
  "ok": true,
  "status": "completed",
  "next_step": 1,
  "next_action": "Call run_step(modality='rna', step=1, config_path='projects/rna/GSE123456/config_GSE123456.yaml')"
}
```

The agent reads `next_action` and executes the suggested command without human intervention.

## 9. FAQ

### Q: Server fails with "FUXI_DATA_ROOT not set"

**A**: Set the environment variable or configure it in `.env`:

```bash
export FUXI_DATA_ROOT=/mnt/e/data
# or
echo "FUXI_DATA_ROOT=/mnt/e/data" >> .env
```

### Q: `download_dataset` takes a long time with no response

**A**: This is expected. GEO datasets can be tens of GB. Try `dry_run=true` first to preview.

### Q: `run_step` timed out. What now?

**A**: Call `run_pipeline(modality=..., config_path=..., resume=true)` to resume from checkpoint. Fuxi's checkpoint mechanism ensures completed work is not repeated.

### Q: How do I verify the MCP server is connected in Claude Desktop?

**A**: Click the 🔌 icon next to the input box in Claude Desktop. You should see "fuxi" server with 10 tools listed.

### Q: Can multiple AI tools connect simultaneously?

**A**: stdio mode supports only one client at a time. For multiple clients, use HTTP mode: `python -m core.ai.mcp_server --http 8080`.

### Q: How do I view server debug logs?

**A**: In Claude Desktop, logs are located at:
- **macOS**: `~/Library/Logs/Claude/mcp*.log`
- **Windows**: `%APPDATA%\Claude\logs\mcp*.log`
- Or modify the server code: `logging.basicConfig(level=logging.DEBUG)`
