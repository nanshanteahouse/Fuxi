# Fuxi MCP 服务器 — 配置与使用指南

> 适用人群：**AI Agent 用户** | 帮助你将 Fuxi 接入 Claude Desktop、VS Code Copilot、Cursor 等支持 MCP（Model Context Protocol）的 AI 工具。

---

## 目录

1. [MCP 是什么？](#1-mcp-是什么)
2. [Fuxi MCP 服务器能做什么？](#2-fuxi-mcp-服务器能做什么)
3. [前置条件](#3-前置条件)
4. [启动 MCP 服务器](#4-启动-mcp-服务器)
5. [配置 Claude Desktop](#5-配置-claude-desktop)
6. [配置 VS Code / Cursor](#6-配置-vs-code--cursor)
7. [配置自定义 AI Agent](#7-配置自定义-ai-agent)
8. [工具参考](#8-工具参考)
9. [常见问题](#9-常见问题)

---

## 1. MCP 是什么？

**MCP**（Model Context Protocol，模型上下文协议）是一种开放标准，让 AI 应用（如 Claude、ChatGPT、VS Code）能够安全地连接到外部工具和数据源。

简单理解：**MCP 就像 AI 的 USB-C 接口**。以前你要手动复制粘贴数据给 AI 看，现在 AI 可以直接"调用"你的分析工具，查询结果，触发计算。

Fuxi 的 MCP 服务器就是一个"翻译官"——把 Fuxi 的注册表、管线状态、数据下载、预处理、运行步骤等功能，全部变成 AI 可以理解和调用的"工具"。

## 2. Fuxi MCP 服务器能做什么？

Fuxi MCP 服务器提供 **10 个工具**，分为三类：

### 查询类（只读，无副作用）

| 工具 | 功能 |
|------|------|
| `registry_status` | 查询某个 GSE 或 PMID 的注册、数据、配置、管线状态 |
| `list_papers` | 按关键词/作者/年份搜索已注册论文 |
| `find_orphans` | 列出所有"孤儿"数据集（有数据但未关联论文）|
| `list_steps` | 列出某个模态下的所有管线步骤 |
| `pipeline_status` | 查看某个数据集的管线进度（哪步完成、下一步是什么）|

### 执行类（有副作用，触发实际操作）

| 工具 | 功能 |
|------|------|
| `download_dataset` | 从 GEO 下载数据集（支持 dry-run 预览）|
| `preprocess_dataset` | 自动检测格式并生成 config YAML |
| `run_step` | 运行单个管线步骤 |
| `run_pipeline` | 运行完整管线（支持从 checkpoint 恢复）|
| `paper_insights` | 通过 AI 提取论文的结构化洞察 |

### 典型工作流

```
你（对 AI 说）："帮我分析 GSE123456"

AI 自动执行：
  registry_status(gse="GSE123456")
  → "已注册，数据已下载，无配置 → 需要生成配置"

  preprocess_dataset(gse="GSE123456")
  → "已生成 config_GSE123456.yaml"

  list_steps(modality="rna")
  → "13 个步骤可用"

  run_pipeline(modality="rna", config="projects/rna/GSE123456/config_GSE123456.yaml")
  → "全部步骤完成！"
```

> **注意**：以上示例中的 `GSE123456` 是占位符，请替换为你自己的 GEO 数据集编号。

## 3. 前置条件

1. **Python 3.14+**（Fuxi 虚拟环境已激活）
2. **FUXI_DATA_ROOT** 环境变量已设置（数据根目录）
3. **LLM_API_KEY** 环境变量已设置（用于 AI 注释功能，可选）
4. **mcp 包已安装**：`pip install "mcp==2.0.0b2"`（位于 Fuxi 虚拟环境中）

验证前置条件：
```bash
source .venv/bin/activate
echo $FUXI_DATA_ROOT        # 应输出数据根目录路径
python -c "import mcp; print(mcp.__version__)"  # 应输出版本号
```

## 4. 启动 MCP 服务器

### 方式一：stdio 模式（推荐，用于本地 AI 工具）

```bash
source .venv/bin/activate
python -m core.ai.mcp_server
```

服务器在后台运行，通过标准输入输出与 AI 宿主通信。按 `Ctrl+C` 退出。

### 方式二：HTTP 模式（用于远程客户端）

```bash
source .venv/bin/activate
python -m core.ai.mcp_server --http 8080
```

服务器监听 `http://127.0.0.1:8080`，接受来自远程 MCP 客户端的请求。

### 方式三：开发调试（使用 MCP Inspector）

```bash
source .venv/bin/activate
uv run mcp dev core/ai/mcp_server.py
```

打开 Inspector 输出的 URL，可在浏览器中逐个测试工具。

## 5. 配置 Claude Desktop

### 5.1 找到配置文件

Claude Desktop 的 MCP 配置位于：
- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux**: `~/.config/Claude/claude_desktop_config.json`

### 5.2 添加 Fuxi 服务器

编辑 `claude_desktop_config.json`，在 `mcpServers` 中添加：

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

> **注意**：把路径换成你自己的。`command` 指向 `.venv/bin/python`，`cwd` 指向项目根目录。

### 5.3 验证

重启 Claude Desktop，在新对话中输入：

> 帮我查一下 Fuxi 注册表里有视网膜相关的论文

Claude 会自动调用 `list_papers(query="retina")` 并返回结果。

## 6. 配置 VS Code / Cursor

VS Code 和 Cursor 通过 MCP 扩展支持 MCP 服务器。推荐使用 **Cline** 或 **Continue** 扩展。

### 使用 Cline

在 Cline 设置中添加 MCP 服务器：

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

### 使用 Continue

在 Continue 配置（`~/.continue/config.json`）中添加：

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

## 7. 配置自定义 AI Agent

任何支持 MCP 协议的 AI Agent 都可以连接到 Fuxi 服务器。

### 通过 stdio 连接（子进程模式）

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
            # 列出所有可用工具
            tools = await client.list_tools()
            for t in tools.tools:
                print(f"  {t.name}: {t.description}")

            # 调用一个工具
            result = await client.call_tool("registry_status", {"gse": "GSE123456"})
            print(result.content[0].text)

asyncio.run(main())
```

### 通过 HTTP 连接（远程模式）

```python
async with Client("http://127.0.0.1:8080/mcp") as client:
    result = await client.call_tool("list_steps", {"modality": "rna"})
```

### 使用 in-memory 客户端（测试模式）

```python
from mcp import Client
from core.ai.mcp_server import _create_server

async def test():
    server = _create_server()
    async with Client(server) as client:
        result = await client.call_tool("find_orphans", {})
        print(result.content[0].text)
```

## 8. 工具参考

### 输入输出格式

所有工具返回 **JSON 字符串**。每个工具的输出都包含 `next_step` 和 `next_action` 字段，AI Agent 可以直接根据这些字段决定下一步操作。

### 完整工具列表

| # | 工具名 | 参数 | 返回内容 |
|---|--------|------|----------|
| 1 | `registry_status` | `gse: str`, `pmid: str`（二选一）| 注册状态 + 数据可用性 + 配置路径 + 下一步建议 |
| 2 | `list_papers` | `query: str`（空=全部，支持 `author:X` 筛选）| 匹配的论文列表，含标题/期刊/年份/关联数据集 |
| 3 | `find_orphans` | 无 | 孤儿数据集列表，含模态和状态 |
| 4 | `list_steps` | `modality: str`（rna/atac/spatial/bulk）| 步骤编号、脚本名、描述 |
| 5 | `pipeline_status` | `config_path: str` 或 `modality + gse` | 每个步骤的完成状态，下一步编号 |
| 6 | `download_dataset` | `gse: str`, `dry_run: bool` | 下载状态、文件数量、下一步建议 |
| 7 | `preprocess_dataset` | `gse: str`, `modality: str`（可选）| 预处理结果、生成的 config 路径 |
| 8 | `run_step` | `modality: str`, `step: int`, `config_path: str` | 步骤状态、耗时、输出路径、下一步编号 |
| 9 | `run_pipeline` | `modality: str`, `config_path: str`, `resume: bool` | 完成状态、摘要 |
| 10 | `paper_insights` | `pmid: str`, `methodology: bool` | 提取状态、insights.yaml 路径、关键发现 |

### Agent 自动决策链

工具的返回值设计为**驱动 Agent 自动决策**：

```json
{
  "ok": true,
  "status": "completed",
  "next_step": 1,
  "next_action": "Call run_step(modality='rna', step=1, config_path='projects/rna/GSE123456/config_GSE123456.yaml')"
}
```

Agent 读到 `next_action` 后可以直接执行建议的命令，无需人工介入。

## 9. 常见问题

### Q: 服务器启动时报 "FUXI_DATA_ROOT not set"

**A**: 设置环境变量或在启动前 `.env` 文件中配置：

```bash
export FUXI_DATA_ROOT=/mnt/e/data
# 或者
echo "FUXI_DATA_ROOT=/mnt/e/data" >> .env
```

### Q: `download_dataset` 运行很久不返回

**A**: 这是正常现象。GEO 数据集可能几十 GB，下载需要时间。可以先用 `dry_run=true` 预览。

### Q: `run_step` 超时了怎么办？

**A**: 调用 `run_pipeline(modality=..., config_path=..., resume=true)` 从断点恢复。Fuxi 的 checkpoint 机制确保不会重复已完成的工作。

### Q: 如何在 Claude Desktop 中确认 MCP 服务器已连接？

**A**: 在 Claude Desktop 中点击输入框旁的 🔌 图标，应能看到 "fuxi" 服务器及其 10 个工具。

### Q: 多个 AI 工具能同时连接吗？

**A**: stdio 模式一次只能有一个客户端。如果需要多个客户端同时访问，使用 HTTP 模式：`python -m core.ai.mcp_server --http 8080`。

### Q: 如何查看服务器的调试日志？

**A**: 在 Claude Desktop 中，日志位于：
- **macOS**: `~/Library/Logs/Claude/mcp*.log`
- **Windows**: `%APPDATA%\Claude\logs\mcp*.log`
- 或修改服务器代码中的 `logging.basicConfig(level=logging.DEBUG)`
