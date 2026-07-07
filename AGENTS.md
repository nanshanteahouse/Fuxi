# Fuxi (伏羲) — Agent Quick Reference

> Full knowledge base: [CLAUDE.md](CLAUDE.md)

## One-liners

```bash
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list
python core/run_pipeline.py --modality rna --config projects/rna/<GSE_ID>/config_<GSE_ID>.py
python core/paper_insights.py --pmid <PMID>       # AI paper interpretation
python core/paper_registry.py --build              # build paper→GSE→config index
python core/paper_registry.py --verify             # check registry consistency
python core/run_reproduce.py --all --dry-run       # preview reproducibility for all papers
python core/run_reproduce.py <paper_dir>           # reproduce a single paper's pipeline
```

> Paper interpretation guide: [docs/paper_insights_zh-CN.md](docs/paper_insights_zh-CN.md)


## Key paths

| Module | Location |
|--------|----------|
| Shared core | `core/` (config, utils, ai_caller, ai_prompts, run_pipeline, preprocess) |
| RNA steps | `rna/steps/` (12 scripts) |
| ATAC steps | `atac/steps/` (10 scripts) |
| Paper insights | `core/paper_insights.py`, `core/paper_converter.py` |
| Paper registry | `core/paper_registry.py`, `core/paper_registry_models.py` |
| Reproduce mode | `core/run_reproduce.py` |
| Paper insights docs | `docs/paper_insights_zh-CN.md` |
| Project configs | `projects/{modality}/{GSE_ID}/config_*.py` |
| Paper index | `projects/papers/paper_index.html`, `projects/papers/registry.yaml` |
## Critical conventions

- Steps run as **subprocesses** via `run_pipeline.py` — never imported directly
- Every step script must add repo root to `sys.path`: `sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))`
- Config loaded dynamically: `CFG = resolve_config(args.config)`
- `data_root()` requires `FUXI_DATA_ROOT` env var (no hardcoded defaults)
- Import pattern: `from core.utils import ...`, `from core.ai_caller import ...`

## Commit message discipline

使用 [Conventional Commits](https://www.conventionalcommits.org/) 规范：

```
<type>(<scope>): <subject>

<body>

<footer>
```

常用 `<type>`：`feat` / `fix` / `docs` / `refactor` / `perf` / `test` / `chore` / `style` / `ci`
`<scope>` 选填，用语义名（如 `enrichment`、`config`），**不要** 用流水线序号（如 `09_enrichment`）
subject 用祈使语气、首字母小写、不超过 72 字符
body 回答「为什么」而非「做了什么」

历史遗留的 `Cx`/`Mx`/`Nx`/`mx`/`Sx` 前缀不再用于新 commit；
若已有声明未实现，在该 commit 的 `git notes` 中标记 `UNFIXED`

## Code organization

**拆分/合并原则**：按代码逻辑与调度边界决策，不按硬性行数指标强行拆。内容高度内聚、单一职责的文件，即使较长也优于强行拆散后相互依赖的碎片。一组总是一起 import 的微型文件可考虑合并。

**按文件类型的行数参考上限**（软性指引，不是硬规则）：

| 文件类型 | 参考上限 | 说明 |
|---------|---------|------|
| 核心模块 (`core/*.py`) | 500 LOC | 跨模态共享逻辑，允许较大但应保持内聚 |
| 步骤脚本 (`steps/*.py`) | 500 LOC | Pipeline 顺序逻辑，天然长，不鼓励拆分 |
| 工具函数集 (`*_utils/*.py`) | 400 LOC | 纯函数集合，按主题分组；超出时拆分独立模块 |
| 测试文件 | 不限 | 覆盖率优先，行数不是质量指标 |
| KB 知识库 (`sources/*.py`) | 不限 | 数据驱动，无逻辑复杂度 |
| Config 文件 | 不限 | 字段声明集中管理更易维护 |

> **算法模块注意**：如果文件路径虽含 `utils` 但内容是完整的算法引擎（如 `rna/utils/marker_scoring.py`、`rna/utils/evidence_fusion.py`），应适用核心模块上限（500 LOC）而非工具函数上限。分类看**职责**，不只看路径。

**数值依据**：基于 194 个 Python 文件（~38,921 LOC）的实际规模分布（中位数 161，P75=284，P90=431）。各类型上限落在对应 P90 附近或略上方，保留合理缓冲，避免对正常代码产生持续噪声。
