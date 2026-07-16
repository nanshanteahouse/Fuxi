# Fuxi 论文解读使用指南

> 适用：**科研人员** | 从论文中自动提取结构化见解，与单细胞分析管线联动

---

## 目录

1. [概述](#1-概述)
2. [安装](#2-安装)
3. [快速开始](#3-快速开始)
4. [三种输入源](#4-三种输入源)
5. [解读结果说明](#5-解读结果说明)
6. [与管线联动](#6-与管线联动)
7. [常见问题（FAQ）](#7-常见问题faq)

---

## 1. 概述

`paper_insights.py` 利用 LLM 从论文中自动提取结构化见解，生成 `insights.yaml`。不再需要手动通读论文——AI 替你完成摘要、关键发现、实验设计、图元数据的提取，并标注每张图是否能基于 GEO 数据复现。

### 三层输入策略

| 层 | 输入 | 质量 | 依赖 |
|----|------|------|------|
| 🥇 PMC XML | `--pmid` / `--xml` | 最佳 | 零 pip 依赖（stdlib） |
| 🥈 PDF | `--pdf` | 良好 | `pymupdf4llm`（可选） |
| 🥉 Markdown | `.md` 文件 | 一般 | 零依赖 |

### 覆盖率

多数视网膜研究论文在 PMC 中有全文 JATS XML，可直接用 `--pmid` 获取，无需任何额外安装。

---

## 2. 安装

### 必须

```bash
# LLM API（用于论文解读）
export LLM_API_KEY=sk-...
export LLM_MODEL=deepseek-v4-flash    # 可选，默认值

# 基础依赖
uv sync  # 或 pip install -r requirements/base.txt
```

### 可选（PDF 回退）

```bash
pip install -r requirements/paper.txt
```

不安装 `pymupdf4llm` 不影响核心功能；PMC XML 路径（占 88% 用例）不需要它。

---

## 3. 快速开始

### 3.1 最基本用法：通过 PMID

```bash
python core/paper_insights.py --pmid 31269016 --force
```

这会自动：
1. 从 NCBI 查询 PMCID
2. 下载 JATS XML 全文
3. 提取章节（摘要、引言、结果、讨论、方法）
4. 解析每张图的标题、类型、基因名
5. 调用 LLM 生成 `insights.yaml`

第一次运行会缓存 PMC XML 到 `projects/papers/{论文名}/`，下次无需重新下载。

### 3.2 通过本地 XML

```bash
python core/paper_insights.py --xml tests/fixtures/pmc6814749.xml --force
```

完全离线，零网络请求。

### 3.3 通过 PDF

```bash
python core/paper_insights.py --pdf paper.pdf --force
```

### 3.4 自动回退（默认行为）

```bash
# 尝试 PMC → 失败则回退到 PDF → 再失败则尝试 .md
python core/paper_insights.py --pmid 31269016 --pdf paper.pdf --source auto
```

### 3.5 .md 文件

```bash
python core/paper_insights.py projects/papers/2019_Menon_NatCommun_Human-Retina-AMD-Atlas.md --force
```

直接传入 `.md` 文件即可。

---

## 4. 三种输入源

### 4.1 PmcXmlSource — PMC XML（推荐）

```bash
python core/paper_insights.py --pmid 31653841    # PubMed ID
python core/paper_insights.py --doi 10.1038/s41467-019-12780-8  # DOI
python core/paper_insights.py --xml local.xml     # 本地 XML 文件
```

**优点**：结构精确（`<sec>` 分节、`<fig>` 分图标签），无文本粘连，不含格式噪声。  
**局限**：需要论文在 PMC 中有全文；约 12% 论文不可用。

### 4.2 Pymupdf4llmSource — PDF 回退

```bash
pip install -r requirements/paper.txt
python core/paper_insights.py --pdf paper.pdf --force
```

**优点**：质量远超 markitdown（54/100 vs 14/100），单 pip 安装。  
**局限**：需要额外依赖；PMC XML 优先。

### 4.3 MarkdownSource — .md 文件

```bash
python core/paper_insights.py paper.md --force
```

直接支持 `.md` 文件输入。

---

## 5. 解读结果说明

### 5.1 输出文件

```
projects/papers/{论文名}/
├── PMC1234567.xml          # PMC XML 缓存（仅 --pmid 模式）
├── paper.md                # 转换后的 markdown（仅 PDF/md 模式）
└── insights.yaml           # AI 结构化见解
```

### 5.2 insights.yaml 结构

```yaml
paper_meta:
  year: "2025"
  first_author: "Zhang"
  journal: "Nature Genetics"

experimental_design:
  species: homo_sapiens
  tissue: retina
  tissue_info: "macular and peripheral retina from postmortem donors"
  models:
    - name: "postmortem human retina"
      description: "Six normal donors"
  conditions:
    - name: "Normal"
      description: "Control retinas"
  modalities:
    - snRNA-seq
  summary: "Single-nucleus RNA-seq on postmortem human retinal samples from 6 donors"

key_findings:
  - "58 transcriptionally distinct cell types identified"
  - "Novel subtypes of amacrine cells discovered"

data_access:
  geo_ids:
    - GSE137537
  sra_ids: []

methods:
  key_methods:
    - "10x Genomics Chromium Single Cell 3' v3"
    - CellRanger
  software_versions:
    CellRanger: "7.0"
  reference_genome: hg38
  sequencing_platforms:
    - "Illumina NovaSeq 6000"

figures:
  - id: Fig_1
    caption: "Single-cell transcriptomic analysis of human retina."
    type: umap
    panels:
      - 1a
      - 1b
    parameters:
      features:
        - PDE6A
      resolution: 0.8
      method: ACTIONet
      conditions:
        - Normal
      n_value: "n=6 donors"
      error_bar_type: SD
    purpose: "Study overview showing all major retinal cell types."
    reproducible: true
    reproducibility_reasoning: "UMAP generated from scRNA-seq data -- reproducible with dataset access."

data_notes:
  - "20,091 cells after QC"
  - "snRNA-seq -- use is_nuclei=True"

reproduction_status:
  pipeline_run: "not_started"
  overall_match: null
  total_figures: 12
  reproducible_count: 9
  verified_figures: []
  notes: ""
```

### reproduction_status
顶层的 `reproduction_status` 字段包含追踪字段和计算聚合字段：
- `pipeline_run`：记录该论文是否已通过下游 QC 流程。
- `overall_match`：可选字段，用于记录复现是否与论文一致（默认 null）。
- `total_figures`：`figures` 数组中的图形条目数量（自动计算）。
- `reproducible_count`：`reproducible: true` 的图形数量（自动计算）。
- `verified_figures`：已手动验证的图形列表。
- `notes`：关于复现的备注。

每个图形条目包含布尔类型 `reproducible` 字段和 `reproducibility_reasoning` 字段用于解释判断依据。逐图的 `reproducible` 值来自 LLM 分类；仅聚合计数为后期计算。

---

## 6. 与管线联动

### 6.1 PaperRegistry — 构建论文索引

`registry.py` 扫描所有论文的 `insights.yaml` 和已有项目配置，自动构建论文→GEO→配置的映射关系：

```bash
python -m core.registry report   # 生成 projects/registry/{papers,datasets,links}.yaml
python -m core.registry verify  # 检查一致性
python -m core.registry report --dry-run  # 预览不写入
```

生成的 `registry.yaml` 为每个 GSE 标记状态：`config_exists`（已有配置）、`not_configured`（待生成）、`data_not_downloaded`（需下载数据）等。

### 6.2 run_reproduce — 论文复现

`run_reproduce.py` 是论文到管线的自动化桥梁：检测每篇论文对应的 GSE，对有配置的直接跑管线，没配置的先预处理再跑：

```bash
# 预览预览全部论文
python core/run_reproduce.py --all --dry-run

# 复现单篇论文
python core/run_reproduce.py projects/papers/2019_Menon_Nature_Com_.../

# 只跑某个指定 GSE
python core/run_reproduce.py projects/papers/.../ --gse GSE107618
```

### 6.3 完整工作流：PMID → 复现

```bash
# Step 1: 解读论文
python core/paper_insights.py --pmid 31269016

# Step 2: 构建注册表
python -m core.registry report

# Step 3: 预览可复现性
python core/run_reproduce.py --all --dry-run

# Step 4: 实际复现（需先下载 GEO 数据到 projects/{modality}/{GSE_ID}/）
python core/run_reproduce.py projects/papers/<paper_dir>/
```

已有手动配置的数据集不会被覆盖（`force=False` 默认行为）。

```

---

## 7. 常见问题（FAQ）

### Q1: --pmid 报错 "PMC full-text not available"

```bash
python core/paper_insights.py --pdf paper.pdf --force  # 回退到 PDF
python core/paper_insights.py paper.md --force          # 或使用旧 .md
```

### Q2: PDF 转换报错 "pymupdf4llm not installed"

```bash
pip install -r requirements/paper.txt
```

或使用 PMC XML（`--pmid` / `--xml`）——零额外依赖。

### Q3: 如何查看缓存了哪些论文？

```bash
ls projects/papers/*/
```

### Q4: 如何强制重新解读已有论文？

```bash
python core/paper_insights.py --pmid 31269016 --force
```

### Q5: AI 解读的可信度如何？

LLM 负责结构化提取（摘要、基因名、图类型），不涉及分析判断。关键发现基于原文，gene list 从论文文本和图表中逐字提取。建议查看 `insights.yaml` 中的 `reproducible` 字段做质量参考。

### Q6: 第一次运行超时？

首次 `--pmid` 需要从 NCBI 下载 XML（~100KB），通常 2-5 秒完成。如果超时，检查网络或使用本地 `--xml` 文件。
