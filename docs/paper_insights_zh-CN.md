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
  year: 2019
  first_author: Menon
  journal: Nature Communications
  doi: 10.1038/s41467-019-12780-8

experimental_design:
  species: Human
  tissue: Retina
  technologies:
    - scRNA-seq
  conditions:
    - Normal
    - AMD

key_findings:
  - 发现 PDGFRA 在视网膜星形胶质细胞中表达
  - 鉴定出 58 种细胞类型
  - ...

data_notes:
  accessions:
    - GSE137537
  cell_count: 149045
  quality: "high"

figures:
  - figure_id: Fig_1
    caption: "Study overview and cell atlas"
    figure_type: overview
    panel_count: 5
    genes: []
    reproduction_feasibility: feasible   # 可基于 GEO 数据复现

  - figure_id: Fig_3
    caption: "Subcluster analysis"
    figure_type: umap
    panel_count: 4
    genes: [PDGFRA, GFAP]
    reproduction_feasibility: feasible

reproduction_status:
  total_figures: 37
  reproducible: 19
  not_reproducible: 18
```

### 5.3 reproduction_feasibility 含义

| 状态 | 条件 |
|------|------|
| `feasible` | 图中有可识别基因名 + scRNA-seq/scATAC-seq 数据 + GEO 编号 |
| `not_feasible` | 图中无可识别基因（纯统计/示意图）或不依赖公开数据 |

---

## 6. 与管线联动

### 6.1 论文基因 → 管线注释

论文中识别的标记基因可自动指导单细胞数据中的细胞类型注释：

```
─pmid ─→ insights.yaml ─→ KB-aware 注释（管线 step 05/06）
```

### 6.2 完整工作流示例

```bash
# Step 1: 解读论文
python core/paper_insights.py --pmid 31269016

# Step 2: 下载对应的 GEO 数据并运行管线
python core/run_pipeline.py --modality rna --config projects/rna/GSE137537/config_GSE137537.py

# Step 3: 对比论文发现与数据结果
# insights.yaml 中的基因/细胞类型 → 管线输出中的 marker_genes.csv / cell_types.json
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

LLM 负责结构化提取（摘要、基因名、图类型），不涉及分析判断。关键发现基于原文，gene list 从论文文本和图表中逐字提取。建议查看 `insights.yaml` 中的 `reproduction_feasibility` 字段做质量参考。

### Q6: 第一次运行超时？

首次 `--pmid` 需要从 NCBI 下载 XML（~100KB），通常 2-5 秒完成。如果超时，检查网络或使用本地 `--xml` 文件。
