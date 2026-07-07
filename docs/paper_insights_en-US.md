# Fuxi Paper Insights Guide

> For: **Researchers** | Extract structured insights from papers automatically, linked with single-cell analysis pipeline

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Quick Start](#3-quick-start)
4. [Three Source Types](#4-three-source-types)
5. [Output Structure](#5-output-structure)
6. [Pipeline Integration](#6-pipeline-integration)
7. [FAQ](#7-faq)

---

## 1. Overview

`paper_insights.py` uses an LLM to extract structured insights from academic papers, producing an `insights.yaml` file. No more manually reading papers—AI handles abstract extraction, key findings, experimental design, figure metadata, and annotates whether each figure can be reproduced from GEO data.

### Three-Tier Source Strategy

| Tier | Input | Quality | Dependencies |
|------|-------|---------|-------------|
| 🥇 PMC XML | `--pmid` / `--xml` | Best | Zero pip deps (stdlib) |
| 🥈 PDF | `--pdf` | Good | `pymupdf4llm` (optional) |
| 🥉 Markdown | `.md` file | Fair | None |

### Coverage

Most retina research papers have full-text JATS XML in PMC, retrievable via `--pmid` with zero extra installation.

---

## 2. Installation

### Required

```bash
# LLM API key (for paper interpretation)
export LLM_API_KEY=sk-...
export LLM_MODEL=deepseek-v4-flash    # optional, this is the default

# Base dependencies
uv sync  # or pip install -r requirements/base.txt
```

### Optional (PDF fallback)

```bash
pip install -r requirements/paper.txt
```

Not installing `pymupdf4llm` does not affect core functionality; the PMC XML path (covering 88% of use cases) does not need it.

---

## 3. Quick Start

### 3.1 Basic: via PMID

```bash
python core/paper_insights.py --pmid 31269016 --force
```

This automatically:
1. Queries NCBI for the PMCID
2. Downloads full-text JATS XML
3. Extracts sections (abstract, introduction, results, discussion, methods)
4. Parses figure captions, types, and gene names
5. Calls LLM to generate `insights.yaml`

The first run caches the PMC XML under `projects/papers/{paper_name}/`; subsequent runs skip re-download.

### 3.2 Via Local XML

```bash
python core/paper_insights.py --xml tests/fixtures/pmc6814749.xml --force
```

Fully offline, zero network requests.

### 3.3 Via PDF

```bash
python core/paper_insights.py --pdf paper.pdf --force
```

### 3.4 Auto Fallback (default behavior)

```bash
# Try PMC → fall back to PDF → fall back to .md
python core/paper_insights.py --pmid 31269016 --pdf paper.pdf --source auto
```

### 3.5 .md files

```bash
python core/paper_insights.py projects/papers/2019_Menon_NatCommun_Human-Retina-AMD-Atlas.md --force
```

Pass the `.md` file directly.

---

## 4. Three Source Types

### 4.1 PmcXmlSource — PMC XML (recommended)

```bash
python core/paper_insights.py --pmid 31653841    # PubMed ID
python core/paper_insights.py --doi 10.1038/s41467-019-12780-8  # DOI
python core/paper_insights.py --xml local.xml     # Local XML file
```

**Pros**: Structurally precise (`<sec>` sections, `<fig>` labels), no text gluing, no formatting noise.  
**Limits**: Requires the paper to have full-text in PMC; ~12% of papers unavailable.

### 4.2 Pymupdf4llmSource — PDF fallback

```bash
pip install -r requirements/paper.txt
python core/paper_insights.py --pdf paper.pdf --force
```

**Pros**: Far better quality than markitdown (54/100 vs 14/100), single pip install.  
**Limits**: Extra dependency; prefer PMC XML when available.

### 4.3 MarkdownSource — .md files

```bash
python core/paper_insights.py paper.md --force
```

Supports `.md` files directly.

---

## 5. Output Structure

### 5.1 Output Files

```
projects/papers/{paper_name}/
├── PMC1234567.xml          # PMC XML cache (--pmid mode only)
├── paper.md                # Converted markdown (PDF/md modes only)
└── insights.yaml           # AI-structured insights
```

### 5.2 insights.yaml Schema

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
  - PDGFRA is expressed in retinal astrocytes
  - 58 cell types identified
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
    reproduction_feasibility: feasible   # reproducible from GEO data

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

### 5.3 reproduction_feasibility Meaning

| Status | Condition |
|--------|-----------|
| `feasible` | Figure has identifiable gene names + scRNA-seq/scATAC-seq data + GEO accession |
| `not_feasible` | No identifiable genes (pure statistics/schematics) or no public data dependency |

---

## 6. Pipeline Integration

### 6.1 Paper Genes → Pipeline Annotation

Marker genes identified from papers can automatically guide cell type annotation in single-cell data:

```
--pmid → insights.yaml → KB-aware annotation (pipeline step 05/06)
```

### 6.2 Complete Workflow Example

```bash
# Step 1: Interpret paper
python core/paper_insights.py --pmid 31269016

# Step 2: Download matching GEO data and run pipeline
python core/run_pipeline.py --modality rna --config projects/rna/GSE137537/config_GSE137537.py

# Step 3: Compare paper findings to data results
# insights.yaml genes/cell types → pipeline output marker_genes.csv / cell_types.json
```

---

## 7. FAQ

### Q1: --pmid fails with "PMC full-text not available"

```bash
python core/paper_insights.py --pdf paper.pdf --force  # fall back to PDF
python core/paper_insights.py paper.md --force          # or use existing .md
```

### Q2: PDF conversion fails with "pymupdf4llm not installed"

```bash
pip install -r requirements/paper.txt
```

Or use PMC XML (`--pmid` / `--xml`) — zero extra dependencies.

### Q3: How do I check which papers are cached?

```bash
ls projects/papers/*/
```

### Q4: How do I force re-interpretation of a paper?

```bash
python core/paper_insights.py --pmid 31269016 --force
```

### Q5: How trustworthy is AI interpretation?

The LLM performs structured extraction (abstracts, gene names, figure types), not analytical judgment. Key findings are based on the original text; gene lists are extracted verbatim from paper text and figures. Check the `reproduction_feasibility` field in `insights.yaml` as a quality reference.

### Q6: First run times out?

The first `--pmid` call downloads XML from NCBI (~100KB), typically completing in 2–5 seconds. If it times out, check your network or use a local `--xml` file.
