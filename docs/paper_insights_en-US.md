# Fuxi Paper Insights Guide

> For: **Researchers** | Extract structured insights from papers automatically, linked with single-cell analysis pipeline

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Quick Start](#3-quick-start)
4. [Four Source Types](#4-four-source-types)
5. [Output Structure](#5-output-structure)
6. [Pipeline Integration](#6-pipeline-integration)
7. [FAQ](#7-faq)

---

## 1. Overview

`paper_insights.py` uses an LLM to extract structured insights from academic papers, producing an `insights.yaml` file. No more manually reading papers—AI handles abstract extraction, key findings, experimental design, figure metadata, and annotates whether each figure can be reproduced from GEO data.

### Four-Tier Source Strategy

| Tier | Input | Quality | Dependencies |
|------|-------|---------|-------------|
| 🥇 PMC XML | `--pmid` / `--xml` | Best | Zero pip deps (stdlib) |
| 🥈 PubMed | `--pmid` | Metadata best | Zero deps (stdlib) |
| 🥉 PDF | `--pdf` | Good | `pymupdf4llm` (optional) |
| 🏅 Markdown | `.md` file | Fair | None |

### Coverage

Most papers are available as full-text JATS XML in PMC. When not in PMC, PubMed metadata (title, journal, year, authors) and abstract are fetched automatically via NCBI E-utilities — reliable registration without any extra dependency. PDF provides full-text analysis.

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

### Optional (PDF full-text)

```bash
pip install -r requirements/paper.txt
```

Not installing `pymupdf4llm` does not affect core functionality; PMC XML and PubMed fallback paths cover the vast majority of use cases.

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
# Try PMC → not in PMC? PubMed → PDF → .md
python core/paper_insights.py --pmid 31269016 --pdf paper.pdf --source auto
```

### 3.5 .md files

```bash
python core/paper_insights.py projects/papers/2019_Menon_NatCommun_Human-Retina-AMD-Atlas.md --force
```

Pass the `.md` file directly.

---

## 4. Four Source Types

### 4.1 PmcXmlSource — PMC XML (recommended)

```bash
python core/paper_insights.py --pmid 31653841    # PubMed ID
python core/paper_insights.py --doi 10.1038/s41467-019-12780-8  # DOI
python core/paper_insights.py --xml local.xml     # Local XML file
```

**Pros**: Structurally precise (`<sec>` sections, `<fig>` labels), no text gluing, no formatting noise.  
**Limits**: Requires full-text JATS XML in PMC; for papers not in PMC, PubMed metadata is used as automatic fallback.

### 4.2 PubmedSource — PubMed metadata (PMC fallback)

```bash
python core/paper_insights.py --pmid 32467236    # Works even if paper is not in PMC
```

Automatically enabled when a paper is not archived in PMC. Fetches authoritative bibliographic metadata (title, journal, year, first author, PMID, DOI) and abstract text via NCBI E-utilities.

**Pros**: 100% reliable metadata (NCBI official), zero extra deps, zero install.
**Limits**: Metadata + abstract only — no body text or figures; pair with `--pdf` for full analysis.
### 4.3 Pymupdf4llmSource — PDF fallback

```bash
pip install -r requirements/paper.txt
python core/paper_insights.py --pdf paper.pdf --force
```

**Pros**: Far better quality than markitdown (54/100 vs 14/100), single pip install.  
**Limits**: Extra dependency; prefer PMC XML when available.

### 4.4 MarkdownSource — .md files


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
The top-level `reproduction_status` field combines tracking and computed aggregate fields:
- `pipeline_run`: tracks whether this paper has been run through the downstream QC pipeline.
- `overall_match`: an optional field for recording whether the reproduction matches the paper (null by default).
- `total_figures`: computed count of figure entries in the `figures` array.
- `reproducible_count`: computed count of figures where `reproducible: true`.
- `verified_figures`: list of figures that have been manually verified.
- `notes`: free-text notes about reproduction.

Each figure entry contains a boolean `reproducible` field and a `reproducibility_reasoning` field explaining the decision. The per-figure `reproducible` value comes from LLM classification; only the aggregate counts are computed post-hoc.

---

## 6. Pipeline Integration

### 6.1 PaperRegistry — Paper Registration & Indexing

`registry.py` manages paper→GEO dataset mappings across five domains. Two registration paths:

#### Auto-register via PMID (recommended)

Downloads full text from NCBI PMC, extracts GSE IDs via AI, registers automatically:

```bash
python -m core.registry register --pmid 31493975
python -m core.registry register --xml local.xml          # Local XML
python -m core.registry register --pdf paper.pdf          # PDF fallback
python -m core.registry register --paper-dir projects/papers/.../  # Existing insights.yaml
```

#### Register via GSE accession

Fetches NCBI GEO SOFT metadata to discover PMID(s) and links to existing papers:

```bash
python -m core.registry register --gse GSE164044            # Register
python -m core.registry register --gse GSE164044 --dry-run  # Preview
```

#### Query & Verify

```bash
python -m core.registry report      # Summary report
python -m core.registry verify      # Consistency check
python -m core.registry find-orphans  # Find datasets with no linked paper
python -m core.registry reset-gse GSE12345  # Reset dataset status
```

### 6.2 run_reproduce — Reproduce papers through the pipeline

`run_reproduce.py` is the automated bridge from paper to pipeline: it detects which GSEs belong to a paper, runs preprocessing for unconfigured datasets, then launches the pipeline for all:

```bash
# Preview reproducibility for all papers
python core/run_reproduce.py --all --dry-run

# Reproduce a single paper
python core/run_reproduce.py projects/papers/2019_Menon_Nature_Com_.../

# Target a specific GSE only
python core/run_reproduce.py projects/papers/.../ --gse GSE107618
```

### 6.3 Complete workflow: PMID → Reproduction

```bash
# Step 1: Interpret paper
python core/paper_insights.py --pmid 31269016

# Step 2: Build registry
python -m core.registry report

# Step 3: Preview reproducibility
python core/run_reproduce.py --all --dry-run

# Step 4: Reproduce (requires GEO data downloaded to projects/{modality}/{GSE_ID}/)
python core/run_reproduce.py projects/papers/<paper_dir>/
```

Existing manually-configured datasets are never overwritten (`force=False` by default).

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

The LLM performs structured extraction (abstracts, gene names, figure types), not analytical judgment. Key findings are based on the original text; gene lists are extracted verbatim from paper text and figures. Check the `reproducible` field in `insights.yaml` as a quality reference.

### Q6: First run times out?

The first `--pmid` call downloads XML from NCBI (~100KB), typically completing in 2–5 seconds. If it times out, check your network or use a local `--xml` file.
