# Fuxi (伏羲) — Unified Single-Cell Multi-Omics Pipeline

> **Fuxi (伏羲)**: painting the eight trigrams, bringing order from chaos — just as this pipeline transforms raw single-cell data into structured biological insights.

## Overview

Fuxi is a unified monorepo for single-cell multi-omics analysis, merging the previously separate `scRNAseq_pipeline` (Scanpy-based) and `ATACseq_pipeline` (Snapatac2-based) into a single codebase with a shared core infrastructure.

### Supported Modalities

| Modality | Engine | Steps | Status |
|----------|--------|:-----:|:------:|
| `rna` | Scanpy 1.10+ | 13 (00-12) | ✅ Production |
| `atac` | Snapatac2 2.9 | 10 (00-09) | ✅ Production |
| `spatial` | Squidpy 1.8+ | 11 (00-10) | ✅ Production |

### Supported Input Formats

| Format | data_format | Modality | Template |
|--------|-------------|----------|----------|
| 10X HDF5 (.h5) | `10X_h5` | RNA | `config_10X_h5.yaml` |
| 10X MTX (matrix.mtx + barcodes + features) | `10X_mtx` | RNA | `config_10X_mtx.yaml` |
| CSV / TSV count matrix | `csv_matrix` | RNA | `config_csv_matrix.yaml` |
| Pre-existing h5ad | `h5ad` | RNA | — |
| Preprocessed TSV (metadata cols + expression) | `preprocessed` | RNA | `config_preprocessed.yaml` |
| 10X Fragments (fragments.tsv.gz) | `10x_fragments` | ATAC | `config_fragments.yaml` |
| 10X Visium (SpaceRanger output) | `visium` | Spatial | `config_visium.yaml` |

**R / Seurat formats (.rds, .qs)** — not natively supported. Use the companion tool [r2h5ad](https://github.com/nanshanteahouse/r2h5ad) to convert RDS/QS files to h5ad before loading with `data_format = "h5ad"`:

### Architecture

```
fuxi/
├── core/              # Shared infrastructure
│   ├── prompts/       # LLM prompt templates (YAML)
│   ├── preprocess/    # Format detection → extraction → config generation
│   ├── ai_prompts.py  # Prompt imports + annotation build functions
│   ├── ai_caller.py   # Unified LLM API with retry + caching
│   ├── config.py      # Unified Config dataclass (Pydantic v2)
│   ├── run_pipeline.py# CLI dispatcher for all modalities
│   ├── paper_insights.py  # AI-assisted paper metadata + methodology extraction
│   ├── paper_converter.py # PMC XML / PDF → structured markdown
│   ├── registry.py    # Paper ↔ dataset unified registry
│   └── methodology_batch.py  # Parallel methodology pattern backfill
├── rna/               # scRNA-seq module (13 steps)
├── atac/              # scATAC-seq module (10 steps)
├── spatial/           # Spatial transcriptomics module (11 steps)
├── projects/          # Project-specific data (gitignored)
│   ├── papers/        # Paper insights + NCBI XML
│   ├── notebook/      # Agent-driven brainstorming notes
│   ├── rna/           # RNA dataset configs
│   ├── atac/          # ATAC dataset configs
│   └── spatial/       # Spatial dataset configs
├── tests/             # Unified test suite
├── templates/         # Config templates + schemas
└── docs/              # Pipeline & architecture docs
```

## Quick Start

### Prerequisites

- Python 3.14+ (WSL2 recommended for ATAC-seq)
- Virtual environment:
  ```bash
  cd <repo_root>
  python -m venv .venv
  source .venv/bin/activate
  pip install -r requirements.txt          # all modalities
  pip install -r requirements/rna.txt        # scRNA-seq only
  pip install -r requirements/atac.txt       # scATAC-seq only
  pip install -r requirements/spatial.txt    # spatial transcriptomics only
  ```

### Running the Pipeline

```bash
# List available steps
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list
python core/run_pipeline.py --modality spatial --list

# Run a full pipeline
python core/run_pipeline.py --modality rna --config projects/rna/<GSE_ID>/config_<GSE_ID>.yaml

# Run a single step
python core/run_pipeline.py --modality atac --step 0 --config projects/atac/<GSE_ID>/config_<GSE_ID>.yaml

# Resume from checkpoint
python core/run_pipeline.py --modality rna --resume --config projects/rna/<GSE_ID>/config_<GSE_ID>.yaml
```

### Data Organization

Raw data lives in a directory configured via the **`FUXI_DATA_ROOT`** environment variable. Each dataset directory contains a `dataset.yaml` metadata file. Pipeline project configs live in `projects/{modality}/{GSE_ID}/`.

```bash
# Required: set data root before running any pipeline
export FUXI_DATA_ROOT=/mnt/e/data              # WSL
# or
set FUXI_DATA_ROOT=E:/data                     # Windows
```

## Project Config Pattern (YAML)

```yaml
# projects/rna/<GSE_ID>/config_<GSE_ID>.yaml
modality: rna
species: human
tissue: retina

# Data input
data_input:
  mtx_prefix: "GSE12345_Sample1_"

# QC thresholds
qc:
  min_genes: 500
  max_genes: 7500
  max_pct_mito: 20.0

# HVG selection
hvg:
  n_top_genes: 4000
  flavor: seurat_v3

# Clustering
clustering:
  cluster_selection_method: multi_metric
  param_grid_n_neighbors: [15, 20, 30]
  param_grid_resolutions: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]

# Annotation
marker:
  marker_dict:
    "Rod Photoreceptor": ["RHO", "GNAT1"]
    "Müller Glia": ["RLBP1", "GLUL"]

# AI annotation
ai:
  enabled: false
```
```

## Key Modules

| Module | Purpose |
|--------|---------|
| `core/utils.py` | I/O, logging, config resolution, AnnData validation, marker loading |
| `core/ai_caller.py` | Unified LLM API with retry, thinking mode, disk caching |
| `core/ai_prompts.py` | Prompt imports from YAML + annotation build functions |
| `core/prompts/` | LLM prompt templates stored as YAML (7 prompt groups) |
| `core/config.py` | Unified Config dataclass (Pydantic v2) for all modalities |
| `core/run_pipeline.py` | CLI with `--modality rna|atac|spatial` dispatch |
| `core/paper_insights.py` | AI-driven paper metadata, figures, methods, and methodology extraction |
| `core/paper_converter.py` | PMC XML / PDF → structured markdown |
| `core/registry.py` | Paper ↔ dataset master registry (44 papers, 69 datasets) |
| `core/methodology_batch.py` | Parallel methodology pattern backfill (ThreadPoolExecutor) |
| `core/dataset_schema.py` | Python model for dataset.yaml |
| `core/preprocess/` | Format detection, archive extraction, config generation |

## Citation

If you use Fuxi in your research, please cite:

> Lun, M. **Fuxi: Unified Single-Cell Multi-Omics Pipeline** (2026).
> GitHub: https://github.com/nanshanteahouse/Fuxi
