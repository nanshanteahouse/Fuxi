# Fuxi (伏羲) — Unified Single-Cell Multi-Omics Pipeline

> **Fuxi (伏羲)**: painting the eight trigrams, bringing order from chaos — just as this pipeline transforms raw single-cell data into structured biological insights.

## Overview

Fuxi is a unified monorepo for single-cell multi-omics analysis, merging the previously separate `scRNAseq_pipeline` (Scanpy-based) and `ATACseq_pipeline` (Snapatac2-based) into a single codebase with a shared core infrastructure.

### Supported Modalities

| Modality | Engine | Steps | Status |
|----------|--------|:-----:|:------:|
| `rna` | Scanpy 1.10+ | 12 (00-11) | ✅ Production |
| `atac` | Snapatac2 2.9 | 10 (00-09) | ✅ Production |
| `spatial` | Squidpy (planned) | TBD | 🔵 Future |

### Architecture

```
fuxi/
├── core/              # Shared infrastructure (utils, ai_caller, config, run_pipeline)
├── rna/               # scRNA-seq module (steps, utils, tissue_ontologies)
├── atac/              # scATAC-seq module (steps)
├── spatial/           # Spatial transcriptomics (future)
├── projects/          # Dataset-specific configs, organized by modality
│   ├── rna/           # 15 RNA projects
│   └── atac/          # 1 ATAC project
├── tests/             # Unified test suite
├── templates/         # Config templates
└── docs/              # Architecture & pipeline documentation
```

## Quick Start

### Prerequisites

- Python 3.14+ (WSL2 recommended for ATAC-seq)
- Virtual environment:
  ```bash
  cd <repo_root>
  python -m venv .venv
  source .venv/bin/activate
  pip install -r requirements.txt
  ```

### Running the Pipeline

```bash
# List available steps
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list

# Run a full pipeline
python core/run_pipeline.py --modality rna --config projects/rna/GSE246169/config_GSE246169.py

# Run a single step
python core/run_pipeline.py --modality atac --step 0 --config projects/atac/GSE246169/config_GSE246169.py

# Resume from checkpoint
python core/run_pipeline.py --modality rna --resume --config projects/rna/GSE107618/config_GSE107618.py
```

### Data Organization

Raw data lives in a directory configured via the **`FUXI_DATA_ROOT`** environment variable. Each dataset directory contains a `dataset.yaml` metadata file. Pipeline project configs live in `projects/{modality}/{GSE_ID}/`.

```bash
# Required: set data root before running any pipeline
export FUXI_DATA_ROOT=/mnt/e/neurobiology    # WSL
# or
set FUXI_DATA_ROOT=E:/neurobiology           # Windows
```

## Project Config Pattern

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
from core.config import Config, AIConfig
from core.utils import data_root

CFG = Config()
CFG.modality = "rna"
CFG.project_dir = os.path.dirname(os.path.abspath(__file__))
CFG.data_dir = os.path.join(data_root(), "GSE107618")
# ... dataset-specific settings ...

CFG.resolve_paths()
```

## Key Modules

| Module | Purpose |
|--------|---------|
| `core/utils.py` | safe_write, safe_plot, setup_logger, resolve_config, validate_adata, monitor_performance |
| `core/ai_caller.py` | Unified LLM calls with retry, thinking mode, disk caching, model auto-discovery |
| `core/ai_prompts.py` | RNA + ATAC annotation prompts, interpretation templates |
| `core/config.py` | Merged Config dataclass with all RNA + ATAC fields |
| `core/run_pipeline.py` | Unified CLI with `--modality rna\|atac` dispatch |
| `core/dataset_schema.py` | Python model for dataset.yaml |
| `core/dataset_detector.py` | Auto-detect modality from file patterns |
