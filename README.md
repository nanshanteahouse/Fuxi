# Fuxi (伏羲) — Unified Single-Cell Multi-Omics Pipeline

> **Fuxi (伏羲)**: painting the eight trigrams, bringing order from chaos — just as this pipeline transforms raw single-cell data into structured biological insights.

## Overview

Fuxi is a unified monorepo for single-cell multi-omics analysis — scRNA-seq (Scanpy), scATAC-seq (Snapatac2), Spatial (Squidpy). Python 3.14+ on WSL2.

### Supported Modalities

| Modality | Engine | Steps | Status |
|----------|--------|:-----:|:------:|
| `rna` | Scanpy 1.10+ | 13 (00-12) | Production |
| `atac` | Snapatac2 2.9 | 14 (00-13) | Production |
| `spatial` | Squidpy 1.8+ | 11 (00-10) | Production |
| `bulk` | PyDESeq2 | 5 (+1) | Beta |

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
| CSV/TSV count matrix | `count_matrix` | Bulk | `config_bulk.yaml` |
| CSV/TSV TPM matrix | `tpm_matrix` | Bulk | `config_bulk.yaml` |

**R / Seurat formats (.rds, .qs)** — not natively supported. Use the companion tool [r2h5ad](https://github.com/nanshanteahouse/r2h5ad) to convert RDS/QS files to h5ad before loading with `data_format = "h5ad"`.

### Architecture

```
fuxi/
├── core/                  # Shared infrastructure (zero modality dependencies)
│   ├── ai/                # LLM caller + prompt templates (YAML)
│   ├── annotation/        # Cell-type annotation engine + standardizer + marker scoring
│   ├── cluster/           # Grid-search clustering + parameter evaluation
│   ├── config/            # Unified Pydantic Config + dataset schema
│   ├── interaction/       # Cell-cell interaction (CCI) utilities
│   ├── kb/                # Tissue knowledge base (markers, adjacency, pathways)
│   ├── paper/             # Paper insights, registry, converter, cross-paper analysis
│   ├── pipeline/          # Pipeline runner, anatomy, enrichment, GRN, reproducibility
│   ├── preprocess/        # Format detection → archive extraction → config generation
│   ├── utils/             # I/O, logging, path resolution, validation, performance
│   ├── run_pipeline.py    # CLI entry point (thin wrapper → pipeline/runner.py)
│   ├── downsample.py      # Cell downsampling (anndata-agnostic)
│   ├── geo_downloader.py  # GEO dataset downloader
│   ├── kb_validator.py    # Marker validation against tissue KB
│   ├── cross_dataset_meta.py  # Cross-dataset meta-analysis
│   └── dataset_detector.py    # Auto-detect modality from file patterns
│   ├── tui/                # Textual-based Terminal UI (v2.x)
│   │   ├── backends/        # Async wrappers for registry/download/pipeline/config
│   │   ├── screens/         # 6 screens (home, registry, pipeline, results, data mgmt, config editor)
│   │   └── widgets/         # 4 reusable widgets (config selector, step selector, progress, log)

│
├── rna/                   # scRNA-seq module (13 steps)
│   ├── steps/             # 00_load → 12_cell_interaction
│   ├── utils/             # RNA-specific utilities (hierarchy, evidence_fusion, etc.)
│   └── ortholog.py        # Cross-species gene mapping
│
├── atac/                  # scATAC-seq module (14 steps)
│   └── steps/             # 00_load → 13_integrate
│
├── spatial/               # Spatial transcriptomics module (11 steps)
│   └── steps/             # 00_load → 10_cell_interaction
|
├── bulk/                  # Bulk RNA-seq module (5+1 steps)
│   └── steps/             # 00_load → 04_exploratory + optional 05_batch
│
├── adhoc/                 # One-off / dataset-specific scripts (use once and discard)
│   ├── migration_scripts/ # Config/kb/methodology migration tools
│   └── ortholog_scripts/  # Ortholog data processing scripts
│
├── projects/              # Project-specific data (gitignored)
│   ├── papers/            # Paper insights + NCBI XML
│   ├── notebook/          # Agent-driven brainstorming notes
│   ├── rna/               # RNA dataset configs
│   ├── atac/              # ATAC dataset configs
│   └── spatial/           # Spatial dataset configs
├── tests/                 # Unified test suite
├── templates/             # Config templates + schemas
└── docs/                  # Pipeline & architecture docs
```

## Quick Start

### Prerequisites

- Python 3.14+ (WSL2 recommended for ATAC-seq)
- Virtual environment:
  ```bash
  cd <repo_root>
  python -m venv .venv
  source .venv/bin/activate
  # 方式 A: requirements 薄包装 (向后兼容)
  pip install -r requirements.txt          # all modalities
  pip install -r requirements/rna.txt        # scRNA-seq only
  pip install -r requirements/atac.txt       # scATAC-seq only
  pip install -r requirements/spatial.txt    # spatial transcriptomics only
  # 方式 B: pyproject extras (推荐, 依赖单一声明源)
  pip install -e ".[all]"              # 等价: 全部模态
  pip install -e ".[rna]"              # 单模态
  pip install -e ".[rna,celltypist]"   # 模态 + 方法学
  ```

可选 extras: `rna` / `atac` / `spatial` / `bulk` / `paper` / `methods` /
`scvi` / `celltypist` / `scvelo` / `sccoda` / `cellbender` / `soupx` / `dev` / `all`.
依赖声明源见 `pyproject.toml` `[project.optional-dependencies]`.
`[all]` 不含方法学 (避免拉入 PyTorch) —— 显式加 `[methods]` 或单个 extra。

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

# ── Multi-sample 10X MTX (each sample in its own subdir) ──
# data_input:
#   mtx_dir: ""               # parent dir containing the sample subdirs
#   mtx_dir_pattern: "GSM*/"    # glob matching sample subdirs; Step 00 merges them
#   mtx_sample_regex: "GSM\d+_(\w+)"   # optional: extract sample name from subdir
#                                  # (empty → use the subdir basename)
#   # Identical gene sets → fast sparse vstack; differing → outer-join concat
#   # mtx_concat_batch: 0        # >0: batched (tree) concat for huge sample sets

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

## Key Modules

| Module | Purpose |
|--------|---------|
| `core/ai/` | LLM API (caller.py) + prompt templates (prompts.py + templates/) |
| `core/annotation/` | Cell-type annotation engine, name standardizer, marker scoring |
| `core/cluster/` | Grid-search clustering + multi-metric parameter evaluation |
| `core/config/` | Unified Pydantic v2 Config (schema.py) + dataset.yaml model |
| `core/interaction/` | Cell-cell interaction via LIANA+ (RNA permutation + Spatial bivariate) |
| `core/kb/` | Tissue knowledge base: markers, adjacency, pathway relevance |
| `core/paper/` | Paper insights, registry, PMC converter, cross-paper analysis |
| `core/pipeline/` | CLI runner, anatomy loading, enrichment, GRN, reproducibility |
| `core/preprocess/` | Format detection, archive extraction, config generation |
| `core/utils/` | I/O, logging, path resolution, config resolution, validation, perf |
| `core/downsample.py` | Cell downsampling (random, stratified, per-sample capping) |
| `core/kb_validator.py` | Empirical marker validation against tissue knowledge base |
| `core/geo_downloader.py` | GEO dataset downloader from NCBI |
| `rna/utils/` | RNA-specific: hierarchy builder, evidence fusion, pseudobulk DE, sex detection |
| `adhoc/` | One-off migration, ortholog processing, dataset-specific analysis |

## Performance

Fuxi is engineered for million-cell datasets on a single workstation (98 GB RAM + RTX 3090).
Measured on retina datasets at 1.05M and 1.98M cells:

| Optimization | Effect |
|--------------|--------|
| Double-write bug fix (Step 03) | write time -50%, total -35% |
| gzip level 1 (configurable) | write time -12%, total -47% vs original |
| zstd compression (hdf5plugin) | write -32%, smaller files on dense data |
| GPU hybrid path (PCA guard + Harmony on X_pca only) | Harmony 560s → 52s (-91%) |
| `stream_raw` (chunked CSR write from 02_qc) | peak RAM -24% (92 → 71 GiB @ 1.98M) |
| Batch diagnosis slim-down | 45 min → 30 s (90x), subsample max 50k cells |
| scVI training defaults (batch 1024, early-stopping, AMP) | up to 40x wall-clock |
| Unified `execution.memory` config + preflight guard | warn/block on estimated peak before running |

Memory policy: `speed` (default) / `balanced` / `memory`; see `core/utils/_memory.py`.


## Citation

If you use Fuxi in your research, please cite:

> Lun, M. **Fuxi: Unified Single-Cell Multi-Omics Pipeline** (2026).
> GitHub: https://github.com/nanshanteahouse/Fuxi
