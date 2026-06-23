# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fuxi (伏羲) is a unified single-cell multi-omics pipeline monorepo, formed by merging previously separate `scRNAseq_pipeline` (Scanpy-based) and `ATACseq_pipeline` (Snapatac2-based) codebases. Python 3.14+, running on WSL2.

## Build / Run / Test

### Environment

```bash
# WSL: create and activate virtual environment
cd <repo_root>
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Required env var before any pipeline run
export FUXI_DATA_ROOT=<path_to_geo_datasets>
```

### Running the Pipeline

```bash
# List available steps for a modality
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list

# Run full pipeline
python core/run_pipeline.py --modality rna --config projects/rna/GSE246169/config_GSE246169.py

# Run a single step (0-indexed)
python core/run_pipeline.py --modality atac --step 0 --config projects/atac/GSE246169/config_GSE246169.py

# Run a range or selection of steps
python core/run_pipeline.py --modality rna --steps 0-2 --config ...
python core/run_pipeline.py --modality rna --steps 1,3,5 --config ...

# Resume from first incomplete checkpoint
python core/run_pipeline.py --modality rna --resume --config ...

# Run subclustering on a specific cell type (RNA Step 07)
python core/run_pipeline.py --modality rna --step 7 --cell-type "Müller Glia" --config ...
```

### Testing

There is no test runner or linting configuration in this repo. The `tests/` directory exists but contains only `__init__.py` files. There is no `pyproject.toml`, `setup.py`, `Makefile`, or CI configuration.

### Adding a Dataset

1. Create `projects/{modality}/{GSE_ID}/` with a `config_GSE_ID.py` that imports from `core.config` and mutates the `CFG` singleton
2. Place raw data in `$FUXI_DATA_ROOT/{GSE_ID}/`
3. Use config templates from `templates/config_templates/` as starting points

## Architecture

### Module Layout

```
core/               Shared infrastructure (no biology libs imported)
  config.py           Unified Config + AIConfig dataclasses; CFG singleton
  run_pipeline.py     CLI orchestrator — dispatches steps via subprocess
  ai_caller.py        LLM client (OpenAI SDK) — retry, caching, model discovery
  ai_prompts.py       RNA + ATAC annotation prompt templates + build_annotation_prompt()
  utils.py            safe_write, safe_plot, setup_logger, resolve_config, validate_adata, data_root()
  dataset_schema.py   Python model for dataset.yaml files
  dataset_detector.py Auto-detect modality from file patterns

rna/                scRNA-seq module (Scanpy 1.10+)
  steps/             12 pipeline steps (00_load → 11_exploratory), each a standalone script
  utils/
    marker_scoring.py  Hypergeometric + cosine scoring of clusters against Knowledge Base
    evidence_fusion.py 5-tier decision engine merging marker scores, expert rules, AI
  annotation_standardizer.py  6-tier name standardization for cell type annotations
  ortholog.py        Cross-species Ensembl→human gene symbol mapping
  tissue_ontologies/ Expert Knowledge Bases — currently retina only, with markers + synonyms

atac/               scATAC-seq module (Snapatac2 2.9+)
  steps/             10 pipeline steps (00_load → 09_integrate)

spatial/            Placeholder for future spatial transcriptomics (Squidpy)

projects/           Dataset-specific configs, organized as projects/{modality}/{GSE_ID}/
templates/          Config templates for different input formats
tests/              Test directory (mostly empty)
```

### Key Design Patterns

**Step dispatch model.** `run_pipeline.py` does NOT import step modules. It runs each step as a separate `subprocess.run()` call, passing `--config=<path>`. Steps self-identify their checkpoint files. This means steps are loosely coupled and can be run independently. However, each step script must add the repo root to `sys.path` manually:
```python
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
```

**Config loading.** `resolve_config()` in `core/utils.py` uses `importlib.util` to dynamically load a config `.py` file as a module, then reads its `CFG` attribute. Config files mutate the global `CFG` singleton imported from `core.config`. The `Config` dataclass contains ALL fields for BOTH modalities in one flat namespace; the `modality` string discriminates.

**Checkpoint system.** Each step reads from a specific checkpoint file and optionally writes one. `run_pipeline.py` maintains step registries (`RNA_STEPS`, `ATAC_STEPS`) and checkpoint file mappings. Steps skip if their output checkpoint already exists. The `--resume` flag scans for the first missing checkpoint.

**Three annotation modes for RNA Step 05:**
1. **Unified KB mode** (if `CFG.tissue_kb` is set): Full pipeline — marker scoring → expert rules → evidence fusion → optional AI fallback for low-confidence clusters
2. **AI mode** (if `CFG.ai.enabled` + `CFG.ai.ai_annotation`): LLM-based annotation with `StandardOntology` normalization
3. **Score_genes mode** (fallback): Simple `sc.tl.score_genes()` with `CFG.marker_dict`

**AI caller.** All LLM calls go through `core/ai_caller.py`'s `ai_query()`. It uses the OpenAI SDK, auto-detects available models on 404, retries empty responses up to 3×, caches responses to disk (SHA-256 keyed), and automatically boosts `max_tokens` to 32768 when thinking mode is enabled.

**Knowledge Base format.** Tissue KBs are Python dicts with per-cell-type marker genes (organized as `confirm`/`add`/`refine` tiers, each gene mapped to PMID references), negative markers, species lists, synonyms, and optional `expert_rules` (priority-ordered deterministic matching rules).

**Path resolution.** `data_root()` reads `FUXI_DATA_ROOT` env var (with `SCRNA_DATA_ROOT` as legacy fallback). WSL paths are auto-detected. `Config.resolve_paths()` resolves all relative paths against `project_dir` and creates output directories.

### RNA Pipeline Steps

| Step | Script | Key Output |
|------|--------|-----------|
| 00 | `00_load.py` | 00_raw.h5ad |
| 01 | `downsample.py` | (optional, overwrites 00_raw) |
| 02 | `01_doublet.py` | 01_doublet.h5ad |
| 03 | `02_qc.py` | 02_qc.h5ad |
| 04 | `03_integrate.py` | 03_integrated.h5ad |
| 05 | `04_cluster_umap.py` | 04_clustered.h5ad |
| 06 | `05_annotate_major.py` | 05_annotated.h5ad |
| 07 | `06_subcluster.py` | (requires --cell-type) |
| 08 | `07_markers_de.py` | marker CSVs |
| 09 | `08_trajectory.py` | trajectory h5ad |
| 10 | `09_enrichment.py` | enrichment CSVs |
| 11 | `06_exploratory.py` | summary figures + CSVs |

### ATAC Pipeline Steps

| Step | Script | Key Output |
|------|--------|-----------|
| 00 | `00_load.py` | 00_raw.h5ad |
| 01 | `01_qc.py` | 01_filtered.h5ad |
| 02 | `02_process.py` | 02_processed.h5ad |
| 03 | `03_cluster.py` | 03_clustered.h5ad |
| 04 | `04_annotate.py` | 04_annotated.h5ad |
| 05 | `05_marker_peaks.py` | marker_peaks.csv |
| 06 | `06_motif.py` | motif_results.csv |
| 07 | `07_trajectory.py` | 07_trajectory.h5ad |
| 08 | `08_enrichment.py` | enrichment CSVs |
| 09 | `09_integrate.py` | 09_integrated.h5ad (RNA+ATAC) |

### Import Path Hack

Because this repo has no `pyproject.toml` or `setup.py`, step scripts cannot do `from fuxi.core import ...`. Instead every step prepends the repo root to `sys.path`:
```python
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
```
This means `from core.utils import ...` works from any step script. The `annotation_standardizer.py` also patches `sys.path` to avoid name collisions with a local `utils/` package.

### Modifying Scripts for a Dataset

Core step scripts (`rna/steps/*.py`, `atac/steps/*.py`) **must not be edited in place**. When a dataset exposes a bug or requires a one-off adaptation:

1. Copy the script into the project directory: `projects/{modality}/{GSE_ID}/`
2. Modify the copy — the original under `rna/steps/` or `atac/steps/` stays untouched
3. Run the copy directly instead of through `run_pipeline.py`
4. After the run completes, write a note to `notes/suggestions/{GSE_ID}.md` describing:
   - What broke and why
   - What the workaround was
   - Whether the root cause should be fixed in the core script

This keeps core scripts reference-stable and builds a searchable record of edge cases that inform future pipeline improvements. See the existing `notes/suggestions/` directory for examples.
