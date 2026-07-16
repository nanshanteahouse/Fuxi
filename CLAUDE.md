# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fuxi (伏羲) is a unified single-cell multi-omics pipeline monorepo, formed by merging previously separate `scRNAseq_pipeline` (Scanpy-based) and `ATACseq_pipeline` (Snapatac2-based) codebases. Python 3.14+, running on WSL2.

### Supported Modalities

| Modality | Engine | Steps | Status |
|----------|--------|:-----:|:------:|
| `rna` (scRNA-seq) | Scanpy 1.10+ | 13 (00–12) | Production |
| `atac` (scATAC-seq) | Snapatac2 2.9 | 10 (00–09) | Production |
| `spatial` | Squidpy 1.8+ | 11 (00–10) | Production |

### Architecture

```
core/              Shared infrastructure (no biology libs)
  config.py          Pydantic Config model + 20 topic sub-models
  run_pipeline.py    CLI orchestrator — subprocess dispatch
  ai_caller.py       LLM client (OpenAI SDK) — retry, caching
  ai_prompts.py      RNA + ATAC annotation prompt templates
  utils/             resolve_config, _perf (monitoring), etc.
  preprocess/        Format detection, archive extract, config gen
  anatomy.py         CCI anatomy constraint system
  enrichment_tissue.py / grn_tissue.py  Tissue-aware scoring
  paper_insights.py / registry.py / run_reproduce.py

rna/               scRNA-seq module
  steps/            13 pipeline steps (00_load → 12_cell_interaction)
  utils/            marker_scoring, evidence_fusion, cluster_evaluation,
                    cell_interaction, sex_detection
  annotation_standardizer.py / ortholog.py / tissue_ontologies/

atac/              scATAC-seq module (Snapatac2 2.9+)
  steps/            10 pipeline steps (00_load → 09_integrate)

spatial/           Spatial module (Squidpy)
  steps/            11 pipeline steps (00_load → 10_cell_interaction)

cross_paper/       Cross-paper pathway comparison
projects/          Dataset configs, projects/{modality}/{GSE_ID}/
templates/         Config templates for different input formats
tests/             Test files (no CI framework)
```

## Development conventions

**Commit message.** Use [Conventional Commits](https://www.conventionalcommits.org/):
```
<type>(<scope>): <subject>
```
Types: `feat` / `fix` / `docs` / `refactor` / `perf` / `test` / `chore`
Scope: semantic name (e.g. `enrichment`), NOT step number (e.g. `09_enrichment`)
Subject: imperative, lowercase, ≤72 chars. Body explains *why*, not *what*.

**Config access.** Always use nested topic paths:
`CFG.hvg.n_top_genes`, `CFG.clustering.cluster_selection_method`, `CFG.qc.min_genes`.
`.py` configs are rejected — use `.yaml`. See `core/config.py` for all 20 topic models.

**Core scripts.** Step scripts under `rna/steps/`, `atac/steps/`, `spatial/steps/` must not be edited in place. Copy to `projects/{modality}/{GSE_ID}/` for dataset-specific changes, then run the copy directly. Write a note to `notes/suggestions/` after.

**Code organization.** 500 LOC soft cap for core modules and step scripts; 400 for utility modules. Algorithm engines under `*_utils/` at 500 LOC.

## Running methods

### Environment

```bash
cd <repo_root>
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt          # all modalities
pip install -r requirements/rna.txt      # scRNA-seq only
pip install -r requirements/atac.txt     # scATAC-seq only

export FUXI_DATA_ROOT=<path_to_geo_datasets>
```

### Pipeline

```bash
# List steps
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list
python core/run_pipeline.py --modality spatial --list

# Run full / single step / range / resume
python core/run_pipeline.py --modality rna --config projects/rna/<GSE_ID>/config_<GSE_ID>.yaml
python core/run_pipeline.py --modality atac --step 0 --config ...
python core/run_pipeline.py --modality rna --steps 0-2 --config ...
python core/run_pipeline.py --modality rna --resume --config ...

# Subclustering (RNA Step 06)
python core/run_pipeline.py --modality rna --step 6 --cell-type "Müller Glia" --config ...

# Subset pipeline (filter cells by sample/obs criteria)
python core/run_pipeline.py --modality rna --config config_pcw8.yaml
# Config: CFG.downsample.sample_keep=["SCR205"] or CFG.downsample.obs_filter="stage=='PCW8'"
```

### Running modes

Both modes use `--step N` to run one step at a time. The difference is whether the Agent pauses for user input:

**Auto mode** — Execute steps sequentially with default settings. After each step, report progress and continue. User can interrupt at any point (Ctrl+C). Suitable for familiar datasets or batch reproduction.

**Interactive mode** — Execute `--step N` one at a time. After each step, present results, ask questions, offer options, and wait for confirmation before proceeding. Suitable for exploratory analysis or new datasets.

### Paper tools

```bash
python core/paper_insights.py --pmid <PMID>       # AI paper interpretation
python -m core.registry report               # print summary
python -m core.registry verify              # check registry consistency
python core/run_reproduce.py --all --dry-run       # preview reproducibility
python core/run_reproduce.py <paper_dir>           # reproduce a single paper
```

### Adding a Dataset

**Automated (recommended):**
```bash
python core/preprocess/preprocessor.py --gse <GSE_ID> --data-root $FUXI_DATA_ROOT
```
**Manual:** Create `projects/{modality}/{GSE_ID}/config_<GSE_ID>.yaml` from `templates/config_templates/`.

### Testing

No CI. Run test files directly:
```bash
python -m pytest tests/test_config_parity.py -v
```

### Key Design Patterns

**Step dispatch.** `run_pipeline.py` runs each step as a separate `subprocess.run()` — never import steps directly. Steps self-identify checkpoint files; steps skip if output checkpoint exists. `--resume` scans for first missing checkpoint.

**Config loading.** `resolve_config()` loads `.yaml` via `yaml.safe_load()` + `Config.model_validate()`. Each call returns a new Config instance. Path resolution via `model_post_init` hook.

**Three annotation modes (RNA Step 05):**
1. **Unified KB** (if `CFG.tissue_kb` set): marker scoring → evidence fusion → optional AI fallback
2. **AI** (if `CFG.ai.enabled` + `CFG.ai.ai_annotation`): LLM with StandardOntology normalization
3. **Score_genes** (fallback): `sc.tl.score_genes()` with `CFG.marker.marker_dict`

**Cluster selection (Step 04).** Controlled by `CFG.clustering.cluster_selection_method`:
- `"multi_metric"` (RNA default, MMACS v2): 5-metric composite. Tissue path (enrichment loop) or Subtype path (DE-gated for homogenous data)
- `"pareto_elbow"` (ATAC/Spatial default): Pareto frontier over n_clusters vs silhouette
- `"silhouette"` / `None` (manual)
- UMAP sweep auto-selects `min_dist`/`spread` via convex hull area.

**snRNA-seq adaptation.** Auto-detected from GEO keywords → `CFG.qc.is_nuclei=True`. Tightens mito threshold to `max_pct_mito_nuclei` (default 3.0%) and MAD multiplier to 1.5×.

**Subset filtering.** Step 00 supports `CFG.downsample.sample_keep` and `CFG.downsample.obs_filter`. Output dirs auto-append `_subset`.

**Import path hack.** No `pyproject.toml`, so every step must prepend repo root to `sys.path`:
```python
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
```

### Pipeline Steps

| Modality | Steps |
|----------|-------|
| RNA | 00_load → 01_doublet → 02_qc → 03_integrate → 04_cluster_umap → 05_annotate_major → 06_subcluster → 07_markers_de → 08_trajectory → 09_enrichment → 10_exploratory → 11_grn → 12_cell_interaction |
| ATAC | 00_load → 01_qc → 02_process → 03_cluster → 04_annotate → 05_marker_peaks → 06_motif → 07_trajectory → 08_enrichment → 09_integrate |
| Spatial | 00_load → 01_qc → 02_image → 03_normalize → 04_cluster → 05_annotate → 06_spatial_de → 07_trajectory → 08_enrichment → 09_exploratory → 10_cell_interaction |

### 更新文档

如果用户明确要求"更新文档"，先总结本次会话已完成的工作，然后通读 `CLAUDE.md`、`README.md`、`docs/`（含 tutorial），评估哪些内容需要更新。按照原有规范格式，删除/修改过时信息，补充新内容。**给出改动摘要让用户确认后再写入。**
