# Fuxi Analysis Pipeline — User Guide

> For: **Single-cell omics researchers** | No programming background required

---

## Table of Contents

1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Prerequisites](#2-prerequisites)
3. [Quick start: running your first pipeline](#3-quick-start-running-your-first-pipeline)
4. [scRNA-seq pipeline in detail](#4-scrna-seq-pipeline-in-detail)
5. [scATAC-seq pipeline in detail](#5-scatac-seq-pipeline-in-detail)
6. [Spatial transcriptomics pipeline in detail](#6-spatial-transcriptomics-pipeline-in-detail)
7. [Output files reference](#7-output-files-reference)
8. [Practical tips](#8-practical-tips)
9. [Configuration file deep-dive](#9-configuration-file-deep-dive)
10. [FAQ](#10-faq)

---

## 1. What the pipeline does

After you've downloaded single-cell data from GEO and run the preprocessor to generate config files, the analysis pipeline automates the **entire computation** from raw data to biological conclusions:

| Stage | What it does | Biological significance |
|-------|-------------|------------------------|
| 🔬 Data loading | Auto-detects and reads 6 common single-cell data formats | Unifies all formats into a single internal representation |
| 🧹 Quality control | Removes doublets, dead cells, and low-quality cells | Ensures downstream analysis is based on reliable data |
| 🔗 Batch integration | Normalization + HVG + PCA + batch correction | Removes technical variation while preserving biological signal |
| 🗺️ Clustering & visualization | Multi-parameter grid search + UMAP | Discovers cell subpopulations, reveals data structure |
| 🏷️ Cell annotation | KB scoring / AI LLM / marker gene — three automatic annotation modes | Converts cluster numbers into biologically meaningful cell types |
| 🔍 Differential analysis | Marker genes + stage-wise comparison + temporal trends — three-layer DE | Identifies characteristic genes and developmental dynamics |
| 🌳 Trajectory inference | PAGA + diffusion pseudotime | Reconstructs cell differentiation/developmental paths |
| 🧬 Pathway enrichment | GO/KEGG over-representation analysis + GSEA | Reveals biological functions of cell types |
| 🧭 GRN analysis | Pseudobulk + decoupler TF activity inference | Identifies transcription factors driving each cell type identity |

**TL;DR: Config file ready → one command → from raw data to publication-quality figures, fully automated.**

---

## 2. Prerequisites

### 2.1 Install the environment

```bash
# Linux / WSL
cd /path/to/Fuxi
python -m venv .venv
source .venv/bin/activate
pip install -r requirements/rna.txt  # or requirements.txt for all modalities
```

### 2.2 Set environment variables

```bash
# Required: data root directory (top-level folder containing all downloaded raw data)
export FUXI_DATA_ROOT=/data/geo_datasets

# WSL users: must disable HDF5 file locking
export HDF5_USE_FILE_LOCKING=FALSE

# Optional: API key for AI-assisted annotation
export LLM_API_KEY=sk-your-api-key-here
```

### 2.3 Ensure config files are ready

Before running the pipeline, you need two files under your dataset directory:

```
projects/{modality}/{dataset_id}/
├── dataset.yaml          # Dataset metadata manifest
└── config_{dataset_id}.py # Pipeline configuration
```

These two files are typically generated automatically by the **preprocessor script**. If you don't have them yet, please refer to the *Fuxi Preprocessor User Guide* first.

---

## 3. Quick start: running your first pipeline

### 3.1 List available steps

```bash
# View all scRNA-seq pipeline steps
python core/run_pipeline.py --modality rna --list

# View all scATAC-seq pipeline steps
python core/run_pipeline.py --modality atac --list

# View all spatial transcriptomics pipeline steps
python core/run_pipeline.py --modality spatial --list
```

You'll see output like this:

```
Fuxi — RNA-seq pipeline step list
============================================================
  [00] Load raw data → 00_raw.h5ad
  ...
  [11] GRN regulatory network analysis (decoupler) → 11_grn.h5ad
```

### 3.2 Run the full pipeline in one command

```bash
# scRNA-seq full workflow (12 steps, from scratch to GRN)
python core/run_pipeline.py --modality rna --config projects/rna/{dataset_id}/config_{dataset_id}.py

# scATAC-seq full workflow (10 steps)
python core/run_pipeline.py --modality atac --config projects/atac/{dataset_id}/config_{dataset_id}.py

# Spatial transcriptomics full workflow (10 steps)
python core/run_pipeline.py --modality spatial --config projects/spatial/{dataset_id}/config_{dataset_id}.py
```

The terminal shows real-time progress with timing for each step:

```
============================================================
[run] [RNA] Step [00]: Load raw data → 00_raw.h5ad
============================================================
[run] Step [00] completed (took 45.2s).

============================================================
[run] [RNA] Step [02]: Scrublet doublet detection (per sample) → 01_doublet.h5ad
============================================================
[run] Step [02] completed (took 120.8s).
...
============================================================
[run] Fuxi RNA-seq pipeline execution finished.
============================================================
[run] Step timing summary:
  [00]    45.2s  Load raw data → 00_raw.h5ad
  ...
  [Total] 1845.3s  10 steps total
```

### 3.3 Checkpoints and resume

The pipeline uses a **checkpoint system**: each step saves intermediate results. If execution is interrupted for any reason, use `--resume` to continue from where it left off — completed steps are automatically skipped:

```bash
python core/run_pipeline.py --modality rna --resume --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

### 3.4 Run individual steps

If you want to run or re-run a specific step:

```bash
# Run only Step 06 (cell annotation)
python core/run_pipeline.py --modality rna --step 5 --config projects/rna/{dataset_id}/config_{dataset_id}.py

# Run Steps 02 through 05
python core/run_pipeline.py --modality rna --steps 2-5 --config projects/rna/{dataset_id}/config_{dataset_id}.py

# Run specific non-consecutive steps
python core/run_pipeline.py --modality rna --steps 0,2,4,11 --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

---

## 4. scRNA-seq pipeline in detail

The scRNA-seq pipeline has 12 steps (numbered 00–11), with data flowing sequentially:

```
Raw data → 00_load → 01_doublet → 02_qc
         → 03_integrate → 04_cluster → 05_annotate
         → 06_subcluster → 07_markers → 08_trajectory
         → 09_enrichment → 10_exploratory → 11_grn
```

### Step 00: Data loading

**Input**: Raw data files (format auto-detected) | **Output**: `00_raw.h5ad`

The pipeline auto-detects and loads one of the following formats:

| Format | Typical files | Common source |
|--------|--------------|---------------|
| `10X_h5` | `*filtered_feature_bc_matrix.h5` | Cell Ranger standard output |
| `10X_mtx` | `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` | Cell Ranger raw output |
| `csv_matrix` | Gene × cell count matrix (CSV/TSV/MTX) | Custom protocols, Smart-seq2, etc. |
| `h5ad` | `*.h5ad` | Pre-processed data |

> **R formats (`.rds` / `.qs`)**: Not natively supported. Use [r2h5ad](https://github.com/nanshanteahouse/r2h5ad) to convert to h5ad before loading.

Automatically handles during loading:
- **Sample/stage mapping**: Assigns each cell to its sample of origin and developmental stage based on barcode suffix (e.g., `-1`, `-2`)
- **Multi-file merging**: If a dataset contains multiple H5 files, concatenates them into a single AnnData object
- **Format compatibility**: Auto-converts legacy 2-column `genes.tsv` to standard 3-column `features.tsv`

### Step 01: Doublet detection (Scrublet)

**Input**: `00_raw.h5ad` | **Output**: `01_doublet.h5ad`

Runs Scrublet independently per sample to detect "doublets" — droplets that captured two cells instead of one.

- Large samples (>15,000 cells) are processed serially to avoid out-of-memory errors
- Small samples are processed in parallel for speed
- **Auto-skip for non-raw-counts data**: When `expression_type` is `TPM`, `FPKM`, `CPM`, or `log1p_counts`, Scrublet is disabled automatically since its negative-binomial assumption is violated for normalized data
- If Scrublet fails for a given sample, all cells in that sample are marked as non-doublet (graceful degradation — the pipeline never blocks)

Output: `doublet_scores` (doublet probability score) and `predicted_doublet` (boolean flag) columns.

### Step 02: Quality control (QC)

**Input**: `01_doublet.h5ad` | **Output**: `02_qc.h5ad`<br>
**Plots**: `{figure_dir}/02_qc/` — `nFeature_distribution.png`, `nCount_vs_nFeature.png`, `pct_mito_distribution.png`

Two modes, controlled by `use_adaptive_thresholds`:

| Mode | Config | Behavior |
|------|--------|----------|
| Hard thresholds (default) | `use_adaptive_thresholds=False` | Uses fixed cutoffs from config |
| MAD adaptive | `use_adaptive_thresholds=True` | Median ± N × MAD per metric, clipped by hard thresholds as safety caps |

Filter dimensions (all expressed as `(lo, hi)` in a thresholds dict):

1. **Remove doublets**: Discard cells flagged `predicted_doublet=True`
2. **Gene count filter**: Remove cells with too few genes (empty droplets) or too many genes (missed doublets). Default range: 500–7,500
3. **Mitochondrial filter**: Remove cells with >20% mitochondrial reads (dead/damaged cells)
4. **nCount filter** (raw_counts only): Upper-bound on `total_counts`. Skipped automatically for TPM/FPKM/CPM
5. **Complexity filter** (raw_counts only): `log10(n_genes) / log10(total_counts)` lower-bound. Skipped automatically for TPM/FPKM/CPM

**3 diagnostic plots are always generated**, annotated with the actual threshold lines (hard or MAD), providing a permanent audit trail without requiring manual inspection.

### Step 03: Normalization & batch integration

**Input**: `02_qc.h5ad` | **Output**: `03_integrated.h5ad`

This is the most critical integration step, aligning data from different samples/batches into a shared analysis space:

1. **Highly variable gene (HVG) selection**: Identifies the most informative genes (default ~4,000). Automatically tries multiple methods (`seurat_v3` → `cell_ranger` → `seurat` → manual variance-based), ensuring success on any dataset
2. **Covariate regression** (optional): Regresses out `total_counts` and `pct_counts_mt` to remove technical noise
3. **Normalization**: Library-size normalization to 10,000 total counts per cell, then log1p transformation
4. **PCA**: Principal component analysis (default 100 dimensions), with elbow plot
5. **Harmony batch correction**: Removes technical variation across samples/batches (optional, controlled via `CFG.use_harmony`)

> 💡 The full gene expression matrix is preserved in `.raw`, which downstream marker gene and differential expression analyses use — ensuring no information loss from HVG filtering.

### Step 04: Clustering & UMAP

**Input**: `03_integrated.h5ad` | **Output**: `04_clustered.h5ad`

Performs a **multi-parameter grid search** to automatically find the optimal clustering:

- Iterates over multiple `n_neighbors` values (default [15, 20, 30]) and Leiden resolutions (default [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
- Computes UMAP and Leiden clustering for each combination
- Uses **Silhouette score** as the objective metric to auto-select the best parameters
- Generates UMAP comparison plots (PDF) for all parameter combinations and a grid search summary table (CSV)

> 💡 This step is computationally heavy because it tries 3×6=18 parameter combinations. But it's worth it — you never need to manually iterate on parameters.

### Step 05: Cell type annotation

**Input**: `04_clustered.h5ad` | **Output**: `05_annotated.h5ad`

This is the core step of the pipeline: turning "Cluster 0, Cluster 1, ..." into biologically meaningful labels like "Rod Photoreceptor, Bipolar Cell, Müller Glia, ...". The pipeline supports three annotation modes, chosen automatically by priority:

#### Mode 1: KB knowledge base mode (highest accuracy)

If your tissue has a pre-built knowledge base (e.g., retina), simply set:
```python
CFG.tissue_kb = "retina"
```

The pipeline runs a sophisticated multi-layer decision process:
1. **Marker gene computation**: Wilcoxon rank-sum test to identify highly-expressed genes per cluster
2. **KB scoring**: Hypergeometric test + cosine similarity dual-scoring of cluster markers against known cell-type markers
3. **Expert rule matching**: Priority-ordered deterministic matching rules (e.g., "co-expresses RHO and PDE6A → Rod Photoreceptor")
4. **Evidence fusion**: 5-tier decision engine combining marker scores, expert rules, and hierarchical structure into a consensus annotation with confidence levels
5. **AI fallback** (optional): Low-confidence clusters are re-annotated by a large language model (LLM) based on their marker genes
6. **Quality control**: Auto-flags mitochondrial/ribosomal-dominated low-quality clusters as "Unknown"; generates an annotation quality report

#### Mode 2: AI LLM mode

If AI annotation is enabled:
```python
CFG.ai.enabled = True
CFG.ai.ai_annotation = True
```

The pipeline sends each cluster's top marker genes to an LLM, which infers cell types from biological knowledge and returns structured annotations (cell type + subtype + state + confidence + reasoning).

Supports multiple LLM backends: OpenAI API, DeepSeek, vLLM, Ollama, etc. (configured via `CFG.ai.api_base`).

#### Mode 3: Score_genes simple scoring (fallback)

If neither KB nor AI is available, the pipeline falls back to classic marker gene scoring — you only need to provide `CFG.marker_dict` (manually curated marker gene lists per cell type).

> 💡 **Annotation output includes**: `cell_type` (primary label), `cell_subtype` (subtype), `cell_state` (state), `annot_confidence` (confidence level), and `annot_reasoning` (rationale).

### Step 06: Subclustering (optional)

**Input**: `05_annotated.h5ad` | **Output**: `05_sub_{cell_type}.h5ad`

Performs fine-grained subtype analysis on a specific cell type (e.g., "Müller Glia"):

```bash
python core/run_pipeline.py --modality rna --step 7 \
    --cell-type "Müller Glia" \
    --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

This step re-runs PCA → neighbor graph → UMAP → Leiden clustering on the subset of the specified cell type, and optionally uses AI for subcluster re-annotation. Results are automatically written back to the main `05_annotated.h5ad` file's `cell_subtype` column.

### Step 07: Differential expression (three-layer DE)

**Input**: `05_annotated.h5ad` | **Output**: CSV tables + heatmaps/dotplots

Three layers of differential expression analysis:

**Layer 1 — Marker genes**: Each cell type vs. all others (Wilcoxon rank-sum test)
- Output: `marker_genes_per_group.csv`
- Identifies identity-defining markers for each cell type

**Layer 2 — Stage-wise pairwise comparison**: Same cell type across adjacent developmental stages (t-test)
- Output: `pairwise_stage_de.csv`
- Tracks transcriptional changes during development (requires `stage` annotation)
- Automatically parallelizes across comparisons

**Layer 3 — Temporal trend genes**: Gene expression Spearman correlation with developmental time
- Output: `temporal_trend_genes.csv`
- Top 20 up- and top 20 down-regulated genes per cell type
- Requires at least 3 developmental stages

> 💡 Also generates: marker gene heatmaps (top 5 per group) and known-marker dotplots.

### Step 08: Trajectory analysis (PAGA + DPT)

**Input**: `04_clustered.h5ad` (typical) | **Output**: `05_final.h5ad`

Reconstructs cellular developmental/differentiation trajectories:

1. **PAGA graph**: Builds a topological connectivity graph between cell types, revealing differentiation relationships
2. **Root cell identification**: Automatically determines the trajectory origin (priority: config-specified root type → highest root-marker expression → earliest developmental stage → first cell in dataset)
3. **Diffusion pseudotime (DPT)**: Computes pseudotime position for each cell along the PAGA topology
4. **Branch analysis**: Identifies lineage-specific genes at branch points
5. **Developmental gene visualization**: Plots known developmental genes (SOX2, PAX6, NEUROD1, etc.) along pseudotime

### Step 09: GO/KEGG pathway enrichment

**Input**: Marker gene table from Step 07 | **Output**: CSV tables + bubble plots + AI interpretation

Two complementary enrichment methods:

| Method | Principle | Output |
|--------|-----------|--------|
| **ORA** (Over-Representation Analysis) | Top marker genes → hypergeometric test for enriched pathways | `enrichment_ora.csv` |
| **GSEA** (Gene Set Enrichment Analysis) | All genes ranked by score → test whether pathway genes cluster at top/bottom | `enrichment_gsea.csv` |

Supports 200+ gene set libraries, commonly used ones include:
- `GO_Biological_Process` — Gene Ontology Biological Process
- `KEGG_2021_Human` — KEGG metabolic/signaling pathways
- `Reactome_2022` — Reactome pathway database
- `MSigDB_Hallmark_2020` — Hallmark gene sets

> 💡 If AI is enabled, enrichment results are automatically accompanied by a biological interpretation report (`ai_interpretation.txt`).

### Step 10: Exploratory analysis

**Input**: `05_annotated.h5ad` | **Output**: CSV tables + various PDF figures

Generates comprehensive summary visualizations to help you quickly understand the data:

- **Cell composition analysis**: Stacked bar charts of cell type proportions per sample/stage
- **QC metric visualization**: Distribution of n_genes, total_counts, and mitochondrial percentage on UMAP
- **Marker gene expression**: Known marker gene expression patterns on UMAP
- **Cluster statistics**: Cell counts and percentages per cluster/cell type

### Step 11: GRN regulatory network analysis (decoupler)

**Input**: `05_annotated.h5ad` | **Output**: `11_grn.h5ad` + CSV tables + heatmap

Transcription factor (TF) activity inference based on pseudobulk aggregation of annotated cell types:

1. **Pseudobulk aggregation**: Computes mean expression per `cell_type` to smooth single-cell dropout noise
2. **Regulon network**: Fetches the CollecTRI database (~1,185 TFs with signed target-gene interactions) via decoupler
3. **TF activity inference**: Runs ULM (Univariate Linear Model) — tests whether a TF's target genes are enriched among highly expressed genes in each cell type
4. **Output**:
   - `11_grn.h5ad` — pseudobulk AnnData (obs = cell types, var = genes) with `obsm['X_tf_activity']` containing TF activity scores
   - `tables/11_grn/tf_activity_per_cell_type.csv` — full TF activity matrix (cell types × TFs)
   - `tables/11_grn/tf_activity_pvals.csv` — corresponding p-values
   - `tables/11_grn/tf_target_edges.csv` — TF→target gene edges for top-variance TFs
   - `tables/11_grn/tf_target_counts.csv` — per-TF target gene count summary
   - `figures/11_grn/tf_activity_heatmap.png` — clustered heatmap of top-N variable TFs across cell types

**Config fields:**

```python
CFG.run_grn = True               # Enable/disable this step
CFG.grn_method = "decoupler"     # Method ('decoupler' only for now; pySCENIC TBD)
CFG.grn_species = "human"        # 'human' | 'mouse'
CFG.grn_n_top_regulons = 50      # Number of top-variance TFs for heatmap
CFG.grn_min_regulon_size = 5     # Minimum target genes per regulon
CFG.grn_confidence_levels = ["A","B","C"]  # DoRothEA confidence levels (if using DoRothEA)
```

> 💡 **No external database downloads required.** decoupler fetches regulon networks online on first use and caches them locally. CollecTRI (default) covers more TFs than DoRothEA and is the recommended network.

---

## 5. scATAC-seq pipeline in detail

The scATAC-seq pipeline has 10 steps (numbered 00–09):

```
Raw data → 00_load → 01_qc → 02_process → 03_cluster
         → 04_annotate → 05_marker_peaks → 06_motif
         → 07_trajectory → 08_enrichment → 09_integrate
```

### Step 00: Data loading

**Input**: Raw ATAC data | **Output**: `00_raw.h5ad`

Supports three formats:

| Format | Typical files | Notes |
|--------|--------------|-------|
| `10x_fragments` | `*fragments.tsv.gz` | ATAC fragment file; stream-imported via SnapATAC2 |
| `10x_peak_h5` | `*filtered_peak_bc_matrix.h5` | 10X peak-by-cell matrix HDF5 |
| `h5ad` | `*.h5ad` | Pre-processed AnnData |

### Step 01: QC + Peak calling

**Input**: `00_raw.h5ad` | **Output**: `01_filtered.h5ad`

1. **Fragment count filter**: Removes cells with abnormal fragment counts (too few = empty droplet, too many = doublet)
2. **MACS3 Peak calling**: Calls peaks on pseudobulk to identify open chromatin regions (standard chromosomes chr1-22, X, Y only)
3. **Peak-by-cell matrix**: Constructs a binary accessibility matrix (open=1, closed=0)
4. **Doublet detection**: Runs Scrublet on the peak matrix

### Step 02: Feature selection & dimensionality reduction

**Input**: `01_filtered.h5ad` | **Output**: `02_processed.h5ad`

1. **Remove doublets**
2. **IDF feature selection**: Selects the most informative peaks (default top 50k), reducing noise
3. **Spectral embedding**: Matrix-free Lanczos spectral decomposition (30 dimensions)
4. **KNN graph**: Builds a cell neighborhood graph in spectral space

### Step 03: Clustering & UMAP

**Input**: `02_processed.h5ad` | **Output**: `03_clustered.h5ad`

Same multi-parameter grid search strategy as the RNA pipeline (iterating `n_neighbors` × `resolution`), with Silhouette score-based auto-selection in spectral embedding space.

### Step 04: AI chromatin state annotation

**Input**: `03_clustered.h5ad` | **Output**: `04_annotated.h5ad`

1. Computes differentially accessible peaks per cluster (marker regions)
2. Sends top peaks and their associated nearby genes to an LLM
3. The LLM infers cell types / chromatin states based on the gene associations of open chromatin regions
4. AI responses are automatically cached to disk (SHA256 deduplication) — re-runs don't incur additional API calls

### Steps 05–08: Downstream analysis

| Step | Content | Output |
|------|---------|--------|
| Step 05 | Differential peak accessibility per cell type | `marker_peaks.csv` |
| Step 06 | TF binding motif enrichment (CIS-BP database) | `motif_enrichment_{cell_type}.csv` |
| Step 07 | ATAC pseudotime trajectory analysis | `07_trajectory.h5ad` |
| Step 08 | GO/KEGG enrichment on peak-associated genes | `enrichment_*.csv` |

### Step 09: RNA+ATAC integration

**Input**: ATAC `04_annotated.h5ad` + RNA h5ad (auto-discovered) | **Output**: `09_integrated.h5ad`

If you have paired multiome data (RNA-seq + ATAC-seq from the same cells):
1. Auto-discovers RNA results under the same dataset
2. Finds common cells by barcode intersection
3. Constructs a MuData multi-modal object (`rna` + `atac` modalities)
4. Runs joint PCA

> 💡 For **multiome datasets**, the preprocessor automatically generates both RNA and ATAC configs. Run both pipelines separately, then run ATAC Step 09 for automatic integration.

---

## 6. Spatial transcriptomics pipeline in detail

The spatial transcriptomics pipeline has 10 steps (numbered 00-09), designed for 10X Visium data (and extensible to other platforms):

```
Raw data → 00_load → 01_qc → 02_image → 03_normalize
         → 04_cluster → 05_annotate → 06_spatial_de
         → 07_trajectory → 08_enrichment → 09_exploratory
```

For a detailed pipeline report including real-world issues and fixes, see `notes/suggestions/spatial_<GSE_ID>.md`.

### Supported platforms

| Platform | Config value | Notes |
|----------|-------------|-------|
| 10X Visium | `"visium"` | SpaceRanger output directory or h5ad with spatial coords |
| Slide-seq | `"slideseq"` | Bead-based spatial barcoding |
| MERFISH | `"merfish"` | Imaging-based, gene panel |
| seqFISH | `"seqfish"` | Imaging-based, gene panel |

### Input format

Two data formats are supported:

| Format | Config | Input |
|--------|--------|-------|
| 10X Visium directory | `CFG.data_format = "visium"` | Directory containing `filtered_feature_bc_matrix.h5` + `spatial/` |
| Pre-built h5ad | `CFG.data_format = "h5ad"` | `.h5ad` with `obsm['spatial']` coordinates and `uns['spatial']` images |

### Step 00: Data loading

**Input**: Raw Visium data or h5ad | **Output**: `00_raw.h5ad`

- For `visium` format: uses `sq.read.visium()` with auto-detection of `library_id`
- For `h5ad` format: reads with `sc.read()` and validates spatial coordinates
- Auto-converts to sparse CSR format and ensures unique observation names
- Adds default `in_tissue` flag for Visium data

### Step 01: Quality control

**Input**: `00_raw.h5ad` | **Output**: `01_qc.h5ad`

Applies QC metrics adapted for spatial data:
- Gene count filter (default 200-7,500)
- Mitochondrial percentage filter (default <25%)
- Gene-UMI complexity filter
- Post-filtering spot count logged with low-count warning

### Step 02: Image processing

**Input**: `01_qc.h5ad` | **Output**: `02_image.h5ad`

Extracts image features from the tissue H&E/IF image:
- Auto-detects the library_id from `uns['spatial']`
- Crops image to tissue region (configurable)
- Extracts basic image features (texture, histogram) via `sq.im.process()`
- Gracefully degrades if no image is present (skips processing, continues pipeline)

### Step 03: Normalization

**Input**: `02_image.h5ad` | **Output**: `03_processed.h5ad`

- Library-size normalization to 10,000 per spot
- Log1p transformation
- Preserves raw counts in `.raw`
- PCA dimensionality reduction

### Step 04: Clustering & UMAP

**Input**: `03_processed.h5ad` | **Output**: `04_clustered.h5ad`

- Multi-resolution Leiden clustering (grid search over resolutions)
- UMAP embedding (2D)
- Spatial neighbor graph construction (`sq.gr.spatial_neighbors`)
- Grid search summary saved to `param_grid_summary.csv`

### Step 05: Cell type annotation

**Input**: `04_clustered.h5ad` | **Output**: `05_annotated.h5ad`

This is the core annotation step. Three annotation modes, chosen by priority:

#### Mode 1: KB knowledge base mode (highest accuracy)

If `CFG.tissue_kb` is set (e.g. `"retina"`), the pipeline reuses the RNA pipeline's full annotation engine:
- Computes marker genes per spatial cluster
- Scores clusters against the tissue Knowledge Base
- Applies expert deterministic rules
- Evidence fusion across scoring tiers
- AI fallback for low-confidence clusters

#### Mode 2: AI LLM mode

If `CFG.ai.enabled` and `CFG.ai.ai_annotation` are set:
- Sends per-cluster marker genes to an LLM for annotation
- Returns structured annotations (cell_type, subtype, state, confidence, reasoning)

#### Mode 3: Score_genes simple scoring (fallback)

Uses `CFG.marker_dict` for per-cluster marker gene scoring. This mode can be enriched by:
- **User-configured markers** in the config file
- **scRNA-derived markers** via Phase 1 marker-list transfer (see below)

#### Phase 1: scRNA marker-list transfer (NEW)

When scRNA-seq has been run on matched samples (same tissue, same timepoints), per-cell-type marker genes can be automatically transferred to spatial annotation:

```python
# In your spatial config:
CFG.rna_ref = "<RNA_dataset_id>"   # e.g. "GSE235585"
CFG.rna_marker_top_n = 10          # top-N markers per cell type
CFG.rna_marker_pval_threshold = 0.05
CFG.rna_marker_logfc_min = 0.0     # 0 = positive LFC only
```

How it works:
1. Auto-discovers the scRNA `marker_genes_per_group_cell_type.csv` from the matched RNA project
2. Extracts top-N significant marker genes per cell type
3. Merges into `CFG.marker_dict` (user-configured entries take priority)
4. Enriches the `score_genes_mode()` fallback without affecting KB or AI mode

> 💡 This is **marker-list transfer**, not cell-to-cell label transfer. It leverages the scRNA pipeline's already-computed differential expression results to inform spatial cluster annotation. Zero new dependencies required.

### Step 06: Spatial DE + SVG

**Input**: `05_annotated.h5ad` | **Output**: `marker_genes_per_group.csv`, `svg_rankings.csv`, `06_svg.h5ad`

- Per-cluster differential expression (Wilcoxon rank-sum)
- Moran's I spatial autocorrelation for spatially variable genes (SVG)
- Top SVG spatial scatter plots
- SVM markers marked in `adata.var['spatially_variable']`

### Step 07: Trajectory analysis

**Input**: `05_annotated.h5ad` | **Output**: `07_trajectory.h5ad`

- PAGA graph construction on spatial clusters
- Diffusion pseudotime (DPT)
- Root cell identification by marker genes or cell type

### Step 08: GO/KEGG enrichment

**Input**: `marker_genes_per_group.csv` from Step 05 | **Output**: Enrichment CSVs + bubble plots

Reuses the RNA pipeline's enrichment engine:
- ORA (over-representation analysis) via Enrichr API
- Pre-ranked GSEA (local computation)
- Supports 200+ gene set libraries (GO, KEGG, Reactome, MSigDB Hallmark, etc.)

### Step 09: Exploratory analysis

**Input**: `05_annotated.h5ad` + `06_svg.h5ad` | **Output**: Figures + CSV tables

- Spatial cell type scatter plots on tissue coordinates
- Gene expression spatial maps (top markers + SVGs)
- Spot composition statistics (cluster / cell type sizes)
- Spatial neighborhood graph summary (edges, average degree)
- UMAP summary plots

### Spatial-specific config fields

```python
CFG.spatial_platform = "visium"          # visium | slideseq | merfish | seqfish
CFG.library_id = ""                      # Visium library ID (auto-detected if empty)
CFG.crop_image = True                    # Crop image to tissue region
CFG.spatial_neighbors_n = 6              # Spatial neighbors count
CFG.spatial_neighbors_radius = 0.0       # Radius mode (0 = use n_neighbors)
CFG.run_spatial_autocorr = True          # Run Moran's I SVG detection
CFG.svg_n_top = 2000                     # Max SVGs for downstream analysis

# Phase 1: scRNA marker-list transfer
CFG.rna_ref = ""                         # scRNA project path or dataset_id
CFG.rna_marker_top_n = 10                # Top-N markers per cell type
CFG.rna_marker_pval_threshold = 0.05     # pvals_adj threshold
CFG.rna_marker_logfc_min = 0.0           # min logfoldchanges
```

---

## 6. Output files reference

After the pipeline completes, all results are organized under the dataset's `results/` directory in three subdirectories:

```
results/
├── h5ad/                          # Intermediate checkpoint files
│   ├── 00_raw.h5ad                # Raw data
│   ├── 01_doublet.h5ad            # After doublet detection
│   ├── 02_qc.h5ad                 # After QC filtering
│   ├── 03_integrated.h5ad         # After batch integration
│   ├── 04_clustered.h5ad          # After clustering + UMAP
│   ├── 05_annotated.h5ad          # After cell annotation ★
│   ├── 05_final.h5ad              # After trajectory analysis
│   └── 11_grn.h5ad               # Pseudobulk + TF activities (GRN) ★
│
├── figures/                       # Visualizations
│   ├── pca_elbow.png              # PCA elbow plot
│   ├── harmony_comparison.png     # Before/after Harmony comparison
│   ├── umap_leiden_resolutions.pdf # Multi-resolution clustering comparison
│   ├── 05_celltype.pdf            # UMAP colored by cell type
│   ├── 07_marker_heatmap.pdf      # Marker gene heatmap
│   ├── 07_dotplot.pdf             # Marker gene dotplot
│   ├── 08_pseudotime.pdf          # Pseudotime UMAP
│   ├── 08_paga_umap.pdf           # PAGA trajectory overlay
│   ├── 08_dev_genes_heatmap.pdf   # Developmental gene heatmap along pseudotime
│   ├── enrichment/                # Enrichment figures
│   │   ├── ora_*_bubble.pdf       # ORA bubble plot
│   │   └── prerank_*_bubble.pdf   # GSEA bubble plot
│   └── 10_exploratory/            # Exploratory analysis atlas
│       ├── composition_by_stage_*.png  # Cell composition stacked bars
│       └── _06_marker_dotplot.pdf      # Known marker dotplot
│
└── tables/                        # Data tables
    ├── marker_genes_per_group.csv # Marker genes (Layer 1)
    ├── pairwise_stage_de.csv      # Stage-wise pairwise DE (Layer 2)
    ├── temporal_trend_genes.csv   # Temporal trend genes (Layer 3)
    ├── branch_deg.csv             # Branch DEG
    ├── cell_type_sizes.csv        # Cell type statistics
    ├── enrichment_ora.csv         # ORA summary
    ├── enrichment_gsea.csv        # GSEA summary
    ├── 11_grn/                    # GRN analysis
    │   ├── tf_activity_per_cell_type.csv  # TF activity matrix
    │   ├── tf_activity_pvals.csv          # TF activity p-values
    │   ├── tf_target_edges.csv            # TF→target gene edges
    │   └── tf_target_counts.csv           # Per-TF target gene counts
    └── enrichment/                # Detailed enrichment results
        ├── ora_*_summary.csv
        ├── prerank_*_summary.csv
        └── ai_interpretation.txt  # AI biological interpretation
```

> 💡 The starred ★ `05_annotated.h5ad` is the most important output — it contains the final annotation labels for every cell and serves as the starting point for most downstream analyses (DE, trajectory, enrichment).

---

## 8. Practical tips

### 7.1 Check pipeline progress

```bash
# List all steps with their checkpoint files
python core/run_pipeline.py --modality rna --list --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

### 7.2 Resume from checkpoint

```bash
# Auto-detect the first incomplete step and continue from there
python core/run_pipeline.py --modality rna --resume --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

The pipeline scans checkpoint files and skips completed steps. Whether the interruption was due to network issues, memory exhaustion, or manual termination, the same command resumes correctly.

### 7.3 Re-run specific steps only

If you're not satisfied with a step's results and want to adjust parameters:

```bash
# Re-run only the annotation step (Step 06)
python core/run_pipeline.py --modality rna --step 5 --config projects/rna/{dataset_id}/config_{dataset_id}.py

# Re-run from annotation onward
python core/run_pipeline.py --modality rna --steps 6-11 --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

### 7.4 Skip time-consuming steps

If you only care about certain analyses:

```bash
# Run only from loading through annotation (first 7 steps)
python core/run_pipeline.py --modality rna --steps 0-6 --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

### 7.5 Clean up intermediate files

The pipeline generates multiple intermediate checkpoint files (each h5ad may be hundreds of MB to several GB). If disk space is tight, auto-delete upstream files after each step:

```bash
python core/run_pipeline.py --modality rna --cleanup --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

### 7.6 Subclustering

Fine-grained subtype analysis on an already-annotated cell type:

```bash
python core/run_pipeline.py --modality rna --step 7 \
    --cell-type "Müller Glia" \
    --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

You can run this for multiple cell types; results are automatically merged back into the main annotation file.

### 7.7 Multiome data workflow

```bash
# Step 1: Run RNA and ATAC separately
python core/run_pipeline.py --modality rna  --config projects/rna/{dataset_id}/config_{dataset_id}.py
python core/run_pipeline.py --modality atac --config projects/atac/{dataset_id}/config_{dataset_id}.py

# Step 2: ATAC Step 09 auto-integrates
# This is the last step of the ATAC full workflow, or run it standalone:
python core/run_pipeline.py --modality atac --step 9 --config projects/atac/{dataset_id}/config_{dataset_id}.py
```

### 7.8 Cross-modality scRNA → spatial marker transfer

If you have scRNA-seq data from matched samples, you can transfer per-cell-type marker genes into spatial annotation:

```bash
# Step 1: Run RNA pipeline first (produces marker_genes_per_group_cell_type.csv)
python core/run_pipeline.py --modality rna --config projects/rna/{rna_dataset_id}/config_{rna_dataset_id}.py

# Step 2: Configure spatial config with rna_ref pointing to the RNA dataset
# In your spatial config:
#   CFG.rna_ref = "{rna_dataset_id}"
# Step 3: Run spatial pipeline — annotation automatically uses scRNA markers
python core/run_pipeline.py --modality spatial --config projects/spatial/{spatial_dataset_id}/config_{spatial_dataset_id}.py
```

The `--list` command will show whether scRNA marker auto-discovery succeeded.

---

## 9. Configuration file deep-dive

The configuration file (`config_{dataset_id}.py`) is a Python script that controls all pipeline behavior by mutating the global `CFG` object. Below are the most commonly adjusted settings:

### 8.1 Data input

```python
# Data format (must match your actual file format)
CFG.data_format = '10X_mtx'     # Options: 10X_h5, 10X_mtx, csv_matrix, h5ad, 10x_fragments, 10x_peak_h5

# For 10X MTX format
CFG.mtx_dir = '.'               # Directory containing MTX files
CFG.mtx_prefix = 'sample_'     # MTX file prefix

# For CSV format
CFG.matrix_file = 'counts.csv.gz'
CFG.barcodes_file = 'barcodes.tsv.gz'
CFG.features_file = 'features.tsv.gz'
```

### 8.2 Sample & stage mapping

```python
# Barcode suffix → sample name
CFG.sample_map = {
    1: 'Control',
    2: 'Treatment',
}

# Sample name → developmental stage
CFG.stage_map = {
    'Control':   'E14.5',
    'Treatment': 'P0',
}
CFG.stage_order = ['E14.5', 'P0']  # Temporal order (used for time-series analysis)
```

### 8.3 Annotation (core)

```python
# Option 1: Use a knowledge base (if your tissue has KB support)
CFG.tissue_kb = "retina"

# Option 2: Manually specify marker genes (when no KB is available)
CFG.marker_dict = {
    'Rod Photoreceptor': ['RHO', 'PDE6A', 'NRL', 'GNAT1'],
    'Bipolar Cell':      ['VSX2', 'GRIK1', 'TRPM1', 'CABP5'],
    'Müller Glia':       ['GLUL', 'RLBP1', 'CLU', 'VIM'],
    # ... add more based on your tissue
}

# Option 3: Enable AI annotation
CFG.ai.enabled = True
CFG.ai.ai_annotation = True
CFG.ai.model = "deepseek-chat"              # Model name
CFG.ai.api_base = "https://api.deepseek.com/v1"
```

### 8.4 Quality control

```python
CFG.expression_type = "raw_counts"     # raw_counts | TPM | FPKM | CPM | log1p_counts
                                       # TPM/FPKM/CPM: total_counts & complexity filters auto-skipped
                                       # Scrublet auto-disabled for non-raw_counts data
CFG.min_genes = 500                    # Minimum detected genes per cell
CFG.max_genes = 7500                   # Maximum detected genes per cell (catches missed doublets)
CFG.max_pct_mito = 20.0                # Maximum mitochondrial percentage
CFG.min_genes_per_umi = 0.70           # Complexity threshold — only for raw_counts
CFG.min_cells_per_gene = 3             # Minimum cells expressing a gene
CFG.use_adaptive_thresholds = False    # True → MAD-based thresholds (auto-adapt to data)
CFG.mad_n_mads = 3.0                   # N × MAD for adaptive thresholds
CFG.qc_ncount_max_mad = 5.0            # Wider MAD multiplier for nCount upper bound
```

### 8.5 Batch correction

```python
CFG.use_harmony = True              # Enable Harmony batch correction
CFG.harmony_batch_key = "sample"    # Column to use as batch key
CFG.use_regress_out = False         # Whether to regress out total_counts and MT%
```

### 8.6 Clustering parameters

```python
CFG.n_neighbors_grid = [15, 20, 30]              # UMAP neighbor count candidates
CFG.resolution_grid = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]  # Leiden resolution candidates
CFG.best_resolution = None                        # Set to a specific value to skip grid search
```

### 8.7 Species & genome

```python
CFG.species = 'human'    # Species (affects MT gene pattern, enrichment database, etc.)
CFG.tissue = 'retina'    # Tissue name
CFG.genome = 'hg38'      # Reference genome (required for ATAC pipeline)
```

### 8.8 Spatial transcriptomics

```python
# Platform and input
CFG.spatial_platform = "visium"          # visium | slideseq | merfish | seqfish
CFG.library_id = ""                      # Visium library ID (auto-detected if empty)
CFG.data_format = "visium"               # visium | h5ad

# Image processing
CFG.crop_image = True                    # Crop image to tissue region

# Spatial graph
CFG.spatial_neighbors_n = 6              # Number of spatial neighbors
CFG.spatial_neighbors_radius = 0.0       # Radius mode (0 = use n_neighbors)

# Spatially variable genes
CFG.run_spatial_autocorr = True          # Run Moran's I
CFG.svg_n_top = 2000                     # Max SVGs for downstream analysis

# Phase 1: scRNA marker-list transfer
CFG.rna_ref = ""                         # scRNA project path or dataset_id
CFG.rna_marker_top_n = 10                # Top-N markers per cell type
CFG.rna_marker_pval_threshold = 0.05     # pvals_adj threshold
CFG.rna_marker_logfc_min = 0.0           # min logfoldchanges
```

---

## 10. FAQ

### Q1: "HDF5 file locking" error

```
OSError: Unable to open file (file locking disabled on this file system)
```

**Cause**: When running under WSL with data on the Windows filesystem (`/mnt/` paths), HDF5 file locking is incompatible.

**Solution**: Set the environment variable before running the pipeline:

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

> 💡 Consider adding this to your `~/.bashrc` or `~/.zshrc` to set it permanently.

### Q2: The generated config file has `# TODO` markers — what should I do?

The preprocessor only auto-fills what it can determine. Sections marked `# TODO` require your input:

- **`CFG.marker_dict`**: Look up known marker genes for your target tissue in the literature. If your tissue has KB support (e.g., retina), set `CFG.tissue_kb` instead
- **`CFG.sample_map`**: Extract barcode → sample mappings from GEO metadata (SRA Run Selector) or the paper's Methods section
- **`CFG.stage_map`**: If studying development, define sample → stage mappings

### Q3: A pipeline step failed — how do I recover?

```bash
# Use --resume; it automatically continues from the failed step
python core/run_pipeline.py --modality rna --resume --config projects/rna/{dataset_id}/config_{dataset_id}.py
```

`--resume` scans checkpoint files and finds the first incomplete step. Completed steps won't re-run, so recovery is fast.

### Q4: My cell annotation results aren't good. What can I try?

First, confirm you've chosen the right annotation mode:

1. **If your tissue has KB support** (e.g., retina) → use KB mode: `CFG.tissue_kb = "retina"`
2. **If you have an API key** → enable AI fallback: `CFG.ai.ai_annotation = True` (AI automatically handles low-confidence clusters in KB mode)
3. **For non-standard tissues** → manually curate `CFG.marker_dict` (reference databases like CellMarker 2.0, PanglaoDB)

> 💡 KB + AI combined gives the best results: KB provides expert-level rule matching, AI handles unknowns and precursor cell types.

### Q5: My data is TPM-formatted and QC filters out ALL cells?

TPM data has each cell's `total_counts` exactly equal to 1,000,000, which makes the complexity metric (`min_genes_per_umi`) degenerate (≈ log10(nFeature)/6).

**Solution**: Set `expression_type = "TPM"` in your config. This automatically:

1. Skips the `total_counts` filter (nCount is non-interpretable on TPM/FPKM/CPM data)
2. Skips the `log_genes_per_umi` complexity filter (equally non-interpretable)
3. Disables Scrublet (negative-binomial assumption is violated for non-count data)

```python
CFG.expression_type = "TPM"  # Auto-adjusts QC for TPM/FPKM/CPM data
```

No need to manually tune `min_genes_per_umi`.

### Q6: ATAC pipeline runs out of memory during peak merging?

In multi-sample ATAC data, independently-called peak sets per sample barely overlap, and their union can reach millions of peaks.

**Recommendation**: If possible, merge fragment files across samples and run MACS3 peak calling once via SnapATAC2 to obtain a consensus peak set.

### Q7: How much do AI API calls cost?

The pipeline has multiple layers of **caching and fallback** to control costs:

- **Disk caching**: Identical inputs (SHA256 match) never trigger repeat API calls — re-runs are free
- **Graceful degradation**: API failures automatically fall back to marker gene scoring — the pipeline never blocks
- **Predictable cost**: A single scRNA-seq annotation run typically makes 1 API call (all clusters sent at once), keeping costs very low

### Q8: My data is from a non-human species. How do I handle it?

```python
CFG.species = 'mouse'    # Supported: human, mouse, macaque, zebrafish, etc.
CFG.genome = 'mm10'      # Corresponding reference genome

# If you have an ortholog mapping file, configure it
CFG.ortholog_map = 'path/to/ortholog_map.csv'
```

> 💡 For cross-species analysis, KB annotation accuracy decreases with evolutionary distance. Consider KB + AI combined mode, letting AI compensate for KB blind spots in non-human species.

### Q9: gseapy installation fails?

GSEApy 0.11.0+ requires a Rust compiler. If `pip install gseapy` fails:

```bash
# Install Rust first
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Then retry installation
pip install gseapy
```

> 💡 If you'd rather not install Rust, use an older version: `pip install gseapy==0.10.8`

### Q10: How do I know what parameters are right for my data?

The pipeline's grid search mechanism already handles most parameter choices (Leiden resolution, UMAP neighbors). The following parameters may need manual adjustment based on your data's characteristics:

| Parameter | How to determine |
|-----------|-----------------|
| `expression_type` | Set to `"TPM"`/`"FPKM"`/`"CPM"` for non-count matrices → auto-adjusts QC |
| `min_genes` / `max_genes` | Check the distribution plots from the QC step; aim for the 1%–99% percentile range |
| `max_pct_mito` | 20% is a safe default; metabolically active tissues (e.g., cardiac muscle) may need a higher threshold |
| `use_adaptive_thresholds` | Enable for auto-adaptive QC (MAD-based); `mad_n_mads=3.0` is a good starting point |
| `n_pcs` | Inspect `pca_elbow.png` and pick the elbow point |
| `harmony_batch_key` | Usually "sample"; if your experiment spans sequencing platforms, use "platform" |

---

## Appendix: Quick command reference

### scRNA-seq

```bash
# Full workflow
python core/run_pipeline.py --modality rna --config projects/rna/{dataset_id}/config_{dataset_id}.py

# Resume from checkpoint
python core/run_pipeline.py --modality rna --resume --config projects/rna/{dataset_id}/config_{dataset_id}.py

# Single step
python core/run_pipeline.py --modality rna --step 5 --config ...

# Range of steps
python core/run_pipeline.py --modality rna --steps 4-8 --config ...

# Subclustering
python core/run_pipeline.py --modality rna --step 7 --cell-type "Cell Type Name" --config ...

# List steps
python core/run_pipeline.py --modality rna --list --config ...
```

### scATAC-seq

```bash
# Full workflow
python core/run_pipeline.py --modality atac --config projects/atac/{dataset_id}/config_{dataset_id}.py

# Resume from checkpoint
python core/run_pipeline.py --modality atac --resume --config projects/atac/{dataset_id}/config_{dataset_id}.py

# Single step
python core/run_pipeline.py --modality atac --step 4 --config ...

# RNA+ATAC integration
python core/run_pipeline.py --modality atac --step 9 --config ...
```

### Spatial transcriptomics

```bash
# Full workflow
python core/run_pipeline.py --modality spatial --config projects/spatial/{dataset_id}/config_{dataset_id}.py

# Resume from checkpoint
python core/run_pipeline.py --modality spatial --resume --config projects/spatial/{dataset_id}/config_{dataset_id}.py

# Single step
python core/run_pipeline.py --modality spatial --step 5 --config ...

# With scRNA marker transfer (Phase 1)
python core/run_pipeline.py --modality spatial --config projects/spatial/{dataset_id}/config_{dataset_id}.py
# Config must set: CFG.rna_ref = "{rna_dataset_id}"
```
