# Fuxi Preprocessor — User Guide

> For: **Single-cell omics researchers** | No programming background required

---

## Table of Contents

1. [What does the preprocessor do?](#1-what-does-the-preprocessor-do)
2. [Prerequisites](#2-prerequisites)
3. [Four usage scenarios](#3-four-usage-scenarios)
4. [What you get: generated files explained](#4-what-you-get-generated-files-explained)
5. [After generation: running the full pipeline](#5-after-generation-running-the-full-pipeline)
6. [Advanced options](#6-advanced-options)
7. [FAQ](#7-faq)

---

## 1. What does the preprocessor do?

After you've downloaded single-cell data from GEO, ArrayExpress, or other sources, the preprocessor automates these tasks:

| Step | What it does | What you used to do manually |
|------|-------------|----------------------------|
| 🔍 File format detection | Auto-detects 10X MTX, HDF5, CSV matrices, ATAC fragments, etc. | Open folders and guess by filename |
| 📦 Archive extraction | Auto-extracts `.tar.gz`, `.zip`, `.gz` files | `tar -xzf` or right-click extract |
| 🧬 Modality detection | Determines scRNA-seq, scATAC-seq, or multiome | Read the paper's Methods section |
| 📋 Generate `dataset.yaml` | Creates a dataset metadata manifest (samples, file paths, formats) | Hand-edit YAML |
| ⚙️ Generate `config_*.py` | Creates a ready-to-run pipeline configuration (format-matched template) | Write ~80 lines of Python config from scratch |
| 🌐 NCBI query (optional) | Fetches title, species, SuperSeries info | Open browser, look up on GEO website |

**TL;DR: Download data → run one command → config files are generated → next step is the actual pipeline.**

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

### 2.2 Set the data root directory

```bash
export FUXI_DATA_ROOT=/data/geo_datasets    # Linux
export FUXI_DATA_ROOT=/mnt/c/geo_datasets  # WSL (adjust to your mount path)
```

> 💡 This is the top-level directory containing all your downloaded GEO folders. The preprocessor looks for datasets at `$FUXI_DATA_ROOT/GSE12345/`.

### 2.3 Prepare API key (for AI annotation in the pipeline)

```bash
export LLM_API_KEY=sk-your-api-key-here
```

> 💡 The API key is NOT embedded into generated configs. This is intentional — you must uncomment the AI section in the config file manually.


## 3. Four usage scenarios

Pick the command that matches your situation:

### Scenario 1: Standard GEO dataset (have accession ID + internet)

```bash
python core/preprocess/preprocessor.py --gse GSE00001 --query-ncbi
```

| Auto-filled | Notes |
|------------|-------|
| File format | ✅ Full auto-detect |
| Archive extraction | ✅ tar.gz / zip etc. |
| Species | ✅ From NCBI API (most accurate) |
| Dataset title | ✅ From NCBI API |
| SuperSeries detection | ✅ NCBI + directory structure + Series Matrix |
| `dataset.yaml` | ✅ Fully generated |
| `config_*.py` | ✅ Fully generated |

**Example output:**
```
============================================================
Fuxi Preprocessing: GSE00001
Data root: /data/geo_datasets
============================================================

[Phase 1] Scanning for archives...
  No archives found.
  Total files: 5

[Phase 2] Checking for SuperSeries structure...
  Not a SuperSeries (single accession).

[Phase 3] Detecting file formats...
  Inferred modality: rna

[Phase 4] Generating dataset.yaml...
  Written: projects/rna/GSE00001/dataset.yaml

[Phase 5] Generating config file...
  Written: projects/rna/GSE00001/config_GSE00001.py

============================================================
[Summary] GSE00001
============================================================
  Type:         SingleAccession
  Modality:     rna
  Data format:  10X_mtx
  Species:      homo_sapiens
  Tissue:       retina
  Elapsed:      0.0s

  Generated:
    projects/rna/GSE00001/dataset.yaml
    projects/rna/GSE00001/config_GSE00001.py

  Next steps:
    1. Review and edit the generated files
    2. Run the pipeline:
       python core/run_pipeline.py --modality rna --config projects/rna/GSE00001/config_GSE00001.py
```

---

### Scenario 2: Have accession ID + no internet (offline / air-gapped)

```bash
python core/preprocess/preprocessor.py --gse GSE00001
# Do NOT pass --query-ncbi
```

| Difference from Scenario 1 | Notes |
|---------------------------|-------|
| Species | ⚠️ Inferred from filenames only (e.g. `_human_`, `_mouse_`). Falls back to `unknown` if no match. |
| Dataset title | ❌ Not available — stays empty in `dataset.yaml` |
| SuperSeries sub-datasets | ⚠️ Detected from directory structure only; may be incomplete |

> ⚠️ **Strong recommendation**: when offline, manually verify the `tissue` and `species` fields. If the preprocessor sets them to `unknown`, edit both `dataset.yaml` and `config_*.py`.

---

### Scenario 3: In-house data + no accession + internet available

```bash
python core/preprocess/preprocessor.py \
    --input-dir /data/my_retina_organoid_experiment \
    --name retina_organoid_batch1
# --query-ncbi has no effect without a GEO accession; you can omit it
```

| Parameter | Purpose |
|-----------|---------|
| `--input-dir` | Direct path to the directory containing your data files |
| `--name` | A human-readable name for the dataset (used for output directories and filenames). Defaults to the directory basename if omitted. |

---

### Scenario 4: Air-gapped + in-house new experimental data

```bash
python core/preprocess/preprocessor.py \
    --input-dir /data/experiment_20260625 \
    --name my_new_data \
    --no-extract
```

| Difference from Scenario 1 | Notes |
|---------------------------|-------|
| Species | ⚠️ Only from filenames / metadata file content; falls back to `unknown` |
| Title | ❌ Empty |
| Accession ID | ❌ None (replaced by `--name`) |
| Archives | If your files are already in standard format, use `--no-extract` to save time |

> ⚠️ **Most critical difference**: in air-gapped environments, species detection is weakest. After generation, manually set `CFG.species` and `CFG.tissue` in the config.


## 4. What you get: generated files explained

After a successful run, two files are created under `projects/{modality}/{dataset_id}/`:

### 4.1 `dataset.yaml` — dataset manifest

```yaml
id: GSE00001
type: SingleAccession
title: ''                           # ← fill in manually
species: homo_sapiens
tissue: retina
modalities:
  - name: scRNA-seq
    status: downloaded
    format: 10X_mtx
    file_count: 5
    total_size_gb: 0.0
samples:
  - id: all
    label: ''
    rna:
      - file: GSE00001_Sample1_filtered_feature_bc_matrix.h5
        format: auto
      - file: GSE00001_Sample1_barcodes.tsv.gz
        format: auto
      ...
```

**Fields you should check / edit:**

| Field | What to do |
|-------|-----------|
| `title` | Paste the paper title from GEO website |
| `species` | If `unknown`, fill in manually (e.g. `mus_musculus`) |
| `tissue` | If `unknown`, fill in manually (e.g. `retina`, `brain`) |

### 4.2 `config_GSE00001.py` — pipeline configuration

```python
from core.config import CFG

CFG.data_format = '10X_mtx'
CFG.mtx_prefix = 'GSE00001_Sample1_'
CFG.mtx_dir = ''               # Leave empty to auto-resolve

CFG.tissue = 'retina'        # ← verify
CFG.species = 'human'

# CFG.sample_map = {        # ← fill in for multi-sample datasets
#     1: 'sample1',
# }

# CFG.marker_dict = {        # ← fill in with known markers
#     'CellTypeA': ['GENE1', 'GENE2'],
# }

# CFG.ai.enabled = True      # ← uncomment to enable AI annotation
# CFG.ai.api_base = 'https://api.deepseek.com/v1'
```

**Sections you must edit (marked with `# TODO`):**

| Section | Priority | Notes |
|---------|---------|-------|
| `CFG.marker_dict` | 🔴 Required | Fill in known marker genes for your tissue. If a KB exists for your tissue, use `CFG.tissue_kb` instead. |
| `CFG.sample_map` | 🟡 Multi-sample only | Map 10X barcode suffixes to sample names |
| `CFG.stage_map` | 🟡 Developmental only | Map samples to developmental stages |
| `CFG.tissue_kb` | 🟢 Recommended | If your tissue (e.g. `retina`, `hypothalamus`) has a KB, set this to skip manual marker curation |
| AI settings | 🟢 Recommended | Uncomment and fill in API key for LLM-assisted annotation |

> 💡 **KB-first mode**: If a tissue knowledge base exists under `rna/tissue_ontologies/`, simply set `CFG.tissue_kb = "retina"` instead of filling in `CFG.marker_dict`. KB mode is more accurate than simple gene scoring.


## 5. After generation: running the full pipeline

### 5.1 Full scRNA-seq workflow

```bash
# Run all 12 steps
python core/run_pipeline.py --modality rna --config projects/rna/GSE00001/config_GSE00001.py

# Optional: subcluster a specific cell type
python core/run_pipeline.py --modality rna --step 7 --cell-type "Müller Glia" \
    --config projects/rna/GSE00001/config_GSE00001.py
```

### 5.2 Full scATAC-seq workflow

```bash
python core/run_pipeline.py --modality atac --config projects/atac/GSE00001/config_GSE00001.py
```

### 5.3 Multiome (paired RNA + ATAC) datasets

If your dataset contains both RNA and ATAC data, the preprocessor auto-generates **two** configs:

```
projects/rna/GSE00001/config_GSE00001.py    ← RNA config
projects/atac/GSE00001/config_GSE00001.py   ← ATAC config
```

Run each modality separately first:
```bash
python core/run_pipeline.py --modality rna  --config projects/rna/GSE00001/config_GSE00001.py
python core/run_pipeline.py --modality atac --config projects/atac/GSE00001/config_GSE00001.py
```

ATAC Step 09 will then auto-discover the RNA results for integration.

### 5.4 Resume from checkpoint

```bash
python core/run_pipeline.py --modality rna --resume --config projects/rna/GSE00001/config_GSE00001.py
```

### 5.5 Run a single step

```bash
# List available steps
python core/run_pipeline.py --modality rna --list

# Run only the annotation step
python core/run_pipeline.py --modality rna --step 6 --config projects/rna/GSE00001/config_GSE00001.py
```

### 5.6 Complete workflow summary

```
1. Download data from GEO
       ↓
2. Preprocess (this guide)
   python core/preprocess/preprocessor.py --gse GSE12345
       ↓
3. Review generated dataset.yaml + config_*.py
   Edit marker_dict / tissue_kb / AI settings
       ↓
4. Run pipeline
   python core/run_pipeline.py --modality rna --config ...
       ↓
5. Analyze results
   results/figures/   → UMAP, heatmaps, trajectory plots
   results/tables/    → annotation tables, DEGs, enrichment
```


## 6. Advanced options

### Full command-line reference

```
python core/preprocess/preprocessor.py --help
```

| Flag | Purpose | When to use |
|------|---------|------------|
| `--gse GSE12345` | GEO accession ID | Scenarios 1, 2 |
| `--input-dir /path/` | Direct input directory path | Scenarios 3, 4 |
| `--name my_label` | Custom dataset identifier | Scenarios 3, 4 |
| `--data-root /path/` | Override `FUXI_DATA_ROOT` | All scenarios |
| `--query-ncbi` | Query NCBI API for metadata | Scenario 1 |
| `--dry-run` | Detect and report, write nothing | Preview before committing |
| `--force` | Overwrite existing files | Re-generation |
| `--no-extract` | Skip archive extraction | Files already extracted |
| `--modality rna\|atac` | Force a specific modality | Multi-modal datasets, separate processing |
| `--output-dir /path/` | Custom output directory | Isolate output from `projects/` tree |
| `--verbose` / `-v` | Show detailed detection info | Troubleshooting |
| `--quiet` / `-q` | Minimal output | Batch processing |

### Preview without writing files

```bash
python core/preprocess/preprocessor.py --gse GSE12345 --dry-run --verbose
```

> 💡 **Recommended**: run `--dry-run` first to confirm detection results before using `--force` for actual generation.

### Batch processing

```bash
# In bash loop
for gse in GSE00001 GSE00002 GSE00003; do
    python core/preprocess/preprocessor.py --gse "$gse" --force
done
```


## 7. FAQ

### Q1: "Directory not found" error

```
[ERROR] Directory not found: /data/geo_datasets/GSE12345
```

**Causes:**
1. `FUXI_DATA_ROOT` environment variable is not set
2. The dataset hasn't been downloaded to `$FUXI_DATA_ROOT/GSE12345/`

**Solutions:**
```bash
echo $FUXI_DATA_ROOT                    # confirm it's set
ls $FUXI_DATA_ROOT/GSE12345/            # confirm data exists
```

Or use `--input-dir` directly:
```bash
python core/preprocess/preprocessor.py --input-dir /path/to/data --name my_data
```

### Q2: The generated files have `# TODO` markers — how do I handle them?

This is **expected**. The preprocessor only auto-fills what it can determine. `# TODO` sections require your input:

- **`CFG.marker_dict`**: Look up known marker genes for your target tissue in the literature.
- **`CFG.sample_map`**: Extract barcode → sample mappings from GEO metadata.
- **`CFG.stage_map`**: If your experiment has a time-series or developmental axis, define the stage mapping.

### Q3: The preprocessor misidentified my file format. What now?

```bash
# First, inspect detailed detection
python core/preprocess/preprocessor.py --gse GSE12345 --dry-run -v

# If the format is wrong, force the correct modality
python core/preprocess/preprocessor.py --gse GSE12345 --modality atac
```

Then manually adjust `CFG.data_format` and file paths in the generated config.

### Q4: My dataset is a SuperSeries — can the preprocessor handle it?

Yes. SuperSeries are auto-detected via directory structure / Series Matrix files / NCBI API.

The generated `dataset.yaml` will list all sub-series, but **config files are not auto-generated for sub-series** — you must run the preprocessor once per sub-series.

### Q5: My data files have non-standard names. Will the preprocessor recognize them?

Partially. Supported custom formats include:
- Prefixed 10X files (e.g. `GSE12345_my_sample_filtered_feature_bc_matrix.h5`)
- Custom count matrices (e.g. `*_counts.mtx.gz`)
- Per-sample individual CSV/TSV files (GSM-style)

If detection still fails:
1. Use the templates in `templates/config_templates/` as a starting point for manual config
2. Create your own `config_*.py` under `projects/{modality}/{dataset_id}/`

### Q6: Will the preprocessor overwrite my existing config files?

**No.** Not unless you pass `--force`.

```bash
# Generate to a temporary location for review
python core/preprocess/preprocessor.py --gse GSE12345 \
    --output-dir /tmp/preprocess_output

# After review, write to the real location
python core/preprocess/preprocessor.py --gse GSE12345 --force
```

### Q7: `FUXI_DATA_ROOT` is not set — what are my options?

```bash
# Option 1: Set the env var (recommended)
export FUXI_DATA_ROOT=/data/geo_datasets

# Option 2: Pass it on every invocation
python core/preprocess/preprocessor.py --gse GSE12345 --data-root /data/geo_datasets

# Option 3: Use --input-dir (no FUXI_DATA_ROOT needed)
python core/preprocess/preprocessor.py --input-dir /data/geo_datasets/GSE12345
```

---

