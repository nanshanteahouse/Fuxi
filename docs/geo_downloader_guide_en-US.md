# Fuxi GEO Downloader — User Guide

> For: **Single-cell omics researchers** | No programming background required

---

## Table of Contents (with anchors)

1. [What does the downloader do?](#1-what-does-the-downloader-do)
2. [Prerequisites](#2-prerequisites)
   - [2.1 Install a download tool (wget or curl)](#21-install-a-download-tool-wget-or-curl)
   - [2.2 Set the data root directory (FUXI_DATA_ROOT)](#22-set-the-data-root-directory-fuxi_data_root)
   - [2.3 Optional: Set NCBI API Key](#23-optional-set-ncbi-api-key)
3. [Three ways to use it](#3-three-ways-to-use-it)
   - [Method 1: Standalone CLI](#method-1-standalone-cli)
   - [Method 2: Preprocessor with --download](#method-2-preprocessor-with---download)
   - [Method 3: Registry register --pmid with --download](#method-3-registry-register---pmid-with---download)
   - [Method 3.1: Selective dataset registration](#method-31-selective-dataset-registration)
   - [Method 3.2: Deregistering (deregister)](#method-32-deregistering-deregister)
4. [Standalone CLI reference](#4-standalone-cli-reference)
   - [4.1 Full flags](#41-full-flags)
   - [4.2 Examples](#42-examples)
5. [Metadata parsing](#5-metadata-parsing)
   - [5.1 What the SOFT file provides](#51-what-the-soft-file-provides)
   - [5.2 Species normalization](#52-species-normalization)
   - [5.3 Metadata cache](#53-metadata-cache)
6. [Supplementary file listing](#6-supplementary-file-listing)
   - [6.1 How file discovery works](#61-how-file-discovery-works)
   - [6.2 RAW.tar marker](#62-rawtar-marker)
   - [6.3 Multi-file series example](#63-multi-file-series-example)
   - [6.4 What if no supplementary files are found?](#64-what-if-no-supplementary-files-are-found)
7. [Resume support](#7-resume-support)
   - [7.1 How resume works](#71-how-resume-works)
   - [7.2 Safe interruption (Ctrl+C)](#72-safe-interruption-ctrlc)
   - [7.3 Timeout](#73-timeout)
   - [7.4 Large file tips](#74-large-file-tips)
8. [Registry integration](#8-registry-integration)
   - [8.1 Automatic integration path](#81-automatic-integration-path)
   - [8.2 Empty-fields-only policy](#82-empty-fields-only-policy)
   - [8.3 CLI trigger](#83-cli-trigger)
   - [8.4 Viewing results](#84-viewing-results)
9. [Edge cases](#9-edge-cases)
   - [9.1 SuperSeries detection](#91-superseries-detection)
   - [9.2 5-digit GSE ID compatibility](#92-5-digit-gse-id-compatibility)
   - [9.3 Non-single-cell data (bulk RNA-seq, STARR-seq, etc.)](#93-non-single-cell-data-bulk-rna-seq-starr-seq-etc)
   - [9.4 Visium spatial data (GSM-level files)](#94-visium-spatial-data-gsm-level-files)
   - [9.5 Multiome data](#95-multiome-data)
   - [9.6 Filtering multi-dataset papers](#96-filtering-multi-dataset-papers)
10. [FAQ](#10-faq)
    - [Q1: Download was interrupted. Do I need to start over?](#q1-download-was-interrupted-do-i-need-to-start-over)
    - [Q2: The downloader says "Neither wget nor curl found". What do I do?](#q2-the-downloader-says-neither-wget-nor-curl-found-what-do-i-do)
    - [Q3: Can I download only specific files from a series?](#q3-can-i-download-only-specific-files-from-a-series)
    - [Q4: My dataset is a SuperSeries. Which GSE should I download?](#q4-my-dataset-is-a-superseries-which-gse-should-i-download)
    - [Q5: What's the next step after download?](#q5-whats-the-next-step-after-download)
    - [Q6: The download says "No supplementary files found". What now?](#q6-the-download-says-no-supplementary-files-found-what-now)
    - [Q7: FUXI_DATA_ROOT is not set. What should I do?](#q7-fuxi_data_root-is-not-set-what-should-i-do)
    - [Q8: Can I batch-download multiple datasets?](#q8-can-i-batch-download-multiple-datasets)
    - [Q9: The tool marked a SuperSeries but it has files. Should I still download sub-series separately?](#q9-the-tool-marked-a-superseries-but-it-has-files-should-i-still-download-sub-series-separately)

---

## 1. What does the downloader do?

The GEO Downloader fetches raw gene expression and epigenomic data from the NCBI Gene Expression Omnibus (GEO) and places it directly into your `$FUXI_DATA_ROOT` directory. It handles three tasks in one run:

| Step | What it does |
|------|-------------|
| Metadata fetch | Downloads and parses `*_family.soft.gz` — extracts title, organism, PMIDs, sample count, SuperSeries status |
| File listing | Scans the NCBI FTP server for all supplementary data files attached to the series |
| Batch download | Downloads every supplementary file using `wget` or `curl` with automatic resume support |

After the download finishes, the downloader can optionally update Fuxi's Master Registry with the fresh metadata (species, sample count, paper links).

**TL;DR: Pick a GSE accession → run one command → data lands in the right folder, ready for the preprocessor.**

---

## 2. Prerequisites

### 2.1 Install a download tool (wget or curl)

The downloader does not use Python's `urllib` for file transfers. It calls an external downloader for better resume support and progress visibility.

```bash
# Check if you have either:
which wget   # recommended
which curl   # fallback

# Install if missing (Linux / WSL)
sudo apt install wget     # recommended
sudo apt install curl
```

The downloader detects the best available tool automatically. It prefers `wget` for its cleaner resume behavior, with `curl` as a fallback. If neither is found, it prints an error and exits.

### 2.2 Set the data root directory (FUXI_DATA_ROOT)

```bash
export FUXI_DATA_ROOT=/data/geo_datasets    # Linux
export FUXI_DATA_ROOT=/mnt/c/geo_datasets  # WSL (adjust to your mount path)
```

The downloader writes files to `$FUXI_DATA_ROOT/GSEXXXXXX/`. If `FUXI_DATA_ROOT` is not set, pass `--data-root` on every invocation (see section 4).

### 2.3 Optional: Set NCBI API Key

Setting an NCBI API key raises the rate limit from 3 requests per second to 10.

```bash
export NCBI_API_KEY=your-ncbi-api-key-here
```

This is optional. The downloader works without it, just slightly slower. You can get a free API key from the NCBI website.

---

## 3. Three ways to use it

### Method 1: Standalone CLI

Directly download any GEO series by its accession. This is the most common path.

```bash
python core/geo_downloader.py --gse GSE107618
```

You get metadata printed to the terminal, then each file downloaded with a progress bar.

### Method 2: Preprocessor with `--download`

When running the preprocessor to generate pipeline configs, let it download the data first.

```bash
python core/preprocess/preprocessor.py --gse GSE107618 --modality rna --download
```

This triggers Phase 0a: the preprocessor calls the downloader, then proceeds directly to archive extraction and config generation. One command from raw accession to runnable pipeline config. If the data is already downloaded, the preprocessor skips the download phase. The registry update happens as part of the same flow, with the dataset marked as downloaded before the preprocessor begins its own phases.

### Method 3: Registry `register --pmid` with `--download`

When adding a new paper to Fuxi's Master Registry, you can auto-download all linked GEO datasets at the same time.

```bash
python -m core.registry register --pmid 31269016 --download
```

This runs the paper insight pipeline (see `paper_insights.py`), registers the paper entry, and then enters an interactive dataset selection mode. It shows each GSE dataset found in the paper's data availability statement with its SOFT metadata and lets you choose which to register.

The `register --pmid` subcommand offers three registration modes:

1. **Interactive (default)**: Shows SOFT metadata for each GSE dataset found, then prompts you to enter comma-separated numbers (e.g. `1,3`) to select which to register. Press Enter with no input to register all.
2. **`--datasets GSE1,GSE2`**: Bypass the interactive prompt and register only the specified GSE accessions. Useful for scripts and automation.
3. **`--all`**: Register every dataset found without prompting. This matches the old behavior of `add-paper`.

### Method 3.1: Selective dataset registration

When a paper mentions multiple GSE datasets, you can choose exactly which ones to register.

**Interactive mode** shows metadata for each GSE before you decide:

```
$ python -m core.registry register --pmid 31493975
Found 2 datasets:
  [1] GSE164044 — scRNA-seq of human retina (Homo sapiens, 64 samples)
  [2] GSE35156 — Hi-C of retinal cells (Homo sapiens, 6 samples)

Enter numbers to register (e.g. 1,2) or leave empty for all: 1
```

Here the user chose only the scRNA-seq dataset (GSE164044) and skipped the Hi-C dataset, since Hi-C is off-topic for a single-cell RNA pipeline.

**Non-interactive mode** with `--datasets` is useful for scripts:

```bash
python -m core.registry register --pmid 31493975 --datasets GSE164044
```

This registers only GSE164044, skipping the selection prompt entirely.

**Register everything** (backward-compatible old behavior):

```bash
python -m core.registry register --pmid 31269016 --all --download
```

### Method 3.2: Deregistering (deregister)

To remove datasets or papers from the registry:

```bash
# Remove a single dataset
python -m core.registry deregister --gse GSE35156

# Skip confirmation (for scripts)
python -m core.registry deregister --gse GSE35156 --force

# Remove a paper + cascade-delete orphan datasets
python -m core.registry deregister --pmid 31493975 --cascade

# Preview without deleting
python -m core.registry deregister --pmid 31493975 --cascade --dry-run
```

The `--cascade` flag only removes datasets that are not shared with other papers. Datasets linked to multiple papers are kept intact.

---

## 4. Standalone CLI reference

### 4.1 Full flags

```bash
python core/geo_downloader.py --gse <GSE_ID> [options]
```

| Flag | Purpose | Default |
|------|---------|---------|
| `--gse GSEXXXXXX` | GEO accession ID (required) | — |
| `--data-root /path/` | Override `FUXI_DATA_ROOT` | `$FUXI_DATA_ROOT` |
| `--dry-run` | Show what would be downloaded, don't download | Off |
| `--skip-soft` | Skip SOFT metadata fetch (use cached `.geo_meta.json`) | Off |
| `--force` | Re-download files even if they already exist | Off |
| `--quiet` / `-q` | Minimal output (errors only) | Off |

Run `python core/geo_downloader.py --help` for a full reference.

### 4.2 Examples

**Dry run (recommended first step):**

```bash
python core/geo_downloader.py --gse GSE107618 --dry-run
```

Output:

```
============================================================
  GEO Download: GSE107618
============================================================

  [METADATA] GSE107618_family.soft.gz
    Title:    Dissecting the transcriptome landscape of human neural retina and retinal
              pigment epithelium by Single-cell RNA sequenci
    Organism: Homo sapiens
    Platform: Illumina HiSeq 4000 (Homo sapiens)
    PMID:     31269016
    Samples:  64

  [SUPPL] Listing files...
    1 file(s) found (76.0 MB total)
       GSE107618_Merge.TPM.csv.gz                                 76.0 MB

  [DRY-RUN] Would download to: /mnt/e/neurobiology/GSE107618/
    → GSE107618_Merge.TPM.csv.gz  (76.0 MB)
```

**Real download:**

```bash
python core/geo_downloader.py --gse GSE107618
```

**Skip metadata (when re-downloading after a failure):**

```bash
python core/geo_downloader.py --gse GSE107618 --skip-soft
```

If the SOFT file was fetched on a previous run, it is cached as `.geo_meta.json`. Skipping re-fetch saves one NCBI request.

**Force re-download of all files:**

```bash
python core/geo_downloader.py --gse GSE107618 --force
```

By default, files that already exist and match the expected size are skipped. `--force` downloads every file from scratch.

**Custom data root:**

```bash
python core/geo_downloader.py --gse GSE107618 --data-root /mnt/e/data
```

**Multimodal dataset (RNA and ATAC in one GSE):**

```bash
python core/geo_downloader.py --gse GSE310245
```

The downloader does not distinguish between modalities. It downloads all supplementary files on the FTP server. Modality detection happens later during preprocessing.

---

## 5. Metadata parsing

Every GEO series has a SOFT format summary file (`GSEXXXXXX_family.soft.gz`). The downloader fetches and parses it automatically during Phase 1.

### 5.1 What the SOFT file provides

| Field | SOFT header | Example |
|-------|-------------|---------|
| GSE ID | From accession | GSE107618 |
| Title | `!Series_title` | "Dissecting the transcriptome landscape of human neural retina..." |
| Organism | `!Sample_organism_ch1` (aggregated) | `Homo sapiens` |
| PMID | `!Series_pubmed_id` | 31269016 (may be multiple) |
| Sample count | Count of `^SAMPLE` blocks | 64 |
| Submission date | `!Series_submission_date` | Jul 17 2024 |
| Platform | `!Platform_title` | Illumina HiSeq 4000 |
| Series type | `!Series_type` | Expression profiling by high throughput sequencing |
| Contributors | `!Series_contributor` | (list) |
| Summary | `!Series_summary` | (free text abstract) |
| SuperSeries | Detected from summary text keywords | True / False |
| Sample list | Per-sample `^SAMPLE` blocks | [{accession, title, organism}, ...] |

### 5.2 Species normalization

If all samples report the same organism, that single value is used. If samples have different organisms (e.g. a cross-species study), the downloader joins them with commas.

When the species is written to the registry (see section 8), the downloader normalizes full scientific names to short identifier slugs. For example, `Homo sapiens` becomes `human`, `Mus musculus` becomes `mouse`, and so on. This normalization ensures consistency across datasets in the registry.

### 5.3 Metadata cache

After a successful fetch, the metadata is saved as `.geo_meta.json` inside the dataset directory:

```bash
$FUXI_DATA_ROOT/GSE107618/.geo_meta.json
```

This cache is used by `--skip-soft` on subsequent runs and by the registry enrichment step (see section 8).

---

## 6. Supplementary file listing

The downloader connects to NCBI's FTP server (`ftp.ncbi.nlm.nih.gov/geo/series/.../suppl/`) to list all supplementary files attached to the series.

### 6.1 How file discovery works

1. The downloader fetches the FTP directory listing as HTML.
2. It parses each `<a href="...">` tag for file names, dates, and sizes.
3. Directory entries, parent links, `.html` files, and symbolic links pointing outside the directory are filtered out.
4. Each file's size is converted from NCBI's format (e.g. `373M`, `1.1G`) to both bytes and a human-readable string.

File sizes are displayed in a compact human-readable format:

| Actual size | Display |
|-------------|---------|
| 647 bytes | 647 B |
| 264,799 bytes | 258.6 KB |
| 1,572,864 bytes | 1.5 MB |
| 65,534,590,976 bytes | 61.0 GB |

RAW.tar files are listed first, sorted before other files. The rest are sorted alphabetically by filename.

### 6.2 RAW.tar marker

Files whose names match the pattern `RAW*.tar` (case-insensitive) are flagged with a star marker in the listing:

```
★ RAW GSE81905_RAW.tar                                           61.0 GB
   filelist.txt                                                 482 B
```

RAW.tar files typically contain the raw per-sample data (FASTQ, BAM, or unprocessed count matrices). They are always listed first in the output, sorted before other files.

### 6.3 Multi-file series example

```bash
python core/geo_downloader.py --gse GSE235583 --dry-run
```

Output:

```
============================================================
  GEO Download: GSE235583
============================================================

  [METADATA] GSE235583_family.soft.gz
    Title:    Deciphering the spatio-temporal transcriptional and chromatin
              accessibility of human retinal organoid development at the
    Organism: Homo sapiens
    Platform: Illumina NovaSeq 6000 (Homo sapiens)
    PMID:     none
    Samples:  24

  [SUPPL] Listing files...
    14 file(s) found (163.7 MB total)
    ★ RAW GSE235583_RAW.tar                                         117.0 MB
       filelist.txt                                                8.0 KB
       GSE235583_AD3_D10_counts.csv.gz                             3.5 MB
       GSE235583_AD3_D10_md.csv.gz                               317.0 KB
       GSE235583_AD3_D150_counts.csv.gz                           11.0 MB
       ...
```

### 6.4 What if no supplementary files are found?

The downloader prints a message and exits cleanly. Some GEO series store data in other repositories (dbGaP, ArrayExpress) or embed it in supplementary tables within the SOFT file itself.

---

## 7. Resume support

Large downloads can be interrupted. The downloader handles this gracefully.

### 7.1 How resume works

Both `wget` and `curl` support resume:

- **wget**: `--continue` flag (used automatically)
- **curl**: `-C -` flag (used automatically)

If a download is interrupted (Ctrl+C, network drop, timeout), simply re-run the same command. The downloader checks each file:

- If the file exists and its size matches the expected size → **skipped** (no re-download).
- If the file exists but is smaller than expected → **resumed** from where it left off.
- If `--force` is passed → **re-downloaded** from scratch regardless of existing files.

Example of a resumed download:

```
    [RESUME] GSE81905_RAW.tar (partial: 4.2 GB / 61.0 GB)
```

### 7.2 Safe interruption (Ctrl+C)

It is safe to press Ctrl+C at any point. Partially downloaded files are left in place so the next run can resume them.

### 7.3 Timeout

The downloader does not impose a download timeout. Both wget and curl have built-in network timeout mechanisms (connect timeout, read timeout) and handle slow or stalled connections on their own. If the network drops for an extended period, wget/curl will retry or exit with an error — simply re-run the command to resume from where it left off.

### 7.4 Large file tips

For datasets with RAW.tar files exceeding 10 GB (such as GSE81905 with its 61 GB RAW.tar), consider using `tmux` or `screen` to keep the download running in the background:

```bash
tmux new -s geo_download
python core/geo_downloader.py --gse GSE81905
# Ctrl+B, D to detach; come back later with:
tmux attach -t geo_download
```

The existing partial download is automatically resumed, so no progress is wasted.

---

## 8. Registry integration

After a successful download, the downloader can update Fuxi's Master Registry with metadata from the SOFT file. This keeps the registry in sync with what was actually downloaded.

### 8.1 Automatic integration path

```
download completes with zero failures
           ↓
update_registry_after_download()
           ↓
1. Dataset status → "data_downloaded"
2. enrich_dataset_from_soft()
   ├─ species       (if empty, filled from SOFT organism)
   ├─ n_samples     (if None, filled from SOFT sample count)
   ├─ paper_pmids   (append unique PubMed IDs from SOFT)
   ├─ links         (auto-create PaperDatasetLink for matching papers)
   └─ notes         (if empty, filled with SOFT summary excerpt, 200 chars max)
```

### 8.2 Empty-fields-only policy

The enrichment step never overwrites manually curated data. It only fills fields that are empty or `None`. If you have already set `species: human` in the registry, the SOFT value is ignored.

### 8.3 CLI trigger

The standalone CLI calls `update_registry_after_download()` automatically after a successful run (zero failures). No extra flags needed.

### 8.4 Viewing results

After the download and registry update complete, inspect the results with the registry CLI:

```bash
python -m core.registry report
```

Or query a specific dataset programmatically:

```bash
python -c "
from core.registry import load_master_registry
reg = load_master_registry()
ds = reg.datasets.get('GSE107618')
print(ds.model_dump_json(indent=2))
"
```

---

## 9. Edge cases

### 9.1 SuperSeries detection

A SuperSeries is a GEO meta-series that groups multiple sub-series. The downloader detects SuperSeries by scanning the SOFT summary text for keywords like "SuperSeries" combined with "composed of".

When a SuperSeries is detected, a warning is printed:

```
    ⚠  SuperSeries detected — individual sub-series should be downloaded separately
```

The parent SuperSeries entry usually contains no data files of its own. You need to download each sub-series individually.

```bash
# Example: GSE81905 is a SuperSeries
python core/geo_downloader.py --gse GSE81905 --dry-run
```

Output:

```
============================================================
  GEO Download: GSE81905
============================================================

  [METADATA] GSE81905_family.soft.gz
    Title:    Single cell RNA-sequencing of retinal bipolar cells
    Organism: Mus musculus
    Platform: Illumina NextSeq 500 (Mus musculus)
    PMID:     27565351
    Samples:  683
    ⚠  SuperSeries detected — individual sub-series should be
       downloaded separately

  [SUPPL] Listing files...
    2 file(s) found (61.0 GB total)
    ★ RAW GSE81905_RAW.tar                                           61.0 GB
       filelist.txt                                                 482 B
```

The parent does have a large RAW.tar (61 GB in this case) containing the raw data for all sub-series combined. For individual sub-series metadata, download each sub-series separately.

### 9.2 5-digit GSE ID compatibility

Older GEO accessions use 5 digits (e.g. GSE81905). Newer ones use 6 digits (e.g. GSE107618). Fuxi's GEO downloader handles both formats correctly.

NCBI's FTP directory groups accessions by prefix — using the first 2 digits for 5-digit IDs and the first 3 digits for 6-digit IDs. The downloader extracts the numeric part and computes the correct directory prefix automatically.

| Accession | Numeric part | Prefix (nnn) |
|-----------|-------------|--------------|
| GSE81905 | 81905 | GSE81nnn |
| GSE107618 | 107618 | GSE107nnn |
| GSE310245 | 310245 | GSE310nnn |

### 9.3 Non-single-cell data (bulk RNA-seq, STARR-seq, etc.)

GEO contains many types of experiments. The downloader fetches all supplementary files regardless of assay type. It does not filter for single-cell data only. After download, the Fuxi preprocessor is responsible for detecting the format and determining whether the data is suitable for single-cell analysis.

Bulk RNA-seq or other non-single-cell data may still download successfully but then fail format detection during preprocessing.

When using `register --pmid`, use interactive mode or `--datasets` to skip non-single-cell datasets during registration.

### 9.4 Visium spatial data (GSM-level files)

For 10X Visium spatial datasets, GEO often stores individual files per GSM sample rather than a single SpaceRanger output directory. The downloader downloads whatever supplementary files NCBI hosts. You may need to reorganize the files after download.

If the downloader finds only per-GSM files (e.g. `GSMXXXXXXX_*.h5`) and no `RAW.tar`, download each sample's files individually and then reconstruct the SpaceRanger output structure manually before running the preprocessor.

### 9.5 Multiome data

For multimodal datasets containing both RNA and ATAC data (such as GSE310245), the downloader fetches all supplementary files without distinguishing between modalities. Download the GSE normally:

```bash
python core/geo_downloader.py --gse GSE310245
```

The Fuxi preprocessor (`preprocessor.py`) automatically detects the multimodal structure during the next phase and generates separate RNA and ATAC pipeline configs.

### 9.6 Filtering multi-dataset papers

A paper may submit scRNA-seq, ATAC-seq, Hi-C, ChIP-seq, and other data types together. Use interactive mode or `--datasets` with `register --pmid` to register only the single-cell-relevant datasets, skipping those unsuitable for the pipeline. For example:

```bash
# Interactive: enter index numbers to select which GSEs to register
python -m core.registry register --pmid 31493975

# Non-interactive: specify directly
python -m core.registry register --pmid 31493975 --datasets GSE164044
```

---

## 10. FAQ

### Q1: Download was interrupted. Do I need to start over?

**No.** Just re-run the same command. The downloader checks existing files:

- Files that completed are **skipped**.
- Partial files are **resumed** from where they stopped.
- If a file was fully downloaded but its size is wrong (e.g. truncated), it is **re-downloaded**.

```bash
# Same command, safe to repeat
python core/geo_downloader.py --gse GSE107618
```

### Q2: The downloader says "Neither wget nor curl found". What do I do?

Install one of them:

```bash
sudo apt update
sudo apt install wget    # recommended (better resume UX)
# or
sudo apt install curl
```

The downloader uses the external tool as a subprocess. Python-only mode is not supported because native resume and progress reporting are much harder to implement robustly.

### Q3: Can I download only specific files from a series?

Not directly. The downloader downloads all supplementary files for the given GSE. If you need only specific files:

1. Run `--dry-run` first to see the full file list.
2. Download individually using `wget` or `curl` with the URLs shown in the dry-run output.
3. Or let the downloader grab everything and delete unwanted files afterward.

The URLs follow this pattern:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{GSE_ID}/suppl/{filename}
```

### Q4: My dataset is a SuperSeries. Which GSE should I download?

Download the parent SuperSeries for the combined RAW.tar (if one exists), then download each sub-series individually for their specific supplementary files and SOFT metadata.

```bash
# Example: GSE81905 is a SuperSeries with 4 sub-series
# Download the parent for the combined data
python core/geo_downloader.py --gse GSE81905

# Then download each sub-series
python core/geo_downloader.py --gse GSE81906
python core/geo_downloader.py --gse GSE81907
# ... etc.
```

Check the paper or the GEO page for the list of sub-series accessions.

### Q5: What's the next step after download?

Run the Fuxi preprocessor to generate pipeline config files:

```bash
python core/preprocess/preprocessor.py --gse GSE107618 --modality rna
```

Or, use the `--download` flag to do both in one step:

```bash
python core/preprocess/preprocessor.py --gse GSE107618 --modality rna --download
```

After the preprocessor finishes, you will have:
- `projects/rna/GSE107618/dataset.yaml` — file manifest
- `projects/rna/GSE107618/config_GSE107618.yaml` — pipeline configuration

Then run the actual pipeline:

```bash
python core/run_pipeline.py --modality rna --config projects/rna/GSE107618/config_GSE107618.yaml
```

### Q6: The download says "No supplementary files found". What now?

Some GEO series do not host supplementary data on the NCBI FTP server. Common reasons:

- Data is in a controlled-access repository (dbGaP, EGA).
- Data is hosted externally (ArrayExpress, Zenodo, the lab's own server).
- Only the processed data tables are available (included in the SOFT file itself, not as separate files).

Check the GEO page for the series (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX) and look for links under "Supplementary data" or "Supplementary file".

### Q7: FUXI_DATA_ROOT is not set. What should I do?

Set the environment variable in your shell profile, or pass `--data-root` on every invocation:

```bash
# Option 1: Set the environment variable (recommended)
export FUXI_DATA_ROOT=/data/geo_datasets

# Option 2: Pass --data-root each time
python core/geo_downloader.py --gse GSE107618 --data-root /data/geo_datasets

# Option 3: Add to .env file and source it
echo "export FUXI_DATA_ROOT=/data/geo_datasets" >> .env
source .env
```

### Q8: Can I batch-download multiple datasets?

Yes. Use a shell loop:

```bash
for gse in GSE107618 GSE118614 GSE235583; do
    python core/geo_downloader.py --gse "$gse"
done
```

To preview all datasets first:

```bash
for gse in GSE107618 GSE118614 GSE235583; do
    echo "----- $gse -----"
    python core/geo_downloader.py --gse "$gse" --dry-run --quiet
done
```

### Q9: The tool marked a SuperSeries but it has files. Should I still download sub-series separately?

Yes. The parent SuperSeries may have a combined RAW.tar containing all raw data, but each sub-series has its own SOFT metadata with sample-specific details (tissue, treatment, batch). The preprocessor and downstream analysis tools work best with per-sub-series metadata.

---

> **Tip**: Always run `--dry-run` before a large download. You see the file sizes and total volume before committing bandwidth and disk space.
