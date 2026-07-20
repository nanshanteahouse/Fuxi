# Fuxi TUI — User Guide

> For: **Single-cell omics researchers** | No programming background required | Replaces multiple CLI commands with a unified dashboard

---

## Table of Contents

1. [What is the TUI?](#1-what-is-the-tui)
2. [Launching the TUI](#2-launching-the-tui)
3. [Screen overview](#3-screen-overview)
4. [Home screen](#4-home-screen)
5. [Registry browser](#5-registry-browser)
6. [Pipeline runner](#6-pipeline-runner)
7. [Config editor](#7-config-editor)
8. [Results viewer](#8-results-viewer)
9. [Data management](#9-data-management)
10. [Keyboard shortcuts](#10-keyboard-shortcuts)
11. [Tips](#11-tips)

---

## 1. What is the TUI?

The Fuxi TUI (Terminal User Interface) is a **unified dashboard** inside your terminal. Instead of remembering and typing half a dozen different commands for registering papers, downloading data, generating configs, running pipelines, and viewing results — you navigate with your keyboard through one consistent interface.

Think of it as a "control center" for your single-cell analysis projects.

**What it replaces:**

| Old way (multiple commands) | New way (TUI) |
|---|---|
| `python -m core.paper.registry status --gse GSE123456` | Registry screen → search or click |
| `python core/geo_downloader.py --gse GSE123456` | Data Management → Download tab |
| `python core/preprocess/preprocessor.py --gse GSE123456` | Data Management → Preprocess tab |
| `python core/run_pipeline.py --modality rna --config ...` | Pipeline screen → select steps → click Run |
| Manually opening CSV files to check results | Results screen → click to view parsed tables |

---

## 2. Launching the TUI

```bash
# Activate environment first
source .venv/bin/activate

# Launch the TUI
python -m core.tui
```

The TUI opens full-screen in your terminal. Press `Ctrl+Q` to exit at any time.

---

## 3. Screen overview

The TUI has **6 screens**, accessible via keyboard shortcuts and the Home screen's quick-launch buttons:

| Screen | What it does | Shortcut |
|--------|-------------|----------|
| 🏠 Home | Landing page; select modality, quick-launch buttons | `Ctrl+H` |
| 📋 Registry | Browse/search registered papers and datasets | `Ctrl+R` |
| ⚙️ Pipeline | Select and run analysis steps | `Ctrl+P` |
| 📊 Results | View QC reports, marker genes, enrichment | `Ctrl+E` |
| 📦 Data Management | Register papers, download GEO data, generate configs | `Ctrl+D` |
| 📝 Config Editor | Edit YAML config files with a form-based editor | `Ctrl+C` |

---

## 4. Home screen

The first screen you see. It shows:

- **Modality selector** (rna / atac / spatial / bulk) — changing this affects which steps and configs are shown throughout the TUI
- **Quick-launch buttons** — jump directly to Registry, Pipeline, or Results

---

## 5. Registry browser

Shows all papers and datasets registered in Fuxi.

- **Search**: Type a GSE ID, PMID, or keyword to filter
- **Status badges**: Color-coded indicators for each dataset (🟢 downloaded, 🔵 pipeline complete, ⚪ not downloaded)
- **Click a row** to see detailed information (title, authors, species, modality, linked papers)

---

## 6. Pipeline runner

This is the **core workflow screen**:

1. **Select modality** at the top (rna/atac/spatial/bulk)
2. **Check steps**: Already-completed steps show green checkmarks
3. **Select steps to run** by checking the boxes
4. **Click "Run Selected Steps"** — the pipeline runs steps sequentially
5. **Watch progress**: The progress bar shows per-step status, the log panel streams real-time output
6. **Stop anytime** with the Stop button

**Dependency warnings**: If you select step 5 but haven't run step 3, a warning indicator appears showing which prerequisite is missing.

---

## 7. Config editor

Edit YAML configuration files without touching raw YAML:

- **Load a config**: Click "Load YAML" and select a config file from the list
- **Form fields**: Each parameter is shown as a form widget (number input, toggle, dropdown) — no YAML syntax to learn
- **Sections**: Parameters grouped by function (QC, HVG, Clustering, etc.)
- **Save**: Click "Save YAML" to write changes
- **Open in Editor**: For advanced users, "Open in $EDITOR" launches vim/nano

---

## 8. Results viewer

After pipeline steps complete, view parsed results:

- **QC Report**: Cell counts before/after filtering, median genes, etc.
- **Marker Genes**: Per-cluster marker gene table (searchable, sortable)
- **Enrichment**: Pathway enrichment results sorted by significance
- **Select config**: Pick which project's results to view

---

## 9. Data management

Three tabs for project setup:

| Tab | Purpose |
|-----|---------|
| **Register** | Add a paper by PMID → fetches metadata from NCBI automatically |
| **Download** | Enter a GSE ID → list files → select and download with progress bar |
| **Preprocess** | Detect data format → generate config file → ready for pipeline |

---

## 10. Keyboard shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+H` | Home screen |
| `Ctrl+R` | Registry browser |
| `Ctrl+P` | Pipeline runner |
| `Ctrl+E` | Results viewer |
| `Ctrl+D` | Data management |
| `Ctrl+C` | Config editor |
| `Ctrl+Q` | Quit |

---

## 11. Tips

- **Quick-launch buttons**: The Home screen's buttons are the fastest mouse route to Registry, Pipeline, and Results
- **Mouse works**: All buttons and selectors are clickable, but keyboard is faster
- **Resume works alongside TUI**: If you stop a pipeline run, you can resume from the CLI with `--resume` — the TUI and CLI share the same checkpoint system
- **Config validation**: The form editor prevents invalid values (e.g., negative cell counts)
