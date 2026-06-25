# Fuxi Environment Setup Guide

> For: **Single-cell omics researchers** | No programming background required

---

## Table of Contents

1. [What you'll need](#1-what-youll-need)
2. [Step 1: Install Python](#2-step-1-install-python)
3. [Step 2: Create a virtual environment](#3-step-2-create-a-virtual-environment)
4. [Step 3: Install dependencies](#4-step-3-install-dependencies)
5. [Step 4: Configure environment variables](#5-step-4-configure-environment-variables)
6. [Verifying the installation](#6-verifying-the-installation)
7. [Platform-specific notes](#7-platform-specific-notes)
8. [FAQ](#8-faq)

---

## 1. What you'll need

Before starting, make sure you have the following:

| Requirement | Details | How to get it |
|-------------|---------|---------------|
| **Operating System** | Linux (recommended), WSL2 (Windows users), or macOS | Windows 10/11 has built-in WSL2 support |
| **Python** | 3.14 or later | [python.org](https://www.python.org/downloads/) or your system package manager |
| **pip** | Python package manager (ships with Python) | Usually no separate installation needed |
| **Git** | For cloning the project repository | [git-scm.com](https://git-scm.com/) |
| **Disk space** | ≥ 50 GB (for raw data and intermediate results) | — |
| **Memory** | ≥ 16 GB (32 GB+ recommended for ATAC-seq) | — |

> 💡 **Heads-up for Windows users**: scATAC-seq analysis depends on Snapatac2, which requires Linux. If you plan to process ATAC data, **use WSL2**. scRNA-seq analysis works fine on native Windows Python.

---

## 2. Step 1: Install Python

### 2.1 Check if Python is already installed

Open a terminal and run:

```bash
python3 --version
```

If you see something like `Python 3.14.x`, you're all set — skip ahead to [Step 2](#3-step-2-create-a-virtual-environment).

If you see `command not found` or a version below 3.14, follow the instructions below.

### 2.2 Linux (Ubuntu / Debian)

```bash
# Add deadsnakes PPA (provides the latest Python releases)
sudo apt update
sudo apt install software-properties-common -y
sudo add-apt-repository ppa:deadsnakes/ppa -y
sudo apt update

# Install Python 3.14
sudo apt install python3.14 python3.14-venv python3.14-dev -y
```

### 2.3 WSL2 (Windows)

In your WSL2 terminal, follow the Linux instructions above. If you haven't installed WSL2 yet:

```powershell
# Run in Windows PowerShell (as Administrator)
wsl --install -d Ubuntu-24.04
```

Restart your computer after installation, then launch the "Ubuntu" app to enter the Linux terminal.

### 2.4 macOS

```bash
# Using Homebrew
brew install python@3.14
```

---

## 3. Step 2: Create a virtual environment

> 💡 **Why a virtual environment?** It isolates this project's dependencies from your system Python and other projects, preventing version conflicts. This is standard practice for all Python projects.

### 3.1 Get the project code

```bash
git clone <repository-url>
cd <project-directory>
```

### 3.2 Create and activate the virtual environment

```bash
# Create the virtual environment in the project root (named .venv)
python3.14 -m venv .venv
```

> ⚠️ Use `python3.14` rather than `python3` to ensure the correct Python version.

**Activate the virtual environment:**

```bash
# Linux / WSL2 / macOS
source .venv/bin/activate
```

Once activated, your terminal prompt will be prefixed with `(.venv)`:

```
(.venv) user@host:~/project$
```

### 3.3 Confirm everything looks right

```bash
# Check Python version
python --version
# Should output: Python 3.14.x

# Confirm pip is inside the virtual environment
which pip
# Should output: .../project-directory/.venv/bin/pip
```

### 3.4 Deactivate and re-enter

```bash
# Exit the virtual environment
deactivate

# Re-activate next time you work on the project
cd <project-directory>
source .venv/bin/activate
```

> 💡 You need to run `source .venv/bin/activate` each time you open a new terminal window. This is normal — your environment hasn't broken.

---

## 4. Step 3: Install dependencies

### 4.1 Basic installation

Make sure the virtual environment is active (your prompt starts with `(.venv)`), then:

```bash
pip install -r requirements.txt
```

This installs all required Python packages:

| Category | Key packages | Purpose |
|----------|-------------|---------|
| **Data processing** | `numpy`, `scipy`, `pandas` | Numerical computation, sparse matrices, data tables |
| **Visualization** | `matplotlib` | UMAP plots, heatmaps, violin plots, etc. |
| **Parallel computing** | `joblib` | Speed up differential expression and other steps |
| **scRNA-seq** | `scanpy`, `scrublet`, `leidenalg` | Single-cell RNA analysis core |
| **Batch correction** | `harmony-pytorch` | Remove batch effects between samples |
| **Functional enrichment** | `gseapy` | GO / KEGG pathway enrichment analysis |
| **scATAC-seq** | `snapatac2`, `macs3` | Single-cell ATAC analysis core |
| **Multi-omics integration** | `muon`, `mudata` | Joint RNA + ATAC analysis |
| **AI annotation** | `openai` | LLM-assisted cell type annotation |

### 4.2 How long will it take?

A fresh install typically takes **5–15 minutes**, depending on your network speed and machine.

### 4.3 Verify the installation

```bash
# Confirm core packages imported successfully
python -c "import scanpy; print('scanpy', scanpy.__version__)"
python -c "import snapatac2; print('snapatac2', snapatac2.__version__)"
python -c "import numpy; print('numpy', numpy.__version__)"
```

If these commands complete without errors, your dependencies are ready.

### 4.4 Keeping up to date

When the project code is updated, you may need to install new dependencies:

```bash
git pull                       # Pull the latest code
pip install -r requirements.txt # Install any newly added packages
```

---

## 5. Step 4: Configure environment variables

### 5.1 Required vs. optional

| Variable | Required? | Purpose |
|----------|:---------:|---------|
| `FUXI_DATA_ROOT` | ✅ Required | Root directory containing your raw downloaded datasets |
| `LLM_API_KEY` | ❌ Optional | API key for AI-assisted annotation |
| `HDF5_USE_FILE_LOCKING` | ⚠️ Recommended | WSL2 users must set this to `FALSE` |

### 5.2 Set the data root directory

```bash
# Point this to where your datasets live (replace with your actual path)
export FUXI_DATA_ROOT=/data/geo_datasets

# Example: if your data is under /home/user/geo_data
export FUXI_DATA_ROOT=/home/user/geo_data
```

The expected directory structure under `FUXI_DATA_ROOT`:

```
$FUXI_DATA_ROOT/
├── GSE12345/          # One folder per dataset (folder name = dataset ID)
│   ├── dataset.yaml   # Dataset metadata file
│   ├── *.h5           # Raw HDF5 files
│   └── *.csv.gz       # Raw CSV files
├── GSE23456/
│   └── ...
└── ...
```

> 💡 The pipeline reads raw data from `$FUXI_DATA_ROOT/<dataset-id>/`.

### 5.3 Set the AI API key (optional)

If you want to use LLM-assisted cell type annotation:

```bash
export LLM_API_KEY=sk-your-api-key-here
```

> 💡 Most analysis steps work fine without this. AI annotation is an optional enhancement.

### 5.4 Extra step for WSL2 users

Accessing the Windows filesystem from WSL2 can trigger HDF5 file-locking errors. Add this:

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

### 5.5 Making the variables permanent

The `export` commands above only last for the current terminal session. To make them apply automatically on every new terminal, append the following to your `~/.bashrc` file:

```bash
# ── Fuxi pipeline environment ──
export FUXI_DATA_ROOT=/data/geo_datasets    # Replace with your actual path
export HDF5_USE_FILE_LOCKING=FALSE           # Needed on WSL2
# export LLM_API_KEY=sk-...                  # Optional: uncomment and fill in your key
```

Save the file, then run `source ~/.bashrc` to apply immediately.

---

## 6. Verifying the installation

Once all the above steps are complete, run this full check to confirm your environment is ready:

```bash
# 1. Make sure you're in the project directory
cd <project-directory>

# 2. Confirm the virtual environment is active
echo $VIRTUAL_ENV
# Should output the .venv path inside your project directory

# 3. Confirm the data root variable is set
echo $FUXI_DATA_ROOT
# Should output the directory path you configured

# 4. Confirm all core packages load correctly
python -c "
import scanpy as sc
import snapatac2 as snap
import anndata
import numpy as np
import pandas as pd
import matplotlib
print('All core packages loaded successfully')
print(f'  Scanpy:    {sc.__version__}')
print(f'  Snapatac2: {snap.__version__}')
print(f'  AnnData:   {anndata.__version__}')
print(f'  NumPy:     {np.__version__}')
print(f'  Pandas:    {pd.__version__}')
print(f'  Matplotlib:{matplotlib.__version__}')
"
```

If you see `All core packages loaded successfully` followed by version numbers, your environment is ready.

---

## 7. Platform-specific notes

### 7.1 Linux (native)

✅ **Best experience**. All features (scRNA-seq + scATAC-seq + multi-omics integration) are fully supported. Ubuntu 22.04 or later is recommended.

### 7.2 WSL2 (Windows users)

✅ **Full functionality**. All analyses work identically to native Linux. Additional considerations:

- Store data inside WSL2 (e.g., `/home/user/data/`) rather than under `/mnt/c/` or `/mnt/e/`. Cross-filesystem HDF5 access is significantly slower.
- Always set `HDF5_USE_FILE_LOCKING=FALSE`.
- If your data is already on a Windows drive, copy it in: `cp -r /mnt/e/geo_data ~/data/`

### 7.3 macOS

⚠️ **Partial functionality**. Snapatac2 depends on Linux-specific features and may not install or run on macOS. scRNA-seq analysis is fully supported. For ATAC analysis, use a Linux VM or cloud server.

### 7.4 Windows (native, without WSL)

⚠️ **scRNA-seq only**. Snapatac2 does not support Windows. If you only need scRNA-seq, you can run on native Windows Python. The installation steps are the same as Linux, except activation uses:

```cmd
.venv\Scripts\activate
```

---

## 8. FAQ

### Q1: `pip install -r requirements.txt` fails — a package won't install

Common package-specific issues and solutions:

**`snapatac2` fails to install:**

```bash
# Confirm Python 3.14+
python --version

# On macOS: snapatac2 is Linux-only. Skip ATAC packages:
pip install scanpy anndata numpy scipy pandas matplotlib joblib scrublet leidenalg gseapy openai python-dotenv harmony-pytorch
```

**`harmony-pytorch` fails to install:**

```bash
# PyTorch may need a separate installation first
pip install torch
pip install harmony-pytorch
```

**`gseapy` fails to install:**

```bash
# gseapy needs a Rust compiler
# Ubuntu / Debian:
sudo apt install cargo rustc -y
pip install gseapy
```

**Other packages fail:**

Check the following:
1. Python development headers installed: `sudo apt install python3.14-dev`
2. Build toolchain available: `sudo apt install build-essential`
3. Network connectivity (packages are downloaded from PyPI)

### Q2: The `python` command is not found after activating the virtual environment

```bash
# Confirm .venv was created
ls -la .venv/bin/python*

# If missing, recreate
python3.14 -m venv .venv --clear
```

### Q3: "Unable to open file" or HDF5 errors in WSL2

```bash
# Cause: HDF5 file locking breaks on cross-filesystem access
# Fix:
export HDF5_USE_FILE_LOCKING=FALSE

# Or move your data inside WSL2
cp -r /mnt/e/geo_data ~/data/
export FUXI_DATA_ROOT=~/data
```

### Q4: After activation, `python` still points to the system Python

```bash
# Check what python resolves to
which python
# Should output: /path/to/project/.venv/bin/python

# If it still points to the system Python:
hash -r                  # Clear bash command cache
source .venv/bin/activate
```

### Q5: How do I completely remove and rebuild the virtual environment?

```bash
deactivate                # Exit first
rm -rf .venv              # Remove the old environment
python3.14 -m venv .venv  # Recreate
source .venv/bin/activate
pip install -r requirements.txt
```

### Q6: Running out of memory during installation or pipeline execution

Recommendations:
- Ensure at least **16 GB of physical memory**; 32 GB+ recommended for ATAC-seq
- Close other memory-intensive applications
- For large scRNA-seq datasets (>50K cells), enable downsampling in the config file
- WSL2 users can adjust the memory cap in `%UserProfile%\.wslconfig`

### Q7: Do I need a GPU?

**No.** All pipeline computations run on CPU. A GPU can accelerate the Harmony batch-correction step (via its PyTorch backend), but it is not required.

### Q8: Can I use conda / Anaconda instead of venv?

Yes. If you're comfortable with conda, replace the virtualenv steps with:

```bash
conda create -n fuxi python=3.14 -y
conda activate fuxi
pip install -r requirements.txt
```

> 💡 These docs use venv as the primary path because conda channels may lag behind Python 3.14, and Snapatac2 recommends pip installation. Both approaches ultimately do `pip install` — pick whichever you're more familiar with.

---

## Environment setup checklist

Before running the pipeline, confirm each item:

- [ ] Python 3.14+ is installed
- [ ] Virtual environment `.venv` is created and activated
- [ ] `pip install -r requirements.txt` completed successfully
- [ ] `FUXI_DATA_ROOT` is set and the directory exists
- [ ] (WSL2) `HDF5_USE_FILE_LOCKING=FALSE` is set
- [ ] (Optional) `LLM_API_KEY` is set if you plan to use AI annotation
- [ ] Verification script runs without errors (see [Section 6](#6-verifying-the-installation))

Once all items are checked off, you're ready to move on to the Preprocessor User Guide and start processing datasets.

---

> 📖 **Next step**: With your environment ready, read the *Fuxi Preprocessor — User Guide* to learn how to generate pipeline configuration files from your raw downloaded data.
