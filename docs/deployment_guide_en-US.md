# Fuxi Production Deployment Guide

> For: **Users who only need to run pipelines, not develop code**
>
> If you need to modify code or run tests, see the [Environment Setup Guide](environment_setup_en-US.md).

---

## Table of Contents

1. [What is this?](#1-what-is-this)
2. [Prerequisites](#2-prerequisites)
3. [First-time Deployment (One Command)](#3-first-time-deployment-one-command)
4. [Configure Environment](#4-configure-environment)
5. [Verify Installation](#5-verify-installation)
6. [Running Pipelines](#6-running-pipelines)
7. [Daily Updates (Sync from Dev Machine)](#7-daily-updates-sync-from-dev-machine)
8. [Directory Layout](#8-directory-layout)
9. [FAQ](#9-faq)

---

## 1. What is this?

Fuxi supports two usage modes:

| Mode | Audience | Description |
|------|----------|-------------|
| **Development** | Code contributors | Full git repo + tests + linters, uses `pip install -e .` |
| **Production** ← this guide | Pipeline users | Runtime code + deps only, uses `bin/bootstrap-prod.sh` |

This guide is for **production** users: you received a copy of the Fuxi source
(via rsync, scp, or USB) and want to get it running with minimal steps.

---

## 2. Prerequisites

| Requirement | Details |
|-------------|---------|
| **OS** | Linux or WSL2 (for Windows users) |
| **Python** | 3.14+ |
| **Disk space** | ≥ 50 GB (raw data + intermediate results) |
| **Memory** | ≥ 16 GB (32 GB+ recommended for ATAC-seq) |
| **Fuxi source** | Obtained via rsync / scp / USB |

> ⚠️ **Windows users**: scATAC-seq depends on Snapatac2, which requires Linux. Use WSL2.

---

## 3. First-time Deployment (One Command)

### 3.1 Default paths

If your layout is standard (source at `/mnt/e/fuxi-prod`, venv at `~/.local/venvs/fuxi-prod`):

```bash
# Enter the source directory (synced from dev machine)
cd /path/to/fuxi-source

# One-command bootstrap
bin/bootstrap-prod.sh
```

The script performs 4 steps automatically:

1. **Sync source** → to production directory (rsync, skips existing files)
2. **Create virtualenv** → on Linux-native filesystem (`~/.local/venvs/fuxi-prod`)
3. **Install runtime deps** → `pip install .[all]` (~3-10 min, reuses pip cache)
4. **Generate config templates** → `.env`, `global.yaml`, `projects/` scaffold

### 3.2 Custom paths

```bash
FUXI_PROD_DIR=/opt/fuxi \
FUXI_PROD_VENV=/opt/fuxi-venv \
bin/bootstrap-prod.sh
```

### 3.3 Expected output

```
[bootstrap] => Production environment ready.

  Source:   /mnt/e/fuxi-prod
  Venv:     /home/user/.local/venvs/fuxi-prod

  Next steps:
    1. Edit  /mnt/e/fuxi-prod/.env          (FUXI_DATA_ROOT, LLM_API_KEY, ...)
    2. Edit  /mnt/e/fuxi-prod/global.yaml   (machine-specific settings)
    ...
```

---

## 4. Configure Environment

### 4.1 Required: `.env`

```bash
vim /mnt/e/fuxi-prod/.env
```

Key fields:

```ini
# Data root directory (required) — where GEO datasets are stored
FUXI_DATA_ROOT=/mnt/e/data

# AI annotation API key (optional)
LLM_API_KEY=sk-your-key-here

# Required for WSL2 users
HDF5_USE_FILE_LOCKING=FALSE
```

### 4.2 Optional: `global.yaml`

```bash
vim /mnt/e/fuxi-prod/global.yaml
```

Contains execution parameters (CPU cores, memory policy), clustering settings,
visualization options, etc. Defaults work for most cases — no edits needed
unless you have specific requirements.

### 4.3 Add dataset configs

Place dataset configs under the appropriate modality directory:

```bash
# RNA dataset
projects/rna/GSE123456/config_GSE123456.yaml

# ATAC dataset
projects/atac/GSE123456/config_GSE123456.yaml
```

> 📖 See the [Preprocessor Guide](preprocessor_guide_en-US.md) and
> [config templates](../templates/config_templates/) for config file format.

---

## 5. Verify Installation

```bash
VENV=~/.local/venvs/fuxi-prod
TARGET=/mnt/e/fuxi-prod

# Check that the pipeline can start
$VENV/bin/python $TARGET/core/run_pipeline.py --modality rna --list
```

If you see the step list (00-12), the environment is ready.

---

## 6. Running Pipelines

### 6.1 Basic commands

```bash
# Define convenience variables (optional, add to ~/.bashrc)
export FUXI_PYTHON=~/.local/venvs/fuxi-prod/bin/python
export FUXI_HOME=/mnt/e/fuxi-prod

# List steps
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py --modality rna --list

# Run full pipeline
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py \
    --modality rna \
    --config $FUXI_HOME/projects/rna/GSE123456/config_GSE123456.yaml

# Run a single step
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py \
    --modality rna --step 3 \
    --config $FUXI_HOME/projects/rna/GSE123456/config_GSE123456.yaml

# Resume from checkpoint
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py \
    --modality rna --resume \
    --config $FUXI_HOME/projects/rna/GSE123456/config_GSE123456.yaml
```

### 6.2 Persist environment variables

Append to `~/.bashrc`:

```bash
# ── Fuxi production environment ──
export FUXI_DATA_ROOT=/mnt/e/data
export HDF5_USE_FILE_LOCKING=FALSE
# export LLM_API_KEY=sk-...

# Convenience alias
alias fuxi='~/.local/venvs/fuxi-prod/bin/python /mnt/e/fuxi-prod/core/run_pipeline.py'
```

Then you can simply:

```bash
fuxi --modality rna --list
```

---

## 7. Daily Updates (Sync from Dev Machine)

When the dev machine has code updates, run **on the dev machine**:

```bash
# On the dev machine
cd /mnt/d/Projects/Fuxi
bin/deploy.sh
```

`deploy.sh` will:

1. Stamp a version marker via `git describe` (`version.txt`)
2. rsync changed source files to the production directory (incremental, seconds)
3. Re-run `pip install` in the production venv (updates the fuxi package, reuses cache)

> 💡 Source-only sync, skip reinstall? Use `bin/deploy.sh --no-reinstall`
>
> 💡 Preview what would change? Use `bin/deploy.sh --dry-run`

### Verify deployed version

```bash
cat /mnt/e/fuxi-prod/version.txt
# Example output: v0.2.0-5-g3a8f2b1
```

---

## 8. Directory Layout

```
Dev Machine                              Production Machine
/mnt/d/Projects/Fuxi/                    /mnt/e/fuxi-prod/
├── .git/                    not synced  ├── core/, rna/, atac/...  ← rsync'd
├── .venv/                   not synced  ├── bin/                   ← rsync'd
├── tests/                   not synced  ├── pyproject.toml         ← rsync'd
├── core/, rna/, atac/...   ── rsync ──→ ├── .env                   ← prod-only
├── bin/deploy.sh            ── rsync ──→ ├── global.yaml            ← prod-only
├── bin/bootstrap-prod.sh   ── rsync ──→ ├── projects/              ← prod-only
└── .rsync-filter                        ├── results/, logs/         ← generated
                                         └── version.txt            ← stamped

                                         ~/.local/venvs/fuxi-prod/  ← Python venv
                                         (Linux-native filesystem)
```

### Why isn't the venv inside the production directory?

WSL2 accesses `/mnt/*` (Windows mounts) via the 9p protocol. Python has
100k+ small files in `site-packages/`, and 9p is extremely slow for many
small files — pip install can take 30+ minutes. On Linux-native ext4
filesystem, it takes 3-5 minutes.

| venv location | pip install time | Daily runtime |
|---------------|:----------------:|:-------------:|
| `/mnt/e/fuxi-prod/.venv/` (Windows drive) | 30+ min | Slow |
| `~/.local/venvs/fuxi-prod/` (Linux fs) | 3-5 min | Normal |

---

## 9. FAQ

### Q1: `pip install` says "No matching distribution found" or version conflict

**Cause**: Python 3.14 is very new. Some scientific packages (e.g. liana, numba)
declare `Requires-Python <3.14` on PyPI, but actually work fine.

**Solution**: `bin/bootstrap-prod.sh` and `bin/deploy.sh` already add
`--ignore-requires-python` automatically. If you run pip install manually,
add this flag:

```bash
pip install --ignore-requires-python .[all]
```

### Q2: rsync or pip install is extremely slow (WSL2)

**Cause**: The venv or production directory is on `/mnt/*` (Windows mount).

**Solution**: Ensure the venv is on Linux-native filesystem. Check:

```bash
echo $FUXI_PROD_VENV
# Should output ~/.local/venvs/fuxi-prod or another Linux path, not /mnt/...
```

If a venv was accidentally created on /mnt, delete it and re-run bootstrap:

```bash
rm -rf /mnt/e/fuxi-prod/.venv
bin/bootstrap-prod.sh
```

### Q3: How to check what version is deployed?

```bash
cat /mnt/e/fuxi-prod/version.txt
```

Output is a git describe string (e.g. `v0.2.0-5-g3a8f2b1-dirty`). Compare
with `git log --oneline -1` on the dev machine.

### Q4: Can I modify code directly in production?

You can, but it's not recommended. Production has no `.git/`, so changes
aren't version-tracked. Best practice: modify on dev → `bin/deploy.sh` →
run in production.

### Q5: How to update all dependencies (not just fuxi itself)?

```bash
# On the dev machine
bin/deploy.sh
# ↑ deploy.sh re-runs pip install by default, which picks up new deps
```

If you only changed code without adding dependencies, use `--no-reinstall`
to skip the install step.

### Q6: Can I manually rsync instead of using deploy.sh?

Yes. deploy.sh is essentially rsync + pip install:

```bash
rsync -av --delete \
    --filter='merge .rsync-filter' \
    /mnt/d/Projects/Fuxi/ /mnt/e/fuxi-prod/

~/.local/venvs/fuxi-prod/bin/pip install --ignore-requires-python /mnt/e/fuxi-prod[all]
```

But deploy.sh also handles version stamping and error checking — the script
is recommended.

### Q7: What extras are available in production?

`bin/bootstrap-prod.sh` installs `.[all]` by default (all modalities,
no methodology packages). To add methodology packages (CellTypist, scVI, etc.)
in production:

```bash
~/.local/venvs/fuxi-prod/bin/pip install --ignore-requires-python /mnt/e/fuxi-prod[methods]
```

---

> 📖 **More docs**:
> - [Environment Setup Guide (for development)](environment_setup_en-US.md)
> - [Pipeline Guide](pipeline_guide_en-US.md)
> - [Preprocessor Guide](preprocessor_guide_en-US.md)
