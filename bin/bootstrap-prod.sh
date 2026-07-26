#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════
# bootstrap-prod.sh — First-time production environment setup
# ═══════════════════════════════════════════════════════════════════════
# Creates:
#   1. Source code at $FUXI_PROD_DIR (Windows-accessible, /mnt/e/...)
#   2. Python venv at $FUXI_PROD_VENV (Linux-native fs, fast I/O)
#   3. .env / global.yaml from templates
#   4. projects/ directory scaffold
#
# WHY venv on Linux fs:  WSL2 /mnt/* mounts use 9p protocol — extremely
# slow for Python's 100k+ small package files.  Linux ext4 is ~50x faster.
#
# Usage:
#   bin/bootstrap-prod.sh
#   FUXI_PROD_DIR=/custom  FUXI_PROD_VENV=/custom/venv  bin/bootstrap-prod.sh
# ═══════════════════════════════════════════════════════════════════════
set -euo pipefail

SOURCE="${FUXI_DEV_DIR:-/mnt/d/Projects/Fuxi}"
TARGET="${FUXI_PROD_DIR:-/mnt/e/fuxi-prod}"
VENV="${FUXI_PROD_VENV:-$HOME/.local/venvs/fuxi-prod}"

# ── Preflight ──────────────────────────────────────────────────────────
[ -f "$SOURCE/.rsync-filter" ] || { echo "[bootstrap] ERROR: $SOURCE/.rsync-filter not found"; exit 1; }
[ -d "$SOURCE/core" ]          || { echo "[bootstrap] ERROR: $SOURCE is not a Fuxi repo"; exit 1; }

if [ -d "$VENV/lib" ] || [ -d "$TARGET/projects" ]; then
    echo "[bootstrap] WARNING: looks already initialized ($VENV or $TARGET/projects exists)."
    echo "             Continue anyway? [y/N]"
    read -r confirm
    [ "$confirm" = "y" ] || [ "$confirm" = "Y" ] || { echo "[bootstrap] Aborted."; exit 0; }
fi

echo "[bootstrap]   source:  $SOURCE"
echo "[bootstrap]   target:  $TARGET"
echo "[bootstrap]   venv:    $VENV  (Linux-native)"

mkdir -p "$TARGET" "$VENV"

# ── Step 1/4: Sync source code ─────────────────────────────────────────
echo "[bootstrap]   [1/4] syncing source code..."
rsync -av --filter="merge $SOURCE/.rsync-filter" "$SOURCE/" "$TARGET/"

# ── Step 2/4: Create venv (Linux filesystem) ───────────────────────────
echo "[bootstrap]   [2/4] creating venv on Linux fs..."
if [ ! -x "$VENV/bin/python" ]; then
    python3 -m venv "$VENV"
fi
"$VENV/bin/pip" install --upgrade pip -q

# ── Step 3/4: Install runtime deps ─────────────────────────────────────
echo "[bootstrap]   [3/4] installing runtime dependencies..."
# --ignore-requires-python: many scientific pkgs declare <3.14 but work fine on 3.14
"$VENV/bin/pip" install --ignore-requires-python "$TARGET[all]"

# ── Step 4/4: Config templates & scaffold ──────────────────────────────
echo "[bootstrap]   [4/4] setting up configuration..."
[ -f "$TARGET/.env" ]          || cp "$TARGET/.env.example" "$TARGET/.env"
[ -f "$TARGET/global.yaml" ]   || cp "$TARGET/global.example.yaml" "$TARGET/global.yaml"
mkdir -p "$TARGET/projects"/{rna,atac,spatial,bulk}

# ── Done ───────────────────────────────────────────────────────────────
echo ""
echo "[bootstrap] => Production environment ready."
echo ""
echo "  Source:   $TARGET"
echo "  Venv:     $VENV"
echo ""
echo "  Next steps:"
echo "    1. Edit  $TARGET/.env          (FUXI_DATA_ROOT, LLM_API_KEY, ...)"
echo "    2. Edit  $TARGET/global.yaml   (machine-specific settings)"
echo "    3. Add configs to  $TARGET/projects/{modality}/{GSE_ID}/"
echo "    4. Run pipeline:"
echo "         $VENV/bin/python $TARGET/core/run_pipeline.py --modality rna --list"
echo ""
echo "  Future updates:  bin/deploy.sh"
