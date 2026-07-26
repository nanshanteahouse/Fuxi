#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════
# deploy.sh — Deploy Fuxi from dev to prod (same-machine, no GitHub)
# ═══════════════════════════════════════════════════════════════════════
# Usage:
#   bin/deploy.sh                  # deploy to default prod
#   bin/deploy.sh --dry-run        # preview what would change
#   bin/deploy.sh --no-reinstall   # skip pip install (source-only sync)
#
# Workflow:
#   1. Stamp version.txt from git
#   2. rsync source code (.rsync-filter controls what syncs)
#   3. pip install in prod venv (non-editable; shared pip cache → fast)
# ═══════════════════════════════════════════════════════════════════════
set -euo pipefail

SOURCE="${FUXI_DEV_DIR:-/mnt/d/Projects/Fuxi}"
TARGET="${FUXI_PROD_DIR:-/mnt/e/fuxi-prod}"
VENV="${FUXI_PROD_VENV:-$HOME/.local/venvs/fuxi-prod}"
DRY_RUN=0
REINSTALL=1

for arg in "$@"; do
    case "$arg" in
        --dry-run|-n)     DRY_RUN=1 ;;
        --no-reinstall)   REINSTALL=0 ;;
        *) echo "[deploy] Unknown arg: $arg"; exit 1 ;;
    esac
done

# ── Preflight ──────────────────────────────────────────────────────────
[ -f "$SOURCE/.rsync-filter" ] || { echo "[deploy] ERROR: $SOURCE/.rsync-filter not found"; exit 1; }
[ -d "$SOURCE/core" ]          || { echo "[deploy] ERROR: $SOURCE is not a Fuxi repo"; exit 1; }

if [ "$DRY_RUN" -eq 0 ]; then
    [ -d "$TARGET" ]       || { echo "[deploy] ERROR: $TARGET not found. Run bootstrap-prod.sh first."; exit 1; }
    [ -x "$VENV/bin/pip" ] || { echo "[deploy] ERROR: $VENV/bin/pip not found. Run bootstrap-prod.sh first."; exit 1; }
fi

echo "[deploy]   source: $SOURCE"
echo "[deploy]   target: $TARGET"
echo "[deploy]   venv:   $VENV"
[ "$DRY_RUN" -eq 1 ] && echo "[deploy]   mode:   DRY RUN (no changes)"

# ── Step 1: Stamp version ──────────────────────────────────────────────
VERSION=$(cd "$SOURCE" && git describe --always --dirty --tags 2>/dev/null || echo "unknown-$(date +%s)")
if [ "$DRY_RUN" -eq 0 ]; then
    echo "$VERSION" > "$SOURCE/version.txt"
fi
echo "[deploy]   rev:    $VERSION"

# ── Step 2: rsync source ───────────────────────────────────────────────
DELETE_FLAG="--delete"
[ "$DRY_RUN" -eq 1 ] && DELETE_FLAG="--dry-run"

echo "[deploy]   syncing source..."
rsync -av $DELETE_FLAG \
    --filter="merge $SOURCE/.rsync-filter" \
    "$SOURCE/" "$TARGET/"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "[deploy] DRY RUN complete. No files changed."
    exit 0
fi

# ── Step 3: Reinstall (non-editable, shared pip cache) ─────────────────
if [ "$REINSTALL" -eq 1 ]; then
    echo "[deploy]   updating packages..."
    "$VENV/bin/pip" install --ignore-requires-python "$TARGET[all]" -q 2>&1 | grep -v '^$' || true
fi

echo "[deploy] => Done. $VERSION deployed to $TARGET"
