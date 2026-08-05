"""Unified memory settings resolution, peak estimation and guard rails.

Steps 01-03 each had their own memory machinery (01: budget-driven worker
capping, 02: live-available block sizing, 03: three-tier policy). This module
provides the single entry point ``execution.memory.{policy, budget, guard}``:

- ``resolve_memory_settings(cfg)`` -> (policy, budget_bytes, guard)
- ``estimate_step_peak(...)``   -> per-step peak RSS estimates (GB)
- ``check_memory_guard(...)``   -> warn / block / off before a run starts

Peak formulas are calibrated against measured runs (see
notes/engineering/2026-08-01_doublet_detection_overhaul.md §8 and
notes/engineering/2026-08-01_qc_step02_performance_optimization.md §8).
"""

from __future__ import annotations

import logging
import math
import re
from typing import Any

logger = logging.getLogger(__name__)

_MEM_UNIT_MULT = {"b": 1, "k": 2**10, "m": 2**20, "g": 2**30, "t": 2**40}


def resolve_memory_budget_bytes(memory_limit: str = "auto") -> int:
    """Resolve a memory budget string to bytes.

    auto -> 80% of physical RAM via psutil; explicit values accept
    64GB / 64 GiB / 512MB etc (case-insensitive). Unparsable -> auto;
    returns 0 only when even psutil fails (callers treat 0 as no budget).
    """
    text = str(memory_limit or "").strip().lower()
    if not text or text in ("auto", "0"):
        text = "auto"
    if text == "auto":
        try:
            import psutil

            return int(psutil.virtual_memory().total * 0.8)
        except Exception:
            return 0
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([a-z]+)?", text)
    if not m:
        return 0
    value = float(m.group(1))
    unit = (m.group(2) or "g").lower()
    if unit == "b":
        pass
    elif unit.endswith("ib"):
        unit = unit[:-2]  # GiB/KiB -> g/k
    elif unit.endswith("b"):
        unit = unit[:-1]  # gb/mb/kb/tb -> g/m/k/t
    mult = _MEM_UNIT_MULT.get(unit, _MEM_UNIT_MULT["g"])
    return int(value * mult)


def resolve_memory_settings(cfg: Any) -> tuple[str, int, str]:
    """Extract (policy, budget_bytes, guard) from a Config-like object.

    Prefers the nested ``execution.memory`` structure; falls back to legacy
    top-level attributes so callers survive both migration states.
    """
    exec_cfg = getattr(cfg, "execution", None)
    mem = getattr(exec_cfg, "memory", None)
    if mem is not None:
        policy = getattr(mem, "policy", "speed")
        budget = getattr(mem, "budget", "auto")
        guard = getattr(mem, "guard", "warn")
    else:
        # legacy flat fields (pre-migration safety)
        policy = getattr(exec_cfg, "memory_policy", "speed")
        budget = getattr(exec_cfg, "memory_limit", "auto")
        guard = "warn"
    budget_bytes = resolve_memory_budget_bytes(str(budget))
    return policy, budget_bytes, guard


# ═══════════════════════════════════════════════════════════════════════
#  Peak estimation (calibrated against measured runs)
# ═══════════════════════════════════════════════════════════════════════


# Step 01 — kNN result matrix is the hidden dominant term (1M cells -> ~18GB).
#   manifold = n_cells * 3 (per-neighbor features)
#   knn_index ~ manifold * 250B;  knn_result = n_cells * k_adj * 4B
#   sparse X = n_cells * n_genes * 4B * density (density = nnz / total)
#   zscore   = n_cells * n_genes * density * dtype
#   zscore   = n_cells * n_genes * density * dtype
#   worker   = base 0.8 GiB + per-worker manifold index share
def _estimate_step01_peak(n_cells: int, n_genes: int, nnz: int, budget_bytes: int) -> float:
    import os

    # kNN runs on the manifold (cells x 3 neighbours) — index ~ manifold x 250B
    # per worker; result matrix = cells x k_adj x 4B. X loads per sample group
    # (serial_threshold ~80k cells/group), never fully resident at once.
    knn_index_pw = n_cells * 3 * 250
    knn_result = n_cells * 30 * 4
    n_groups = max(1, math.ceil(n_cells / 80_000))
    group_x = nnz * 12 / n_groups  # one group resident (CSR, 12B/nnz)
    base = (knn_result + group_x) / 1e9

    n_workers = 1
    if budget_bytes:
        per_worker = 2**30 * 1.2
        n_workers = max(1, int((budget_bytes * 0.95 - 2 * 2**30) / per_worker))
    n_workers = min(n_workers, os.cpu_count() or 1)  # never exceed real cores
    if n_cells < 80_000:  # serial path below scrublet.serial_threshold
        n_workers = 1
    return base + n_workers * (knn_index_pw + 0.8 * 2**30) / 1e9


# Step 02 — peak ~ 0.85 * full-gene CSR resident + buffer; live-adaptive
#   block sizing keeps the working set at O(block) on top of X resident.
#   Measured: StressTest nnz 3.54B -> 35.1 GB peak.
def _estimate_step02_peak(nnz: int) -> float:
    return max(4.0, nnz * 12 / 1e9 * 0.85 + 3.0)


# Step 03 — speed: dense PCA (n_cells x n_genes x 4B) + .raw CSR (12B/nnz);
#   balanced/memory: CSR X + regress_out skipped -> much lower.
def _estimate_step03_peak(n_cells: int, n_genes: int, nnz: int, policy: str) -> float:
    raw_csr = nnz * 12 / 1e9  # float32 data + int64 indices
    if policy in ("balanced", "memory"):
        return max(8.0, raw_csr + n_cells * n_genes * 4 * 0.10 / 1e9 + 2.0)
    return max(10.0, raw_csr + n_cells * n_genes * 4 / 1e9 + 2.0)


def estimate_step_peak(
    step: int,
    n_cells: int,
    n_genes: int,
    nnz: int = 0,
    policy: str = "speed",
    budget_bytes: int = 0,
    approximation: str = "exact",
) -> float:
    """Estimated peak RSS (GB) for a pipeline step."""
    if step == 1:
        return _estimate_step01_peak(n_cells, n_genes, nnz, budget_bytes)
    if step == 2:
        return _estimate_step02_peak(nnz)
    if step == 3:
        return _estimate_step03_peak(n_cells, n_genes, nnz, policy)
    if step == 5:
        return _estimate_step05_peak(n_cells, n_genes, nnz, approximation)
    return 0.0


def _estimate_step05_peak(
    n_cells: int, n_genes: int, nnz: int, approximation: str = "exact"
) -> float:
    """Estimated peak RSS (GB) for step 05 (KB annotation).

    exact: raw CSR (~12B/nnz) + transposed copy + adata overhead.
    fast : cluster-downsampled workloads shrink with the sample; the
    n_cells (sampled nnz is bounded by sampling*n_genes)."""
    raw_csr = nnz * 12 / 1e9
    if approximation == "fast":
        # downsampled workloads shrink transposed/rank buffers, but the
        # raw + X views stay resident and the 05 h5ad write is unchanged
        return max(5.0, raw_csr * 1.4 + 5.0)
    # exact: raw + X CSR stay resident (~24B/nnz when X is full-gene, as in
    # stress atlases) plus the transposed copy; calibrated on the 1.05M and
    # 1.93M cell runs (formula previously underestimated peak by ~10%).
    return max(8.0, raw_csr * 2.25 + 6.0)


# ═══════════════════════════════════════════════════════════════════════
#  Guard rails
# ═══════════════════════════════════════════════════════════════════════


def check_memory_guard(
    estimates: dict[int, float],
    budget_bytes: int,
    guard: str,
    *,
    logger_obj=None,
) -> bool:
    """Compare per-step peak estimates against the budget.

    guard="warn"  -> log warning, continue (returns True)
    guard="block" -> raise RuntimeError when any estimate > 1.2*budget
    guard="off"   -> skip entirely (returns True)
    Returns True when the run may proceed.
    """
    log = logger_obj or logger
    if guard == "off" or budget_bytes <= 0:
        return True
    budget_gb = budget_bytes / 1e9
    over = {s: g for s, g in estimates.items() if g > budget_gb}
    if not over:
        if guard == "block":
            log.info("[memory-guard] all steps within budget (%.1f GB)", budget_gb)
        return True
    lines = [f"[memory-guard] estimated peak exceeds budget {budget_gb:.1f} GB:"]
    for s, g in sorted(over.items()):
        lines.append(f"  step {s:02d}: ~{g:.1f} GB")
    msg = "\n".join(lines)
    if guard == "block":
        raise RuntimeError(
            msg + "\nLower CFG.downsample.sample_keep / obs_filter or set "
            "execution.memory.policy=balanced, then retry."
        )
    log.warning(
        msg + "\n  Consider downsample (CFG.downsample) or execution.memory.policy=balanced."
    )
    return True
