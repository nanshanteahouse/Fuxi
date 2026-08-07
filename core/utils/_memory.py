"""Unified memory settings resolution, peak estimation and guard rails.

Steps 01-03 each had their own memory machinery (01: budget-driven worker
capping, 02: live-available block sizing, 03: three-tier policy). This module
provides the single entry point ``execution.memory.{policy, budget, guard}``:

- ``resolve_memory_settings(cfg)`` -> (policy, budget_bytes, guard)
- ``estimate_step_peak(...)``   -> per-step peak RSS estimates (GB)
- ``check_memory_guard(...)``   -> warn / block / off before a run starts

Peak formulas are calibrated against measured runs (step 00: 5 datasets
 1.2M cells @ 34.11 GiB, 70.5k @ 7.49 GiB, 51k @ 3.42 GiB, 32k @ 2.14 GiB,
 137.5k @ 3.76 GiB; step 01: 155k cells @ 2.8 GB; step 02: 1.05M cells
 @ 33 GB, 1.68M cells @ 35 GB; step 04: 110k @ 12 GB, 620k @ 21 GB,
 1.05M @ 40.5 GB, 1.676M @ <100 GB (run OK on a 100 GB budget); step
 06: subcluster subsets 2.5k-61k cells; step 07: 9 runs from 2.3k to
 1.05M cells, RSS 1.8-16.3 GB).
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
#   manifold = n_cells * 3 (per-neighbor features);  k_adj = 1.5 * sqrt(n_cells)
#   knn_index ~ manifold * 250B (25B/point/tree x 10 trees; 411k -> 98 MiB)
#   knn_result = manifold * k_adj * 4B — O(cells^1.5), the hidden big one
#   zscore_main = max_group * n_genes * 4B * 2 (sparse_zscore densifies the
#     largest sample group: result + temp, float32; float64 doubles it)
#     max_group ~ 2 x mean group size (StressTest 1.97M -> 157k measured)
#   worker = base 0.8 GiB + per-worker manifold index share
#   Peak ~ max(kNN stage, zscore stage) + one group X + worker shares
#   Anchors (2026-08-01 overhaul): StressTest 1.97M -> 52 GB (f32) / 71 GB
#   (f64); 155k-cell runs ~3 GB (small-data floor).
def _estimate_step01_peak(n_cells: int, n_genes: int, nnz: int, budget_bytes: int) -> float:
    import os

    manifold = n_cells * 3
    k_adj = 1.5 * math.sqrt(n_cells)
    knn_index_pw = manifold * 250
    knn_result = manifold * k_adj * 4

    n_groups = max(1, math.ceil(n_cells / 80_000))
    max_group = 2 * n_cells / n_groups  # largest sample group ~ 2x the mean
    group_x = nnz * 12 / n_groups  # one group resident (CSR, 12B/nnz)

    # sparse_zscore densifies the largest group: dense result + temp, float32.
    # 157k x 36,601 -> 46 GB (f32) / 92 GB (f64) measured.
    zscore_main = max_group * n_genes * 4 * 2
    stage_gb = max(knn_result, zscore_main) / 1e9  # kNN and zscore barely overlap

    n_workers = 1
    if budget_bytes:
        per_worker = 2**30 * 1.2
        # mem_cap = (budget - main-process peak) / per-worker peak (report
        # model); the scheduler caps ~12 parallel sample groups in practice
        # (round 7 used 8 buckets, Lobe_Neurons measured ~13 workers).
        main_gb = stage_gb + group_x / 1e9 + 0.8
        n_workers = max(1, int((budget_bytes * 0.95 - main_gb * 1e9) / per_worker))
    n_workers = min(n_workers, os.cpu_count() or 1, 12)  # never exceed real cores
    if n_cells < 80_000:  # serial path below scrublet.serial_threshold
        n_workers = 1
    return stage_gb + group_x / 1e9 + 0.8 + n_workers * (knn_index_pw + 0.8 * 2**30) / 1e9


# Step 02 — streaming (backed) QC keeps X out of RAM, so peak is budget-
#   driven, not matrix-sized: ~0.4 x the memory budget (write staging) +
#   536 B/cell + 1.0 GiB base.  Calibrated 2026-08-05 on 63 runs (496k ->
#   42.4 GiB on 128G, Lobe_Neurons 1.05M -> 32.9, StressTest 1.97M -> 34.4,
#   all +-7%).  Fallback (no budget): old nnz-resident model.
def _estimate_step02_peak(n_cells: int, nnz: int, budget_bytes: int = 0) -> float:
    if budget_bytes > 0:
        return max(4.0, budget_bytes / 1e9 * 0.4 + n_cells * 536 / 1e9 + 1.0)
    return max(4.0, nnz * 12 / 1e9 * 0.85 + 3.0)


# Step 03 — speed: dense HVG PCA (n_cells x n_genes x 4B) + .raw CSR (12B/nnz);
#   balanced/memory: CSR X + regress_out skipped -> much lower.  stream_raw
#   keeps the full-gene matrix off-RAM (chunked reads), so the raw term is
#   scaled 0.95 as a conservative residual (StressTest 1.97M: 71.8 est vs
#   71.1 measured; Lobe_Neurons 1.05M: 43.2 vs 46.9).
def _estimate_step03_peak(n_cells: int, n_genes: int, nnz: int, policy: str) -> float:
    raw_csr = nnz * 12 / 1e9  # float32 data + int64 indices
    if policy in ("balanced", "memory"):
        return max(8.0, raw_csr + n_cells * n_genes * 4 * 0.10 / 1e9 + 2.0)
    return max(10.0, raw_csr * 0.95 + n_cells * n_genes * 4 / 1e9 + 2.0)


def _estimate_step03_peak(n_cells: int, n_genes: int, nnz: int, policy: str) -> float:
    raw_csr = nnz * 12 / 1e9  # float32 data + int64 indices
    if policy in ("balanced", "memory"):
        return max(8.0, raw_csr + n_cells * n_genes * 4 * 0.10 / 1e9 + 2.0)
    return max(10.0, raw_csr * 0.95 + n_cells * n_genes * 4 / 1e9 + 2.0)


# Step 04 — neighbors + UMAP + grid Leiden + stability on the integrated
#   object (X_integrated is 100-dim PCA, so the full-gene / HVG X terms
#   do not apply here).  Two-segment linear fit over measured runs:
#   <620k: 12 GB fixed + ~17.6 GB/M (110k -> 12, 620k -> 21);
#   >=620k: 21 GB + ~45 GB/M (1.05M -> 40.5).  1.676M projects to
#   ~68.5 GB and ran fine on a 100 GB budget.  The 12 GB floor covers
#   scanpy/GPU import + h5ad load + grid/stability working set.
def _estimate_step04_peak(n_cells: int) -> float:
    """Estimated peak RSS (GB) for step 04 (cluster + UMAP).

    X_integrated is 100-dim, so the estimate depends on cell count only.
    Two-segment linear fit over 110k/620k/1.05M measured runs; 1.676M
    projected 68.5 GB and completed on a 100 GB WSL budget."""
    if n_cells < 620_000:
        return max(12.0, 12.0 + 17.6 * (n_cells - 110_000) / 1e6)
    return 21.0 + 45.0 * (n_cells - 620_000) / 1e6


# Step 05 — KB annotation on the annotated object.
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


# Step 06 — subcluster on one cell type.  Dual-variable: the parent
# 05_annotated full read materialises raw (full-gene dense lognorm, ~4B/gene)
# and the subset work materialises sub_raw from 02_qc + PCA dense:
#   parent_read = parent_cells * parent_genes * 4B        (transient peak)
#   sub_dense   = n_cells * n_genes * 4B                 (to_memory + PCA)
#   sub_sparse  = nnz * 12B                              (resident CSR copy)
# Calibrated: 61k subset @ 27k genes -> ~8-9 GB peak; small subsets floor
# at ~2.5 GB.  Wall-clock: t = 5 + 0.30ms*parent_cells + 2.4ms*n_cells
# (GSE118546 61k: 5 + 18 + 146 = 169s vs 156.7s measured).
def _estimate_step06_wall(n_cells: int, parent_cells: int = 0) -> float:
    """Estimated wall-clock (s) for step 06 (GPU subcluster path).

    Fixed ~5s startup + full parent read (~0.30ms/cell of parent 05, sc.read
    materialising raw) + subset compute (neighbors/UMAP/leiden/plot ~2.4ms/cell)
    + h5ad write (~0.1ms/cell).  Calibrated on GSE118546 61k (156.7s) and
    four other post-optimisation runs (2.5k-52k cells, +-20%)."""
    return 5.0 + 0.30e-3 * parent_cells + 2.4e-3 * n_cells + 0.1e-3 * n_cells


def _estimate_step06_peak(
    n_cells: int,
    n_genes: int,
    nnz: int,
    parent_cells: int = 0,
    parent_genes: int = 0,
) -> float:
    """Estimated peak RSS (GB) for step 06 (subcluster).

    ``n_cells/n_genes/nnz`` describe the subset; ``parent_cells/parent_genes``
    the parent 05_annotated (defaults to the subset when omitted)."""
    parent_read_gb = parent_cells * (parent_genes or n_genes) * 4 / 1e9
    sub_dense_gb = n_cells * n_genes * 4 / 1e9
    sub_sparse_gb = nnz * 12 / 1e9
    return max(2.5, max(parent_read_gb, sub_dense_gb + sub_sparse_gb) + 1.5)


# Step 07 — DE on the annotated object.  Two regimes:
#
#   use_raw=True (default):  the whole raw gene set is DE'd; peak is dominated
#     by the in-memory raw CSR + the wilcoxon transpose (xt) working set, i.e.
#     ~ file_size x 6-8 (zstd decompress + sparse->dense + xt).  Anchors:
#       71.5k x 4002 (nnzRaw 158M):  12.1 GB   (wall 85.7s)
#       136k  x 4002 (nnzRaw 229M):  16.3 GB   (wall 158.6s)
#       234k  x 33770 (nnzRaw 658M): 25.0 GB   (wall ~300s, re-run estimate)
#
#   use_raw=False:  backed load + raw dropped -> only the HVG X dense slice
#     stays resident; peak ~ 1.5 x X dense + base.  Anchors:
#       1.05M x 4000: 13.3 GB (wall 182.3s)
#       1.68M x 4000: 21.4 GB (wall 295.3s)
#
#   Wall-clock: two regimes split at the 50k fast-path threshold
#     (numba patches active above it):
#       slow (<50k cells):    wall = 5.0e-4 * (n*g)^0.666
#       fast (>=50k cells):   wall = 1.2 * (n*g)^0.229 * (1.0 if raw else 0.55)


def _estimate_step07_peak(
    n_cells: int, n_genes: int, nnz: int, approximation: str = "exact"
) -> float:
    """Estimated peak RSS (GB) for step 07 (multi-layer DE).

    ``n_genes`` is the DE gene count (HVG ~4000 when ``use_raw=False``, else
    the raw gene set ~34k); ``nnz`` should be the *raw* matrix nnz when
    available (CSR 12B/nnz on disk, ~6-8x decompressed in memory).
    """
    if approximation == "fast":
        # use_raw=False / backed path: HVG dense X only (n x 4000 x 4B)
        x_mib = n_cells * 4000 * 4 / 2**20
        return max(2.0, x_mib * 0.9 / 1e3 + 1.7)
    # exact: raw CSR resident + wilcoxon xt transpose.
    # Calibrated on 9 measured step-07 runs; the per-nnz cost shrinks with
    # matrix size (transpose chunking amortises), so a plain linear fit
    # over-predicts at >600M nnz.  Two-segment model: A≈5.5 below 60M nnz,
    # A≈3.2 above (GSE116106/241268/249004 anchors; Zuo2024 658M).
    raw_gb = nnz * 12 / 1e9
    a_coef = 5.5 if nnz < 60_000_000 else 3.2
    x_gb = n_cells * 4000 * 4 / 1e9
    return max(3.0, raw_gb * a_coef + x_gb + 1.7)


def _estimate_step07_wall(n_cells: int, n_genes: int, use_raw: bool = True) -> float:
    """Estimated wall-clock (s) for step 07; split at the 50k fast-path gate."""
    work = n_cells * n_genes
    if n_cells < 50_000:
        return 5.0e-4 * work**0.666
    return 1.2 * work**0.229 * (1.0 if use_raw else 0.55)


# Step 08 — trajectory (PAGA/DPT) on the annotated object.  Peak is
#   dominated by the 05_final raw (full-gene *dense* lognorm, 4B/gene —
#   unlike step 07's sparse raw), read once and held with the HVG X +
#   dense diffmap (n x 15 x 8B) + per-branch X[mask] DE copies.  Anchors
#   (post-optimisation runs):
#     2.3k  x 58.8k (GSE107618):  1.9 GB   (wall 11.8s)
#     17.8k x 32.5k (GSE202212):  2.7 GB   (wall 51.6s)
#     33.9k x 29.0k (GSE122680):  4.9 GB   (wall 153.0s)
#     39.3k x 33.7k (GSE116106):  7.9 GB   (wall 257.9s)
#     71.5k x 32.3k (GSE241268): 18.9 GB   (wall 138.4s — kendall patched,
#                                        few branches; upper bound 3.5x here)
def _estimate_step08_peak(n_cells: int, n_genes: int, nnz: int, parent_genes: int = 0) -> float:
    """Estimated peak RSS (GB) for step 08 (trajectory).

    ``n_genes`` is the *full* raw gene count (the dense raw dominates);
    ``nnz`` is only used as a floor via the X sparse term.  Calibrated so
    every anchor is an upper bound: raw_dense x 1.7 fits GSE241268
    (18.9 GB) and over-covers the four smaller anchors."""
    g_full = parent_genes or n_genes
    raw_dense = n_cells * g_full * 4 / 1e9
    x_hvg = n_cells * n_genes * 4 * 0.10 / 1e9
    diffmap = n_cells * 15 * 8 / 1e9
    return max(2.5, 1.5 + raw_dense * 1.7 + x_hvg + diffmap + 1.0)


def _estimate_step08_wall(n_cells: int, n_genes: int, n_branches: int = 10) -> float:
    """Estimated wall-clock (s) for step 08 — conservative upper bound.

    Linear load/DPT term + O(n x g) vectorised Spearman (HVG scale,
    capped at 4000 genes) + per-branch DE.  Not strictly monotonic in n:
    GSE116106 (39.3k, 9 stages) runs 257.9s while GSE241268 (71.5k, 2
    groups) runs 138.4s because branch count and KB trend-gene count
    dominate."""
    g_eff = min(n_genes, 4000)
    return 5.0 + 1.5e-3 * n_cells + 1.3e-6 * n_cells * g_eff + 1.2e-4 * n_branches * g_eff


# Step 09 — enrichment.  Reads 05_annotated.h5ad (zstd) for the marker-
#   validation PASS-rate check: raw CSR resident + HVG X sparse + buffers.
#   Order-of-magnitude anchors (not precisely RSS-calibrated yet):
#     GSE107618 2.3k x 58.8k (67M h5ad):  ~4 GB
#     GSE241268 71.5k x 32.3k (1.4G h5ad): ~8 GB
def _estimate_step09_peak(n_cells: int, n_genes: int, nnz: int) -> float:
    """Estimated peak RSS (GB) for step 09 (enrichment).

    Peak is the 05_annotated read: full-gene raw CSR (12B/nnz) + HVG X
    sparse slice + python/pandas buffers.  Formula is an upper bound; no
    precise RSS calibration yet."""
    raw_csr = nnz * 12 / 1e9
    x_hvg = n_cells * 4000 * 4 * 0.10 / 1e9
    return max(4.0, raw_csr + x_hvg + 2.0)


# Step 10 — exploratory plots on the annotated object.  Peak is the
#   full 05_annotated read (raw CSR 12B/nnz + HVG X CSR resident) plus
#   the light plot AnnData (stratified downsample) and matplotlib
#   transient; the plot term is decoupled from n once sampling kicks in.
#   Anchors (post-OOM-fix runs, plot_max_cells=20000):
#     GSE234963 166.8k x 4.0k (raw nnz 352M):  7.1 GB   (wall 47.3s)
#     GSE116106 39.3k x ~4.0k:                  2.6 GB   (wall 22.3s)
def _estimate_step10_peak(
    n_cells: int, n_genes: int, nnz: int, plot_max_cells: int = 20_000
) -> float:
    """Estimated peak RSS (GB) for step 10 (exploratory).

    ``nnz`` is the *raw* matrix nnz (12B/nnz CSR, resident for the whole
    step); ``n_genes`` is the HVG count for the X sparse slice.  The plot
    term covers plot_adata (X slice + marker-only raw) and matplotlib
    transient — constant once n exceeds plot_max_cells."""
    raw_csr = nnz * 12 / 1e9
    x_hvg = n_cells * n_genes * 0.10 * 12 / 1e9
    obs = 0.2 * n_cells / 1e6
    plot = 0.6 * min(1.0, n_cells / plot_max_cells)
    return max(2.0, 1.5 + raw_csr + x_hvg + obs + plot)


def _estimate_step10_wall(file_bytes: int, plot_max_cells: int = 20_000) -> float:
    """Estimated wall-clock (s) for step 10 — IO-dominated load.

    The full 05_annotated read at ~60-90 MB/s (9p + gzip) plus a ~8s
    constant for plotting (sampling keeps it n-invariant).  Calibrated:
    GSE234963 3.56 GB -> 55s (measured 47.3s); GSE116106 0.9 GB -> 20s
    (measured 22.3s)."""
    return file_bytes / 75e6 + 8.0


# Step 11 — GRN.  Streaming pseudobulk keeps peak constant regardless of
#   scale: one 20k-cell chunk (6B/nnz) + sums (n_groups x n_genes
#   float64) + pseudo DataFrame; ULM/CollecTRI are tiny matrices.
#   Anchors (2026-08-06, streaming rewrite): GSE138002 110k and
#   GSE137398 76k — both <1.5 GB peak.
def _estimate_step11_peak() -> float:
    return 1.2


def _estimate_step11_wall(nnz: int) -> float:
    """Estimated wall-clock (s) for step 11 — IO-dominated load.

    Fixed ~3.5s (import + CollecTRI parse + ULM/export/plot) plus the raw
    read at ~10s per 1e9 nnz (h5py, page-cache dependent).  Calibrated:
    GSE138002 250M nnz -> 3.78s; GSE137398 175M nnz -> 5.76s (warm cache).
    """
    return 3.5 + nnz / 1e9 * 10.0


# Step 12 — cell-cell interaction (LIANA).  Reads 05_annotated in full:
#   raw CSR (12B/nnz, int64 indices) stays resident for the whole step;
#   the LR-filtered slice (~10% of nnz, int32) plus its normcounts copy
#   (logfc always triggers) add a second term; the permutation cube is a
#   sparse matvec X.T @ S (MB-scale, ignored).  Calibrated on the 1.05M-
#   cell Li2026_Lobe_Neurons run (peak 51.9 GB measured under a 60 GB
#   RLIMIT; raw CSR 27.1 GB at 2.257e9 nnz):
#     peak = 1.6 * (nnz*12e-9 + 2 * 0.1*nnz*8e-9 + 3)
def _estimate_step12_peak(nnz: int) -> float:
    """Estimated peak RSS (GB) for step 12 (cell-cell interaction).

    The full-gene raw CSR (12B/nnz) dominates; LIANA's filtered slice
    (~10% of nnz) plus the normcounts layer copy are the second term.
    ``nnz`` should be the *raw* matrix nnz.  Calibrated: Lobe_Neurons
    2.257e9 nnz -> 53.9 GB est vs 51.9 GB measured."""
    return max(6.0, 1.6 * (nnz * 12e-9 + 2 * 0.1 * nnz * 8e-9 + 3.0))


def _estimate_step12_wall(
    nnz: int,
    n_groups: int = 10,
    permutations: int = 100,
    n_jobs: int = 24,
) -> float:
    """Estimated wall-clock (s) for step 12 — three additive phases.

    Load: h5py direct sparse read + CSR build at ~2e7 elements/s (Lobe
    2.257e9 nnz -> 113s).  Permutation cube: one sparse matvec X.T @ S
    per permutation, thread-parallel at ~3.5e8 nnz*s/s/thread.  _get_lr:
    per-group mean/zscore passes are *single-threaded* and do not shrink
    with n_jobs.  Calibrated: Lobe_Neurons 100 perms -> 267.6s measured
    (267 est); GSE241268 1000 perms -> 49.5s measured (39 est, -20%)."""
    nnz_f = 0.1 * nnz  # LR-filtered slice (~10% of raw nnz)
    load = nnz / 2e7
    perm = nnz_f * n_groups * permutations / (n_jobs * 3.5e8)
    lr = nnz_f * n_groups * 12 / 2e8
    return load + perm + lr + 8.0


# Step 00 — raw load (post-T8/T1b rework).  The load builds ONE output CSR
#   (nnz×12B: float32 data + int64 indices) plus the chunked-write staging
#   (T2 block-bounded: one block ≲ 500k rows) plus the obs DataFrame
#   (~536 B/cell) and a 1.0 GiB base.  ``concat_factor`` covers union-var
#   growth: 1.0 for single-file loads, ~1.3 for multi-file merges
#   (identical gene sets take the in-place fast path, so 1.3 is a
#   conservative bound on near-1.0 reality — metis G11).  Calibrated on
#   the 5 T8/T1b measured datasets (runner peak_rss_mib / time -v):
#     Lobe_Neurons 28×10X_h5 1.204M c / 2.801e9 nnz -> 34.11 GiB (est 43.6)
#     Multiome     10×10X_h5  70.5k c / 0.426e9 nnz ->  7.49 GiB (est  8.6)
#     GSE173180    csv_table  50.9k c / 0.108e9 nnz ->  3.42 GiB (est  3.6)
#     GSE202735    preproc    32.1k c / 0.053e9 nnz ->  2.14 GiB (est  2.9)
#     GSE239410    MTX-mmread 137.5k c / 0.156e9 nnz -> 3.76 GiB (est  4.1)
#   StressTest 83×10X_h5 1.973M c / 4.468e9 nnz -> est 68.2 GiB (<100 GB,
#   metis G8).  write_staging = 1.5 GB keeps the 5 anchors in bracket.
def _estimate_step00_peak(
    n_cells: int,
    n_genes: int,
    nnz: int,
    budget_bytes: int = 0,
    concat_factor: float = 1.0,
) -> float:
    """Estimated peak RSS (GB) for step 00 (raw load).

    ``nnz`` is the *raw input* non-zero count (summed over all input files);
    ``concat_factor`` is 1.0 for single-file loads and ~1.3 for multi-file
    merges (union-var growth).  The output CSR (12 B/nnz) is the dominant
    term; the T2 block-bounded write adds a fixed ~1.5 GB staging; the obs
    frame adds ~536 B/cell; 1.0 GiB covers scanpy/pandas/python base."""
    raw_csr = nnz * 12 * concat_factor / 1e9  # resident merged CSR
    write_staging = 1.5  # one chunked-write block + h5py buffers
    obs = n_cells * 536 / 1e9
    return raw_csr + write_staging + obs + 1.0


def estimate_step_peak(
    step: int,
    n_cells: int,
    n_genes: int,
    nnz: int = 0,
    policy: str = "speed",
    budget_bytes: int = 0,
    approximation: str = "exact",
    parent_cells: int = 0,
    parent_genes: int = 0,
    plot_max_cells: int = 20_000,
    concat_factor: float = 1.0,
) -> float:
    """Estimated peak RSS (GB) for a pipeline step.

    For step 06 the first three arguments describe the *subset* (cells /
    full gene count / subset nnz); ``parent_cells`` / ``parent_genes``
    describe the parent 05_annotated object whose full read dominates the
    loading phase."""
    if step == 0:
        return _estimate_step00_peak(n_cells, n_genes, nnz, budget_bytes, concat_factor)
    if step == 1:
        return _estimate_step01_peak(n_cells, n_genes, nnz, budget_bytes)
    if step == 2:
        return _estimate_step02_peak(n_cells, nnz, budget_bytes)
    if step == 3:
        return _estimate_step03_peak(n_cells, n_genes, nnz, policy)
    if step == 4:
        return _estimate_step04_peak(n_cells)
    if step == 5:
        return _estimate_step05_peak(n_cells, n_genes, nnz, approximation)
    if step == 6:
        return _estimate_step06_peak(
            n_cells,
            n_genes,
            nnz,
            parent_cells=parent_cells,
            parent_genes=parent_genes,
        )
    if step == 7:
        return _estimate_step07_peak(n_cells, n_genes, nnz, approximation=approximation)
    if step == 8:
        return _estimate_step08_peak(
            n_cells,
            n_genes,
            nnz,
            parent_genes=parent_genes,
        )
    if step == 9:
        return _estimate_step09_peak(n_cells, n_genes, nnz)
    if step == 10:
        return _estimate_step10_peak(
            n_cells,
            n_genes,
            nnz,
            plot_max_cells=plot_max_cells,
        )
    if step == 11:
        return _estimate_step11_peak()
    if step == 12:
        return _estimate_step12_peak(nnz)
    return 0.0


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
