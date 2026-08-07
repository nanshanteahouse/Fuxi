#!/usr/bin/env python3
"""
Step 01a: Scrublet 双细胞检测 (per sample, joblib 并行)
=========================================================
从 02_qc.py 中独立出的 Scrublet 步骤，不含 QC 指标计算或过滤。

输入: 00_raw.h5ad
输出: 01_doublet.h5ad (含 doublet_scores / predicted_doublet 列)
"""

import argparse
import math
import multiprocessing
import os
import sys
import tempfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import scanpy as sc
from joblib import Parallel, delayed

from core.utils import resolve_config, resolve_memory_settings, setup_logger

# ── Multiprocess ANN kNN patch (Scrublet approx mode) ──
# Annoy's get_nns_by_item holds the GIL, so threads cannot parallelize it.
# We monkey-patch scrublet's get_knn_graph: the main process builds the index
# once, saves it to a temp file, and fork workers each load a copy to query
# a row block. Results are identical because index + query params are
# deterministic. Parallelism only kicks in for large manifolds (>50k rows)
# and only outside loky workers (no nested fork).

_MP_KNN_NJOBS = 1
_MP_KNN_ORIGINAL = None
_MP_MEM_BUDGET_BYTES = 0  # 0 = no budget (unlimited)

_ZSCORE_CHUNK_ROWS = 20000
# Output dtype of the chunked sparse_zscore. float64 = bit-exact parity with
# scrublet's original implementation; float32 halves zscore peak memory (~20
# GiB per 157k-cell group) at the cost of ~+9% doublet-label drift. Set via
# scrublet.zscore_float32.
_ZSCORE_DTYPE = np.float64


def _sparse_zscore_chunked(e, gene_mean=None, gene_stdev=None):
    # Numerically equivalent to scrublet's sparse_zscore, which materializes
    # (e - gene_mean) as a dense matrix (~20-23 GiB transient per large
    # group). Computing row-chunked keeps the transient allocation at chunk
    # size while producing the same dense result.
    from scrublet.helper_functions import sparse_var

    if gene_mean is None:
        gene_mean = e.mean(0)
    if gene_stdev is None:
        gene_stdev = np.sqrt(sparse_var(e))
    gm = np.asarray(gene_mean).ravel()
    gs = np.asarray(gene_stdev).ravel()
    nrow, ncol = e.shape
    out = np.empty((nrow, ncol), dtype=_ZSCORE_DTYPE)
    for i in range(0, nrow, _ZSCORE_CHUNK_ROWS):
        j = min(i + _ZSCORE_CHUNK_ROWS, nrow)
        out[i:j] = (e[i:j].toarray() - gm) / gs
    return out


def _greedy_buckets(peaks, n_buckets):
    """Bin-pack group peaks into n_buckets with a largest-first greedy.",
    Returns list of bucket -> list of group indices (bucket-internal order
    is irrelevant for peak estimation: groups run serially inside a bucket)."""
    order = sorted(range(len(peaks)), key=lambda i: peaks[i], reverse=True)
    buckets = [[] for _ in range(n_buckets)]
    acc = [0.0] * n_buckets
    for gi in order:
        b = min(range(n_buckets), key=lambda j: acc[j])
        buckets[b].append(gi)
        acc[b] += peaks[gi]
    return buckets


def _run_bucket(raw_path, small_names, small_idxs, bucket, cfg):
    """Run one greedy bucket serially (module-level: pickled by loky).
    Workers open the backed h5ad themselves — the main process's AnnData
    holds an h5py handle and cannot be pickled."""
    import anndata as ad

    adata = ad.read_h5ad(raw_path, backed="r")
    out = []
    for gi in bucket:
        scores, pred = run_scrublet_sample(
            _extract_subset(adata, small_idxs[gi]), small_names[gi], cfg
        )
        out.append((gi, scores, pred))
    return out


def _run_small_parallel(raw_path, small_names, small_idxs, buckets, cfg):
    """Joblib-backed small-group runner (executed in a background thread so the
    main process can process large groups serially at the same time).

    Each bucket is one job: groups inside a bucket run serially (peak memory
    adds up), buckets run in parallel. Buckets are produced by _greedy_buckets
    so the per-worker peak tracks the *sum* of the bucket's groups instead of
    the largest single group."""
    return Parallel(n_jobs=len(buckets), initializer=_install_zscore_patch)(
        delayed(_run_bucket)(raw_path, small_names, small_idxs, b, cfg) for b in buckets if b
    )


def _set_zscore_dtype(cfg):
    """Pick the chunked zscore output dtype from scrublet.zscore_float32."""
    global _ZSCORE_DTYPE
    _ZSCORE_DTYPE = np.float32 if cfg.scrublet.zscore_float32 else np.float64


def _install_zscore_patch():
    """Replace scrublet's dense-materializing sparse_zscore with the chunked
    version. Also passed as the joblib Parallel initializer so loky worker
    processes (which re-import scrublet fresh) get the patched version."""
    import scrublet.helper_functions as _hf
    import scrublet.scrublet as _scr

    _hf.sparse_zscore = _sparse_zscore_chunked
    _scr.sparse_zscore = _sparse_zscore_chunked


def _annoy_query_block(path, npc, metric, k, i0, i1):
    from annoy import AnnoyIndex

    idx = AnnoyIndex(npc, metric=metric)
    idx.load(path)
    knn = np.zeros((i1 - i0, k), dtype=np.int32)
    for i in range(i0, i1):
        knn[i - i0] = idx.get_nns_by_item(i, k + 1)[1:]
    return knn


def _get_knn_graph_mp(
    x, k=5, dist_metric="euclidean", approx=False, return_edges=True, random_seed=0
):
    if not approx:
        return _MP_KNN_ORIGINAL(x, k, dist_metric, approx, return_edges, random_seed)
    if dist_metric == "cosine":
        dist_metric = "angular"
    n, npc = x.shape
    if (
        _MP_KNN_NJOBS <= 1
        or n <= 50000
        or multiprocessing.current_process().name.startswith("SpawnProcess")
    ):
        return _MP_KNN_ORIGINAL(x, k, dist_metric, approx, return_edges, random_seed)

    n_jobs = min(_MP_KNN_NJOBS, max(1, math.ceil(n / 25000)))
    if _MP_MEM_BUDGET_BYTES > 0:
        # ~25 B/point/tree with 10 trees; + kNN result matrix + 2 GiB headroom
        index_bytes = n * 10 * 25
        knn_bytes = n * (k + 1) * 4
        headroom = 2 * 2**30
        mem_cap = max(1, (_MP_MEM_BUDGET_BYTES - knn_bytes - headroom) // max(index_bytes, 1))
        if mem_cap < n_jobs:
            warnings.warn(
                f"kNN parallelism capped by memory budget: {n_jobs} -> {mem_cap} workers "
                f"(index ~{index_bytes / 2**20:.0f} MiB x workers, knn matrix ~{knn_bytes / 2**20:.0f} MiB)"
            )
        n_jobs = min(n_jobs, mem_cap)
    if n_jobs <= 1:
        return _MP_KNN_ORIGINAL(x, k, dist_metric, approx, return_edges, random_seed)

    fd, tmp_path = tempfile.mkstemp(suffix=".annoy")
    os.close(fd)
    try:
        from annoy import AnnoyIndex

        idx = AnnoyIndex(npc, metric=dist_metric)
        idx.set_seed(random_seed)
        for i in range(n):
            idx.add_item(i, list(x[i]))
        idx.build(10)
        idx.save(tmp_path)
        del idx

        n_jobs = min(_MP_KNN_NJOBS, max(1, math.ceil(n / 25000)))
        edges = np.linspace(0, n, n_jobs + 1).astype(int)
        blocks = [
            (tmp_path, npc, dist_metric, k, int(edges[i]), int(edges[i + 1]))
            for i in range(n_jobs)
        ]
        ctx = multiprocessing.get_context("fork")
        with ctx.Pool(n_jobs) as pool:
            parts = pool.starmap(_annoy_query_block, blocks)
        knn = np.vstack(parts).astype(int)
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    if return_edges:
        links = set()
        for i in range(knn.shape[0]):
            for j in knn[i]:
                links.add(tuple(sorted((i, j))))
        return links, knn
    return knn


def _install_knn_mp_patch(cfg):
    """Replace scrublet's approx kNN query with the multiprocess version."""
    global _MP_KNN_NJOBS, _MP_KNN_ORIGINAL, _MP_MEM_BUDGET_BYTES
    import scrublet.helper_functions as _hf
    import scrublet.scrublet as _scr

    if _MP_KNN_ORIGINAL is None:
        _MP_KNN_ORIGINAL = _hf.get_knn_graph
    _MP_KNN_NJOBS = min(cfg.execution.n_jobs or os.cpu_count() or 1, 16)
    _MP_MEM_BUDGET_BYTES = resolve_memory_settings(cfg)[1]
    _hf.get_knn_graph = _get_knn_graph_mp
    _scr.get_knn_graph = _get_knn_graph_mp


def _resolve_doublet_rate(cfg, n_cells: int) -> float:
    """返回用于 Scrublet 的 expected_doublet_rate。

    如果 config 显式设置了值，直接使用；否则按 10X 官方拟合：
      y (%) = 0.000759 * x + 0.053    (x = recovered cell 数)
    钳位在 [0.004, 0.15]。
    """
    if cfg.scrublet.expected_doublet_rate is not None:
        return cfg.scrublet.expected_doublet_rate
    rate = 0.00000759 * n_cells + 0.00000053
    return min(max(rate, 0.004), 0.15)


def _stable_threshold(scrub, expected_rate):
    """Stable, machine-independent fallback threshold: the quantile that
    labels ~expected_rate of observed cells (unlike scrublet's fixed-score
    fallback which can flip entire groups when the score distribution sits
    around the constant)."""
    obs = np.asarray(scrub.doublet_scores_obs_).ravel()
    return float(np.percentile(obs, 100 * (1 - expected_rate)))


def _threshold_is_stable(scrub, rel_tol=0.05, max_abs_change=0.05):
    """Scrublet's auto-threshold (skimage threshold_minimum on the sim
    histogram) can sit in the densest part of the observed score
    distribution when there is no clear bimodality. There, tiny input
    perturbations (BLAS non-determinism) move the threshold and flip whole
    groups between runs. Check the robustness of the *observed* crossing
    fraction against a ±rel_tol threshold shift: if the fraction changes by
    more than max_abs_change (absolute), the auto-threshold is not
    trustworthy and callers should use _stable_threshold instead. Groups
    with <1% or >99% crossing are exempt: the threshold sits on a sparse
    tail, so wobble cannot flip a meaningful number of cells."""
    obs = np.asarray(scrub.doublet_scores_obs_).ravel()
    thr = float(scrub.threshold_)
    base = np.mean(obs > thr)
    if base <= 0.01 or base >= 0.99:
        return True  # sparse tail: threshold wobble is harmless
    for f in (1 - rel_tol, 1 + rel_tol):
        frac = np.mean(obs > thr * f)
        if abs(frac - base) > max_abs_change:
            return False
    return True


def run_scrublet_sample(adata_sub, sample_name, cfg):
    try:
        import scipy.sparse as sp
        import scrublet as scr

        expected_rate = _resolve_doublet_rate(cfg, adata_sub.n_obs)
        x_mat = adata_sub.X
        if not isinstance(x_mat, sp.spmatrix):
            if hasattr(x_mat, "to_memory"):
                # Backed AnnData (e.g. "all" branch): _CSRDataset is not a scipy
                # matrix and np.asarray() on it yields an object array, which
                # crashes sp.csr_matrix() → Scrublet silently returned all zeros.
                x_mat = x_mat.to_memory()
            else:
                x_mat = sp.csr_matrix(x_mat)
        scrub = scr.Scrublet(
            x_mat,
            expected_doublet_rate=expected_rate,
            random_state=cfg.execution.random_seed,
        )
        scores, predicted = scrub.scrub_doublets(
            min_counts=cfg.scrublet.min_counts,
            min_cells=cfg.scrublet.min_cells,
            min_gene_variability_pctl=cfg.scrublet.min_gene_var_pctl,
            n_prin_comps=cfg.scrublet.n_prin_comps,
            svd_solver=cfg.scrublet.svd_solver,
        )
        if predicted is None:
            fallback = _stable_threshold(scrub, expected_rate)
            warnings.warn(
                f"Scrublet auto-threshold failed for {sample_name}, "
                f"falling back to stable threshold={fallback:.4f}"
            )
            predicted = scrub.call_doublets(threshold=fallback)
            if predicted is None:
                warnings.warn(
                    f"Scrublet threshold fallback failed for {sample_name}, assuming no doublets"
                )
                predicted = np.zeros(adata_sub.n_obs, dtype=bool)
        elif not _threshold_is_stable(scrub):
            fallback = _stable_threshold(scrub, expected_rate)
            warnings.warn(
                f"Auto-threshold unstable for {sample_name} (thr={scrub.threshold_:.4f}), "
                f"using stable fallback={fallback:.4f}"
            )
            predicted = scrub.call_doublets(threshold=fallback)
            if predicted is None:
                predicted = np.zeros(adata_sub.n_obs, dtype=bool)
        return scores, predicted
    except Exception as e:
        warnings.warn(f"Scrublet failed for {sample_name}: {e}")
        return np.zeros(adata_sub.n_obs), np.zeros(adata_sub.n_obs, dtype=bool)


def _extract_subset(adata, idx):
    """Extract a fully-materialized subset from a backed AnnData.

    Backed AnnData views (produced by ``adata[idx].to_memory()``) retain
    ``ElementRef``/``SparseCSRMatrixView`` references to the HDF5 backing
    file, making them unpicklable for joblib parallel dispatch.  This
    constructs a clean AnnData from raw slices, breaking all back references.
    """
    import scipy.sparse as sp

    sub = adata[idx]
    x_mat = sub.X
    if sp.issparse(x_mat):
        x_mat = sp.csr_matrix(x_mat)  # materialize view to concrete CSR
    else:
        x_mat = np.array(x_mat)
    return sc.AnnData(x_mat, obs=sub.obs.copy(), var=sub.var.copy())


def detect_doublets_parallel(adata, cfg, log):
    global _MP_KNN_NJOBS
    if not cfg.scrublet.run:
        adata.obs["doublet_scores"] = 0.0
        adata.obs["predicted_doublet"] = False
        return adata.obs["doublet_scores"].values, adata.obs["predicted_doublet"].values

    if cfg.expression_type != "raw_counts":
        _policy = cfg.scrublet.on_non_counts
        if _policy == "abort":
            log.error(
                "Scrublet requires raw UMI counts but expression_type='%s'. "
                "Set scrublet.on_non_counts to 'skip_warn' or 'skip_silent' to proceed, "
                "or use raw UMI counts as input.",
                cfg.expression_type,
            )
            sys.exit(1)
        elif _policy == "skip_warn":
            log.warning(
                "Scrublet is designed for raw UMI counts. "
                "expression_type='%s' violates the negative-binomial assumption. "
                "Disabling Scrublet.",
                cfg.expression_type,
            )
        # skip_silent: no log message
        adata.obs["doublet_scores"] = 0.0
        adata.obs["predicted_doublet"] = False
        return adata.obs["doublet_scores"].values, adata.obs["predicted_doublet"].values

    log.info("Running Scrublet (per sample, parallel)...")
    configured_key = cfg.scrublet.batch_key
    if configured_key in adata.obs:
        groupby_col = configured_key
        log.info("  Using configured batch column: %s", groupby_col)
    else:
        log.warning(
            "Configured batch column '%s' not in adata.obs, falling back to 'sample'/'stage'",
            configured_key,
        )
        groupby_col = "sample" if "sample" in adata.obs else "stage"
    if groupby_col not in adata.obs:
        log.warning("Group column (%s) not found, running Scrublet on all data.", groupby_col)
        scores, pred = run_scrublet_sample(adata, "all", cfg)
        adata.obs["doublet_scores"] = scores
        adata.obs["predicted_doublet"] = pred
        log.info(
            "  Predicted doublets: %d / %d (%.1f%%)", pred.sum(), adata.n_obs, 100 * pred.mean()
        )
        return scores, pred

    sample_groups = adata.obs.groupby(groupby_col, observed=True)

    # Memory-aware scheduling: large samples (>serial_threshold cells) serially,
    # small samples (<=serial_threshold) in parallel. For backed AnnData, subsets are
    # extracted via .to_memory() to avoid pulling the full sparse matrix.
    memory_threshold = getattr(cfg.scrublet, "serial_threshold", 15000)
    large_names, large_idxs = [], []
    small_names, small_idxs = [], []
    for name, idx in sample_groups.indices.items():
        if len(idx) > memory_threshold:
            large_names.append(name)
            large_idxs.append(idx)
        else:
            small_names.append(name)
            small_idxs.append(idx)

    results = [None] * (len(large_names) + len(small_names))
    n_cpu = cfg.execution.n_jobs or os.cpu_count() or 1
    budget = _MP_MEM_BUDGET_BYTES or resolve_memory_settings(cfg)[1]

    def _extract_bytes(n_cells):
        # In-memory float32 AnnData for one group (happens in the main process).
        return int(n_cells * adata.n_vars * 4)

    def _run_peak(n_cells):
        # Per-worker residency: normalized sparse matrix (~10% density) + zscore
        # result (~10% of genes survive the variability filter) at _ZSCORE_DTYPE,
        # plus a fixed 0.8 GiB base.
        sparse = n_cells * adata.n_vars * 4 * 0.10
        zscore = n_cells * adata.n_vars * 0.1 * np.dtype(_ZSCORE_DTYPE).itemsize
        return int(sparse + zscore + int(0.8 * 2**30))

    small_n_jobs = min(n_cpu - 1, len(small_names)) if small_names else 0
    avail = None
    if small_n_jobs > 0 and budget > 0:
        max_large = max(large_idxs, key=len, default=[])
        max_small = max(small_idxs, key=len, default=[])
        # Main process: one large extraction + one small extraction happen
        # concurrently (large in this thread, small in the background thread),
        # plus the large group's run-time peak (zscore + sparse, no base).
        main_peak = (
            _extract_bytes(len(max_large))
            + _extract_bytes(len(max_small))
            + (_run_peak(len(max_large)) - int(0.8 * 2**30))
            + int(0.8 * 2**30)
        )
        worker_peak = _run_peak(len(max_small))
        avail = max(0, int(budget * 0.95) - main_peak)
        mem_cap = max(1, avail // worker_peak)
        if mem_cap < small_n_jobs:
            log.warning(
                "  Parallel sample workers capped by memory budget: %d -> %d "
                "(~%.1f GiB/worker, main-process peak ~%.1f GiB). "
                "Adjust scrublet.serial_threshold or execution.memory.budget if needed.",
                small_n_jobs,
                mem_cap,
                worker_peak / 2**30,
                main_peak / 2**30,
            )
        # Reserve cores for the large groups' kNN fork pool (manifold ≈ 3×
        # cells: n obs + 2n simulated doublets, ≥25k rows per worker, cap 16
        # workers). This is machine-independent: derived from the largest
        # large group and the CPU count, not hard-coded.
        largest_cells = max((len(i) for i in large_idxs), default=0)
        knn_reserve = max(1, min(16, math.ceil(largest_cells * 3 / 25000)))
        small_n_jobs = min(small_n_jobs, mem_cap, max(1, n_cpu - knn_reserve))
    buckets = None
    if small_names:
        # Greedy bin-packing: bucket the small groups by estimated peak so the
        # per-worker budget tracks the bucket sum, then shrink the bucket count
        # until every bucket fits the available budget.
        peaks = [_run_peak(len(i)) for i in small_idxs]
        while True:
            buckets = _greedy_buckets(peaks, small_n_jobs)
            if avail is None or small_n_jobs <= 1:
                break
            max_bucket = max((sum(peaks[gi] for gi in b) for b in buckets if b), default=0)
            if max_bucket <= avail:
                break
            small_n_jobs -= 1
            log.warning(
                "  Re-bucketing: largest bucket (~%.1f GiB) exceeds budget, "
                "reducing workers %d -> %d",
                max_bucket / 2**30,
                small_n_jobs + 1,
                small_n_jobs,
            )

    small_future = None
    if small_names:
        # Overlap: small groups run in the background (loky via a helper thread)
        # while large groups are processed serially in this (main) process, so
        # the single-threaded PCA/zscore phases of large groups no longer idle
        # the remaining cores.
        small_future = ThreadPoolExecutor(max_workers=1).submit(
            _run_small_parallel, cfg.raw_h5ad, small_names, small_idxs, buckets, cfg
        )
    if large_names:
        log.info(
            "  Large groups (%s) — processing serially, small groups in parallel (n_jobs=%d)",
            ", ".join(f"{n}({len(i)} cells)" for n, i in zip(large_names, large_idxs)),
            small_n_jobs,
        )
        # Hand remaining cores to the small-group workers; the kNN fork pool is
        # capped below the CPU count while the background pool is active.
        _MP_KNN_NJOBS = max(1, min(n_cpu - small_n_jobs, 16))
        for i, (name, idx) in enumerate(zip(large_names, large_idxs)):
            sub = _extract_subset(adata, idx)
            results[i] = run_scrublet_sample(sub, name, cfg)
            del sub  # release early
    if small_future is not None:
        for bucket_res in small_future.result():
            for gi, scores, pred in bucket_res:
                results[len(large_names) + gi] = (scores, pred)
    all_scores = np.zeros(adata.n_obs)
    all_pred = np.zeros(adata.n_obs, dtype=bool)
    all_names = large_names + small_names
    all_idxs = list(large_idxs) + list(small_idxs)
    for (scores, pred), name, idx in zip(results, all_names, all_idxs):
        all_scores[idx] = scores
        all_pred[idx] = pred
        used_rate = _resolve_doublet_rate(cfg, len(idx))
        log.info(
            "  Sample %s: %d / %d doublets (%.1f%%) [expected_rate=%.4f]",
            name,
            pred.sum(),
            len(idx),
            100 * pred.mean(),
            used_rate,
        )

    adata.obs["doublet_scores"] = all_scores
    adata.obs["predicted_doublet"] = all_pred
    log.info(
        "  Total predicted doublets: %d / %d (%.1f%%)",
        all_pred.sum(),
        adata.n_obs,
        100 * all_pred.mean(),
    )
    return all_scores, all_pred


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("01_doublet", os.path.join(cfg.log_dir, "01_doublet.log"))
    log.info("Step 01a: Scrublet doublet detection")
    _install_knn_mp_patch(cfg)
    _set_zscore_dtype(cfg)
    _install_zscore_patch()

    # ── Pre-read memory guard: use load_meta persisted by step 00 ──
    # Step 00 writes real nnz/n_cells into perf_report.json's pipeline layer.
    # Read it BEFORE opening the h5ad (zero-copy) so step 01/02/03 peaks are
    # checked against the true matrix shape, not an estimate.
    import json as _json

    _pre_est = None
    _perf_report_path = os.path.join(cfg.results_dir, "perf_report.json")
    if os.path.exists(_perf_report_path):
        try:
            with open(_perf_report_path) as _prf:
                _pr_data = _json.load(_prf)
            _lm = (_pr_data.get("pipeline", {}) or {}).get("load_meta")
            if _lm and _lm.get("nnz"):
                from core.utils import (
                    check_memory_guard,
                    estimate_step_peak,
                    resolve_memory_settings,
                )

                _mpolicy, _mbudget, _mguard = resolve_memory_settings(cfg)
                _pre_est = {
                    s: estimate_step_peak(
                        s,
                        int(_lm["n_cells"]),
                        int(_lm["n_genes"]),
                        int(_lm["nnz"]),
                        policy=_mpolicy,
                        budget_bytes=_mbudget,
                    )
                    for s in (1, 2, 3)
                }
                if _mbudget > 0:
                    log.info(
                        "[memory-guard] pre-read peaks (from load_meta): "
                        + ", ".join(f"step{s} ~{g:.0f}GB" for s, g in _pre_est.items())
                    )
                check_memory_guard(_pre_est, _mbudget, _mguard, logger_obj=log)
        except Exception:
            log.warning("Failed to read load_meta from %s", _perf_report_path, exc_info=True)

    # Use backed mode — only load one sample group into memory at a time.
    # The raw_h5ad for large datasets (1M+ cells) occupies ~56 GiB uncompressed;
    # backed='r' keeps it on disk and reads only the requested indices.
    adata = sc.read(cfg.raw_h5ad, backed="r")
    log.info(
        "Opened in backed mode: %s — %d cells × %d genes", cfg.raw_h5ad, adata.n_obs, adata.n_vars
    )

    # ── Memory guard: estimate step 01/02/03 peaks vs budget before heavy work ──
    from core.utils import check_memory_guard, estimate_step_peak, resolve_memory_settings

    _mem_policy, _mem_budget, _mem_guard = resolve_memory_settings(cfg)
    # backed _CSRDataset has no .nnz — read from h5ad X/data directly (zero-copy)
    _nnz = 0
    try:
        import h5py

        with h5py.File(cfg.raw_h5ad, "r") as _h5:
            _nnz = int(_h5["X/data"].shape[0])
    except Exception:
        _nnz = getattr(adata.X, "nnz", 0)
    _est = {
        s: estimate_step_peak(
            s, adata.n_obs, adata.n_vars, _nnz, policy=_mem_policy, budget_bytes=_mem_budget
        )
        for s in (1, 2, 3)
    }
    if _mem_budget > 0:
        log.info(
            "[memory-guard] estimated peaks: "
            + ", ".join(f"step{s} ~{g:.0f}GB" for s, g in _est.items())
        )
    check_memory_guard(_est, _mem_budget, _mem_guard, logger_obj=log)

    doublet_scores, doublet_pred = detect_doublets_parallel(adata, cfg, log)

    out_path = os.path.join(cfg.h5ad_dir, "01_doublet.h5ad")
    # Lightweight write: copy raw h5ad + append only the 2 new obs columns,
    # avoiding a full ~9 GB AnnData rewrite for just metadata.
    import shutil

    shutil.copy2(cfg.raw_h5ad, out_path)
    log.info("Copied raw h5ad → %s (%.1f GiB)", out_path, os.path.getsize(out_path) / 2**30)
    import pandas as pd

    from core.utils import write_obs_columns_lightweight

    obs_to_write = pd.DataFrame(
        {"doublet_scores": doublet_scores, "predicted_doublet": doublet_pred},
        index=adata.obs_names,
    )
    write_obs_columns_lightweight(out_path, obs_to_write, logger=log)
    log.info("Step 01a complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
