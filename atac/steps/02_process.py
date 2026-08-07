#!/usr/bin/env python3
"""
Step 02: Feature selection + spectral + KNN
======================================================================================
  - Remove doublets (predicted_doublet column)
  - Select top features (IDF-weighted)
  - Matrix-free spectral embedding (SnapATAC2 Lanczos algorithm)
  - Harmony batch correction (optional)
  - KNN graph construction

Input:  01_filtered.h5ad
Output: 02_processed.h5ad
"""

import argparse
import gc
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import scipy.sparse as sp
import snapatac2 as snap

from core.utils import (
    check_memory_guard,
    estimate_step_peak,
    monitor_performance,
    resolve_config,
    resolve_memory_settings,
    safe_write,
    setup_logger,
    validate_adata,
)


def _run_step(cfg, log, t0):
    """Core step 02 body — extracted for the perf wrapper."""
    # ── Memory guard: estimate step-02 peak (in-memory float64 X) before the
    #    heavy load, mirroring RNA 00_load.  Pre-read uses the configured
    #    n_features as the gene dimension; a second estimate with the real
    #    shape runs once the file is open.  Estimation failure never blocks.
    try:
        _policy, _budget, _guard = resolve_memory_settings(cfg)
        _n_peaks = getattr(cfg.atac, "n_features", 50_000)
        _est = {
            2: estimate_step_peak(
                2, 0, _n_peaks, modality="atac", policy=_policy, budget_bytes=_budget
            )
        }
        if _budget > 0:
            log.info("[memory-guard] estimated step 02 peak: ~%.0f GB", _est[2])
        check_memory_guard(_est, _budget, _guard, logger_obj=log)
    except Exception as e:
        log.warning("[memory-guard] step 02 estimation skipped: %s", e)

    # Load to memory (SnapATAC2 backed mode does not support subscript/copy)
    data = snap.read(cfg.filtered_h5ad, backed=None)
    log.info("Loaded: %d cells, %d peaks (in-memory)", data.n_obs, data.n_vars)

    # Re-estimate with the real cell count now that the file is open.
    try:
        _policy, _budget, _guard = resolve_memory_settings(cfg)
        _est2 = {
            2: estimate_step_peak(
                2,
                data.n_obs,
                data.n_vars,
                modality="atac",
                policy=_policy,
                budget_bytes=_budget,
            )
        }
        if _budget > 0:
            log.info("[memory-guard] estimated step 02 peak (real shape): ~%.0f GB", _est2[2])
        check_memory_guard(_est2, _budget, _guard, logger_obj=log)
    except Exception as e:
        log.warning("[memory-guard] step 02 re-estimate skipped: %s", e)

    # ── Remove predicted doublets ──
    pred_dbl = data.obs["predicted_doublet"]
    d = int(pred_dbl.sum())
    if d > 0:
        keep = ~pred_dbl.values.astype(bool)
        data = data[keep].copy()
        gc.collect()
        log.info("Removed %d doublets → %d cells", d, data.n_obs)

    # ── Feature selection (out-of-core, works on backed data) ──
    snap.pp.select_features(data, n_features=cfg.atac.n_features)

    # ── Ensure float64 for SnapATAC2 spectral (Rust backend requires it) ──
    if sp.issparse(data.X) and data.X.dtype != np.float64:
        data.X = data.X.astype(np.float64, copy=False)
        log.info("X converted to float64 for spectral embedding")

    # ── Spectral embedding (matrix-free Lanczos) ──
    # Use sample_size for large datasets to enable Nyström approximation
    spectral_kwargs = dict(
        n_comps=cfg.atac.n_spectral,
        random_state=cfg.execution.random_seed,
    )
    sample_size = cfg.atac.spectral_sample_size
    if sample_size and data.n_obs > sample_size:
        spectral_kwargs["sample_size"] = sample_size
        log.info("Spectral with Nyström (sample_size=%s)", sample_size)
    snap.tl.spectral(data, **spectral_kwargs)  # type: ignore[reportArgumentType]

    # ── KNN graph (with optional Harmony batch correction) ──
    if cfg.atac.harmony_use_harmony and cfg.atac.harmony_batch_key in data.obs:
        n_b = data.obs[cfg.atac.harmony_batch_key].nunique()
        if n_b >= 2:
            log.info("Harmony (batch=%s, %d batches)...", cfg.atac.harmony_batch_key, n_b)
            snap.pp.harmony(data, batch=cfg.atac.harmony_batch_key)
            snap.pp.knn(data, n_neighbors=cfg.clustering.n_neighbors, use_rep="X_spectral_harmony")
        else:
            log.info("Harmony skipped: only %d batch(es)", n_b)
            snap.pp.knn(data, n_neighbors=cfg.clustering.n_neighbors)
    else:
        snap.pp.knn(data, n_neighbors=cfg.clustering.n_neighbors)

    validate_adata(data, stage_name="02_process", logger=log)
    safe_write(data, cfg.processed_h5ad, cfg=cfg, compression_override=None)
    gc.collect()
    log.info("Step 02 complete, took %.1fs", time.time() - t0)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("02_process", os.path.join(cfg.log_dir, "02_process.log"))
    log.info("Step 02: Feature selection + spectral + KNN")

    if os.path.exists(cfg.processed_h5ad):
        log.info("Skip: %s exists.", cfg.processed_h5ad)
        return

    with monitor_performance("02_process", log=log):
        _run_step(cfg, log, t0)


if __name__ == "__main__":
    main()
