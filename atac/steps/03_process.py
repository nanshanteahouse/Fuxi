#!/usr/bin/env python3
"""
Step 03: Feature selection + spectral + KNN
======================================================================================
  - Remove doublets (predicted_doublet column)
  - Select top features (IDF-weighted)
  - Matrix-free spectral embedding (SnapATAC2 Lanczos algorithm)
  - Harmony batch correction (optional)
  - KNN graph construction

Input:  01_doublet.h5ad
Output: 03_processed.h5ad
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

from core.utils import resolve_config, safe_write, setup_logger, validate_adata


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("03_process", os.path.join(cfg.log_dir, "03_process.log"))
    log.info("Step 03: Feature selection + spectral + KNN")

    if os.path.exists(cfg.processed_h5ad):
        log.info("Skip: %s exists.", cfg.processed_h5ad)
        return

    # Load to memory (SnapATAC2 backed mode does not support subscript/copy)
    data = snap.read(cfg.filtered_h5ad, backed=None)
    log.info("Loaded: %d cells, %d peaks (in-memory)", data.n_obs, data.n_vars)

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
    sample_size = getattr(cfg.atac, "spectral_sample_size", None)
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

    validate_adata(data, stage_name="03_process", logger=log)
    safe_write(data, cfg.processed_h5ad, cfg=cfg, compression_override=None)
    gc.collect()
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
