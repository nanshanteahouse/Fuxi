#!/usr/bin/env python3
"""
Step 02: Feature selection + spectral + KNN
=============================================
  - Remove doublets (predicted_doublet column)
  - Select top features (IDF-weighted)
  - Matrix-free spectral embedding (SnapATAC2 Lanczos algorithm)
  - KNN graph construction

Input:  01_filtered.h5ad
Output: 02_processed.h5ad
"""

import sys, os, time, argparse, gc
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, validate_adata
import snapatac2 as snap
import numpy as np
import scipy.sparse as sp


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("02_process", os.path.join(CFG.log_dir, "02_process.log"))
    log.info("Step 02: Feature selection + spectral + KNN")

    if os.path.exists(CFG.processed_h5ad):
        log.info("Skip: %s exists.", CFG.processed_h5ad)
        return

    # Load to memory (SnapATAC2 backed mode does not support subscript/copy)
    data = snap.read(CFG.filtered_h5ad, backed=None)
    log.info("Loaded: %d cells, %d peaks (in-memory)", data.n_obs, data.n_vars)

    # ── Remove predicted doublets ──
    pred_dbl = data.obs['predicted_doublet']
    d = int(pred_dbl.sum())
    if d > 0:
        keep = ~pred_dbl.values.astype(bool)
        data = data[keep].copy()
        gc.collect()
        log.info("Removed %d doublets → %d cells", d, data.n_obs)

    # ── Feature selection (out-of-core, works on backed data) ──
    snap.pp.select_features(data, n_features=CFG.n_features)

    # ── Ensure float64 for SnapATAC2 spectral (Rust backend requires it) ──
    if sp.issparse(data.X) and data.X.dtype != np.float64:
        data.X = data.X.astype(np.float64, copy=False)
        log.info("X converted to float64 for spectral embedding")

    # ── Spectral embedding (matrix-free Lanczos) ──
    # Use sample_size for large datasets to enable Nyström approximation
    spectral_kwargs = dict(
        n_comps=CFG.n_spectral,
        random_state=CFG.random_seed,
    )
    sample_size = getattr(CFG, 'spectral_sample_size', None)
    if sample_size and data.n_obs > sample_size:
        spectral_kwargs['sample_size'] = sample_size
        log.info("Spectral with Nyström (sample_size=%s)", sample_size)
    snap.tl.spectral(data, **spectral_kwargs)

    # ── KNN graph ──
    snap.pp.knn(data, n_neighbors=CFG.n_neighbors)

    validate_adata(data, stage_name="02_process", logger=log)
    safe_write(data, CFG.processed_h5ad, cfg=CFG, compression_override=None)
    gc.collect()
    log.info("Step 02 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
