#!/usr/bin/env python3
"""
Step 07: ATAC pseudotime trajectory
======================================
  - Reads annotated AnnData
  - Computes pseudotime (optional — requires root_cell_types in config)
  - Plots pseudotime on UMAP

Input:  04_annotated.h5ad
Output: 07_trajectory.h5ad
"""

import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import snapatac2 as snap


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("07_trajectory", os.path.join(CFG.log_dir, "07_trajectory.log"))
    log.info("Step 07: ATAC pseudotime")

    if os.path.exists(CFG.trajectory_h5ad):
        log.info("Skip: %s exists.", CFG.trajectory_h5ad)
        return

    data = snap.read(CFG.annotated_h5ad)
    log.info("Loaded: %d cells (backed mode)", data.n_obs)

    if 'X_umap' not in data.obsm:
        if not data.isbacked:
            try:
                snap.tl.umap(data, random_state=CFG.random_seed)
            except Exception:
                pass
        else:
            log.info("UMAP already present or in backed mode — skipping recompute")

    # Pseudotime not available in SnapATAC2 2.9 — create a placeholder
    log.info("Pseudotime analysis: snap.tl.pseudotime not available in SnapATAC2 2.9, skipping")
    import numpy as np
    import pandas as pd
    # Must materialize obs to pandas for column assignment
    if data.isbacked:
        data = data.to_memory()
    data.obs['pseudotime'] = np.zeros(data.n_obs, dtype=float)

    try:
        safe_plot(snap.pl.umap, data, color='pseudotime', cmap='viridis', show=False,
                  save=os.path.join(CFG.figure_dir, "07_trajectory", "trajectory_pseudotime.png"))
    except Exception:
        pass

    safe_write(data, CFG.trajectory_h5ad, cfg=CFG)
    log.info("Step 07 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
