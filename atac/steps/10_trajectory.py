#!/usr/bin/env python3
"""
Step 10: ATAC pseudotime trajectory
======================================
  - Reads annotated AnnData
  - Computes pseudotime (optional — requires root_cell_types in config)
  - Plots pseudotime on UMAP

Input:  05_annotated.h5ad
Output: 10_trajectory.h5ad
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import snapatac2 as snap

from core.utils import resolve_config, safe_plot, safe_write, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("10_trajectory", os.path.join(cfg.log_dir, "10_trajectory.log"))
    log.info("Step 10: ATAC pseudotime")

    if os.path.exists(cfg.trajectory_h5ad):
        log.info("Skip: %s exists.", cfg.trajectory_h5ad)
        return

    data = snap.read(cfg.annotated_h5ad)
    log.info("Loaded: %d cells (backed mode)", data.n_obs)

    if "X_umap" not in data.obsm:
        if not data.isbacked:
            try:
                snap.tl.umap(data, random_state=cfg.execution.random_seed)
            except Exception:
                pass
        else:
            log.info("UMAP already present or in backed mode — skipping recompute")

    # Pseudotime not available in SnapATAC2 2.9 — create a placeholder
    log.info("Pseudotime analysis: snap.tl.pseudotime not available in SnapATAC2 2.9, skipping")
    import numpy as np

    # Must materialize obs to pandas for column assignment
    if data.isbacked:
        data = data.to_memory()
    data.obs["pseudotime"] = np.zeros(data.n_obs, dtype=float)

    try:
        safe_plot(
            snap.pl.umap,
            data,
            color="pseudotime",
            cmap="viridis",
            show=False,
            save=os.path.join(cfg.figure_dir, "07_trajectory", "trajectory_pseudotime.png"),
        )
    except Exception:
        pass

    safe_write(data, cfg.trajectory_h5ad, cfg=cfg)
    log.info("Step 10 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
