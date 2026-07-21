#!/usr/bin/env python3
"""Step 01: Scrublet doublet detection (ATAC)"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import snapatac2 as snap

from core.utils import resolve_config, safe_write, setup_logger, validate_adata


def main():
    t0 = time.time()
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    a = p.parse_args()
    cfg = resolve_config(a.config)
    log = setup_logger("01_doublet", os.path.join(cfg.log_dir, "01_doublet.log"))
    log.info("Step 01: Scrublet doublet detection")
    data = snap.read(cfg.raw_h5ad, backed="r")
    log.info("Loaded: %d cells", data.n_obs)
    try:
        if data.var is not None and data.n_vars > 0:
            snap.pp.scrublet(data, random_state=cfg.execution.random_seed)
            n = int(data.obs["predicted_doublet"].sum())
            log.info("Doublets: %d / %d (%.1f%%)", n, data.n_obs, 100 * n / max(data.n_obs, 1))
        else:
            log.info("Raw fragment data (no peak matrix) — deferring scrublet to QC step")
            import numpy as np

            data.obs["predicted_doublet"] = np.full(data.n_obs, False)
            data.obs["doublet_scores"] = np.full(data.n_obs, 0.0)
    except Exception as e:
        log.warning("Scrublet failed: %s", e)
        import numpy as np

        data.obs["predicted_doublet"] = np.full(data.n_obs, False)
        data.obs["doublet_scores"] = np.full(data.n_obs, 0.0)
    validate_adata(data, stage_name="01_doublet", logger=log)
    safe_write(data, cfg.doublet_h5ad, cfg=cfg)
    log.info("Done %.1fs", time.time() - t0)
