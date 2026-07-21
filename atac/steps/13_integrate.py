#!/usr/bin/env python3
"""
Step 13: RNA+ATAC integration via muon (optional)
===================================================
  - Reads ATAC annotated AnnData + RNA AnnData
  - Matches common cells
  - Builds MuData object with muon

Input:  05_annotated.h5ad + cfg.rna_h5ad
Output: 13_integrated.h5ad
"""

import argparse
import gc
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import scanpy as sc

from core.utils import resolve_config, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("13_integrate", os.path.join(cfg.log_dir, "13_integrate.log"))
    log.info("Step 13: RNA+ATAC integration")

    if os.path.exists(cfg.integrated_h5ad):
        log.info("Skip: %s exists.", cfg.integrated_h5ad)
        return

    if not cfg.has_rna_data():
        # Last-resort auto-discovery (run_pipeline.py may already have set it)
        from core.utils import find_rna_h5ad

        auto_rna = find_rna_h5ad(cfg=cfg, log=log)
        if auto_rna:
            cfg.rna_h5ad = auto_rna
            log.info("Auto-discovered RNA h5ad: %s", auto_rna)
        else:
            log.info("No RNA data. Skipping.")
            return

    try:
        import snapatac2 as snap

        atac = snap.read(cfg.annotated_h5ad)
    except Exception:
        atac = sc.read(cfg.annotated_h5ad)
    log.info("ATAC: %d cells", atac.n_obs)

    rna = sc.read(cfg.rna_h5ad)
    log.info("RNA: %d cells", rna.n_obs)

    common = rna.obs_names.intersection(atac.obs_names)
    log.info("Common cells: %d", len(common))

    if len(common) == 0:
        log.warning("No common cells between RNA and ATAC.")
        return

    try:
        from mudata import MuData

        mdata = MuData({"rna": rna[common].copy(), "atac": atac[common].copy()})
        log.info("MuData: %s", str(mdata))

        import muon as mu

        mu.pp.pca(mdata, n_comps=min(30, mdata.n_obs - 1))

        mdata.write(cfg.integrated_h5ad)
        log.info("Saved integrated MuData.")
    except Exception as e:
        log.warning("Integration failed: %s", e)
    finally:
        gc.collect()

    log.info("Step 13 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
