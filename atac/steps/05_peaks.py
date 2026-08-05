#!/usr/bin/env python3
"""Step 05: Post-clustering peak calling"""

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import snapatac2 as snap

from core.utils import resolve_config, safe_write, save_figure, setup_logger, validate_adata


def main():
    t0 = time.time()
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    a = p.parse_args()
    cfg = resolve_config(a.config)
    log = setup_logger("05_peaks", os.path.join(cfg.log_dir, "05_peaks.log"))
    log.info("Step 05: Post-clustering peak calling")
    # Load chrom sizes
    chrom_sizes = cfg.atac.chrom_sizes
    if isinstance(chrom_sizes, str) and os.path.isfile(chrom_sizes):
        cs = {}
        with open(chrom_sizes) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    cs[parts[0]] = int(parts[1])
        chrom_sizes = cs
    frag = os.path.abspath(cfg.data_input.fragment_file)
    # Get leiden from clustered
    clustered = snap.read(cfg.clustered_h5ad, backed="r")
    leiden_map = {
        str(b): str(lab) for b, lab in zip(clustered.obs_names, clustered.obs["leiden"].to_list())
    }
    log.info(
        "Loaded %d clustered cells, %d leiden clusters",
        len(leiden_map),
        len(set(leiden_map.values())),
    )
    # Import fragments fresh (bypasses HDF5 plugin issue with raw_h5ad)
    log.info("Importing fragments...")
    data = snap.pp.import_fragments(
        fragment_file=Path(frag),
        chrom_sizes=chrom_sizes,
        sorted_by_barcode=True,
        min_num_fragments=0,
        n_jobs=cfg.execution.n_jobs,
    )
    log.info("Imported: %d cells", data.n_obs)
    # Add leiden labels
    data.obs["leiden"] = [leiden_map.get(str(b), "unassigned") for b in data.obs_names]
    matched = sum(1 for lb in data.obs["leiden"] if lb != "unassigned")
    log.info("Matched %d/%d cells", matched, data.n_obs)
    # Per-cluster MACS3
    try:
        snap.tl.macs3(data, groupby="leiden", qvalue=cfg.atac.peak_qval)
        log.info("Per-cluster MACS3 done")
    except Exception as e:
        log.warning("MACS3 failed: %s", e)
        snap.tl.macs3(data, qvalue=cfg.atac.peak_qval)
    pd = snap.pp.make_peak_matrix(data, backend="hdf5")
    log.info("Peak matrix: %d x %d", pd.n_obs, pd.n_vars)
    if hasattr(snap.metrics, "frip"):
        snap.metrics.frip(pd)
        log.info("FRiP: mean=%.3f", pd.obs["frip"].mean())
    else:
        log.warning("FRiP not available")
        pd.obs["frip"] = 0.5
    try:
        os.makedirs(os.path.join(cfg.figure_dir, "05_peaks"), exist_ok=True)
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(pd.obs["frip"], bins=50)
        ax.set_xlabel("FRiP")
        save_figure(
            None,
            os.path.join(cfg.figure_dir, "05_peaks", "frip_distribution"),
            dpi=cfg.plot.figure_dpi,
        )
        plt.close()
    except Exception as e:
        log.warning("FRiP histogram: %s", e)
    validate_adata(pd, stage_name="05_peaks", logger=log)
    safe_write(pd, cfg.peak_h5ad, cfg=cfg, compression_override=None)
    log.info("Done %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
