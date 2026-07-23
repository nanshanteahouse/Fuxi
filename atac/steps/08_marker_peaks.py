#!/usr/bin/env python3
"""
Step 08: Marker peaks — differential accessibility
=====================================================
  - Computes marker peaks per cell_type (or leiden cluster)
  - Summarizes result as CSV

Input:  05_annotated.h5ad
Output: marker_peaks.csv
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import pandas as pd
import snapatac2 as snap

from core.utils import resolve_config, safe_plot, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("08_marker_peaks", os.path.join(cfg.log_dir, "08_marker_peaks.log"))
    log.info("Step 08: Marker peaks")

    data = snap.read(cfg.annotated_h5ad)
    log.info("Loaded: %d cells, %d peaks (backed mode)", data.n_obs, data.n_vars)

    groupby = "cell_type" if "cell_type" in data.obs else None
    if groupby is None:
        for c in data.obs.columns:
            if c.startswith("leiden"):
                groupby = c
                break
    if groupby is None:
        log.error("No clustering column found.")
        sys.exit(1)

    # SnapATAC2 2.9 marker_regions returns Dict[str, Index], not stored in uns
    markers = snap.tl.marker_regions(data, groupby=groupby, pvalue=cfg.atac.marker_peaks_fdr)
    # markers is Dict[str, pd.Index] — group -> peak names
    rows = []
    for grp, peaks in markers.items():
        for i, p in enumerate(peaks):
            rows.append({"group": grp, "peak": p, "rank": i + 1})
    if rows:
        markers_df = pd.DataFrame(rows)
        markers_df.to_csv(os.path.join(cfg.table_dir, "marker_peaks.csv"), index=False)
        log.info("Saved marker_peaks.csv (%d rows)", len(markers_df))

    try:
        safe_plot(
            snap.pl.heatmap,
            data,
            groupby=groupby,
            show=False,
            save=os.path.join(cfg.figure_dir, "05_marker_peaks", "marker_peaks_heatmap.png"),
            cfg=cfg,
        )
    except Exception as e:
        log.warning("Heatmap failed: %s", e)

    log.info("Step 08 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
