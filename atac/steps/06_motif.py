#!/usr/bin/env python3
"""
Step 06: Motif enrichment (Snapatac2 2.9 standalone API)
==========================================================
  - Load JASPAR motifs
  - Run motif enrichment per cluster (top significant peaks, sorted by p-value)
  - Store results as CSV

Note: chromVAR not available in Snapatac2 2.9. Motif enrichment
      uses standalone API with genome_fasta requirement.
"""

import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_plot
import numpy as np
import pandas as pd
import snapatac2 as snap


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("06_motif", os.path.join(CFG.log_dir, "06_motif.log"))
    log.info("Step 06: Motif enrichment")

    data = snap.read(CFG.annotated_h5ad)
    log.info("Loaded: %d cells, %d peaks (backed mode)", data.n_obs, data.n_vars)

    groupby = 'cell_type' if 'cell_type' in data.obs else None
    if groupby is None:
        for c in data.obs.columns:
            if c.startswith('leiden'):
                groupby = c; break

    try:
        log.info("Loading motifs (CIS-BP human)...")
        motifs = snap.datasets.cis_bp(unique=True)
        log.info("Loaded %d motifs", len(motifs))

        # Build region sets per cluster from marker_regions dict
        # marker_regions returns Dict[str, pd.Index] in SnapATAC2 2.9
        # If clusters already re-run after step 4 fix, markers may not be in data.uns
        # Instead, re-run marker_regions here if needed
        log.info("Computing marker regions for motif enrichment...")
        markers = snap.tl.marker_regions(data, groupby=groupby, pvalue=CFG.marker_peaks_fdr)
        # markers is Dict[str, pd.Index]
        regions_by_group = {}
        for grp in data.obs[groupby].unique():
            sg = str(grp)
            if sg in markers and len(markers[sg]) >= 5:
                n_top = min(500, len(markers[sg]))
                regions_by_group[sg] = markers[sg][:n_top].tolist()
                log.info("  Group %s: %d marker peaks", sg, len(regions_by_group[sg]))
            else:
                log.warning("  Group %s: insufficient marker peaks", sg)

        # Use pre-built Genome object for hg38 (auto-downloads FASTA)
        from snapatac2.genome import hg38
        genome = hg38

        if regions_by_group:
            log.info("Running motif enrichment for %d groups...", len(regions_by_group))
            try:
                result = snap.tl.motif_enrichment(
                    motifs=motifs,
                    regions=regions_by_group,
                    genome_fasta=genome,
                )
                for grp, df in result.items():
                    csv_out = os.path.join(CFG.table_dir, f"motif_enrichment_{grp}.csv")
                    df.write_csv(csv_out)
                    log.info("  Saved: %s", csv_out)
            except Exception as e:
                log.warning("Motif enrichment failed: %s", e)

    except Exception as e:
        log.warning("Motif analysis failed: %s", e)

    log.info("Step 06 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
