#!/usr/bin/env python3
"""
Step 08: GO/KEGG enrichment on marker peak-associated regions
===============================================================
  - Reads marker_peaks.csv from step 05
  - Peak-to-gene mapping via nearest TSS (pybedtools)
  - Enrichr over-representation analysis per group

Input:  marker_peaks.csv (from step 05)
Output: enrichment_results.csv
"""

import sys, os, time, argparse, gc
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_plot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def peak_to_gene(peak_df, genome="hg38"):
    """Map peak coordinates to nearest gene symbols using pybedtools.
    Falls back to using peak names directly if mapping is unavailable.
    """
    try:
        import pybedtools
        # Build a BED from the first column (assumed chr:start-end format)
        if peak_df.empty:
            return []
        peaks = peak_df.iloc[:, 0].astype(str)
        bed_lines = []
        for p in peaks:
            p = p.replace(':', '\t').replace('-', '\t')
            bed_lines.append(p)
        if not bed_lines:
            return peaks.tolist()
        bed = pybedtools.BedTool('\n'.join(bed_lines), from_string=True)
        nearest = bed.closest(
            pybedtools.BedTool(os.path.join(pybedtools.helpers.get_tempdir(),
                                            f"{genome}_genes.bed")),
            d=True)
        genes = []
        for interval in nearest:
            name = str(interval).split('\t')[-2] if len(str(interval).split('\t')) > 1 else str(interval)
            if name.startswith('chr'):
                name = interval.name if hasattr(interval, 'name') else name
            genes.append(name)
        return genes
    except Exception:
        # Fallback: return peak names directly (caller should handle this)
        return peak_df.iloc[:, 0].astype(str).tolist()


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("08_enrichment", os.path.join(CFG.log_dir, "08_enrichment.log"))
    log.info("Step 08: GO/KEGG enrichment")

    marker_csv = os.path.join(CFG.table_dir, "marker_peaks.csv")
    if not os.path.exists(marker_csv):
        log.warning("marker_peaks.csv not found. Run Step 05 first.")
        return

    markers_df = pd.read_csv(marker_csv)
    log.info("Loaded: %d rows", len(markers_df))

    import gseapy as gp

    all_results = []
    group_col = 'group' if 'group' in markers_df.columns else None

    for grp in (markers_df[group_col].unique() if group_col else ['all']):
        if group_col:
            sub = markers_df[markers_df[group_col] == grp]
            genes = peak_to_gene(sub, genome=CFG.genome)
        else:
            genes = peak_to_gene(markers_df, genome=CFG.genome)

        # Deduplicate and filter — strip invalid names
        # peak_to_gene returns peak names (chr:start-end) as fallback — filter them
        genes_clean = []
        for g in genes:
            if isinstance(g, str) and g and len(g) > 2:
                if not g.startswith('chr'):
                    genes_clean.append(g)
        genes = genes_clean

        if len(genes) < 5:
            log.info("Enrichment: %s — too few genes (%d), skipping", str(grp), len(genes))
            continue
        log.info("Enrichment: %s (%d unique genes)", str(grp), len(genes))
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=CFG.enrichment_gene_sets,
                organism=CFG.enrichment_organism,
                outdir=None, no_plot=True,
            )
            enr.results['group'] = str(grp)
            all_results.append(enr.results)
        except Exception as e:
            log.debug("Enrichment failed for %s: %s", str(grp), e)

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(os.path.join(CFG.table_dir, "enrichment_results.csv"), index=False)
        log.info("Saved enrichment_results.csv (%d rows)", len(combined))

        try:
            top = combined.sort_values('Adjusted P-value').head(20)
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.barh(range(len(top)), -np.log10(top['Adjusted P-value'].values + 1e-10))
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(top['Term'].str[:50])
            ax.set_xlabel('-log10(Adjusted P-value)')
            plt.tight_layout()
            plt.savefig(os.path.join(CFG.figure_dir, "enrichment_barplot.png"), dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            log.warning("Barplot failed: %s", e)
    else:
        pd.DataFrame().to_csv(os.path.join(CFG.table_dir, "enrichment_results.csv"), index=False)

    gc.collect()
    log.info("Step 08 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
