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

import sys, os, time, argparse, gc, bisect
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_plot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed


def peak_to_gene(peak_df, genome="hg38", gene_bed=None, max_distance=100000):
    """Map peak coordinates to nearest gene symbols.

    Tries in order:
      1. pybedtools with genome-specific gene BED
      2. Pure-Python nearest-TSS using a configured gene annotation BED
      3. Returns empty list on failure (caller skips enrichment)

    Parameters:
        peak_df: DataFrame whose first column contains peak names (chr:start-end)
        genome: reference genome (used for pybedtools gene BED lookup)
        gene_bed: path to gene TSS annotation BED (chr, start, end, gene_name, strand)
        max_distance: maximum TSS distance (bp) to associate a peak with a gene

    Returns:
        list of gene symbols (one per peak, empty string for unmapped)
    """
    if peak_df.empty:
        return []

    peaks = peak_df.iloc[:, 0].astype(str)

    # ── Strategy 1: pybedtools ────────────────────────────────────────
    try:
        import pybedtools
        bed_lines = []
        for p in peaks:
            p = p.replace(':', '\t').replace('-', '\t')
            bed_lines.append(p)
        if not bed_lines:
            return peaks.tolist()
        bed = pybedtools.BedTool('\n'.join(bed_lines), from_string=True)
        gb = os.path.join(pybedtools.helpers.get_tempdir(), f"{genome}_genes.bed")
        nearest = bed.closest(pybedtools.BedTool(gb), d=True)
        genes = []
        for interval in nearest:
            parts = str(interval).split('\t')
            name = parts[-2] if len(parts) > 4 else str(interval)
            if name.startswith('chr'):
                name = interval.name if hasattr(interval, 'name') else name
            if name and not name.startswith('chr'):
                genes.append(name)
        if genes and len(genes) == len(peaks):
            return genes
    except Exception:
        pass

    # ── Strategy 2: Pure-Python nearest-TSS ───────────────────────────
    if gene_bed and os.path.exists(gene_bed):
        try:
            # Parse gene BED: chrom, start, end, gene_name, [score, strand]
            gene_coords = {}
            with open(gene_bed) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    chrom = parts[0]
                    # TSS: start for + strand, end for - strand
                    strand = parts[5] if len(parts) > 5 else '+'
                    tss = int(parts[1]) if strand != '-' else int(parts[2])
                    name = parts[3]
                    gene_coords.setdefault(chrom, []).append((tss, name))

            # Sort each chromosome's genes by TSS for bisect
            for chrom in gene_coords:
                gene_coords[chrom].sort(key=lambda x: x[0])

            genes = []
            for p in peaks:
                try:
                    chrom, rest = p.split(':', 1) if ':' in p else (p, '')
                    mid = int(rest.split('-')[0]) if rest else 0
                except (ValueError, IndexError):
                    genes.append('')
                    continue

                if chrom not in gene_coords:
                    genes.append('')
                    continue

                tss_list = gene_coords[chrom]
                tss_positions = [t[0] for t in tss_list]
                idx = bisect.bisect_left(tss_positions, mid)

                best = (float('inf'), '')
                if idx < len(tss_list):
                    d = abs(tss_list[idx][0] - mid)
                    if d < best[0]:
                        best = (d, tss_list[idx][1])
                if idx > 0:
                    d = abs(tss_list[idx - 1][0] - mid)
                    if d < best[0]:
                        best = (d, tss_list[idx - 1][1])

                genes.append(best[1] if best[0] <= max_distance else '')

            return genes
        except Exception:
            pass

    # ── Strategy 3: Cannot map ────────────────────────────────────────
    return []


def _enrichr_one_group(grp, genes, CFG, log):
    """Run Enrichr ORA for a single group (used by ThreadPoolExecutor)."""
    import gseapy as gp
    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=CFG.enrichment_gene_sets,
            organism=CFG.enrichment_organism,
            outdir=None, no_plot=True,
        )
        enr.results['group'] = str(grp)
        return (grp, enr.results)
    except Exception as e:
        log.debug("Enrichment failed for %s: %s", str(grp), e)
        return (grp, None)


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

    all_results = []
    group_col = 'group' if 'group' in markers_df.columns else None

    # ── Build task list: gene list per group ──
    tasks = []
    for grp in (markers_df[group_col].unique() if group_col else ['all']):
        if group_col:
            sub = markers_df[markers_df[group_col] == grp]
            genes = peak_to_gene(sub, genome=CFG.genome,
                                 gene_bed=getattr(CFG, 'gene_annotation_bed', ''),
                                 max_distance=CFG.peak_gene_distance)
        else:
            genes = peak_to_gene(markers_df, genome=CFG.genome,
                                 gene_bed=getattr(CFG, 'gene_annotation_bed', ''),
                                 max_distance=CFG.peak_gene_distance)

        # Deduplicate and filter — strip empty/invalid names
        genes_clean = [g for g in genes if isinstance(g, str) and g and len(g) > 2]

        if len(genes_clean) < 5:
            if len(genes_clean) == 0 and genes:
                log.warning("Enrichment: %s — no genes mapped (check pybedtools or gene_annotation_bed config)", str(grp))
            else:
                log.info("Enrichment: %s — too few genes (%d), skipping", str(grp), len(genes_clean))
            continue
        genes = genes_clean
        log.info("Enrichment: %s (%d unique genes)", str(grp), len(genes))
        tasks.append((grp, genes))

    # ── Parallel Enrichr calls (like RNA-09) ──
    if tasks:
        max_workers = min(5, getattr(CFG, 'n_jobs', 4))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_grp = {
                executor.submit(_enrichr_one_group, grp, genes, CFG, log): grp
                for grp, genes in tasks
            }
            for future in as_completed(future_to_grp):
                grp, res = future.result()
                if res is not None:
                    all_results.append(res)

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
