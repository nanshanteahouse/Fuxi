#!/usr/bin/env python3
"""
Step 01: QC filtering + MACS3 peak calling + peak matrix + doublet detection
==============================================================================
  - Filter cells by fragment counts (Snapatac2 2.9: no TSS add_tsse avail)
  - Call peaks via MACS3 → export BED → create peak-by-cell matrix
  - Scrublet doublet detection (requires .X, i.e. peak matrix)

Input:  00_raw.h5ad (fragment-level AnnData)
Output: 01_filtered.h5ad (peak-by-cell matrix with qc flags)
"""

import sys, os, time, argparse, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, validate_adata
import snapatac2 as snap


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("01_qc", os.path.join(CFG.log_dir, "01_qc.log"))
    log.info("Step 01: QC + MACS3 + peak matrix + doublet")

    if os.path.exists(CFG.filtered_h5ad):
        log.info("Skip: %s exists.", CFG.filtered_h5ad)
        return

    # SnapATAC2 default backed mode — lazy loading, near-zero memory
    data = snap.read(CFG.raw_h5ad)
    log.info("Loaded: %d cells (backed mode)", data.n_obs)

    # ── Filter cells (counts only; TSS not available in 2.9) ──
    n0 = data.n_obs
    snap.pp.filter_cells(data, min_counts=CFG.min_fragments,
                         max_counts=CFG.max_fragments, min_tsse=None)
    log.info("Filtered: %d → %d cells (-%.1f%%)",
             n0, data.n_obs, 100 * (n0 - data.n_obs) / max(n0, 1))

    # ── MACS3 peak calling (SnapATAC2 stores result in uns) ──
    log.info("MACS3 (qval=%.2f)...", CFG.peak_qval)
    snap.tl.macs3(data, qvalue=CFG.peak_qval, n_jobs=CFG.n_jobs)
    # In SnapATAC2 2.9 backed mode, uns is PyElemCollection — use subscript, not .get()
    import polars as pl
    try:
        peaks = data.uns['macs3_pseudobulk']
    except KeyError:
        log.error("MACS3 returned no peaks (macs3_pseudobulk not in uns).")
        sys.exit(1)
    if not isinstance(peaks, pl.DataFrame) or len(peaks) == 0:
        log.error("MACS3 returned empty peaks.")
        sys.exit(1)
    # Filter to standard chromosomes only (avoid alt/haplotype contigs)
    standard_chroms = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}
    peaks = peaks.filter(pl.col('chrom').is_in(standard_chroms))
    if len(peaks) == 0:
        log.error("MACS3 returned no peaks on standard chromosomes.")
        sys.exit(1)
    log.info("  Peaks: %d (on standard chromosomes)", len(peaks))

    bed = os.path.join(tempfile.gettempdir(), "atac_peaks.bed")
    peaks.select(['chrom', 'start', 'end']).write_csv(bed, separator='\t', include_header=False)

    # ── Peak-by-cell matrix (SnapATAC2 2.9 make_peak_matrix has no n_jobs) ──
    log.info("Creating peak matrix...")
    peak_data = snap.pp.make_peak_matrix(data, peak_file=bed, backend='hdf5')
    log.info("  Matrix: %d cells × %d peaks", peak_data.n_obs, peak_data.n_vars)

    # ── Scrublet ──
    log.info("Scrublet doublet detection...")
    try:
        snap.pp.scrublet(peak_data, features=None, random_state=CFG.random_seed)
        peak_data.obs['predicted_doublet'] = peak_data.obs['doublet_probability'] > 0.5
        n_dbl = peak_data.obs['predicted_doublet'].sum()
        log.info("  Doublets: %d (%.1f%%)", int(n_dbl),
                 100 * n_dbl / max(peak_data.n_obs, 1))
    except Exception as e:
        log.warning("Scrublet failed (likely OOM), marking all cells as non-doublets: %s", e)
        peak_data.obs['predicted_doublet'] = False

    validate_adata(peak_data, stage_name="01_qc", logger=log)
    safe_write(peak_data, CFG.filtered_h5ad, cfg=CFG, compression_override=None)
    log.info("Step 01 complete (%.1fs)", time.time() - t0)


if __name__ == '__main__':
    main()
