#!/usr/bin/env python3
"""
Step 01: QC filtering + MACS3 peak calling + peak matrix
================================
  - Filter cells by fragment counts + TSS enrichment
  - Call peaks via MACS3 → export BED → create peak-by-cell matrix
  - Scrublet doublet detection (requires .X, i.e. peak matrix)

Input:  00_raw.h5ad (fragment-level AnnData)
Output: 01_filtered.h5ad (peak-by-cell matrix with qc flags)
"""

import argparse
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import snapatac2 as snap

from core.utils import (
    check_memory_guard,
    estimate_step_peak,
    monitor_performance,
    resolve_config,
    resolve_memory_settings,
    safe_write,
    setup_logger,
    validate_adata,
)

# Standard chromosome sets per species (autosomes + X + Y)
CHROMOSOME_SETS = {
    "human": {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"},
    "mouse": {f"chr{i}" for i in range(1, 20)} | {"chrX", "chrY"},
    "rat": {f"chr{i}" for i in range(1, 21)} | {"chrX", "chrY"},
}


def _run_step(cfg, log, t0):
    """Core step 01 body — extracted for the perf wrapper."""
    # SnapATAC2 default backed mode — lazy loading, near-zero memory
    data = snap.read(cfg.raw_h5ad)
    log.info("Loaded: %d cells (backed mode)", data.n_obs)

    # ── Memory guard: estimate step-01 peak (peak-matrix build + MACS3
    #    workspace) before the heavy compute, mirroring RNA 00_load.
    try:
        _policy, _budget, _guard = resolve_memory_settings(cfg)
        _est = {
            1: estimate_step_peak(
                1, data.n_obs, 0, modality="atac", policy=_policy, budget_bytes=_budget
            )
        }
        if _budget > 0:
            log.info("[memory-guard] estimated step 01 peak: ~%.0f GB", _est[1])
        check_memory_guard(_est, _budget, _guard, logger_obj=log)
    except Exception as e:
        log.warning("[memory-guard] step 01 estimation skipped: %s", e)

    # ── TSS enrichment score
    try:
        genome_obj = getattr(snap.genome, cfg.atac.genome)  # type: ignore[reportAttributeAccessIssue]
        snap.metrics.tsse(data, genome_obj, inplace=True)
        log.info("TSS enrichment: mean=%.2f", data.obs["tsse"].mean())
        data.uns["library_tsse"] = float(data.obs["tsse"].median())
    except Exception as e:
        log.warning("TSS enrichment failed: %s - setting NaN", e)
        data.obs["tsse"] = float("nan")
        data.uns["library_tsse"] = 0.0

    # ── Fragment size distribution
    try:
        from core.utils import safe_plot

        os.makedirs(os.path.join(cfg.figure_dir, "01_qc"), exist_ok=True)
        safe_plot(
            snap.pl.frag_size_distr,
            data,
            show=False,
            save=os.path.join(cfg.figure_dir, "01_qc", "fragment_size_distribution"),
            cfg=cfg,
        )
    except Exception as e:
        log.warning("Frag size plot failed: %s", e)
    # ── Filter cells ──
    n0 = data.n_obs
    snap.pp.filter_cells(
        data,
        min_counts=cfg.atac.min_fragments,
        max_counts=cfg.atac.max_fragments,
        min_tsse=cfg.atac.min_tsse,
    )
    log.info(
        "Filtered: %d → %d cells (-%.1f%%)", n0, data.n_obs, 100 * (n0 - data.n_obs) / max(n0, 1)
    )

    # ── MACS3 peak calling (SnapATAC2 stores result in uns) ──
    blacklist = Path(cfg.atac.blacklist_bed) if cfg.atac.blacklist_bed else None
    log.info("MACS3 (qval=%.2f)...", cfg.atac.peak_qval)
    if blacklist is not None:
        log.info("  blacklist: %s", blacklist)
    snap.tl.macs3(
        data, qvalue=cfg.atac.peak_qval, n_jobs=cfg.execution.n_jobs, blacklist=blacklist
    )
    # In SnapATAC2 2.9 backed mode, uns is PyElemCollection — use subscript, not .get()
    import polars as pl

    try:
        peaks = data.uns["macs3_pseudobulk"]
    except KeyError:
        log.error("MACS3 returned no peaks (macs3_pseudobulk not in uns).")
        sys.exit(1)
    if not isinstance(peaks, pl.DataFrame) or len(peaks) == 0:
        log.error("MACS3 returned empty peaks.")
        sys.exit(1)
    # Filter to standard chromosomes only (avoid alt/haplotype contigs)
    species_lower = (cfg.species or "").strip().lower()
    standard_chroms = CHROMOSOME_SETS.get(
        species_lower,
        # Default fallback: 22 autosomes + XY (human-like)
        {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"},
    )
    peaks = peaks.filter(pl.col("chrom").is_in(standard_chroms))
    if len(peaks) == 0:
        log.error("MACS3 returned no peaks on standard chromosomes.")
        sys.exit(1)
    log.info("  Peaks: %d (on standard chromosomes)", len(peaks))

    bed = os.path.join(tempfile.gettempdir(), "atac_peaks.bed")
    peaks.select(["chrom", "start", "end"]).write_csv(bed, separator="\t", include_header=False)

    # ── Peak-by-cell matrix (SnapATAC2 2.9 make_peak_matrix has no n_jobs) ──
    log.info("Creating peak matrix...")
    peak_data = snap.pp.make_peak_matrix(data, peak_file=bed, backend="hdf5")
    log.info("  Matrix: %d cells × %d peaks", peak_data.n_obs, peak_data.n_vars)

    # ── Scrublet ──
    log.info("Scrublet doublet detection...")
    try:
        snap.pp.scrublet(peak_data, features=None, random_state=cfg.execution.random_seed)
        peak_data.obs["predicted_doublet"] = peak_data.obs["doublet_probability"] > 0.5
        n_dbl = peak_data.obs["predicted_doublet"].sum()
        log.info("  Doublets: %d (%.1f%%)", int(n_dbl), 100 * n_dbl / max(peak_data.n_obs, 1))
    except Exception as e:
        log.error("Scrublet failed: %s", e)
        log.error("Traceback:\n%s", traceback.format_exc())
        log.warning("Scrublet failed, marking all cells as non-doublets")
        peak_data.obs["predicted_doublet"] = False

    validate_adata(peak_data, stage_name="01_qc", logger=log)
    safe_write(peak_data, cfg.filtered_h5ad, cfg=cfg, compression_override=None)
    log.info("Step 01 complete (%.1fs)", time.time() - t0)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("01_qc", os.path.join(cfg.log_dir, "01_qc.log"))
    log.info("Step 01: QC + MACS3 + peak matrix + doublet")

    if os.path.exists(cfg.filtered_h5ad):
        log.info("Skip: %s exists.", cfg.filtered_h5ad)
        return

    with monitor_performance("01_qc", log=log):
        _run_step(cfg, log, t0)


if __name__ == "__main__":
    main()
