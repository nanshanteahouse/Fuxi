#!/usr/bin/env python3
"""Step 04: Post-clustering peak calling"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import snapatac2 as snap

from core.utils import (
    check_memory_guard,
    estimate_step_peak,
    monitor_performance,
    resolve_config,
    resolve_memory_settings,
    safe_write,
    save_figure,
    setup_logger,
    validate_adata,
)


def _run_step(cfg, log, t0):
    """Core step 04 body — extracted for the perf wrapper."""
    # Load chrom sizes (file) or auto-detect from the fragment file like
    # 00_load — mm10/other genomes get observed-max sizes, not hg38.
    chrom_sizes = cfg.atac.chrom_sizes
    if isinstance(chrom_sizes, str) and os.path.isfile(chrom_sizes):
        cs = {}
        with open(chrom_sizes) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    cs[parts[0]] = int(parts[1])
        chrom_sizes = cs
    elif not chrom_sizes:
        from core.atac_utils.chrom_sizes import auto_chrom_sizes

        chrom_sizes = auto_chrom_sizes(os.path.abspath(cfg.data_input.fragment_file))
        log.info("Auto-detected %d chromosome sizes", len(chrom_sizes))
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

    # ── Memory guard: estimate step-04 peak (fragment import + larger
    #    per-cluster peak matrix) before the heavy compute, mirroring RNA
    #    00_load.  Estimation failure never blocks.
    try:
        _policy, _budget, _guard = resolve_memory_settings(cfg)
        _n_cells = len(leiden_map)
        _est = {
            4: estimate_step_peak(
                4, _n_cells, 0, modality="atac", policy=_policy, budget_bytes=_budget
            )
        }
        if _budget > 0:
            log.info("[memory-guard] estimated step 04 peak: ~%.0f GB", _est[4])
        check_memory_guard(_est, _budget, _guard, logger_obj=log)
    except Exception as e:
        log.warning("[memory-guard] step 04 estimation skipped: %s", e)

    # Import fragments fresh (bypasses HDF5 plugin issue with raw_h5ad)
    log.info("Importing fragments...")
    data = snap.pp.import_fragments(
        fragment_file=Path(frag),
        chrom_sizes=chrom_sizes,
        sorted_by_barcode=getattr(cfg.data_input, "sorted_by_barcode", True),
        min_num_fragments=0,
        n_jobs=cfg.execution.n_jobs,
    )
    log.info("Imported: %d cells", data.n_obs)
    # Add leiden labels
    data.obs["leiden"] = [leiden_map.get(str(b), "unassigned") for b in data.obs_names]
    matched = sum(1 for lb in data.obs["leiden"] if lb != "unassigned")
    log.info("Matched %d/%d cells", matched, data.n_obs)
    # blacklist: ENCODE BED (Path object required by macs3); None disables
    blacklist = Path(cfg.atac.blacklist_bed) if cfg.atac.blacklist_bed else None
    if blacklist is not None:
        log.info("  blacklist: %s", blacklist)
    # Pseudo-replicate split for reproducible-peak calling (macs3 `replicate`).
    # Split each leiden cluster in half deterministically. Values MUST be strings
    # (macs3 casts the obs column to str internally — int values raise).
    use_rep = getattr(cfg.atac, "use_pseudo_replicates", True)
    if use_rep:
        data.obs["pseudo_replicate"] = "r0"
        for cl in sorted(set(data.obs["leiden"]) - {"unassigned"}):
            idx = [i for i, lb in enumerate(data.obs["leiden"]) if lb == cl]
            rng = np.random.RandomState(cfg.execution.random_seed)
            data.obs["pseudo_replicate"].iloc[idx] = [f"r{v}" for v in rng.randint(0, 2, len(idx))]
        log.info("  pseudo_replicate: enabled (split-half per cluster)")
    # Per-cluster MACS3
    try:
        if use_rep:
            snap.tl.macs3(
                data,
                groupby="leiden",
                qvalue=cfg.atac.peak_qval,
                replicate="pseudo_replicate",
                replicate_qvalue=cfg.atac.peak_qval,
                blacklist=blacklist,
            )
        else:
            snap.tl.macs3(
                data,
                groupby="leiden",
                qvalue=cfg.atac.peak_qval,
                blacklist=blacklist,
            )
        log.info("Per-cluster MACS3 done")
    except Exception as e:
        log.warning("MACS3 failed: %s", e)
        snap.tl.macs3(data, qvalue=cfg.atac.peak_qval, blacklist=blacklist)
    # ── Collect peaks into a BED file (per-cluster macs3 stores a dict under
    #    uns['macs3']; the pooled fallback under uns['macs3_pseudobulk']). ──
    _peaks = None
    if "macs3" in data.uns:
        _peaks = data.uns["macs3"]
    elif "macs3_pseudobulk" in data.uns:
        _peaks = data.uns["macs3_pseudobulk"]
    if _peaks is None:
        log.error("MACS3 produced no peak table (uns['macs3'/'macs3_pseudobulk'] missing).")
        sys.exit(1)
    # dict[group -> DataFrame] (per-cluster) or a single DataFrame (pooled)
    _frames = _peaks.values() if isinstance(_peaks, dict) else [_peaks]
    _dfs = [df for df in _frames if df is not None and len(df) > 0]
    if not _dfs:
        log.error("MACS3 returned empty peak tables.")
        sys.exit(1)
    import tempfile

    import polars as pl

    bed = os.path.join(tempfile.gettempdir(), "atac_percluster_peaks.bed")
    _cols = ["chrom", "start", "end"]
    if hasattr(_dfs[0], "select"):  # polars
        _all = (
            pl.concat([df.select(_cols) for df in _dfs])
            if len(_dfs) > 1
            else _dfs[0].select(_cols)
        )
        _all.write_csv(bed, separator="\t", include_header=False)
    else:  # pandas
        import pandas as _pd

        _all = _pd.concat([df[_cols] for df in _dfs]) if len(_dfs) > 1 else _dfs[0][_cols]
        _all.to_csv(bed, sep="\t", header=False, index=False)
    log.info("Peak BED written: %d rows across %d group(s)", len(_dfs[0]) * len(_dfs), len(_dfs))
    pd = snap.pp.make_peak_matrix(data, peak_file=Path(bed), backend="hdf5")
    log.info("Peak matrix: %d x %d", pd.n_obs, pd.n_vars)
    # FRiP (fraction of reads in peaks).  SnapATAC2 2.9 requires a regions
    # dict {group -> peak regions}; fall back to a count-based proxy when
    # the fragment-level API is unavailable.
    try:
        regions = {"peaks": list(pd.var_names)}
        snap.metrics.frip(pd, regions)
        log.info("FRiP: mean=%.3f", pd.obs["frip"].mean())
    except Exception as e:
        log.warning("FRiP unavailable (%s) — using count-based proxy", e)
        import scipy.sparse as _sp

        if _sp.issparse(pd.X):
            total = np.asarray(pd.X.sum(axis=1)).ravel()
            pd.obs["frip"] = total / max(total.max(), 1.0)
        else:
            pd.obs["frip"] = 0.5

    # ── Merge clustered metadata (obs + obsm) into peak matrix ──
    # 04_peaks.h5ad becomes the downstream superset consumed by 05_annotate:
    # per-cluster peak matrix + full obs metadata + spectral/UMAP embeddings
    # carried over from 03_clustered.h5ad. Cells imported by import_fragments
    # but filtered out upstream keep leiden="unassigned" and NaN embeddings.
    if clustered.isbacked:
        clustered = clustered.to_memory()
    clustered_names = [str(b) for b in clustered.obs_names]
    peak_names = [str(b) for b in pd.obs_names]
    idx_map = {name: i for i, name in enumerate(clustered_names)}
    match_idx = [idx_map.get(b, -1) for b in peak_names]
    n_matched = sum(1 for i in match_idx if i >= 0)
    n_unassigned = pd.n_obs - n_matched
    # obs: copy clustered columns not already present (leiden already set from
    # leiden_map above — never overwrite it; existing cols like frip stay put)
    for col in clustered.obs.columns:
        if col not in pd.obs.columns:
            vals = clustered.obs[col].to_numpy()
            # Preserve dtype: bool/int/float columns get a dtype-safe NA instead
            # of np.nan (which would promote to object and break h5py writes).
            if vals.dtype.kind in "biufc":
                out = np.full(len(match_idx), np.nan, dtype=float)
                out[[i for i in match_idx if i >= 0]] = [vals[i] for i in match_idx if i >= 0]
                pd.obs[col] = out
            else:
                pd.obs[col] = [vals[i] if i >= 0 else np.nan for i in match_idx]
    # obsm: spectral + UMAP embeddings, NaN rows for unmatched cells
    for emb in ("X_spectral", "X_umap"):
        if emb in clustered.obsm:
            mat = clustered.obsm[emb]
            out = np.full((pd.n_obs, mat.shape[1]), np.nan, dtype=np.float64)
            for r, i in enumerate(match_idx):
                if i >= 0:
                    out[r] = mat[i]
            pd.obsm[emb] = out
    log.info(
        "Merged clustered metadata: %d matched, %d unassigned",
        n_matched,
        n_unassigned,
    )
    try:
        os.makedirs(os.path.join(cfg.figure_dir, "04_peaks"), exist_ok=True)
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(pd.obs["frip"], bins=50)
        ax.set_xlabel("FRiP")
        save_figure(
            None,
            os.path.join(cfg.figure_dir, "04_peaks", "frip_distribution"),
            dpi=cfg.plot.figure_dpi,
        )
        plt.close()
    except Exception as e:
        log.warning("FRiP histogram: %s", e)
    validate_adata(pd, stage_name="04_peaks", logger=log)
    safe_write(pd, cfg.peak_h5ad, cfg=cfg, compression_override=None)
    log.info("Done %.1fs", time.time() - t0)


def main():
    t0 = time.time()
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    a = p.parse_args()
    cfg = resolve_config(a.config)
    log = setup_logger("04_peaks", os.path.join(cfg.log_dir, "04_peaks.log"))
    log.info("Step 04: Post-clustering peak calling")

    with monitor_performance("04_peaks", log=log):
        _run_step(cfg, log, t0)


if __name__ == "__main__":
    main()
