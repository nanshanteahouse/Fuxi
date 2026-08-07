#!/usr/bin/env python3
"""
Step 00: Load ATAC fragments → AnnData
=========================================
  - fragments.tsv.gz → Snapatac2 AnnData (via import_fragments)
  - Auto-detect chrom_sizes from fragment file (with hg38 reference fallback)
  - Optional max_cells downsampling (view-based, no copy)
  - Streaming write via file= parameter to avoid full RAM load

Output: 00_raw.h5ad
"""

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import snapatac2 as snap

from core.atac_utils.chrom_sizes import HG38_CHROM_SIZES, auto_chrom_sizes  # noqa: E402
from core.utils import resolve_config, safe_write, setup_logger, validate_adata

_N_STANDARD_CHROMS = len(HG38_CHROM_SIZES)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(cfg.log_dir, "00_load.log"))
    log.info("Step 00: Load raw ATAC data | Format: %s", cfg.data_format)

    if os.path.exists(cfg.raw_h5ad):
        log.info("Skip: %s exists.", cfg.raw_h5ad)
        return

    if cfg.data_format == "10x_fragments":
        frag_path = os.path.abspath(cfg.data_input.fragment_file)
        if not os.path.exists(frag_path):
            log.error("Fragment file not found: %s", frag_path)
            sys.exit(1)

        chrom_sizes = cfg.atac.chrom_sizes
        if isinstance(chrom_sizes, str):
            chrom_sizes = (
                None
                if not os.path.isfile(chrom_sizes)
                else {k: int(v) for k, v in (line.strip().split() for line in open(chrom_sizes))}
            )
        if not chrom_sizes:
            log.info("Auto-detecting chrom_sizes...")
            chrom_sizes = auto_chrom_sizes(frag_path)
            log.info("  Found %d chromosomes", len(chrom_sizes))

        whitelist = None
        if cfg.data_input.barcodes_file and os.path.exists(cfg.data_input.barcodes_file):
            with open(cfg.data_input.barcodes_file) as f:
                whitelist = [ln.strip() for ln in f if ln.strip()]
            log.info("Loaded %d whitelist barcodes", len(whitelist))

        # Determine sorted_by_barcode from config (default True = faster)
        sorted_by_bc = getattr(cfg.data_input, "sorted_by_barcode", True)

        log.info("Importing fragments (sorted_by_barcode=%s)...", sorted_by_bc)
        data = snap.pp.import_fragments(
            fragment_file=Path(frag_path),
            chrom_sizes=chrom_sizes,
            whitelist=whitelist,
            sorted_by_barcode=sorted_by_bc,
            min_num_fragments=0,
            n_jobs=cfg.execution.n_jobs,
            file=Path(cfg.raw_h5ad),
        )
        log.info("Imported: %d cells → %s", data.n_obs, cfg.raw_h5ad)
        log.info("Written directly to checkpoint by SnapATAC2.")
        return

    elif cfg.data_format == "h5ad":
        h5ad_path = os.path.abspath(cfg.data_input.input_h5ad)
        if not os.path.exists(h5ad_path):
            log.error("h5ad not found: %s", h5ad_path)
            sys.exit(1)
        backed = getattr(cfg.data_input, "backed", "") or None
        if backed:
            data = snap.read(h5ad_path, backed=backed)
        else:
            import scanpy as sc

            data = sc.read(h5ad_path)
        log.info("Loaded h5ad: %d cells", data.n_obs)

    elif cfg.data_format == "10x_peak_h5":
        h5_path = os.path.abspath(cfg.data_input.input_h5ad)
        if not os.path.exists(h5_path):
            log.error("10x peak h5 not found: %s", h5_path)
            sys.exit(1)
        try:
            import scanpy as sc

            data = sc.read_10x_h5(h5_path)
            # 10x peak h5 stores features in 'gene_ids' column — rename for clarity
            if "gene_ids" in data.var.columns:
                data.var_names = data.var["gene_ids"].astype(str)
            # Ensure feature names are peak-like (chr:start-end)
            # 10x ATAC uses 'interval' or 'feature_type' column if available
            if "feature_types" in data.var.columns and "Peaks" in data.var["feature_types"].values:
                log.info(
                    "Detected 10x ATAC peak matrix: %d cells, %d peaks", data.n_obs, data.n_vars
                )
            else:
                log.info("Loaded 10x h5: %d cells, %d features", data.n_obs, data.n_vars)
        except Exception as e:
            log.error(
                "Failed to read 10x peak h5 (%s). "
                "Try converting to h5ad format first with: "
                "import scanpy as sc; sc.read_10x_h5(path).write('out.h5ad')",
                e,
            )
            sys.exit(1)
    else:
        log.error("Unknown data_format: %s", cfg.data_format)
        sys.exit(1)

    data.uns["config"] = {"genome": cfg.atac.genome, "data_format": cfg.data_format}

    # ── Downsampling (view-based, no .copy()) ──
    max_cells = getattr(cfg.atac, "max_cells", None)
    if max_cells and data.n_obs > max_cells:
        rng = np.random.RandomState(cfg.execution.random_seed)
        idx = rng.choice(data.n_obs, size=max_cells, replace=False)
        idx.sort()
        data = data[idx]  # view-based indexing — no full memory copy
        log.info("Downsampled to %d cells (view)", data.n_obs)

    # ── Unified sparse format & precision (skip if X is None, as in raw fragments) ──
    if data.X is not None:
        if getattr(cfg.execution, "force_csr", True):
            import scipy.sparse as sp

            if hasattr(data, "X") and sp.issparse(data.X):
                if not sp.isspmatrix_csr(data.X):
                    data.X = data.X.tocsr()
                    log.info("X format converted to CSR")
        if getattr(cfg.execution, "use_float32", False):
            import scipy.sparse as sp

            if sp.issparse(data.X):
                data.X = data.X.astype("float32", copy=False)
            else:
                data.X = data.X.astype("float32")
            log.info("X precision converted to float32")
    else:
        log.info("X is None (raw fragment data) — skipping format conversion")

    if data.X is not None:
        validate_adata(data, stage_name="00_load", logger=log)
    safe_write(data, cfg.raw_h5ad, cfg=cfg)
    log.info("Step 00 complete (%.1fs)", time.time() - t0)


if __name__ == "__main__":
    main()
