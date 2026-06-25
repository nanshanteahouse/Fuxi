#!/usr/bin/env python3
"""
Step 00: Load raw spatial transcriptomics data
================================================
Supports:
  1. 10X Visium (SpaceRanger output) — sq.read.visium()
  2. Generic h5ad with spatial coordinates in obsm['spatial']

Input:  Raw data directory or .h5ad file
Output: 00_raw.h5ad (with spatial coordinates + image in uns)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import scanpy as sc
import numpy as np
import scipy.sparse as sp


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(CFG.log_dir, "00_load.log"))
    log.info("Step 00: Load raw spatial transcriptomics data")
    log.info("Format: %s, Platform: %s", CFG.data_format, CFG.spatial_platform)

    if os.path.exists(CFG.raw_h5ad):
        log.info("Skip: %s already exists. Delete it to force reload.", CFG.raw_h5ad)
        return

    # ── Load by data_format ──────────────────────────────────────────────
    if CFG.data_format == "visium":
        import squidpy as sq

        # If library_id is not set, auto-detect the first Visium directory
        if not CFG.library_id:
            candidates = [
                d for d in sorted(os.listdir(CFG.data_dir))
                if os.path.isdir(os.path.join(CFG.data_dir, d))
                and os.path.exists(os.path.join(CFG.data_dir, d, 'filtered_feature_bc_matrix.h5'))
            ]
            if not candidates:
                # Try the data_dir itself as a visium directory
                if os.path.exists(os.path.join(CFG.data_dir, 'filtered_feature_bc_matrix.h5')):
                    candidates = ['']
            if not candidates:
                log.error(
                    "No Visium directory found in %s. "
                    "Set CFG.library_id to the directory name containing filtered_feature_bc_matrix.h5.",
                    CFG.data_dir,
                )
                sys.exit(1)
            CFG.library_id = candidates[0] if candidates[0] else os.path.basename(CFG.data_dir)
            log.info("Auto-detected library_id: '%s'", CFG.library_id)

        visium_dir = os.path.join(CFG.data_dir, CFG.library_id) if CFG.library_id else CFG.data_dir
        log.info("Loading Visium data from: %s", visium_dir)

        adata = sq.read.visium(
            visium_dir,
            library_id=CFG.library_id or None,
            load_images=True,
        )

        if adata is None:
            log.error("sq.read.visium() returned None — check data directory structure")
            sys.exit(1)

        log.info("Visium data loaded: %d spots × %d genes", adata.n_obs, adata.n_vars)

        # Verify spatial coordinates exist
        if 'spatial' not in adata.obsm:
            log.error("No spatial coordinates (obsm['spatial']) found in loaded data")
            sys.exit(1)
        log.info("  Spatial coordinates: shape=%s", adata.obsm['spatial'].shape)

        # Log library_ids stored in uns
        if 'spatial' in adata.uns:
            log.info("  Library IDs in uns['spatial']: %s", list(adata.uns['spatial'].keys()))

    elif CFG.data_format == "h5ad":
        log.info("Loading from h5ad: %s", CFG.input_h5ad)
        backed = getattr(CFG, 'backed', None) or None
        adata = sc.read(CFG.input_h5ad, backed=backed) if backed else sc.read(CFG.input_h5ad)
        log.info("Loaded: %d cells/spots × %d genes", adata.n_obs, adata.n_vars)

        # If spatial coords are missing, try to infer from common keys
        if 'spatial' not in adata.obsm:
            # Check for common coordinate keys
            coord_keys = [k for k in adata.obsm if 'spatial' in k.lower() or 'coord' in k.lower()]
            if coord_keys and CFG.spatial_platform != "visium":
                adata.obsm['spatial'] = adata.obsm[coord_keys[0]]
                log.info("Mapped '%s' → obsm['spatial']", coord_keys[0])
            else:
                log.warning(
                    "No spatial coordinates in obsm. "
                    "Downstream spatial analysis will be limited. "
                    "Set CFG.spatial_platform appropriately."
                )

    else:
        log.error("Unknown data_format for spatial: '%s'. Supported: 'visium', 'h5ad'", CFG.data_format)
        sys.exit(1)

    # ── Ensure CSR format ────────────────────────────────────────────────
    if getattr(CFG, 'force_csr', True) and sp.issparse(adata.X):
        if not sp.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
            log.info("X format converted to CSR")

    # ── Ensure observation names are unique ──────────────────────────────
    if not adata.obs_names.is_unique:
        log.warning("Observation names not unique, calling make_unique()")
        adata.obs_names_make_unique()

    # ── Add in_tissue flag if missing ────────────────────────────────────
    if 'in_tissue' not in adata.obs and CFG.spatial_platform == "visium":
        adata.obs['in_tissue'] = 1
        log.info("Added default 'in_tissue' column (all spots marked as tissue)")

    # ── Save ─────────────────────────────────────────────────────────────
    safe_write(adata, CFG.raw_h5ad, cfg=CFG)
    log.info("Step 00 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
