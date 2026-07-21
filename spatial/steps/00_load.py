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

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import scanpy as sc
import scipy.sparse as sp

from core.utils import resolve_config, safe_write, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(cfg.log_dir, "00_load.log"))
    log.info("Step 00: Load raw spatial transcriptomics data")
    log.info("Format: %s, Platform: %s", cfg.data_format, cfg.spatial.platform)

    if os.path.exists(cfg.raw_h5ad):
        log.info("Skip: %s already exists. Delete it to force reload.", cfg.raw_h5ad)
        return

    # ── Load by data_format ──────────────────────────────────────────────
    if cfg.data_format == "visium":
        import squidpy as sq

        # If library_id is not set, auto-detect the first Visium directory
        if not cfg.spatial.library_id:
            candidates = [
                d
                for d in sorted(os.listdir(cfg.data_dir))
                if os.path.isdir(os.path.join(cfg.data_dir, d))
                and os.path.exists(os.path.join(cfg.data_dir, d, "filtered_feature_bc_matrix.h5"))
            ]
            if not candidates:
                # Try the data_dir itself as a visium directory
                if os.path.exists(os.path.join(cfg.data_dir, "filtered_feature_bc_matrix.h5")):
                    candidates = [""]
            if not candidates:
                log.error(
                    "No Visium directory found in %s. "
                    "Set CFG.spatial.library_id to the directory name containing filtered_feature_bc_matrix.h5.",
                    cfg.data_dir,
                )
                sys.exit(1)
            cfg.spatial.library_id = (
                candidates[0] if candidates[0] else os.path.basename(cfg.data_dir)
            )
            log.info("Auto-detected library_id: '%s'", cfg.spatial.library_id)

        visium_dir = (
            os.path.join(cfg.data_dir, cfg.spatial.library_id)
            if cfg.spatial.library_id
            else cfg.data_dir
        )
        log.info("Loading Visium data from: %s", visium_dir)

        adata = sq.read.visium(
            visium_dir,
            library_id=cfg.spatial.library_id or None,
            load_images=True,
        )

        if adata is None:
            log.error("sq.read.visium() returned None — check data directory structure")
            sys.exit(1)

        log.info("Visium data loaded: %d spots × %d genes", adata.n_obs, adata.n_vars)

        # Verify spatial coordinates exist
        if "spatial" not in adata.obsm:
            log.error("No spatial coordinates (obsm['spatial']) found in loaded data")
            sys.exit(1)
        log.info("  Spatial coordinates: shape=%s", adata.obsm["spatial"].shape)

        # Log library_ids stored in uns
        if "spatial" in adata.uns:
            log.info("  Library IDs in uns['spatial']: %s", list(adata.uns["spatial"].keys()))

    elif cfg.data_format == "h5ad":
        log.info("Loading from h5ad: %s", cfg.data_input.input_h5ad)
        backed = getattr(cfg.data_input, "backed", None) or None
        adata = (
            sc.read(cfg.data_input.input_h5ad, backed=backed)
            if backed
            else sc.read(cfg.data_input.input_h5ad)
        )
        log.info("Loaded: %d cells/spots × %d genes", adata.n_obs, adata.n_vars)

        # If spatial coords are missing, try to infer from common keys
        if "spatial" not in adata.obsm:
            # Check for common coordinate keys
            coord_keys = [k for k in adata.obsm if "spatial" in k.lower() or "coord" in k.lower()]
            if coord_keys and cfg.spatial.platform != "visium":
                adata.obsm["spatial"] = adata.obsm[coord_keys[0]]
                log.info("Mapped '%s' → obsm['spatial']", coord_keys[0])
            else:
                log.warning(
                    "No spatial coordinates in obsm. "
                    "Downstream spatial analysis will be limited. "
                    "Set CFG.spatial.platform appropriately."
                )

    else:
        log.error(
            "Unknown data_format for spatial: '%s'. Supported: 'visium', 'h5ad'", cfg.data_format
        )
        sys.exit(1)

    # ── Ensure CSR format ────────────────────────────────────────────────
    if getattr(cfg.execution, "force_csr", True) and sp.issparse(adata.X):
        if not sp.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
            log.info("X format converted to CSR")

    # ── Ensure observation names are unique ──────────────────────────────
    if not adata.obs_names.is_unique:
        log.warning("Observation names not unique, calling make_unique()")
        adata.obs_names_make_unique()

    # ── Add in_tissue flag if missing ────────────────────────────────────
    if "in_tissue" not in adata.obs and cfg.spatial.platform == "visium":
        adata.obs["in_tissue"] = 1
        log.info("Added default 'in_tissue' column (all spots marked as tissue)")

    # ── Save ─────────────────────────────────────────────────────────────
    safe_write(adata, cfg.raw_h5ad, cfg=cfg)
    log.info("Step 00 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
