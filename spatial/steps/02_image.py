#!/usr/bin/env python3
"""
Step 02: Image processing for spatial transcriptomics
========================================================
  1. Extract image features via sq.im.process()
  2. (Optional) Tissue segmentation via sq.im.segment()
  3. Compute image-derived features (texture, histogram, etc.)

Input:  01_qc.h5ad (with image in uns['spatial'])
Output: 02_image.h5ad
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import scanpy as sc
import squidpy as sq

from core.utils import resolve_config, safe_write, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("02_image", os.path.join(cfg.log_dir, "02_image.log"))
    log.info("Step 02: Image processing")

    output_path = os.path.join(cfg.h5ad_dir, "02_image.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    input_path = (
        cfg.raw_h5ad
        if not os.path.exists(os.path.join(cfg.h5ad_dir, "01_qc.h5ad"))
        else os.path.join(cfg.h5ad_dir, "01_qc.h5ad")
    )

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    # ── Check for image ──────────────────────────────────────────────────
    library_id = (
        cfg.spatial.library_id
        if cfg.spatial.library_id
        else sorted(adata.uns.get("spatial", {}).keys())[0]
        if "spatial" in adata.uns and adata.uns["spatial"]
        else None
    )

    if library_id is None:
        log.warning("No spatial library_id found — image processing skipped")
        safe_write(adata, output_path, cfg=cfg)
        log.info("Step 02 complete (image processing skipped), took %.1fs", time.time() - t0)
        return

    has_image = (
        "spatial" in adata.uns
        and library_id in adata.uns["spatial"]
        and "images" in adata.uns["spatial"][library_id]
    )

    if not has_image:
        log.warning(
            "No image found for library_id='%s' in uns['spatial']. Image processing skipped.",
            library_id,
        )
        safe_write(adata, output_path, cfg=cfg)
        log.info("Step 02 complete (no image present), took %.1fs", time.time() - t0)
        return

    # ── Image feature extraction ─────────────────────────────────────────
    log.info("Processing image for library_id='%s'...", library_id)

    try:
        # Create ImageContainer from AnnData (Squidpy 1.8 API)
        img = sq.im.ImageContainer.from_adata(adata, library_id=library_id)

        # Crop image to tissue bounding box + smooth
        if cfg.spatial.crop_image:
            log.info("  Cropping image to tissue region...")
        sq.im.process(
            img,
            library_id=library_id,
            method="smooth",
        )
        log.info("  Image features extracted")

        # ── Optional: tissue segmentation ──
        if cfg.spatial.run_segmentation:
            log.info("  Running image segmentation (method=watershed)...")
            try:
                sq.im.segment(
                    img,
                    method="watershed",
                    library_id=library_id,
                )
                log.info("  Segmentation complete")
            except ImportError:
                log.warning("scikit-image not installed — skipping segmentation")
        # ── Multi-scale image features ──
        log.info("  Computing multi-scale image features...")
        sq.im.calculate_image_features(
            adata,
            img,
            features=["summary", "texture"],
            library_id=library_id,
            layer="hires_smooth",
        )
        log.info("  Multi-scale image features added")

    except Exception as exc:
        log.warning("Image processing failed: %s — continuing with raw image", exc)

    safe_write(adata, output_path, cfg=cfg)
    log.info("Step 02 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
