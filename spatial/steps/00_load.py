#!/usr/bin/env python3
"""
Step 00: Load raw spatial transcriptomics data
================================================
Supports:
  1. 10X Visium (SpaceRanger output) — sq.read.visium()
  2. Generic h5ad with spatial coordinates in obsm['spatial']
  3. Seurat wide CSV (counts.csv.gz genes×spots + md.csv.gz) — seurat_csv

Input:  Raw data directory or .h5ad file
Output: 00_raw.h5ad (with spatial coordinates + image in uns)
"""

import argparse
import logging
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from core.utils import resolve_config, safe_write, setup_logger


def _widen_csr_index_dtype(x_csr: sp.csr_matrix) -> sp.csr_matrix:
    """Keep CSR index arrays dtype-consistent (both int64).

    scipy >= 1.18 rejects mixed index dtypes (e.g. int64 ``indptr`` with int32
    ``indices``) in in-place ops like ``eliminate_zeros``.  The seurat_csv
    loader widens ``indptr`` for large merged slides, so ``indices`` must follow
    to keep the matrix internally consistent."""
    x_csr.indices = x_csr.indices.astype(np.int64)
    x_csr.indptr = x_csr.indptr.astype(np.int64)
    return x_csr


def _load_seurat_csv(counts_file: str, md_file: str, log: logging.Logger) -> sc.AnnData:
    """Load a Seurat wide-format counts table + metadata CSV into AnnData.

    ``counts_file`` is a genes×spots wide CSV (rows = ENSG gene IDs, columns = spot
    barcodes ``cell_N_M``) and ``md_file`` the Seurat metadata table (index = spot
    IDs, including ``pixel_x``/``pixel_y``).  Returns an AnnData with X = spots×genes
    CSR float32 (int64 indptr), ``obsm['spatial']`` from the pixel coordinates,
    obs = the full md table (column order/dtypes preserved), and var_names left as
    the raw ENSG IDs (symbol mapping happens downstream in 06_deconvolve).

    Behaviour is pinned by ``tests/test_spatial/test_00_load_seurat_csv.py`` (oracle
    parity): NaN preserved, exact zeros dropped (csr dense semantics), obs equal to
    ``pd.read_csv(md_file, index_col=0)`` under ``assert_frame_equal``.
    """
    log.info("Loading Seurat counts CSV: %s", os.path.basename(counts_file))
    counts = pd.read_csv(counts_file, index_col=0)
    log.info("  counts: %d genes × %d spots", counts.shape[0], counts.shape[1])
    x_csr = sp.csr_matrix(counts.values.T.astype(np.float32))
    x_csr = _widen_csr_index_dtype(x_csr)
    adata = sc.AnnData(X=x_csr)
    adata.var_names = counts.index.astype(str)
    adata.obs_names = counts.columns.astype(str)
    md = pd.read_csv(md_file, index_col=0)
    log.info("  metadata: %d columns × %d spots", md.shape[1], md.shape[0])
    adata.obs = md
    adata.obsm["spatial"] = md[["pixel_x", "pixel_y"]].to_numpy()
    return adata


def _load_multi_slide_seurat(cfg, log: logging.Logger) -> sc.AnnData:
    """Load one or more Seurat CSV slides and merge into a single AnnData.

    Multi-slide: ``CFG.spatial.samples`` (back-compat ``getattr`` — the schema field
    lands in Wave 3) lists per-sample directories under ``cfg.data_dir``, each holding
    ``counts.csv.gz`` + ``md.csv.gz``; an obs ``sample`` column distinguishes slides.
    Without ``samples`` a single slide is auto-discovered from ``counts.csv.gz`` +
    ``md.csv.gz`` directly in ``data_dir`` (no sample column — obs stays exactly the
    md table).  Missing files abort with a clear error.
    """
    data_dir = cfg.data_dir
    samples = list(getattr(cfg.spatial, "samples", None) or [])
    if samples:
        pairs = [
            (
                name,
                os.path.join(data_dir, name, "counts.csv.gz"),
                os.path.join(data_dir, name, "md.csv.gz"),
            )
            for name in samples
        ]
        log.info("Seurat CSV multi-slide: %d sample(s) under %s", len(pairs), data_dir)
    else:
        pairs = [
            (
                os.path.basename(os.path.normpath(data_dir)),
                os.path.join(data_dir, "counts.csv.gz"),
                os.path.join(data_dir, "md.csv.gz"),
            )
        ]
        log.info("Seurat CSV single-slide auto-discovered under %s", data_dir)

    for _name, counts_file, md_file in pairs:
        for f in (counts_file, md_file):
            if not os.path.exists(f):
                log.error("Seurat CSV file not found: %s", f)
                sys.exit(1)

    if len(pairs) == 1 and not samples:
        # Single slide, no sample column — obs is exactly the md table.
        adata = _load_seurat_csv(pairs[0][1], pairs[0][2], log)
        log.info("Loaded: %d spots × %d genes", adata.n_obs, adata.n_vars)
        return adata

    adatas = []
    for i, (name, counts_file, md_file) in enumerate(pairs):
        log.info("  [%d/%d] sample '%s'", i + 1, len(pairs), name)
        a = _load_seurat_csv(counts_file, md_file, log)
        a.obs["sample"] = name
        adatas.append(a)

    gene_sets = [list(a.var_names) for a in adatas]
    identical = all(gs == gene_sets[0] for gs in gene_sets[1:])
    log.info(
        "Gene sets %s across %d slides",
        "identical — vstack fast path" if identical else "differ — outer-join concat",
        len(adatas),
    )
    if identical:
        x_stack = sp.vstack([a.X for a in adatas], format="csr")
        x_stack = _widen_csr_index_dtype(x_stack)
        adata = sc.AnnData(
            X=x_stack,
            obs=pd.concat([a.obs for a in adatas], axis=0),
            var=adatas[0].var.copy(),
        )
    else:
        adata = sc.concat(adatas, join="outer", fill_value=0)
        if sp.issparse(adata.X):
            adata.X = _widen_csr_index_dtype(adata.X.tocsr())
    if all("spatial" in a.obsm for a in adatas):
        adata.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adatas])
    log.info("Loaded: %d spots × %d genes", adata.n_obs, adata.n_vars)
    return adata


def _maybe_map_symbols(adata: sc.AnnData, cfg, log: logging.Logger) -> sc.AnnData:
    """ENSG → gene-symbol var_names when ``spatial.symbol_mapping`` is set.

    The retina KB (and most downstream annotation/marker logic) matches on
    HGNC symbols, so GSE235583-style ENSG counts must be converted to be
    annotatable.  Original IDs are preserved in ``var['ensembl_id']``;
    unmapped/duplicate genes are dropped (mirrors ``ensure_gene_symbols``)."""
    if not getattr(cfg.spatial, "symbol_mapping", False):
        return adata
    if not any(str(v).startswith("ENSG") for v in adata.var_names):
        log.info("symbol_mapping enabled but var_names already symbols — skip")
        return adata
    from core.interaction import ensure_gene_symbols

    log.info("symbol_mapping: converting ENSG var_names → gene symbols (mygene)...")
    adata.var["ensembl_id"] = adata.var_names.astype(str)
    adata = ensure_gene_symbols(adata, log=log, species=cfg.species)
    log.info("  %d genes after symbol mapping", adata.n_vars)
    return adata


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

    elif cfg.data_format == "seurat_csv":
        log.info("Loading Seurat CSV (wide counts + md table)")
        adata = _load_multi_slide_seurat(cfg, log)
        adata = _maybe_map_symbols(adata, cfg, log)

    else:
        log.error(
            "Unknown data_format for spatial: '%s'. Supported: 'visium', 'h5ad', 'seurat_csv'",
            cfg.data_format,
        )
        sys.exit(1)

    # ── Ensure CSR format ────────────────────────────────────────────────
    if getattr(cfg.execution, "force_csr", True) and sp.issparse(adata.X) and adata.X is not None:
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
