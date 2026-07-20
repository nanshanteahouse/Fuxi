#!/usr/bin/env python3
"""
Step 00: Load bulk RNA-seq count matrix

Supports three input formats:
  1. count_matrix — CSV/TSV count matrix (genes x samples)
  2. tpm_matrix   — TPM/FPKM expression matrix
  3. h5ad         — Pre-existing h5ad (must be samples x genes)

Output: 00_raw.h5ad
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp

def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(CFG.log_dir, "00_load.log"))
    log.info("Step 00: Load bulk RNA-seq data")
    log.info("Format: %s", CFG.data_format)

    if os.path.exists(CFG.raw_h5ad):
        log.info("Skip: %s already exists. Delete it to force reload.", CFG.raw_h5ad)
        return

    # --- count_matrix format ---
    if CFG.data_format == "count_matrix":
        log.info("Loading from count matrix: %s", CFG.data_input.matrix_file)
        sep = getattr(CFG.data_input, "csv_sep", None)
        if not sep:
            # Auto-detect: try comma first
            import io
            with open(CFG.data_input.matrix_file, 'r') as f:
                first_line = f.readline()
            sep = "," if "," in first_line and "\t" not in first_line else "\t"
        df = pd.read_csv(CFG.data_input.matrix_file, index_col=0, sep=sep)
        log.info("Matrix shape (genes x samples): %s", df.shape)
        # Transpose to AnnData convention: samples x genes
        X = df.values.T
        if not sp.issparse(X):
            X = sp.csr_matrix(X.astype(np.float32))
        adata = sc.AnnData(X=X)
        adata.obs_names = df.columns.astype(str)
        adata.var_names = df.index.astype(str)
        log.info("Loading complete: %d samples x %d genes", adata.n_obs, adata.n_vars)

    # --- tpm_matrix format ---
    elif CFG.data_format == "tpm_matrix":
        log.info("Loading TPM matrix: %s", CFG.data_input.matrix_file)
        sep = getattr(CFG.data_input, "csv_sep", None)
        if not sep:
            with open(CFG.data_input.matrix_file, 'r') as f:
                first_line = f.readline()
            sep = "," if "," in first_line and "\t" not in first_line else "\t"
        df = pd.read_csv(CFG.data_input.matrix_file, index_col=0, sep=sep)
        log.info("Matrix shape (genes x samples): %s", df.shape)
        X = df.values.T
        if not sp.issparse(X):
            X = sp.csr_matrix(X.astype(np.float32))
        adata = sc.AnnData(X=X)
        adata.obs_names = df.columns.astype(str)
        adata.var_names = df.index.astype(str)
        adata.uns["expression_type"] = "tpm"
        log.info("Loading complete: %d samples x %d genes", adata.n_obs, adata.n_vars)

    # --- h5ad format ---
    elif CFG.data_format == "h5ad":
        log.info("Loading from h5ad: %s", CFG.data_input.input_h5ad)
        adata = sc.read(CFG.data_input.input_h5ad)
        # Validate: must be samples x genes (few obs, many vars)
        if adata.n_obs > adata.n_vars:
            log.warning(
                "adata shape (%d x %d) looks transposed (more obs than vars). "
                "Bulk RNA-seq data should be samples x genes. "
                "If your matrix is genes x samples, transpose it before saving as h5ad.",
                adata.n_obs, adata.n_vars,
            )
        log.info("Loading complete: %d samples x %d genes", adata.n_obs, adata.n_vars)

    else:
        log.error("Unknown data_format: %s. Supported: count_matrix, tpm_matrix, h5ad", CFG.data_format)
        sys.exit(1)

    # --- Uniform CSR format ---
    force_csr = getattr(CFG.execution, "force_csr", True)
    if force_csr:
        if sp.issparse(adata.X):
            if not sp.isspmatrix_csr(adata.X):
                adata.X = adata.X.tocsr()
        else:
            adata.X = sp.csr_matrix(adata.X)

    # --- Save ---
    log.info("Saving to %s...", CFG.raw_h5ad)
    from core.utils import safe_write
    if not adata.obs_names.is_unique:
        log.warning("Observation names not unique, calling make_unique()")
        adata.obs_names_make_unique()

    # --- Load sample metadata if provided ---
    meta_file = getattr(CFG.data_input, "metadata_file", "")
    if meta_file and os.path.exists(meta_file):
        log.info("Loading sample metadata from: %s", meta_file)
        meta_df = pd.read_csv(meta_file)
        meta_col = meta_df.columns[0]
        meta_df = meta_df.set_index(meta_col)
        for col in meta_df.columns:
            adata.obs[col] = meta_df[col].values
        log.info("Added metadata columns: %s", list(meta_df.columns))
        log.info("adata.obs columns: %s", list(adata.obs.columns))
    elif meta_file:
        log.warning("Metadata file not found: %s", meta_file)
    else:
        log.info("No metadata file provided")
    safe_write(adata, CFG.raw_h5ad, cfg=CFG)
    log.info("Step 00 complete, took %.1fs", time.time() - t0)

if __name__ == "__main__":
    main()
