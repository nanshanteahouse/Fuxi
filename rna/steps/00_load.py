#!/usr/bin/env python3
"""
Step 00: 加载原始 scRNA-seq 数据
===================================
支持五种输入格式:
  1. 10X MTX (CellRanger 输出): sc.read_10x_mtx()
  2. CSV 矩阵 + 元数据文件:     mmread() + pandas
  3. 已有 h5ad:                sc.read()
  4. 10X HDF5 (.h5):           sc.read_10x_h5()
  5. Preprocessed TSV:         pd.read_csv() + auto split metadata/expression

输出: 00_raw.h5ad
"""

import argparse
import gzip
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmread

from core.utils import resolve_config, setup_logger


def _read_features_with_header_detection(features_path: str, sep=None) -> pd.DataFrame:
    """Read features file, auto-detecting whether it has a header row.

    Standard features files (e.g., from CellRanger) have lowercase column names
    like 'id', 'gene_short_name', 'feature_type'. Headerless files have a gene
    symbol or ID as the first column name (typically starting with uppercase).
    """
    try:
        peek = pd.read_csv(features_path, nrows=0, sep=sep)
        first_col = peek.columns[0]
        # Lowercase first char → standard header; otherwise → headerless data
        has_header = bool(first_col) and first_col[0].islower()
    except (pd.errors.EmptyDataError, IndexError):
        has_header = False
    if has_header:
        return pd.read_csv(features_path, sep=sep)
    return pd.read_csv(features_path, header=None, names=["gene_symbol"], sep=sep)


def _parse_barcodes(adata, cfg, log):
    """Parse barcode names using configurable regex patterns.

    Supports single string or list of patterns (first match wins).
    Extracted groups are added as obs columns per CFG.sample_meta.barcode_parse_groups.
    """
    if not cfg.sample_meta.barcode_parse_regex:
        return

    log.info("Using barcode regex parsing: %s", cfg.sample_meta.barcode_parse_regex)

    regex_patterns = cfg.sample_meta.barcode_parse_regex
    if isinstance(regex_patterns, str):
        regex_patterns = [regex_patterns]

    parsed = None
    for i, pattern in enumerate(regex_patterns):
        candidates = adata.obs_names.to_series().str.extract(pattern)
        if candidates.iloc[:, 0].notna().any():
            log.info("  Regex pattern #%d matched: %s", i + 1, pattern)
            parsed = candidates
            break

    if parsed is None:
        log.warning("  No barcode regex pattern matched; skipping barcode_parse_groups")
    else:
        for obs_col, group_key in cfg.sample_meta.barcode_parse_groups.items():
            if group_key in parsed.columns or (
                isinstance(group_key, int) and group_key < len(parsed.columns)
            ):
                adata.obs[obs_col] = parsed[group_key].values
                log.info("  Extracted %s from barcode", obs_col)


def _run_ambient_correction(adata, cfg, log):
    """Run ambient RNA correction (CellBender or SoupX) on raw count data.

    Writes an intermediate ``ambient_removed.h5ad`` before the main
    raw h5ad is written.  On missing dependency the step is skipped
    with a logged error (does not crash).
    """
    from core.utils import safe_write
    from core.utils._optional import require_cellbender, require_soupx

    method = cfg.ambient.method
    out_dir = os.path.dirname(cfg.raw_h5ad)
    ambient_pth = os.path.join(out_dir, "ambient_removed.h5ad")

    if method == "cellbender":
        try:
            require_cellbender()
        except ImportError:
            log.error("CellBender not installed — skipping ambient correction")
            return
        log.info("Running CellBender ambient RNA removal...")
        log.warning("CellBender GPU training not yet implemented — Phase 2 scope")
        # Phase 1: write pass-through until GPU training is wired
        safe_write(ambient_pth, adata=adata, file_type="h5ad")
        log.info("Ambient-corrected data written to %s", ambient_pth)

    elif method == "soupx":
        try:
            require_soupx()
        except ImportError:
            log.error("SoupX not installed — skipping ambient correction")
            return
        log.info("Running SoupX ambient RNA removal...")
        log.warning("SoupX ambient estimation not yet implemented — Phase 2 scope")
        # Phase 1: write pass-through until SoupX estimation is wired
        safe_write(ambient_pth, adata=adata, file_type="h5ad")
        log.info("Ambient-corrected data written to %s", ambient_pth)

    else:
        log.warning("Unknown ambient method: %s — skipping", method)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(cfg.log_dir, "00_load.log"))
    log.info("Step 00: Load raw data")
    log.info("Format: %s", cfg.data_format)

    if os.path.exists(cfg.raw_h5ad):
        log.info("Skip: %s already exists. Delete it to force reload.", cfg.raw_h5ad)
        return

    # ── 3 种加载方式 ──────────────────────────────────────────────
    if cfg.data_format == "10X_mtx":
        # Legacy 2-column genes.tsv.gz → 3-column features.tsv.gz
        genes_path = os.path.join(
            cfg.data_input.mtx_dir, cfg.data_input.mtx_prefix + "genes.tsv.gz"
        )
        features_path = os.path.join(
            cfg.data_input.mtx_dir, cfg.data_input.mtx_prefix + "features.tsv.gz"
        )
        if not os.path.exists(features_path) and os.path.exists(genes_path):
            log.info("Detected legacy 2-column genes.tsv.gz — converting to features.tsv.gz...")
            with gzip.open(genes_path, "rt") as f_in:
                with gzip.open(features_path, "wt") as f_out:
                    for line in f_in:
                        f_out.write(line.rstrip("\n") + "\tGene Expression\n")
            log.info("  features.tsv.gz created")

        log.info("Loading from MTX (prefix='%s') ...", cfg.data_input.mtx_prefix)
        adata = sc.read_10x_mtx(
            cfg.data_input.mtx_dir,
            var_names="gene_symbols",
            prefix=cfg.data_input.mtx_prefix,
            cache=True,
            gex_only=False,
        )
        log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

        # 解析 barcode 后缀 → 样本/阶段映射
        if cfg.has_sample_mapping() or cfg.has_stage_mapping():
            bc_suffix = adata.obs_names.to_series().str.extract(r"-(\d+)$")[0].astype(int)
            if cfg.has_sample_mapping():
                adata.obs["sample"] = bc_suffix.map(cfg.sample_meta.sample_map).values
            if cfg.has_stage_mapping():
                adata.obs["stage"] = bc_suffix.map(cfg.sample_meta.stage_map).values
                if cfg.sample_meta.stage_order:
                    adata.obs["stage"] = pd.Categorical(
                        adata.obs["stage"],
                        categories=cfg.sample_meta.stage_order,
                        ordered=True,
                    )
            log.info("Sample mapping applied. Sample distribution:")
            if "sample" in adata.obs:
                for s, cnt in adata.obs["sample"].value_counts().items():
                    log.info("  %-20s %5d cells", s, cnt)

        # 可配置 barcode 正则解析
        _parse_barcodes(adata, cfg, log)

        # 清理 gene_ids 列（如果有）
        if "gene_ids" in adata.var:
            adata.var.drop(columns=["gene_ids"], inplace=True)

    elif cfg.data_format == "csv_matrix":
        base = (
            cfg.data_input.matrix_file[:-3]
            if cfg.data_input.matrix_file.endswith(".gz")
            else cfg.data_input.matrix_file
        )
        matrix_ext = os.path.splitext(base)[1].lower()
        if matrix_ext in (".csv",):
            # True CSV format: gene × cell, first column = gene names
            log.info("Loading from CSV: %s", cfg.data_input.matrix_file)
            sep = getattr(cfg.data_input, "csv_sep", ",")
            decimal = getattr(cfg.data_input, "csv_decimal", ".")
            df = pd.read_csv(cfg.data_input.matrix_file, index_col=0, sep=sep)
            if decimal != ".":
                df = pd.read_csv(cfg.data_input.matrix_file, index_col=0, sep=sep, decimal=decimal)
            log.info("CSV shape: %s", df.shape)
            # Transpose to AnnData convention: cells × genes
            adata = sc.AnnData(X=df.values.T.astype(np.float32))
            adata.var_names = df.index.astype(str)
            adata.obs_names = df.columns.astype(str)
            # Load metadata if barcodes/features files provided
            if cfg.data_input.barcodes_file and os.path.exists(cfg.data_input.barcodes_file):
                metadata = pd.read_csv(cfg.data_input.barcodes_file, index_col=0, sep=sep)
                # Apply meta_columns renaming (same as MTX branch below)
                if cfg.sample_meta.meta_columns:
                    rename_map = {}
                    for target_col, source_col in cfg.sample_meta.meta_columns.items():
                        if source_col in metadata.columns:
                            rename_map[source_col] = target_col
                    if rename_map:
                        metadata.rename(columns=rename_map, inplace=True)
                adata.obs = adata.obs.join(metadata, how="left")
            if cfg.data_input.features_file and os.path.exists(cfg.data_input.features_file):
                genes = _read_features_with_header_detection(cfg.data_input.features_file, sep=sep)
                if len(genes) == adata.n_vars:
                    gene_symbol_col = getattr(cfg.data_input, "gene_symbol_column", "")
                    if gene_symbol_col and gene_symbol_col in genes.columns:
                        adata.var_names = genes[gene_symbol_col].values.astype(str)
                        genes = genes.drop(columns=[gene_symbol_col])
                    elif "gene_short_name" in genes.columns:
                        adata.var_names = genes["gene_short_name"].values.astype(str)
                        genes = genes.drop(columns=["gene_short_name"])
                    elif "symbol" in genes.columns:
                        adata.var_names = genes["symbol"].values.astype(str)
                        genes = genes.drop(columns=["symbol"])
                    else:
                        adata.var_names = genes.iloc[:, 0].values.astype(str)
                        genes = genes.drop(columns=[genes.columns[0]])
                    adata.var = genes
        else:
            # Original MTX path (mmread)
            matrix_path = cfg.data_input.matrix_file
            # Auto-decompress .gz files (scipy mmread cannot read .gz directly)
            if matrix_path.endswith(".gz"):
                decompressed_path = matrix_path.rstrip(".gz")
                if not os.path.exists(decompressed_path) or os.path.getmtime(
                    matrix_path
                ) > os.path.getmtime(decompressed_path):
                    log.info(
                        "Decompressing %s → %s ...",
                        os.path.basename(matrix_path),
                        os.path.basename(decompressed_path),
                    )
                    with gzip.open(matrix_path, "rb") as f_in:
                        with open(decompressed_path, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    log.info("  Decompression complete")
                matrix_path = decompressed_path

            log.info("Loading from MTX matrix: %s", matrix_path)
            mtx = mmread(matrix_path)
            log.info("Matrix shape: %s, nnz=%d", mtx.shape, mtx.nnz)
            mtx.data = mtx.data.astype(np.float32)
            mtx = mtx.T.tocsr()

            genes = _read_features_with_header_detection(cfg.data_input.features_file)
            gene_symbol_col = getattr(cfg.data_input, "gene_symbol_column", "")
            if gene_symbol_col and gene_symbol_col in genes.columns:
                gene_names = genes[gene_symbol_col].values.astype(str)
            elif "gene_short_name" in genes.columns:
                gene_names = genes["gene_short_name"].values.astype(str)
            elif "symbol" in genes.columns:
                gene_names = genes["symbol"].values.astype(str)
            else:
                gene_names = genes.iloc[:, 0].values.astype(str)
            gene_names = pd.Index(gene_names)
            is_dup = gene_names.duplicated(keep=False)
            if is_dup.any():
                log.warning("Duplicate gene names found, adding suffixes to deduplicate")
                gene_names_series = gene_names.to_series().astype(str)
                gene_names_series[is_dup] = (
                    gene_names_series[is_dup]
                    + "_"
                    + gene_names_series.groupby(gene_names_series).cumcount().astype(str)[is_dup]
                )
                gene_names = gene_names_series.values

            metadata = pd.read_csv(cfg.data_input.barcodes_file, index_col=0)
            if cfg.sample_meta.meta_columns:
                rename_map = {}
                for target_col, source_col in cfg.sample_meta.meta_columns.items():
                    if source_col in metadata.columns:
                        rename_map[source_col] = target_col
                if rename_map:
                    metadata.rename(columns=rename_map, inplace=True)

            adata = sc.AnnData(X=mtx, obs=metadata, var=pd.DataFrame(index=gene_names))

        _parse_barcodes(adata, cfg, log)
        log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

    elif cfg.data_format == "h5ad":
        log.info("Loading from h5ad: %s", cfg.data_input.input_h5ad)
        backed = getattr(cfg.data_input, "backed", None) or None
        adata = (
            sc.read(cfg.data_input.input_h5ad, backed=backed)
            if backed
            else sc.read(cfg.data_input.input_h5ad)
        )
        log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)
        _parse_barcodes(adata, cfg, log)

    elif cfg.data_format == "10X_h5":
        import glob as glob_mod

        h5_dir = getattr(cfg.data_input, "h5_dir", "") or cfg.data_dir
        pattern = os.path.join(h5_dir, cfg.data_input.h5_file_pattern)
        h5_files = sorted(glob_mod.glob(pattern))

        if not h5_files:
            log.error(
                "No .h5 files matching %s found (directory: %s)",
                cfg.data_input.h5_file_pattern,
                h5_dir,
            )
            sys.exit(1)

        suffix = cfg.data_input.h5_file_pattern.lstrip("*")

        if len(h5_files) == 1:
            log.info("Loading from 10X HDF5 (single file): %s", h5_files[0])
            adata = sc.read_10x_h5(h5_files[0])
            sample_name = os.path.basename(h5_files[0])
            if suffix and sample_name.endswith(suffix):
                sample_name = sample_name[: -len(suffix)].rstrip("_")
            elif suffix:
                alt = suffix.lstrip("_")
                if alt and sample_name.endswith(alt):
                    sample_name = sample_name[: -len(alt)].rstrip("_")
            else:
                sample_name = os.path.splitext(sample_name)[0]
            adata.obs["sample"] = sample_name
            log.info("  Sample: %s, %d cells × %d genes", sample_name, adata.n_obs, adata.n_vars)
        else:
            adatas = []
            for f in h5_files:
                log.info("Loading from 10X HDF5: %s", f)
                a = sc.read_10x_h5(f)
                sample_name = os.path.basename(f)
                if suffix and sample_name.endswith(suffix):
                    sample_name = sample_name[: -len(suffix)].rstrip("_")
                elif suffix:
                    alt = suffix.lstrip("_")
                    if alt and sample_name.endswith(alt):
                        sample_name = sample_name[: -len(alt)].rstrip("_")
                else:
                    sample_name = os.path.splitext(sample_name)[0]
                a.obs["sample"] = sample_name
                adatas.append(a)
                log.info("  %s: %d cells", sample_name, a.n_obs)
            for a in adatas:
                if a.var_names.duplicated().any():
                    a.var_names_make_unique()
            adata = sc.concat(adatas, index_unique="-")
            log.info("Merge complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

        if "gene_ids" in adata.var:
            adata.var.drop(columns=["gene_ids"], inplace=True)
        _parse_barcodes(adata, cfg, log)

    elif cfg.data_format == "preprocessed":
        import glob as glob_mod

        sep = cfg.data_input.separator if cfg.data_input.separator else ""
        pattern = getattr(cfg.data_input, "file_pattern", "") or "*.tsv.gz"

        file_list = sorted(glob_mod.glob(os.path.join(cfg.data_dir, pattern)))
        if not file_list:
            log.error("No files matching '%s' found in %s", pattern, cfg.data_dir)
            sys.exit(1)

        log.info("Found %d files matching '%s'", len(file_list), pattern)

        # ── Auto-detect separator if not configured ──
        if not sep:
            log.info("Auto-detecting separator from %s ...", os.path.basename(file_list[0]))
            try:
                if file_list[0].endswith(".gz"):
                    import gzip as _gzip

                    with _gzip.open(file_list[0], "rt", encoding="utf-8", errors="replace") as fh:
                        peek = fh.readline()
                else:
                    with open(file_list[0], "r", encoding="utf-8", errors="replace") as fh:
                        peek = fh.readline()
                if "," in peek and "\t" not in peek:
                    sep = ","
                else:
                    sep = "\t"
                log.info("  Separator: %s", repr(sep))
            except Exception:
                log.info("  Using default tab separator")
                sep = "\t"

        # ── Auto-detect metadata/expression boundary ──
        log.info("Auto-detecting meta/expr boundary from %s ...", os.path.basename(file_list[0]))
        try:
            sample = pd.read_csv(file_list[0], sep=sep, nrows=100)
        except Exception as e:
            log.error("Failed to read sample from %s: %s", file_list[0], e)
            sys.exit(1)

        n_sampled = len(sample)
        if n_sampled < 2:
            log.error("Sample has %d rows — too few to detect boundary", n_sampled)
            sys.exit(1)

        # Three-metric classification: numeric ratio → sparsity → cardinality
        classifications = []
        for col in sample.columns:
            numeric = pd.to_numeric(sample[col], errors="coerce")
            numeric_ratio = numeric.notna().sum() / n_sampled
            if numeric_ratio < 0.5:
                classifications.append("M")  # mostly non-numeric → metadata
            else:
                non_na = numeric.dropna()
                # Step 1: Sparsity check — if >80% of values are zero, it's expression
                zero_frac = (non_na == 0).sum() / len(non_na) if len(non_na) > 0 else 0
                if zero_frac > 0.8:
                    classifications.append("E")  # sparse numeric → expression
                    continue
                # Step 2: Small-integer categorical metadata (e.g. cluster labels)
                is_small_int = False
                if len(non_na) > 0:
                    if all(v == int(v) for v in non_na):
                        rng = non_na.max() - non_na.min()
                        if rng < 50:
                            is_small_int = True
                if is_small_int:
                    classifications.append("M")  # integer-coded meta
                else:
                    # Step 3: Cardinality — high-cardinality numeric = expression
                    unique_ratio = numeric.nunique() / n_sampled
                    if unique_ratio < 0.5:
                        classifications.append("M")  # low-cardinality meta
                    else:
                        classifications.append("E")  # expression
        meta_cols = 0
        for i, c in enumerate(classifications):
            if c == "E":
                meta_cols = i
                break

        log.info(
            "Column classification: %d meta + %d expression",
            meta_cols,
            len(classifications) - meta_cols,
        )

        if meta_cols < 1:
            # Check if truly no metadata or all metadata
            if all(c == "E" for c in classifications):
                log.error(
                    "All columns classified as expression — this looks like a pure count matrix. "
                    "Use data_format='csv_matrix' instead."
                )
            else:
                log.error(
                    "Auto-detection: first column is expression data. "
                    "Use data_format='csv_matrix' for count matrices."
                )
            sys.exit(1)

        if meta_cols >= sample.shape[1] - 2:
            log.error(
                "%d metadata columns leaves only %d expression columns — doesn't look ",
                "like a gene expression matrix.",
                meta_cols,
                sample.shape[1] - meta_cols,
            )
            sys.exit(1)

        log.info(
            "Detected: %d metadata cols, %d expression cols (~%d genes)",
            meta_cols,
            sample.shape[1] - meta_cols,
            sample.shape[1] - meta_cols,
        )

        # ── Load and concat all files ──
        all_dfs = []
        for fpath in file_list:
            log.info("Loading %s ...", os.path.basename(fpath))
            all_dfs.append(pd.read_csv(fpath, sep=sep))

        combined = pd.concat(all_dfs, axis=0, ignore_index=True)
        log.info("Combined shape: %s", combined.shape)
        del all_dfs

        # ── Separate metadata and expression ──
        meta = combined.iloc[:, :meta_cols].copy()
        expr = combined.iloc[:, meta_cols:]

        # Make gene names unique
        if expr.columns.duplicated().any():
            dup_count = expr.columns.duplicated().sum()
            log.warning("Duplicate gene names: %d — keeping first occurrence", dup_count)
            expr = expr.loc[:, ~expr.columns.duplicated(keep="first")]

        # Build AnnData
        barcodes = meta.iloc[:, 0].values.astype(str)
        if not pd.Index(barcodes).is_unique:
            log.warning("Barcodes not unique — appending row index")
            barcodes = [f"{bc}_{i}" for i, bc in enumerate(barcodes)]

        x = expr.values.astype(np.float32)
        adata = sc.AnnData(X=sp.csr_matrix(x))
        adata.obs_names = barcodes
        adata.var_names = expr.columns.astype(str)

        # Remaining metadata columns → obs
        for col_idx in range(1, meta_cols):
            col_name = str(meta.columns[col_idx]).strip()
            adata.obs[col_name] = meta.iloc[:, col_idx].values

        # Apply meta_columns renaming if configured
        if cfg.sample_meta.meta_columns:
            rename_map = {}
            for target_col, source_col in cfg.sample_meta.meta_columns.items():
                if source_col in adata.obs.columns:
                    rename_map[source_col] = target_col
            if rename_map:
                log.info("Obs column renaming via meta_columns: %s", rename_map)
                adata.obs.rename(columns=rename_map, inplace=True)

        # Auto-categorize string columns
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object:
                adata.obs[col] = adata.obs[col].astype("category")

        log.info(
            "Loading complete: %d cells × %d genes, %d obs columns",
            adata.n_obs,
            adata.n_vars,
            len(adata.obs.columns),
        )
        _parse_barcodes(adata, cfg, log)

    else:
        log.error("Unknown data_format: %s", cfg.data_format)
        sys.exit(1)

    # ── 统一稀疏格式: CSR (行优先) ──
    # Handle backed mode before CSR conversion
    if adata.isbacked:
        log.info("Backed mode detected — loading fully into memory for processing")
        adata = adata.to_memory()
    force_csr = getattr(cfg.execution, "force_csr", True)
    if force_csr and adata.X is not None:
        if sp.issparse(adata.X):
            if not sp.isspmatrix_csr(adata.X):
                adata.X = adata.X.tocsr()
                log.info("X format converted to CSR")
        else:
            log.info("Converting dense X (shape=%s) to CSR sparse...", adata.X.shape)
            adata.X = sp.csr_matrix(adata.X)
            log.info("  CSR conversion complete")

    # ── 可选 float32 精度 ──
    if getattr(cfg.execution, "use_float32", False) and adata.X is not None:
        adata.X = adata.X.astype("float32", copy=False) if sp.issparse(adata.X) else adata.X
        log.info("X precision converted to float32")

    # ── 可选 ambient RNA 校正 ──
    if cfg.ambient.run and cfg.data_format in ("10X_h5", "10X_mtx"):
        _run_ambient_correction(adata, cfg, log)

    # ── 可选细胞过滤 + 降采样（config-driven） ──
    from core.downsample import downsample_by_config, filter_by_config

    adata = filter_by_config(adata, cfg, log)
    adata = downsample_by_config(adata, cfg, log)

    # ── 保存 ──
    log.info("Saving to %s...", cfg.raw_h5ad)
    from core.utils import safe_write

    if not adata.obs_names.is_unique:
        log.warning("Observation names not unique, calling make_unique()")
        adata.obs_names_make_unique()
    safe_write(adata, cfg.raw_h5ad, cfg=cfg)
    log.info("Step 00 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
