#!/usr/bin/env python3
"""
Step 00: 加载原始 scRNA-seq 数据
===================================
支持四种输入格式:
  1. 10X MTX (CellRanger 输出): sc.read_10x_mtx()
  2. CSV 矩阵 + 元数据文件:     mmread() + pandas
  3. 已有 h5ad:                sc.read()
  4. 10X HDF5 (.h5):           sc.read_10x_h5()

输出: 00_raw.h5ad
"""
import sys, os, time, argparse, gzip, shutil
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp


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
    return pd.read_csv(features_path, header=None, names=['gene_symbol'], sep=sep)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("00_load", os.path.join(CFG.log_dir, "00_load.log"))
    log.info("Step 00: Load raw data")
    log.info("Format: %s", CFG.data_format)

    if os.path.exists(CFG.raw_h5ad):
        log.info("Skip: %s already exists. Delete it to force reload.", CFG.raw_h5ad)
        return

    # ── 3 种加载方式 ──────────────────────────────────────────────
    if CFG.data_format == "10X_mtx":
        # Legacy 2-column genes.tsv.gz → 3-column features.tsv.gz
        genes_path = os.path.join(CFG.mtx_dir, CFG.mtx_prefix + 'genes.tsv.gz')
        features_path = os.path.join(CFG.mtx_dir, CFG.mtx_prefix + 'features.tsv.gz')
        if not os.path.exists(features_path) and os.path.exists(genes_path):
            log.info("Detected legacy 2-column genes.tsv.gz — converting to features.tsv.gz...")
            with gzip.open(genes_path, 'rt') as f_in:
                with gzip.open(features_path, 'wt') as f_out:
                    for line in f_in:
                        f_out.write(line.rstrip('\n') + '\tGene Expression\n')
            log.info("  features.tsv.gz created")

        log.info("Loading from MTX (prefix='%s') ...", CFG.mtx_prefix)
        adata = sc.read_10x_mtx(
            CFG.mtx_dir,
            var_names='gene_symbols',
            prefix=CFG.mtx_prefix,
            cache=True,
            gex_only=False,
        )
        log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

        # 解析 barcode 后缀 → 样本/阶段映射
        if CFG.has_sample_mapping() or CFG.has_stage_mapping():
            bc_suffix = (
                adata.obs_names.to_series()
                .str.extract(r'-(\d+)$')[0]
                .astype(int)
            )
            if CFG.has_sample_mapping():
                adata.obs['sample'] = bc_suffix.map(CFG.sample_map).values
            if CFG.has_stage_mapping():
                adata.obs['stage'] = bc_suffix.map(CFG.stage_map).values
                if CFG.stage_order:
                    adata.obs['stage'] = pd.Categorical(
                        adata.obs['stage'],
                        categories=CFG.stage_order,
                        ordered=True,
                    )
            log.info("Sample mapping applied. Sample distribution:")
            if 'sample' in adata.obs:
                for s, cnt in adata.obs['sample'].value_counts().items():
                    log.info("  %-20s %5d cells", s, cnt)

        # 可配置 barcode 正则解析（非 10X 格式用，如 CSV/per_sample_csv）
        # 支持单字符串或列表（多级回退：first match wins）
        if CFG.barcode_parse_regex:
            log.info("Using barcode regex parsing: %s", CFG.barcode_parse_regex)

            regex_patterns = CFG.barcode_parse_regex
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
                for obs_col, group_key in CFG.barcode_parse_groups.items():
                    if group_key in parsed.columns or (isinstance(group_key, int) and group_key < len(parsed.columns)):
                        adata.obs[obs_col] = parsed[group_key].values
                        log.info("  Extracted %s from barcode", obs_col)

        # 清理 gene_ids 列（如果有）
        if 'gene_ids' in adata.var:
            adata.var.drop(columns=['gene_ids'], inplace=True)

    elif CFG.data_format == "csv_matrix":
        base = CFG.matrix_file[:-3] if CFG.matrix_file.endswith('.gz') else CFG.matrix_file
        matrix_ext = os.path.splitext(base)[1].lower()
        if matrix_ext in ('.csv',):
            # True CSV format: gene × cell, first column = gene names
            log.info("Loading from CSV: %s", CFG.matrix_file)
            sep = getattr(CFG, 'csv_sep', None)
            decimal = getattr(CFG, 'csv_decimal', '.')
            df = pd.read_csv(CFG.matrix_file, index_col=0, sep=sep)
            if decimal != '.':
                df = pd.read_csv(CFG.matrix_file, index_col=0, sep=sep, decimal=decimal)
            log.info("CSV shape: %s", df.shape)
            # Transpose to AnnData convention: cells × genes
            adata = sc.AnnData(X=df.values.T.astype(np.float32))
            adata.var_names = df.index.astype(str)
            adata.obs_names = df.columns.astype(str)
            # Load metadata if barcodes/features files provided
            if CFG.barcodes_file and os.path.exists(CFG.barcodes_file):
                metadata = pd.read_csv(CFG.barcodes_file, index_col=0, sep=sep)
                adata.obs = adata.obs.join(metadata, how='left')
            if CFG.features_file and os.path.exists(CFG.features_file):
                genes = _read_features_with_header_detection(CFG.features_file, sep=sep)
                if len(genes) == adata.n_vars:
                    gene_symbol_col = getattr(CFG, 'gene_symbol_column', '')
                    if gene_symbol_col and gene_symbol_col in genes.columns:
                        adata.var_names = genes[gene_symbol_col].values.astype(str)
                        genes = genes.drop(columns=[gene_symbol_col])
                    elif 'gene_short_name' in genes.columns:
                        adata.var_names = genes['gene_short_name'].values.astype(str)
                        genes = genes.drop(columns=['gene_short_name'])
                    elif 'symbol' in genes.columns:
                        adata.var_names = genes['symbol'].values.astype(str)
                        genes = genes.drop(columns=['symbol'])
                    else:
                        adata.var_names = genes.iloc[:, 0].values.astype(str)
                        genes = genes.drop(columns=[genes.columns[0]])
                    adata.var = genes
        else:
            # Original MTX path (mmread)
            matrix_path = CFG.matrix_file
            # Auto-decompress .gz files (scipy mmread cannot read .gz directly)
            if matrix_path.endswith('.gz'):
                decompressed_path = matrix_path.rstrip('.gz')
                if not os.path.exists(decompressed_path) or os.path.getmtime(matrix_path) > os.path.getmtime(decompressed_path):
                    log.info("Decompressing %s → %s ...", os.path.basename(matrix_path), os.path.basename(decompressed_path))
                    with gzip.open(matrix_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    log.info("  Decompression complete")
                matrix_path = decompressed_path

            log.info("Loading from MTX matrix: %s", matrix_path)
            mtx = mmread(matrix_path)
            log.info("Matrix shape: %s, nnz=%d", mtx.shape, mtx.nnz)
            mtx.data = mtx.data.astype(np.float32)
            mtx = mtx.T.tocsr()

            genes = _read_features_with_header_detection(CFG.features_file)
            gene_symbol_col = getattr(CFG, 'gene_symbol_column', '')
            if gene_symbol_col and gene_symbol_col in genes.columns:
                gene_names = genes[gene_symbol_col].values.astype(str)
            elif 'gene_short_name' in genes.columns:
                gene_names = genes['gene_short_name'].values.astype(str)
            elif 'symbol' in genes.columns:
                gene_names = genes['symbol'].values.astype(str)
            else:
                gene_names = genes.iloc[:, 0].values.astype(str)
            gene_names = pd.Index(gene_names)
            is_dup = gene_names.duplicated(keep=False)
            if is_dup.any():
                log.warning("Duplicate gene names found, adding suffixes to deduplicate")
                gene_names_series = gene_names.to_series().astype(str)
                gene_names_series[is_dup] = (
                    gene_names_series[is_dup]
                    + '_'
                    + gene_names_series.groupby(gene_names_series).cumcount().astype(str)[is_dup]
                )
                gene_names = gene_names_series.values

            metadata = pd.read_csv(CFG.barcodes_file, index_col=0)
            if CFG.meta_columns:
                rename_map = {}
                for target_col, source_col in CFG.meta_columns.items():
                    if source_col in metadata.columns:
                        rename_map[source_col] = target_col
                if rename_map:
                    metadata.rename(columns=rename_map, inplace=True)

            adata = sc.AnnData(X=mtx, obs=metadata, var=pd.DataFrame(index=gene_names))

        log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

    elif CFG.data_format == "h5ad":
        log.info("Loading from h5ad: %s", CFG.input_h5ad)
        backed = getattr(CFG, 'backed', None) or None
        adata = sc.read(CFG.input_h5ad, backed=backed) if backed else sc.read(CFG.input_h5ad)
        log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

    elif CFG.data_format == "10X_h5":
        import glob as glob_mod
        h5_dir = getattr(CFG, 'h5_dir', '') or CFG.data_dir
        pattern = os.path.join(h5_dir, CFG.h5_file_pattern)
        h5_files = sorted(glob_mod.glob(pattern))

        if not h5_files:
            log.error("No .h5 files matching %s found (directory: %s)", CFG.h5_file_pattern, h5_dir)
            sys.exit(1)

        suffix = CFG.h5_file_pattern.lstrip('*')

        if len(h5_files) == 1:
            log.info("Loading from 10X HDF5 (single file): %s", h5_files[0])
            adata = sc.read_10x_h5(h5_files[0])
            sample_name = os.path.basename(h5_files[0])
            if suffix and sample_name.endswith(suffix):
                sample_name = sample_name[:-len(suffix)].rstrip('_')
            elif suffix:
                alt = suffix.lstrip('_')
                if alt and sample_name.endswith(alt):
                    sample_name = sample_name[:-len(alt)].rstrip('_')
            else:
                sample_name = os.path.splitext(sample_name)[0]
            adata.obs['sample'] = sample_name
            log.info("  Sample: %s, %d cells × %d genes", sample_name, adata.n_obs, adata.n_vars)
        else:
            adatas = []
            for f in h5_files:
                log.info("Loading from 10X HDF5: %s", f)
                a = sc.read_10x_h5(f)
                sample_name = os.path.basename(f)
                if suffix and sample_name.endswith(suffix):
                    sample_name = sample_name[:-len(suffix)].rstrip('_')
                elif suffix:
                    alt = suffix.lstrip('_')
                    if alt and sample_name.endswith(alt):
                        sample_name = sample_name[:-len(alt)].rstrip('_')
                else:
                    sample_name = os.path.splitext(sample_name)[0]
                a.obs['sample'] = sample_name
                adatas.append(a)
                log.info("  %s: %d cells", sample_name, a.n_obs)
            for a in adatas:
                if a.var_names.duplicated().any():
                    a.var_names_make_unique()
            adata = sc.concat(adatas, index_unique='-')
            log.info("Merge complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

        if 'gene_ids' in adata.var:
            adata.var.drop(columns=['gene_ids'], inplace=True)

    else:
        log.error("Unknown data_format: %s", CFG.data_format)
        sys.exit(1)

    # ── 统一稀疏格式: CSR (行优先) ──
    if getattr(CFG, 'force_csr', True) and sp.issparse(adata.X):
        if not sp.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
            log.info("X format converted to CSR")

    # ── 可选 float32 精度 ──
    if getattr(CFG, 'use_float32', False):
        adata.X = adata.X.astype('float32', copy=False) if sp.issparse(adata.X) else adata.X
        log.info("X precision converted to float32")

    # ── 保存 ──
    log.info("Saving to %s...", CFG.raw_h5ad)
    from core.utils import safe_write
    if not adata.obs_names.is_unique:
        log.warning("Observation names not unique, calling make_unique()")
        adata.obs_names_make_unique()
    safe_write(adata, CFG.raw_h5ad, cfg=CFG)
    log.info("Step 00 complete, took %.1fs", time.time() - t0)

if __name__ == '__main__':
    main()
