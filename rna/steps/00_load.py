#!/usr/bin/env python3
"""
Step 00: 加载原始 scRNA-seq 数据
===================================
支持五种输入格式:
  #  1. 10X MTX (CellRanger 输出): sc.read_10x_mtx()
#     多样本: mtx_dir_pattern 匹配样本子目录 (hstack 快路径 / outer 合并)
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
    if sep is not None:
        try:
            peek = pd.read_csv(features_path, nrows=0, sep=sep)
            first_col = peek.columns[0]
            has_header = bool(first_col) and first_col[0].islower()
        except (pd.errors.EmptyDataError, IndexError, TypeError):
            has_header = False
        if has_header:
            return pd.read_csv(features_path, sep=sep)
        return pd.read_csv(features_path, header=None, names=["gene_symbol"], sep=sep)

    import gzip

    opener = gzip.open if features_path.endswith(".gz") else open
    mode = "rt" if features_path.endswith(".gz") else "r"
    with opener(features_path, mode) as f:
        lines = [line.rstrip("\n\r") for line in f if line.rstrip("\n\r")]
    if not lines:
        return pd.DataFrame(columns=["gene_symbol"])

    ncols = len(lines[0].split("\t"))
    if ncols > 1:
        import io

        buf = io.StringIO("\n".join(lines))
        try:
            peek = pd.read_csv(buf, nrows=0, sep="\t")
            first_col = peek.columns[0]
            has_header = bool(first_col) and first_col[0].islower()
        except (pd.errors.EmptyDataError, IndexError):
            has_header = False
        buf.seek(0)
        if has_header:
            return pd.read_csv(buf, sep="\t")
        buf.seek(0)
        return pd.read_csv(buf, header=None, names=["gene_symbol"], sep="\t")

    return pd.DataFrame({"gene_symbol": lines})


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
        safe_write(adata, ambient_pth, cfg=cfg)
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
        safe_write(adata, ambient_pth, cfg=cfg)
        log.info("Ambient-corrected data written to %s", ambient_pth)

    else:
        log.warning("Unknown ambient method: %s — skipping", method)


def _detect_mtx_prefix(mtx_dir: str) -> str:
    """Detect the 10X MTX file prefix inside a sample directory."""
    for f in sorted(os.listdir(mtx_dir)):
        if f.endswith(".mtx.gz") or f.endswith(".mtx"):
            for suffix in ("matrix.mtx.gz", "matrix.mtx", ".mtx.gz", ".mtx"):
                if f.endswith(suffix):
                    return f[: -len(suffix)]
    return ""


def _ensure_10x_features(mtx_dir: str, prefix: str, log) -> str:
    """Return the features file path, converting legacy genes.tsv.gz in place."""
    features = os.path.join(mtx_dir, prefix + "features.tsv.gz")
    if os.path.exists(features):
        return features
    genes = os.path.join(mtx_dir, prefix + "genes.tsv.gz")
    if not os.path.exists(genes):
        raise FileNotFoundError(
            f"Neither {prefix}features.tsv.gz nor {prefix}genes.tsv.gz in {mtx_dir}"
        )
    log.info("  legacy genes.tsv.gz in %s — converting to features.tsv.gz", mtx_dir)
    with gzip.open(genes, "rt") as f_in, gzip.open(features, "wt") as f_out:
        for line in f_in:
            f_out.write(line.rstrip("\n") + "\tGene Expression\n")
    return features


def _read_10x_features(mtx_dir: str, prefix: str) -> pd.DataFrame:
    for name in ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"):
        p = os.path.join(mtx_dir, prefix + name)
        if os.path.exists(p):
            return _read_features_with_header_detection(p, sep="\t")
    raise FileNotFoundError(f"No features/genes file found in {mtx_dir}")


def _features_gene_names(features: pd.DataFrame, gene_symbol_column: str = "") -> list:
    if gene_symbol_column and gene_symbol_column in features.columns:
        return features[gene_symbol_column].astype(str).tolist()
    if "gene_short_name" in features.columns:
        return features["gene_short_name"].astype(str).tolist()
    if "symbol" in features.columns:
        return features["symbol"].astype(str).tolist()
    return features.iloc[:, 0].astype(str).tolist()


def _concat_mtx_batched(adatas: list, batch: int, log):
    """Tree-merge adatas in batches of ``batch`` (bounded peak memory)."""
    import gc

    while len(adatas) > 1:
        groups = [adatas[i : i + batch] for i in range(0, len(adatas), batch)]
        adatas = [sc.concat(g, join="outer", fill_value=0) for g in groups]
        _ = gc.collect()
        log.info("  batched concat: %d group(s) remaining", len(adatas))
    return adatas[0]


def _unique_columns_keep_first(columns, log=None):
    """Dedup column names keeping the first occurrence.

    Mirrors ``expr.columns.duplicated(keep="first")`` in the preprocessed
    branch — the first occurrence wins, later duplicates are dropped.
    """
    seen = set()
    out = []
    n_dup = 0
    for c in columns:
        if c in seen:
            n_dup += 1
            continue
        seen.add(c)
        out.append(c)
    if n_dup and log is not None:
        log.warning("Duplicate gene names: %d — keeping first occurrence", n_dup)
    return out


def _classify_preprocessed_columns(sample: pd.DataFrame) -> list:
    """Three-metric metadata/expression column classification.

    Returns one of ``"M"`` (metadata) / ``"E"`` (expression) per column of
    ``sample``: numeric ratio → sparsity → small-int cardinality → unique
    ratio.  Logic is byte-for-byte the preprocessed branch's auto-detection
    (numeric ratio < 0.5 → meta; >80% zeros → expression; small-int
    categorical → meta; low-cardinality numeric → meta).
    """
    n_sampled = len(sample)
    classifications = []
    for col in sample.columns:
        numeric = pd.to_numeric(sample[col], errors="coerce")
        numeric_ratio = numeric.notna().sum() / n_sampled
        if numeric_ratio < 0.5:
            classifications.append("M")
        else:
            non_na = numeric.dropna()
            zero_frac = (non_na == 0).sum() / len(non_na) if len(non_na) > 0 else 0
            if zero_frac > 0.8:
                classifications.append("E")
                continue
            is_small_int = False
            if len(non_na) > 0:
                if all(v == int(v) for v in non_na):
                    rng = non_na.max() - non_na.min()
                    if rng < 50:
                        is_small_int = True
            if is_small_int:
                classifications.append("M")
            else:
                unique_ratio = numeric.nunique() / n_sampled
                if unique_ratio < 0.5:
                    classifications.append("M")
                else:
                    classifications.append("E")
    return classifications


def _detect_preprocessed_boundary(file_list, sep, log) -> int:
    """Auto-detect the metadata/expression boundary from the first file.

    Reads only the first 100 rows for classification, then derives ``meta_cols``
    (the first expression-column position).  Error paths are identical to the
    previous inline logic.
    """
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

    classifications = _classify_preprocessed_columns(sample)
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
    return meta_cols


def _build_preprocessed_sparse(file_list, sep, meta_cols, log, chunksize=None):
    """Memory-bounded sparse build for the preprocessed branch.

    Replaces the old dense build (concat of all expression frames → dense
    ``expr.values.astype(float32)`` → ``csr_matrix``; measured 29.4 GiB peak
    at 32k cells) with a chunked by-name assembly:

    * the union gene order is read from file headers only (``nrows=0``),
      reproducing ``pd.concat(axis=0)`` order-of-first-appearance;
    * each file is streamed in row-chunks; every chunk is re-aligned to the
      final gene order BY NAME (missing columns → NaN, metis G6) and turned
      into a float32 CSR slice;
    * per-file chunk CSRs are vstacked, then the file CSRs — only the output
      sparse matrix plus a bounded chunk block ever live in memory.

    NaN is preserved and exact-0.0 dropped (``csr_matrix(dense)`` semantics).
    The old ``ignore_index=True`` row numbering is reproduced by streaming
    files in order, so duplicate barcodes get the same ``_{i}`` suffixes.
    """
    # ── header-only column union (pd.concat order-of-first-appearance) ──
    first_columns = list(pd.read_csv(file_list[0], sep=sep, nrows=0).columns)
    combined_columns = list(first_columns)
    seen = set(first_columns)
    for fpath in file_list[1:]:
        for c in pd.read_csv(fpath, sep=sep, nrows=0).columns:
            if c not in seen:
                seen.add(c)
                combined_columns.append(c)

    meta_names = combined_columns[:meta_cols]
    expr_final = _unique_columns_keep_first(combined_columns[meta_cols:], log)
    n_genes = len(expr_final)

    # ~10M cells×genes (~200 MB) of dense per chunk → chunksize by width.
    if chunksize is None:
        chunksize = max(256, min(100_000, int(10_000_000 // max(n_genes, 1))))

    file_csrs = []
    meta_blocks = []
    for fpath in file_list:
        log.info("Loading %s ...", os.path.basename(fpath))
        chunk_csrs = []
        for chunk in pd.read_csv(fpath, sep=sep, chunksize=chunksize):
            expr_block = chunk.reindex(columns=expr_final)
            chunk_csrs.append(sp.csr_matrix(expr_block.values.astype(np.float32)))
            meta_blocks.append(chunk.reindex(columns=meta_names))
        if len(chunk_csrs) == 1:
            file_csrs.append(chunk_csrs[0])
        elif len(chunk_csrs) > 1:
            file_csrs.append(sp.vstack(chunk_csrs, format="csr"))
        else:
            file_csrs.append(sp.csr_matrix((0, n_genes), dtype=np.float32))

    if len(file_csrs) == 1:
        x = file_csrs[0]
    else:
        x = sp.vstack(file_csrs, format="csr")
    log.info("Combined shape: %s", (x.shape[0], x.shape[1]))

    meta = (
        pd.concat(meta_blocks, axis=0, ignore_index=True)
        if meta_blocks
        else pd.DataFrame(columns=meta_names)
    )
    del meta_blocks

    barcodes = meta.iloc[:, 0].values.astype(str)
    if not pd.Index(barcodes).is_unique:
        log.warning("Barcodes not unique — appending row index")
        barcodes = [f"{bc}_{i}" for i, bc in enumerate(barcodes)]

    adata = sc.AnnData(X=x)
    adata.obs_names = barcodes
    adata.var_names = [str(c) for c in expr_final]
    for col_idx in range(1, meta_cols):
        col_name = str(meta.columns[col_idx]).strip()
        adata.obs[col_name] = meta.iloc[:, col_idx].values
    return adata


def _csv_chunk_to_float32(chunk: pd.DataFrame) -> np.ndarray:
    """Per-block float32 coercion of a genes×cells chunk, transposed.

    Mirrors the old branch's whole-file ``astype(np.float32)`` conversion on a
    chunk of the genes×cells table: fully-numeric blocks take the exact fast
    path.  Object blocks (a non-numeric cell anywhere) are coerced per column
    with ``pd.to_numeric(errors="coerce")`` so the offending cell becomes NaN
    — the old whole-file ``astype`` raised ValueError on such cells; the plan
    mandates NaN instead (no crash).
    """

    vals = chunk.values
    if vals.dtype == object:
        vals = chunk.apply(pd.to_numeric, errors="coerce").values
    return vals.T.astype(np.float32)


def _build_csv_matrix_sparse(matrix_file, sep, decimal, log, chunksize=None):
    """Memory-bounded chunked sparse build for the csv_matrix table branch.

    The table file is genes×cells (first column = gene names, header row = cell
    barcodes).  The old build read the WHOLE file into one dense DataFrame and
    materialized a full cells×genes dense transpose before ``force_csr``
    (measured 24.7 GiB peak at 50,954 cells / GSE173180).  This build streams
    gene-blocks instead:

    * a header-only peek (``nrows=0``) sizes the block so a dense block holds
      ≲10M elements (mirrors ``_build_preprocessed_sparse``);
    * each ``pd.read_csv(..., chunksize=...)`` block is a genes×cells slice
      transposed per-chunk into a cells×genes float32 CSR slice;
    * slices are ``sp.hstack``-ed — peak memory is one block's dense transpose
      plus the accumulated sparse output, never a full dense transpose.

    Semantics replicated from the old branch:
    * header/meta handling — ``index_col=0``, default header row (cell
      barcodes become obs_names, the gene-name column becomes var_names);
    * the ``csv_decimal != "."`` re-read is applied by passing ``decimal=`` to
      every chunk read (equivalent to the old branch's second read, which
      overwrote the first);
    * NaN preserved, exact-0.0 dropped (``csr_matrix(dense)`` semantics);
    * obs_names/var_names order and dtype float32 identical to the old build.
    """
    if chunksize is None:
        header = pd.read_csv(matrix_file, index_col=0, sep=sep, nrows=0)
        n_cells = max(len(header.columns), 1)
        chunksize = max(256, min(100_000, int(10_000_000 // n_cells)))

    blocks = []
    var_names = []
    obs_names = None
    for chunk in pd.read_csv(
        matrix_file, index_col=0, sep=sep, decimal=decimal, chunksize=chunksize
    ):
        if obs_names is None:
            obs_names = chunk.columns
        var_names.extend(chunk.index.astype(str))
        blocks.append(sp.csr_matrix(_csv_chunk_to_float32(chunk)))

    if obs_names is None:
        # header-only file: the chunked reader yields a 0-row chunk with columns
        obs_names = pd.read_csv(matrix_file, index_col=0, sep=sep, nrows=0).columns

    x = blocks[0] if len(blocks) == 1 else sp.hstack(blocks, format="csr")
    log.info("CSV shape: %s", (x.shape[1], x.shape[0]))
    adata = sc.AnnData(X=x)
    adata.var_names = list(var_names)
    adata.obs_names = obs_names.astype(str)
    return adata


def _load_multi_sample_10x_mtx(cfg, log):
    """Load multiple 10X MTX sample dirs matched by ``mtx_dir_pattern``.

    Fast path (identical gene sets) uses sparse hstack — O(nnz), no concat
    re-alignment.  Otherwise falls back to one-shot ``sc.concat(join="outer",
    fill_value=0)`` (batched tree merge when ``mtx_concat_batch > 0``).

    Fail-fast: a sample that cannot be loaded aborts with a clear message
    naming the sample and its files.
    """
    import glob as glob_mod
    import re

    mtx_dir = cfg.data_input.mtx_dir
    candidates = sorted(glob_mod.glob(os.path.join(mtx_dir, cfg.data_input.mtx_dir_pattern)))

    sample_dirs = []
    for d in candidates:
        if not os.path.isdir(d):
            continue
        # Any 10X file (matrix/features/barcodes) marks it as a sample dir;
        # corrupt dirs (e.g. missing matrix) are kept so loading fails loudly.
        is_10x = any(
            glob_mod.glob(os.path.join(d, pat))
            for pat in (
                "*matrix.mtx*",
                "*features.tsv*",
                "*genes.tsv*",
                "*barcodes.tsv*",
            )
        )
        if is_10x:
            sample_dirs.append(d)

    if not sample_dirs:
        log.error(
            "mtx_dir_pattern=%r matched no 10X MTX directories under %s",
            cfg.data_input.mtx_dir_pattern,
            mtx_dir,
        )
        sys.exit(1)

    log.info("Multi-sample MTX: %d sample directories", len(sample_dirs))

    sample_names = []
    prefixes = []
    for d in sample_dirs:
        name = os.path.basename(os.path.normpath(d))
        if cfg.data_input.mtx_sample_regex:
            m = re.search(cfg.data_input.mtx_sample_regex, name)
            if m:
                name = m.group(1) if m.groups() else m.group(0)
        sample_names.append(name)
        prefixes.append(_detect_mtx_prefix(d))
        log.info("  %s → sample '%s' (prefix='%s')", os.path.basename(d), name, prefixes[-1])

    for d, prefix in zip(sample_dirs, prefixes):
        try:
            _ensure_10x_features(d, prefix, log)
        except FileNotFoundError as e:
            log.error("Sample dir %s: %s", d, e)
            sys.exit(1)

    gene_sets = []
    for d, prefix in zip(sample_dirs, prefixes):
        try:
            feats = _read_10x_features(d, prefix)
        except FileNotFoundError as e:
            log.error("Sample dir %s: %s", d, e)
            sys.exit(1)
        gene_sets.append(_features_gene_names(feats, cfg.data_input.gene_symbol_column))
        del feats

    first_genes = gene_sets[0]
    identical = all(gs == first_genes for gs in gene_sets[1:])
    log.info(
        "Gene sets %s across %d samples",
        "identical — fast hstack path" if identical else "differ — outer-join concat path",
        len(sample_dirs),
    )

    adatas = []
    for i, d in enumerate(sample_dirs):
        prefix = prefixes[i]
        sname = sample_names[i]
        log.info("  [%d/%d] %s — loading...", i + 1, len(sample_dirs), sname)
        try:
            a = sc.read_10x_mtx(
                d, var_names="gene_symbols", prefix=prefix, cache=False, gex_only=False
            )
        except Exception as e:
            log.error(
                "Failed to load sample '%s' (dir: %s, matrix: %smatrix.mtx[.gz]): %s",
                sname,
                d,
                prefix,
                e,
            )
            sys.exit(1)
        a.X = a.X.tocsr()
        a.obs_names = [f"{bc}-{i}" for bc in a.obs_names]
        a.obs["sample"] = sname
        log.info(
            "  [%d/%d] %s — %d cells × %d genes",
            i + 1,
            len(sample_dirs),
            sname,
            a.n_obs,
            a.n_vars,
        )
        adatas.append(a)

    if identical:
        x_stack = sp.vstack([a.X for a in adatas], format="csr")
        adata = sc.AnnData(
            X=x_stack, obs=pd.concat([a.obs for a in adatas], axis=0), var=adatas[0].var.copy()
        )
        log.info("Merge complete (vstack): %d cells × %d genes", adata.n_obs, adata.n_vars)
    else:
        batch = cfg.data_input.mtx_concat_batch
        if batch and batch > 0:
            log.info("Batched outer-join concat (batch=%d)...", batch)
            adata = _concat_mtx_batched(adatas, batch, log)
        else:
            adata = sc.concat(adatas, join="outer", fill_value=0)
        log.info("Merge complete (outer join): %d cells × %d genes", adata.n_obs, adata.n_vars)

    del adatas

    # sample_map 重命名 (目录名/正则提取名 → 自定义名)
    if cfg.has_sample_mapping() and "sample" in adata.obs:
        map_str = {str(k): v for k, v in cfg.sample_meta.sample_map.items()}
        if map_str:
            mapped = adata.obs["sample"].astype(str).map(map_str)
            n_mapped = int(mapped.notna().sum())
            if n_mapped:
                log.info("sample_map remapped %d/%d cells", n_mapped, len(mapped))
                adata.obs["sample"] = mapped.fillna(adata.obs["sample"]).astype(str)

    _parse_barcodes(adata, cfg, log)
    if "gene_ids" in adata.var:
        adata.var.drop(columns=["gene_ids"], inplace=True)
    return adata


def _extract_sample_name(fname, suffix):
    """Replicate the step-00 sample-name extraction from the pattern suffix."""
    sample_name = os.path.basename(fname)
    if suffix and sample_name.endswith(suffix):
        sample_name = sample_name[: -len(suffix)].rstrip("_")
    elif suffix:
        alt = suffix.lstrip("_")
        if alt and sample_name.endswith(alt):
            sample_name = sample_name[: -len(alt)].rstrip("_")
    else:
        sample_name = os.path.splitext(sample_name)[0]
    return sample_name


def _read_10x_h5_meta(h5_file):
    """Cheap metadata read of a 10x v3 feature-barcode h5 — no X materialized.

    Mirrors scanpy's read (gex_only=True + var_names_make_unique) so the fast
    path assembles X in place with byte-identical obs/var.
    """
    import h5py

    with h5py.File(h5_file, "r") as f:
        m = f["matrix"]
        nnz = int(m["indptr"][-1])
        obs_index = m["barcodes"][:].astype(str)
        names = m["features"]["name"][:].astype(str)
        gene_ids = m["features"]["id"][:].astype(str)
        feature_types = m["features"]["feature_type"][:].astype(str)
        genomes = m["features"]["genome"][:].astype(str)
        extras = {
            k: v[:].astype(bool if v.dtype.kind == "b" else str)
            for k, v in m["features"].items()
            if isinstance(v, h5py.Dataset)
            and k not in ("name", "feature_type", "id", "gene_id", "_all_tag_keys")
        }

    n_obs = len(obs_index)
    gex = feature_types == "Gene Expression"  # gex_only=True (scanpy default)
    names, gene_ids, feature_types, genomes = (
        a[gex] for a in (names, gene_ids, feature_types, genomes)
    )
    extras = {k: v[gex] for k, v in extras.items()}
    if len(names) != len(set(names)):  # anndata var_names_make_unique (join='-')
        counts = {}
        uniq = []
        for n in names:
            c = counts.get(n, 0)
            counts[n] = c + 1
            uniq.append(n if c == 0 else f"{n}-{c}")
        names = uniq

    var = {"gene_ids": gene_ids, "feature_types": feature_types, "genome": genomes, **extras}
    return n_obs, nnz, list(names), obs_index, pd.DataFrame(var, index=names)


def _assemble_csr_inplace(meta, n_vars):
    """Stream per-file X into one preallocated CSR — no vstack transient.

    Only ONE per-file matrix + the preallocated output live at any moment →
    peak ≈ 1× output + 1× largest file. The returned csr shares the arrays.
    """
    import gc
    import warnings as _warn

    total_nnz = sum(nnz for _, _, nnz in meta)
    n_obs = sum(n for _, n, _ in meta)
    idx_dtype = sp.get_index_dtype(maxval=max(total_nnz, n_vars))  # vstack parity
    data = np.empty(total_nnz, dtype=np.float32)
    indices = np.empty(total_nnz, dtype=idx_dtype)
    indptr = np.zeros(n_obs + 1, dtype=idx_dtype)

    row = 0
    off = 0
    for h5_file, n_obs_i, _ in meta:
        with _warn.catch_warnings():
            _warn.simplefilter("ignore", UserWarning)
            a = sc.read_10x_h5(h5_file)
        x = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
        end = off + x.nnz
        data[off:end] = x.data
        indices[off:end] = x.indices
        # cast to int64: per-file indptr is int32, but the global offset can
        # exceed int32 max at 2.8e9 nnz (numpy raises on int32+big-int add)
        indptr[row : row + n_obs_i] = x.indptr[:-1].astype(np.int64) + off
        row += n_obs_i
        off = end
        del a, x
        _ = gc.collect()
    indptr[-1] = off
    return sp.csr_matrix((data, indices, indptr), shape=(n_obs, n_vars))


def _load_multi_sample_10x_h5(cfg, log):
    """Load multiple 10X_h5 files matched by ``h5_file_pattern``.

    PASS 1 reads cheap metadata (shape/nnz/barcodes/features) via h5py — no X
    is materialized.  Identical gene sets (EXACT order) take PASS 2: each
    file's X is streamed into one preallocated CSR so peak ≈ 1× merged matrix
    (vs ~2× for ``sp.vstack``); differing sets use a batched outer-join tree.
    Per-file obs_names get their ``-{i}`` suffix BEFORE any grouping because
    ``sc.concat`` without ``index_unique`` errors on duplicate obs_names.
    """
    import gc
    import glob as glob_mod
    import warnings as _warn

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
    n_total = len(h5_files)
    log.info("Loading from 10X HDF5 (multi-file): %d files", n_total)

    gene_sets, meta, obs_frames, var0 = [], [], [], None
    for i, h5_file in enumerate(h5_files):
        fname = os.path.basename(h5_file)
        log.info("  [%d/%d] %s — loading...", i + 1, n_total, fname)
        try:
            n_obs, nnz, var_names, obs_index, var_df = _read_10x_h5_meta(h5_file)
        except Exception as e:
            log.error("Failed to load 10X_h5 file '%s': %s", fname, e)
            sys.exit(1)
        sample_name = _extract_sample_name(fname, suffix)
        obs = pd.DataFrame(
            {"sample": [sample_name] * n_obs}, index=[f"{b}-{i}" for b in obs_index]
        )
        var0 = var_df if i == 0 else var0
        gene_sets.append(var_names)
        log.info(
            "  [%d/%d] %s — %d cells × %d genes", i + 1, n_total, fname, n_obs, len(var_names)
        )
        meta.append((h5_file, n_obs, nnz))
        obs_frames.append(obs)

    # ORDERED list equality — set equality would silently misalign the vstack
    identical = all(gs == gene_sets[0] for gs in gene_sets[1:])
    log.info(
        "Gene sets %s across %d files",
        "identical — vstack fast path" if identical else "differ — batched outer-join",
        n_total,
    )

    if identical:
        x_stack = _assemble_csr_inplace(meta, len(gene_sets[0]))
        adata = sc.AnnData(X=x_stack, obs=pd.concat(obs_frames, axis=0), var=var0.copy())
        log.info("Merge complete (vstack): %d cells × %d genes", adata.n_obs, adata.n_vars)
    else:
        batch = cfg.data_input.mtx_concat_batch
        if not batch or batch <= 0:
            batch = 10
        log.info("Batched outer-join concat (batch=%d)...", batch)
        adatas = []
        for i, h5_file in enumerate(h5_files):
            with _warn.catch_warnings():
                _warn.simplefilter("ignore", UserWarning)
                a = sc.read_10x_h5(h5_file)
            if a.var_names.duplicated().any():
                a.var_names_make_unique()
            a.obs_names = [f"{bc}-{i}" for bc in a.obs_names]
            a.obs["sample"] = _extract_sample_name(os.path.basename(h5_file), suffix)
            adatas.append(a)
        adata = _concat_mtx_batched(adatas, batch, log)
        log.info("Merge complete (outer join): %d cells × %d genes", adata.n_obs, adata.n_vars)
        del adatas

    del meta, obs_frames, gene_sets
    _ = gc.collect()
    return adata


def _raw_h5ad_compression(cfg):
    """Resolve the step-00 raw-write compression — mirrors safe_write's
    ``step_alias="raw"`` resolution: ``per_step_h5ad_compression["raw"]``
    ("zstd" by schema default) > ``cfg.h5ad_compression`` > "gzip".
    Returns ``(compression, opts)``; opts only apply to the gzip family.
    """
    compression = getattr(cfg, "per_step_h5ad_compression", {}).get("raw")
    if compression is None:
        compression = getattr(cfg, "h5ad_compression", "gzip")
    opts = getattr(cfg, "h5ad_compression_opts", None)
    if opts is not None and not compression.startswith("gzip"):
        opts = None
    return compression, opts


def _filter_downsample_mask(obs, cfg, log) -> np.ndarray:
    """Combined filter + downsample cell mask for the final on-write path.

    Replicates the EXACT selection semantics — incl. the RNG draw order — of
    ``core.downsample.filter_by_config`` (L163-193) then ``downsample_by_config``
    (L196-234): filter first (``sample_keep`` isin on "sample", then ``obs_filter``
    via ``obs.eval`` with the per-row eval fallback), then the downsample RNG over
    the FILTERED cell count with a fresh ``RandomState(seed)``. Returns a bool mask
    of length ``len(obs)`` (kept cells in original order) so the chunked writer can
    drop cells during the block-wise write (metis G4).
    """
    n = len(obs)
    m1 = np.ones(n, dtype=bool)
    # ── filter mask ── (core/downsample.py filter_by_config L163-193)
    sample_keep = getattr(cfg.downsample, "sample_keep", None) or []
    obs_filter = getattr(cfg.downsample, "obs_filter", None) or ""
    if sample_keep:
        m1 &= obs["sample"].isin(sample_keep).to_numpy()
        log.info("sample_keep filter: %d → %d cells", n, int(m1.sum()))
    if obs_filter:
        n_before = int(m1.sum())
        sub = obs[m1]  # eval must see the already-filtered population (order!)
        mask = sub.eval(obs_filter)
        if not isinstance(mask, pd.Series):
            mask = sub.apply(lambda row: eval(obs_filter, {"__builtins__": {}}, dict(row)), axis=1)
        m1[m1] &= np.asarray(mask.values, dtype=bool)
        log.info("obs_filter filter: %d → %d cells", n_before, int(m1.sum()))
    # ── downsample mask ── (core/downsample.py downsample_by_config L196-234)
    target = getattr(cfg, "downsample_target", None)
    n_f = int(m1.sum())
    if target is None or target >= n_f:
        return m1
    strategy = getattr(cfg, "downsample_strategy", "stratified")
    seed = getattr(cfg, "downsample_random_seed", 42)
    max_per = getattr(cfg, "downsample_max_per_sample", 5000)
    sample_key = "sample" if "sample" in obs.columns else None
    if sample_key is None:  # core/downsample._check_sample_col candidate order
        for cand in ["Sample", "samples", "batch", "Batch", "stage", "Stage"]:
            if cand in obs.columns:
                sample_key = cand
                log.info("Using '%s' as group column ('sample' not found)", cand)
                break
    if sample_key is None:
        sample_key = "sample"
    rng = np.random.RandomState(seed)
    f_pos = np.where(m1)[0]
    f_obs = obs.iloc[f_pos]
    counts = f_obs[sample_key].value_counts()
    log.info("Downsampling: target=%d, strategy=%s, seed=%d", target, strategy, seed)
    if strategy == "random":  # core/downsample.downsample_random L30-41
        idx = rng.choice(n_f, size=target, replace=False)
        idx.sort()
    elif strategy in ("stratified", "max_per_sample"):  # L44-85 / L88-122
        per_sample = None
        if strategy == "stratified":
            fractions = counts / n_f
            per_sample = (fractions * target).astype(int)
            remainder = target - int(per_sample.sum())
            if remainder > 0:
                order = np.argsort((fractions * target) - per_sample)[::-1]
                per_sample.iloc[order[:remainder]] += 1
        parts = []
        for sample_name in counts.index:
            s_idx = np.where(f_obs[sample_key].values == sample_name)[0]
            cap = (
                min(int(per_sample[sample_name]), len(s_idx))
                if per_sample is not None
                else min(max_per, len(s_idx))
            )
            if cap < len(s_idx):
                parts.append(rng.choice(s_idx, size=cap, replace=False))
            else:
                parts.append(s_idx)
        idx = np.concatenate(parts)
        idx.sort()
    else:
        log.warning("Unknown downsample strategy '%s', skipping downsample", strategy)
        return m1
    m2 = np.zeros(n_f, dtype=bool)
    m2[idx] = True
    combined = np.zeros(n, dtype=bool)
    combined[f_pos[m2]] = True
    log.info(
        "Downsampled: %d → %d cells (%.1f%%)",
        n_f,
        int(combined.sum()),
        100 * int(combined.sum()) / n_f,
    )
    return combined


def _preflight_step00_meta(cfg, log):
    """Metadata-only pre-scan feeding the step-00 memory guard.

    Reads ONLY file headers / h5 dataset shapes — never the expression
    matrix — so the guard can run before the heavy load (metis G7).
    Returns ``(n_cells, n_genes, nnz, concat_factor)`` or ``None`` when the
    format has no cheap metadata probe (the guard is silently skipped, so
    existing runs are unaffected unless ``execution.memory.guard=block``)."""
    import glob as glob_mod

    fmt = cfg.data_format

    if fmt == "10X_h5":
        import h5py

        h5_dir = getattr(cfg.data_input, "h5_dir", "") or cfg.data_dir
        files = sorted(glob_mod.glob(os.path.join(h5_dir, cfg.data_input.h5_file_pattern)))
        if not files:
            return None  # the load branch reports the missing-file error
        n_cells = nnz = n_genes = 0
        for fp in files:
            try:
                with h5py.File(fp, "r") as f:
                    m = f["matrix"]
                    nnz += int(m["data"].shape[0])  # zero-copy shape read
                    n_cells += int(m["indptr"].shape[0]) - 1
                    n_genes = max(n_genes, int(m["features"]["name"].shape[0]))
            except Exception:
                return None  # let the real load report the corrupt file
        return n_cells, n_genes, nnz, 1.3 if len(files) > 1 else 1.0

    if fmt == "10X_mtx":
        if getattr(cfg.data_input, "mtx_dir_pattern", None):
            return None  # multi-sample MTX: no single cheap probe (dead branch)
        for name in (
            cfg.data_input.mtx_prefix + "matrix.mtx",
            cfg.data_input.mtx_prefix + "matrix.mtx.gz",
        ):
            path = os.path.join(cfg.data_input.mtx_dir, name)
            if os.path.exists(path):
                info = _mtx_header_info(path)
                if info:
                    n_genes, n_cells, nnz = info  # MTX is genes×cells
                    return n_cells, n_genes, nnz, 1.0
        return None

    if fmt == "csv_matrix":
        matrix_file = cfg.data_input.matrix_file
        base = matrix_file[:-3] if matrix_file.endswith(".gz") else matrix_file
        ext = os.path.splitext(base)[1].lower()
        if ext in (".csv", ".tsv", ".txt"):
            # table: header-only read → n_cells; line count → n_genes;
            # nnz estimated from a bounded density sample (cheap).
            sep = getattr(cfg.data_input, "csv_sep", None)
            if not sep:
                try:
                    peek = pd.read_csv(matrix_file, sep=None, engine="python", nrows=1)
                    sep = "\t" if len(peek.columns) > 1 else ","
                except Exception:
                    sep = "\t"
            try:
                header = pd.read_csv(matrix_file, index_col=0, sep=sep, nrows=0)
                n_cells = max(len(header.columns), 1)
                n_genes = max(_count_text_lines(matrix_file) - 1, 1)
                density = _sample_table_density(matrix_file, sep)
                return n_cells, n_genes, int(n_cells * n_genes * density), 1.0
            except Exception:
                return None
        info = _mtx_header_info(matrix_file)  # MTX-mmread sub-branch
        if info:
            n_genes, n_cells, nnz = info
            return n_cells, n_genes, nnz, 1.0
        return None

    if fmt == "preprocessed":
        sep = cfg.data_input.separator or "\t"
        pattern = getattr(cfg.data_input, "file_pattern", "") or "*.tsv.gz"
        files = sorted(glob_mod.glob(os.path.join(cfg.data_dir, pattern)))
        if not files:
            return None
        n_cells = n_genes = 0
        density_total = 0.0
        for fp in files:
            try:
                header = pd.read_csv(fp, sep=sep, nrows=0)
                n_genes = max(n_genes, len(header.columns) - 1)  # minus barcode col
                rows = max(_count_text_lines(fp) - 1, 0)
                n_cells += rows
                density_total += _sample_table_density(fp, sep) * rows
            except Exception:
                return None
        density = density_total / n_cells if n_cells else 0.0
        return n_cells, n_genes, int(n_cells * n_genes * density), 1.0

    return None  # h5ad / unknown formats: no metadata probe


def _mtx_header_info(path):
    """Parse rows/cols/nnz from a (possibly .gz) MatrixMarket header."""
    import gzip as _gzip

    opener = _gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    return int(parts[0]), int(parts[1]), int(parts[2])
    except Exception:
        return None
    return None


def _count_text_lines(path):
    """Count lines of a (possibly .gz) text file; used for cheap row counts."""
    import gzip as _gzip

    opener = _gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as f:
            return sum(1 for _ in f)
    except Exception:
        return 0


def _sample_table_density(path, sep, max_rows: int = 200):
    """Non-zero fraction over a bounded table sample (upper-bound nnz)."""
    try:
        df = pd.read_csv(path, index_col=0, sep=sep, nrows=max_rows)
        vals = df.values
        if vals.dtype == object:
            vals = df.apply(pd.to_numeric, errors="coerce").values
        v32 = vals.astype(np.float32)
        n = v32.size
        if n == 0:
            return 0.0
        return float(np.count_nonzero(~np.isnan(v32) & (v32 != 0))) / n
    except Exception:
        return 0.0


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

    # ── Memory guard: metadata-only pre-scan → estimate step-00 peak ──
    # Reads ONLY file headers / h5 shapes (never X) so the guard runs before
    # the heavy load, mirroring 01_doublet.  Default guard is "warn" —
    # existing runs are unaffected unless execution.memory.guard=block.
    from core.utils import check_memory_guard, estimate_step_peak, resolve_memory_settings

    _mem_policy, _mem_budget, _mem_guard = resolve_memory_settings(cfg)
    _meta = _preflight_step00_meta(cfg, log)
    if _meta is not None:
        _n_cells, _n_genes, _nnz, _concat = _meta
        _est = {
            0: estimate_step_peak(
                0,
                _n_cells,
                _n_genes,
                _nnz,
                policy=_mem_policy,
                budget_bytes=_mem_budget,
                concat_factor=_concat,
            )
        }
        if _mem_budget > 0:
            log.info(
                "[memory-guard] estimated peaks: "
                + ", ".join(f"step{s} ~{g:.0f}GB" for s, g in _est.items())
            )
        check_memory_guard(_est, _mem_budget, _mem_guard, logger_obj=log)

    # ── 3 种加载方式 ──────────────────────────────────────────────
    if cfg.data_format == "10X_mtx":
        if cfg.data_input.mtx_dir_pattern:
            # ── 多样本 MTX (mtx_dir_pattern): 每个样本一个子目录 ──
            adata = _load_multi_sample_10x_mtx(cfg, log)
        else:
            # Legacy 2-column genes.tsv.gz → 3-column features.tsv.gz
            genes_path = os.path.join(
                cfg.data_input.mtx_dir, cfg.data_input.mtx_prefix + "genes.tsv.gz"
            )
            features_path = os.path.join(
                cfg.data_input.mtx_dir, cfg.data_input.mtx_prefix + "features.tsv.gz"
            )
            if not os.path.exists(features_path) and os.path.exists(genes_path):
                log.info(
                    "Detected legacy 2-column genes.tsv.gz — converting to features.tsv.gz..."
                )
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
                cache=False,
                gex_only=False,
            )
            log.info("Loading complete: %d cells × %d genes", adata.n_obs, adata.n_vars)

            # 解析 barcode 后缀 → 样本/阶段映射
            if cfg.has_sample_mapping() or cfg.has_stage_mapping():
                bc_suffix = adata.obs_names.to_series().str.extract(r"-(\d+)$")[0].astype(str)
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
        if matrix_ext in (".csv", ".tsv", ".txt"):
            # True table format (CSV/TSV/TXT): gene x cell, first column = gene names
            log.info("Loading from table matrix: %s", cfg.data_input.matrix_file)
            # Auto-detect separator if not explicitly configured
            csv_sep_cfg = getattr(cfg.data_input, "csv_sep", None)
            if csv_sep_cfg:
                sep = csv_sep_cfg
            else:
                try:
                    peek = pd.read_csv(
                        cfg.data_input.matrix_file, sep=None, engine="python", nrows=1
                    )
                    sep = "	" if len(peek.columns) > 1 else ","
                    log.info("Auto-detected separator: %r", sep)
                except Exception:
                    sep = "	"
                    log.info("Fallback to tab separator")
            decimal = getattr(cfg.data_input, "csv_decimal", ".")
            # Chunked sparse build (memory-bounded): stream gene-blocks, per-block
            # float32 transpose into cells×genes CSR, then hstack — never a full
            # dense transpose (old build peaked 24.7 GiB @ 50,954 cells).
            adata = _build_csv_matrix_sparse(cfg.data_input.matrix_file, sep, decimal, log)
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
            adata = _load_multi_sample_10x_h5(cfg, log)

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
        meta_cols = _detect_preprocessed_boundary(file_list, sep, log)

        # ── Load + build sparse (chunked, memory-bounded) ──
        adata = _build_preprocessed_sparse(file_list, sep, meta_cols, log)

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

    # ── Apply stage_map from 'sample' column (all formats, non-destructive) ──
    # The 10X MTX branch already creates 'stage' from barcode suffix;
    # this covers h5ad, csv_matrix, preprocessed, and 10X_h5 where
    # 'stage' doesn't exist yet but 'sample' column + stage_map are available.
    if cfg.has_stage_mapping() and "stage" not in adata.obs.columns:
        if "sample" in adata.obs.columns:
            adata.obs["stage"] = adata.obs["sample"].astype(str).map(cfg.sample_meta.stage_map)
            n_mapped = adata.obs["stage"].notna().sum()
            log.info(
                "Stage mapping applied: %d/%d cells assigned a stage",
                n_mapped,
                adata.n_obs,
            )
            if n_mapped == 0:
                log.warning(
                    "0 cells matched stage_map — check stage_map keys match 'sample' column values"
                )
            if cfg.sample_meta.stage_order:
                adata.obs["stage"] = pd.Categorical(
                    adata.obs["stage"],
                    categories=cfg.sample_meta.stage_order,
                    ordered=True,
                )
        else:
            log.warning("has_stage_mapping but no 'sample' column in obs; stage not assigned")

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

    # ── 可选细胞过滤 + 降采样（config-driven, 写时过滤, 无整矩阵拷贝）──
    # Selection semantics replicate core/downsample.py exactly; the chunked
    # writer drops cells during the block-wise write (metis G4).
    from anndata.utils import make_index_unique

    from core.utils import safe_write, write_csr_chunked

    cell_mask = _filter_downsample_mask(adata.obs, cfg, log)
    n_keep = int(cell_mask.sum())

    # ── 保存 ──
    log.info("Saving to %s...", cfg.raw_h5ad)
    raw_compression, raw_compression_opts = _raw_h5ad_compression(cfg)

    if n_keep == adata.n_obs:
        # No-op filter/downsample → unchanged write path (identical output).
        if not adata.obs_names.is_unique:
            log.warning("Observation names not unique, calling make_unique()")
            adata.obs_names_make_unique()
        safe_write(adata, cfg.raw_h5ad, cfg=cfg, step_alias="raw")
    else:
        # Filter/downsample applied ON WRITE via the cell mask — no full-matrix
        # copies; X and obs are both sliced, index-aligned.
        out_obs = adata.obs[cell_mask].copy()
        if not out_obs.index.is_unique:
            log.warning("Observation names not unique, calling make_unique()")
            out_obs.index = make_index_unique(out_obs.index)
        write_csr_chunked(
            cfg.raw_h5ad,
            adata.X,
            cell_mask=cell_mask,
            n_obs=n_keep,
            obs=out_obs,
            var=adata.var.copy(),
            compression=raw_compression,
            compression_opts=raw_compression_opts,
        )
    log.info("Step 00 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
