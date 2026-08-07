"""Shared chunked CSR → h5ad writer, extracted from the 02_qc filter-on-write pattern.

This module is the write-side workhorse for memory-bounded sparse output across
the step-00 load paths. It exposes two pure, type-annotated functions:

* ``resolve_block_size`` — density-derived streaming block size. Ported from
  ``rna/steps/02_qc.py::_resolve_block_size`` (L638-658) but as a pure function
  over ``(shape, nnz, avail_bytes)`` — no AnnData needed, unit-testable.
* ``write_csr_chunked`` — write a CSR matrix to a brand-new h5ad in blocks, with
  optional on-write row filtering (``cell_mask``), a hidden temp file and an
  atomic ``os.replace`` on success. Ported from
  ``rna/steps/02_qc.py::_write_qc_h5ad`` (L678-789).

Deliberate differences from the 02_qc originals:

* compression is passed **explicitly** — the ``FUXI_QC_COMPR`` env mechanism is
  02_qc-specific and is NOT copied here;
* an all-zero block still emits its selected (empty) rows into ``indptr`` so the
  output row count always matches the declared shape (the 02_qc loop skipped
  zero-nnz blocks, which silently dropped rows);
* no ThreadPool prefetch producer-consumer — blocks are processed sequentially
  (``prefetch=1``); add a pool only when profiling shows a need.

Consumers import the module directly::

    from core.utils._write_chunked import resolve_block_size, write_csr_chunked

(A ``core.utils`` package re-export lands together with the 00_load rewrite that
consumes this module.)
"""

import logging
import os
import time
from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp

from core.utils._io import _resolve_compression_kwargs

log = logging.getLogger(__name__)

__all__ = ["resolve_block_size", "write_csr_chunked"]

# h5py chunk shapes for the resizeable CSR datasets (mirrors 02_qc).
_DATA_CHUNKS = (65536,)
_INDICES_CHUNKS = (65536,)
_IND_PTR_CHUNKS = (4096,)


def resolve_block_size(
    shape: Tuple[int, int],
    nnz: int,
    avail_bytes: Optional[int] = None,
    prefetch: int = 1,
    min_block: int = 50_000,
    max_block: int = 500_000,
    dtype_bytes: int = 12,
) -> int:
    """Density-derived streaming block size.

    Memory model (mirrors ``02_qc._resolve_block_size``): the write loop holds at
    most ``(prefetch + 1)`` blocks, each costing ``block × nnz-per-row ×
    dtype_bytes`` bytes, where ``dtype_bytes`` covers float32 data (4 B) + int64
    indices (8 B) = 12 B. Wall time is insensitive to block size in the
    100k-500k plateau, so memory is the binding constraint::

        block = avail_bytes × 0.4 / (density × dtype_bytes × (prefetch + 1)),
        clamped to [min_block, max_block].

    ``density`` is the average nnz per cell row: ``nnz / shape[0]`` (02_qc samples
    the first 20k rows; here the caller passes the numbers directly so the
    function stays pure). A zero-density matrix returns ``max_block``.

    Args:
        shape: matrix shape ``(n_cells, n_genes)``.
        nnz: total non-zero entries of the matrix.
        avail_bytes: available memory to budget against (``None`` → psutil).
        prefetch: number of additionally-buffered blocks (1 → sequential loop).
        min_block: lower clamp.
        max_block: upper clamp.
        dtype_bytes: bytes per nnz entry (float32 data + int64 indices).

    Returns:
        int block size within ``[min_block, max_block]``.
    """
    if avail_bytes is None:
        import psutil

        avail_bytes = psutil.virtual_memory().available
    density = nnz / max(shape[0], 1)
    bytes_per_row = density * dtype_bytes * (prefetch + 1)
    if bytes_per_row <= 0:
        return max_block
    bs = int(avail_bytes * 0.4 / bytes_per_row)
    return int(min(max(bs, min_block), max_block))


def _masked_blocks(x, mask, block_size):
    """Yield the per-block CSR slice after optional boolean row filtering.

    Blocks whose mask slice is entirely False are skipped — they contribute no
    rows. Rows that are selected but empty (zero nnz) are still yielded so the
    caller can emit their (empty) ``indptr`` entries.
    """
    n = x.shape[0]
    for i in range(0, n, block_size):
        xr = x[i : i + block_size]
        if mask is not None:
            xr = xr[mask[i : i + block_size]]
        if xr.shape[0] == 0:
            continue
        yield xr


def _empty_df(n: int):
    """Minimal empty pandas DataFrame with ``n`` string-index rows."""
    import pandas as pd

    # anndata requires a str index on obs/var (a RangeIndex triggers an
    # ImplicitModificationWarning on read, promoted to an error in this repo).
    return pd.DataFrame(index=pd.Index(np.arange(n).astype(str)))


def write_csr_chunked(
    h5ad_path: str,
    X_csr,  # noqa: N803 (mandated API name)
    cell_mask=None,
    n_obs: Optional[int] = None,
    obs=None,
    var=None,
    compression: str = "gzip",
    compression_opts: Optional[int] = None,
    block_size: Optional[int] = None,
    verify: bool = True,
) -> None:
    """Write a CSR matrix to a brand-new h5ad in blocks, optionally filtering rows.

    ``X`` is written block-by-block into resizeable h5py datasets
    (chunks=(65536,)) with no intermediate dense copy, so peak memory is bounded
    by one block plus the caller's already-resident ``X_csr``. Output goes to a
    hidden temp file next to ``h5ad_path`` and is atomically renamed via
    ``os.replace`` only on success.

    Args:
        h5ad_path: destination path.
        X_csr: scipy sparse matrix (converted to CSR internally).
        cell_mask: optional boolean array of length ``X_csr.shape[0]``; True rows
            are kept in order, False rows are dropped on write. ``None`` keeps
            all rows.
        n_obs: expected output row count. Defaults to the mask sum (or full row
            count when no mask). When given, it sizes ``/X`` shape/``indptr`` and
            a mismatch with the rows actually written raises.
        obs: optional pandas DataFrame written to ``/obs``.
        var: optional pandas DataFrame written to ``/var``.
        compression: ``"gzip"`` | ``"lzf"`` | ``"zstd"`` (zstd requires
            hdf5plugin and falls back to gzip exactly like ``core.utils._io``).
        compression_opts: h5py compression level (gzip family only).
        block_size: fixed int block size, or ``None`` to auto-resolve from the
            matrix density and available memory (see :func:`resolve_block_size`).
        verify: after ``os.replace``, re-read the file block-wise and assert the
            CSR round-trips to ``X_csr[cell_mask]``. Set ``False`` only when even
            block-wise verification is undesirable.

    Returns:
        None. The file is written atomically at ``h5ad_path``.
    """
    import h5py
    from anndata._io.h5ad import write_elem

    if not sp.issparse(X_csr):
        raise TypeError(f"X_csr must be a scipy sparse matrix, got {type(X_csr).__name__}")
    x = X_csr.tocsr()
    n, m = x.shape

    mask = None
    if cell_mask is not None:
        mask = np.asarray(cell_mask, dtype=bool)
        if mask.shape[0] != n:
            raise ValueError(f"cell_mask length {mask.shape[0]} != matrix rows {n}")

    n_keep = n_obs if n_obs is not None else (int(mask.sum()) if mask is not None else n)

    if block_size is None:
        block_size = resolve_block_size((n, m), x.nnz)
    elif not isinstance(block_size, int) or block_size <= 0:
        raise ValueError("block_size must be a positive int or None for auto-resolution")

    comp_kwargs = _resolve_compression_kwargs(compression, compression_opts)

    target_dir = os.path.dirname(h5ad_path) or "."
    os.makedirs(target_dir, exist_ok=True)
    tmp_path = os.path.join(target_dir, f".{os.path.basename(h5ad_path)}.tmp.{os.getpid()}")

    t0 = time.time()
    try:
        with h5py.File(tmp_path, "w") as f:
            f.attrs["encoding-type"] = "anndata"
            f.attrs["encoding-version"] = "0.1.0"
            for key in ("layers", "obsm", "obsp", "varm", "varp"):
                write_elem(f, key, {})
            write_elem(f, "uns", {})
            # read_h5ad's backwards-compat check unconditionally accesses /obs
            # (and read_dispatched walks /var) — write empty DataFrames when the
            # caller has no metadata, so the file is always anndata-readable.
            write_elem(f, "obs", obs if obs is not None else _empty_df(n_keep))
            write_elem(f, "var", var if var is not None else _empty_df(m))

            xg = f.create_group("X")
            xg.attrs["encoding-type"] = "csr_matrix"
            xg.attrs["encoding-version"] = "0.1.0"
            xg.attrs["shape"] = (n_keep, m)
            d_data = xg.create_dataset(
                "data",
                (0,),
                maxshape=(None,),
                dtype=x.dtype,
                chunks=_DATA_CHUNKS,
                **comp_kwargs,
            )
            d_idx = xg.create_dataset(
                "indices",
                (0,),
                maxshape=(None,),
                dtype=np.int64,
                chunks=_INDICES_CHUNKS,
                **comp_kwargs,
            )
            d_iptr = xg.create_dataset(
                "indptr",
                (0,),
                maxshape=(None,),
                dtype=np.int64,
                chunks=_IND_PTR_CHUNKS,
                **comp_kwargs,
            )

            new_indptr = np.zeros(n_keep + 1, dtype=np.int64)
            pos = 0  # rows written so far
            w = 0  # nnz written so far
            for xr in _masked_blocks(x, mask, block_size):
                k = xr.shape[0]
                # Per-row nnz accumulates onto the running total; empty rows get
                # delta 0 so the row count stays exact even for all-zero blocks.
                new_indptr[pos + 1 : pos + k + 1] = new_indptr[pos] + xr.indptr[1:]
                pos += k
                n_w = xr.nnz
                if n_w:
                    d_data.resize(w + n_w, axis=0)
                    d_idx.resize(w + n_w, axis=0)
                    d_data[w : w + n_w] = xr.data
                    d_idx[w : w + n_w] = xr.indices
                    w += n_w
            if pos != n_keep:
                raise ValueError(f"row count mismatch: wrote {pos} rows, expected {n_keep}")
            d_iptr.resize(n_keep + 1, axis=0)
            d_iptr[:] = new_indptr
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise

    os.replace(tmp_path, h5ad_path)
    log.info(
        "  Saved %s (%.1f MB, chunked %.1fs)",
        os.path.basename(h5ad_path),
        os.path.getsize(h5ad_path) / 1e6,
        time.time() - t0,
    )
    if verify:
        _verify_roundtrip(h5ad_path, x, mask, n_keep, block_size)


def _verify_roundtrip(path: str, x, mask, n_keep: int, block_size: int) -> None:
    """Block-wise integrity check: re-derive the expected CSR from the source.

    Memory-bounded: only one block plus the full ``indptr`` (n_keep+1 int64) is
    materialized at a time. Raises IOError on any mismatch.
    """
    import h5py

    _, m = x.shape
    with h5py.File(path, "r") as f:
        xg = f["X"]
        if tuple(xg.attrs["shape"]) != (n_keep, m):
            raise IOError(
                f"shape mismatch: wrote {tuple(xg.attrs['shape'])}, expected {(n_keep, m)}"
            )
        optr = xg["indptr"][...]
        if len(optr) != n_keep + 1:
            raise IOError(f"indptr length {len(optr)} != expected {n_keep + 1}")
        # Expected indptr, re-derived from the source + mask (same math as write).
        eptr = np.zeros(n_keep + 1, dtype=np.int64)
        pos = 0
        for xr in _masked_blocks(x, mask, block_size):
            k = xr.shape[0]
            eptr[pos + 1 : pos + k + 1] = eptr[pos] + xr.indptr[1:]
            pos += k
        if not np.array_equal(optr, eptr):
            raise IOError("indptr does not round-trip to X_csr[cell_mask]")
        # data / indices in blocks (exact equality, NaN-aware on data).
        w = 0
        for xr in _masked_blocks(x, mask, block_size):
            n_w = xr.nnz
            if n_w:
                if not np.array_equal(xg["data"][w : w + n_w], xr.data, equal_nan=True):
                    raise IOError("data does not round-trip to X_csr[cell_mask]")
                if not np.array_equal(xg["indices"][w : w + n_w], xr.indices):
                    raise IOError("indices do not round-trip to X_csr[cell_mask]")
                w += n_w
    log.info("Integrity check: %s verified OK", os.path.basename(path))
