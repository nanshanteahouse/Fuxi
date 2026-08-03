"""Incremental h5ad write engine — in-place obs/obsm/obsp/uns appends.

Backed by :func:`anndata.io.write_elem` / :func:`anndata.io.read_elem`
(anndata ≥ 0.13). Appends to an existing h5ad file without rewriting X;
``X`` / ``layers`` / ``raw`` are never touched.

Design notes
------------
- **Overwrite semantics**: keys passed in are always overwritten, so re-running
  a step replaces stale results (e.g. an old ``X_umap`` / ``leiden_*``). This
  deliberately rejects an "only write missing keys" scheme, which would leave
  stale results behind on re-runs.
- **Compression consistency**: appends use the same compression as the target
  file's ``X`` dataset (gzip for RNA, uncompressed for ATAC) so the file keeps
  a uniform layout. Probe failures fall back to uncompressed appends.
- **No WAL**: HDF5 appends mutate B-tree metadata in place — a crash mid-write
  can corrupt the whole file. Every append is therefore re-read and verified
  (:func:`verify_incremental_write`). The failure *policy* (delete corrupt
  copy vs. restore a ``.bak``) is the caller's choice — see ``_io.py``.
"""

from __future__ import annotations

import logging
import os
from typing import Mapping, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

__all__ = ["write_h5ad_incremental", "verify_incremental_write"]


_ZSTD_FILTER_ID = 32015  # HDF5 registered filter id for Zstandard (hdf5plugin)


def _detect_x_compression(h5ad_path: str) -> Optional[str]:
    """Return the compression of the target file's ``X`` dataset.

    ``None`` means the file is uncompressed (ATAC writes with
    ``compression_override=None``). Any probe failure falls back to
    uncompressed appends.

    h5py only knows gzip/lzf/szip by name — third-party filters (zstd
    #32015) surface as the placeholder ``"unknown"``. Probe the dataset's
    create plist for the registered filter id so zstd files append with
    the same compression instead of crashing on ``compression="unknown"``.
    """
    import h5py

    try:
        with h5py.File(h5ad_path, "r") as f:
            if "X" not in f:
                return None
            x = f["X"]
            # Sparse X is stored as a Group (csr_matrix encoding); only
            # Dataset objects expose .compression.
            if isinstance(x, h5py.Dataset):
                comp = x.compression
                if comp in ("gzip", "lzf"):
                    return comp
                if comp == "unknown":
                    # Third-party filter: read the filter id from the plist.
                    try:
                        plist = x.id.get_create_plist()
                        code, _flags, _values = plist.get_filter(0)
                        if code == _ZSTD_FILTER_ID:
                            return "zstd"
                    except Exception:
                        pass
                return None
            return None
    except (OSError, AttributeError):
        return None


def _write_kwargs(compression: Optional[str]) -> dict:
    """dataset_kwargs for write_elem — only pass compression when set.

    zstd is not a h5py string name: pass the filter the same way
    ``core.utils._io._resolve_compression_kwargs`` does (dict() yields
    ``compression=32015`` + opts, which h5py accepts once hdf5plugin has
    registered the filter). Any registration failure falls back to
    uncompressed appends (safe).
    """
    if not compression:
        return {}
    if compression == "zstd":
        try:
            import hdf5plugin  # registers filter #32015 with the HDF5 lib

            return dict(hdf5plugin.Zstd(clevel=1))
        except Exception:
            return {}
    return {"compression": compression}


def _to_write_elem(series: pd.Series):
    """Convert a pandas Series into a write_elem-compatible element.

    - categorical → the ``Categorical`` itself (write_elem picks the minimal
      codes dtype — int8/int16/int32 by n_categories; fixes int8 wraparound
      at 127 categories).
    - other pandas extension arrays (string / Int64 / boolean) → the array
      (written as nullable mask+values; fixes NaN→"nan" string corruption).
    - plain numpy dtypes → ``.to_numpy()`` preserving the source dtype.
    - object/mixed → promoted to nullable string (missing kept as ``<NA>``).
    """
    dtype = series.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        return series.array
    if isinstance(dtype, pd.api.extensions.ExtensionDtype):
        return series.array
    if pd.api.types.is_bool_dtype(dtype):
        return series.to_numpy(dtype=bool)
    if pd.api.types.is_integer_dtype(dtype) or pd.api.types.is_float_dtype(dtype):
        return series.to_numpy()
    return series.astype(pd.StringDtype()).array


def _write_obs_columns(
    obs_group, obs_df: pd.DataFrame, dataset_kwargs: Mapping[str, object]
) -> None:
    """Write every obs column via write_elem and refresh column-order.

    write_elem does NOT maintain the ``/obs/column-order`` attribute; without
    it anndata's reader drops the columns, so new names must be appended here.
    """
    from anndata.io import write_elem

    for col in obs_df.columns:
        write_elem(obs_group, col, _to_write_elem(obs_df[col]), dataset_kwargs=dataset_kwargs)
    existing = [str(x) for x in obs_group.attrs.get("column-order", [])]
    merged = existing + [c for c in obs_df.columns if c not in existing]
    obs_group.attrs["column-order"] = merged


# ── verification helpers ─────────────────────────────────────────────────


def _obs_elem_equal(actual, series: pd.Series) -> bool:
    """Compare a read-back obs element with its source pandas Series."""
    dtype = series.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        if not isinstance(actual, pd.Categorical):
            return False
        return bool(
            np.array_equal(actual.codes, series.cat.codes.to_numpy())
            and np.array_equal(
                actual.categories.to_numpy(dtype=object),
                series.cat.categories.to_numpy(dtype=object),
            )
        )
    if isinstance(dtype, pd.api.extensions.ExtensionDtype):
        # StringArray / IntegerArray / BooleanArray have NA-aware .equals()
        try:
            return bool(pd.array(actual).equals(series.array))
        except (TypeError, ValueError):
            return False
    if pd.api.types.is_object_dtype(dtype):
        # object/mixed columns are promoted to nullable string on write
        expected = series.astype(pd.StringDtype()).array
        return bool(pd.array(actual, dtype="string").equals(expected))
    expected = series.to_numpy()
    actual_np = np.asarray(actual)
    if actual_np.shape != expected.shape or actual_np.dtype != expected.dtype:
        return False
    return bool(np.array_equal(actual_np, expected, equal_nan=True))


def _arrays_equal(actual, expected) -> bool:
    """Shape/dtype/value equality for dense arrays (NaN-tolerant)."""
    a, b = np.asarray(actual), np.asarray(expected)
    return bool(a.shape == b.shape and a.dtype == b.dtype and np.array_equal(a, b, equal_nan=True))


def _sparse_equal(actual, expected) -> bool:
    """Exact CSR-layout equality for sparse matrices (NaN-tolerant)."""
    a, b = actual.tocsr(), expected.tocsr()
    return bool(
        a.shape == b.shape
        and a.dtype == b.dtype
        and np.array_equal(a.indices, b.indices)
        and np.array_equal(a.indptr, b.indptr)
        and np.array_equal(a.data, b.data, equal_nan=True)
    )


def _uns_values_equal(actual, expected) -> bool:
    """Recursive equality for arbitrary ``uns`` values."""
    if isinstance(expected, dict) and isinstance(actual, dict):
        return set(actual.keys()) == set(expected.keys()) and all(
            _uns_values_equal(actual[k], expected[k]) for k in expected
        )
    if sp.issparse(actual) or sp.issparse(expected):
        return sp.issparse(actual) and sp.issparse(expected) and _sparse_equal(actual, expected)
    if isinstance(expected, pd.DataFrame) and isinstance(actual, pd.DataFrame):
        return expected.equals(actual)
    if isinstance(expected, (list, tuple)) or isinstance(actual, (list, tuple)):
        return _arrays_equal(actual, expected)
    if isinstance(actual, np.ndarray) or isinstance(expected, np.ndarray):
        return _arrays_equal(actual, expected)
    if expected is None or actual is None:
        return expected is None and actual is None
    if isinstance(expected, float):
        return bool((np.isnan(expected) and np.isnan(actual)) or expected == actual)
    return actual == expected


def verify_incremental_write(
    h5ad_path: str,
    obs: Optional[pd.DataFrame] = None,
    obsm: Optional[Mapping[str, object]] = None,
    obsp: Optional[Mapping[str, object]] = None,
    uns: Optional[Mapping[str, object]] = None,
) -> None:
    """Re-open the file and verify every written key round-trips.

    Raises ``RuntimeError`` on the first mismatch. This is the integrity gate
    for incremental appends — HDF5 has no WAL, so a silent corruption here
    would poison every downstream reader.
    """
    import h5py
    from anndata.io import read_elem

    def _fail(what: str) -> RuntimeError:
        return RuntimeError(f"Incremental write integrity FAILED: {what}")

    with h5py.File(h5ad_path, "r") as f:
        if obs is not None:
            obs_group = f["obs"]
            for col in obs.columns:
                if col not in obs_group:
                    raise _fail(f"obs column '{col}' missing from /obs")
                if not _obs_elem_equal(read_elem(obs_group[col]), obs[col]):
                    raise _fail(f"obs column '{col}' does not round-trip")
        for label, group_name, mapping in (
            ("obsm", "obsm", obsm),
            ("obsp", "obsp", obsp),
        ):
            if not mapping:
                continue
            group = f[group_name]
            for k, expected in mapping.items():
                if k not in group:
                    raise _fail(f"{label} key '{k}' missing from /{group_name}")
                actual = read_elem(group[k])
                ok = (
                    _sparse_equal(actual, expected)
                    if sp.issparse(expected)
                    else _arrays_equal(actual, expected)
                )
                if not ok:
                    raise _fail(f"{label} key '{k}' does not round-trip")
        if uns:
            uns_group = f["uns"]
            for k, expected in uns.items():
                if k not in uns_group:
                    raise _fail(f"uns key '{k}' missing from /uns")
                if not _uns_values_equal(read_elem(uns_group[k]), expected):
                    raise _fail(f"uns key '{k}' does not round-trip")


# ── engine ───────────────────────────────────────────────────────────────


def write_h5ad_incremental(
    h5ad_path: str,
    obs: Optional[pd.DataFrame] = None,
    obsm: Optional[Mapping[str, object]] = None,
    obsp: Optional[Mapping[str, object]] = None,
    uns: Optional[Mapping[str, object]] = None,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Append obs / obsm / obsp / uns to an existing h5ad file in place.

    Uses :func:`anndata.io.write_elem` — never touches ``X`` / ``layers`` /
    ``raw``. Keys passed in are always overwritten (re-runs replace stale
    results). New datasets use the same compression as the file's ``X``.
    The whole append is read back and verified afterwards.

    Parameters
    ----------
    h5ad_path : str
        Existing h5ad file to append to (created files are rejected — an
        append to a fresh file would produce an unreadable h5ad).
    obs : DataFrame, optional
        Obs columns to write. Must have the same index length as the file.
    obsm, obsp, uns : mappings, optional
        Arrays / sparse matrices / arbitrary values written under the
        matching top-level group.
    logger : logging.Logger, optional
        Progress logger (no-op when None).

    Raises
    ------
    FileNotFoundError
        If ``h5ad_path`` does not exist.
    ValueError
        If ``obs`` length does not match the file's ``n_obs``.
    RuntimeError
        If the read-back verification fails.
    """
    import anndata as ad
    import h5py
    from anndata.io import write_elem

    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"Cannot append to non-existent h5ad: {h5ad_path}")

    ad.settings.allow_write_nullable_strings = True
    _log = logger.info if logger else (lambda msg, *a: None)
    for label, mapping in (("obs", obs), ("obsm", obsm), ("obsp", obsp), ("uns", uns)):
        if mapping is not None and len(mapping):
            _log(
                "Incremental %s: %d key(s) → %s", label, len(mapping), os.path.basename(h5ad_path)
            )

    compression = _detect_x_compression(h5ad_path)
    dkw = _write_kwargs(compression)

    with h5py.File(h5ad_path, "a") as f:
        if obs is not None and len(obs.columns):
            if "X" in f:
                # Sparse X is a Group whose shape lives in attrs; Dataset has .shape
                x_node = f["X"]
                if isinstance(x_node, h5py.Dataset):
                    n_x_obs = x_node.shape[0]
                elif "shape" in x_node.attrs:
                    n_x_obs = x_node.attrs["shape"][0]
                else:
                    n_x_obs = None
                if n_x_obs is not None and len(obs.index) != n_x_obs:
                    raise ValueError(
                        f"Incremental obs length {len(obs.index)} does not match "
                        f"file n_obs {n_x_obs}"
                    )
            _write_obs_columns(f.require_group("obs"), obs, dkw)
        for group_name, mapping in (("obsm", obsm), ("obsp", obsp), ("uns", uns)):
            if not mapping:
                continue
            group = f.require_group(group_name)
            for k, v in mapping.items():
                write_elem(group, k, v, dataset_kwargs=dkw)

    verify_incremental_write(h5ad_path, obs=obs, obsm=obsm, obsp=obsp, uns=uns)
    _log("Incremental write verified OK → %s", os.path.basename(h5ad_path))
