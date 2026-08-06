#!/usr/bin/env python3
"""grn_tissue.py — Tissue-aware TF activity annotation via KB marker overlap

``core/grn_tissue.py`` enriches the TF activity matrix produced by GRN
inference (ULM enrichment) with knowledge-base (KB) marker overlap
statistics.  The single public function ``compute_tf_relevance()`` is
called by downstream steps (T6 export, T7 gating) to compute per-TF
overlap ratios and produce an activity matrix whose values are weighted
by tissue relevance.

Functions
---------
compute_tf_relevance
    Compute KB-overlap statistics per TF and return a weighted activity
    matrix plus a separate annotation table.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

# ── Module-level logger ────────────────────────────────────────────────
_log = logging.getLogger(__name__)


# ======================================================================
#  Public API
# ======================================================================


def compute_tf_relevance(
    activity_df: pd.DataFrame,
    net: pd.DataFrame,
    kb_markers: Optional[set[str]] = None,
    log: Optional[logging.Logger] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Annotate each TF in the activity matrix with its overlap against KB markers.

    For every TF column in *activity_df*, the regulon target genes from
    *net* are intersected with *kb_markers*.  Two outputs are returned:

    1. **annotated_activity_df** — same shape as *activity_df* (rows = cell
       types, cols = TFs) with each cell replaced by the **weighted score**
       ``abs(activity) * (1 + kb_overlap_ratio)``.  TFs with no matching
       regulon in *net* are weighted by the raw ``abs(activity)`` only.

    2. **tf_annotation_table** — one row per TF with:

       - ``tf``                 — TF symbol
       - ``n_targets``          — number of target genes in the regulon
       - ``kb_overlap_count``   — target genes also present in *kb_markers*
       - ``kb_overlap_ratio``   — ``kb_overlap_count / n_targets``

    Parameters
    ----------
    activity_df : pd.DataFrame
        TF activity matrix from ULM enrichment.  Rows are cell types /
        groups, columns are TF symbols.  Values are numeric activity
        scores (may be positive or negative).
    net : pd.DataFrame
        Regulon network with at least two columns: ``source`` (TF
        regulator) and ``target`` (regulated gene).  This is the same
        ``net`` used by ``decoupler.mt.ulm()``.
    kb_markers : set of str or None
        Knowledge-base marker gene symbols (UPPERCASE).  When ``None``
        or empty, all overlap ratios are 0.0 and the weighted score
        degrades to ``abs(activity)``.
    log : logging.Logger or None
        Optional logger.  If ``None`` the module-level logger is used.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        ``(annotated_activity_df, tf_annotation_table)`` — see description
        above.

    Examples
    --------
    >>> markers = {"RHO", "NRL", "GNAT1", "VSX2"}
    >>> act_w, tf_ann = compute_tf_relevance(estimates_df, net, markers)
    """
    logger = log or _log

    # ── Normalise kb_markers ──────────────────────────────────────────
    markers: set[str] | None = kb_markers
    if markers is not None and len(markers) == 0:
        markers = None  # treat empty set as "not provided"

    # ── Handle missing/empty net ───────────────────────────────────────
    has_net: bool = not net.empty
    if not has_net:
        logger.warning("net is empty. All overlap statistics will be zero.")

    # ── Collect all TFs present in the activity matrix ─────────────────
    tf_list: list[str] = list(activity_df.columns)
    if not tf_list:
        logger.warning("activity_df has no columns (no TFs).")
        empty_ann = pd.DataFrame(
            {
                "tf": pd.Series(dtype=str),
                "n_targets": pd.Series(dtype=int),
                "kb_overlap_count": pd.Series(dtype=int),
                "kb_overlap_ratio": pd.Series(dtype=float),
            }
        )
        return activity_df.copy(), empty_ann

    # ── Uppercase net targets (case normalisation) ─────────────────────
    tf_target_map: dict[str, list[str]] = {}
    if has_net:
        net_upper = net.copy()
        net_upper["target"] = net_upper["target"].astype(str).str.upper()

        for tf in tf_list:
            targets = net_upper.loc[net_upper["source"] == tf, "target"].tolist()
            tf_target_map[tf] = targets
    else:
        for tf in tf_list:
            tf_target_map[tf] = []

    # ── Compute per-TF overlap statistics ──────────────────────────────
    n_targets_list: list[int] = []
    kb_overlap_count_list: list[int] = []
    kb_overlap_ratio_list: list[float] = []

    for tf in tf_list:
        targets = tf_target_map[tf]
        n_t = len(targets)
        n_targets_list.append(n_t)

        if markers and n_t > 0:
            overlap = sum(1 for g in targets if g in markers)
            ratio = overlap / n_t
        else:
            overlap = 0
            ratio = 0.0
        kb_overlap_count_list.append(overlap)
        kb_overlap_ratio_list.append(ratio)

    # ── Build annotation table ─────────────────────────────────────────
    tf_annotation_table = pd.DataFrame(
        {
            "tf": tf_list,
            "n_targets": n_targets_list,
            "kb_overlap_count": kb_overlap_count_list,
            "kb_overlap_ratio": kb_overlap_ratio_list,
        }
    )

    # ── Build weighted activity matrix ─────────────────────────────────
    #  annotated_df[tf] = abs(activity_df[tf]) * (1 + kb_overlap_ratio[tf])
    ratio_map: dict[str, float] = dict(zip(tf_list, kb_overlap_ratio_list, strict=False))
    annotated_activity_df = activity_df.copy()
    for tf in tf_list:
        r = ratio_map.get(tf, 0.0)
        annotated_activity_df[tf] = annotated_activity_df[tf].abs() * (1.0 + r)

    # ── Log summary ────────────────────────────────────────────────────
    n_with_overlap = sum(1 for c in kb_overlap_count_list if c > 0)
    logger.debug(
        "compute_tf_relevance: %d/%d TFs have KB overlap (markers=%s). Net edges: %s.",
        n_with_overlap,
        len(tf_list),
        "provided" if markers else "None/empty",
        "provided" if has_net else "empty",
    )

    return annotated_activity_df, tf_annotation_table


# ======================================================================
#  Streaming pseudobulk (bypasses full anndata load)
# ======================================================================


def streaming_pseudobulk(
    h5ad_path: str,
    group_col: str = "cell_type",
    log: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Per-group mean-expression pseudobulk by streaming the h5ad matrix.

    Reads ``raw/X`` (fallback ``X``) chunk-wise with h5py and accumulates
    per-group sums via CSR construction, never materialising the full matrix.
    This avoids loading the dense HVG ``X`` that anndata would otherwise
    read in full.  Returns a log1p-transformed group x gene DataFrame.
    """
    logger = log or _log
    import h5py
    import numpy as np

    try:
        import hdf5plugin  # noqa: F401  (registers the zstd filter)
    except ImportError:
        pass

    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]

        # ---- group codes: categorical (codes+categories) or plain array ----
        col = group_col if group_col in obs else "leiden"
        if col != group_col:
            logger.warning("%s not in obs - using 'leiden'", group_col)
        if f"{col}/categories" in obs:
            codes = np.asarray(obs[f"{col}/codes"][...])
            cats = obs[f"{col}/categories"][...]
            categories = [str(c, "utf-8") if isinstance(c, bytes) else str(c) for c in cats]
        else:
            raw_col = obs[col][...]
            dec = [c.decode() if isinstance(c, bytes) else c for c in raw_col]
            codes, cats = pd.factorize(dec)
            categories = [str(c) for c in cats]

        # ---- expression source: raw/X sparse (preferred), else X ----
        if "raw" in f and "X" in f["raw"]:
            xgrp = f["raw/X"]
            var_index = f["raw/var/_index"]
            use_raw = True
        else:
            xgrp = f["X"]
            var_index = f["var/_index"]
            use_raw = False

        var_names = [v.decode() if isinstance(v, bytes) else v for v in var_index[...]]
        n_groups = len(categories)
        n_genes = len(var_names)

        if isinstance(xgrp, h5py.Group):  # sparse CSR storage
            data = xgrp["data"]
            indices = xgrp["indices"]
            indptr_arr = xgrp["indptr"][...]
            n_cells = int(indptr_arr.shape[0] - 1)
            dense = None
        else:  # dense dataset
            data = indices = None
            indptr_arr = None
            n_cells = int(xgrp.shape[0])
            dense = xgrp

        sums = np.zeros((n_groups, n_genes), dtype=np.float64)
        chunk = 20_000
        if data is not None:
            for start in range(0, n_cells, chunk):
                end = min(start + chunk, n_cells)
                p0, p1 = int(indptr_arr[start]), int(indptr_arr[end])
                d = data[p0:p1]
                idx = indices[p0:p1]
                seg = np.diff(indptr_arr[start : end + 1])
                rows = np.repeat(np.arange(end - start), seg)
                gr = codes[start:end][rows]
                lin = gr * n_genes + idx
                bc = np.bincount(lin, weights=d, minlength=n_groups * n_genes)
                sums += bc.reshape(n_groups, n_genes)
        else:
            onehot = np.arange(n_groups)
            for start in range(0, n_cells, chunk):
                end = min(start + chunk, n_cells)
                block = dense[start:end]
                gc = codes[start:end]
                sums += (gc[:, None] == onehot).T.astype(np.float64) @ block

        counts = np.bincount(codes, minlength=n_groups)
        pseudo = np.log1p(sums / np.maximum(counts, 1)[:, None])
        df = pd.DataFrame(pseudo, index=categories, columns=var_names)

    logger.info("Pseudobulk: %d cells -> %d groups", n_cells, n_groups)
    logger.info("  Pseudobulk matrix: %d x %d (raw=%s)", n_groups, n_genes, use_raw)
    return df


# ======================================================================
#  Decoupler network download cache
# ======================================================================


def patch_decoupler_cache(log: Optional[logging.Logger] = None) -> None:
    """Cache decoupler network downloads under ``~/.cache/fuxi/decoupler/``.

    decoupler's ``_download`` (Zenodo CollecTRI, HCOP ortholog maps, ...)
    has no disk cache and re-fetches on every process; this wrapper makes
    each resource a single-fetch.  Idempotent - safe to call repeatedly.
    """
    logger = log or _log
    import hashlib
    import importlib
    import io
    import os

    from core.utils import fuxi_cache_dir

    cache_dir = fuxi_cache_dir("decoupler")
    os.makedirs(cache_dir, exist_ok=True)
    import decoupler._download as _dl

    _orig = _dl._download
    if getattr(_dl, "_fuxi_cached", False):
        return

    def _cached(
        url: str, verbose: bool = False, retries: int = 5, wait_time: int = 20
    ) -> io.BytesIO:
        key = hashlib.sha1(url.encode()).hexdigest()[:16]
        fname = url.rsplit("/", 1)[-1].split("?", 1)[0]
        path = os.path.join(cache_dir, f"{key}_{fname}")
        if os.path.exists(path):
            with open(path, "rb") as fh:
                data = io.BytesIO(fh.read())
            logger.info("decoupler cache hit: %s", fname)
            return data
        data = _orig(url, verbose=verbose, retries=retries, wait_time=wait_time)
        with open(path, "wb") as fh:
            fh.write(data.getvalue())
        data.seek(0)
        logger.info("decoupler cached: %s", fname)
        return data

    # Module-level ``from decoupler._download import _download`` bindings
    _dl._download = _cached
    _dl._fuxi_cached = True
    for _mod_name in (
        "decoupler.op._collectri",
        "decoupler.op._translate",
        "decoupler.op._dorothea",
    ):
        try:
            _mod = importlib.import_module(_mod_name)
        except ImportError:
            continue
        if getattr(_mod, "_download", None) is _orig:
            _mod._download = _cached
        _mod._fuxi_cached = True
    logger.info("decoupler download cache enabled -> %s", cache_dir)
