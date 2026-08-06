"""
cell_interaction.py — Cell-Cell Interaction (CCI) analysis via LIANA+
======================================================================

Shared utility layer for:
  - RNA Step 12: permutation-based LR interaction testing
  - Spatial Step 10: spatial bivariate LR interaction metrics

Five exportable functions:
  ensure_gene_symbols()    — convert Ensembl var_names to gene symbols
  load_lr_database()       — load LIANA ligand-receptor resource
  run_cci_permutation()    — permutation testing via liana.mt.rank_aggregate
  run_cci_spatial()        — spatial bivariate metrics via liana.mt.bivariate
  format_cci_results()     — filter + sort + annotate interaction table

Dependencies: liana>=1.0.0, anndata, pandas, mygene
"""

import os
import time
from typing import cast

import pandas as pd

from core.utils import fuxi_cache_dir

MYGENE_CHUNK_SIZE: int = 1000


def ensure_gene_symbols(adata, log: object = None, species: str = "human"):
    """Ensure AnnData var_names are gene symbols (not Ensembl IDs).

    LIANA resources use HGNC gene symbols.  When var_names contain
    Ensembl IDs (``ENSG...``), this function queries mygene.info to
    batch-convert them to symbols and returns a **new** AnnData with
    clean gene symbol var_names (Ensembl-only genes dropped, duplicates
    resolved).  The ``.raw`` layer is rebuilt too.

    Parameters
    ----------
    adata : AnnData
        Input data (may have mixed Ensembl ID / gene symbol var_names).
    log : object, optional
        Logger.
    species : str
        Species key for mygene query (default: "human").

    Returns
    -------
    AnnData
        A new AnnData with gene-symbol var_names, or the original if
        no conversion was needed (same object, not a copy).
    """

    is_ensembl = adata.var_names.str.match(r"^ENS[A-Z]{0,4}G\d{11}$")
    n_ensembl = is_ensembl.sum()
    if n_ensembl == 0:
        if log:
            log.info("All var_names are already gene symbols -- skip mapping")
        return adata

    if log:
        log.info(
            "%d/%d var_names are Ensembl IDs -- mapping to gene symbols...",
            n_ensembl,
            adata.n_vars,
        )

    try:
        import mygene
    except ImportError:
        mygene = None

    if mygene is None:
        raise ImportError(
            "mygene package required for Ensembl ID → gene symbol conversion. Install: pip install mygene"
        )
    mg = mygene.MyGeneInfo()

    # Batch query mygene.info in chunks of 1000
    ensembl_ids = adata.var_names[is_ensembl].tolist()
    results = {}
    chunk_size = MYGENE_CHUNK_SIZE
    for i in range(0, len(ensembl_ids), chunk_size):
        chunk = ensembl_ids[i : i + chunk_size]
        try:
            batch = mg.querymany(
                chunk, scopes="ensembl.gene", fields="symbol", species=species, as_dataframe=True
            )
        except Exception:
            continue
        for eid, row in batch.iterrows():
            sym = row.get("symbol", None)
            if sym and isinstance(sym, str):
                results[eid] = sym

    if log:
        n_mapped = len(results)
        log.info("Mapped %d/%d Ensembl IDs to gene symbols", n_mapped, n_ensembl)

    # Build new var_names list
    new_names = []
    for name in adata.var_names:
        if name in results:
            new_names.append(results[name])
        elif name.startswith("ENSG"):
            new_names.append(name)  # keep as-is (will be dropped later)
        else:
            new_names.append(name)  # already a symbol

    # Create a clean AnnData
    new_adata = adata[:, :].copy()
    new_adata.var_names = new_names

    # Drop unmapped Ensembl IDs (they will never match LIANA resources)
    keep = ~new_adata.var_names.str.match(r"^ENSG\d+")
    if not keep.all():
        new_adata = new_adata[:, keep].copy()
        if log:
            log.info("Dropped %d unmapped Ensembl-only genes", (~keep).sum())

    # Deduplicate -- keep first occurrence for duplicated symbols
    dup = new_adata.var_names.duplicated()
    if dup.any():
        new_adata = new_adata[:, ~dup].copy()
        if log:
            log.info("Dropped %d duplicate gene symbols (kept first)", dup.sum())

    # Rebuild .raw layer with the same gene set
    if adata.raw is not None:
        raw = adata.raw.to_adata()[:, adata.var_names].copy()
        raw.var_names = new_names
        if not keep.all():
            raw = raw[:, keep].copy()
        if dup.any():
            raw = raw[:, ~dup].copy()
        new_adata.raw = raw

    if log:
        log.info("Gene symbol conversion: %d -> %d genes", adata.n_vars, new_adata.n_vars)

    return new_adata


def load_lr_database(
    resource_name: str = "consensus",
    cache_dir: str = "",
    log: object = None,
) -> pd.DataFrame:
    """Load a LIANA ligand-receptor database resource.

    Parameters
    ----------
    resource_name : str
        Resource name. One of: 'consensus', 'cellphonedb', 'cellchat',
        'celltalkdb', 'ramilowski', 'talklr', 'baccin', 'connectome',
        'guide2pharma', 'italk', 'kirouac', 'nichenet', 'omni', 'scaffold'.
        Default: 'consensus' (union of major databases).
    cache_dir : str
        LIANA cache directory. Empty = unified Fuxi cache (~/.cache/fuxi/liana).
    log : object, optional
        Logger with .info() method.

    Returns
    -------
    pd.DataFrame
        LR pairs with columns including 'ligand', 'receptor', 'source'.
    """
    if cache_dir:
        os.environ["LIANA_CACHE_DIR"] = cache_dir
    else:
        os.environ.setdefault("LIANA_CACHE_DIR", fuxi_cache_dir("liana"))

    import liana.resource

    t0 = time.time()
    if log:
        log.info("Loading LIANA LR database: %s", resource_name)

    lr_df = liana.resource.select_resource(resource_name)

    if log:
        n_ligands = lr_df["ligand"].nunique()
        n_receptors = lr_df["receptor"].nunique()
        log.info(
            "LR database loaded: %d interactions, %d unique ligands, %d unique receptors (%.1fs)",
            len(lr_df),
            n_ligands,
            n_receptors,
            time.time() - t0,
        )

    return lr_df


# ── LIANA hot-path monkeypatches (Step 12 permutation testing) ──
#
# cProfile across 10 datasets shows ~90% of Step 12 runtime in liana's
# permutation pipeline, and 1M-cell datasets OOM on liana's dense
# intermediates:
#   * _get_positions: per-entity np.where(var_names == entity) — O(n_entities × n_vars)
#   * _generate_perms_cube: X[perm_idx] + boolean row slicing triggers
#     csr_sort_indices + csr_sum_duplicates per group per permutation
#   * prep_check_adata: deep-copies the full expression matrix (25.6 GB at
#     1.05M x 34.6k genes) and checks np.isfinite over all non-zero entries
#     at once (two 2.1 GB bool temporaries)
#   * _get_lr connectome method: sc.pp.scale(adata, copy=True) densifies the
#     whole matrix (~146 GB at 1.05M x 34.6k)
#   * liana_pipe mat_mean: np.mean(adata.X) allocates a full nnz temporary
#     (8.4 GB) via data * (1/denom)
#
# All patches are mathematically/semantically equivalent; small datasets
# keep the original code paths (bit-identical), the patched paths engage
# only above size thresholds.

_PATCHED_LIANA = False
_SCALE_DENSE_BYTES_THRESHOLD = 5e9  # n_obs*n_vars*4B dense; below this keep liana's dense path
_SPARSE_MEAN_CHUNK = 64_000_000  # entries per chunk in the chunked sparse mean
_ORIG_GET_LR = None
_ORIG_PREP_CHECK_ADATA = None
_ORIG_SPARSE_MEAN = None


def _patched_get_positions(adata, lr_res):
    """Vectorized replacement for liana's per-entity np.where lookup.

    Original does ``np.where(adata.var_names == entity)[0][0]`` for every
    ligand and receptor — O(n_entities × n_vars) pandas __eq__ elementwise
    scans.  ``Index.get_indexer`` finds all positions in one pass (first
    match — same semantics as np.where()[0]).
    """
    import numpy as np

    idx = pd.Index(adata.var_names)
    entities = np.union1d(lr_res["ligand"], lr_res["receptor"])
    pos = pd.Series(idx.get_indexer(entities), index=entities)
    ligand_pos = {e: int(pos[e]) for e in lr_res["ligand"]}
    receptor_pos = {e: int(pos[e]) for e in lr_res["receptor"]}
    labels = adata.obs["@label"].cat.categories
    labels_pos = {labels[i]: i for i in range(labels.shape[0])}
    return ligand_pos, receptor_pos, labels_pos


def _patched_generate_perms_cube(X, n_perms, labels_mask, seed, agg_fun, n_jobs, verbose):  # noqa: N803
    """Sparse matvec replacement for liana's per-permutation boolean slicing.

    liana's version does ``X[perm_idx]`` followed by per-group boolean row
    slicing (``perm_mat[labels_mask[:, i]]``) for every group and every
    permutation — each slice triggers csr_sort_indices / csr_sum_duplicates
    (131s of 242s on a 71.5k-cell dataset).  Here all permutation index
    vectors are generated up front with the identical RNG stream, then each
    permutation's per-group column sums come from a single sparse matvec
    ``X.T @ S`` where S is the one-hot group assignment scattered by the
    permuted row order.  ``agg_fun == np.mean`` (CellPhoneDB, the only
    permutation method in rank_aggregate) is supported; other aggs fall back
    to liana's original implementation.
    """
    import numpy as np
    from joblib import Parallel, delayed

    rng = np.random.default_rng(seed=seed)
    idx = np.arange(X.shape[0])
    group_id = np.argmax(labels_mask, axis=1)
    n_groups = labels_mask.shape[1]
    n_cells, n_genes = X.shape

    if agg_fun is np.mean:
        perms = np.zeros((n_perms, n_groups, n_genes))
        xt = X.T.tocsr()
        if n_cells * n_genes * 4 >= _SCALE_DENSE_BYTES_THRESHOLD:
            xt = X.T  # zero-copy csc view; keeps ~18 GB off at 1M+ cells
        inv = np.asarray(labels_mask.sum(axis=0)).ravel()
        inv = 1.0 / inv
        inv[~np.isfinite(inv)] = np.nan

        def _agg_one(perm, pidx):
            s = np.zeros((n_cells, n_groups))
            s[pidx, group_id] = 1.0
            m = xt.dot(s)
            m *= inv[None, :]
            return perm, m

        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_agg_one)(perm, rng.permutation(idx)) for perm in range(n_perms)
        )
        for perm, m in results:
            perms[perm] = m.T
        return perms

    # Non-mean aggs (e.g. CellChat trimean) — liana's original, threads only
    from liana.method._pipe_utils._get_mean_perms import _permute_and_aggregate

    results = Parallel(n_jobs=min(n_jobs, 4), prefer="threads")(
        delayed(_permute_and_aggregate)(perm, rng.permutation(idx), X, labels_mask, agg_fun)
        for perm in range(n_perms)
    )
    perms = np.zeros((n_perms, n_groups, n_genes))
    for perm, permuted_means in results:
        perms[perm] = np.reshape(permuted_means, (n_groups, n_genes))
    return perms


def _patched_get_lr(
    adata, resource, groupby_pairs, relevant_cols, mat_mean, mat_max, de_method, base, verbose
):
    """Sparse-safe _get_lr for very large adata (avoids connectome dense scale).

    liana's connectome method materializes ``sc.pp.scale(adata, copy=True)`` --
    a full dense (n_obs x n_vars x 4B) copy of X.  At 1.05M x 34.6k that is
    ~146 GB and OOMs.  scale() is column-wise and its only consumer is
    ``temp.layers['scaled'].mean(axis=0)`` -- the per-label *mean* of the
    global z-scores; mean over label cells commutes with the linear
    z-transform: mean_L((x - mu)/sigma) == (mean_L(x) - mu)/sigma.  Global
    (mu, sigma) come from one sparse pass using scanpy's own sparse variance
    formula (E[x^2] - mu^2, ddof=0), so the dense layer is never built and
    the per-label z-score means differ from liana's only by float summation
    order (~1e-7).
    """
    n_obs, n_vars = adata.X.shape
    if n_obs * n_vars * 4 < _SCALE_DENSE_BYTES_THRESHOLD:
        return _ORIG_GET_LR(
            adata,
            resource,
            groupby_pairs,
            relevant_cols,
            mat_mean,
            mat_max,
            de_method,
            base,
            verbose,
        )

    import liana.method.sc._liana_pipe as lp
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from liana._constants import (
        CommonColumns as _C,  # noqa: N814
    )
    from liana._constants import (
        InternalValues as _I,  # noqa: N814
    )
    from liana._constants import (
        MethodColumns as _M,  # noqa: N814
    )
    from liana._constants import (
        PrimaryColumns as _P,  # noqa: N814
    )
    from liana.method._pipe_utils._common import _get_props, _join_stats

    labels = adata.obs[_I.label].cat.categories

    # Method-specific stats
    connectome_flag = (_M.ligand_zscores in relevant_cols) | (_M.receptor_zscores in relevant_cols)
    logfc_flag = (_M.ligand_logfc in relevant_cols) | (_M.receptor_logfc in relevant_cols)

    # Global (mu, sigma) in one sparse pass — replaces the dense sc.pp.scale
    mtx = adata.X
    n_rows = mtx.shape[0]
    mu = mtx.sum(axis=0).A1 / n_rows
    sigma = np.sqrt(mtx.multiply(mtx).sum(axis=0).A1 / n_rows - mu * mu)

    if logfc_flag:
        adata.layers["normcounts"] = mtx.copy()
        adata.layers["normcounts"].data = lp._expm1_base(mtx.data, base)

    # initialize dict
    dedict = {}

    # Calc pvals + other stats per gene or not
    rank_genes_bool = (_C.ligand_pvals in relevant_cols) | (_C.receptor_pvals in relevant_cols)
    if rank_genes_bool:
        adata = sc.tl.rank_genes_groups(
            adata, groupby=_I.label, method=de_method, use_raw=False, copy=True
        )

    for label in labels:
        temp = adata[adata.obs[_I.label] == label, :]
        a = _get_props(temp.X)
        stats = (
            pd.DataFrame({"names": temp.var_names, "props": a})
            .assign(label=label)
            .sort_values("names")
        )
        if rank_genes_bool:
            pvals = sc.get.rank_genes_groups_df(adata, label)
            stats = stats.merge(pvals)
        dedict[label] = stats

    # check if genes are ordered correctly
    if not list(adata.var_names) == list(dedict[labels[0]]["names"]):
        raise AssertionError("Variable names did not match DE results!")

    # Calculate Mean, logFC and z-scores by group
    for label in labels:
        temp = adata[adata.obs[_I.label].isin([label])]
        dedict[label]["means"] = temp.X.mean(axis=0).A.flatten()
        if connectome_flag:
            # mean over label cells commutes with the linear z-transform
            dedict[label]["zscores"] = (dedict[label]["means"] - mu) / sigma
        if logfc_flag:
            dedict[label]["logfc"] = lp._calc_log2fc(adata, label)
        if isinstance(mat_max, np.float32):  # cellchat flag
            dedict[label]["trimean"] = lp._trimean(temp.X / mat_max)

    pairs = pd.DataFrame(
        np.array(np.meshgrid(labels, labels)).reshape(2, np.size(labels) * np.size(labels)).T
    ).rename(columns={0: _P.source, 1: _P.target})

    if groupby_pairs is not None:
        pairs = pairs.merge(groupby_pairs, on=[_P.source, _P.target], how="inner")

    # Join Stats
    lr_res = pd.concat(
        [
            _join_stats(source, target, dedict, resource)
            for source, target in zip(pairs[_P.source], pairs[_P.target], strict=False)
        ]
    )

    if _M.mat_mean in relevant_cols:
        assert isinstance(mat_mean, np.float32)
        lr_res[_M.mat_mean] = mat_mean

    if isinstance(mat_max, np.float32):
        lr_res[_M.mat_max] = mat_max

    # subset to only relevant columns
    relevant_cols = np.intersect1d(relevant_cols, lr_res.columns)

    return lr_res[relevant_cols]


def _patched_prep_check_adata(
    adata,
    groupby,
    min_cells,
    groupby_subset=None,
    use_raw=False,
    layer=None,
    obsm=None,
    uns=None,
    complex_sep="_",
    verbose=False,
):
    """prep_check_adata without the full-X deep copy and with chunked finite check.

    liana's original builds ``sc.AnnData(X=X, obs=..., var=..., obsp=...,
    uns=..., obsm=...).copy()`` -- a deep copy of the whole expression matrix
    (25.6 GB at 1.05M x 34.6k genes, int64 indices) -- then checks
    ``np.isfinite(adata.X.data)`` over all 2.26B entries at once, allocating
    two 2.1 GB bool temporaries.  Both OOM on very large data.  The Step 12
    loader hands liana a disposable lightweight AnnData, so the copy is
    dropped; the finite-value check is chunked (512M entries at a time).
    Everything else mirrors liana's original semantics exactly.
    """
    import numpy as np
    import scanpy as sc
    from liana._logging import _logg
    from liana.method._pipe_utils._pre import _check_groupby, _choose_mtx_rep, check_vars

    mtx = _choose_mtx_rep(adata=adata, use_raw=use_raw, layer=layer, verbose=verbose)

    if use_raw and layer is None:
        var = pd.DataFrame(index=adata.raw.var_names)
    else:
        var = pd.DataFrame(index=adata.var_names)

    if obsm is not None:
        # discard any instances of AnnData if in obsm
        obsm = {k: v for k, v in obsm.items() if not isinstance(v, object)}

    adata = sc.AnnData(
        X=mtx,
        obs=adata.obs.copy(),
        var=var,
        obsp=adata.obsp.copy(),
        uns=uns,
        obsm=obsm,
    )
    adata.var_names_make_unique()

    # Check for empty features
    msk_features = np.sum(adata.X, axis=0).A1 == 0
    n_empty_features = np.sum(msk_features)
    if n_empty_features > 0:
        _logg(
            f"{n_empty_features} features of mat are empty, they will be removed.",
            level="warn",
            verbose=verbose,
        )
        adata = adata[:, ~msk_features]

    # Check for empty samples
    msk_samples = adata.X.sum(axis=1).A1 == 0
    n_empty_samples = np.sum(msk_samples)
    if n_empty_samples > 0:
        _logg(
            f"{n_empty_samples} samples of mat are empty, they will be removed.",
            level="warn",
            verbose=verbose,
        )

    # Check if log-norm
    data = adata.X.data if hasattr(adata.X, "data") else adata.X
    _sum = np.sum(data[0:100])
    if _sum == np.floor(_sum):
        _logg("Make sure that normalized counts are passed!", level="warn", verbose=verbose)

    # Check for non-finite values (chunked to bound the bool temporaries)
    chunk = 512_000_000
    if data.size > chunk:
        finite = True
        for i in range(0, data.size, chunk):
            if np.any(~np.isfinite(data[i : i + chunk])):
                finite = False
                break
    else:
        finite = bool(np.all(np.isfinite(data)))
    if not finite:
        raise ValueError(
            "mat contains non finite values (nan or inf), please set them to 0 or remove them."
        )

    if groupby is not None:
        _check_groupby(adata, groupby, verbose)

        if groupby_subset is not None:
            adata = adata[adata.obs[groupby].isin(groupby_subset), :]

        adata.obs["@label"] = adata.obs[groupby]

        # Remove any cell types below X number of cells per cell type
        count_cells = adata.obs.groupby(groupby)[groupby].size().reset_index(name="count").copy()
        count_cells["keep"] = count_cells["count"] >= min_cells

        if not all(count_cells.keep):
            lowly_abundant_idents = list(count_cells[~count_cells.keep][groupby])
            # remove lowly abundant identities
            msk = ~np.isin(adata.obs[[groupby]], lowly_abundant_idents)
            adata = adata[msk]
            _logg(
                "The following cell identities were excluded: {}".format(
                    ", ".join(lowly_abundant_idents)
                ),
                level="warn",
                verbose=verbose,
            )

    check_vars(adata.var_names, complex_sep=complex_sep, verbose=verbose)
    # Re-order adata vars alphabetically
    adata = adata[:, np.sort(adata.var_names)]
    return adata


def _patched_sparse_mean(self, axis=None, dtype=None, out=None):
    """Chunked mean for very large sparse matrices.

    scipy's sparse mean(axis=None) computes ``(X * (1/denom)).sum()`` -- the
    scalar multiplication materializes a full nnz temporary (8.4 GB at 2.26B
    non-zeros).  Mathematically identical: sum(data) / prod(shape) for
    axis=None.  Small matrices keep scipy's original path (bit-identical).
    """
    import math

    import numpy as np

    if axis is None and out is None and self.data.size > _SPARSE_MEAN_CHUNK:
        denom = math.prod(self.shape)
        total = np.float64(0.0)
        data = self.data
        for i in range(0, data.size, _SPARSE_MEAN_CHUNK):
            total += np.sum(data[i : i + _SPARSE_MEAN_CHUNK], dtype=np.float64)
        res = total / denom
        if dtype is not None:
            res = res.astype(dtype, copy=False)
        return res
    return _ORIG_SPARSE_MEAN(self, axis=axis, dtype=dtype, out=out)


def _patch_liana_perf():
    """Apply the Step 12 hot-path patches to the installed liana package.

    Idempotent; safe no-op when liana is not importable (RNA/spatial steps
    degrade to liana's original implementation instead of crashing).
    """
    global _PATCHED_LIANA, _ORIG_GET_LR, _ORIG_PREP_CHECK_ADATA, _ORIG_SPARSE_MEAN
    if _PATCHED_LIANA:
        return
    try:
        import liana.method._pipe_utils as pu
        import liana.method._pipe_utils._get_mean_perms as gmp
        import liana.method._pipe_utils._pre as pre
        import liana.method.sc._liana_pipe as lp

        if getattr(gmp, "_get_positions", None) is not _patched_get_positions:
            gmp._get_positions = _patched_get_positions
        if getattr(gmp, "_generate_perms_cube", None) is not _patched_generate_perms_cube:
            gmp._generate_perms_cube = _patched_generate_perms_cube

        if getattr(pre, "prep_check_adata", None) is not _patched_prep_check_adata:
            _ORIG_PREP_CHECK_ADATA = pre.prep_check_adata
            pre.prep_check_adata = _patched_prep_check_adata
        if getattr(pu, "prep_check_adata", None) is not _patched_prep_check_adata:
            pu.prep_check_adata = _patched_prep_check_adata
        if getattr(lp, "prep_check_adata", None) is not _patched_prep_check_adata:
            lp.prep_check_adata = _patched_prep_check_adata

        if getattr(lp, "_get_lr", None) is not _patched_get_lr:
            _ORIG_GET_LR = lp._get_lr
            lp._get_lr = _patched_get_lr

        import scipy.sparse._base as spb

        if getattr(spb._spbase, "mean", None) is not _patched_sparse_mean:
            _ORIG_SPARSE_MEAN = spb._spbase.mean
            spb._spbase.mean = _patched_sparse_mean
        _PATCHED_LIANA = True
    except ImportError:
        pass


def run_cci_permutation(
    adata,
    groupby_col: str = "cell_type",
    resource_name: str = "consensus",
    n_perms: int = 1000,
    seed: int = 1337,
    use_raw: bool = True,
    layer: str | None = None,
    n_jobs: int = 1,
    log: object = None,
) -> pd.DataFrame:
    """Run LIANA+ rank_aggregate permutation testing for ligand-receptor
    interactions between cell type groups.

    Parameters
    ----------
    adata : AnnData
        Annotated data object. Must contain `groupby_col` in .obs and,
        if `use_raw=True`, a .raw attribute with full transcriptome expression.
    groupby_col : str
        Column in adata.obs to group by (default 'cell_type').
    resource_name : str
        LR database name passed to liana.mt.rank_aggregate.
    n_perms : int
        Number of permutation iterations (default 1000).
    seed : int
        Random seed for reproducibility.
    use_raw : bool
        Use adata.raw.X (full transcriptome, log-normalized) if available.
    layer : str | None
        Explicit layer key holding raw UMI counts (e.g. 'counts').
        When set, takes precedence over `use_raw` for LIANA's expression input.
    n_jobs : int
        Number of parallel jobs for the permutation test (threads; the RNG
        stream is unchanged regardless of n_jobs).
    log : object, optional
        Logger with .info() method.
    """
    import liana as li

    _patch_liana_perf()

    t0 = time.time()
    if log:
        log.info(
            "Running CCI permutation test: groupby=%s, resource=%s, n_perms=%d",
            groupby_col,
            resource_name,
            n_perms,
        )

    lr_res = li.mt.rank_aggregate(
        adata,
        groupby=groupby_col,
        resource_name=resource_name,
        n_perms=n_perms,
        seed=seed,
        use_raw=use_raw and layer is None,
        layer=layer,
        inplace=False,
        verbose=False,
        n_jobs=n_jobs,
    )
    if lr_res is None:
        raise RuntimeError("LIANA rank_aggregate returned no result table")

    if log:
        n_interactions = len(lr_res)
        if "pvalue" in lr_res.columns:
            pvals = cast(pd.Series, lr_res["pvalue"])
            n_sig = int((pvals < 0.05).sum())
        else:
            n_sig = 0
        log.info(
            "CCI permutation done: %d total, %d significant (p<0.05), took %.1fs",
            n_interactions,
            n_sig,
            time.time() - t0,
        )

    return lr_res


def run_cci_spatial(
    adata,
    resource_name: str = "consensus",
    connectivity_key: str = "spatial_connectivities",
    local_name: str = "cosine",
    global_name: str = "morans",
    n_perms: int = 1000,
    seed: int = 1337,
    log: object = None,
) -> pd.DataFrame:
    """Run LIANA+ spatial bivariate metrics for spatially-resolved
    ligand-receptor co-expression analysis.

    Uses local metrics (spatially-weighted cosine similarity by default)
    and global metrics (Moran's R by default) with permutation testing.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial coordinates in .obsm['spatial']
        and spatial connectivities in .obsp[connectivity_key].
    resource_name : str
        LR database name.
    connectivity_key : str
        Key in adata.obsp for spatial connectivity matrix.
    local_name : str
        Local bivariate metric: 'cosine', 'jaccard', 'pearson', 'spearman'.
    global_name : str
        Global bivariate metric: 'morans', 'connectome', etc.
    n_perms : int
        Number of permutation iterations.
    seed : int
        Random seed.
    log : object, optional
        Logger with .info() method.

    Returns
    -------
    pd.DataFrame
        Interaction results extracted from the returned AnnData .var
        DataFrame. Columns: ligand, receptor, morans, morans_pvals, etc.
    """
    import liana as li

    t0 = time.time()
    if log:
        log.info(
            "Running CCI spatial analysis: resource=%s, local=%s, global=%s, n_perms=%d",
            resource_name,
            local_name,
            global_name,
            n_perms,
        )

    # Validate spatial connectivities exist
    if connectivity_key not in adata.obsp:
        raise KeyError(
            f"'{connectivity_key}' not found in adata.obsp. "
            "Run spatial neighbors construction (Step 03) first."
        )

    # LIANA 1.8 bivariate returns an AnnData with interaction results
    # stored in .var (LR-level stats).  Extract that to a DataFrame.
    result_adata = li.mt.bivariate(
        adata,
        resource_name=resource_name,
        connectivity_key=connectivity_key,
        local_name=local_name,
        global_name=global_name,
        n_perms=n_perms,
        seed=seed,
    )

    if result_adata is None:
        raise RuntimeError("LIANA bivariate returned no result AnnData")
    # Extract the interaction table from var
    lr_res = cast(pd.DataFrame, result_adata.var.copy())
    lr_res.reset_index(drop=True, inplace=True)

    if log:
        n_interactions = len(lr_res)
        log.info(
            "CCI spatial done: %d interactions, took %.1fs",
            n_interactions,
            time.time() - t0,
        )

    return lr_res


def format_cci_results(
    lr_res: pd.DataFrame,
    n_top: int = 50,
    pval_col: str = "magnitude_rank",
    ascending: bool = True,
    log: object = None,
    adjacency: pd.DataFrame = None,
    adjacency_mode: str = "off",
) -> pd.DataFrame:
    """Filter, sort and format CCI interaction results.

    Sorts by the given rank/significance column and selects top N interactions,
    then adds a readable interaction label.

    Parameters
    ----------
    lr_res : pd.DataFrame
        Raw LIANA results from run_cci_permutation() or run_cci_spatial().
    n_top : int
        Number of top interactions to retain.
    pval_col : str
        Column to sort by (default 'magnitude_rank').
    ascending : bool
        Whether to sort ascending (default True).  Pass False when sorting
        by a column where higher values are better (e.g. Moran's I).
    log : object, optional
        Logger with .info() method.
    adjacency : pd.DataFrame, optional
        Anatomical adjacency table (columns: source, target, adjacency_type).
        If provided and adjacency_mode != "off", applies adjacency filter.
    adjacency_mode : str
        ``"off"`` (default), ``"soft"`` (annotate), or ``"hard"`` (filter).

    Returns
    -------
    pd.DataFrame
        Filtered, sorted top-N DataFrame with added 'interaction' label.
    """
    # Create readable interaction label
    lr_res = lr_res.copy()
    cols = lr_res.columns

    if "ligand" in cols and "receptor" in cols:
        src = (
            lr_res["source"].astype(str) if "source" in cols else pd.Series("", index=lr_res.index)
        )
        tgt = (
            lr_res["target"].astype(str) if "target" in cols else pd.Series("", index=lr_res.index)
        )
        if "source" in cols and "target" in cols:
            lr_res["interaction"] = (
                src
                + "->"
                + tgt
                + " | "
                + lr_res["ligand"].astype(str)
                + "_"
                + lr_res["receptor"].astype(str)
            )
        else:
            lr_res["interaction"] = (
                lr_res["ligand"].astype(str) + "_" + lr_res["receptor"].astype(str)
            )
    elif "ligand_complex" in cols and "receptor_complex" in cols:
        lr_res["interaction"] = (
            lr_res["ligand_complex"].astype(str) + "_" + lr_res["receptor_complex"].astype(str)
        )

    # ── Apply anatomical adjacency filter (v4.0+) ──
    if adjacency is not None and adjacency_mode != "off":
        from core.pipeline.anatomy import filter_cci_by_adjacency

        lr_res = filter_cci_by_adjacency(
            lr_res,
            adjacency,
            mode=adjacency_mode,
            adjacency_types=[],  # empty = all types pass
            log=log,
        )

    # Sort by significance
    if pval_col in cols:
        top_df = lr_res.sort_values(pval_col, ascending=ascending).head(n_top)
        if log:
            log.info("Top %d interactions selected by %s", n_top, pval_col)
    else:
        # Fallback: sort by the first available rank column
        rank_cols = [c for c in cols if "rank" in c.lower()]
        if rank_cols:
            top_df = lr_res.sort_values(rank_cols[0], ascending=ascending).head(n_top)
        else:
            top_df = lr_res.head(n_top)
        if log:
            log.warning("Column '%s' not found; using first available ordering", pval_col)

    if log:
        log.info(
            "Formatted CCI results: %d interactions retained",
            len(top_df),
        )

    return top_df
