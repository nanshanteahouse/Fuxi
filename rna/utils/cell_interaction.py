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
import pandas as pd
from typing import Optional


def ensure_gene_symbols(adata, log: object = None):
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

    Returns
    -------
    AnnData
        A new AnnData with gene-symbol var_names, or the original if
        no conversion was needed (same object, not a copy).
    """
    import numpy as np

    is_ensembl = adata.var_names.str.match(r"^ENSG\d+")
    n_ensembl = is_ensembl.sum()
    if n_ensembl == 0:
        if log:
            log.info("All var_names are already gene symbols -- skip mapping")
        return adata

    if log:
        log.info("%d/%d var_names are Ensembl IDs -- mapping to gene symbols...",
                 n_ensembl, adata.n_vars)

    import mygene
    mg = mygene.MyGeneInfo()

    # Batch query mygene.info in chunks of 1000
    ensembl_ids = adata.var_names[is_ensembl].tolist()
    results = {}
    chunk_size = 1000
    for i in range(0, len(ensembl_ids), chunk_size):
        chunk = ensembl_ids[i:i + chunk_size]
        try:
            batch = mg.querymany(chunk, scopes="ensembl.gene", fields="symbol",
                                 species="human", as_dataframe=True)
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
            new_names.append(name)   # keep as-is (will be dropped later)
        else:
            new_names.append(name)   # already a symbol

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
        raw = adata.raw.to_adata()[:, :].copy()
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
        LIANA cache directory. Empty = auto (~/.cache/liana).
    log : object, optional
        Logger with .info() method.

    Returns
    -------
    pd.DataFrame
        LR pairs with columns including 'ligand', 'receptor', 'source'.
    """
    if cache_dir:
        os.environ["LIANA_CACHE_DIR"] = cache_dir

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
            len(lr_df), n_ligands, n_receptors, time.time() - t0,
        )

    return lr_df


def run_cci_permutation(
    adata,
    groupby_col: str = "cell_type",
    resource_name: str = "consensus",
    n_perms: int = 1000,
    seed: int = 1337,
    use_raw: bool = True,
    n_jobs: int = 1,
    log: object = None,
) -> pd.DataFrame:
    """Run LIANA+ rank_aggregate permutation testing for ligand-receptor
    interactions between cell type groups.

    Parameters
    ----------
    adata : AnnData
        Annotated data object. Must contain `groupby_col` in .obs and,
        if `use_raw=True`, a .raw attribute with raw counts.
    groupby_col : str
        Column in adata.obs to group by (default 'cell_type').
    resource_name : str
        LR database name passed to liana.mt.rank_aggregate.
    n_perms : int
        Number of permutation iterations (default 1000).
    seed : int
        Random seed for reproducibility.
    use_raw : bool
        Use adata.raw.X (raw counts) if available.
    n_jobs : int
        Number of parallel jobs for the permutation test.
    log : object, optional
        Logger with .info() method.

    Returns
    -------
    pd.DataFrame
        LIANA results with columns: source, target, ligand, receptor,
        magnitude_rank, specificity_rank, pvalue, etc.
    """
    import liana as li

    t0 = time.time()
    if log:
        log.info(
            "Running CCI permutation test: groupby=%s, resource=%s, n_perms=%d",
            groupby_col, resource_name, n_perms,
        )

    lr_res = li.mt.rank_aggregate(
        adata,
        groupby=groupby_col,
        resource_name=resource_name,
        n_perms=n_perms,
        seed=seed,
        use_raw=use_raw,
        n_jobs=n_jobs,
        inplace=False,
        verbose=False,
    )

    if log:
        n_interactions = len(lr_res)
        n_sig = (lr_res.get("pvalue", 1.0) < 0.05).sum() if "pvalue" in lr_res.columns else 0
        log.info(
            "CCI permutation done: %d total, %d significant (p<0.05), took %.1fs",
            n_interactions, n_sig, time.time() - t0,
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
            resource_name, local_name, global_name, n_perms,
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

    # Extract the interaction table from var
    lr_res = result_adata.var.copy()
    lr_res.reset_index(drop=True, inplace=True)

    if log:
        n_interactions = len(lr_res)
        log.info(
            "CCI spatial done: %d interactions, took %.1fs",
            n_interactions, time.time() - t0,
        )

    return lr_res


def format_cci_results(
    lr_res: pd.DataFrame,
    n_top: int = 50,
    pval_col: str = "magnitude_rank",
    log: object = None,
) -> pd.DataFrame:
    """Filter, sort and format CCI interaction results.

    Sorts by the given rank/significance column (lower = more significant),
    selects top N interactions, and adds a readable interaction label.

    Parameters
    ----------
    lr_res : pd.DataFrame
        Raw LIANA results from run_cci_permutation() or run_cci_spatial().
    n_top : int
        Number of top interactions to retain.
    pval_col : str
        Column to sort by (default 'magnitude_rank').
    log : object, optional
        Logger with .info() method.

    Returns
    -------
    pd.DataFrame
        Filtered, sorted top-N DataFrame with added 'interaction' label.
    """
    # Create readable interaction label
    lr_res = lr_res.copy()
    cols = lr_res.columns

    if "ligand" in cols and "receptor" in cols:
        src = lr_res["source"].astype(str) if "source" in cols else pd.Series("", index=lr_res.index)
        tgt = lr_res["target"].astype(str) if "target" in cols else pd.Series("", index=lr_res.index)
        if "source" in cols and "target" in cols:
            lr_res["interaction"] = src + "->" + tgt + " | " + lr_res["ligand"].astype(str) + "_" + lr_res["receptor"].astype(str)
        else:
            lr_res["interaction"] = lr_res["ligand"].astype(str) + "_" + lr_res["receptor"].astype(str)
    elif "ligand_complex" in cols and "receptor_complex" in cols:
        lr_res["interaction"] = (
            lr_res["ligand_complex"].astype(str) + "_" +
            lr_res["receptor_complex"].astype(str)
        )

    # Sort by significance
    if pval_col in cols:
        top_df = lr_res.sort_values(pval_col, ascending=True).head(n_top)
        if log:
            log.info("Top %d interactions selected by %s", n_top, pval_col)
    else:
        # Fallback: sort by the first available rank column
        rank_cols = [c for c in cols if "rank" in c.lower()]
        if rank_cols:
            top_df = lr_res.sort_values(rank_cols[0], ascending=True).head(n_top)
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
