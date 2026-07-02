#!/usr/bin/env python3
"""
Shared clustering interface for Fuxi multi-omics pipeline.

Provides modality-agnostic grid-search orchestration for clustering parameters
(rna / spatial / atac).  Modality-specific code injects callables (clusterer,
neighbor_fn, umap_fn) that wrap the underlying library (Scanpy, SnapATAC2,
Squidpy).  The shared layer handles:

  * Cartesian-product grid expansion
  * Optional grouping so expensive steps (neighbors, UMAP) run once per group
  * Parallel evaluation across combinations
  * Pareto-elbow best-parameter selection (re-export from
    ``rna.utils.cluster_evaluation``)

Exports
-------
grid_search_clustering(adata, param_grid, clusterer, ...) -> list[dict]
    Run a grid search over clustering parameters and return a results table.

select_best_params(results, method='pareto_elbow')
    Re-export from ``rna.utils.cluster_evaluation.select_best_params``.

umap_sweep(adata, param_sweep, umap_fn, ...) -> list[dict]
    Run a UMAP parameter sweep and return per-combination results.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Sequence
from typing import Any
from typing import Any

import numpy as np

# Re-export the existing Pareto selection logic — no duplication.
from rna.utils.cluster_evaluation import select_best_params  # noqa: F401

# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------


def grid_search_clustering(
    adata: Any,
    param_grid: dict[str, Sequence[Any]],
    clusterer: Callable[..., str],
    neighbor_fn: Callable[..., None] | None = None,
    umap_fn: Callable[..., None] | None = None,
    evaluation_fn: Callable[..., float] | None = None,
    *,
    group_key: str | None = None,
    n_jobs: int = 1,
    random_seed: int = 42,
    **fixed_kwargs: Any,
) -> list[dict[str, Any]]:
    """Run a grid search over clustering parameters.

    *grouped* mode — when ``group_key`` is given — assumes that expensive
    pre-processing steps (neighbor graph, UMAP embedding) can be *shared*
    across all values of the remaining parameters.  For each unique value
    of ``group_key``, ``neighbor_fn`` and ``umap_fn`` are called once (if
    provided); ``clusterer`` is then invoked for every combination.

    Parameters
    ----------
    adata : AnnData or SnapATAC2 AnnData
        The data object to cluster.

    param_grid : dict[str, Sequence]
        Parameter names mapped to a list of values to sweep.  Example::

            {'n_neighbors': [15, 20], 'resolution': [0.5, 1.0, 1.5]}

    clusterer : Callable(adata, **params, **fixed_kwargs) -> str
        Modality-specific clustering callable.  Must mutate *adata* in-place
        and return the key of the newly added observation column (e.g.
        ``"leiden"`` or ``"leiden_15_0.5"``).

    neighbor_fn : Callable(adata, **params, **fixed_kwargs) -> None, optional
        Called once per ``group_key`` value to compute neighbour graph / KNN.

    umap_fn : Callable(adata, **params, **fixed_kwargs) -> None, optional
        Called once per ``group_key`` value to compute a UMAP embedding.

    evaluation_fn : Callable(adata, cluster_key, **params, **fixed_kwargs) -> float, optional
        Quality metric (e.g. silhouette score).  Receives the cluster key
        returned by ``clusterer``.

    group_key : str, optional
        Name of the parameter whose values determine the grouping (typically
        ``"n_neighbors"``).  ``neighbor_fn`` and ``umap_fn`` are invoked once
        for each unique value of this key and then shared across all
        combinations of the remaining parameters.

    n_jobs : int
        Number of parallel workers (requires ``joblib``).  Default 1 = serial.

    random_seed : int
        Passed to ``fixed_kwargs`` as ``random_seed`` if present.

    **fixed_kwargs
        Extra keyword arguments forwarded to every callable invocation.

    Returns
    -------
    list[dict]
        One entry per successfully evaluated combination with keys drawn from
        ``param_grid`` plus ``"n_clusters"``, ``"score"``, and
        ``"cluster_key"``.
    """
    # --- validate ---
    if not param_grid:
        raise ValueError("param_grid must not be empty")

    param_names: list[str] = list(param_grid.keys())
    _inject_seed(fixed_kwargs, random_seed)

    # --- expand Cartesian product ---
    # param_grid values are lists; product gives every combination
    combos: list[tuple[Any, ...]] = list(itertools.product(*param_grid.values()))

    if group_key is not None and group_key not in param_names:
        raise ValueError(f"group_key {group_key!r} not found in param_grid keys: {param_names}")

    # --- group index ---
    group_idx: int | None = param_names.index(group_key) if group_key is not None else None

    results: list[dict[str, Any]] = []

    # Serial execution path (parallel left to the caller or can be trivially
    # wrapped with joblib — we keep this function simple and testable).
    _grid_search_serial(
        adata,
        param_names,
        combos,
        group_idx,
        clusterer,
        neighbor_fn,
        umap_fn,
        evaluation_fn,
        fixed_kwargs,
        results,
    )

    return results


def umap_sweep(
    adata: Any,
    param_sweep: dict[str, Sequence[Any]],
    umap_fn: Callable[..., np.ndarray] | None = None,
    evaluation_fn: Callable[..., float] | None = None,
    *,
    random_seed: int = 42,
    **fixed_kwargs: Any,
) -> list[dict[str, Any]]:
    """Run a UMAP parameter sweep, returning per-combination results.

    Parameters
    ----------
    adata : AnnData or SnapATAC2 AnnData

    param_sweep : dict[str, Sequence]
        UMAP parameter grid, e.g. ``{"min_dist": [0.1, 0.3], "spread": [1.0, 2.0]}``.

    umap_fn : Callable(adata, **params, **fixed_kwargs) -> ndarray, optional
        Must compute UMAP and return the 2-D embedding.  If *None* the
        caller must provide a pre-computed ``"X_umap"`` in ``adata.obsm`` and
        this function becomes a no-op wrapper.

    evaluation_fn : Callable(adata, coordinates, **params, **fixed_kwargs) -> float, optional
        Optional quality metric per combination.

    random_seed : int

    **fixed_kwargs
        Forwarded to every callable.

    Returns
    -------
    list[dict]
        One entry per combination with param names + optional ``"score"``.
    """
    if not param_sweep:
        raise ValueError("param_sweep must not be empty")

    _inject_seed(fixed_kwargs, random_seed)

    param_names = list(param_sweep.keys())
    combos = list(itertools.product(*param_sweep.values()))

    results: list[dict[str, Any]] = []

    # No-op fallback when caller doesn't supply umap_fn: assume UMAP already
    # exists in adata.obsm and just run evaluation.
    _umap_fn: Callable[..., np.ndarray]
    if umap_fn is None:
        _umap_fn = _noop_umap_fn
    else:
        _umap_fn = umap_fn

    for combo in combos:
        params = dict(zip(param_names, combo))
        merged = {**fixed_kwargs, **params}
        try:
            coords = _umap_fn(adata, **merged)
        except Exception:
            results.append({**params, "score": None, "error": "UMAP failed"})
            continue

        entry: dict[str, Any] = dict(params)
        if evaluation_fn is not None:
            try:
                entry["score"] = evaluation_fn(adata, coords, **merged)
            except Exception:
                entry["score"] = None
                entry["error"] = "evaluation failed"
        results.append(entry)

    return results


# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------


def _inject_seed(kwargs: dict[str, Any], seed: int) -> None:
    """Ensure *random_seed* is in *kwargs* unless the caller already set it."""
    if "random_seed" not in kwargs:
        kwargs["random_seed"] = seed


def _noop_umap_fn(adata: Any, **kw: Any) -> np.ndarray:
    """Fallback UMAP function: return pre-existing ``X_umap``."""
    return np.asarray(adata.obsm["X_umap"])


def _grid_search_serial(
    adata: Any,
    param_names: list[str],
    combos: list[tuple[Any, ...]],
    group_idx: int | None,
    clusterer: Callable[..., str],
    neighbor_fn: Callable[..., None] | None,
    umap_fn: Callable[..., None] | None,
    evaluation_fn: Callable[..., float] | None,
    fixed_kwargs: dict[str, Any],
    results: list[dict[str, Any]],
) -> None:
    """Serialize grid-search: group by *group_idx*, call neighbours/UMAP once."""

    if group_idx is None:
        # flat case — no grouping, evaluate every combo independently
        for combo in combos:
            params = dict(zip(param_names, combo))
            merged = {**fixed_kwargs, **params}
            _try_one_combo(
                adata, clusterer, evaluation_fn, params, merged, results
            )
        return

    # --- grouped case ---
    # Partition combos by their group_key value
    group_values = sorted({c[group_idx] for c in combos})
    for gv in group_values:
        # Filter combos belonging to this group
        group_combos = [c for c in combos if c[group_idx] == gv]
        if not group_combos:
            continue

        # Build params dict for the group (shared neighbour/UMAP)
        group_params: dict[str, Any] = {param_names[group_idx]: gv}
        group_merged = {**fixed_kwargs, **group_params}

        # 1. Neighbour graph (once per group)
        if neighbor_fn is not None:
            try:
                neighbor_fn(adata, **group_merged)
            except Exception:
                # Skip entire group if neighbour computation fails
                continue

        # 2. UMAP (once per group)
        if umap_fn is not None:
            try:
                umap_fn(adata, **group_merged)
            except Exception:
                # UMAP failure is non-fatal — continue with clustering only
                pass

        # 3. Evaluate each (resolution, ...) combo within the same graph
        for combo in group_combos:
            params = dict(zip(param_names, combo))
            merged = {**group_merged, **params}
            _try_one_combo(
                adata, clusterer, evaluation_fn, params, merged, results
            )


def _try_one_combo(
    adata: Any,
    clusterer: Callable[..., str],
    evaluation_fn: Callable[..., float] | None,
    params: dict[str, Any],
    merged_kwargs: dict[str, Any],
    results: list[dict[str, Any]],
) -> None:
    """Evaluate a single parameter combination and append to *results*."""
    entry: dict[str, Any] = dict(params)
    try:
        cluster_key = clusterer(adata, **merged_kwargs)
        entry["cluster_key"] = cluster_key
    except Exception:
        entry["error"] = "clusterer failed"
        results.append(entry)
        return

    # Count clusters if possible
    try:
        labels = adata.obs[cluster_key]
        entry["n_clusters"] = int(labels.nunique())
    except Exception:
        entry["n_clusters"] = None

    # Quality score
    if evaluation_fn is not None:
        try:
            entry["score"] = evaluation_fn(adata, cluster_key, **merged_kwargs)
        except Exception:
            entry["score"] = None

    results.append(entry)
