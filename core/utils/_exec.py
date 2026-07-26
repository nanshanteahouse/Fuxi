"""Execution / subprocess helpers for the Fuxi pipeline."""

from __future__ import annotations

import os

__all__ = ["_set_blas_env"]


_BLAS_ENV_VARS: tuple[str, ...] = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "PYTORCH_ENABLE_MPS_FALLBACK",
    "TORCH_NUM_THREADS",
)


def _set_blas_env(n_jobs: int, *, overwrite: bool = True) -> dict[str, str]:
    """Set BLAS / OpenMP thread-limit environment variables.

    Prevents numerical-library thread oversubscription when running
    multiple pipeline steps concurrently on a shared node.

    Parameters
    ----------
    n_jobs:
        Number of threads to cap each BLAS/OpenMP library to.
        Must be >= 1.  Values <= 0 are silently ignored (no-op).
    overwrite:
        If True, force-set each variable (``os.environ[key] = val``).
        If False, only set variables not already present
        (``os.environ.setdefault(key, val)``).

    Returns
    -------
    dict[str, str]
        The subset of ``os.environ`` entries that were set or updated
        by this call (empty dict when *n_jobs* <= 0).
    """
    if n_jobs < 1:
        return {}

    value = str(n_jobs)
    changed: dict[str, str] = {}

    if overwrite:
        for var in _BLAS_ENV_VARS:
            os.environ[var] = value
            changed[var] = value
    else:
        for var in _BLAS_ENV_VARS:
            if var not in os.environ:
                os.environ[var] = value
                changed[var] = value

    return changed
