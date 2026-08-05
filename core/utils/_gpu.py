"""GPU/CPU dispatch helpers for scanpy-equivalent operations.

This module is the single integration point between Fuxi and rapids-singlecell.
Each helper:

1. Reads the active ``cfg.execution.device`` setting (``auto`` / ``gpu`` / ``cpu``)
2. Routes to ``rapids_singlecell`` when GPU is requested and available
3. Falls back to the equivalent ``scanpy`` call when GPU is unavailable
4. Logs which path was taken (once per operation type, cached after first emit)

The dispatch decision is cached per-process after the first call to avoid
repeated ``nvidia-smi`` probes.

Public API (scanpy-shaped, single-implementation import pattern):

    from core.utils._gpu import gpu_neighbors, gpu_umap, gpu_leiden, gpu_pca, gpu_harmony

Designed so that step scripts can call these unconditionally — the runtime
decision happens here, not at call sites.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# ── Per-process GPU-availability cache ───────────────────────────────────
# Resolved on first auto-detect; reset to None only if caller explicitly
# invalidates (e.g. for tests). NOTE: This caches the *GPU probe* result,
# NOT the per-call decision — explicit device=cpu / device=gpu overrides
# must always be respected regardless of cache state.
_gpu_available_cache: bool | None = None


# ── Per-process log-once cache ─────────────────────────────────────────
def _fmt_gb(n_bytes: float) -> str:
    """Format bytes as a GB string with one decimal (e.g. 30.9GB)."""
    return f"{n_bytes / 1e9:.1f}GB"


# ── Per-process log-once cache ─────────────────────────────────────────
# Prevents redundant dispatch logging during grid searches that call
# gpu_* functions dozens of times without the device decision changing.
_auto_device_logged: str | None = None  # last auto-resolve path logged
_dispatched_ops_logged: set[str] = set()  # op names already emitted


def _auto_detect_gpu() -> bool:
    """Cached gpu_available_rapids() probe. Set/reset by reset_device_cache()."""
    global _gpu_available_cache
    if _gpu_available_cache is None:
        from core.utils._optional import gpu_available_rapids

        _gpu_available_cache = gpu_available_rapids()
    return _gpu_available_cache


def resolve_device(device: str = "auto", log: logging.Logger | None = None) -> bool:
    """Resolve the effective device mode for THIS call.

    Args:
        device: ``"auto"`` (detect), ``"gpu"" (force), or ``"cpu"" (force).
        log: optional logger for the resolution decision.

    Returns:
        True if GPU path should be used, False if CPU.

    Raises:
        RuntimeError: if device="gpu" but RAPIDS/CUDA unavailable.

    Note:
        Explicit device="cpu" / device="gpu" overrides are honored on EVERY call;
        only the underlying gpu_available_rapids() probe is cached. This avoids the
        footgun where the first call (often "auto") poisons the cache and silences
        later force-cpu requests.
    """
    if device == "cpu":
        if log is not None:
            log.info("[device] CPU forced (device=cpu)")
        return False

    gpu_ok = _auto_detect_gpu()
    if device == "gpu":
        if not gpu_ok:
            raise RuntimeError(
                "device=gpu but rapids-singlecell + CUDA are not available. "
                "Install with: pip install fuxi[rapids] --extra-index-url=https://pypi.nvidia.com"
            )
        if log is not None:
            log.info("[device] GPU forced (device=gpu)")
        return True

    # device == "auto"
    global _auto_device_logged
    path = "GPU (rapids-singlecell)" if gpu_ok else "CPU (scanpy)"
    if log is not None and _auto_device_logged != path:
        log.info("[device] auto-resolved to %s", path)
        _auto_device_logged = path
    return gpu_ok


def reset_device_cache() -> None:
    """Reset the cached GPU probe result. Mainly for tests."""
    global _gpu_available_cache, _auto_device_logged, _dispatched_ops_logged
    _gpu_available_cache = None
    _auto_device_logged = None
    _dispatched_ops_logged.clear()


def is_gpu_active() -> bool:
    """Read-only probe of cached GPU availability. False before any resolve_device call."""
    return bool(_gpu_available_cache)


# ═══════════════════════════════════════════════════════════════════════
# Per-operation dispatchers
# ═══════════════════════════════════════════════════════════════════════
# Each function:
#   - Mirrors the scanpy signature for the equivalent op (so call sites stay
#     unchanged)
#   - Calls rapids_singlecell when GPU mode is active
#   - Falls through to scanpy otherwise
#   - First GPU call in a process triggers a one-time adata → GPU transfer
#     via rsc.get.anndata_to_GPU(); subsequent ops reuse the pinned buffer.


def gpu_neighbors(
    adata,
    log: logging.Logger | None = None,
    device: str = "auto",
    **kwargs: Any,
):
    """scanpy.pp.neighbors equivalent with GPU dispatch.

    GPU path uses ``rapids_singlecell.pp.neighbors`` (cuml brute_force / ivfflat).
    """
    if resolve_device(device, log):
        import rapids_singlecell as rsc

        rsc.get.anndata_to_GPU(adata)
        if log is not None and "neighbors" not in _dispatched_ops_logged:
            log.info("[device] sc.pp.neighbors → rsc.pp.neighbors (GPU)")
            _dispatched_ops_logged.add("neighbors")
        return rsc.pp.neighbors(adata, **kwargs)
    import scanpy as sc

    return sc.pp.neighbors(adata, **kwargs)


def gpu_umap(
    adata,
    log: logging.Logger | None = None,
    device: str = "auto",
    **kwargs: Any,
):
    """scanpy.tl.umap equivalent with GPU dispatch.

    GPU path uses ``rapids_singlecell.tl.umap`` (cuml UMAP). Speed-up on 1M
    cells is typically 10-30x vs CPU scanpy+numba.
    """
    if resolve_device(device, log):
        import rapids_singlecell as rsc

        rsc.get.anndata_to_GPU(adata)
        # rsc.tl.umap accepts maxiter but NOT n_epochs — cuml derives n_epochs internally
        gpu_kwargs = dict(kwargs)
        gpu_kwargs.pop("n_epochs", None)
        if log is not None and "umap" not in _dispatched_ops_logged:
            log.info("[device] sc.tl.umap → rsc.tl.umap (GPU)")
            _dispatched_ops_logged.add("umap")
        return rsc.tl.umap(adata, **gpu_kwargs)
    import scanpy as sc

    return sc.tl.umap(adata, **kwargs)


def gpu_leiden(
    adata,
    log: logging.Logger | None = None,
    device: str = "auto",
    **kwargs: Any,
):
    """scanpy.tl.leiden equivalent with GPU dispatch.

    GPU path uses ``rapids_singlecell.tl.leiden`` (cuGraph Leiden).
    Note: cuGraph Leiden does not support all the flavor= variants scanpy
    does; only ``flavor="igraph"`` (the default RAPIDS backend) is honored.
    """
    if resolve_device(device, log):
        import rapids_singlecell as rsc

        rsc.get.anndata_to_GPU(adata)
        # cuGraph ignores flavor; force-remove incompatible kwargs for parity
        gpu_kwargs = dict(kwargs)
        gpu_kwargs.pop("flavor", None)
        gpu_kwargs.pop("directed", None)
        gpu_kwargs.pop("n_iterations", None)  # cuGraph uses max_iter
        if log is not None and "leiden" not in _dispatched_ops_logged:
            log.info("[device] sc.tl.leiden → rsc.tl.leiden (GPU)")
            _dispatched_ops_logged.add("leiden")
        return rsc.tl.leiden(adata, **gpu_kwargs)
    import scanpy as sc

    return sc.tl.leiden(adata, **kwargs)


def gpu_pca(
    adata,
    log: logging.Logger | None = None,
    device: str = "auto",
    cfg: Any = None,
    step: str = "03_integrate",
    **kwargs: Any,
):
    """scanpy.pp.pca equivalent with GPU dispatch.

    GPU path uses ``rapids_singlecell.pp.pca`` (cuVS / cuSOLVER SVD).
    Speed-up is most visible for n_obs > 100k.

    ``cfg`` (optional): resolved config — when provided, a VRAM-triggered
    GPU→CPU fallback is recorded to ``<results_dir>/memory_skips.jsonl``
    via :func:`core.utils._perf.record_memory_skip`, mirroring the
    memory_policy skip audit.
    """
    if resolve_device(device, log):
        # NOTE(2026-08-01): 显存 guard —— rsc.pp.pca 需要 X dense 化进显存
        # (anndata_to_GPU 把 CSR → dense)。估算 dense 体积 (n_obs × n_vars × dtype)
        # 超过可用显存时自动降级 CPU，避免 2M 细胞 × 4k 基因 (≈30GB) 直接 OOM 显存。
        import cupy as cp
        import rapids_singlecell as rsc

        n_obs, n_vars = adata.n_obs, adata.n_vars
        dtype_size = np.dtype(adata.X.dtype).itemsize if hasattr(adata.X, "dtype") else 4
        dense_bytes = n_obs * n_vars * dtype_size
        free_bytes, _ = cp.cuda.runtime.memGetInfo()
        if dense_bytes > int(free_bytes * 0.9):
            if log is not None:
                log.warning(
                    "[device] PCA dense %s exceeds free VRAM %.1fGB — falling back to CPU (arpack)",
                    _fmt_gb(dense_bytes),
                    free_bytes / 1e9,
                )
            from core.utils._perf import record_memory_skip  # lazy: _perf imports _gpu

            record_memory_skip(
                step=step,
                operation="pca GPU→CPU fallback",
                reason=(
                    f"dense X {_fmt_gb(dense_bytes)} exceeds free VRAM "
                    f"{free_bytes / 1e9:.1f}GB (0.9x guard) — n_obs={n_obs}, n_vars={n_vars}"
                ),
                cfg=cfg,
                log=log,
            )
        else:
            rsc.get.anndata_to_GPU(adata)
            if log is not None and "pca" not in _dispatched_ops_logged:
                log.info("[device] sc.pp.pca → rsc.pp.pca (GPU)")
                _dispatched_ops_logged.add("pca")
            return rsc.pp.pca(adata, **kwargs)
    import scanpy as sc

    return sc.pp.pca(adata, **kwargs)


def gpu_harmony(
    adata,
    key: str | list[str],
    *,
    output_key: str = "X_integrated",
    log: logging.Logger | None = None,
    device: str = "auto",
    **kwargs: Any,
):
    """harmonypy.run_harmony equivalent with GPU dispatch.

    Both modes normalize the output to ``adata.obsm[output_key]`` (default
    ``'X_integrated'``) so callers don't need to branch on the backend.

    GPU path uses ``rapids_singlecell.pp.harmony_integrate``, which writes
    to ``adata.obsm[adjusted_basis]`` (default ``'X_pca_harmony'``); this
    wrapper sets ``adjusted_basis=output_key`` so callers don't need to branch
    on the backend.

    Returns:
        ``None`` — caller reads the corrected embedding from
        ``adata.obsm[output_key]``.
    """
    keys_list = [key] if isinstance(key, str) else list(key)
    if resolve_device(device, log):
        import rapids_singlecell as rsc

        # NOTE(2026-08-01): 不调用 anndata_to_GPU —— harmony 只需 obsm[basis]
        # (X_pca)，rsc.pp.harmony_integrate 内部自行 cp.array() 搬 GPU 并写回 CPU
        # (harmony_out.get())。之前无脑 anndata_to_GPU 会把整个 X dense 化搬进
        # 显存 (2M 细胞 × 4k 基因 ≈ 30GB > 24GB 显存墙)，导致大数据无法走 GPU
        # harmony。只搬 X_pca (n×100 float32 ≈ 0.7GB) 即可。
        if log is not None and "harmony" not in _dispatched_ops_logged:
            log.info("[device] harmonypy → rsc.pp.harmony_integrate (GPU, basis=X_pca only)")
            _dispatched_ops_logged.add("harmony")
        # rsc.pp.harmony_integrate writes to adjusted_basis (default X_pca_harmony).
        # Pass through any caller-provided kwargs (random_state, max_iter_harmony, etc.)
        # but let adjusted_basis be controlled by output_key (with caller override possible).
        kwargs.setdefault("adjusted_basis", output_key)
        rsc.pp.harmony_integrate(adata, key, **kwargs)
        return None

    import harmonypy as hm

    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        vars_use=keys_list,
        **kwargs,
    )
    adata.obsm[output_key] = ho.Z_corr
    return None


# ── Transfer-out helper ────────────────────────────────────────────────
def sync_to_cpu(adata, log: logging.Logger | None = None) -> None:
    """Pull GPU-resident arrays back to host memory.

    Call before any h5ad write or before passing adata to CPU-only code paths
    (matplotlib plotting, joblib workers without RAPIDS, etc.).
    """
    if not is_gpu_active():
        return
    import rapids_singlecell as rsc

    rsc.get.anndata_to_CPU(adata)
    if log is not None:
        log.info("[device] adata transferred GPU → CPU")
