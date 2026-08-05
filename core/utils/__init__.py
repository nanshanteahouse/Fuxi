#!/usr/bin/env python3
# ruff: noqa: E402
"""
utils — Fuxi pipeline shared utilities
========================================

Sub-modules (private, import through this package):
  _path           — WSL detection, data_root, repo_root
  _io             — safe_write, safe_plot
  _logging        — setup_logger
  _config         — resolve_config, species validation, dataset.yaml helpers
  _cross_modality — find_rna_h5ad, find_rna_marker_csv, load_scRNA_markers
  _perf           — PerformanceReport, PerformanceSummary, monitor_performance, timed_substep, record_memory_skip
  _validation     — validate_adata, validate_pipeline_state
"""

import os
import platform

# ── WSL h5py file locking auto-detection ────────────────────────────
if "microsoft" in platform.release().lower():
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# ── Auto-load .env ──────────────────────────────────────────────────
try:
    from dotenv import load_dotenv

    _dotenv_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        ".env",
    )
    if os.path.exists(_dotenv_path):
        load_dotenv(_dotenv_path, override=True)
except ImportError:
    pass


# ── NVIDIA CUDA runtime auto-discovery (.venv-installed nvidia-* wheels) ──
# pip-installed RAPIDS / cupy wheels ship .so files under nvidia/*/lib but
# ld.so does not search those paths by default. Without this preload, any
# `import rapids_singlecell` (or cuml/cugraph) fails with:
#   ImportError: libcudart.so.12: cannot open shared object file
# We preload the critical libs by absolute path AND add their dirs to
# LD_LIBRARY_PATH so subsequent dlopen calls (incl. child subprocesses via
# run_pipeline.py) succeed. Must run BEFORE any `import rapids_singlecell`.
def _preload_nvidia_cuda() -> None:
    import ctypes
    import glob
    import sys

    site_pkgs = next((p for p in sys.path if p.endswith("site-packages")), None)
    if not site_pkgs:
        return
    lib_dirs = glob.glob(os.path.join(site_pkgs, "nvidia", "*", "lib"))
    if not lib_dirs:
        return  # nvidia-* wheels not installed (CPU-only env), no-op

    # Append (do not overwrite) to LD_LIBRARY_PATH for child processes
    existing = os.environ.get("LD_LIBRARY_PATH", "")
    existing_set = set(existing.split(":")) if existing else set()
    new_dirs = [d for d in lib_dirs if d not in existing_set]
    if new_dirs:
        merged = ":".join(new_dirs + ([existing] if existing else []))
        os.environ["LD_LIBRARY_PATH"] = merged

    # Preload critical CUDA libs by absolute path (bypasses ld.so startup cache).
    # RTLD_GLOBAL makes symbols visible to subsequently-loaded libs.
    # Try both .so.12 and .so.13 — pip wheels exist in both CUDA major versions.
    _critical = (
        "libcudart.so.12",
        "libcudart.so.13",
        "libcublas.so.12",
        "libcublas.so.13",
        "libcusparse.so.12",
        "libcusparse.so.13",
        "libcusolver.so.11",
        "libcusolver.so.12",
        "libcufft.so.11",
        "libcufft.so.12",
        "libcurand.so.10",
        "libcudnn.so.12",
        "libnccl.so.2",
        "libnvrtc.so.12",
        "libnvrtc.so.13",
    )
    for _d in lib_dirs:
        for _lib in _critical:
            _p = os.path.join(_d, _lib)
            if os.path.exists(_p):
                try:
                    ctypes.CDLL(_p, mode=ctypes.RTLD_GLOBAL)
                except OSError:
                    pass


_preload_nvidia_cuda()

# ── Re-exports ──────────────────────────────────────────────────────
from core.utils._config import (  # noqa: F401
    _KNOWN_SPECIES_KEYS,
    _find_dataset_yaml,
    _validate_species,
    resolve_config,
)
from core.utils._cross_modality import (  # noqa: F401
    find_rna_h5ad,
    find_rna_marker_csv,
    load_scRNA_markers,
)
from core.utils._exec import _set_blas_env  # noqa: F401
from core.utils._gpu import (  # noqa: F401
    gpu_harmony,
    gpu_leiden,
    gpu_neighbors,
    gpu_pca,
    gpu_umap,
    is_gpu_active,
    reset_device_cache,
    resolve_device,
    sync_to_cpu,
)
from core.utils._io import (  # noqa: F401
    safe_plot,
    safe_write,
    save_figure,
    stream_write_raw,
    write_obs_columns_inplace,
    write_obs_columns_lightweight,
)
from core.utils._io_incremental import write_h5ad_incremental  # noqa: F401
from core.utils._logging import setup_logger  # noqa: F401
from core.utils._memory import (  # noqa: F401
    check_memory_guard,
    estimate_step_peak,
    resolve_memory_budget_bytes,
    resolve_memory_settings,
)
from core.utils._optional import (  # noqa: F401
    gpu_available_nvidia_smi,
    gpu_available_rapids,
    gpu_available_torch,
    require_celltypist,
    require_rapids,
    require_scvi,
)
from core.utils._path import (  # noqa: F401
    _DATA_ROOT_CACHE,
    _REPO_ROOT_CACHE,
    data_root,
    is_wsl,
    repo_root,
    wsl_to_win,
)
from core.utils._perf import (  # noqa: F401
    PerformanceReport,
    PerformanceSummary,
    monitor_performance,
    record_memory_skip,
    timed_substep,
)
from core.utils._validation import (  # noqa: F401
    _STEP_REQUIREMENTS,
    validate_adata,
    validate_pipeline_state,
)
