"""Optional dependency guards and GPU detection for scVI integration."""

import importlib.util
import logging
import shutil
import subprocess

logger = logging.getLogger(__name__)

# Module-level cache for scvi availability
_scvi_available: bool | None = None


def require_scvi(feature: str = "scVI integration") -> None:
    """Lazily check if scvi-tools is available.

    Results are cached at module level after the first check.
    Raises ``ImportError`` with a clear install hint when scvi-tools is missing.
    """
    global _scvi_available
    if _scvi_available is None:
        _scvi_available = importlib.util.find_spec("scvi") is not None
    if not _scvi_available:
        raise ImportError(
            f"scvi-tools is required for {feature}. Install with: pip install fuxi[scvi]"
        )


# Module-level cache for celltypist availability
_celltypist_available: bool | None = None


def require_celltypist(feature: str = "CellTypist annotation") -> None:
    """Lazily check if celltypist is available.

    Results are cached at module level after the first check.
    Raises ``ImportError`` with a clear install hint when celltypist is missing.
    """
    global _celltypist_available
    if _celltypist_available is None:
        _celltypist_available = importlib.util.find_spec("celltypist") is not None
    if not _celltypist_available:
        raise ImportError(
            f"celltypist is required for {feature}. Install with: pip install fuxi[celltypist]"
        )


# Module-level cache for sccoda availability
_sccoda_available: bool | None = None


def require_sccoda(feature: str = "scCODA compositional analysis") -> None:
    """Lazily check if sccoda is available.

    Results are cached at module level after the first check.
    Raises ``ImportError`` with a clear install hint when sccoda is missing.
    """
    global _sccoda_available
    if _sccoda_available is None:
        _sccoda_available = importlib.util.find_spec("sccoda") is not None
    if not _sccoda_available:
        raise ImportError(
            f"sccoda is required for {feature}. Install with: pip install fuxi[sccoda]"
        )


# Module-level cache for scvelo availability
_scvelo_available: bool | None = None


def require_scvelo(feature: str = "scVelo RNA velocity") -> None:
    """Lazily check if scvelo is available.

    Results are cached at module level after the first check.
    Raises ``ImportError`` with a clear install hint when scvelo is missing.
    """
    global _scvelo_available
    if _scvelo_available is None:
        _scvelo_available = importlib.util.find_spec("scvelo") is not None
    if not _scvelo_available:
        raise ImportError(
            f"scvelo is required for {feature}. Install with: pip install fuxi[scvelo]"
        )


# Module-level cache for cellbender availability
_cellbender_available: bool | None = None


def require_cellbender(feature: str = "CellBender ambient RNA removal") -> None:
    """Lazily check if cellbender is available.

    Results are cached at module level after the first check.
    Raises ``ImportError`` with a clear install hint when cellbender is missing.
    """
    global _cellbender_available
    if _cellbender_available is None:
        _cellbender_available = importlib.util.find_spec("cellbender") is not None
    if not _cellbender_available:
        raise ImportError(
            f"cellbender is required for {feature}. Install with: pip install fuxi[cellbender]"
        )


# Module-level cache for soupx availability
_soupx_available: bool | None = None


def require_soupx(feature: str = "SoupX ambient RNA removal") -> None:
    """Lazily check if soupx is available.

    Results are cached at module level after the first check.
    Raises ``ImportError`` with a clear install hint when soupx is missing.
    """
    global _soupx_available
    if _soupx_available is None:
        _soupx_available = importlib.util.find_spec("soupx") is not None
    if not _soupx_available:
        raise ImportError(
            f"soupx is required for {feature}. "
            "Install with: pip install fuxi[soupx] "
            "(third-party Python port, soupx-python; "
            "for the canonical R implementation use the SoupX R package)."
        )


def gpu_available_nvidia_smi() -> bool:
    """Check for an NVIDIA GPU by probing ``nvidia-smi``.

    Returns ``True`` only when the binary is found and exits with code 0.
    Returns ``False`` if the binary is missing, the call times out (5 s),
    exits non-zero, or any ``OSError`` occurs.
    """
    if shutil.which("nvidia-smi") is None:
        return False
    try:
        result = subprocess.run(
            ["nvidia-smi"],
            timeout=5,
            capture_output=True,
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError, OSError):
        return False


def gpu_available_torch() -> bool:
    """Authoritative GPU check via ``torch.cuda.is_available()``.

    Only runs if ``nvidia-smi`` pre-check passes, then lazily imports
    torch. Returns ``False`` gracefully if torch is not installed.
    """
    if not gpu_available_nvidia_smi():
        return False
    try:
        import torch

        return torch.cuda.is_available()
    except ImportError:
        return False
