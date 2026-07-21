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
