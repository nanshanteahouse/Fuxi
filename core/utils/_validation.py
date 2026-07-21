"""Data validation — AnnData integrity checks and pipeline state validation."""

import logging
from typing import Any

import numpy as np
import scipy.sparse as sp


def validate_adata(adata, stage_name="", logger=None, fix_nan_inf=True) -> bool:
    """检查 AnnData X 矩阵完整性，自动修复 NaN/Inf。

    在后续步骤开始前调用，避免因前一步意外产生的 NaN/Inf
    导致下游（PCA、UMAP 等）崩溃。

    参数:
        adata: AnnData 对象
        stage_name: 当前步骤名称（仅用于日志标记）
        logger: 日志记录器（None 则自动获取）
        fix_nan_inf: 是否修复（替换为 0）

    返回:
        True — 发现并修复了问题
        False — X 矩阵清洁
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    x_data = adata.X.data if sp.issparse(adata.X) else adata.X
    n_nan = int(np.isnan(x_data).sum())
    n_inf = int(np.isinf(x_data).sum())
    total = n_nan + n_inf

    if total == 0:
        logger.info("[%s] X matrix clean — no NaN/Inf values", stage_name or "validate")
        return False

    logger.warning(
        "[%s] Found %d NaN and %d Inf values in X matrix — fixing",
        stage_name or "validate",
        n_nan,
        n_inf,
    )

    if fix_nan_inf:
        if sp.issparse(adata.X):
            adata.X.data = np.nan_to_num(adata.X.data, nan=0, posinf=0, neginf=0)
        else:
            adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)

    return True


# ── Cross-step data flow validation ──────────────────────────────────

_STEP_REQUIREMENTS = {
    "rna": {
        "01": {"obs": ["doublet_scores", "predicted_doublet"]},
        "02": {"obs": ["doublet_scores", "predicted_doublet"]},
        "04": {"obsm": ["X_pca"]},
        "05": {"obs": ["leiden"], "obsm": ["X_umap", "X_pca"]},
        "06": {"obs": ["cell_type"], "obsm": ["X_umap"]},
        "07": {"obs": ["cell_type", "leiden"]},
        "08": {"obs": ["cell_type"], "obsm": ["X_pca", "X_umap"]},
        "10": {"obs": ["cell_type"], "obsm": ["X_umap"]},
        "11": {"obs": ["cell_type"], "obsm": ["X_umap"]},
        "12": {"obs": ["cell_type"], "obsm": ["X_umap"]},
    },
    "atac": {
        "02": {"obs": ["predicted_doublet"]},
        "04": {"obs": ["leiden"], "obsm": ["X_umap"]},
        "05": {"obs": ["cell_type"], "obsm": ["X_umap"]},
        "06": {"obs": ["cell_type"]},
        "07": {"obs": ["cell_type"]},
        "08": {"obs": ["cell_type"]},
        "09": {"obs": ["cell_type"]},
    },
    "spatial": {
        "04": {"obsm": ["X_pca"]},
        "05": {"obs": ["leiden"], "obsm": ["X_umap"]},
        "06": {"obs": ["cell_type"], "obsp": ["spatial_connectivities"]},
        "07": {"obs": ["cell_type"], "obsm": ["X_umap"]},
        "09": {"obs": ["cell_type"], "obsm": ["X_umap"]},
        "10": {"obs": ["cell_type"]},
    },
}


def validate_pipeline_state(adata: Any, step: str, modality: str = "rna") -> None:
    """Assert that required obs/obsm/obsp columns exist at a step boundary.

    Parameters
    ----------
    adata : AnnData
        The current AnnData object to validate.
    step : str
        Step number (e.g. "00", "01", "02").
    modality : str
        One of "rna", "atac", or "spatial".

    Raises
    ------
    AssertionError
        If any required column is missing, with a clear message.
    """
    logger = logging.getLogger(__name__)

    # Normalise step key: try as-is first, then zero-padded to two digits.
    mod_reqs = _STEP_REQUIREMENTS.get(modality, {})
    step_key = step if step in mod_reqs else step.zfill(2)

    req = mod_reqs.get(step_key)
    if req is None:
        logger.debug(
            "[step=%s modality=%s] No requirements defined — skipping validation",
            step,
            modality,
        )
        return

    missing: list[str] = []

    for col in req.get("obs", []):
        if col not in adata.obs:
            missing.append(f"obs['{col}']")

    for col in req.get("obsm", []):
        if col not in adata.obsm:
            missing.append(f"obsm['{col}']")

    for col in req.get("obsp", []):
        if col not in adata.obsp:
            missing.append(f"obsp['{col}']")

    if missing:
        msg = f"Pipeline state validation failed at step {step} ({modality}):\n" + "\n".join(
            f"  - Missing: {m}" for m in missing
        )
        raise AssertionError(msg)

    logger.debug(
        "[step=%s modality=%s] All required columns present: obs=%s obsm=%s obsp=%s",
        step,
        modality,
        sorted(req.get("obs", [])),
        sorted(req.get("obsm", [])),
        sorted(req.get("obsp", [])),
    )
