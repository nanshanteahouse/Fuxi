"""
downsample.py — 细胞降采样核心逻辑 (config-driven + 手动 CLI 共同使用)
======================================================================

导出:
  - downsample_by_config(adata, cfg, logger) → AnnData   # 主要入口
  - downsample_random / downsample_stratified / downsample_max_per_sample  # 可直接调用
  - estimate_memory_gb
  - _check_sample_col
"""

import numpy as np
import scanpy as sc
import scipy.sparse as sp


def _check_sample_col(adata: sc.AnnData, sample_key: str, log) -> str | None:
    """查找可用的样本分组列。返回实际使用的列名或 None。"""
    if sample_key and sample_key in adata.obs:
        return sample_key
    for candidate in ['sample', 'Sample', 'samples', 'batch', 'Batch', 'stage', 'Stage']:
        if candidate in adata.obs:
            log.info("Using '%s' as group column ('%s' not found)", candidate, sample_key)
            return candidate
    return None


def downsample_random(adata: sc.AnnData, target: int, rng: np.random.RandomState,
                      log) -> sc.AnnData:
    """完全随机采样 target 个细胞。"""
    n_cells = adata.n_obs
    if target >= n_cells:
        log.info("target_total (%d) >= current cell count (%d), no downsampling needed", target, n_cells)
        return adata
    idx = rng.choice(n_cells, size=target, replace=False)
    idx.sort()
    log.info("Random sampling: %d → %d cells (%.1f%%)", n_cells, target, 100 * target / n_cells)
    return adata[idx].copy()


def downsample_stratified(adata: sc.AnnData, target: int, sample_key: str,
                          rng: np.random.RandomState, log) -> sc.AnnData:
    """按样本分层采样，保持各样本比例。"""
    n_cells = adata.n_obs
    if target >= n_cells:
        log.info("target_total (%d) >= current cell count (%d), no downsampling needed", target, n_cells)
        return adata

    counts = adata.obs[sample_key].value_counts()
    log.info("Stratified sampling, group=%s, target_total=%d", sample_key, target)
    for s, c in counts.items():
        log.info("  Sample %s: %d cells (%.1f%%)", s, c, 100 * c / n_cells)

    fractions = counts / n_cells
    per_sample_targets = (fractions * target).astype(int)
    remainder = target - per_sample_targets.sum()
    if remainder > 0:
        sorted_idx = np.argsort((fractions * target) - per_sample_targets)[::-1]
        for i in range(remainder):
            per_sample_targets.iloc[int(sorted_idx[i])] += 1

    indices = []
    for sample_name in counts.index:
        mask = adata.obs[sample_key] == sample_name
        sample_idx = np.where(mask)[0]
        n_sample = len(sample_idx)
        t = min(per_sample_targets[sample_name], n_sample)
        if t < n_sample:
            chosen = rng.choice(sample_idx, size=t, replace=False)
        else:
            chosen = sample_idx
        indices.append(chosen)

    idx = np.concatenate(indices)
    idx.sort()
    log.info("Stratified sampling: %d → %d cells (%.1f%%)", n_cells, len(idx), 100 * len(idx) / n_cells)
    return adata[idx].copy()


def downsample_max_per_sample(adata: sc.AnnData, max_per: int, sample_key: str,
                              rng: np.random.RandomState, log) -> sc.AnnData:
    """每个样本最多保留 max_per 个细胞。"""
    counts = adata.obs[sample_key].value_counts()
    log.info("Capping per sample, max %d cells per sample", max_per)

    indices = []
    for sample_name in counts.index:
        mask = adata.obs[sample_key] == sample_name
        sample_idx = np.where(mask)[0]
        n_sample = len(sample_idx)
        if n_sample > max_per:
            chosen = rng.choice(sample_idx, size=max_per, replace=False)
            log.info("  Sample %s: %d → %d (truncated %d)", sample_name, n_sample, max_per, n_sample - max_per)
        else:
            chosen = sample_idx
            log.info("  Sample %s: %d (unchanged)", sample_name, n_sample)
        indices.append(chosen)

    idx = np.concatenate(indices)
    idx.sort()
    log.info("Capped sampling: %d → %d cells (%.1f%%)", adata.n_obs, len(idx), 100 * len(idx) / adata.n_obs)
    return adata[idx].copy()


def estimate_memory_gb(adata: sc.AnnData) -> float:
    """粗略估计 AnnData 在内存中的大小 (GB)。"""
    total = 0.0
    if hasattr(adata, 'X') and adata.X is not None:
        if sp.issparse(adata.X):
            total += adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
        else:
            total += adata.X.nbytes
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        if dtype == object:
            continue
        total += adata.obs[col].values.nbytes if hasattr(adata.obs[col].values, 'nbytes') else 0
    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        if dtype == object:
            continue
        total += adata.var[col].values.nbytes if hasattr(adata.var[col].values, 'nbytes') else 0
    if hasattr(adata, 'layers'):
        for layer_name in adata.layers.keys():
            layer = adata.layers[layer_name]
            if sp.issparse(layer):
                total += layer.data.nbytes + layer.indices.nbytes + layer.indptr.nbytes
            elif layer is not None:
                total += layer.nbytes
    if hasattr(adata, 'obsm'):
        for key in adata.obsm.keys():
            arr = adata.obsm[key]
            if hasattr(arr, 'nbytes'):
                total += arr.nbytes
    if hasattr(adata, 'varm'):
        for key in adata.varm.keys():
            arr = adata.varm[key]
            if hasattr(arr, 'nbytes'):
                total += arr.nbytes
    return total / (1024 ** 3)


def downsample_by_config(adata: sc.AnnData, cfg, logger) -> sc.AnnData:
    """根据 config 设置对 adata 降采样。作为 pipeline 内联调用入口。

    Args:
        adata: 需要降采样的 AnnData。
        cfg: Fuxi Config 对象（读取 downsample_target / downsample_strategy 等）。
        logger: logging.Logger 实例。
    Returns:
        降采样后的 AnnData（如果 target 未设置或无需采样则返回原始对象）。
    """
    target = getattr(cfg, 'downsample_target', None)
    if target is None or target >= adata.n_obs:
        return adata

    strategy = getattr(cfg, 'downsample_strategy', 'stratified')
    seed = getattr(cfg, 'downsample_random_seed', 42)
    rng = np.random.RandomState(seed)
    sample_key = _check_sample_col(adata, 'sample', logger)
    n_before = adata.n_obs

    logger.info("Downsampling: target=%d, strategy=%s, seed=%d", target, strategy, seed)
    if strategy == 'random':
        adata = downsample_random(adata, target, rng, logger)
    elif strategy == 'stratified':
        adata = downsample_stratified(adata, target, sample_key or 'sample', rng, logger)
    elif strategy == 'max_per_sample':
        max_per = getattr(cfg, 'downsample_max_per_sample', 5000)
        adata = downsample_max_per_sample(adata, max_per, sample_key or 'sample', rng, logger)

    # 可选 float32 节省内存
    if getattr(cfg, 'use_float32', False):
        if sp.issparse(adata.X):
            adata.X = adata.X.astype('float32', copy=False)
        else:
            adata.X = adata.X.astype('float32')
        logger.info("X precision converted to float32")

    logger.info("Downsampled: %d → %d cells", n_before, adata.n_obs)
    return adata
