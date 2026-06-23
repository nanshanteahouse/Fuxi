#!/usr/bin/env python3
"""
Downsampling — 大型 scRNA-seq 数据集的细胞降采样
===================================================
当数据集过大导致下游步骤 OOM 时，在任意 h5ad checkpoint 之间插入此脚本，
减少细胞数以控制内存使用。

支持三种降采样策略:
  1. stratified  (默认): 按样本分层采样，保持各样本比例 → 适合有 sample 列的数据
  2. random:             完全随机采样 → 适合无分组的简单降采样
  3. max_per_sample:     每个样本最多保留 N 个细胞 → 适合样本大小极不均衡的数据

使用方法:
  # 插入在 pipeline 步骤之间:
  python run_pipeline.py --step 0 --config config_large.py
  ./venv/bin/python scripts/downsample.py --config config_large.py --target-total 50000
  python run_pipeline.py --step 1 --config config_large.py

  # 指定输入输出 checkpoint (路径从 config 读取):
  ./venv/bin/python scripts/downsample.py \\
      --config config_myproject.py \\
      --strategy stratified --target-total 30000

  # 超大文件使用 backed 模式（低内存读取）:
  ./venv/bin/python scripts/downsample.py \\
      --config config_large.py --backed --target-total 50000

  # 每个样本最多 3000 细胞:
  ./venv/bin/python scripts/downsample.py \\
      --config config_large.py --strategy max_per_sample --max-per-sample 3000
"""
import sys, os, time, argparse
import numpy as np
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write

# 与 run_pipeline.py CHECKPOINT_FILES 同步，副本避免跨目录导入
CHECKPOINT_FILES = [
    "00_raw.h5ad", "01_doublet.h5ad", "02_qc.h5ad",
    "03_integrated.h5ad", "04_clustered.h5ad", "05_annotated.h5ad",
    "05_annotated.h5ad", "05_annotated.h5ad", "04_clustered.h5ad",
    "marker_genes_per_group.csv", "05_annotated.h5ad",
]


def _check_sample_col(adata: sc.AnnData, sample_key: str, log) -> str:
    """查找可用的样本分组列。返回实际使用的列名或 None。"""
    if sample_key and sample_key in adata.obs:
        return sample_key
    # 尝试常见列名
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

    # 按比例分配 target
    fractions = counts / n_cells
    per_sample_targets = (fractions * target).astype(int)
    # 处理余数 — 从余数最大的样本补 1
    remainder = target - per_sample_targets.sum()
    if remainder > 0:
        sorted_idx = np.argsort((fractions * target) - per_sample_targets)[::-1]
        for i in range(remainder):
            per_sample_targets.iloc[int(sorted_idx[i])] += 1

    # 每个样本分别采样
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
    # X matrix
    if hasattr(adata, 'X') and adata.X is not None:
        if sp.issparse(adata.X):
            # CSR: data + indices + indptr
            total += adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
        else:
            total += adata.X.nbytes
    # obs
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        if dtype == object:
            continue  # string columns are harder to estimate
        total += adata.obs[col].values.nbytes if hasattr(adata.obs[col].values, 'nbytes') else 0
    # var
    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        if dtype == object:
            continue
        total += adata.var[col].values.nbytes if hasattr(adata.var[col].values, 'nbytes') else 0
    # layers
    if hasattr(adata, 'layers'):
        for layer_name in adata.layers.keys():
            layer = adata.layers[layer_name]
            if sp.issparse(layer):
                total += layer.data.nbytes + layer.indices.nbytes + layer.indptr.nbytes
            elif layer is not None:
                total += layer.nbytes
    # obsm
    if hasattr(adata, 'obsm'):
        for key in adata.obsm.keys():
            arr = adata.obsm[key]
            if hasattr(arr, 'nbytes'):
                total += arr.nbytes
    # varm
    if hasattr(adata, 'varm'):
        for key in adata.varm.keys():
            arr = adata.varm[key]
            if hasattr(arr, 'nbytes'):
                total += arr.nbytes
    # uns (skip, too heterogeneous)
    return total / (1024 ** 3)


def main():
    t0 = time.time()

    parser = argparse.ArgumentParser(
        description="scRNA-seq cell downsampling — for pipeline use",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--config", default=None,
                        help="Config file path (auto-reads h5ad_dir etc.)")
    parser.add_argument("--input", "--input-h5ad", dest="input_h5ad", default=None,
                        help="Input h5ad path (default: checkpoint from config)")
    parser.add_argument("--output", "--output-h5ad", dest="output_h5ad", default=None,
                        help="Output h5ad path (default: overwrite input)")
    parser.add_argument("--step", type=int, default=0,
                        help="Which checkpoint step to read (default 0 = 00_raw.h5ad)")

    strategy_group = parser.add_mutually_exclusive_group()
    strategy_group.add_argument("--strategy", default="stratified",
                                choices=["random", "stratified", "max_per_sample"],
                                help="Downsampling strategy (default: stratified)")
    strategy_group.add_argument("--random", action="store_true",
                                help="Equivalent to --strategy random")

    target_group = parser.add_mutually_exclusive_group()
    target_group.add_argument("--target-total", type=int, default=None,
                              help="Target total cell count (pipeline mode reads config.downsample_target)")
    target_group.add_argument("--target-fraction", type=float, default=None,
                              help="Retention fraction (0.0-1.0)")
    target_group.add_argument("--max-per-sample", type=int, default=None,
                              help="Max cells per sample (only for max_per_sample strategy)")

    parser.add_argument("--sample-key", type=str, default="sample",
                        help="Sample name column in obs (default: sample)")
    parser.add_argument("--random-seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--backed", action="store_true",
                        help="Use backed='r' mode to read h5ad (low memory, for very large datasets)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview statistics only, do not downsample")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite input file directly (instead of writing new file)")
    parser.add_argument("--to-config", action="store_true",
                        help="Use downsampled result as input for subsequent steps (update config checkpoint references)")

    args = parser.parse_args()

    # ── 解析策略别名 ──
    strategy = args.strategy
    if args.random:
        strategy = "random"
    if strategy == "max_per_sample":
        if args.max_per_sample is None:
            parser.error("--strategy max_per_sample requires --max-per-sample")
        target_per_sample = args.max_per_sample
    elif args.max_per_sample is not None:
        parser.error("--max-per-sample only for --strategy max_per_sample")

    # ── 加载配置 ──
    CFG = None
    if args.config:
        CFG = resolve_config(args.config)
        CFG.resolve_paths()

    # ── 确定输入/输出路径 ──
    if args.input_h5ad:
        input_path = args.input_h5ad
    elif CFG is not None:
        step_idx = min(args.step, len(CHECKPOINT_FILES) - 1)
        input_path = os.path.join(CFG.h5ad_dir, CHECKPOINT_FILES[step_idx])
    else:
        parser.error("Please provide --config or --input")

    if args.output_h5ad:
        output_path = args.output_h5ad
    elif args.overwrite or args.config is None:
        output_path = input_path
    elif CFG is not None:
        # 默认生成 *_downsampled.h5ad 同级文件
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_downsampled{ext}"
    else:
        output_path = input_path

    # ── 从 config 读取降采样参数（仅当未通过 CLI 指定时） ──
    if CFG is not None:
        if args.target_total is None and args.target_fraction is None and args.max_per_sample is None:
            if CFG.downsample_target is not None:
                args.target_total = CFG.downsample_target
                if CFG.downsample_strategy:
                    strategy = CFG.downsample_strategy
                if strategy == "max_per_sample" and CFG.downsample_max_per_sample is not None:
                    args.max_per_sample = CFG.downsample_max_per_sample
                args.random_seed = CFG.downsample_random_seed
            else:
                # 配置未启用降采样 — 跳过（管道模式无操作）
                print("[downsample] Skipped: downsample_target not configured")
                return

    # ── 设置日志 ──
    if CFG is not None:
        log_dir = CFG.log_dir
    else:
        # 无 config 时，在输出文件同目录下写日志
        log_dir = os.path.join(os.path.dirname(os.path.abspath(output_path)), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log = setup_logger("downsample", os.path.join(log_dir, "downsample.log"))
    log.info("=" * 60)
    log.info("Downsampling")
    log.info("=" * 60)
    log.info("Input: %s", input_path)
    log.info("Output: %s", output_path)
    log.info("Strategy: %s", strategy)
    if args.target_total:
        log.info("Target total cells: %d", args.target_total)
    elif args.target_fraction:
        log.info("Retention fraction: %.2f", args.target_fraction)
    elif strategy == "max_per_sample":
        log.info("Max per sample: %d", args.max_per_sample)

    if not os.path.exists(input_path):
        log.error("Input file not found: %s", input_path)
        sys.exit(1)
    if args.dry_run:
        log.info("[DRY RUN] Preview only, no downsampling")

    # ── 读取数据 ──
    log.info("Reading data...")
    read_t = time.time()
    try:
        adata = sc.read(input_path, backed='r' if args.backed else None)
    except Exception as e:
        log.error("Read failed: %s", e)
        sys.exit(1)
    log.info("  Cells: %d × Genes: %d", adata.n_obs, adata.n_vars)
    log.info("  Read time: %.1fs", time.time() - read_t)

    # 如果 backed，加载到内存
    if args.backed:
        log.info("  backed mode — converting to in-memory AnnData...")
        adata = adata.to_memory()

    # 预估内存
    est_gb = estimate_memory_gb(adata)
    log.info("  Estimated memory: %.2f GB (this object only)", est_gb)

    # ── 确定采样目标 ──
    n_cells = adata.n_obs
    if args.target_fraction is not None:
        target_total = int(n_cells * args.target_fraction)
        log.info("Retention fraction %.2f → target %d cells", args.target_fraction, target_total)
    elif args.target_total is not None:
        target_total = args.target_total
    else:
        target_total = n_cells  # max_per_sample handles itself

    # ── 预览统计 ──
    if strategy != "max_per_sample":
        log.info("Current: %d cells → target: %d cells (reduction %.1f%%)",
                  n_cells, target_total,
                  (1 - target_total / n_cells) * 100 if target_total < n_cells else 0)
    else:
        log.info("Current: %d cells → max per sample: %d", n_cells, args.max_per_sample)

    sample_col = _check_sample_col(adata, args.sample_key, log)
    if sample_col:
        log.info("Cells per group:")
        counts = adata.obs[sample_col].value_counts()
        for s, c in counts.items():
            log.info("  %-30s %6d cells", s, c)

    if args.dry_run:
        log.info("[DRY RUN] Preview complete. Remove --dry-run to execute downsampling.")
        log.info("Tip: if estimated memory is insufficient, suggest target_total ≤ %d",
                 int(n_cells * 0.3))
        return

    # ── 执行降采样 ──
    rng = np.random.RandomState(args.random_seed)
    log.info("Executing downsampling (strategy=%s)...", strategy)

    if strategy == "random":
        adata = downsample_random(adata, target_total, rng, log)
    elif strategy == "stratified":
        if sample_col is None:
            log.warning("Group column not found, falling back to random sampling")
            adata = downsample_random(adata, target_total, rng, log)
        else:
            adata = downsample_stratified(adata, target_total, sample_col, rng, log)
    elif strategy == "max_per_sample":
        if sample_col is None:
            log.warning("Group column not found, falling back to random sampling")
            adata = downsample_random(adata, target_total, rng, log)
        else:
            adata = downsample_max_per_sample(adata, args.max_per_sample, sample_col, rng, log)

    # 可选 float32 节省内存
    if CFG and getattr(CFG, 'use_float32', False):
        if sp.issparse(adata.X):
            adata.X = adata.X.astype('float32', copy=False)
        else:
            adata.X = adata.X.astype('float32')
        log.info("X precision converted to float32")

    # ── 写入 ──
    log.info("Saving to: %s", output_path)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    if output_path == input_path and CFG:
        # 覆盖原始 checkpoint — 使用 safe_write
        safe_write(adata, output_path, cfg=CFG)
    else:
        # 新文件 — safe_write 是安全的
        safe_write(adata, output_path, cfg=CFG)

    # ── 摘要 ──
    elapsed = time.time() - t0
    new_est = estimate_memory_gb(adata)
    log.info("=" * 60)
    log.info("Downsampling complete!")
    log.info("  %d → %d cells (%d genes)", n_cells, adata.n_obs, adata.n_vars)
    log.info("  Estimated memory: %.2f GB → %.2f GB", est_gb, new_est)
    log.info("  Time: %.1fs (%.1fmin)", elapsed, elapsed / 60)
    log.info("  Output: %s", output_path)
    log.info("=" * 60)

    if CFG and args.to_config:
        log.info("Tip: rename downsampled file to original checkpoint name for seamless continuation")
        log.info("  mv %s %s", output_path, input_path)

    print(f"\n✅ Downsampling complete: {n_cells} → {adata.n_obs} cells")
    print(f"   Output: {output_path}")
    print(f"   Time: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
