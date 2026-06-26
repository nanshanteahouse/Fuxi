#!/usr/bin/env python3
"""
Downsampling — CLI 包装器
===============================
核心逻辑在 core/downsample.py，本文件仅处理 CLI 参数解析和文件 I/O。

使用方式:
  # 插入在 pipeline 步骤之间:
  python rna/steps/downsample.py --config config_large.py --target-total 50000

  # 指定输入输出 checkpoint (路径从 config 读取):
  python rna/steps/downsample.py \\
      --config config_myproject.py \\
      --strategy stratified --target-total 30000

  # 超大文件使用 backed 模式（低内存读取）:
  python rna/steps/downsample.py \\
      --config config_large.py --backed --target-total 50000

  # 每个样本最多 3000 细胞:
  python rna/steps/downsample.py \\
      --config config_large.py --strategy max_per_sample --max-per-sample 3000
"""
import sys, os, time, argparse
import numpy as np
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
from core.downsample import (
    downsample_random,
    downsample_stratified,
    downsample_max_per_sample,
    estimate_memory_gb,
    _check_sample_col,
)

# 与 run_pipeline.py CHECKPOINT_FILES 同步，副本避免跨目录导入
CHECKPOINT_FILES = [
    "00_raw.h5ad", "01_doublet.h5ad", "02_qc.h5ad",
    "03_integrated.h5ad", "04_clustered.h5ad", "05_annotated.h5ad",
    "05_annotated.h5ad", "05_annotated.h5ad", "04_clustered.h5ad",
    "marker_genes_per_group.csv", "05_annotated.h5ad",
]


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
                print("[downsample] Skipped: downsample_target not configured")
                return

    # ── 设置日志 ──
    if CFG is not None:
        log_dir = CFG.log_dir
    else:
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

    if args.backed:
        log.info("  backed mode — converting to in-memory AnnData...")
        adata = adata.to_memory()

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
