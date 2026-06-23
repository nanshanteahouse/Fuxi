#!/usr/bin/env python3
"""
run_pipeline.py — Fuxi (伏羲) 统一管线主控
===========================================

支持多组学类型:
    python run_pipeline.py --modality rna                    # scRNA-seq 全流程
    python run_pipeline.py --modality atac                   # scATAC-seq 全流程
    python run_pipeline.py --modality rna --step 3           # RNA 单步
    python run_pipeline.py --modality atac --steps 0-2       # ATAC 步骤范围
    python run_pipeline.py --modality rna --resume           # 从断点恢复
    python run_pipeline.py --modality rna --list             # 列出 RNA 步骤
    python run_pipeline.py --config my_config.py             # 自定义配置
    python run_pipeline.py --modality atac --cleanup         # 清理中间 checkpoint

用法:
    python run_pipeline.py --modality rna                      # 全部顺序执行
    python run_pipeline.py --modality atac --step 3            # 只跑第 3 步
    python run_pipeline.py --modality rna --steps 0-2          # 跑步骤 0~2
    python run_pipeline.py --modality rna --steps 1,3,5        # 跑步骤 1, 3, 5
    python run_pipeline.py --modality rna --resume             # 从第一个未完成步骤继续
    python run_pipeline.py --modality rna --list               # 列出所有步骤
    python run_pipeline.py --config my_config.py               # 使用自定义配置
"""

import sys
import os
import subprocess
import argparse
import logging
import time

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# ── Performance monitor (optional) ──────────────────────────────────
try:
    from core.utils import monitor_performance
    _HAVE_MONITOR = True
except ImportError:
    _HAVE_MONITOR = False
    from contextlib import nullcontext as _nullcontext
    def monitor_performance(*a, **kw): return _nullcontext()


# ═══════════════════════════════════════════════════════════════════════
#  RNA step registry
# ═══════════════════════════════════════════════════════════════════════
RNA_STEPS = [
    ("00", "00_load.py",                "Load raw data → 00_raw.h5ad"),
    ("01", "downsample.py",             "Downsampling (optional, config: downsample_target)"),
    ("02", "01_doublet.py",             "Scrublet doublet detection (per sample) → 01_doublet.h5ad"),
    ("03", "02_qc.py",                  "QC filtering (doublets removed) → 02_qc.h5ad"),
    ("04", "03_integrate.py",           "Normalize + HVG + PCA + Harmony → 03_integrated.h5ad"),
    ("05", "04_cluster_umap.py",        "Multi-param UMAP + multi-resolution Leiden"),
    ("06", "05_annotate_major.py",      "AI-assisted major cell type annotation (dual mode)"),
    ("07", "06_subcluster.py",          "Interactive subtype analysis (requires --cell-type)"),
    ("08", "07_markers_de.py",          "Differential expression (multi-layer)"),
    ("09", "08_trajectory.py",          "PAGA + DPT trajectory analysis"),
    ("10", "09_enrichment.py",          "GO/KEGG enrichment + AI interpretation"),
    ("11", "06_exploratory.py",         "Exploratory analysis (composition/QC/marker)"),
]

RNA_CHECKPOINT_FILES = [
    "00_raw.h5ad",               # step 00
    "00_raw.h5ad",               # step 01 (downsample overwrites)
    "01_doublet.h5ad",           # step 02
    "02_qc.h5ad",                # step 03
    "03_integrated.h5ad",        # step 04
    "04_clustered.h5ad",         # step 05
    "05_annotated.h5ad",         # step 06
    "05_annotated.h5ad",         # step 07 (reads 05_annotated)
    "05_annotated.h5ad",         # step 08 (reads 05_annotated)
    "04_clustered.h5ad",         # step 09 (reads 04_clustered)
    "marker_genes_per_group.csv",# step 10 (reads CSV from tables/)
    "05_annotated.h5ad",         # step 11 (reads 05_annotated)
]

RNA_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5, 6}


# ═══════════════════════════════════════════════════════════════════════
#  ATAC step registry
# ═══════════════════════════════════════════════════════════════════════
ATAC_STEPS = [
    ("00", "00_load.py",          "Load fragments.tsv.gz → 00_raw.h5ad"),
    ("01", "01_qc.py",            "QC + TSS enrichment + doublet detection → 01_filtered.h5ad"),
    ("02", "02_process.py",       "Peak calling + TF-IDF + spectral + KNN → 02_processed.h5ad"),
    ("03", "03_cluster.py",       "Multi-param Leiden + UMAP → 03_clustered.h5ad"),
    ("04", "04_annotate.py",      "AI-assisted chromatin state annotation → 04_annotated.h5ad"),
    ("05", "05_marker_peaks.py",  "Differential peak accessibility → marker_peaks.csv"),
    ("06", "06_motif.py",         "Motif enrichment + chromVAR → motif_results.csv"),
    ("07", "07_trajectory.py",    "ATAC pseudotime trajectory → 07_trajectory.h5ad"),
    ("08", "08_enrichment.py",    "GO/KEGG enrichment on peak-associated genes → enrichment_*.csv"),
    ("09", "09_integrate.py",     "RNA+ATAC integration via muon → 09_integrated.h5ad"),
]

ATAC_CHECKPOINT_FILES = [
    "00_raw.h5ad",           # step 00
    "01_filtered.h5ad",      # step 01
    "02_processed.h5ad",     # step 02
    "03_clustered.h5ad",     # step 03
    "04_annotated.h5ad",     # step 04
    "marker_peaks.csv",      # step 05
    "motif_results.csv",     # step 06
    "07_trajectory.h5ad",    # step 07
    "enrichment_*.csv",      # step 08 (glob)
    "09_integrated.h5ad",    # step 09
]

ATAC_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}


# ═══════════════════════════════════════════════════════════════════════
#  Modality dispatch
# ═══════════════════════════════════════════════════════════════════════
MODALITY_MAP = {
    "rna": {
        "steps": RNA_STEPS,
        "checkpoints": RNA_CHECKPOINT_FILES,
        "write_checkpoints": RNA_STEPS_WRITE_CHECKPOINT,
        "dir": "rna",
    },
    "atac": {
        "steps": ATAC_STEPS,
        "checkpoints": ATAC_CHECKPOINT_FILES,
        "write_checkpoints": ATAC_STEPS_WRITE_CHECKPOINT,
        "dir": "atac",
    },
}


def find_first_incomplete(h5ad_dir: str, steps, checkpoints, write_checkpoints, cfg=None) -> int:
    """扫描 checkpoint 目录，找到第一个未完成的步骤。"""
    if not h5ad_dir:
        logging.getLogger("run_pipeline").warning(
            "find_first_incomplete: h5ad_dir is empty (%r), falling back to current working directory",
            h5ad_dir,
        )
        h5ad_dir = "."

    for i in range(len(steps)):
        if i not in write_checkpoints:
            continue
        ckpt = os.path.join(h5ad_dir, checkpoints[i])
        if '*' in ckpt:
            import glob as glob_mod
            if not glob_mod.glob(ckpt):
                return i
        elif not os.path.exists(ckpt) or os.path.getsize(ckpt) == 0:
            return i

        # ── RNA Step 01 special: downsample may not be done ──────────
        if i == 1 and cfg is not None and hasattr(cfg, 'downsample_target') and cfg.downsample_target is not None:
            try:
                import scanpy as sc
                adata = sc.read(ckpt, backed='r')
                n_cells = adata.n_obs
                adata.file.close()
                if n_cells > cfg.downsample_target:
                    return i
            except Exception:
                pass
    return len(steps)


def _get_step_dependency(step: int, steps, checkpoints) -> str:
    """Return the checkpoint file that step `step` reads from."""
    # ATAC dependencies
    if len(steps) == 10:  # ATAC
        deps = {
            5: checkpoints[4],
            6: checkpoints[4],
            7: checkpoints[4],
            8: checkpoints[5],
            9: checkpoints[4],
        }
        return deps.get(step, checkpoints[step - 1] if step > 0 else "")
    # RNA dependencies
    deps = {
        5: checkpoints[3],
        6: checkpoints[4],
        7: checkpoints[6],
        8: checkpoints[6],
        9: checkpoints[4],
        10: checkpoints[6],
        11: checkpoints[6],
    }
    return deps.get(step, checkpoints[step - 1] if step > 0 else "")


def parse_step_range(spec: str) -> list:
    """解析步骤范围: "3-6" → [3,4,5,6], "1,3,5" → [1,3,5]"""
    if "-" in spec:
        a, b = map(int, spec.split("-"))
        return list(range(a, b + 1))
    else:
        return [int(s) for s in spec.split(",")]


def main():
    parser = argparse.ArgumentParser(
        description="Fuxi (伏羲) — Unified single-cell multi-omics pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--modality", type=str, choices=["rna", "atac", "spatial"],
                        default="rna",
                        help="Modality: rna (default), atac, spatial")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--steps", type=str,
                       help="Step range (e.g. 0-2) or list (e.g. 1,3,5)")
    group.add_argument("--step", type=int,
                       help="Run a single step (0-based)")
    group.add_argument("--resume", action="store_true",
                       help="Resume from first incomplete checkpoint")
    parser.add_argument("--list", action="store_true",
                        help="List all steps")
    parser.add_argument("--config", type=str, default="config.py",
                        help="Config file path (default: config.py)")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove upstream intermediate checkpoint files after each step")
    parser.add_argument("--cell-type", type=str,
                        help="(RNA only) Cell type to subcluster (Step 07)")
    parser.add_argument("--annotate-method", type=str,
                        choices=["auto", "unified"], default="auto",
                        help="(RNA only) Annotation method: auto=AI, unified=KB-based")
    args = parser.parse_args()

    # ── Get modality config ──────────────────────────────────────────
    if args.modality not in MODALITY_MAP:
        print(f"[run] Error: unknown modality '{args.modality}'. Supported: {list(MODALITY_MAP.keys())}")
        sys.exit(1)

    mod = MODALITY_MAP[args.modality]
    STEPS = mod["steps"]
    CHECKPOINT_FILES = mod["checkpoints"]
    STEPS_WRITE_CHECKPOINT = mod["write_checkpoints"]

    # ── --list mode ──────────────────────────────────────────────────
    if args.list:
        print(f"Fuxi — {args.modality.upper()}-seq pipeline step list")
        print("=" * 60)
        for num, script, desc in STEPS:
            ckpt = CHECKPOINT_FILES[STEPS.index((num, script, desc))]
            print(f"  [{num}] {desc}")
            print(f"        script: {script}  |  checkpoint: {ckpt}")
        print(f"\nUsage: python {os.path.basename(__file__)} --modality {args.modality} --step 3")
        return

    # ── Load config ──────────────────────────────────────────────────
    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        print(f"[run] Error: config file not found: {config_path}")
        sys.exit(1)

    import importlib.util
    spec = importlib.util.spec_from_file_location("pipeline_config", config_path)
    cfg_module = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_config"] = cfg_module
    sys.modules["config"] = cfg_module
    spec.loader.exec_module(cfg_module)
    CFG = cfg_module.CFG
    CFG.resolve_paths()

    # ── Resolve n_jobs ───────────────────────────────────────────────
    _nc = getattr(CFG, 'n_jobs', 0)
    if _nc == 0:
        _nc = os.cpu_count() or 1
        CFG.n_jobs = _nc
    print(f"[run] Using {_nc} CPU core(s)")

    # ── BLAS / OpenMP thread limits ──────────────────────────────────
    if CFG.n_jobs > 0 and getattr(CFG, 'limit_blas_threads', True):
        for var in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                     "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                     "PYTORCH_ENABLE_MPS_FALLBACK", "TORCH_NUM_THREADS"]:
            if var not in os.environ:
                os.environ[var] = str(CFG.n_jobs)
        print(f"[run] Set BLAS/OpenMP threads to {CFG.n_jobs} via env vars")

    # ── Resolve paths ────────────────────────────────────────────────
    scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', mod["dir"], 'steps')

    # ── ATAC: auto-discover RNA h5ad for Step 09 integration ─────────
    if args.modality == "atac" and not getattr(CFG, 'rna_h5ad', ''):
        from core.utils import find_rna_h5ad
        auto_rna = find_rna_h5ad(cfg=CFG)
        if auto_rna:
            CFG.rna_h5ad = auto_rna
            print(f"[run] Auto-discovered RNA h5ad for integration: {auto_rna}")
        else:
            print("[run] No RNA h5ad auto-discovered — Step 09 will be skipped.")

    python_exe = sys.executable

    # ── Parse step range ─────────────────────────────────────────────
    if args.resume:
        start = find_first_incomplete(CFG.h5ad_dir, STEPS, CHECKPOINT_FILES, STEPS_WRITE_CHECKPOINT, cfg=CFG)
        if start >= len(STEPS):
            print("[run] All steps completed.")
            return
        dep = _get_step_dependency(start, STEPS, CHECKPOINT_FILES)
        if dep:
            dep_path = os.path.join(CFG.h5ad_dir, dep)
            if '*' in dep:
                import glob as glob_mod
                if not glob_mod.glob(dep_path):
                    print(f"[run] Step [{STEPS[start][0]}] dependency missing: {dep_path}")
            elif not os.path.exists(dep_path):
                print(f"[run] Step [{STEPS[start][0]}] dependency missing: {dep_path}")
        step_indices = list(range(start, len(STEPS)))
        print(f"[run] Resuming from step [{STEPS[start][0]}]")
    elif args.steps:
        step_indices = parse_step_range(args.steps)
        for i in step_indices:
            if i < 0 or i >= len(STEPS):
                print(f"[run] Error: invalid step number {i} (valid range: 0-{len(STEPS) - 1})")
                sys.exit(1)
    elif args.step is not None:
        if args.step < 0 or args.step >= len(STEPS):
            print(f"[run] Error: step number {args.step} out of range (0-{len(STEPS) - 1})")
            sys.exit(1)
        step_indices = [args.step]
    else:
        step_indices = list(range(len(STEPS)))

    # ── Execute steps ────────────────────────────────────────────────
    step_times = []
    for i in step_indices:
        num, script, desc = STEPS[i]
        script_path = os.path.join(scripts_dir, script)

        if not os.path.exists(script_path):
            print(f"[run] Error: script not found: {script_path}")
            sys.exit(1)

        print(f"\n{'=' * 60}")
        cell_info = ""
        if args.modality == "rna" and i == 7 and args.cell_type:
            cell_info = f" (cell-type: {args.cell_type})"
        print(f"[run] [{args.modality.upper()}] Step [{num}]: {desc}{cell_info}")
        print(f"{'=' * 60}")

        extra_args = [f"--config={config_path}"]
        if args.modality == "rna" and i == 1:
            extra_args.append("--overwrite")
        if args.modality == "rna" and i == 7 and args.cell_type:
            extra_args.extend(["--cell-type", args.cell_type])

        step_t0 = time.time()
        if _HAVE_MONITOR:
            with monitor_performance(f"Step[{num}]", log=None) as perf:
                result = subprocess.run(
                    [python_exe, script_path] + extra_args,
                )
        else:
            result = subprocess.run(
                [python_exe, script_path] + extra_args,
            )
        elapsed = time.time() - step_t0

        if result.returncode != 0:
            if args.modality == "rna" and i == 7 and not args.cell_type:
                print(f"\n[run] Step [{num}] skipped (no --cell-type for pipeline mode)")
                step_times.append((num, desc, 0))
                continue
            print(f"\n[run] Step [{num}] failed (exit code={result.returncode})")
            print(f"[run] To continue after fixing the issue:")
            print(f"      python {__file__} --modality {args.modality} --resume --config {args.config}")
            sys.exit(1)

        # ── Optional checkpoint cleanup ──────────────────────────────
        if args.cleanup or getattr(CFG, 'cleanup_intermediates', False):
            dep = _get_step_dependency(i, STEPS, CHECKPOINT_FILES)
            if dep and i in STEPS_WRITE_CHECKPOINT:
                dep_path = os.path.join(CFG.h5ad_dir, dep)
                if '*' not in dep_path and os.path.exists(dep_path):
                    try:
                        os.remove(dep_path)
                        print(f"[run]   Cleaned up: {dep}")
                    except OSError:
                        pass

        print(f"[run] Step [{num}] completed (took {elapsed:.1f}s).")
        step_times.append((num, desc, elapsed))

    total_elapsed = sum(t for _, _, t in step_times)
    print(f"\n{'=' * 60}")
    print(f"[run] Fuxi {args.modality.upper()}-seq pipeline execution finished.")
    print(f"{'=' * 60}")
    print(f"[run] Step timing summary:")
    for num, desc, elapsed in step_times:
        print(f"  [{num}] {elapsed:7.1f}s  {desc}")
    print(f"  {'─' * 50}")
    print(f"  [Total] {total_elapsed:7.1f}s  {len(step_times)} steps total")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
