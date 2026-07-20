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

# Ensure repo root is on sys.path for core package imports
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

# ── Performance monitor (optional) ──────────────────────────────────
try:
    from core.utils import monitor_performance
    _HAVE_MONITOR = True
except ImportError:
    _HAVE_MONITOR = False
    from contextlib import nullcontext as _nullcontext
    def monitor_performance(step_name: str = "", log=None, child_pid=None): return _nullcontext()


# ═══════════════════════════════════════════════════════════════════════
#  RNA step registry
# ═══════════════════════════════════════════════════════════════════════
RNA_STEPS = [
    ("00", "00_load.py",                "Load raw data → 00_raw.h5ad"),
    ("01", "01_doublet.py",             "Scrublet doublet detection (per sample) → 01_doublet.h5ad"),
    ("02", "02_qc.py",                  "QC filtering (doublets removed) → 02_qc.h5ad"),
    ("03", "03_integrate.py",           "Normalize + HVG + PCA + Harmony → 03_integrated.h5ad"),
    ("04", "04_cluster_umap.py",        "Multi-param UMAP + multi-resolution Leiden"),
    ("05", "05_annotate_major.py",      "AI-assisted major cell type annotation (dual mode)"),
    ("06", "06_subcluster.py",          "Interactive subtype analysis (requires --cell-type)"),
    ("07", "07_markers_de.py",          "Differential expression (multi-layer)"),
    ("08", "08_trajectory.py",          "PAGA + DPT trajectory analysis"),
    ("09", "09_enrichment.py",          "GO/KEGG enrichment + AI interpretation"),
    ("10", "10_exploratory.py",         "Exploratory analysis (composition/QC/marker)"),
    ("11", "11_grn.py",                 "GRN regulatory network analysis (decoupler) → 11_grn.h5ad"),
    ("12", "12_cell_interaction.py",   "CCI cell-cell interaction (LIANA+) → tables + figures"),
]

RNA_CHECKPOINT_FILES = [
    "00_raw.h5ad",               # step 00
    "01_doublet.h5ad",           # step 01
    "02_qc.h5ad",                # step 02
    "03_integrated.h5ad",        # step 03
    "04_clustered.h5ad",         # step 04
    "05_annotated.h5ad",         # step 05
    "05_annotated.h5ad",         # step 06 (reads 05_annotated)
    "05_annotated.h5ad",         # step 07 (reads 05_annotated)
    "04_clustered.h5ad",         # step 08 (reads 04_clustered)
    "marker_genes_per_group.csv",# step 09 (reads CSV from tables/)
    "05_annotated.h5ad",         # step 10 (reads 05_annotated)
    "11_grn.h5ad",               # step 11
    "05_annotated.h5ad",         # step 12 (reads 05_annotated)
]

RNA_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5, 11}


# ═══════════════════════════════════════════════════════════════════════
#  ATAC step registry
# ═══════════════════════════════════════════════════════════════════════
ATAC_STEPS = [
    ("00", "00_load.py",            "Load fragments.tsv.gz → 00_raw.h5ad"),
    ("01", "01_doublet.py",          "Scrublet doublet detection → 01_doublet.h5ad"),
    ("02", "02_qc.py",               "QC filtering + TSS + peak calling + peak matrix → 02_filtered.h5ad"),
    ("03", "03_process.py",           "Feature selection + spectral + Harmony + KNN → 03_processed.h5ad"),
    ("04", "04_cluster.py",           "Multi-param Leiden + UMAP → 04_clustered.h5ad"),
    ("05", "05_peaks.py",             "Post-clustering peak calling → 05_peaks.h5ad"),
    ("06", "06_annotate.py",          "AI-assisted chromatin state annotation → 05_annotated.h5ad"),
    ("07", "07_subcluster.py",        "Subcluster analysis (placeholder)"),
    ("08", "08_marker_peaks.py",      "Differential peak accessibility → marker_peaks.csv"),
    ("09", "09_motif.py",             "Motif enrichment → motif_results.csv"),
    ("10", "10_trajectory.py",        "ATAC pseudotime trajectory → 10_trajectory.h5ad"),
    ("11", "11_enrichment.py",        "GO/KEGG enrichment on peak-associated genes → enrichment_*.csv"),
    ("12", "12_exploratory.py",       "Exploratory analysis (placeholder)"),
    ("13", "13_integrate.py",         "RNA+ATAC integration via muon → 13_integrated.h5ad"),
]

ATAC_CHECKPOINT_FILES = [
    "00_raw.h5ad",           # step 00
    "01_doublet.h5ad",       # step 01
    "02_filtered.h5ad",      # step 02
    "03_processed.h5ad",     # step 03
    "04_clustered.h5ad",     # step 04
    "05_peaks.h5ad",         # step 05
    "05_annotated.h5ad",     # step 06
    "",                      # step 07 (placeholder)
    "marker_peaks.csv",      # step 08
    "motif_results.csv",     # step 09
    "10_trajectory.h5ad",    # step 10
    "enrichment_*.csv",      # step 11
    "",                      # step 12 (placeholder)
    "13_integrated.h5ad",    # step 13
]

ATAC_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13}


# ═══════════════════════════════════════════════════════════════════════
#  Spatial step registry
# ═══════════════════════════════════════════════════════════════════════
SPATIAL_STEPS = [
    ("00", "00_load.py",           "Load spatial data -> 00_raw.h5ad (coords + image)"),
    ("01", "01_qc.py",             "QC filtering (spots + tissue detection) -> 01_qc.h5ad"),
    ("02", "02_image.py",          "Image processing (sq.im.process) -> 02_image.h5ad"),
    ("03", "03_normalize.py",      "Normalize + HVG + spatial graph -> 03_processed.h5ad"),
    ("04", "04_cluster.py",        "PCA + UMAP + Leiden clustering -> 04_clustered.h5ad"),
    ("05", "05_annotate.py",       "Cell type annotation (AI / score_genes) -> 05_annotated.h5ad"),
    ("06", "06_spatial_stats.py",     "DE + SVG + nhood enrichment + co-occurrence -> CSVs + figures"),
    ("07", "07_trajectory.py",     "Pseudotime analysis -> 07_trajectory.h5ad"),
    ("08", "08_enrichment.py",     "GO/KEGG enrichment -> enrichment CSVs"),
    ("09", "09_exploratory.py",    "Spatial visualization -> figures + CSVs"),
    ("10", "10_cell_interaction.py",  "CCI spatial cell-cell interaction (LIANA+) -> tables + figures"),
    ("11", "subcluster.py",          "Conditional subclustering per cell type -> 05_sub_{type}.h5ad"),
    ("12", "grn.py",                  "Conditional GRN analysis via decoupler -> TF activity CSV + heatmap"),
]

SPATIAL_CHECKPOINT_FILES = [
    "00_raw.h5ad",           # step 00
    "01_qc.h5ad",            # step 01
    "02_image.h5ad",         # step 02
    "03_processed.h5ad",     # step 03
    "04_clustered.h5ad",     # step 04
    "05_annotated.h5ad",     # step 05
    "05_annotated.h5ad",     # step 06
    "05_annotated.h5ad",     # step 07
    "05_annotated.h5ad",     # step 08
    "05_annotated.h5ad",     # step 09
    "05_annotated.h5ad",     # step 10
    "05_annotated.h5ad",     # step 11 (subcluster reads 05_annotated)
    "05_annotated.h5ad",     # step 12 (GRN reads 05_annotated)
]

SPATIAL_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5}


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
    "spatial": {
        "steps": SPATIAL_STEPS,
        "checkpoints": SPATIAL_CHECKPOINT_FILES,
        "write_checkpoints": SPATIAL_STEPS_WRITE_CHECKPOINT,
        "dir": "spatial",
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
    return len(steps)


def _get_step_dependency(step: int, steps, checkpoints, modality: str = "rna") -> str:
    """Return the checkpoint file that step `step` reads from."""
    if modality == "atac":
        deps = {
            2: checkpoints[0],    # 02_qc reads raw_h5ad
            3: checkpoints[1],    # 03_process reads doublet_h5ad
            6: checkpoints[4],    # 06_annotate reads clustered_h5ad
            8: checkpoints[6],    # 08_marker_peaks reads annotated_h5ad
            9: checkpoints[6],    # 09_motif reads annotated_h5ad
            10: checkpoints[6],   # 10_trajectory reads annotated_h5ad
            11: checkpoints[8],   # 11_enrichment reads marker_peaks.csv
            12: checkpoints[6],   # 12_exploratory reads annotated_h5ad
            13: checkpoints[6],   # 13_integrate reads annotated_h5ad
        }
        return deps.get(step, checkpoints[step - 1] if step > 0 else "")
    if modality == "spatial":
        deps = {
            5: checkpoints[3],   # annotate reads clustered
            6: checkpoints[5],   # spatial_de reads annotated
            7: checkpoints[5],   # trajectory reads annotated
            8: checkpoints[6],   # enrichment reads DE CSVs
            9: checkpoints[5],   # exploratory reads annotated
            10: checkpoints[5],  # CCI reads 05_annotated
            11: checkpoints[5],  # subcluster reads 05_annotated
            12: checkpoints[5],  # GRN reads 05_annotated
        }
        return deps.get(step, checkpoints[step - 1] if step > 0 else "")
    # RNA dependencies
    deps = {
        4: checkpoints[3],
        5: checkpoints[4],
        6: checkpoints[5],
        7: checkpoints[5],
        8: checkpoints[4],
        9: checkpoints[5],
        10: checkpoints[5],
        11: checkpoints[5],
        12: checkpoints[5],   # CCI reads 05_annotated
    }
    return deps.get(step, checkpoints[step - 1] if step > 0 else "")


def parse_step_range(spec: str) -> list:
    """解析步骤范围: "3-6" → [3,4,5,6], "1,3,5" → [1,3,5]"""
    if "-" in spec:
        a, b = map(int, spec.split("-"))
        return list(range(a, b + 1))
    else:
        return [int(s) for s in spec.split(",")]


def _get_checkpoint_shape(ckpt_path: str) -> tuple:
    """Fast-read h5ad shape via backed mode or CSV row count.

    Returns (n_obs, n_vars) for h5ad, (n_rows, 0) for CSV,
    or (0, 0) on any failure.
    """
    if not ckpt_path:
        return (0, 0)
    if '*' in ckpt_path:
        import glob
        files = glob.glob(ckpt_path)
        if not files:
            return (0, 0)
        ckpt_path = files[0]
    if not os.path.exists(ckpt_path):
        return (0, 0)
    try:
        if ckpt_path.endswith('.h5ad'):
            import anndata
            ad = anndata.read_h5ad(ckpt_path, backed='r')
            return (ad.n_obs, ad.n_vars)
        elif ckpt_path.endswith('.csv'):
            with open(ckpt_path) as f:
                n_rows = sum(1 for _ in f) - 1  # minus header
            return (n_rows, 0)
    except Exception:
        pass
    return (0, 0)

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
    from core.utils._config import resolve_config
    CFG = resolve_config(config_path)
    print(f"[run] Using {CFG.execution.n_jobs} CPU core(s)")

    # ── BLAS / OpenMP thread limits ──────────────────────────────────
    if CFG.execution.n_jobs > 0 and getattr(CFG.execution, 'limit_blas_threads', True):
        for var in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                     "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                     "PYTORCH_ENABLE_MPS_FALLBACK", "TORCH_NUM_THREADS"]:
            if var not in os.environ:
                os.environ[var] = str(CFG.execution.n_jobs)
        print(f"[run] Set BLAS/OpenMP threads to {CFG.execution.n_jobs} via env vars")

    # ── Resolve paths ────────────────────────────────────────────────
    scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', mod["dir"], 'steps')

    # ── ATAC: auto-discover RNA h5ad for Step 09 integration ─────────
    if args.modality == "atac" and not getattr(CFG, 'rna_h5ad', ''):
        from core.utils import find_rna_h5ad
        auto_rna = find_rna_h5ad(cfg=CFG)
        if auto_rna:
            CFG.rna_h5ad = auto_rna
            print(f"[run] Auto-discovered RNA h5ad for integration: {auto_rna}")
        else:
            print("[run] No RNA h5ad auto-discovered — Step 09 will be skipped.")

    # ── Spatial: auto-discover scRNA marker CSV for Phase 1 transfer ─
    if args.modality == "spatial" and not getattr(CFG, 'rna_ref', ''):
        from core.utils import find_rna_marker_csv
        auto_csv = find_rna_marker_csv(cfg=CFG)
        if auto_csv:
            print(f"[run] Auto-discovered scRNA marker CSV: {auto_csv}")
        else:
            print("[run] No scRNA marker CSV auto-discovered.")

    python_exe = sys.executable

    # ── Parse step range ─────────────────────────────────────────────
    if args.resume:
        start = find_first_incomplete(CFG.h5ad_dir, STEPS, CHECKPOINT_FILES, STEPS_WRITE_CHECKPOINT, cfg=CFG)
        if start >= len(STEPS):
            print("[run] All steps completed.")
            return
        dep = _get_step_dependency(start, STEPS, CHECKPOINT_FILES, modality=args.modality)
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

    if getattr(CFG, 'perf_monitoring', True):
        from core.utils import PerformanceSummary
        pipeline_summary = PerformanceSummary()
        pipeline_summary.pipeline_info = {
            "modality": args.modality,
            "config_path": config_path,
            "n_jobs": CFG.execution.n_jobs,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
    else:
        pipeline_summary = None
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
        if args.modality == "rna" and i == 6 and args.cell_type:
            cell_info = f" (cell-type: {args.cell_type})"
        print(f"[run] [{args.modality.upper()}] Step [{num}]: {desc}{cell_info}")
        print(f"{'=' * 60}")

        extra_args = [f"--config={config_path}"]
        if args.modality == "rna" and i == 6 and args.cell_type:
            extra_args.extend(["--cell-type", args.cell_type])

        step_t0 = time.time()
        step_proc = subprocess.Popen(
            [python_exe, script_path] + extra_args,
            stdout=None, stderr=None,
        )
        _perf_report = None
        if _HAVE_MONITOR and getattr(CFG, 'perf_monitoring', True):
            with monitor_performance(f"Step[{num}]", child_pid=step_proc.pid) as perf:
                step_proc.wait()
                _perf_report = perf
        else:
            step_proc.wait()
        elapsed = time.time() - step_t0
        result = step_proc

        if result.returncode == 2:
            print(f"[run] Step [{num}] skipped (no --cell-type for pipeline mode)")
            step_times.append((num, desc, 0))
            continue

        if result.returncode != 0:
            print(f"\n[run] Step [{num}] failed (exit code={result.returncode})")
            print(f"[run] To continue after fixing the issue:")
            print(f"      python {__file__} --modality {args.modality} --resume --config {args.config}")
            sys.exit(1)
        if pipeline_summary is not None and _HAVE_MONITOR and _perf_report is not None:
            # Read output checkpoint shape
            ckpt_file = CHECKPOINT_FILES[i]
            ckpt_path = os.path.join(CFG.h5ad_dir, ckpt_file) if ckpt_file else ""
            n_cells, n_genes = _get_checkpoint_shape(ckpt_path)
            if n_cells == 0 and n_genes == 0:
                dep_file = _get_step_dependency(i, STEPS, CHECKPOINT_FILES, modality=args.modality)
                if dep_file:
                    dep_path = os.path.join(CFG.h5ad_dir, dep_file)
                    if dep_path != ckpt_path:
                        n_cells, n_genes = _get_checkpoint_shape(dep_path)

            _perf_report.n_cells = n_cells
            _perf_report.n_genes = n_genes
            ckpt_size = os.path.getsize(ckpt_path) / (1024 * 1024) if ckpt_path and os.path.exists(ckpt_path) and '*' not in ckpt_path else 0.0
            _perf_report.checkpoint_mib = round(ckpt_size, 1)

            pipeline_summary.add_step(num, desc, _perf_report)

        # ── Optional checkpoint cleanup ──────────────────────────────
        if args.cleanup or getattr(CFG, 'cleanup_intermediates', False):
            dep = _get_step_dependency(i, STEPS, CHECKPOINT_FILES, modality=args.modality)
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
    if pipeline_summary is not None:
        pipeline_summary.pipeline_info["total_wall_sec"] = total_elapsed
    print(f"\n{'=' * 60}")
    print(f"[run] Fuxi {args.modality.upper()}-seq pipeline execution finished.")
    print(f"{'=' * 60}")
    print(f"[run] Step timing summary:")
    for num, desc, elapsed in step_times:
        print(f"  [{num}] {elapsed:7.1f}s  {desc}")
    print(f"  {'─' * 50}")
    print(f"  [Total] {total_elapsed:7.1f}s  {len(step_times)} steps total")
    print(f"{'=' * 60}")
    if pipeline_summary is not None:
        pipeline_summary.print_terminal_summary(
            n_jobs=CFG.execution.n_jobs,
            modality=args.modality.upper(),
            config_path=config_path,
        )
        report_path = os.path.join(CFG.results_dir, "perf_report.json")
        pipeline_summary.save_json(report_path)
        print(f"[run] Performance report saved: {report_path}")


if __name__ == "__main__":
    main()
