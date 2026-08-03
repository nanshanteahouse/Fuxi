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

import argparse
import cProfile
import logging
import os
import runpy
import subprocess
import sys
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
from core.utils import _set_blas_env  # noqa: E402

# ── Performance monitor (optional) ──────────────────────────────────
try:
    from core.utils import monitor_performance

    _HAVE_MONITOR = True
except ImportError:
    _HAVE_MONITOR = False
    from contextlib import nullcontext as _nullcontext

    def monitor_performance(step_name: str = "", log=None, child_pid=None):
        return _nullcontext()


# ═══════════════════════════════════════════════════════════════════════
#  RNA step registry
# ═══════════════════════════════════════════════════════════════════════
RNA_STEPS = [
    ("00", "00_load.py", "Load raw data → 00_raw.h5ad"),
    ("01", "01_doublet.py", "Scrublet doublet detection (per sample) → 01_doublet.h5ad"),
    ("02", "02_qc.py", "QC filtering (doublets removed) → 02_qc.h5ad"),
    ("03", "03_integrate.py", "Normalize + HVG + PCA + batch integration → 03_integrated.h5ad"),
    ("04", "04_cluster_umap.py", "Multi-param UMAP + multi-resolution Leiden"),
    ("05", "05_annotate_major.py", "AI-assisted major cell type annotation (dual mode)"),
    ("06", "06_subcluster.py", "Interactive subtype analysis (requires --cell-type)"),
    ("07", "07_markers_de.py", "Differential expression (multi-layer)"),
    ("08", "08_trajectory.py", "Trajectory analysis (PAGA/DPT or scVelo)"),
    ("09", "09_enrichment.py", "GO/KEGG enrichment + AI interpretation"),
    ("10", "10_exploratory.py", "Exploratory analysis (composition/QC/marker)"),
    ("11", "11_grn.py", "GRN regulatory network analysis (decoupler) → 11_grn.h5ad"),
    ("12", "12_cell_interaction.py", "CCI cell-cell interaction (LIANA+) → tables + figures"),
]

RNA_CHECKPOINT_FILES = [
    "00_raw.h5ad",  # step 00
    "01_doublet.h5ad",  # step 01
    "02_qc.h5ad",  # step 02
    "03_integrated.h5ad",  # step 03
    "04_clustered.h5ad",  # step 04
    "05_annotated.h5ad",  # step 05
    "05_annotated.h5ad",  # step 06 (reads 05_annotated)
    "05_annotated.h5ad",  # step 07 (reads 05_annotated)
    "04_clustered.h5ad",  # step 08 (reads 04_clustered)
    "marker_genes_per_group.csv",  # step 09 (reads CSV from tables/)
    "05_annotated.h5ad",  # step 10 (reads 05_annotated)
    "11_grn.h5ad",  # step 11
    "05_annotated.h5ad",  # step 12 (reads 05_annotated)
]

RNA_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5, 11}

# ── Sentinel 完成度标记 (plan h5ad-incremental-io Item 1.6) ─────────
# 原地写回 / 可能无产物 的步骤不产出新 checkpoint 文件，其锚定文件（如
# 05_annotated.h5ad）在更早步骤就已存在，不能作为完成标志。改用独立
# sentinel 文件标记步骤完成。
#
# 命名约定: "<base_checkpoint>.step{NN}_done"，与 base checkpoint 放在
# 同一 h5ad 目录（如 05_annotated.h5ad.step06_done）。
#
# 写入方: 步骤脚本自身。runner 以 subprocess 运行步骤、无法访问步骤内部
# 状态，因此 sentinel 由步骤脚本在成功路径的最后一步创建（内容必须非空，
# 建议 JSON 或 "done" 文本；空文件视为未完成）。T5/T6 按此约定在
# 06_subcluster.py / 08_trajectory.py 中写 sentinel。
#
# 读取方: find_first_incomplete() 对 sentinel 步骤只看 sentinel 存在性
# （替代 checkpoint 检查）。
# 清理方: --cleanup 删除锚定文件时同步删除其 sentinel（_remove_anchored_sentinels）。
RNA_SENTINEL_FILES = {
    # step 06 (06_subcluster.py): 原地写回 05_annotated.h5ad，不产新文件
    # → 以 sentinel 证明完成。注意：无 --cell-type 且未配置 subcluster_types
    # 时步骤以 exit 2 跳过、不写 sentinel（resume 会从 step 06 重新开始并
    # 再次跳过，属预期行为）。
    6: "05_annotated.h5ad.step06_done",
    # step 08 (08_trajectory.py): CFG.trajectory.save_final_h5ad=False 时不产
    # 05_final.h5ad → 以 sentinel 证明完成（PAGA/scVelo 分支都要写）。
    8: "05_final.h5ad.step08_done",
}


# ═══════════════════════════════════════════════════════════════════════
#  ATAC step registry
# ═══════════════════════════════════════════════════════════════════════
ATAC_STEPS = [
    ("00", "00_load.py", "Load fragments.tsv.gz → 00_raw.h5ad"),
    ("01", "01_doublet.py", "Scrublet doublet detection → 01_doublet.h5ad"),
    ("02", "02_qc.py", "QC filtering + TSS + peak calling + peak matrix → 02_filtered.h5ad"),
    (
        "03",
        "03_process.py",
        "Feature selection + spectral + batch correction + KNN → 03_processed.h5ad",
    ),
    ("04", "04_cluster.py", "Multi-param Leiden + UMAP → 04_clustered.h5ad"),
    ("05", "05_peaks.py", "Post-clustering peak calling → 05_peaks.h5ad"),
    ("06", "06_annotate.py", "AI-assisted chromatin state annotation → 05_annotated.h5ad"),
    ("07", "07_subcluster.py", "Subcluster analysis (placeholder)"),
    ("08", "08_marker_peaks.py", "Differential peak accessibility → marker_peaks.csv"),
    ("09", "09_motif.py", "Motif enrichment → motif_results.csv"),
    ("10", "10_trajectory.py", "ATAC pseudotime trajectory → 10_trajectory.h5ad"),
    ("11", "11_enrichment.py", "GO/KEGG enrichment on peak-associated genes → enrichment_*.csv"),
    ("12", "12_exploratory.py", "Exploratory analysis (placeholder)"),
    ("13", "13_integrate.py", "RNA+ATAC integration via muon → 13_integrated.h5ad"),
]

ATAC_CHECKPOINT_FILES = [
    "00_raw.h5ad",  # step 00
    "01_doublet.h5ad",  # step 01
    "02_filtered.h5ad",  # step 02
    "03_processed.h5ad",  # step 03
    "04_clustered.h5ad",  # step 04
    "05_peaks.h5ad",  # step 05
    "05_annotated.h5ad",  # step 06
    "",  # step 07 (placeholder)
    "marker_peaks.csv",  # step 08
    "motif_results.csv",  # step 09
    "10_trajectory.h5ad",  # step 10
    "enrichment_*.csv",  # step 11
    "",  # step 12 (placeholder)
    "13_integrated.h5ad",  # step 13
]

ATAC_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13}

# Sentinel 完成度标记 — ATAC 暂无原地写回/无产物步骤（后续推广时按 RNA 约定补充）。
ATAC_SENTINEL_FILES: dict[int, str] = {}


# ═══════════════════════════════════════════════════════════════════════
#  Spatial step registry
# ═══════════════════════════════════════════════════════════════════════
SPATIAL_STEPS = [
    ("00", "00_load.py", "Load spatial data -> 00_raw.h5ad (coords + image)"),
    ("01", "01_qc.py", "QC filtering (spots + tissue detection) -> 01_qc.h5ad"),
    ("02", "02_image.py", "Image processing (sq.im.process) -> 02_image.h5ad"),
    ("03", "03_normalize.py", "Normalize + HVG + spatial graph -> 03_processed.h5ad"),
    ("04", "04_cluster.py", "PCA + UMAP + Leiden clustering -> 04_clustered.h5ad"),
    ("05", "05_annotate.py", "Cell type annotation (AI / score_genes) -> 05_annotated.h5ad"),
    ("06", "06_spatial_stats.py", "DE + SVG + nhood enrichment + co-occurrence -> CSVs + figures"),
    ("07", "07_trajectory.py", "Pseudotime analysis -> 07_trajectory.h5ad"),
    ("08", "08_enrichment.py", "GO/KEGG enrichment -> enrichment CSVs"),
    ("09", "09_exploratory.py", "Spatial visualization -> figures + CSVs"),
    (
        "10",
        "10_cell_interaction.py",
        "CCI spatial cell-cell interaction (LIANA+) -> tables + figures",
    ),
    ("11", "subcluster.py", "Conditional subclustering per cell type -> 05_sub_{type}.h5ad"),
    ("12", "grn.py", "Conditional GRN analysis via decoupler -> TF activity CSV + heatmap"),
]

SPATIAL_CHECKPOINT_FILES = [
    "00_raw.h5ad",  # step 00
    "01_qc.h5ad",  # step 01
    "02_image.h5ad",  # step 02
    "03_processed.h5ad",  # step 03
    "04_clustered.h5ad",  # step 04
    "05_annotated.h5ad",  # step 05
    "05_annotated.h5ad",  # step 06
    "05_annotated.h5ad",  # step 07
    "05_annotated.h5ad",  # step 08
    "05_annotated.h5ad",  # step 09
    "05_annotated.h5ad",  # step 10
    "05_annotated.h5ad",  # step 11 (subcluster reads 05_annotated)
    "05_annotated.h5ad",  # step 12 (GRN reads 05_annotated)
]

SPATIAL_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 3, 4, 5}

# Sentinel 完成度标记 — Spatial 暂无（后续推广时按 RNA 约定补充）。
SPATIAL_SENTINEL_FILES: dict[int, str] = {}

# ═══════════════════════════════════════════════════════════════════════
#  Bulk step registry
# ═══════════════════════════════════════════════════════════════════════
BULK_STEPS = [
    ("00", "00_load.py", "Load count matrix (CSV/TSV/h5ad) -> 00_raw.h5ad"),
    ("01", "01_qc.py", "Sample QC (library size, gene detection) -> 01_qc.h5ad"),
    ("02", "02_de.py", "DESeq2 normalization + DE -> 02_de.h5ad + CSVs + figures"),
    ("03", "03_enrichment.py", "GO/KEGG enrichment (GSEApy) -> tables/"),
    ("04", "04_exploratory.py", "PCA, heatmaps, volcano plots -> figures/"),
    ("05", "05_batch.py", "Batch correction (optional, pycombat) -> 05_batch_corrected.h5ad"),
]

BULK_CHECKPOINT_FILES = [
    "00_raw.h5ad",  # step 00
    "01_qc.h5ad",  # step 01
    "02_de.h5ad",  # step 02
    "",  # step 03 (CSV output, no h5ad checkpoint)
    "",  # step 04 (figures output)
    "05_batch_corrected.h5ad",  # step 05 (optional)
]

BULK_STEPS_WRITE_CHECKPOINT = {0, 1, 2, 5}

# Sentinel 完成度标记 — Bulk 暂无（后续推广时按 RNA 约定补充）。
BULK_SENTINEL_FILES: dict[int, str] = {}

# ═══════════════════════════════════════════════════════════════════════
#  Modality dispatch
# ═══════════════════════════════════════════════════════════════════════
MODALITY_MAP = {
    "rna": {
        "steps": RNA_STEPS,
        "checkpoints": RNA_CHECKPOINT_FILES,
        "write_checkpoints": RNA_STEPS_WRITE_CHECKPOINT,
        "sentinels": RNA_SENTINEL_FILES,
        "dir": "rna",
    },
    "atac": {
        "steps": ATAC_STEPS,
        "checkpoints": ATAC_CHECKPOINT_FILES,
        "write_checkpoints": ATAC_STEPS_WRITE_CHECKPOINT,
        "sentinels": ATAC_SENTINEL_FILES,
        "dir": "atac",
    },
    "spatial": {
        "steps": SPATIAL_STEPS,
        "checkpoints": SPATIAL_CHECKPOINT_FILES,
        "write_checkpoints": SPATIAL_STEPS_WRITE_CHECKPOINT,
        "sentinels": SPATIAL_SENTINEL_FILES,
        "dir": "spatial",
    },
    "bulk": {
        "steps": BULK_STEPS,
        "checkpoints": BULK_CHECKPOINT_FILES,
        "write_checkpoints": BULK_STEPS_WRITE_CHECKPOINT,
        "sentinels": BULK_SENTINEL_FILES,
        "dir": "bulk",
    },
}


def find_first_incomplete(
    h5ad_dir: str,
    steps,
    checkpoints,
    write_checkpoints,
    sentinels=None,
    cfg=None,
) -> int:
    """扫描 checkpoint 目录，找到第一个未完成的步骤。

    sentinel 步骤（见 ``*_SENTINEL_FILES``，如 RNA 的 step 06/08）的完成度
    只看 sentinel 文件是否存在且非空——其锚定 checkpoint（如 05_annotated.h5ad）
    在更早步骤就已存在，不能作为完成标志。其余步骤保持原有 checkpoint 判定。
    """
    if not h5ad_dir:
        logging.getLogger("run_pipeline").warning(
            "find_first_incomplete: h5ad_dir is empty (%r), falling back to current working directory",
            h5ad_dir,
        )
        h5ad_dir = "."

    sentinels = sentinels or {}

    for i in range(len(steps)):
        if i in sentinels:
            sent = os.path.join(h5ad_dir, sentinels[i])
            if not os.path.exists(sent) or os.path.getsize(sent) == 0:
                return i
            continue
        if i not in write_checkpoints:
            continue
        ckpt = os.path.join(h5ad_dir, checkpoints[i])
        if "*" in ckpt:
            import glob as glob_mod

            if not glob_mod.glob(ckpt):
                return i
        elif not os.path.exists(ckpt) or os.path.getsize(ckpt) == 0:
            return i
    return len(steps)


def _sentinel_base(sentinel_file: str) -> str:
    """从 sentinel 文件名还原其锚定的 base checkpoint 文件名。

    命名约定: ``"<base_checkpoint>.step{NN}_done"`` → ``"<base_checkpoint>"``
    （如 ``05_annotated.h5ad.step06_done`` → ``05_annotated.h5ad``）。
    """
    return sentinel_file.rsplit(".step", 1)[0]


def _remove_anchored_sentinels(h5ad_dir: str, base_file: str, sentinels) -> None:
    """删除锚定在 ``base_file`` 上的 sentinel 文件（``--cleanup`` 用）。

    base checkpoint 被清理后其 sentinel 失去意义，必须同步删除，避免
    ``--resume`` 依据 stale sentinel 误判步骤已完成。
    """
    if not sentinels:
        return
    for sent_file in sentinels.values():
        if _sentinel_base(sent_file) != base_file:
            continue
        sent_path = os.path.join(h5ad_dir, sent_file)
        try:
            if os.path.exists(sent_path):
                os.remove(sent_path)
                print(f"[run]   Cleaned up: {sent_file}")
        except OSError:
            pass


def _get_step_dependency(step: int, steps, checkpoints, modality: str = "rna") -> str:
    """Return the checkpoint file that step `step` reads from."""
    if modality == "atac":
        deps = {
            2: checkpoints[0],  # 02_qc reads raw_h5ad
            3: checkpoints[1],  # 03_process reads doublet_h5ad
            6: checkpoints[4],  # 06_annotate reads clustered_h5ad
            8: checkpoints[6],  # 08_marker_peaks reads annotated_h5ad
            9: checkpoints[6],  # 09_motif reads annotated_h5ad
            10: checkpoints[6],  # 10_trajectory reads annotated_h5ad
            11: checkpoints[8],  # 11_enrichment reads marker_peaks.csv
            12: checkpoints[6],  # 12_exploratory reads annotated_h5ad
            13: checkpoints[6],  # 13_integrate reads annotated_h5ad
        }
        return deps.get(step, checkpoints[step - 1] if step > 0 else "")
    if modality == "spatial":
        deps = {
            5: checkpoints[3],  # annotate reads clustered
            6: checkpoints[5],  # spatial_de reads annotated
            7: checkpoints[5],  # trajectory reads annotated
            8: checkpoints[6],  # enrichment reads DE CSVs
            9: checkpoints[5],  # exploratory reads annotated
            10: checkpoints[5],  # CCI reads 05_annotated
            11: checkpoints[5],  # subcluster reads 05_annotated
            12: checkpoints[5],  # GRN reads 05_annotated
        }
        return deps.get(step, checkpoints[step - 1] if step > 0 else "")
    if modality == "bulk":
        deps = {
            1: checkpoints[0],  # 01_qc reads raw
            2: checkpoints[1],  # 02_de reads qc
            3: checkpoints[2],  # 03_enrichment reads DE output
            4: checkpoints[2],  # 04_exploratory reads DE output
            5: checkpoints[1],  # 05_batch reads qc
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
        12: checkpoints[5],  # CCI reads 05_annotated
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
    if "*" in ckpt_path:
        import glob

        files = glob.glob(ckpt_path)
        if not files:
            return (0, 0)
        ckpt_path = files[0]
    if not os.path.exists(ckpt_path):
        return (0, 0)
    try:
        if ckpt_path.endswith(".h5ad"):
            import anndata

            ad = anndata.read_h5ad(ckpt_path, backed="r")
            return (ad.n_obs, ad.n_vars)
        elif ckpt_path.endswith(".csv"):
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
    parser.add_argument(
        "--modality",
        type=str,
        choices=["rna", "atac", "spatial", "bulk"],
        default="rna",
        help="Modality: rna (default), atac, spatial, bulk",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--steps", type=str, help="Step range (e.g. 0-2) or list (e.g. 1,3,5)")
    group.add_argument("--step", type=int, help="Run a single step (0-based)")
    group.add_argument(
        "--resume", action="store_true", help="Resume from first incomplete checkpoint"
    )
    parser.add_argument("--list", action="store_true", help="List all steps")
    parser.add_argument(
        "--config", type=str, default="config.py", help="Config file path (default: config.py)"
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove upstream intermediate checkpoint files after each step",
    )
    parser.add_argument(
        "--cell-type", type=str, help="(RNA only) Cell type to subcluster (Step 07)"
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        default=False,
        help="(Step 04 only) Re-render summary figures from the saved checkpoint without re-clustering",
    )
    parser.add_argument(
        "--annotate-method",
        type=str,
        choices=["auto", "unified"],
        default="auto",
        help="(RNA only) Annotation method: auto=AI, unified=KB-based",
    )
    parser.add_argument(
        "--in-process",
        action="store_true",
        default=False,
        help="Run steps in-process via runpy (debug mode: avoids Python/scanpy re-import overhead)",
    )
    parser.add_argument(
        "--profile",
        type=str,
        default=None,
        help="Enable cProfile profiling; path to .prof file or directory for per-step profiles",
    )
    args = parser.parse_args()

    # ── Get modality config ──────────────────────────────────────────
    if args.modality not in MODALITY_MAP:
        print(
            f"[run] Error: unknown modality '{args.modality}'. Supported: {list(MODALITY_MAP.keys())}"
        )
        sys.exit(1)

    mod = MODALITY_MAP[args.modality]
    STEPS = mod["steps"]  # noqa: N806
    CHECKPOINT_FILES = mod["checkpoints"]  # noqa: N806
    STEPS_WRITE_CHECKPOINT = mod["write_checkpoints"]  # noqa: N806

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

    CFG = resolve_config(config_path)  # noqa: N806
    print(f"[run] Using {CFG.execution.n_jobs} CPU core(s)")

    # ── BLAS / OpenMP thread limits ──────────────────────────────────
    if CFG.execution.n_jobs > 0 and getattr(CFG.execution, "limit_blas_threads", True):
        _set_blas_env(CFG.execution.n_jobs, overwrite=False)
        print(f"[run] Set BLAS/OpenMP threads to {CFG.execution.n_jobs} via env vars")

    # ── Resolve paths ────────────────────────────────────────────────
    scripts_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", mod["dir"], "steps"
    )

    # ── ATAC: auto-discover RNA h5ad for Step 09 integration ─────────
    if args.modality == "atac" and not getattr(CFG, "rna_h5ad", ""):
        from core.utils import find_rna_h5ad

        auto_rna = find_rna_h5ad(cfg=CFG)
        if auto_rna:
            CFG.rna_h5ad = auto_rna
            print(f"[run] Auto-discovered RNA h5ad for integration: {auto_rna}")
        else:
            print("[run] No RNA h5ad auto-discovered — Step 09 will be skipped.")

    # ── Spatial: auto-discover scRNA marker CSV for Phase 1 transfer ─
    if args.modality == "spatial" and not getattr(CFG, "rna_ref", ""):
        from core.utils import find_rna_marker_csv

        auto_csv = find_rna_marker_csv(cfg=CFG)
        if auto_csv:
            print(f"[run] Auto-discovered scRNA marker CSV: {auto_csv}")
        else:
            print("[run] No scRNA marker CSV auto-discovered.")

    python_exe = sys.executable

    # ── Parse step range ─────────────────────────────────────────────
    if args.resume:
        start = find_first_incomplete(
            CFG.h5ad_dir,
            STEPS,
            CHECKPOINT_FILES,
            STEPS_WRITE_CHECKPOINT,
            sentinels=mod["sentinels"],
            cfg=CFG,
        )
        if start >= len(STEPS):
            print("[run] All steps completed.")
            return
        dep = _get_step_dependency(start, STEPS, CHECKPOINT_FILES, modality=args.modality)
        if dep:
            dep_path = os.path.join(CFG.h5ad_dir, dep)
            if "*" in dep:
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

    if getattr(CFG, "perf_monitoring", True):
        from core.utils import PerformanceSummary

        # Preserve historical perf data across --step N / --resume runs:
        # load any existing perf_report.json so subsequent add_step() calls
        # upsert into history instead of overwriting it.
        _perf_report_path = os.path.join(CFG.results_dir, "perf_report.json")
        pipeline_summary = PerformanceSummary.load_existing(_perf_report_path)
        if pipeline_summary is None:
            pipeline_summary = PerformanceSummary()
        _prev_info = pipeline_summary.pipeline_info or {}
        pipeline_summary.pipeline_info = {
            "modality": args.modality,
            "config_path": config_path,
            "n_jobs": CFG.execution.n_jobs,
            "first_run_timestamp": _prev_info.get("first_run_timestamp")
            or _prev_info.get("timestamp")
            or time.strftime("%Y-%m-%dT%H:%M:%S"),
            "last_run_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "partial": True,  # set to False at end of full run
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
        if getattr(args, "plot_only", False) and num == "04":
            extra_args.append("--plot-only")

        step_t0 = time.time()

        # ── Profiling setup ─────────────────────────────────────────
        # The main process only waits on subprocesses, so a profiler here
        # would record nothing but subprocess.wait. For subprocess steps we
        # inject `python -m cProfile -o <path>` into the child command so the
        # step's own code is profiled; the in-process (debug) path keeps an
        # in-process profiler instead.
        _prof_path = None
        if args.profile:
            if args.profile.endswith(".prof"):
                _prof_path = args.profile
            else:
                _prof_path = os.path.join(args.profile, f"step-{num}.prof")
                os.makedirs(args.profile, exist_ok=True)
        if args.profile and args.in_process:
            _profiler = cProfile.Profile()
            _profiler.enable()
        if args.in_process:
            # ── In-process execution (debug mode) ──────────────────
            # Uses runpy.run_path to avoid Python+scanpy re-import
            # overhead during development. Not for production: the
            # step script shares the same process namespace, which
            # can cause state leakage between steps.
            _perf_report = None
            _was_interrupted = False
            _saved_argv = sys.argv
            sys.argv = [script_path] + extra_args
            try:
                runpy.run_path(script_path, run_name="__main__")
                _exit_code = 0
            except SystemExit as e:
                _exit_code = e.code if e.code is not None else 0
            except KeyboardInterrupt:
                _was_interrupted = True
                _exit_code = 130
            finally:
                sys.argv = _saved_argv
            elapsed = time.time() - step_t0
            result = type("_ProcResult", (), {"returncode": _exit_code})()
        else:
            _child_cmd = [python_exe, script_path] + extra_args
            if _prof_path:
                _child_cmd = [
                    python_exe,
                    "-m",
                    "cProfile",
                    "-o",
                    _prof_path,
                    script_path,
                ] + extra_args
            step_proc = subprocess.Popen(
                _child_cmd,
                stdout=None,
                stderr=None,
            )
            _perf_report = None
            _was_interrupted = False
            try:
                if _HAVE_MONITOR and getattr(CFG, "perf_monitoring", True):
                    with monitor_performance(f"Step[{num}]", child_pid=step_proc.pid) as perf:
                        step_proc.wait()
                        _perf_report = perf
                else:
                    step_proc.wait()
            except (KeyboardInterrupt, SystemExit):
                # SIGTERM/SIGINT/Ctrl+C — kill child if still alive, then fall through
                # to the normal add_step path so partial perf is recorded.
                _was_interrupted = True
                try:
                    step_proc.kill()
                    step_proc.wait(timeout=5)
                except Exception:
                    pass
            elapsed = time.time() - step_t0
            result = step_proc

        # ── Profiling teardown ──────────────────────────────────────
        if args.profile and args.in_process:
            _profiler.disable()
            if args.profile.endswith(".prof"):
                _prof_path = args.profile
            else:
                _prof_path = os.path.join(args.profile, f"step-{num}.prof")
            _profiler.dump_stats(_prof_path)

        # Determine exit status for perf_report
        if _was_interrupted:
            _exit_status = "killed"
        elif result.returncode == 2:
            _exit_status = "skipped"
        elif result.returncode != 0:
            _exit_status = "failed"
        else:
            _exit_status = "completed"

        # Always record + persist immediately (even on kill/failure) so perf_report
        # survives interruption — single-step / SIGTERM no longer lose history.
        if pipeline_summary is not None and _HAVE_MONITOR and _perf_report is not None:
            _perf_report.exit_status = _exit_status
            # Read checkpoint shape (only meaningful when step completed)
            if _exit_status == "completed":
                ckpt_file = CHECKPOINT_FILES[i]
                ckpt_path = os.path.join(CFG.h5ad_dir, ckpt_file) if ckpt_file else ""
                n_cells, n_genes = _get_checkpoint_shape(ckpt_path)
                if n_cells == 0 and n_genes == 0:
                    dep_file = _get_step_dependency(
                        i, STEPS, CHECKPOINT_FILES, modality=args.modality
                    )
                    if dep_file:
                        dep_path = os.path.join(CFG.h5ad_dir, dep_file)
                        if dep_path != ckpt_path:
                            n_cells, n_genes = _get_checkpoint_shape(dep_path)
                _perf_report.n_cells = n_cells
                _perf_report.n_genes = n_genes
                ckpt_size = (
                    os.path.getsize(ckpt_path) / (1024 * 1024)
                    if ckpt_path and os.path.exists(ckpt_path) and "*" not in ckpt_path
                    else 0.0
                )
                _perf_report.checkpoint_mib = round(ckpt_size, 1)

            pipeline_summary.add_step(num, desc, _perf_report)
            # Per-step persistence: save perf_report after every step so
            # single-step runs and interrupted runs both leave a usable trail.
            try:
                pipeline_summary.save_json(_perf_report_path)
            except Exception as _save_err:
                print(f"[run] Warning: failed to save perf_report: {_save_err}")

        # ── Exit handling per status ───────────────────────────────
        if _exit_status == "skipped":
            print(f"[run] Step [{num}] skipped (no --cell-type for pipeline mode)")
            step_times.append((num, desc, 0))
            continue
        if _was_interrupted:
            print(f"\n[run] Step [{num}] killed (interrupted by user)")
            sys.exit(130)
        if result.returncode != 0:
            print(f"\n[run] Step [{num}] failed (exit code={result.returncode})")
            print("[run] To continue after fixing the issue:")
            print(
                f"      python {__file__} --modality {args.modality} --resume --config {args.config}"
            )
            sys.exit(1)

        # ── Optional checkpoint cleanup ──────────────────────────────
        if args.cleanup or getattr(CFG, "cleanup_intermediates", False):
            dep = _get_step_dependency(i, STEPS, CHECKPOINT_FILES, modality=args.modality)
            if dep and i in STEPS_WRITE_CHECKPOINT:
                dep_path = os.path.join(CFG.h5ad_dir, dep)
                if "*" not in dep_path and os.path.exists(dep_path):
                    try:
                        os.remove(dep_path)
                        print(f"[run]   Cleaned up: {dep}")
                    except OSError:
                        pass
                    # 锚定文件被清理 → 同步删除其 sentinel（如有），
                    # 避免 stale 完成标记误导 --resume。
                    _remove_anchored_sentinels(CFG.h5ad_dir, dep, mod["sentinels"])

        print(f"[run] Step [{num}] completed (took {elapsed:.1f}s).")
        step_times.append((num, desc, elapsed))

    total_elapsed = sum(t for _, _, t in step_times)
    if pipeline_summary is not None:
        pipeline_summary.pipeline_info["total_wall_sec"] = total_elapsed
        # Mark pipeline as complete (not partial) only if we reached the end normally.
        # Single-step / interrupted runs keep partial=True (set at init).
        pipeline_summary.pipeline_info["partial"] = False
    print(f"\n{'=' * 60}")
    print(f"[run] Fuxi {args.modality.upper()}-seq pipeline execution finished.")
    print(f"{'=' * 60}")
    print("[run] Step timing summary:")
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
