#!/usr/bin/env python3
"""
utils.py — Fuxi 管线通用工具函数
===================================

合并 scRNAseq_pipeline 和 ATACseq_pipeline 的共享工具函数:
  - safe_write:       WSL 兼容的 h5ad 安全保存
  - safe_plot:        容错的 matplotlib 绘图包装
  - setup_logger:     统一日志配置
  - resolve_config:   动态加载项目配置文件
  - validate_adata:   检查 AnnData X 矩阵完整性
  - monitor_performance: 性能监控上下文管理器 (来自 ATACseq)
  - is_wsl / data_root / repo_root / wsl_to_win: 跨平台路径工具 (来自 scRNAseq)
"""

import os
import sys
import shutil
import logging
import warnings
import platform
import time as _time
import threading
import subprocess as _sp
from pathlib import Path
from typing import Optional
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
import scanpy as sc
import scipy.sparse as sp


# ── WSL h5py file locking auto-detection ────────────────────────────
# WSL /mnt mounts don't support fcntl locking; silence h5py errors.
if 'microsoft' in platform.release().lower():
    os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')


# ── Auto-load .env ──────────────────────────────────────────────────
# Standalone step scripts (python rna/steps/00_load.py --config=...)
# call resolve_config() but don't load the .env file beforehand.
# run_pipeline.py already does this, but steps spawned via subprocess
# inherit the parent env anyway — the real fix is for direct runs.
try:
    from dotenv import load_dotenv
    _dotenv_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        '.env',
    )
    if os.path.exists(_dotenv_path):
        load_dotenv(_dotenv_path, override=True)
except ImportError:
    pass


# ── Cross-platform path helpers ─────────────────────────────────────
_DATA_ROOT_CACHE: Optional[str] = None
_REPO_ROOT_CACHE: Optional[str] = None


def is_wsl() -> bool:
    """True when running inside Windows Subsystem for Linux."""
    return ('microsoft' in platform.uname().release.lower()
            and os.path.exists('/mnt/c'))


def data_root() -> str:
    """Absolute root of the raw data tree.

    Resolved in order of precedence:
      1. FUXI_DATA_ROOT  env var (canonical name)
      2. SCRNA_DATA_ROOT env var (legacy name, backward compat)

    If neither is set the function raises RuntimeError with setup
    instructions — this is intentional: every machine has its own
    data layout, so the path must be configured explicitly.

    Cached after first resolution.
    """
    global _DATA_ROOT_CACHE
    if _DATA_ROOT_CACHE is None:
        _DATA_ROOT_CACHE = (
            os.environ.get('FUXI_DATA_ROOT')
            or os.environ.get('SCRNA_DATA_ROOT')
        )
        if not _DATA_ROOT_CACHE:
            raise RuntimeError(
                "Data root not configured.\n"
                "  Set the FUXI_DATA_ROOT environment variable to the\n"
                "  directory containing your GEO dataset folders, e.g.:\n"
                '    export FUXI_DATA_ROOT=/mnt/e/data   # WSL\n'
                '    set FUXI_DATA_ROOT=E:/data          # Windows'
            )
    return _DATA_ROOT_CACHE


def repo_root() -> str:
    """Absolute path to this repository's root, located from __file__.

    Override (rare) via SCRNA_REPO_ROOT. Cached at import.
    """
    global _REPO_ROOT_CACHE
    if _REPO_ROOT_CACHE is None:
        _REPO_ROOT_CACHE = os.environ.get('SCRNA_REPO_ROOT') or os.path.abspath(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
    return _REPO_ROOT_CACHE


def wsl_to_win(path: str) -> str:
    """Translate /mnt/X/... -> X:/...; pass-through if not /mnt-prefixed."""
    if path.startswith('/mnt/') and len(path) > 5:
        return f"{path[5]}:/{path[6:]}"
    return path


# ── Core utilities ──────────────────────────────────────────────────

def safe_write(adata, target: str,
               tmpdir: str = "/tmp/Fuxi",
               compression: str = "gzip", cfg=None,
               compression_override: Optional[str] = None) -> None:
    """
    安全写入 h5ad 文件，避免 WSL /mnt 挂载的文件锁定问题。

    策略: 先写入 tmpdir，再 mv 到目标路径。
    mv 是原子操作（在同一文件系统内），确保不会留下损坏的中间文件。

    参数:
        adata: AnnData 对象
        target: 目标 .h5ad 路径
        tmpdir: 临时目录（cfg 传入时优先使用 cfg.h5ad_tempdir）
        compression: 默认 h5py 压缩方式 ('gzip' | 'lzf' | 'zstd')
        cfg: 可选的 Config 对象 — 传入后优先使用 cfg.h5ad_compression
        compression_override: 显式覆盖 — 优先级高于 cfg.h5ad_compression。
            用于 SnapATAC2 兼容写（compression_override=None 写未压缩文件）。
    """
    # Resolution order: compression_override > cfg.h5ad_compression > compression default
    if compression_override is not None:
        compression = compression_override
    elif cfg is not None:
        compression = getattr(cfg, 'h5ad_compression', compression)
    # Respect cfg.h5ad_tempdir (from ATACseq config)
    if cfg is not None:
        tmpdir = getattr(cfg, 'h5ad_tempdir', tmpdir)

    # WSL /mnt mounts require tmp+mv to avoid h5py file locking issues.
    _wsl = (target.startswith("/mnt/")
            and os.environ.get("HDF5_USE_FILE_LOCKING", "").upper() != "FALSE")
    logging.getLogger(__name__).info("Writing %s ...", os.path.basename(target))
    if _wsl:
        os.makedirs(tmpdir, exist_ok=True)
        tmp_path = os.path.join(tmpdir, os.path.basename(target))
        adata.write(tmp_path, compression=compression)
        shutil.move(tmp_path, target)
    else:
        adata.write(target, compression=compression)

    size_mb = os.path.getsize(target) / 1e6
    logger = logging.getLogger(__name__)
    logger.info("Saved %s (%.1f MB)", os.path.basename(target), size_mb)


def safe_plot(func, *args, **kwargs):
    """
    容错的 scanpy 绘图包装。

    某些 scanpy 绘图函数在某些版本组合下可能因 matplotlib 兼容性崩溃。
    本函数捕获异常并记录警告，避免整个步骤因此中断。
    自动处理已弃用的 save 参数 — 拦截并改用 plt.savefig。

    用法:
        safe_plot(sc.pl.umap, adata, color='stage', show=False, save='_stage.pdf')
    """
    logger = logging.getLogger(__name__)
    save_path = kwargs.pop('save', None)
    if save_path:
        kwargs.setdefault('show', False)
    try:
        result = func(*args, **kwargs)
        if save_path:
            import matplotlib.pyplot as plt
            if not os.path.isabs(save_path):
                save_path = os.path.join(sc.settings.figdir, save_path)
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            plt.close()
        return result
    except Exception as e:
        logger.warning("Plot failed (skipped): %s", e)
        return None


def setup_logger(name: str, log_file: str,
                 level: int = logging.INFO,
                 force: bool = False) -> logging.Logger:
    """
    统一配置日志: 同时输出到 stdout 和文件。

    格式:
        14:30:00 | INFO    | 消息内容

    参数:
        name: logger 名称
        log_file: 日志文件路径
        level: 日志级别
        force: 是否强制重建 handler（清除已有 handler）

    返回:
        配置好的 logger 实例
    """
    import pandas as pd
    warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
    warnings.filterwarnings("ignore", message=".*fragmented.*")

    logger = logging.getLogger(name)
    logger.setLevel(level)

    if force and logger.handlers:
        logger.handlers.clear()

    if logger.handlers:
        return logger

    formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-7s | %(message)s',
        datefmt='%H:%M:%S',
    )

    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setFormatter(formatter)
    logger.addHandler(stdout_h)

    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    file_h = logging.FileHandler(log_file, mode='w')
    file_h.setFormatter(formatter)
    logger.addHandler(file_h)

    return logger


def resolve_config(config_path: Optional[str] = None):
    """
    解析 --config 参数，返回配置模块的 CFG 对象。

    所有步骤脚本统一使用本函数加载配置。
    """
    if config_path is None:
        # 默认寻找父目录的 config.py
        config_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "config.py",
        )

    config_path = os.path.abspath(config_path)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"配置文件不存在: {config_path}")

    import importlib.util
    spec = importlib.util.spec_from_file_location("pipeline_config", config_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_config"] = mod
    # Also register as "config" for backward compatibility with scripts that do
    # "from config import CFG" after calling resolve_config()
    sys.modules["config"] = mod
    spec.loader.exec_module(mod)

    # Auto-detect project_dir from config file location if not explicitly set.
    if not mod.CFG.project_dir:
        mod.CFG.project_dir = os.path.dirname(config_path)

    # ── Resolve n_jobs ───────────────────────────────────────────────
    # 0 means "auto-detect" but joblib.Parallel rejects 0 outright.
    # Resolve here so both standalone step scripts and subprocess steps
    # get a usable value (run_pipeline.py also does this, but steps
    # shouldn't depend on the launcher).
    _nc = getattr(mod.CFG, 'n_jobs', 0)
    if _nc == 0:
        mod.CFG.n_jobs = os.cpu_count() or 1

    mod.CFG.resolve_paths()

    # ── Species sanity check ───────────────────────────────────────────
    # Warn early when the species string won't match any KB entry, rather
    # than silently failing downstream (Fisher scores all zero, etc.).
    _validate_species(mod.CFG)

    return mod.CFG


# Species pipeline-keys known to rna/ortholog.py and rna/utils/marker_scoring.py.
# Keep in sync when adding new species.
_KNOWN_SPECIES_KEYS = frozenset({
    "human", "mouse", "macaque", "cynomolgus", "marmoset",
    "zebrafish", "chicken", "lamprey", "frog",
    "pig", "cow", "sheep", "ferret",
    "squirrel", "opossum", "tree_shrew", "treeshrew",
    "lizard", "anolis", "peromyscus", "deer_mouse",
    "rhabdomys", "striped_mouse",
    # mus_musculus is used by some legacy mouse configs
    "mus_musculus",
})


def _validate_species(cfg) -> None:
    """Log a warning if *cfg.species* is not a recognised pipeline key."""
    species = getattr(cfg, 'species', '')
    if not species:
        return  # empty species = user hasn't set it yet; not necessarily wrong
    if species not in _KNOWN_SPECIES_KEYS:
        logging.getLogger("core").warning(
            "CFG.species=%r is not a recognised species key. "
            "KB marker scoring and phylogenetic weighting may silently fail. "
            "Add this species to _SPECIES_SYNONYMS in rna/utils/marker_scoring.py "
            "if it is a valid species.",
            species,
        )


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

    X_data = adata.X.data if sp.issparse(adata.X) else adata.X
    n_nan = int(np.isnan(X_data).sum())
    n_inf = int(np.isinf(X_data).sum())
    total = n_nan + n_inf

    if total == 0:
        logger.info("[%s] X matrix clean — no NaN/Inf values", stage_name or "validate")
        return False

    logger.warning(
        "[%s] Found %d NaN and %d Inf values in X matrix — fixing",
        stage_name or "validate", n_nan, n_inf,
    )

    if fix_nan_inf:
        if sp.issparse(adata.X):
            adata.X.data = np.nan_to_num(adata.X.data, nan=0, posinf=0, neginf=0)
        else:
            adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)

    return True


# ── RNA result auto-discovery for ATAC integration ──────────────────

def find_rna_h5ad(cfg=None, dataset_id: str = None, log=None) -> Optional[str]:
    """Auto-discover the RNA annotated h5ad for ATAC Step 09 integration.

    Invoked automatically by run_pipeline.py when --modality atac and
    cfg.rna_h5ad is empty.  Also callable from 09_integrate.py as a
    last-resort fallback.

    Search order (first existing file wins):
      1. cfg.rna_h5ad (explicitly set by user — always honoured first)
      2. projects/rna/{dataset_id}/results/h5ad/05_annotated.h5ad
      3. projects/rna/{dataset_id}/results/h5ad/04_clustered.h5ad
      4. projects/rna/{dataset_id}/results/h5ad/03_integrated.h5ad

    ``dataset_id`` is inferred from ``cfg.project_dir`` when omitted.
    """
    if log is None:
        log = logging.getLogger(__name__)

    # ── 1. Explicit setting always wins ──────────────────────────────
    if cfg is not None and getattr(cfg, 'rna_h5ad', ''):
        return cfg.rna_h5ad

    # ── 2. Derive dataset ID from config ─────────────────────────────
    if dataset_id is None and cfg is not None:
        proj = getattr(cfg, 'project_dir', '')
        if proj:
            dataset_id = os.path.basename(os.path.normpath(proj))
    if not dataset_id:
        return None

    # ── 3. Locate repo root ──────────────────────────────────────────
    repo = repo_root()  # e.g. /mnt/d/Projects/Fuxi or D:\Projects\Fuxi

    # Candidate h5ad files (most complete first)
    candidates = [
        "05_annotated.h5ad",
        "04_clustered.h5ad",
        "03_integrated.h5ad",
    ]

    rna_project_dir = os.path.join(repo, "projects", "rna", dataset_id)
    if not os.path.isdir(rna_project_dir):
        log.debug("find_rna_h5ad: no RNA project dir at %s", rna_project_dir)
        return None

    h5ad_dir = os.path.join(rna_project_dir, "results", "h5ad")
    if not os.path.isdir(h5ad_dir):
        # Fallback: try legacy structure where h5ad is in results/ directly
        h5ad_dir = os.path.join(rna_project_dir, "results")

    for fname in candidates:
        full = os.path.join(h5ad_dir, fname)
        if os.path.isfile(full) and os.path.getsize(full) > 0:
            log.info("Auto-discovered RNA h5ad: %s", full)
            return full

    log.debug("find_rna_h5ad: no RNA h5ad found in %s", h5ad_dir)
    return None


# ── scRNA marker CSV auto-discovery for spatial annotation ──────────

def find_rna_marker_csv(cfg=None, dataset_id: str = None, log=None) -> Optional[str]:
    """Auto-discover the scRNA marker CSV for spatial annotation (Phase 1).

    Invoked from ``spatial/steps/05_annotate.py`` when ``cfg.rna_ref``
    is set.  Returns the path to a ``marker_genes_per_group_cell_type.csv``
    file produced by RNA Step 07.

    Search order (first existing file wins):
      1. ``cfg.rna_ref`` — if path to an existing CSV file, use directly
      2. ``projects/rna/{id}/results/tables/marker_genes_per_group_cell_type.csv``
      3. ``projects/rna/{id}/results/tables/marker_genes_per_group.csv``  (fallback)

    ``dataset_id`` is inferred from ``cfg.rna_ref`` when that is a bare
    ID string, otherwise from ``cfg.project_dir``.
    """
    if log is None:
        log = logging.getLogger(__name__)

    # ── 1. Explicit path to CSV — use directly ────────────────────────
    if cfg is not None and getattr(cfg, 'rna_ref', ''):
        ref = cfg.rna_ref
        if ref.endswith('.csv') and os.path.isfile(ref):
            return ref
        # If rna_ref is a directory path, use its basename as dataset_id
        if os.path.isdir(ref):
            dataset_id = os.path.basename(os.path.normpath(ref))

    # ── 2. Resolve dataset_id from rna_ref string or project_dir ─────
    if dataset_id is None and cfg is not None:
        ref = getattr(cfg, 'rna_ref', '')
        if ref and not os.path.sep in ref and '/' not in ref:
            dataset_id = ref   # bare ID like "GSE235585"
    if dataset_id is None and cfg is not None:
        proj = getattr(cfg, 'project_dir', '')
        if proj:
            dataset_id = os.path.basename(os.path.normpath(proj))
    if not dataset_id:
        log.debug("find_rna_marker_csv: cannot determine dataset_id")
        return None

    # ── 3. Locate repo root ───────────────────────────────────────────
    repo = repo_root()

    # Candidate tables (cell-type markers preferred over leiden markers)
    candidates = [
        "marker_genes_per_group_cell_type.csv",
        "marker_genes_per_group.csv",
    ]

    rna_project_dir = os.path.join(repo, "projects", "rna", dataset_id)
    if not os.path.isdir(rna_project_dir):
        log.debug("find_rna_marker_csv: no RNA project dir at %s", rna_project_dir)
        return None

    tables_dir = os.path.join(rna_project_dir, "results", "tables")
    for fname in candidates:
        full = os.path.join(tables_dir, fname)
        if os.path.isfile(full) and os.path.getsize(full) > 0:
            log.info("Auto-discovered RNA marker CSV: %s", full)
            return full

    log.debug("find_rna_marker_csv: no marker CSV found in %s", tables_dir)
    return None


def load_scRNA_markers(
    csv_path: str,
    top_n: int = 10,
    pval_threshold: float = 0.05,
    logfc_min: float = 0.0,
    log=None,
) -> Dict[str, List[str]]:
    """Extract per-cell-type top-N marker genes from scRNA DE results.

    Reads a ``marker_genes_per_group_cell_type.csv`` file (produced by
    RNA Step 07) and converts it to ``CFG.marker_dict`` format.

    Parameters
    ----------
    csv_path : str
        Path to the RNA marker CSV (columns: ``group``, ``names``,
        ``scores``, ``logfoldchanges``, ``pvals``, ``pvals_adj``,
        ``pct_nz_group``, ``pct_nz_reference``).
    top_n : int
        Number of top marker genes to extract per cell type.
    pval_threshold : float
        Only consider genes with ``pvals_adj`` < this value.
    logfc_min : float
        Minimum ``logfoldchanges`` to include (0.0 = positive LFC only).

    Returns
    -------
    Dict[str, List[str]]
        ``{cell_type: [gene1, gene2, ..., geneN]}`` suitable for
        ``CFG.marker_dict``.
    """
    import pandas as pd

    if log is None:
        log = logging.getLogger(__name__)

    df = pd.read_csv(csv_path)
    log.info(
        "load_scRNA_markers: read %d rows from %s",
        len(df), os.path.basename(csv_path),
    )

    # Filter: significant + positively enriched
    mask = (
        (df['pvals_adj'] < pval_threshold)
        & (df['logfoldchanges'] > logfc_min)
    )
    df_filt = df[mask]
    log.info(
        "  after filtering (pval<%.0e, logfc>%.1f): %d rows retained",
        pval_threshold, logfc_min, len(df_filt),
    )

    # Group by cell type, sort by score descending, take top-N
    marker_dict: Dict[str, List[str]] = {}
    for group, group_df in df_filt.groupby('group', observed=True):
        top_genes = (
            group_df
            .sort_values('scores', ascending=False)
            .head(top_n)['names']
            .tolist()
        )
        if top_genes:
            marker_dict[group] = top_genes
            log.debug("  %s: %d markers → %s", group, len(top_genes),
                      top_genes[:5])

    log.info(
        "load_scRNA_markers: extracted markers for %d cell types "
        "(top_n=%d, pval<%.0e, logfc>%.1f)",
        len(marker_dict), top_n, pval_threshold, logfc_min,
    )
    return marker_dict


@dataclass
class PerformanceReport:
    """Performance metrics for a pipeline step."""
    step: str = ""
    wall_sec: float = 0.0
    cpu_sec: float = 0.0
    peak_rss_mb: float = 0.0
    avg_cpu_pct: float = 0.0
    gpu_mem_mb: float = -1.0


@contextmanager
def monitor_performance(step_name: str = "", log=None):
    """Time CPU, memory, and GPU usage for a code block."""
    report = PerformanceReport(step=step_name)
    t0 = _time.time(); tcpu0 = _time.process_time()
    import psutil
    proc = psutil.Process()
    cpu_samples: list[float] = []
    peak_rss = 0
    stop = threading.Event()

    def _sample():
        nonlocal peak_rss
        while not stop.is_set():
            try:
                m = proc.memory_info().rss
                if m > peak_rss: peak_rss = m
                cpu_samples.append(proc.cpu_percent())
            except Exception:
                pass
            stop.wait(1.0)

    sampler = threading.Thread(target=_sample, daemon=True)
    sampler.start()
    try:
        yield report
    finally:
        stop.set(); sampler.join(timeout=5)
        dt = _time.time() - t0
        report.wall_sec = round(dt, 1)
        report.cpu_sec = round(_time.process_time() - tcpu0, 1)
        report.peak_rss_mb = round(peak_rss / 1e6, 1)
        report.avg_cpu_pct = round(sum(cpu_samples) / max(len(cpu_samples), 1), 1)
        try:
            out = _sp.check_output(["nvidia-smi","--query-gpu=memory.used",
                "--format=csv,noheader,nounits"]).decode().strip().split("\n")
            report.gpu_mem_mb = sum(float(m) for m in out)
        except Exception:
            report.gpu_mem_mb = -1.0
        if log:
            log.info("[perf] wall=%.1fs cpu=%.1fs mem=%.1fMB cpu%%=%.1f%% gpu=%.0fMB",
                     report.wall_sec, report.cpu_sec, report.peak_rss_mb,
                     report.avg_cpu_pct, report.gpu_mem_mb)
