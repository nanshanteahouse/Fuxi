"""Logging utilities — unified log configuration."""

import logging
import os
import sys
import warnings
from datetime import datetime

# Module-level guard: root logger config happens once per process.
_root_configured = False


def setup_logger(
    name: str, log_file: str, level: int = logging.INFO, force: bool = False
) -> logging.Logger:
    """
    统一配置日志: 同时输出到 stdout 和文件。

    格式（终端）:
        14:30:00 | INFO    | step_name | 消息内容
    格式（文件）:
        2026-07-29 14:30:00 | INFO    | step_name | 消息内容

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

    global _root_configured

    logger = logging.getLogger(name)
    logger.setLevel(level)

    if force:
        logger.handlers.clear()
        # 强制模式也重建 root logger,保证格式一致
        root = logging.getLogger()
        root.handlers.clear()
        _root_configured = False

    if logger.handlers:
        return logger

    # ── 两套 Formatter ────────────────────────────────────────────────
    fmt_console = logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        datefmt="%H:%M:%S",
    )
    fmt_file = logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # ── Step logger 的 handler ─────────────────────────────────────────
    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setFormatter(fmt_console)
    stdout_h.setLevel(logging.INFO)
    logger.addHandler(stdout_h)

    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # 文件追加模式,启动时写一条分隔横幅
    run_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as f:
        f.write(f"\n===== run [{name}] {run_ts} =====\n")

    file_h = logging.FileHandler(log_file, mode="a")
    file_h.setFormatter(fmt_file)
    file_h.setLevel(logging.DEBUG)
    logger.addHandler(file_h)

    # P0: 不再向 root 冒泡,消除库 basicConfig 造成的格式混杂 + 重复行
    logger.propagate = False

    # ── Root logger 接管（整个进程只做一次）────────────────────────────
    if not _root_configured:
        root = logging.getLogger()
        root.setLevel(logging.INFO)

        # Suppress verbose DEBUG from numba JIT (floods 1000s of lines in grid search)
        logging.getLogger("numba").setLevel(logging.WARNING)
        # 清掉库 import 时 basicConfig 装的 handler
        root.handlers.clear()

        root_stdout = logging.StreamHandler(sys.stdout)
        root_stdout.setFormatter(fmt_console)
        root_stdout.setLevel(logging.INFO)
        root.addHandler(root_stdout)

        root_file = logging.FileHandler(log_file, mode="a")
        root_file.setFormatter(fmt_file)
        root_file.setLevel(logging.DEBUG)
        root.addHandler(root_file)

        _root_configured = True

    return logger
