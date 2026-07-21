"""Logging utilities — unified log configuration."""

import logging
import os
import sys
import warnings


def setup_logger(
    name: str, log_file: str, level: int = logging.INFO, force: bool = False
) -> logging.Logger:
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
        "%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%H:%M:%S",
    )

    stdout_h = logging.StreamHandler(sys.stdout)
    stdout_h.setFormatter(formatter)
    logger.addHandler(stdout_h)

    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    file_h = logging.FileHandler(log_file, mode="w")
    file_h.setFormatter(formatter)
    logger.addHandler(file_h)

    return logger
