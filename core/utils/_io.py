"""I/O utilities — safe h5ad write and plot wrappers."""

import logging
import os
from typing import Optional

import anndata

from core.utils._path import is_wsl


def safe_write(
    adata,
    target: str,
    tmpdir: str = "/tmp/Fuxi",
    compression: str = "gzip",
    cfg=None,
    compression_override: Optional[str] = None,
    step_alias: str | None = None,
) -> None:
    """
    安全写入 h5ad 文件，避免 WSL /mnt 挂载的文件锁定问题。

    策略: 先写入 tmpdir，再 mv 到目标路径。
    mv 是原子操作（在同一文件系统内），确保不会留下损坏的中间文件。

    参数:
        adata: AnnData 或 MuData 对象
        tmpdir: 临时目录（cfg 传入时优先使用 cfg.h5ad_tempdir）
        compression: 默认 h5py 压缩方式 ('gzip' | 'lzf' | 'zstd')
        cfg: 可选的 Config 对象 — 传入后优先使用 cfg.h5ad_compression
        compression_override: 显式覆盖 — 优先级高于 cfg.h5ad_compression。
            用于 SnapATAC2 兼容写（compression_override=None 写未压缩文件）。
        step_alias: 步别名（如 "integrated"）— 在 cfg.per_step_h5ad_compression 中查找压缩配置。
            优先级高于 compression_override 和 cfg.h5ad_compression。
    """
    # Resolution order: per_step_h5ad_compression > compression_override > cfg.h5ad_compression > default
    if step_alias is not None and cfg is not None:
        per_step_cfg = getattr(cfg, "per_step_h5ad_compression", {})
        step_compression = per_step_cfg.get(step_alias)
        if step_compression is not None:
            compression = step_compression
    if compression_override is not None:
        compression = compression_override
    elif cfg is not None:
        compression = getattr(cfg, "h5ad_compression", compression)
    # Respect cfg.h5ad_tempdir (from ATACseq config)
    if cfg is not None:
        tmpdir = getattr(cfg, "h5ad_tempdir", tmpdir)
    anndata.settings.allow_write_nullable_strings = True

    # WSL /mnt mounts: prefer writing to a tmp file then atomic rename.
    # If the configured tmpdir is on a different filesystem than target
    # (e.g. tmpfs -> 9p DrvFs), shutil.copy2 corrupts h5 metadata on 4GB+
    # files. Fall back to a hidden tmp file in target's directory so the
    # final rename is a same-fs atomic os.replace.
    _wsl = is_wsl() and target.startswith("/mnt/")
    logging.getLogger(__name__).info("Writing %s ...", os.path.basename(target))
    if _wsl:
        target_dir = os.path.dirname(target) or "."
        try:
            target_dev = os.stat(target_dir).st_dev
        except OSError:
            target_dev = None
        tmp_dev = None
        if os.path.exists(tmpdir):
            try:
                tmp_dev = os.stat(tmpdir).st_dev
            except OSError:
                tmp_dev = None

        # Cross-device -> use hidden tmp file in target dir for atomic rename.
        if tmp_dev is not None and target_dev is not None and tmp_dev != target_dev:
            os.makedirs(target_dir, exist_ok=True)
            tmp_path = os.path.join(target_dir, f".{os.path.basename(target)}.tmp.{os.getpid()}")
        else:
            os.makedirs(tmpdir, exist_ok=True)
            tmp_path = os.path.join(tmpdir, os.path.basename(target))

        try:
            adata.write(tmp_path, compression=compression)
            os.replace(tmp_path, target)  # atomic same-fs rename
        except Exception:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            raise
    else:
        adata.write(target, compression=compression)

    size_mb = os.path.getsize(target) / 1e6
    logger = logging.getLogger(__name__)
    logger.info("Saved %s (%.1f MB)", os.path.basename(target), size_mb)
    # Verify file integrity (gated by cfg.verify_write_integrity)
    if cfg and cfg.verify_write_integrity:
        if hasattr(adata, "mod"):
            # MuData — scanpy's sc.read does not support .h5mu
            try:
                import muon as mu

                _verify = mu.read_h5mu(target)
                logger.info(
                    "Integrity check (MuData): %s verified OK",
                    os.path.basename(target),
                )
            except ImportError:
                logger.info(
                    "Integrity check skipped for MuData %s (muon not available)",
                    os.path.basename(target),
                )
            except Exception as e:
                logger.error(
                    "Integrity check FAILED for %s: %s — file may be corrupted!",
                    os.path.basename(target),
                    e,
                )
        else:
            try:
                import scanpy as sc

                _verify = sc.read(target, backed="r")
                logger.info("Integrity check: %s verified OK", os.path.basename(target))
            except Exception as e:
                logger.error(
                    "Integrity check FAILED for %s: %s — file may be corrupted!",
                    os.path.basename(target),
                    e,
                )


def safe_plot(
    func,
    *args,
    cfg=None,
    dpi=None,
    fmt=None,
    **kwargs,
):
    """
    容错的 scanpy 绘图包装。

    某些 scanpy 绘图函数在某些版本组合下可能因 matplotlib 兼容性崩溃。
    本函数捕获异常并记录警告，避免整个步骤因此中断。
    自动处理已弃用的 save 参数 — 拦截并改用 plt.savefig。

    参数:
        func: 绘图函数 (e.g. sc.pl.umap)
        *args: 传递给 func 的位置参数
        cfg: 可选的 Config 对象 — 使用 cfg.plot.figure_dpi / figure_format / figure_transparent
        dpi: 显式 DPI — 优先级高于 cfg.plot.figure_dpi，高于 150
        fmt: 显式图片格式 — 优先级高于 cfg.plot.figure_format，高于 "pdf"
        **kwargs: 传递给 func 的关键字参数

    用法:
        safe_plot(sc.pl.umap, adata, color='stage', show=False, save='_stage.pdf')
        safe_plot(sc.pl.umap, adata, cfg=cfg, color='stage', show=False, save='_stage')
    """

    import scanpy as sc

    logger = logging.getLogger(__name__)
    save_path = kwargs.pop("save", None)
    if save_path:
        kwargs.setdefault("show", False)
    try:
        result = func(*args, **kwargs)
        if save_path:
            import matplotlib.pyplot as plt

            # Resolve DPI: explicit dpi > cfg.plot.figure_dpi > 150
            final_dpi = dpi or (cfg.plot.figure_dpi if cfg else 150)
            # Resolve format: explicit fmt > cfg.plot.figure_format > "pdf"
            final_fmt = fmt or (cfg.plot.figure_format if cfg else "pdf")
            # Resolve transparency
            transparent = cfg.plot.figure_transparent if cfg else True

            if not os.path.isabs(save_path):
                save_path = os.path.join(sc.settings.figdir, save_path)
            # If save_path has no extension, append final_fmt
            if not os.path.splitext(save_path)[1]:
                save_path = f"{save_path}.{final_fmt}"
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            plt.savefig(save_path, dpi=final_dpi, bbox_inches="tight", transparent=transparent)
            plt.close()
        return result
    except Exception as e:
        logger.warning("Plot failed (skipped): %s", e)
        return None
