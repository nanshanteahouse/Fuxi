"""I/O utilities — safe h5ad write and plot wrappers."""

import logging
import os
import shutil
from typing import Optional

import anndata


def safe_write(
    adata,
    target: str,
    tmpdir: str = "/tmp/Fuxi",
    compression: str = "gzip",
    cfg=None,
    compression_override: Optional[str] = None,
) -> None:
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
        compression = getattr(cfg, "h5ad_compression", compression)
    # Respect cfg.h5ad_tempdir (from ATACseq config)
    if cfg is not None:
        tmpdir = getattr(cfg, "h5ad_tempdir", tmpdir)
    anndata.settings.allow_write_nullable_strings = True

    # WSL /mnt mounts: use explicit copy+unlink instead of shutil.move
    # to avoid "Invalid cross-device link" on DrvFs (Windows mounts).
    _wsl = target.startswith("/mnt/")
    logging.getLogger(__name__).info("Writing %s ...", os.path.basename(target))
    if _wsl:
        os.makedirs(tmpdir, exist_ok=True)
        tmp_path = os.path.join(tmpdir, os.path.basename(target))
        adata.write(tmp_path, compression=compression)
        shutil.copy2(tmp_path, target)
        os.unlink(tmp_path)
    else:
        adata.write(target, compression=compression)

    size_mb = os.path.getsize(target) / 1e6
    logger = logging.getLogger(__name__)
    logger.info("Saved %s (%.1f MB)", os.path.basename(target), size_mb)

    # Verify file integrity
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
