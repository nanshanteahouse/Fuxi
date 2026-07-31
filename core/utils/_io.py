"""I/O utilities — safe h5ad write and plot wrappers."""

import logging
import os
import shutil
from typing import Optional

import anndata

from core.utils._io_incremental import write_h5ad_incremental
from core.utils._path import is_wsl


def safe_write(
    adata,
    target: str,
    tmpdir: str = "/tmp/Fuxi",
    compression: str = "gzip",
    cfg=None,
    compression_override: Optional[str] = None,
    step_alias: str | None = None,
    *,
    delta_only: bool = False,
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
        delta_only: 若为 True 且目标文件已存在（且 adata 是 AnnData），则改走
            write_h5ad_incremental in-place 追加 obs/obsm/obsp/uns——不重写 X，
            追加的 key 一律覆盖。目标不存在或 MuData 时回退全量路径。
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

    # Incremental append path (Item 1.5): append obs/obsm/obsp/uns in place
    # when the target exists. MuData (.h5mu) is not supported by the engine
    # → falls through to the full-write path below.
    if delta_only and os.path.exists(target) and not hasattr(adata, "mod"):
        write_h5ad_incremental(
            target,
            obs=adata.obs,
            obsm=adata.obsm,
            obsp=adata.obsp,
            uns=adata.uns,
            logger=logging.getLogger(__name__),
        )
        return
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


# FICLONE ioctl request (linux/fs.h): clone the src fd's extents into the dst fd.
FICLONE = 0x40049409


def copy_h5ad(src: str, dst: str) -> None:
    """Copy an h5ad file, opportunistically using a reflink (CoW clone).

    On filesystems that support FICLONE (btrfs, xfs, …) the copy is
    near-instant and shares the underlying extents: subsequent in-place
    writes to either side trigger copy-on-write, so the ``.bak`` backup
    never changes when the live file is mutated later. On filesystems
    without reflink support (ext4, WSL DrvFs) — or when the ioctl fails
    for any other reason — this silently falls back to
    :func:`shutil.copy2`, preserving copy2's metadata/permission
    semantics. Only a copy2 failure propagates.

    Parameters
    ----------
    src, dst : str
        Source and destination paths. ``dst`` is created/overwritten.
    """
    try:
        import fcntl  # non-POSIX (e.g. Windows) → ImportError → copy2
    except ImportError:
        pass
    else:
        try:
            with open(src, "rb") as src_file, open(dst, "wb") as dst_file:
                fcntl.ioctl(dst_file, FICLONE, src_file.fileno())
            shutil.copystat(src, dst)  # copy2-equivalent metadata
            return
        except Exception:
            pass  # EINVAL/ENOTTY/EXDEV/EPERM → reflink unsupported → copy2
    shutil.copy2(src, dst)


def write_obs_columns_lightweight(
    h5ad_path: str,
    obs_df,  # pandas DataFrame with only the new columns
    logger=None,
):
    """Append obs columns to an existing h5ad file (NO X rewrite).

    Delegates to :func:`write_h5ad_incremental` (anndata.io.write_elem based)
    — categorical codes get a minimal dtype by n_categories (fixes int8
    wraparound at 127 categories) and string columns keep missing values as
    missing (fixes NaN→"nan" corruption).

    Copy+append mode: the h5ad is derived from an intact source file, so a
    failed append deletes the corrupt copy — re-running the step regenerates
    it (commit d80836b's intentional design, kept).

    Parameters
    ----------
    h5ad_path : str
        Path to the existing h5ad file.
    obs_df : pd.DataFrame
        DataFrame containing ONLY the new columns to write (not all obs).
        Must have the same index as the h5ad's /obs index.
    logger : logging.Logger or None
    """
    _log = logger.info if logger else (lambda msg, *a: None)
    _log("Lightweight obs write: %d column(s) → %s", len(obs_df.columns), h5ad_path)
    try:
        write_h5ad_incremental(h5ad_path, obs=obs_df, logger=logger)
    except Exception:
        if logger:
            logger.error(
                "Lightweight write integrity FAILED — removing corrupt file: %s",
                h5ad_path,
            )
        if os.path.exists(h5ad_path):
            os.remove(h5ad_path)
        raise


def write_obs_columns_inplace(
    h5ad_path: str,
    obs_df,  # pandas DataFrame with the columns to write back
    logger=None,
):
    """Write obs columns back into a SHARED checkpoint, with .bak backup/restore.

    Unlike copy+append mode (delete-on-failure), this mutates a checkpoint the
    rest of the pipeline depends on — deleting it would force expensive
    upstream re-runs (e.g. losing 05_annotated means redoing step 05's AI
    annotation calls). Strategy:

      1. ``copy_h5ad`` to ``<file>.bak`` (same directory; reflink when possible,
         CoW guarantees the backup stays pristine across the in-place write);
      2. append via :func:`write_h5ad_incremental` (which verifies);
      3. success → ``unlink`` the backup; failure → ``os.replace`` the backup
         back over the file (atomic same-dir restore) and re-raise.

    Parameters
    ----------
    h5ad_path : str
        Path to the shared checkpoint h5ad file.
    obs_df : pd.DataFrame
        DataFrame with the columns to write back.
    logger : logging.Logger or None
    """
    _log = logger.info if logger else (lambda msg, *a: None)
    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"Cannot write back to non-existent h5ad: {h5ad_path}")
    bak_path = h5ad_path + ".bak"
    _log("In-place obs writeback: backup → %s", os.path.basename(bak_path))
    copy_h5ad(h5ad_path, bak_path)
    try:
        write_h5ad_incremental(h5ad_path, obs=obs_df, logger=logger)
    except Exception:
        try:
            os.replace(bak_path, h5ad_path)  # atomic same-dir restore
        except OSError as restore_err:
            if logger:
                logger.critical(
                    "In-place writeback FAILED and backup restore FAILED: %s",
                    restore_err,
                )
        if logger:
            logger.error(
                "In-place writeback FAILED — restored %s from %s",
                h5ad_path,
                os.path.basename(bak_path),
            )
        raise
    os.unlink(bak_path)
    _log("In-place writeback OK — backup removed")


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
