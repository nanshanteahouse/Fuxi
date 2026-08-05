"""I/O utilities — safe h5ad write and plot wrappers."""

import logging
import os
import shutil
from typing import Optional

try:
    import hdf5plugin  # 注册第三方 HDF5 滤镜（zstd 等）；读 zstd 文件必须

    # 让不 import core 的子进程/独立脚本也能读 zstd 文件（H5PL 动态插件路径）
    os.environ.setdefault(
        "HDF5_PLUGIN_PATH",
        os.path.join(os.path.dirname(hdf5plugin.__file__), "plugins"),
    )
except ImportError:
    hdf5plugin = None
    logging.getLogger(__name__).warning(
        "hdf5plugin missing: zstd-compressed h5ad files (e.g. obsm/X_integrated)",
        " cannot be read. Run: pip install hdf5plugin",
    )

import anndata

from core.utils._io_incremental import write_h5ad_incremental
from core.utils._path import is_wsl


def _resolve_compression_kwargs(compression: str, compression_opts: Optional[int]):
    """压缩名 → h5py 写参数；zstd 需 hdf5plugin 滤镜对象（h5py 不识别字符串）。"""
    if compression == "zstd" and hdf5plugin is not None:
        return dict(hdf5plugin.Zstd(clevel=1))
    if compression == "zstd":
        return {"compression": "gzip"}
    kwargs = {"compression": compression}
    if compression_opts is not None:
        kwargs["compression_opts"] = compression_opts
    return kwargs


def safe_write(
    adata,
    target: str,
    tmpdir: str = "/tmp/Fuxi",
    compression: str = "gzip",
    cfg=None,
    compression_override: Optional[str] = None,
    step_alias: str | None = None,
    compression_opts: Optional[int] = None,
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
            append 的 key 一律覆盖。目标不存在或 MuData 时回退全量路径。
        compression_opts: h5py 压缩级别（如 gzip level 1）— 非 None 时透传给
            adata.write(compression_opts=...)。仅当 compression 为 gzip 系列时生效。
            None 时回退 cfg.h5ad_compression_opts（与 compression 同源解析）。
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
    # Compression level: explicit arg > cfg.h5ad_compression_opts (gzip only)
    if compression_opts is None and cfg is not None:
        compression_opts = getattr(cfg, "h5ad_compression_opts", None)
    if compression_opts is not None and not compression.startswith("gzip"):
        compression_opts = None  # opts only valid for gzip family
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
            write_kwargs = _resolve_compression_kwargs(compression, compression_opts)
            adata.write(tmp_path, **write_kwargs)
            os.replace(tmp_path, target)  # atomic same-fs rename
        except Exception:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            raise
    else:
        write_kwargs = _resolve_compression_kwargs(compression, compression_opts)
        adata.write(target, **write_kwargs)

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
                import gc

                import scanpy as sc

                _verify = sc.read(target, backed="r")
                logger.info("Integrity check: %s verified OK", os.path.basename(target))
                # 关键：backed 读持有的 h5py 句柄因循环引用不会立即回收，
                # 若后续 r+ 打开同一文件（如 stream_write_raw 写 /raw）会冲突：
                # gzip 时偶发 OSError，zstd 滤镜下触发 C 层 double free (SIGABRT)。
                # 显式关闭 + gc 确保文件释放。
                try:
                    _verify.file.close()
                except Exception:
                    pass
                del _verify
                gc.collect()
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


# ── 流式写 .raw ──────────────────────────────────────────────────────
def stream_write_raw(
    target: str,
    source: str,
    *,
    target_sum: float = 1e4,
    compression: str = "gzip",
    compression_opts: Optional[int] = None,
    chunk_size: int = 200_000,
    var_source: str | None = None,
    var_df=None,
    transform: bool = True,
    logger=None,
) -> int:
    """Stream-write the .raw (full-gene) group into an existing h5ad file.

    从 source (02_qc.h5ad) 分块读取 counts，逐块 normalize_total + log1p，
    直写 target 的 /raw 组。nnz 在变换前后严格不变（行缩放零保持零、
    log1p(0)=0）→ 一次性 create_dataset，零 resize，避免 HDF5 反复扩展。

    参数:
        target: 已存在的 03 输出 h5ad（r+ 追加 /raw 组）
        source: 02_qc.h5ad（counts 源）
        target_sum: normalize_total 的 target_sum
        chunk_size: 每次分块读取的细胞数
        var_source: 若提供，则从该文件拷贝 /var 到 /raw/var（否则拷贝 source 的 /var）

    返回:
        写入的 nnz 总数。
    """
    import h5py
    import numpy as np
    import scipy.sparse as sp

    log = logger or logging.getLogger(__name__)
    log.info("Stream-writing .raw (full genes) from %s → %s ...", source, target)

    # 压缩解析与 safe_write 一致：zstd 需 hdf5plugin filter 对象（h5py 不识别字符串）
    comp_kwargs = _resolve_compression_kwargs(compression, compression_opts)

    src_h5 = h5py.File(source, "r")
    dst_h5 = h5py.File(target, "r+")
    try:
        src_x = src_h5["X"]
        nnz = src_x["data"].shape[0]
        shape = tuple(src_x.attrs["shape"])
        n_obs = shape[0]
        n_genes = shape[1]
        var_src = var_source or source

        # ── /raw 组 ──
        if "raw" in dst_h5:
            del dst_h5["raw"]
        rg = dst_h5.create_group("raw")

        # raw/var: var_df 优先（anndata 编码，保留 03 算出的 highly_variable 列）
        if var_df is not None:
            import tempfile

            _vtmp = os.path.join(tempfile.gettempdir(), f"fuxi_rawvar_{os.getpid()}.h5ad")
            try:
                ad_tmp = anndata.AnnData(var=var_df)
                # NOTE(2026-08-01): var 仅 ~1-2MB，永远用 gzip 写——
                # 实测 anndata 写 zstd 临时文件在 GPU 上下文 (cupy) 共存时非确定性触发
                # C 层堆损坏 (double free / free(): invalid size，~2/3 概率)，gzip 稳定。
                # 大数组 data/indices/indptr 仍走主压缩 (zstd)，性能不受影响。
                var_kwargs = {"compression": "gzip"}
                if isinstance(comp_kwargs.get("compression_opts"), int):
                    var_kwargs["compression_opts"] = 1
                ad_tmp.write(_vtmp, **var_kwargs)
                with h5py.File(_vtmp, "r") as vh:
                    vh.copy("var", rg, name="var")
            finally:
                if os.path.exists(_vtmp):
                    os.unlink(_vtmp)
        elif var_src == source:
            src_h5.copy("var", rg, name="var")
        else:
            with h5py.File(var_src, "r") as vh:
                vh.copy("var", rg, name="var")
        # raw/X: csr_matrix 编码
        rg.create_group("X")
        rg["X"].attrs["encoding-type"] = "csr_matrix"
        rg["X"].attrs["encoding-version"] = "0.1.0"
        rg["X"].attrs["shape"] = np.array(shape, dtype=np.int64)
        d = rg["X"].create_dataset("data", shape=(nnz,), dtype="f4", **comp_kwargs)
        rg["X"].create_dataset("indices", data=src_x["indices"][...], **comp_kwargs)
        rg["X"].create_dataset("indptr", data=src_x["indptr"][...], **comp_kwargs)

        # ── 分块变换并写入 data ──
        indptr = src_x["indptr"][...]
        pos = 0
        for i in range(0, n_obs, chunk_size):
            j = min(i + chunk_size, n_obs)
            lo, hi = int(indptr[i]), int(indptr[j])
            if hi <= lo:
                continue
            cdata = src_x["data"][lo:hi]
            cidx = src_x["indices"][lo:hi]
            cptr = indptr[i : j + 1] - indptr[i]
            c = sp.csr_matrix((cdata, cidx, cptr), shape=(j - i, n_genes))
            if transform:
                # normalize_total (与 scanpy bit 一致):
                #   counts_per_cell = 每行 float64 求和 → float32 → /target_sum
                #   X = X / counts_per_cell (除法, 非乘倒数)
                rs = np.zeros(c.shape[0], dtype=np.float64)
                for k in range(c.shape[0]):
                    rs[k] = c.data[c.indptr[k] : c.indptr[k + 1]].sum(dtype=np.float64)
                cpc = (rs.astype(np.float32) / np.float32(target_sum)).astype(np.float32)
                c = c.copy()
                for k in range(c.shape[0]):
                    lo2, hi2 = c.indptr[k], c.indptr[k + 1]
                    if hi2 > lo2:
                        c.data[lo2:hi2] = c.data[lo2:hi2] / cpc[k]
                # log1p 原地
                c.data = np.log1p(c.data).astype(np.float32)
            d[pos : pos + c.nnz] = c.data
            pos += c.nnz
        log.info(".raw streamed: %d nnz (%d × %d)", pos, n_obs, n_genes)
        return pos
    finally:
        src_h5.close()
        dst_h5.close()
