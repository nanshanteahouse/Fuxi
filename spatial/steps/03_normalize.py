#!/usr/bin/env python3
"""
Step 03: Normalization + HVG selection + spatial neighbor graph
==================================================================
Phase 1: Store raw counts in adata.raw (BEFORE any transformation)
Phase 2: Library-size normalize + log1p
Phase 3: HVG selection with fallback chain + forced gene inclusion
Phase 4: Build spatial neighbor graph via sq.gr.spatial_neighbors()
Phase 5: PCA on HVGs
Phase 6: Optional multi-slide Harmony batch integration (method='harmony')
Phase 7: Checkpoint write (03_processed.h5ad) — BEFORE any plotting
Phase 8: Harmony comparison plot (after checkpoint)

Input:  02_image.h5ad (or 01_qc.h5ad if Step 02 was skipped)
Output: 03_processed.h5ad

Integration is opt-in and off by default for spatial (read via
``getattr(cfg.spatial, 'integration_method', ...)`` — spatial-scoped enum,
default 'none'). It runs only when method=='harmony' AND the obs batch_key
(default 'sample', the column added by spatial 00_load for merged slides)
has >= 2 unique values (multi-slide). Single-slide data skips.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import scanpy as sc
import scipy.sparse
import squidpy as sq

from core.utils import (
    gpu_harmony,
    resolve_config,
    resolve_device,
    safe_write,
    save_figure,
    setup_logger,
    timed_substep,
)


def _scan_collinearity_warnings(harmony_result) -> list[str]:
    """Extract perfectly-collinear warnings from a harmony run result.

    harmonypy itself does not expose a warnings channel in this version, so
    the scan is attribute-safe: a ``report`` object carrying a ``warnings``
    list is scanned when present; otherwise the guard no-ops.

    Ref: notes/engineering/2026-07-15_cross_batch_critical_fixes.md (T3).
    """
    report = getattr(harmony_result, "report", None)
    warnings_list = getattr(report, "warnings", None) or []
    return [w for w in warnings_list if "perfectly collinear" in str(w).lower()]


def _effective_n_pcs(n_pcs_use: int, n_embedding_dims: int) -> int:
    """Cap a requested n_pcs at the actual embedding width (min() guard).

    Integration backends may emit fewer dims than ``cfg.pca.n_pcs_use`` (e.g. a
    reduced-dim corrected embedding). Downstream PCA/neighbor consumers must
    never slice beyond ``obsm[...].shape[1]``.

    Ref: notes/engineering/2026-07-28_n_pcs_use_scvi_dimension_mismatch.md.
    """
    return min(n_pcs_use, n_embedding_dims)


def _log_integrated_shape(adata, cfg, log, backend: str) -> None:
    """Log the corrected embedding shape and apply the n_pcs min() guard."""
    n_integrated = adata.obsm["X_integrated"].shape[1]
    n_pcs_use = getattr(cfg.pca, "n_pcs_use", 30)
    if n_integrated < n_pcs_use:
        log.warning(
            "X_integrated has %d dims < cfg.pca.n_pcs_use=%d — downstream PCA/neighbor ",
            "n_pcs capped to min(%d, %d)=%d",
            n_integrated,
            n_pcs_use,
            n_pcs_use,
            n_integrated,
            _effective_n_pcs(n_pcs_use, n_integrated),
        )
    log.info("Harmony complete (%s), output shape: %s", backend, adata.obsm["X_integrated"].shape)


def _integrate_harmony(adata, cfg, log) -> None:
    """GPU-first multi-slide Harmony batch correction with collinearity guard.

    Runs only when ``cfg.integration.method == 'harmony'`` AND the obs
    ``batch_key`` has >= 2 unique values (multi-slide). Single-slide data
    skips with a log.info and X_pca remains the primary representation.

    Collinearity guard (CPU path): the harmonypy run result's report is
    scanned for "perfectly collinear" warnings; when found the corrected
    embedding is NOT applied and ``uns['harmony_skipped']`` is recorded.
    """
    batch_key = getattr(cfg.integration, "batch_key", "sample")
    if batch_key not in adata.obs:
        log.info("batch_key '%s' not in obs — skipping Harmony integration", batch_key)
        return
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.info(
            "Single slide (obs['%s'] has %d unique value) — skipping Harmony",
            batch_key,
            n_batches,
        )
        return
    if "X_pca" not in adata.obsm:
        log.warning("obsm['X_pca'] missing — skipping Harmony integration")
        return

    device = getattr(cfg.execution, "device", "auto")
    random_state = getattr(cfg.execution, "random_seed", 0)
    max_iter = getattr(cfg.integration, "max_iter", 20)
    # Harmony runs on the PCA embedding (obsm['X_pca']) — never on raw counts
    # (use_raw guard: adata.raw / adata.X stay untouched).

    # ── GPU path (rapids-singlecell) — no run report exposed ──
    if resolve_device(device, log):
        try:
            with timed_substep("Harmony (GPU)", log=log):
                gpu_harmony(
                    adata,
                    key=batch_key,
                    output_key="X_integrated",
                    log=log,
                    device=device,
                    random_state=random_state,
                    max_iter_harmony=max_iter,
                )
            _log_integrated_shape(adata, cfg, log, backend="GPU")
            return
        except Exception as e:
            log.warning("GPU Harmony failed (%s) — falling back to CPU", e)

    # ── CPU path (harmonypy) — run report available → collinearity guard ──
    try:
        import harmonypy as hm
    except ImportError:
        log.warning(
            "harmonypy not installed — skipping Harmony integration (install fuxi[rna])",
        )
        return
    try:
        with timed_substep("Harmony (CPU)", log=log):
            ho = hm.run_harmony(
                adata.obsm["X_pca"],
                adata.obs,
                vars_use=[batch_key],
                random_state=random_state,
                max_iter_harmony=max_iter,
            )
        collinear = _scan_collinearity_warnings(ho)
        if collinear:
            log.error(
                "[collinearity-guard] Harmony ABORTED — batch_key '%s' perfectly ",
                "collinear with biology",
                batch_key,
            )
            for w in collinear[:3]:
                log.error("  %s", w)
            adata.uns["harmony_skipped"] = {
                "reason": "collinearity",
                "warnings": collinear,
            }
            return
        adata.obsm["X_integrated"] = ho.Z_corr
        _log_integrated_shape(adata, cfg, log, backend="CPU")
    except Exception as e:
        log.warning("Harmony correction failed (%s) — continuing with raw PCA", e)
        adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()


def _plot_harmony_comparison(adata, cfg, log) -> None:
    """Draw PCA-before vs Harmony-after comparison (checkpoint already written)."""
    import matplotlib.pyplot as plt

    batch_key = getattr(cfg.integration, "batch_key", "sample")
    fig_dir = os.path.join(cfg.figure_dir, "03_normalize")
    os.makedirs(fig_dir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.embedding(
        adata,
        basis="X_pca",
        color=batch_key,
        ax=axes[0],
        show=False,
        title="PCA (before Harmony)",
    )
    sc.pl.embedding(
        adata,
        basis="X_integrated",
        color=batch_key,
        ax=axes[1],
        show=False,
        title="Harmony-corrected",
    )
    fig.tight_layout()
    save_figure(
        fig,
        os.path.join(fig_dir, "harmony_comparison"),
        cfg=cfg,
        dpi=cfg.plot.figure_dpi,
    )
    plt.close(fig)
    log.info("  Harmony comparison plot saved")


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("03_normalize", os.path.join(cfg.log_dir, "03_normalize.log"))
    log.info("Step 03: Normalization + HVG + spatial graph + PCA")

    output_path = os.path.join(cfg.h5ad_dir, "03_processed.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    # ── Determine input ─────────────────────────────────────────────────
    image_path = os.path.join(cfg.h5ad_dir, "02_image.h5ad")
    qc_path = os.path.join(cfg.h5ad_dir, "01_qc.h5ad")
    if os.path.exists(image_path):
        input_path = image_path
    elif os.path.exists(qc_path):
        input_path = qc_path
    else:
        log.error("Neither %s nor %s found. Run Steps 01–02 first.", image_path, qc_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    # ═══ Phase 1: Store raw counts (before any transformation) ═══
    adata.raw = adata.copy()
    log.info("  Raw counts stored in adata.raw")

    # ═══ Phase 2: Normalization ═══
    log.info("Normalizing to target sum=%.0f...", cfg.normalization.normalize_target_sum)
    sc.pp.normalize_total(adata, target_sum=cfg.normalization.normalize_target_sum)
    log.info("  Normalization complete")
    sc.pp.log1p(adata)
    log.info("  log1p transformation applied")

    # ═══ Phase 3: Highly variable genes ═══
    flavors = []
    for f in [cfg.hvg.flavor, "seurat_v3", "cell_ranger", "seurat"]:
        if f not in flavors:
            flavors.append(f)

    hvg_selected = False
    for flavor in flavors:
        try:
            log.info("Selecting %d HVGs (flavor=%s)...", cfg.hvg.n_top_genes, flavor)
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=cfg.hvg.n_top_genes,
                flavor=flavor,
                batch_key=cfg.hvg.batch_key if cfg.has_sample_mapping() else None,
            )
            n_hvg = adata.var["highly_variable"].sum()
            log.info("  Selected %d HVGs (flavor=%s)", n_hvg, flavor)
            hvg_selected = True
            break
        except Exception as e:
            log.warning("  HVG flavor '%s' failed: %s", flavor, e)

    if not hvg_selected:
        log.info("Falling back to manual variance-based HVG selection...")
        x = adata.X
        if x is None:
            raise ValueError("adata.X is None — cannot compute manual variance HVG")
        if scipy.sparse.issparse(x):
            mean = np.array(x.mean(axis=0)).flatten()
            var = np.array(x.multiply(x).mean(axis=0)).flatten() - mean**2
        else:
            var = x.var(axis=0)
        top_idx = np.argsort(var)[::-1][: cfg.hvg.n_top_genes]
        adata.var["highly_variable"] = False
        adata.var.iloc[top_idx, adata.var.columns.get_loc("highly_variable")] = True
        n_hvg = len(top_idx)
        log.info("  Selected %d HVGs (manual variance fallback)", n_hvg)
        hvg_selected = True

    # Force-include genes from CFG.hvg.forced_genes and CFG.marker.marker_dict
    forced_genes = set(cfg.hvg.forced_genes)
    if cfg.marker.marker_dict:
        for genes in cfg.marker.marker_dict.values():
            forced_genes.update(genes)
    genes_in_data = set(adata.var_names)
    force_set = forced_genes & genes_in_data
    if force_set:
        adata.var.loc[list(force_set), "highly_variable"] = True
        log.info("  Force-included %d marker/forced genes as HVG", len(force_set))

    n_hvg_final = adata.var["highly_variable"].sum()
    log.info("  Total HVGs after force-include: %d", n_hvg_final)

    if n_hvg_final == 0:
        log.error("No HVGs selected — check your data quality")
        sys.exit(1)

    # ═══ Phase 4: Spatial neighbor graph ═══
    log.info("Building spatial neighbor graph...")
    sq.gr.spatial_neighbors(
        adata,
        n_neighs=cfg.spatial.neighbors_n,
        radius=None if cfg.spatial.neighbors_radius == 0 else cfg.spatial.neighbors_radius,
        coord_type="generic",
    )
    log.info(
        "  Spatial graph: n_neighs=%d, radius=%.1f",
        cfg.spatial.neighbors_n,
        cfg.spatial.neighbors_radius,
    )

    if "spatial_connectivities" not in adata.obsp:
        log.error("Spatial neighbor graph NOT created — spatial_connectivities missing")
        sys.exit(1)

    # ═══ Phase 5: PCA ═══
    log.info("Computing PCA (n_comps=%d)...", cfg.pca.n_pcs_use)
    sc.pp.pca(
        adata,
        n_comps=cfg.pca.n_pcs_use,
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=cfg.execution.random_seed,
    )
    log.info("  PCA complete: %d components stored", cfg.pca.n_pcs_use)

    # ═══ Phase 6: Multi-slide batch integration (optional, Harmony) ═══
    # Spatial-scoped read first (schema default 'none'); falls back to the
    # shared RNA-style enum for backward compatibility — a missing spatial
    # method means no integration, keeping existing behaviour unchanged.
    integration_method = getattr(
        cfg.spatial, "integration_method", getattr(cfg.integration, "method", "none")
    )
    if integration_method == "harmony":
        _integrate_harmony(adata, cfg, log)
    elif integration_method != "none":
        log.info(
            "Integration method '%s' not supported by spatial Step 03 — skipping",
            integration_method,
        )

    # ═══ Phase 7: Checkpoint-before-plot (hard convention) ═══
    # safe_write MUST precede any comparison plotting so that all computed
    # results survive a plot-stage crash (OOM / rendering failure).
    # Ref: notes/engineering/2026-07-30_checkpoint_before_plot.md
    safe_write(adata, output_path, cfg=cfg)

    # ═══ Phase 8: Harmony comparison plot (after checkpoint) ═══
    if integration_method == "harmony" and "X_integrated" in adata.obsm:
        _plot_harmony_comparison(adata, cfg, log)

    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
