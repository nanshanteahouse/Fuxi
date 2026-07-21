#!/usr/bin/env python3
"""
Step 03: Normalization + HVG selection + spatial neighbor graph
==================================================================
Phase 1: Library-size normalize + log1p
Phase 2: Store raw counts in adata.raw
Phase 3: HVG selection with fallback chain + forced gene inclusion
Phase 4: Build spatial neighbor graph via sq.gr.spatial_neighbors()
Phase 5: PCA on HVGs

Input:  02_image.h5ad (or 01_qc.h5ad if Step 02 was skipped)
Output: 03_processed.h5ad
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

from core.utils import resolve_config, safe_write, setup_logger


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

    # ═══ Phase 1: Normalization ═══
    log.info("Normalizing to target sum=%.0f...", cfg.normalization.normalize_target_sum)
    sc.pp.normalize_total(adata, target_sum=cfg.normalization.normalize_target_sum)
    log.info("  Normalization complete")
    sc.pp.log1p(adata)
    log.info("  log1p transformation applied")

    # ═══ Phase 2: Store raw counts ═══
    adata.raw = adata.copy()
    log.info("  Raw counts stored in adata.raw")

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

    # ── Save ────────────────────────────────────────────────────────────
    safe_write(adata, output_path, cfg=cfg)
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
