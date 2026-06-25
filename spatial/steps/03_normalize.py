#!/usr/bin/env python3
"""
Step 03: Normalization + HVG selection + spatial neighbor graph
==================================================================
  1. Library-size normalize + log1p
  2. Identify highly variable genes (HVG)
  3. Build spatial neighbor graph via sq.gr.spatial_neighbors()
  4. PCA on HVGs

Input:  02_image.h5ad (or 01_qc.h5ad if Step 02 was skipped)
Output: 03_processed.h5ad
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import scanpy as sc
import squidpy as sq


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("03_normalize", os.path.join(CFG.log_dir, "03_normalize.log"))
    log.info("Step 03: Normalization + HVG + spatial graph + PCA")

    output_path = os.path.join(CFG.h5ad_dir, "03_processed.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    # ── Determine input ─────────────────────────────────────────────────
    image_path = os.path.join(CFG.h5ad_dir, "02_image.h5ad")
    qc_path = os.path.join(CFG.h5ad_dir, "01_qc.h5ad")
    if os.path.exists(image_path):
        input_path = image_path
    elif os.path.exists(qc_path):
        input_path = qc_path
    else:
        log.error("Neither %s nor %s found. Run Steps 01–02 first.", image_path, qc_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    # ── 1. Normalize ────────────────────────────────────────────────────
    log.info("Normalizing to target sum=%.0f...", CFG.normalize_target_sum)
    sc.pp.normalize_total(adata, target_sum=CFG.normalize_target_sum)
    log.info("  Normalization complete")

    # ── 2. Log1p ─────────────────────────────────────────────────────────
    sc.pp.log1p(adata)
    log.info("  log1p transformation applied")

    # Save raw counts for later use (e.g., marker gene scoring)
    adata.raw = adata.copy()
    log.info("  Raw counts stored in adata.raw")

    # ── 3. Highly variable genes ────────────────────────────────────────
    log.info("Selecting %d HVGs (flavor=%s)...", CFG.n_top_genes, CFG.hvg_flavor)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=CFG.n_top_genes,
        flavor=CFG.hvg_flavor,
        batch_key=CFG.hvg_batch_key if CFG.has_sample_mapping() else None,
    )
    n_hvg = adata.var['highly_variable'].sum()
    log.info("  Selected %d HVGs", n_hvg)

    if n_hvg == 0:
        log.error("No HVGs selected — check your data quality")
        sys.exit(1)

    # ── 4. Spatial neighbor graph ────────────────────────────────────────
    log.info("Building spatial neighbor graph...")
    # Determine connectivity mode
    if CFG.spatial_neighbors_radius > 0:
        sq.gr.spatial_neighbors(
            adata,
            radius=CFG.spatial_neighbors_radius,
            library_key=None,
            coord_type='generic',
        )
        log.info("  Spatial graph: radius=%.1f", CFG.spatial_neighbors_radius)
    else:
        sq.gr.spatial_neighbors(
            adata,
            n_neighs=CFG.spatial_neighbors_n,
            library_key=None,
            coord_type='generic',
        )
        log.info("  Spatial graph: n_neighs=%d", CFG.spatial_neighbors_n)

    # Verify graph exists
    if 'spatial_connectivities' not in adata.obsp:
        log.error("Spatial neighbor graph NOT created — spatial_connectivities missing")
        sys.exit(1)

    # ── 5. PCA ────────────────────────────────────────────────────────────
    log.info("Computing PCA (n_comps=%d)...", CFG.n_pcs_full)
    sc.pp.pca(
        adata,
        n_comps=CFG.n_pcs_full,
        use_highly_variable=True,
        svd_solver='arpack',
        random_state=CFG.random_seed,
    )
    n_pcs_used = min(CFG.n_pcs_use, adata.obsm['X_pca'].shape[1])
    log.info("  PCA complete: %d components stored", n_pcs_used)

    # ── Save ────────────────────────────────────────────────────────────
    safe_write(adata, output_path, cfg=CFG)
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
