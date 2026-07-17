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
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import scanpy as sc
import squidpy as sq
import numpy as np
import scipy.sparse


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config, modality="spatial")
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

    # ═══ Phase 1: Normalization ═══
    log.info("Normalizing to target sum=%.0f...", CFG.normalization.normalize_target_sum)
    sc.pp.normalize_total(adata, target_sum=CFG.normalization.normalize_target_sum)
    log.info("  Normalization complete")
    sc.pp.log1p(adata)
    log.info("  log1p transformation applied")

    # ═══ Phase 2: Store raw counts ═══
    adata.raw = adata.copy()
    log.info("  Raw counts stored in adata.raw")

    # ═══ Phase 3: Highly variable genes ═══
    flavors = []
    for f in [CFG.hvg.flavor, "seurat_v3", "cell_ranger", "seurat"]:
        if f not in flavors:
            flavors.append(f)

    hvg_selected = False
    last_error = None
    for flavor in flavors:
        try:
            log.info("Selecting %d HVGs (flavor=%s)...", CFG.hvg.n_top_genes, flavor)
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=CFG.hvg.n_top_genes,
                flavor=flavor,
                batch_key=CFG.hvg.batch_key if CFG.has_sample_mapping() else None,
            )
            n_hvg = adata.var['highly_variable'].sum()
            log.info("  Selected %d HVGs (flavor=%s)", n_hvg, flavor)
            hvg_selected = True
            break
        except Exception as e:
            log.warning("  HVG flavor '%s' failed: %s", flavor, e)
            last_error = e

    if not hvg_selected:
        log.info("Falling back to manual variance-based HVG selection...")
        X = adata.X
        if scipy.sparse.issparse(X):
            mean = np.array(X.mean(axis=0)).flatten()
            var = np.array(X.multiply(X).mean(axis=0)).flatten() - mean ** 2
        else:
            var = X.var(axis=0)
        top_idx = np.argsort(var)[::-1][:CFG.hvg.n_top_genes]
        adata.var['highly_variable'] = False
        adata.var.iloc[top_idx, adata.var.columns.get_loc('highly_variable')] = True
        n_hvg = len(top_idx)
        log.info("  Selected %d HVGs (manual variance fallback)", n_hvg)
        hvg_selected = True

    # Force-include genes from CFG.hvg.forced_genes and CFG.marker.marker_dict
    forced_genes = set(CFG.hvg.forced_genes)
    if CFG.marker.marker_dict:
        for genes in CFG.marker.marker_dict.values():
            forced_genes.update(genes)
    genes_in_data = set(adata.var_names)
    force_set = forced_genes & genes_in_data
    if force_set:
        adata.var.loc[list(force_set), 'highly_variable'] = True
        log.info("  Force-included %d marker/forced genes as HVG", len(force_set))

    n_hvg_final = adata.var['highly_variable'].sum()
    log.info("  Total HVGs after force-include: %d", n_hvg_final)

    if n_hvg_final == 0:
        log.error("No HVGs selected — check your data quality")
        sys.exit(1)

    # ═══ Phase 4: Spatial neighbor graph ═══
    log.info("Building spatial neighbor graph...")
    sq.gr.spatial_neighbors(
        adata,
        n_neighs=CFG.spatial.neighbors_n,
        radius=CFG.spatial.neighbors_radius,
        coord_type='generic',
    )
    log.info("  Spatial graph: n_neighs=%d, radius=%.1f", CFG.spatial.neighbors_n, CFG.spatial.neighbors_radius)

    if 'spatial_connectivities' not in adata.obsp:
        log.error("Spatial neighbor graph NOT created — spatial_connectivities missing")
        sys.exit(1)

    # ═══ Phase 5: PCA ═══
    log.info("Computing PCA (n_comps=%d)...", CFG.pca.n_pcs_use)
    sc.pp.pca(
        adata,
        n_comps=CFG.pca.n_pcs_use,
        use_highly_variable=True,
        svd_solver='arpack',
        random_state=CFG.execution.random_seed,
    )
    log.info("  PCA complete: %d components stored", CFG.pca.n_pcs_use)

    # ── Save ────────────────────────────────────────────────────────────
    safe_write(adata, output_path, cfg=CFG)
    log.info("Step 03 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
