#!/usr/bin/env python3
"""
Step 07: Pseudotime trajectory analysis for spatial transcriptomics
======================================================================
  - PAGA graph abstraction
  - Diffusion pseudotime (DPT)
  - Optional: CellRank for fate mapping

Input:  05_annotated.h5ad
Output: 07_trajectory.h5ad
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("07_trajectory", os.path.join(CFG.log_dir, "07_trajectory.log"))
    log.info("Step 07: Trajectory analysis (PAGA + DPT)")

    output_path = os.path.join(CFG.h5ad_dir, "07_trajectory.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return

    input_path = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_path):
        log.error("Input not found: %s. Run Step 05 first.", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots × %d genes", input_path, adata.n_obs, adata.n_vars)

    # Check we have PCA and neighbors
    if 'X_pca' not in adata.obsm:
        log.error("No PCA found — run Step 03 first")
        sys.exit(1)

    if 'neighbors' not in adata.uns:
        log.warning("No neighbor graph found, building default neighbors...")
        sc.pp.neighbors(adata, n_pcs=CFG.n_pcs_use, random_state=CFG.random_seed)

    # ── 1. PAGA graph ────────────────────────────────────────────────────
    log.info("Computing PAGA graph...")
    group_col = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
    log.info("  Grouping by '%s'", group_col)

    try:
        sc.tl.paga(adata, groups=group_col)
        log.info("  PAGA complete")

        # PAGA plot
        safe_plot(sc.pl.paga, adata, show=False,
                  save='_07_paga.png', threshold=0.1)
        safe_plot(sc.pl.paga_compare, adata, show=False,
                  save='_07_paga_compare.png',
                  legend_fontsize=8)
        log.info("  PAGA plots saved")
    except Exception as e:
        log.warning("PAGA failed: %s — skipping", e)

    # ── 2. Diffusion pseudotime (DPT) ────────────────────────────────────
    log.info("Computing diffusion pseudotime...")

    # Try root cell type selection
    root_cell_types = CFG.root_cell_types
    try:
        if root_cell_types:
            adata.uns['iroot'] = np.flatnonzero(
                adata.obs[group_col].isin(root_cell_types)
            )[0]
            log.info("  Root set from root_cell_types: %s", root_cell_types)
        elif CFG.root_markers:
            # Find cells with highest expression of root markers
            root_markers_present = [g for g in CFG.root_markers if g in adata.var_names]
            if root_markers_present:
                root_score = adata[:, root_markers_present].X.mean(axis=1)
                if scipy.sparse.issparse(root_score):
                    root_score = root_score.toarray()
                adata.uns['iroot'] = int(np.argmax(root_score))
                log.info("  Root auto-detected from markers: %s → cell %d",
                         root_markers_present[:3], adata.uns['iroot'])
            else:
                log.warning("No root markers found in data")
        else:
            # Auto-select root: first cluster (typically stem/progenitor)
            root_cluster = sorted(adata.obs[group_col].unique())[0]
            root_cells = adata.obs[group_col] == root_cluster
            adata.uns['iroot'] = int(np.flatnonzero(root_cells.values)[0])
            log.info("  Root auto-selected from first cluster '%s'", root_cluster)
    except (ValueError, IndexError) as e:
        log.warning("Could not set root cell: %s — using first cell", e)
        adata.uns['iroot'] = 0

    try:
        sc.tl.diffmap(adata, n_comps=CFG.n_diffmap_comps, random_state=CFG.random_seed)
        log.info("  Diffusion map computed (%d components)", CFG.n_diffmap_comps)
        sc.tl.dpt(adata, n_branchings=CFG.n_branchings)
        log.info("  DPT computed")
    except Exception as e:
        log.warning("DPT failed: %s", e)

    # ── 3. Visualizations ────────────────────────────────────────────────
    sc.settings.figdir = CFG.figure_dir
    sc.settings.autoshow = False

    plot_vars = ['dpt_pseudotime'] if 'dpt_pseudotime' in adata.obs else []
    if group_col in adata.obs:
        plot_vars.insert(0, group_col)

    for var in plot_vars:
        try:
            safe_plot(sc.pl.umap, adata, color=var, show=False,
                      save=f'_07_{var}.pdf')
        except Exception:
            pass

    # ── Save ─────────────────────────────────────────────────────────────
    safe_write(adata, output_path, cfg=CFG)
    log.info("Step 07 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
