#!/usr/bin/env python3
"""
Step 09: ATAC pseudotime trajectory
======================================
  - Reads annotated AnnData (05_annotated.h5ad)
  - Path A (no cfg.atac.terminal_cell_types): explicit skip — no fake
    pseudotime, no misleading UMAP; checkpoint carries a status annotation.
  - Path B (terminal_cell_types set): scanpy diffusion pseudotime (diffmap +
    dpt) rooted at cells farthest from all terminal-type centroids.

Input:  05_annotated.h5ad
Output: 09_trajectory.h5ad
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import snapatac2 as snap

from core.utils import resolve_config, safe_plot, safe_write, setup_logger


def _path_a(cfg, log, data, reason="no terminal_cell_types configured"):
    """Explicit skip: annotate the checkpoint, do not fabricate pseudotime."""
    log.info("Pseudotime skipped: %s", reason)
    # Materialize so safe_write's anndata compression kwargs are accepted.
    if data.isbacked:
        data = data.to_memory()
    data.uns["trajectory"] = {"status": "skipped", "reason": reason}
    return data


def _path_b(cfg, log, data):
    """Diffusion pseudotime rooted away from terminal cell types."""
    terminal = list(cfg.atac.terminal_cell_types)
    log.info(
        "Pseudotime: diffusion DPT, terminal_cell_types=%s",
        terminal,
    )
    if "cell_type" not in data.obs:
        return _path_a(cfg, log, data, reason="diffusion pseudotime failed: no cell_type column")
    cell_types = data.obs["cell_type"].astype(str)
    present = [t for t in terminal if (cell_types == t).any()]
    if not present:
        return _path_a(
            cfg,
            log,
            data,
            reason=(
                "diffusion pseudotime failed: no terminal_cell_types match "
                f"obs['cell_type'] (available: {sorted(set(cell_types))[:10]})"
            ),
        )

    import scanpy as sc

    if data.isbacked:
        data = data.to_memory()
    if "X_spectral" not in data.obsm:
        return _path_a(
            cfg, log, data, reason="diffusion pseudotime failed: no X_spectral embedding"
        )
    # snap.pp.knn only writes obsp['distances']; scanpy diffmap/dpt need
    # connectivities, so rebuild the graph on the spectral embedding.
    if "connectivities" not in data.obsp:
        sc.pp.neighbors(data, n_neighbors=cfg.clustering.n_neighbors, use_rep="X_spectral")
    sc.tl.diffmap(data, n_comps=15)

    # Root: cells farthest from ALL terminal-type centroids in diffmap space.
    xd = data.obsm["X_diffmap"]
    centroids = np.stack([xd[cell_types.values == t].mean(axis=0) for t in present])
    dist = np.min(np.linalg.norm(xd[:, None, :] - centroids[None, :, :], axis=2), axis=1)
    k = max(1, int(0.01 * data.n_obs))
    far_idx = np.argsort(-dist)[:k]
    root_idx = int(far_idx[np.argsort(dist[far_idx])[len(far_idx) // 2]])

    # scanpy 1.12 removed root_user — set uns['iroot'] instead.
    data.uns["iroot"] = root_idx
    sc.tl.dpt(data, n_dcs=10)
    # dpt emits inf for cells unreachable from the root in the diffusion graph
    # (typically disconnected components); sanitise to NaN.
    pseudo = data.obs["dpt_pseudotime"].astype(float).replace([np.inf, -np.inf], np.nan)
    data.obs["pseudotime"] = pseudo.values
    n_nan = int(pseudo.isna().sum())
    log.info(
        "DPT done: pseudotime range [%.3f, %.3f], %d NaN (unreachable)",
        pseudo.min() if not pseudo.isna().all() else float("nan"),
        pseudo.max() if not pseudo.isna().all() else float("nan"),
        n_nan,
    )

    try:
        safe_plot(
            snap.pl.umap,
            data,
            color="pseudotime",
            cmap=cfg.plot.palette.pseudotime,
            show=False,
            save=os.path.join(cfg.figure_dir, "09_trajectory", "trajectory_pseudotime"),
            cfg=cfg,
        )
    except Exception as e:
        log.warning("Pseudotime UMAP plot failed: %s", e)

    data.uns["trajectory"] = {
        "status": "computed",
        "root": int(root_idx),
        "terminal_cell_types": present,
    }
    return data


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("09_trajectory", os.path.join(cfg.log_dir, "09_trajectory.log"))
    log.info("Step 09: ATAC pseudotime")

    if os.path.exists(cfg.trajectory_h5ad):
        log.info("Skip: %s exists.", cfg.trajectory_h5ad)
        return

    data = snap.read(cfg.annotated_h5ad)
    log.info("Loaded: %d cells (backed mode)", data.n_obs)

    terminal = list(getattr(cfg.atac, "terminal_cell_types", []) or [])
    if not terminal:
        data = _path_a(cfg, log, data)
    else:
        try:
            data = _path_b(cfg, log, data)
        except Exception as e:
            # Scanpy / connectivity failure -> degrade to explicit skip.
            log.warning("Diffusion pseudotime failed: %s", e)
            data = _path_a(cfg, log, data, reason=f"diffusion pseudotime failed: {e}")

    safe_write(data, cfg.trajectory_h5ad, cfg=cfg)
    log.info("Step 09 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
