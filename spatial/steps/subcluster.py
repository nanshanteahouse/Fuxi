#!/usr/bin/env python3
"""
Spatial subcluster: conditional subclustering of a selected cell type
======================================================================
  Ported from rna/steps/06_subcluster.py to spatial transcriptomics.

  Extracts a user-selected cell type, re-runs HVG → PCA → Leiden at
  multiple resolutions, optionally annotates subclusters via AI, and
  optionally loads scRNA markers via CFG.rna_ref.

Input:  05_annotated.h5ad  (from Step 05, or downstream)
Output: 05_sub_{cell_type}.h5ad  (per cell type, in CFG.h5ad_dir)
"""
import sys, os, time, argparse, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write, safe_plot
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def resolve_cell_type(cell_types, requested, cutoff=0.6):
    if requested in cell_types:
        return requested
    casefold_map = {ct.casefold(): ct for ct in cell_types}
    folded = requested.casefold()
    if folded in casefold_map:
        return casefold_map[folded]
    import difflib
    matches = difflib.get_close_matches(requested, cell_types, n=1, cutoff=cutoff)
    return matches[0] if matches else None


def auto_writeback(sub, cell_type, main_path, log=None, cfg=None):
    if log is None:
        log = setup_logger(
            "subcluster_writeback",
            os.path.join(os.path.dirname(main_path) or ".",
                         "subcluster_writeback.log"),
        )
    if 'sub_ai_label' not in sub.obs:
        log.info("No sub_ai_label on subcluster — skipping writeback.")
        return 0
    if not os.path.exists(main_path):
        log.warning("Main annotation file not found, skipping writeback: %s", main_path)
        return 0
    main = sc.read(main_path)
    if 'cell_type' not in main.obs:
        log.warning("'cell_type' missing in %s — skipping writeback.", main_path)
        return 0
    if 'cell_subtype' not in main.obs:
        main.obs['cell_subtype'] = main.obs['cell_type'].astype(str)
    main.obs['cell_subtype'] = main.obs['cell_subtype'].astype(str)
    sub_labels = sub.obs['sub_ai_label'].astype(str)
    n_updated = 0
    for bc, label in sub_labels.items():
        if bc not in main.obs_names:
            continue
        if main.obs.at[bc, 'cell_type'] != cell_type:
            continue
        main.obs.at[bc, 'cell_subtype'] = label
        n_updated += 1
    if n_updated == 0:
        log.warning("No cells in subcluster matched %s in main file — "
                    "skipping writeback.", cell_type)
        return 0
    main.obs['cell_subtype'] = main.obs['cell_subtype'].astype('category')
    safe_write(main, main_path, cfg=cfg)
    log.info("Wrote back %d subcluster labels to %s", n_updated, main_path)
    return n_updated


def main():
    t0 = time.time()

    parser = argparse.ArgumentParser(
        description="Subcluster a selected cell type from spatial annotated data."
    )
    parser.add_argument("--config", default="../config.py",
                        help="Path to config.py")
    parser.add_argument("--cell-type", default=None,
                        help="Cell type to extract and subcluster")
    args = parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("subcluster",
                       os.path.join(CFG.log_dir, "subcluster.log"))
    log.info("Spatial subcluster step")
    log.info("Cell type: %s", args.cell_type)

    if args.cell_type is None and not CFG.marker.subcluster_types:
        log.info("subcluster_types not configured, skipping.")
        sys.exit(2)

    input_path = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_path):
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d spots, %d genes",
             input_path, adata.n_obs, adata.n_vars)

    if 'cell_type' not in adata.obs:
        log.error("'cell_type' column not found in adata.obs. "
                  "Run annotation step first.")
        sys.exit(1)

    available_types = sorted(adata.obs['cell_type'].unique().tolist())
    resolved_type = resolve_cell_type(available_types, args.cell_type)
    if resolved_type is None:
        log.error(
            "Cell type '%s' not found in adata.obs['cell_type']. "
            "Available types: %s",
            args.cell_type, available_types,
        )
        sys.exit(1)
    if resolved_type != args.cell_type:
        log.info("Resolved cell type: '%s' -> '%s'", args.cell_type, resolved_type)
    args.cell_type = resolved_type

    mask = adata.obs['cell_type'] == args.cell_type
    n_cells = mask.sum()
    sub = adata[mask].copy()
    log.info("Filtered: %d spots for cell type '%s'", n_cells, args.cell_type)

    min_cells = CFG.marker.min_cells_subcluster
    if n_cells < min_cells:
        log.warning("Too few spots (%d < %d). Skipping subclustering.",
                    n_cells, min_cells)
        safe_cell_type = args.cell_type.replace(" ", "_").replace("/", "_")
        output_path = os.path.join(
            CFG.h5ad_dir, f"05_sub_{safe_cell_type}.h5ad")
        safe_write(sub, output_path, cfg=CFG)
        log.info("Saved subset (no subclustering performed): %s", output_path)
        return

    processed_path = os.path.join(CFG.h5ad_dir, "03_processed.h5ad")
    if os.path.exists(processed_path):
        raw_adata = sc.read(processed_path)
        common_bc = list(sub.obs_names.intersection(raw_adata.obs_names))
        n_common = len(common_bc)
        if n_common < n_cells:
            log.info("Matched %d / %d spots in 03_processed.h5ad", n_common, n_cells)

        if n_common >= min_cells:
            sub = sub[common_bc].copy()
            sub_raw = raw_adata[common_bc].copy()
            del raw_adata
            log.info("Raw subset loaded: %d spots, %d genes",
                     sub_raw.n_obs, sub_raw.n_vars)

            sc.pp.normalize_total(sub_raw, target_sum=CFG.normalization.normalize_target_sum)
            sc.pp.log1p(sub_raw)

            hvg_flavors = [CFG.hvg.flavor, 'seurat_v3', 'seurat']
            hvg_ok = False
            for flavor in hvg_flavors:
                try:
                    sc.pp.highly_variable_genes(
                        sub_raw, n_top_genes=CFG.hvg.n_top_genes, flavor=flavor,
                    )
                    hvg_ok = True
                    log.info("Subset HVG: flavor=%s, %d genes selected",
                             flavor, sub_raw.var['highly_variable'].sum())
                    break
                except (ValueError, ImportError):
                    continue

            if hvg_ok and sub_raw.var['highly_variable'].sum() > 0:
                sub_raw = sub_raw[:, sub_raw.var['highly_variable']].copy()
                n_comps_sub = min(50, sub_raw.n_obs - 2)
                sc.pp.pca(sub_raw, n_comps=n_comps_sub, svd_solver='randomized',
                          random_state=CFG.execution.random_seed)
                log.info("Subset PCA: n_comps=%d", n_comps_sub)
                sub.obsm['X_pca'] = sub_raw.obsm['X_pca']
            else:
                log.warning("Subset HVG failed, falling back to parent-level PCA")
                n_comps_sub = min(50, sub.n_obs - 2)
                sc.pp.pca(sub, n_comps=n_comps_sub, svd_solver='arpack')
        else:
            log.warning("Too few overlapping spots (%d < %d), "
                        "falling back to parent-level PCA",
                        n_common, min_cells)
            n_comps_sub = min(50, sub.n_obs - 2)
            sc.pp.pca(sub, n_comps=n_comps_sub, svd_solver='arpack')
    else:
        log.warning("03_processed.h5ad not found at %s — falling back to "
                    "parent-level PCA", processed_path)
        n_comps_sub = min(50, sub.n_obs - 2)
        sc.pp.pca(sub, n_comps=n_comps_sub, svd_solver='arpack')

    n_pcs_use = min(CFG.pca.n_pcs_use, n_comps_sub)
    log.info("Computing neighbor graph (n_pcs=%d)...", n_pcs_use)
    sc.pp.neighbors(sub, n_pcs=n_pcs_use,
                    random_state=CFG.execution.random_seed)

    log.info("Computing UMAP...")
    sc.tl.umap(sub, random_state=CFG.execution.random_seed)

    log.info("Leiden subclustering, resolutions: %s", CFG.clustering.leiden_resolutions)
    for res in CFG.clustering.leiden_resolutions:
        key = f'sub_leiden_r{res}'
        sc.tl.leiden(sub, resolution=res, key_added=key,
                     random_state=CFG.execution.random_seed, flavor=CFG.clustering.leiden_flavor,
                     directed=False, n_iterations=2)
        n_cl = sub.obs[key].nunique()
        log.info("  r=%.1f -> %d subclusters", res, n_cl)

    best_key = f'sub_leiden_r{CFG.marker.subcluster_resolution}'
    if best_key in sub.obs:
        sub.obs['leiden'] = sub.obs[best_key].copy()
        log.info("  Subcluster resolution: sub_leiden_r%.1f -> leiden set",
                 CFG.marker.subcluster_resolution)
    else:
        fallback_key = f'sub_leiden_r{CFG.clustering.best_resolution}'
        if fallback_key in sub.obs:
            sub.obs['leiden'] = sub.obs[fallback_key].copy()
            log.info("  Fallback to best_resolution: %.1f -> leiden set",
                     CFG.clustering.best_resolution)
        else:
            avail = [k for k in sub.obs if k.startswith('sub_leiden_')]
            if avail:
                sub.obs['leiden'] = sub.obs[avail[-1]].copy()
                log.info("  Fallback to %s for 'leiden'", avail[-1])
            else:
                sub.obs['leiden'] = "0"

    fig_dir = os.path.join(CFG.figure_dir, 'subcluster')
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    sc.settings.autoshow = False
    safe_cell_type = args.cell_type.replace(" ", "_").replace("/", "_")
    safe_filename = f"05_sub_{safe_cell_type}.h5ad"
    output_path = os.path.join(CFG.h5ad_dir, safe_filename)

    safe_plot(sc.pl.umap, sub, color='leiden', show=False,
              save=f'_sub_{safe_cell_type}_leiden.pdf',
              title=f'{args.cell_type} — leiden')

    res_keys = [f'sub_leiden_r{r}' for r in CFG.clustering.leiden_resolutions
                if f'sub_leiden_r{r}' in sub.obs]
    i = -1
    if res_keys:
        n_res = len(res_keys)
        n_cols = min(3, n_res)
        n_rows = int(np.ceil(n_res / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(6 * n_cols, 5 * n_rows))
        axes = axes.ravel() if n_res > 1 else [axes]
        for i, key in enumerate(res_keys):
            sc.pl.umap(sub, color=key, ax=axes[i], show=False,
                       legend_fontsize=8,
                       title=f'{safe_cell_type} — {key}')
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir,
                                 f'umap_sub_{safe_cell_type}_resolutions.pdf'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info("Multi-resolution UMAP saved for %s", safe_cell_type)

    safe_write(sub, output_path, cfg=CFG)
    log.info("Intermediate results saved (pre-AI): %s", output_path)

    if CFG.ai.enabled and CFG.ai.subcluster:
        try:
            from core.ai.prompts import build_annotation_prompt
            from core.ai.caller import ai_query

            log.info("Running AI subcluster re-annotation...")
            sys_prompt, user_prompt = build_annotation_prompt(
                sub, tissue=args.cell_type, species="unknown",
                precomputed_rank=True,
            )

            result = ai_query(sys_prompt, user_prompt, CFG.ai)

            if not result:
                log.warning("AI returned empty response — "
                            "falling back to numeric labels")
                sub.obs['sub_ai_label'] = (
                    "Subcluster_" + sub.obs['leiden'].astype(str)
                ).astype('category')
            else:
                log.info("AI response received (%d chars)", len(result))
                cleaned = result.strip()
                if cleaned.startswith("```"):
                    lines = cleaned.split("\n")
                    start = 0
                    for i, line in enumerate(lines):
                        if line.strip().startswith("```"):
                            start = i + 1
                            break
                    end = len(lines)
                    for i in range(len(lines) - 1, start - 1, -1):
                        if lines[i].strip().startswith("```"):
                            end = i
                            break
                    cleaned = "\n".join(lines[start:end]).strip()

                parsed = json.loads(cleaned)
                ai_labels = {}
                for cluster_id, info in parsed.items():
                    cell_type = info.get("cell_type", "Unknown")
                    subtype = info.get("subtype", "N/A")
                    if subtype and subtype.upper() != "N/A":
                        label = f"{cell_type}_{subtype}"
                    else:
                        label = cell_type
                    label = label.replace(" ", "_").replace("/", "_")
                    ai_labels[cluster_id] = label

                sub.obs['leiden'] = sub.obs['leiden'].astype(str)
                sub.obs['sub_ai_label'] = (
                    sub.obs['leiden'].map(ai_labels)
                ).astype('category')
                n_ai_types = sub.obs['sub_ai_label'].nunique()
                log.info("AI annotation: %d subcluster types identified",
                         n_ai_types)

                for cluster_id in sorted(sub.obs['leiden'].unique(),
                                         key=lambda x: int(x)):
                    label = ai_labels.get(str(cluster_id), "Unmapped")
                    count = (sub.obs['leiden'] == cluster_id).sum()
                    log.info("  Subcluster %s -> %s (%d spots)",
                             cluster_id, label, count)

                safe_plot(sc.pl.umap, sub, color='sub_ai_label', show=False,
                          save=f'_sub_{safe_cell_type}_ai.pdf',
                          title=f'{args.cell_type} — AI subcluster')

        except Exception as e:
            log.warning("AI subcluster annotation failed: %s", e)
            log.info("Falling back to numeric subcluster labels.")
            sub.obs['sub_ai_label'] = (
                "Subcluster_" + sub.obs['leiden'].astype(str)
            ).astype('category')
    else:
        log.info("AI subcluster annotation disabled "
                 "(CFG.ai.enabled=%s, CFG.ai.subcluster=%s)",
                 CFG.ai.enabled, CFG.ai.subcluster)

    sub.obs['leiden'] = sub.obs['leiden'].astype('category')
    safe_write(sub, output_path, cfg=CFG)
    log.info("Saved: %s", output_path)

    main_path = os.path.join(CFG.h5ad_dir, "05_annotated.h5ad")
    auto_writeback(
        sub=sub,
        cell_type=args.cell_type,
        main_path=main_path,
        log=log,
        cfg=CFG,
    )

    n_clusters = sub.obs['leiden'].nunique()
    log.info("=" * 50)
    log.info("Subcluster Summary")
    log.info("  Cell type:       %s", args.cell_type)
    log.info("  Spots:           %d", n_cells)
    log.info("  Subclusters:     %d", n_clusters)
    log.info("  Resolution:      %.1f", CFG.marker.subcluster_resolution)
    log.info("  Per-cluster counts:")
    for cluster_id in sorted(sub.obs['leiden'].unique(),
                             key=lambda x: int(x)):
        count = (sub.obs['leiden'] == cluster_id).sum()
        pct = 100.0 * count / n_cells
        log.info("    Cluster %s: %d spots (%.1f%%)",
                 cluster_id, count, pct)
    log.info("=" * 50)
    log.info("Spatial subcluster done, %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
