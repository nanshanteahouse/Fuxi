#!/usr/bin/env python3
"""
Step 06: Interactive subclustering of a selected cell type
===========================================================
  Extract a user-selected cell type, re-run PCA + neighbors + UMAP + Leiden,
  and optionally use AI to re-annotate subclusters.

  Designed to be called in a loop (once per cell type) for iterative
  refinement of subpopulations within major cell types.

Input:  05_annotated.h5ad  (from Step 05)
Output: 05_sub_{cell_type}.h5ad  (per cell type, in CFG.h5ad_dir)
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import json

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from core.ai.json_extract import extract_json_block
from core.utils import resolve_config, safe_plot, safe_write, setup_logger


def resolve_cell_type(cell_types, requested, cutoff=0.6):
    """
    Resolve a user-supplied cell-type name against an available list,
    with a forgiving three-tier fallback:

        1. exact membership              — "T cell"   → "T cell"
        2. case-insensitive (casefold)   — "t cell"   → "T cell"
        3. fuzzy (difflib, cutoff=0.6)   — "Müller Glia" → "Müller glia"

    ``str.casefold()`` is used (not ``str.lower()``) so non-ASCII
    case-folding behaves correctly on ß, İ, ı, etc.

    Returns the matched canonical cell-type string, or ``None`` if no
    tier succeeded. ``difflib`` is imported lazily so the cold-start
    cost of this module is unaffected for callers who never hit the
    fuzzy tier.
    """
    if requested in cell_types:
        return requested

    casefold_map = {ct.casefold(): ct for ct in cell_types}
    folded = requested.casefold()
    if folded in casefold_map:
        return casefold_map[folded]

    import difflib  # local import — fuzzy tier is the slow path

    matches = difflib.get_close_matches(requested, cell_types, n=1, cutoff=cutoff)
    return matches[0] if matches else None


def auto_writeback(sub, cell_type, main_path, log=None, cfg=None):
    """
    Merge ``sub.obs['sub_ai_label']`` back into ``main_path`` as the
    ``cell_subtype`` column, scoped to cells whose ``cell_type`` matches
    ``cell_type``.

    No-op (returns 0) when:
      - ``sub_ai_label`` is absent (AI annotation was not run, or
        ``CFG.ai.ai_subcluster`` was disabled)
      - ``main_path`` does not exist on disk
      - ``main_path`` is missing ``cell_type`` / ``cell_subtype`` columns
      - no cells in ``sub`` are present in the main file

    Returns the number of cells whose ``cell_subtype`` was updated.
    """
    if log is None:
        log = setup_logger(
            "06_subcluster_writeback",
            os.path.join(os.path.dirname(main_path) or ".", "06_subcluster_writeback.log"),
        )

    if "sub_ai_label" not in sub.obs:
        log.info("No sub_ai_label on subcluster — skipping writeback.")
        return 0

    if not os.path.exists(main_path):
        log.warning("Main annotation file not found, skipping writeback: %s", main_path)
        return 0

    main = sc.read(main_path)
    if "cell_type" not in main.obs:
        log.warning("'cell_type' missing in %s — skipping writeback.", main_path)
        return 0

    if "cell_subtype" not in main.obs:
        # Match the score_genes fallback: cell_subtype starts as a copy
        # of cell_type (see 05_annotate_major.py score_genes_mode).
        main.obs["cell_subtype"] = main.obs["cell_type"].astype(str)

    # h5ad round-trips object columns through categorical, so the column
    # on disk has a fixed set of categories. Casting to str first lets
    # us add new AI-derived labels (e.g. "Tcell_Naive") without having
    # to enumerate them in advance via cat.add_categories.
    main.obs["cell_subtype"] = main.obs["cell_subtype"].astype(str)

    sub_labels = sub.obs["sub_ai_label"].astype(str)
    n_updated = 0
    for bc, label in sub_labels.items():
        if bc not in main.obs_names:
            continue
        if main.obs.at[bc, "cell_type"] != cell_type:
            continue
        main.obs.at[bc, "cell_subtype"] = label
        n_updated += 1

    if n_updated == 0:
        log.warning(
            "No cells in subcluster matched %s in main file — skipping writeback.", cell_type
        )
        return 0

    main.obs["cell_subtype"] = main.obs["cell_subtype"].astype("category")
    safe_write(main, main_path, cfg=cfg)
    log.info("Wrote back %d subcluster labels to %s", n_updated, main_path)
    return n_updated


def main():
    t0 = time.time()

    # ── Argument parsing ──────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        description="Subcluster a selected cell type from annotated data."
    )
    parser.add_argument("--config", default="../config.py", help="Path to config.py")
    parser.add_argument(
        "--cell-type",
        default=None,
        help="Cell type to extract and subcluster (e.g. 'T cell', 'Microglia')",
    )
    args = parser.parse_args()

    # ── Config & logger ───────────────────────────────────────────────
    cfg = resolve_config(args.config)
    log = setup_logger("06_subcluster", os.path.join(cfg.log_dir, "06_subcluster.log"))
    log.info("Step 06: Interactive subclustering")
    log.info("Cell type: %s", args.cell_type)

    # Early exit when neither --cell-type nor CFG.marker.subcluster_types is configured
    if args.cell_type is None and not cfg.marker.subcluster_types:
        log.info("subcluster_types not configured, skipping.")
        sys.exit(2)

    # Fallback: use first configured subcluster type when --cell-type is not provided
    if args.cell_type is None and cfg.marker.subcluster_types:
        args.cell_type = cfg.marker.subcluster_types[0]
        log.info("Using configured subcluster type: %s", args.cell_type)

    # ── (a) Load annotated data ───────────────────────────────────────
    input_path = os.path.join(cfg.h5ad_dir, "05_annotated.h5ad")
    if not os.path.exists(input_path):
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %s — %d cells, %d genes", input_path, adata.n_obs, adata.n_vars)

    if "cell_type" not in adata.obs:
        log.error("'cell_type' column not found in adata.obs. Run Step 05 (annotate) first.")
        sys.exit(1)

    # ── (b) Filter to selected cell type ──────────────────────────────
    # Strict `==` would silently yield 0 cells on minor typos, so we
    # resolve via resolve_cell_type() and update args in place.
    available_types = sorted(adata.obs["cell_type"].unique().tolist())
    resolved_type = resolve_cell_type(available_types, args.cell_type)
    if resolved_type is None:
        log.error(
            "Cell type '%s' not found in adata.obs['cell_type']. Available types: %s",
            args.cell_type,
            available_types,
        )
        sys.exit(1)
    if resolved_type != args.cell_type:
        log.info("Resolved cell type: '%s' -> '%s'", args.cell_type, resolved_type)
    args.cell_type = resolved_type

    mask = adata.obs["cell_type"] == args.cell_type
    n_cells = mask.sum()
    sub = adata[mask].copy()
    log.info("Filtered: %d cells for cell type '%s'", n_cells, args.cell_type)

    # ── (c) Minimum cell check ────────────────────────────────────────
    min_cells = 50
    if n_cells < min_cells:
        log.warning("Too few cells (%d < %d). Skipping subclustering.", n_cells, min_cells)
        # Still save the subset (without subcluster fields)
        safe_cell_type = args.cell_type.replace(" ", "_").replace("/", "_")
        output_path = os.path.join(cfg.h5ad_dir, f"05_sub_{safe_cell_type}.h5ad")
        safe_write(sub, output_path, cfg=cfg)
        log.info("Saved subset (no subclustering performed): %s", output_path)
        return

    # ── (d) Re-process from raw counts on subset ──────────────────────
    # Subsetting changes which genes are informative across the selected
    # cell type.  Re-select HVGs, re-normalize, re-run PCA + Harmony so
    # subcluster resolution is driven by relevant variation, not by the
    # full-dataset HVG set.
    raw_path = os.path.join(cfg.h5ad_dir, "02_qc.h5ad")
    if os.path.exists(raw_path):
        raw_adata = sc.read(raw_path)
        common_bc = list(sub.obs_names.intersection(raw_adata.obs_names))
        n_common = len(common_bc)
        if n_common < n_cells:
            log.info("Matched %d / %d cells in 02_qc.h5ad", n_common, n_cells)

        if n_common >= min_cells:
            # Align both objects to the same barcode order
            sub = sub[common_bc].copy()
            sub_raw = raw_adata[common_bc].copy()
            del raw_adata
            log.info("Raw subset loaded: %d cells, %d genes", sub_raw.n_obs, sub_raw.n_vars)

            # Re-select HVGs within the subset
            hvg_flavors = [cfg.hvg.flavor, "seurat_v3", "seurat"]
            hvg_ok = False
            for flavor in hvg_flavors:
                try:
                    sc.pp.highly_variable_genes(
                        sub_raw,
                        n_top_genes=cfg.hvg.n_top_genes,
                        flavor=flavor,
                    )
                    hvg_ok = True
                    log.info(
                        "Subset HVG: flavor=%s, %d genes selected",
                        flavor,
                        sub_raw.var["highly_variable"].sum(),
                    )
                    break
                except (ValueError, ImportError):
                    continue

            if hvg_ok and sub_raw.var["highly_variable"].sum() > 0:
                # Normalize + log1p on HVG subset
                sub_raw = sub_raw[:, sub_raw.var["highly_variable"]].copy()
                sc.pp.normalize_total(sub_raw, target_sum=cfg.normalization.normalize_target_sum)
                sc.pp.log1p(sub_raw)

                # PCA
                n_comps_sub = min(50, sub_raw.n_obs - 2)
                sc.pp.pca(
                    sub_raw,
                    n_comps=n_comps_sub,
                    svd_solver="randomized",
                    random_state=cfg.execution.random_seed,
                )
                log.info("Subset PCA: n_comps=%d", n_comps_sub)

                # Harmony batch correction
                if (
                    cfg.integration.method == "harmony"
                    and cfg.integration.batch_key in sub_raw.obs.columns
                ):
                    import harmonypy as hm

                    n_pcs_use = min(cfg.pca.n_pcs_use, n_comps_sub)
                    log.info(
                        "Subset Harmony (batch_key=%s, n_pcs_use=%d)...",
                        cfg.integration.batch_key,
                        n_pcs_use,
                    )
                    try:
                        ho = hm.run_harmony(
                            sub_raw.obsm["X_pca"][:, :n_pcs_use],
                            sub_raw.obs,
                            vars_use=cfg.integration.batch_key,
                            random_state=cfg.execution.random_seed,
                            max_iter_harmony=cfg.integration.max_iter,
                        )
                        sub_raw.obsm["X_integrated"] = ho.Z_corr
                    except Exception as e:
                        log.warning("Subset Harmony failed (%s), using raw PCA", e)
                        sub_raw.obsm["X_integrated"] = sub_raw.obsm["X_pca"].copy()
                else:
                    log.warning(
                        "Subcluster Harmony skipped — integration method is '%s', not 'harmony'. "
                        "Subcluster-level batch effects may persist.",
                        cfg.integration.method,
                    )
                    sub_raw.obsm["X_integrated"] = sub_raw.obsm["X_pca"].copy()

                # Copy embeddings back to sub (already index-aligned)
                sub.obsm["X_pca"] = sub_raw.obsm["X_pca"]
                sub.obsm["X_integrated"] = sub_raw.obsm["X_integrated"]
            else:
                log.warning("Subset HVG selection failed, falling back to parent-level PCA")
                n_comps_sub = min(50, sub.n_obs - 2)
                sc.pp.pca(sub, n_comps=n_comps_sub, svd_solver="arpack")
        else:
            log.warning(
                "Too few overlapping cells in 02_qc.h5ad (%d < %d), "
                "falling back to parent-level PCA",
                n_common,
                min_cells,
            )
            n_comps_sub = min(50, sub.n_obs - 2)
            sc.pp.pca(sub, n_comps=n_comps_sub, svd_solver="arpack")
    else:
        log.warning("02_qc.h5ad not found at %s — falling back to parent-level PCA", raw_path)
        n_comps_sub = min(50, sub.n_obs - 2)
        sc.pp.pca(sub, n_comps=n_comps_sub, svd_solver="arpack")

    # ── (e) Neighbors (use Harmony-corrected PCA when available) ─────
    n_pcs_use = min(cfg.pca.n_pcs_use, n_comps_sub)
    use_rep = "X_integrated" if "X_integrated" in sub.obsm else "X_pca"
    if use_rep in sub.obsm:
        n_pcs_use = min(n_pcs_use, sub.obsm[use_rep].shape[1])
    use_rep = "X_integrated" if "X_integrated" in sub.obsm else "X_pca"
    log.info("Computing neighbor graph (use_rep=%s, n_pcs=%d)...", use_rep, n_pcs_use)
    sc.pp.neighbors(sub, n_pcs=n_pcs_use, use_rep=use_rep, random_state=cfg.execution.random_seed)

    # ── (f) UMAP ──────────────────────────────────────────────────────
    log.info("Computing UMAP...")
    sc.tl.umap(sub, random_state=cfg.execution.random_seed)

    # ── (g) Multi-resolution Leiden ───────────────────────────────────
    log.info("Leiden subclustering, resolutions: %s", cfg.clustering.param_grid_resolutions)
    for res in cfg.clustering.param_grid_resolutions:
        key = f"sub_leiden_r{res}"
        sc.tl.leiden(
            sub,
            resolution=res,
            key_added=key,
            random_state=cfg.execution.random_seed,
            flavor=cfg.clustering.leiden_flavor,
            directed=False,
            n_iterations=2,
        )
        n_cl = sub.obs[key].nunique()
        log.info("  r=%.1f → %d subclusters", res, n_cl)

    # ── (h) Set best-resolution leiden ─────────────────────────────────
    best_key = f"sub_leiden_r{cfg.marker.subcluster_resolution}"
    if best_key in sub.obs:
        sub.obs["leiden"] = sub.obs[best_key].copy()
        log.info(
            "  Subcluster resolution: sub_leiden_r%.1f → leiden set",
            cfg.marker.subcluster_resolution,
        )
    else:
        # Fallback: try best_resolution, then last available
        fallback_key = f"sub_leiden_r{cfg.clustering.best_resolution}"
        if fallback_key in sub.obs:
            sub.obs["leiden"] = sub.obs[fallback_key].copy()
            log.info(
                "  Fallback to best_resolution: %.1f → leiden set", cfg.clustering.best_resolution
            )
        else:
            avail = [k for k in sub.obs if k.startswith("sub_leiden_")]
            if avail:
                sub.obs["leiden"] = sub.obs[avail[-1]].copy()
                log.info("  Fallback to %s for 'leiden'", avail[-1])
            else:
                log.warning("No Leiden results available — skipping cluster label.")
                sub.obs["leiden"] = "0"

    # ── (i) Save UMAP plots ───────────────────────────────────────────
    fig_dir = os.path.join(cfg.figure_dir, "06_subcluster")
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    sc.settings.autoshow = False
    safe_cell_type = args.cell_type.replace(" ", "_").replace("/", "_")
    safe_filename = f"05_sub_{safe_cell_type}.h5ad"
    output_path = os.path.join(cfg.h5ad_dir, safe_filename)

    # Leiden at best resolution
    safe_plot(
        sc.pl.umap,
        sub,
        color="leiden",
        show=False,
        save=f"sub_{safe_cell_type}_leiden_umap.pdf",
        title=f"{args.cell_type} — leiden",
        cfg=cfg,
    )

    # Multi-resolution comparison
    res_keys = [
        f"sub_leiden_r{r}"
        for r in cfg.clustering.param_grid_resolutions
        if f"sub_leiden_r{r}" in sub.obs
    ]
    i = -1
    if res_keys:
        n_res = len(res_keys)
        n_cols = min(3, n_res)
        n_rows = int(np.ceil(n_res / n_cols))
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(cfg.plot.umap_panel_size[0] * n_cols, cfg.plot.umap_panel_size[1] * n_rows),
        )
        axes = axes.ravel() if n_res > 1 else [axes]
        for i, key in enumerate(res_keys):
            sc.pl.umap(
                sub,
                color=key,
                ax=axes[i],
                show=False,
                legend_fontsize=8,
                title=f"{safe_cell_type} — {key}",
            )
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")
        fig.tight_layout()
        fig.savefig(
            os.path.join(fig_dir, f"sub_{safe_cell_type}_multires.pdf"),
            dpi=cfg.plot.figure_dpi,
            bbox_inches="tight",
        )
        plt.close(fig)
        log.info("Multi-resolution UMAP saved for %s", safe_cell_type)

    # ── Save intermediate results before AI annotation ─────────────────
    safe_write(sub, output_path, cfg=cfg)
    log.info("Intermediate results saved (pre-AI): %s", output_path)

    # ── (j) AI-based subcluster annotation ────────────────────────────
    if cfg.ai.enabled and cfg.ai.subcluster:
        try:
            from core.ai.caller import ai_query
            from core.ai.prompts import build_annotation_prompt

            log.info("Running AI subcluster re-annotation...")

            # build_annotation_prompt runs rank_genes_groups internally
            # and returns (system_prompt, user_prompt)
            sys_prompt, user_prompt = build_annotation_prompt(
                sub,
                tissue=args.cell_type,
                species="unknown",
                precomputed_rank=True,
            )

            result = ai_query(sys_prompt, user_prompt, cfg.ai)

            if not result:
                log.warning("AI returned empty response — falling back to numeric labels")
                sub.obs["sub_ai_label"] = ("Subcluster_" + sub.obs["leiden"].astype(str)).astype(
                    "category"
                )
            else:
                log.info("AI response received (%d chars)", len(result))

                # ── Parse JSON from AI response ──
                # Strip potential markdown code fences
                cleaned = extract_json_block(result)

                parsed = json.loads(cleaned)

                # Build per-cluster labels with cell_type + subtype
                ai_labels = {}
                for cluster_id, info in parsed.items():
                    cell_type = info.get("cell_type", "Unknown")
                    subtype = info.get("subtype", "N/A")
                    if subtype and subtype.upper() != "N/A":
                        label = f"{cell_type}_{subtype}"
                    else:
                        label = cell_type
                    # Sanitize for categorical use
                    label = label.replace(" ", "_").replace("/", "_")
                    ai_labels[cluster_id] = label

                # Convert to string first to avoid Categorical restrictions
                # when adding new values via map()
                sub.obs["leiden"] = sub.obs["leiden"].astype(str)

                # Map to sub.obs as categorical
                sub.obs["sub_ai_label"] = (sub.obs["leiden"].map(ai_labels)).astype("category")
                n_ai_types = sub.obs["sub_ai_label"].nunique()
                log.info("AI annotation: %d subcluster types identified", n_ai_types)

                # Log per-cluster AI mapping
                for cluster_id in sorted(sub.obs["leiden"].unique(), key=lambda x: int(x)):
                    label = ai_labels.get(str(cluster_id), "Unmapped")
                    count = (sub.obs["leiden"] == cluster_id).sum()
                    log.info("  Subcluster %s → %s (%d cells)", cluster_id, label, count)

                # Save AI-annotated UMAP
                safe_plot(
                    sc.pl.umap,
                    sub,
                    color="sub_ai_label",
                    show=False,
                    save=f"sub_{safe_cell_type}_umap_ai.pdf",
                    title=f"{args.cell_type} — AI subcluster",
                    cfg=cfg,
                )

        except Exception as e:
            log.warning("AI subcluster annotation failed: %s", e)
            log.info("Falling back to numeric subcluster labels.")
            sub.obs["sub_ai_label"] = ("Subcluster_" + sub.obs["leiden"].astype(str)).astype(
                "category"
            )
    else:
        log.info(
            "AI subcluster annotation disabled (CFG.ai.enabled=%s, CFG.ai.subcluster=%s)",
            cfg.ai.enabled,
            cfg.ai.subcluster,
        )

    # Restore categorical dtype for consistency
    sub.obs["leiden"] = sub.obs["leiden"].astype("category")

    # ── (k) Save subcluster result ────────────────────────────────────
    safe_write(sub, output_path, cfg=cfg)
    log.info("Saved: %s", output_path)

    # ── (k2) Auto-writeback into 05_annotated.h5ad ───────────────────
    # Merge sub_ai_label back into the main annotation file as the
    # ``cell_subtype`` column, so downstream steps see a single source
    # of truth. No-op when AI annotation was not actually run.
    main_path = os.path.join(cfg.h5ad_dir, "05_annotated.h5ad")
    auto_writeback(
        sub=sub,
        cell_type=args.cell_type,
        main_path=main_path,
        log=log,
        cfg=cfg,
    )

    # ── (l) Summary ───────────────────────────────────────────────────
    n_clusters = sub.obs["leiden"].nunique()
    log.info("=" * 50)
    log.info("Subcluster Summary")
    log.info("  Cell type:       %s", args.cell_type)
    log.info("  Cells:           %d", n_cells)
    log.info("  Subclusters:     %d", n_clusters)
    log.info("  Resolution:      %.1f", cfg.marker.subcluster_resolution)
    log.info("  Per-cluster counts:")
    for cluster_id in sorted(sub.obs["leiden"].unique(), key=lambda x: int(x)):
        count = (sub.obs["leiden"] == cluster_id).sum()
        pct = 100.0 * count / n_cells
        log.info("    Cluster %s: %d cells (%.1f%%)", cluster_id, count, pct)
    log.info("=" * 50)
    log.info("Step 06 done, %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
