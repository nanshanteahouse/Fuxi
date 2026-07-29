#!/usr/bin/env python3
"""
Step 08: 轨迹分析 — PAGA/DPT 或 scVelo + CellRank
===================================================
继承深度轨迹分析最佳实践:
  1. PAGA (在子聚类或聚类级别上)
  2. 根细胞自动识别 (ROI 类型或标记基因)
  3. 扩散伪时间 (DPT)
  4. 分支间差异表达
  5. 基因沿伪时间表达趋势

输入: 04_clustered.h5ad (需要 Stage 05 注释结果)
输出: 05_final.h5ad (含 PAGA, DPT, 分支结果) + tables + figures
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from core.kb import load_kb
from core.utils import resolve_config, safe_plot, safe_write, setup_logger


def recompute_neighbors(adata, cfg, log):
    """确保邻居图和 UMAP 存在"""
    if "neighbors" not in adata.uns:
        log.info("Recomputing neighbors...")
        use_rep = "X_integrated" if "X_integrated" in adata.obsm else "X_pca"
        sc.pp.neighbors(
            adata,
            n_pcs=min(cfg.pca.n_pcs_use, adata.obsm[use_rep].shape[1]),
            n_neighbors=cfg.clustering.n_neighbors,
            use_rep=use_rep,
            random_state=cfg.execution.random_seed,
        )
    if "X_umap" not in adata.obsm:
        log.info("Recomputing UMAP...")
        sc.tl.umap(adata, random_state=cfg.execution.random_seed)
    log.info("Neighbors + UMAP ready")


def run_paga(adata, cfg, log):
    """PAGA 轨迹拓扑"""
    group_col = (
        "cell_type_sub"
        if "cell_type_sub" in adata.obs
        else ("cell_type" if "cell_type" in adata.obs else "leiden")
    )
    log.info("PAGA (groupby=%s)...", group_col)
    sc.tl.paga(adata, groups=group_col)
    n_edges = np.sum(adata.uns["paga"]["connectivities"].data > 0)
    log.info("  PAGA edges: %d", n_edges)
    safe_plot(
        sc.pl.paga,
        adata,
        color=group_col,
        show=False,
        save="paga_graph.pdf",
        title="PAGA trajectory",
        cfg=cfg,
    )
    safe_plot(
        sc.pl.paga_compare,
        adata,
        basis="umap",
        color=group_col,
        show=False,
        save="paga_umap.pdf",
        edge_width_scale=0.5,
        title="PAGA on UMAP",
        cfg=cfg,
    )


def find_root_cells(adata, cfg, log):
    """自动识别根细胞"""
    # 方法 1: 指定根细胞类型
    if cfg.trajectory.root_cell_types:
        log.info("Root cells: type %s + earliest stage", cfg.trajectory.root_cell_types)
        if "stage" in adata.obs and cfg.sample_meta.stage_order:
            root_mask = adata.obs["cell_type"].isin(cfg.trajectory.root_cell_types) & (
                adata.obs["stage"] == cfg.sample_meta.stage_order[0]
            )
        else:
            root_mask = adata.obs["cell_type"].isin(cfg.trajectory.root_cell_types)
        if root_mask.sum() > 0:
            log.info("  Root cells: %d", root_mask.sum())
            return root_mask.values
        else:
            log.warning("  Root cells of specified type not found, trying marker gene method")

    # 方法 2: 标记基因自动检测
    if cfg.trajectory.root_markers:
        log.info("Root cells: marker gene method %s", cfg.trajectory.root_markers)
        markers_present = [g for g in cfg.trajectory.root_markers if g in adata.raw.var_names]
        if markers_present:
            group_col = "cell_type" if "cell_type" in adata.obs else "leiden"
            cluster_scores = []
            for cl in adata.obs[group_col].cat.categories:
                mask = adata.obs[group_col] == cl
                sub = adata.raw[mask]
                gene_indices = [list(adata.raw.var_names).index(g) for g in markers_present]
                gene_exprs = sub.X[:, gene_indices]
                if issparse(gene_exprs):
                    scores = gene_exprs.mean(axis=0).A1.tolist()
                else:
                    scores = gene_exprs.mean(axis=0).tolist()
                cluster_scores.append((cl, np.mean(scores)))
            cluster_scores.sort(key=lambda x: -x[1])
            best_cl = cluster_scores[0][0]
            root_mask = adata.obs[group_col] == best_cl
            log.info("  High-score cluster: %s (score=%.4f)", best_cl, cluster_scores[0][1])
            log.info("  Root cells: %d", root_mask.sum())
            return root_mask.values

    # 方法 3: 回退到最早阶段的细胞
    log.warning("  Cannot auto-determine root, using earliest stage cells.")
    if "stage" in adata.obs and cfg.sample_meta.stage_order:
        root_mask = adata.obs["stage"] == cfg.sample_meta.stage_order[0]
        log.info(
            "  Root cells: %d (earliest stage %s)", root_mask.sum(), cfg.sample_meta.stage_order[0]
        )
        return root_mask.values
    # 最终回退: 第一个细胞
    log.warning("  Final fallback: using first cell as root.")
    root_mask = np.zeros(adata.n_obs, dtype=bool)
    root_mask[0] = True
    return root_mask


def compute_dpt(adata, root_mask, cfg, log):
    """扩散图 + 扩散伪时间"""
    log.info("Diffusion map (n_comps=%d)...", cfg.trajectory.n_diffmap_comps)
    sc.tl.diffmap(adata, n_comps=cfg.trajectory.n_diffmap_comps)

    log.info("Diffusion pseudotime...")
    adata.uns["iroot"] = np.flatnonzero(root_mask)[0]
    for nb in [cfg.trajectory.n_branchings, 1, 0]:
        try:
            sc.tl.dpt(adata, n_branchings=nb)
            log.info("  DPT complete (n_branchings=%d)", nb)
            break
        except ValueError:
            log.warning(
                "DPT n_branchings=%d failed, trying n_branchings=%d", nb, 1 if nb > 1 else 0
            )
            continue
    log.info(
        "  DPT range: %.3f – %.3f",
        adata.obs["dpt_pseudotime"].min(),
        adata.obs["dpt_pseudotime"].max(),
    )

    safe_plot(
        sc.pl.umap,
        adata,
        color="dpt_pseudotime",
        show=False,
        save="pseudotime_umap.pdf",
        cmap=cfg.plot.palette.pseudotime,
        cfg=cfg,
    )
    safe_plot(
        sc.pl.diffmap,
        adata,
        color="dpt_pseudotime",
        show=False,
        save="pseudotime_diffmap.pdf",
        cmap=cfg.plot.palette.pseudotime,
        cfg=cfg,
    )


def branch_analysis(adata, cfg, log, table_dir) -> Optional[pd.DataFrame]:
    """分支间差异表达 (分支间配对比较策略)"""
    if "cell_type" not in adata.obs:
        log.info("No cell_type annotation, skipping branch analysis.")
        return

    if hasattr(cfg, "trajectory_branches") and cfg.trajectory_branches:
        branches = cfg.trajectory_branches
    else:
        # Auto-detect: use cell type pairs from PAGA graph
        if "cell_type" in adata.obs:
            avail_types = list(adata.obs["cell_type"].cat.categories)
            branches = []
            for i in range(len(avail_types) - 1):
                branches.append((avail_types[i], avail_types[i + 1]))
        else:
            branches = []
    # 仅保留数据中存在的分支
    avail_types = set(adata.obs["cell_type"].cat.categories)
    branches = [(p, c) for p, c in branches if p in avail_types and c in avail_types]

    if not branches:
        log.info("No matching branch pairs, skipping.")
        return

    log.info("Branch differential expression analysis...")
    branch_results = []
    for parent, child in branches:
        mask = adata.obs["cell_type"].isin([parent, child])
        sub = adata[mask].copy()
        if sub.obs["cell_type"].value_counts().min() < 10:
            log.info("  %s → %s: insufficient cells", parent, child)
            continue
        try:
            sc.tl.rank_genes_groups(
                sub,
                groupby="cell_type",
                groups=[child],
                reference=parent,
                method=cfg.de.branch_method,
                n_genes=50,
                use_raw=True,
                random_state=cfg.execution.random_seed,
            )
            de_df = sc.get.rank_genes_groups_df(sub, group=child)
            if cfg.de.pval_cutoff is not None:
                de_df = de_df[de_df["pvals_adj"] < cfg.de.pval_cutoff].copy()
            de_df["branch"] = f"{child}_vs_{parent}"
            branch_results.append(de_df)
            n_up = (de_df["logfoldchanges"] > 0).sum()
            n_down = (de_df["logfoldchanges"] < 0).sum()
            log.info(
                "  %s → %s: %d DEGs (%d up, %d down)", parent, child, len(de_df), n_up, n_down
            )
        except (KeyError, ValueError) as e:
            log.warning("  %s → %s branch DE failed: %s", parent, child, e)

    if not branch_results:
        log.warning(
            "branch_deg.csv NOT generated: all branch DE pairs failed or were filtered. "
            "Check cfg.de.branch_method and pval_cutoff."
        )
        return None
    if branch_results:
        combined = pd.concat(branch_results, ignore_index=True)
        out_path = os.path.join(table_dir, "branch_deg.csv")
        combined.to_csv(out_path, index=False)
        log.info("  Branch DEG exported: %s (%d rows)", out_path, len(combined))
        return combined

    return None


def _select_pseudotime_correlated(adata, cfg) -> Tuple[List[str], pd.DataFrame]:
    """Select top genes correlated with dpt_pseudotime via Spearman correlation.

    Pre-filters to expressed genes and top HVGs, applies BH correction,
    balances positive and negative correlations.
    """
    if adata.raw is None:
        return [], pd.DataFrame(columns=("gene", "rho", "pval_raw", "pval_adj"))  # type: ignore[arg-type]

    # 1. Pre-filter: expressed in >=1% of cells
    if issparse(adata.raw.X):
        expr_frac = np.array((adata.raw.X > 0).mean(axis=0)).ravel()
    else:
        expr_frac = (adata.raw.X > 0).mean(axis=0)

    expressed_mask = expr_frac >= 0.01

    # Restrict to HVGs if available, else top 3000 by mean expression
    if "highly_variable" in adata.raw.var:
        hvgs = adata.raw.var["highly_variable"].values
        candidate_mask = expressed_mask & hvgs
    else:
        expr_mean = np.array(adata.raw.X.mean(axis=0)).ravel()
        top_n = min(3000, int(expressed_mask.sum()))
        if top_n == 0:
            return [], pd.DataFrame(columns=("gene", "rho", "pval_raw", "pval_adj"))  # type: ignore[arg-type]
        expr_mean_sorted_idx = np.argsort(-expr_mean)
        top_idx = set(expr_mean_sorted_idx[: top_n * 3])
        candidate_mask = np.array([i in top_idx for i in range(adata.raw.n_vars)]) & expressed_mask

    candidate_indices = np.where(candidate_mask)[0]
    if len(candidate_indices) == 0:
        return [], pd.DataFrame(columns=("gene", "rho", "pval_raw", "pval_adj"))  # type: ignore[arg-type]

    # 2. Extract pseudotime (drop NaN)
    pseudotime = adata.obs["dpt_pseudotime"].values
    pt_mask = ~np.isnan(pseudotime)
    if pt_mask.sum() < 2:
        return [], pd.DataFrame(columns=("gene", "rho", "pval_raw", "pval_adj"))  # type: ignore[arg-type]
    pseudotime_clean = pseudotime[pt_mask]

    # 3. Compute Spearman correlation per candidate gene
    rhos: List[float] = []
    pvals: List[float] = []
    gene_names: List[str] = []
    for idx in candidate_indices:
        # Extract expression values as float array
        if issparse(adata.raw.X):
            raw_values: np.ndarray = adata.raw.X[:, idx].toarray().ravel()
        else:
            raw_values = np.asarray(adata.raw.X[:, idx], dtype=float)
        expr = np.asarray(raw_values, dtype=float)[pt_mask]

        if np.var(expr) == 0:
            continue

        result = spearmanr(expr, pseudotime_clean)
        rho: float = result.statistic  # type: ignore[assignment]
        pval: float = result.pvalue  # type: ignore[assignment]
        if np.isnan(rho) or np.isnan(pval):
            continue

        rhos.append(rho)
        pvals.append(pval)
        gene_names.append(adata.raw.var_names[idx])

    if len(rhos) == 0:
        return [], pd.DataFrame(columns=("gene", "rho", "pval_raw", "pval_adj"))  # type: ignore[arg-type]

    # 4. BH correction
    _pvals_adj = multipletests(pvals, method="fdr_bh")  # type: ignore[var-annotated]
    _, pvals_adj, _, _ = _pvals_adj

    full_df = pd.DataFrame(
        {
            "gene": gene_names,
            "rho": rhos,
            "pval_raw": pvals,
            "pval_adj": pvals_adj if pvals_adj is not None else [1.0] * len(gene_names),
        }
    )

    # 5. Filter by adjusted p-value and correlation strength
    selected = []
    for i in range(len(gene_names)):
        if (
            pvals_adj is not None
            and pvals_adj[i] < cfg.trajectory.pseudotime_cor_pval
            and abs(rhos[i]) > 0.2
        ):
            selected.append((gene_names[i], rhos[i]))

    # 6. Sort by |rho| descending
    selected.sort(key=lambda x: abs(x[1]), reverse=True)

    # 7. Balanced split: up to half from each sign
    n = cfg.trajectory.pseudotime_n_correlated
    half_n = n // 2
    pos = [(g, r) for g, r in selected if r > 0]
    neg = [(g, r) for g, r in selected if r < 0]

    result = []
    result.extend(g for g, _ in pos[:half_n])
    result.extend(g for g, _ in neg[:half_n])
    return result, full_df


def gene_trends(adata, cfg, log, table_dir, branch_results: Optional[pd.DataFrame] = None):
    """基因表达沿伪时间趋势——四源数据驱动选择"""
    # Guard A: DPT exists and has variance
    if "dpt_pseudotime" not in adata.obs:
        log.info("No DPT, skipping gene trends.")
        return
    if adata.obs["dpt_pseudotime"].dropna().nunique() < 2:
        log.info("DPT has insufficient variance, skipping gene trends.")
        return

    # Guard B: raw data exists
    if adata.raw is None:
        log.info("No raw data, skipping gene trends.")
        return

    selected_genes = []
    source_counts = {}
    branch_top: List[str] = []

    # Source 1: Branch DE
    if branch_results is not None and not branch_results.empty:
        try:
            if all(c in branch_results.columns for c in ["names", "scores", "pvals_adj"]):
                branch_de = branch_results.sort_values("scores", ascending=False)
                branch_top = (
                    branch_de["names"]
                    .drop_duplicates()
                    .head(cfg.trajectory.pseudotime_n_branch_de)
                    .tolist()
                )
                branch_top = [g for g in branch_top if g in adata.raw.var_names]
                selected_genes.extend(branch_top)
                source_counts["branch_DE"] = len(branch_top)
            else:
                log.warning("branch_results has unexpected columns, skipping branch DE source")
        except Exception as e:
            log.warning("Failed to extract branch DE genes: %s", e)

    corr_genes: List[str] = []
    corr_df: pd.DataFrame = pd.DataFrame(columns=("gene", "rho", "pval_raw", "pval_adj"))  # type: ignore[arg-type]
    # Source 2: Pseudotime correlation
    try:
        corr_genes, corr_df = _select_pseudotime_correlated(adata, cfg)
        selected_genes.extend(corr_genes)
        source_counts["pseudotime_correlation"] = len(corr_genes)
    except Exception as e:
        log.warning("Failed to compute pseudotime correlation: %s", e)

    # Source 3: CFG override
    if cfg.trajectory.pseudotime_genes:
        override = [g for g in cfg.trajectory.pseudotime_genes if g in adata.raw.var_names]
        excluded = set(cfg.trajectory.pseudotime_genes) - set(override)
        if excluded:
            log.warning(
                "CFG.trajectory.pseudotime_genes excluded (not in var_names): %s", excluded
            )
        selected_genes.extend(override)
        source_counts["CFG_override"] = len(override)

    kb_genes = []
    # Source 4: KB markers
    if cfg.tissue_kb:
        try:
            kb = load_kb(cfg.tissue_kb)
            kb_genes = set()
            for cell_type, info in kb.items():
                if cell_type in ("expert_rules", "_meta"):
                    continue
                if "markers" in info:
                    for cat in ("confirm", "add"):
                        if cat in info["markers"]:
                            kb_genes.update(info["markers"][cat].keys())
            kb_genes = [g for g in kb_genes if g in adata.raw.var_names]
            selected_genes.extend(kb_genes)
            source_counts["KB_markers"] = len(kb_genes)
        except ValueError as e:
            log.warning("Unsupported tissue_kb '%s': %s. Skipping KB markers.", cfg.tissue_kb, e)
        except Exception as e:
            log.warning("Failed to load KB markers: %s", e)

    # Union: deduplicate preserving first occurrence (priority order 1->2->3->4)
    seen = set()
    union_genes = []
    for g in selected_genes:
        if g not in seen:
            seen.add(g)
            union_genes.append(g)

    # Cap total
    max_genes = cfg.trajectory.pseudotime_n_correlated * 2
    union_genes = union_genes[:max_genes]

    # Export Spearman correlation full results
    if not corr_df.empty:
        corr_csv = os.path.join(table_dir, "pseudotime_trend_genes.csv")
        corr_df.to_csv(corr_csv, index=False)
        log.info("  Pseudotime trend genes exported: %s (%d rows)", corr_csv, len(corr_df))

    # Export selected union with source annotation
    source_map = {}
    for g in union_genes:
        if g in branch_top:
            source_map.setdefault(g, set()).add("branch_DE")
        if g in corr_genes:
            source_map.setdefault(g, set()).add("pseudotime_correlation")
        if g in cfg.trajectory.pseudotime_genes:
            source_map.setdefault(g, set()).add("CFG_override")
        if g in kb_genes:
            source_map.setdefault(g, set()).add("KB_markers")
    sel_df = pd.DataFrame(
        {
            "gene": list(source_map.keys()),
            "source": ["+".join(sorted(v)) for v in source_map.values()],
        },
    )
    sel_csv = os.path.join(table_dir, "pseudotime_trend_genes_selected.csv")
    sel_df.to_csv(sel_csv, index=False)
    log.info("  Selected pseudotime trend genes exported: %s (%d rows)", sel_csv, len(sel_df))
    log.info(
        "Gene trends along pseudotime (%d unique genes from %s)...",
        len(union_genes),
        source_counts,
    )

    if not union_genes:
        log.warning("No genes selected for pseudotime trends from any source.")
        return

    # Scatter plots: first 6 genes
    for gene in union_genes[:6]:
        safe_plot(
            sc.pl.scatter,
            adata,
            x="dpt_pseudotime",
            y=gene,
            use_raw=True,
            show=False,
            save=f"trend_{gene}.pdf",
            cfg=cfg,
        )

    # Heatmap: if >=5 genes, with binned pseudotime
    if len(union_genes) >= 5:
        n_bins = min(10, int(adata.obs["dpt_pseudotime"].dropna().nunique() - 1))
        if n_bins >= 2:
            try:
                adata_sub = adata[adata.obs["dpt_pseudotime"].notna()].copy()
                adata_sub.obs["dpt_pseudotime_bin"] = pd.qcut(
                    adata_sub.obs["dpt_pseudotime"], q=n_bins, duplicates="drop"
                ).astype(str)  # type: ignore[union-attr]
                safe_plot(
                    sc.pl.heatmap,
                    adata_sub,
                    var_names=union_genes,
                    groupby="dpt_pseudotime_bin",
                    use_raw=True,
                    show=False,
                    save="dev_gene_heatmap.pdf",
                    cfg=cfg,
                )
            except (ValueError, TypeError) as e:
                log.warning("Gene trend heatmap skipped: qcut failed (%s)", e)
        else:
            log.info("Not enough unique pseudotime values (%d) for heatmap binning.", n_bins)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    cfg = resolve_config(args.config)
    log = setup_logger("08_trajectory", os.path.join(cfg.log_dir, "08_trajectory.log"))
    log.info("Step 08: Trajectory analysis")

    input_path = cfg.annotated_h5ad if os.path.exists(cfg.annotated_h5ad) else cfg.cluster_h5ad
    adata = sc.read(input_path)
    log.info("Loaded: %s — %d cells", input_path, adata.n_obs)

    # 设置图输出目录（必须在 plot 调用之前，否则 scanpy save= 默认写到 ./figures/）
    sc.settings.figdir = os.path.join(cfg.figure_dir, "08_trajectory")
    os.makedirs(sc.settings.figdir, exist_ok=True)

    # ── Table output subdirectory ──────────────────────────
    table_dir = os.path.join(cfg.table_dir, "08_trajectory")
    os.makedirs(table_dir, exist_ok=True)

    # 当 marker_validation PASS 率极低时，退回到 leiden 聚类
    if "marker_validation" in adata.obs and adata.n_obs > 0:
        pass_cells = (adata.obs["marker_validation"] == "PASS").sum()
        pass_rate = pass_cells / adata.n_obs
        pass_rate_min = getattr(cfg.marker, "validation_pass_rate_min", 0.1)
        if pass_rate < pass_rate_min:
            if getattr(cfg, "interactive", False):
                print(f"\n⚠  Annotation validation PASS rate = {pass_rate * 100:.1f}%")
                try:
                    choice = (
                        input(
                            "Trajectory analysis options:\n"
                            "  [l] Use leiden clusters (safe fallback)\n"
                            "  [c] Use cell_type labels anyway\n"
                            "Choice> "
                        )
                        .strip()
                        .lower()
                    )
                except (EOFError, KeyboardInterrupt):
                    choice = "l"
                if choice == "c":
                    log.warning(
                        "User chose cell_type labels despite %.1f%% PASS rate",
                        pass_rate * 100,
                    )
                else:
                    adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
                    log.info("Falling back to leiden clusters for trajectory")
            else:
                log.warning(
                    "marker_validation PASS rate %.1f%% (<%.0f%%) — "
                    "cell_type labels are unreliable, falling back to "
                    "leiden clusters",
                    pass_rate * 100,
                    pass_rate_min * 100,
                )
                adata.obs["cell_type"] = adata.obs["leiden"].astype(str)

    recompute_neighbors(adata, cfg, log)

    # ── scVelo RNA velocity dispatch ──────────────────────────
    if cfg.trajectory.method == "scvelo_cellrank":
        from core.utils._optional import require_scvelo

        require_scvelo()
        import scvelo as scv

        if "spliced" not in adata.layers or "unspliced" not in adata.layers:
            log.error(
                "scVelo requires spliced/unspliced layers in adata. "
                "Did you preprocess with velocyto (e.g., velocyto run10x / run.py) "
                "and load the resulting loom file? "
                "See scVelo documentation for loom preparation."
            )
            sys.exit(1)

        log.info("scVelo RNA velocity (mode=%s)...", cfg.trajectory.scvelo.mode)
        scv.tl.velocity(adata, mode=cfg.trajectory.scvelo.mode)
        scv.tl.velocity_graph(adata)
        log.info("scVelo velocity + velocity_graph computed")

        safe_write(adata, cfg.final_h5ad, cfg=cfg)
        log.info("Step 08 complete (scVelo), took %.1fs", time.time() - t0)
        return
    run_paga(adata, cfg, log)
    root_mask = find_root_cells(adata, cfg, log)
    compute_dpt(adata, root_mask, cfg, log)
    branch_results = branch_analysis(adata, cfg, log, table_dir)
    gene_trends(adata, cfg, log, table_dir, branch_results=branch_results)

    # 最终可视化 (figdir 已在上面设置)
    for color in ["stage", "cell_type", "cell_type_sub", "dpt_pseudotime"]:
        if color in adata.obs or color in adata.obsm:
            safe_plot(
                sc.pl.umap, adata, color=color, show=False, save=f"final_umap_{color}.pdf", cfg=cfg
            )

    safe_write(adata, cfg.final_h5ad, cfg=cfg)
    log.info("Step 08 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
