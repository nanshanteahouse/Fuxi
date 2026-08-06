#!/usr/bin/env python3
"""
Step 03: GO/KEGG enrichment analysis for bulk DEGs via GSEApy.

Performs:
  - ORA (over-representation analysis) on up/down-regulated DEGs separately
  - Pre-ranked GSEA on all genes ranked by -log10(padj) × sign(logFC)
  - Dot plots for top significant terms

Input:  tables/02_de_significant.csv  (from step 02)
        tables/02_de_results.csv      (for GSEA ranking)
Output: tables/03_enrichment_up.csv
        tables/03_enrichment_down.csv
        tables/03_gsea.csv
        figures/03_enrichment_up_dot.png
        figures/03_enrichment_down_dot.png
        figures/03_gsea_dot.png
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from core.utils import resolve_config, save_figure, setup_logger


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("03_enrichment", os.path.join(cfg.log_dir, "03_enrichment.log"))
    log.info("Step 03: GO/KEGG enrichment for bulk DEGs")

    # ── Check if enrichment is enabled ──
    if not cfg.enrichment.run:
        log.info("Enrichment disabled (enrichment.run=False) — skipping")
        return

    # ── Import gseapy ──
    try:
        import gseapy as gp
    except ImportError:
        log.error("gseapy not installed — run: pip install gseapy")
        sys.exit(1)

    # ── Load significant DEGs ──
    sig_path = os.path.join(cfg.table_dir, "02_de_significant.csv")
    if not os.path.exists(sig_path):
        log.warning("Skip: %s not found. Run step 02 first.", sig_path)
        return

    sig_df = pd.read_csv(sig_path)
    if sig_df.empty:
        log.warning("Skip: %s is empty — no significant DEGs found.", sig_path)
        return
    log.info("Loaded %d significant DEGs from %s", len(sig_df), sig_path)

    gene_col = "gene" if "gene" in sig_df.columns else sig_df.columns[0]
    lfc_col = "log2FoldChange"

    # ── Create output dirs ──
    os.makedirs(cfg.table_dir, exist_ok=True)
    os.makedirs(cfg.figure_dir, exist_ok=True)

    # ── Split genes by direction ──
    if lfc_col not in sig_df.columns:
        log.warning("Column '%s' not found; can't split by direction — skipping ORA.", lfc_col)
        up_genes, down_genes = [], []
    else:
        up_genes = (
            sig_df.loc[sig_df[lfc_col] > 0, gene_col]
            .dropna()
            .astype(str)
            .str.upper()
            .unique()
            .tolist()
        )
        down_genes = (
            sig_df.loc[sig_df[lfc_col] < 0, gene_col]
            .dropna()
            .astype(str)
            .str.upper()
            .unique()
            .tolist()
        )
        log.info(
            "Up-regulated DEGs: %d  |  Down-regulated DEGs: %d", len(up_genes), len(down_genes)
        )

    # ── ORA: up-regulated ──
    if up_genes and len(up_genes) >= cfg.enrichment.min_size:
        log.info("ORA: %d up-regulated genes → %s", len(up_genes), cfg.enrichment.gene_sets)
        try:
            ora_up = gp.enrichr(
                gene_list=up_genes,
                gene_sets=cfg.enrichment.gene_sets,
                organism=cfg.enrichment.organism,
                outdir=None,
            )
            if ora_up.results is not None and not ora_up.results.empty:
                ora_up.results.to_csv(
                    os.path.join(cfg.table_dir, "03_enrichment_up.csv"), index=False
                )
                n = (ora_up.results["Adjusted P-value"] < cfg.enrichment.pval_cutoff).sum()
                log.info("  → %d/%d significant", n, len(ora_up.results))
                _dot_ora(ora_up.results, "up", cfg, log)
        except Exception as e:
            log.warning("Up ORA failed: %s", e)
    else:
        log.info("Skip up ORA (%d genes, need ≥ %d)", len(up_genes), cfg.enrichment.min_size)

    # ── ORA: down-regulated ──
    if down_genes and len(down_genes) >= cfg.enrichment.min_size:
        log.info("ORA: %d down-regulated genes → %s", len(down_genes), cfg.enrichment.gene_sets)
        try:
            ora_down = gp.enrichr(
                gene_list=down_genes,
                gene_sets=cfg.enrichment.gene_sets,
                organism=cfg.enrichment.organism,
                outdir=None,
            )
            if ora_down.results is not None and not ora_down.results.empty:
                ora_down.results.to_csv(
                    os.path.join(cfg.table_dir, "03_enrichment_down.csv"), index=False
                )
                n = (ora_down.results["Adjusted P-value"] < cfg.enrichment.pval_cutoff).sum()
                log.info("  → %d/%d significant", n, len(ora_down.results))
                _dot_ora(ora_down.results, "down", cfg, log)
        except Exception as e:
            log.warning("Down ORA failed: %s", e)
    else:
        log.info("Skip down ORA (%d genes, need ≥ %d)", len(down_genes), cfg.enrichment.min_size)

    # ── Pre-ranked GSEA ──
    all_de_path = os.path.join(cfg.table_dir, "02_de_results.csv")
    if not os.path.exists(all_de_path):
        log.warning("Skip GSEA: %s not found", all_de_path)
    elif "padj" not in sig_df.columns:
        log.warning("Skip GSEA: 'padj' column missing from significant results")
    else:
        log.info("GSEA: preranked on all %d results", len(sig_df))
        all_df = pd.read_csv(all_de_path).dropna(subset=["padj", lfc_col]).copy()
        all_df["padj_clip"] = all_df["padj"].clip(lower=1e-300)
        all_df["rank_metric"] = -np.log10(all_df["padj_clip"]) * np.sign(all_df[lfc_col])
        all_df = all_df.sort_values("rank_metric", ascending=False)

        gene_col_all = "gene" if "gene" in all_df.columns else all_df.columns[0]
        rank_series = all_df.set_index(gene_col_all)["rank_metric"]
        rank_series = rank_series[~rank_series.index.duplicated(keep="first")]

        if len(rank_series) < cfg.enrichment.min_size:
            log.info(
                "Skip GSEA: %d ranked genes < min_size=%d",
                len(rank_series),
                cfg.enrichment.min_size,
            )
        else:
            try:
                gsea_res = gp.prerank(
                    rnk=rank_series,
                    gene_sets=cfg.enrichment.gene_sets,
                    min_size=cfg.enrichment.min_size,
                    max_size=cfg.enrichment.max_size,
                    permutation_num=cfg.enrichment.permutations,
                    outdir=None,
                    seed=cfg.execution.random_seed,
                    verbose=False,
                    no_plot=True,
                )
                if gsea_res.res2d is not None and not gsea_res.res2d.empty:
                    gsea_res.res2d.to_csv(os.path.join(cfg.table_dir, "03_gsea.csv"), index=False)
                    n = (gsea_res.res2d["FDR q-val"] < cfg.enrichment.pval_cutoff).sum()
                    log.info("  → %d/%d significant", n, len(gsea_res.res2d))
                    _dot_gsea(gsea_res.res2d, cfg, log)
            except Exception as e:
                log.warning("GSEA failed: %s", e)

    log.info("Step 03 complete, took %.1fs", time.time() - t0)


# ── Plot helpers ──────────────────────────────────────────────────────────


def _dot_ora(df: pd.DataFrame, direction: str, cfg, log):
    """Horizontal dot plot: x = -log10(padj), color = overlap fraction."""
    sig = df[df["Adjusted P-value"] < cfg.enrichment.pval_cutoff].copy()
    if sig.empty:
        sig = df.head(10)
    if len(sig) < 3:
        log.info("  Skip dot plot (%s): <3 terms", direction)
        return
    sig = sig.sort_values("Adjusted P-value").head(20)
    overlap_parts = sig["Overlap"].astype(str).str.split("/")
    overlap_ratio = overlap_parts.str[0].astype(float) / overlap_parts.str[1].astype(float).clip(
        lower=1
    )

    # Shorten term names: strip trailing '(GO:NNNNNNN)'
    terms = sig["Term"].str.replace(r"\s*\(GO:\d+\)$", "", regex=True).str[:60]

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.3 * len(sig))))
    x = -np.log10(sig["Adjusted P-value"].clip(lower=1e-300))
    ax.scatter(
        x,
        range(len(sig)),
        s=overlap_ratio * 200 + 20,
        c=overlap_ratio,
        cmap=cfg.plot.palette.dotplot_fill,
        edgecolors=cfg.plot.palette.significance_edge,
        linewidths=0.5,
    )
    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(terms, fontsize=8)
    ax.set_xlabel("-log10(Adjusted P-value)")
    ax.set_title(f"Enrichment ({direction}-regulated), top {len(sig)} terms")
    fig.tight_layout()
    path = os.path.join(cfg.figure_dir, f"03_enrichment_{direction}_dot")
    save_figure(fig, path, cfg=cfg, bbox_inches="tight")
    plt.close(fig)
    log.info("  Dot plot: %s", path)


def _dot_gsea(df: pd.DataFrame, cfg, log):
    """Horizontal dot plot: x = -log10(FDR), color = NES."""
    sig = df[df["FDR q-val"] < cfg.enrichment.pval_cutoff].copy()
    if sig.empty:
        sig = df.head(10)
    if len(sig) < 3:
        log.info("  Skip GSEA dot plot: <3 terms")
        return
    sig = sig.sort_values("FDR q-val").head(20)
    terms = sig["Term"].str[:60]
    log10_fdr = -np.log10(sig["FDR q-val"].clip(lower=1e-300))
    vmax = max(abs(sig["NES"].min()), abs(sig["NES"].max()))

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.3 * len(sig))))
    sc = ax.scatter(
        log10_fdr,
        range(len(sig)),
        s=log10_fdr * 20,
        c=sig["NES"],
        cmap=cfg.plot.palette.heatmap,
        norm=plt.Normalize(-vmax, vmax),
        edgecolors=cfg.plot.palette.significance_edge,
        linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label="NES")
    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(terms, fontsize=8)
    ax.set_xlabel("-log10(FDR q-val)")
    ax.set_title(f"GSEA (preranked), top {len(sig)} terms")
    fig.tight_layout()
    path = os.path.join(cfg.figure_dir, "03_gsea_dot")
    save_figure(fig, path, cfg=cfg, bbox_inches="tight")
    plt.close(fig)
    log.info("  GSEA dot plot: %s", path)


if __name__ == "__main__":
    main()
