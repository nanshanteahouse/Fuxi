#!/usr/bin/env python3
"""
Pseudobulk differential expression via decoupler + PyDESeq2.
"""

import os

import numpy as np
import pandas as pd


def run_pseudobulk_de(adata, cfg, log):
    """Run pseudobulk DE analysis using decoupler + PyDESeq2.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts.
    cfg : Config
        Resolved Fuxi config object. Reads cfg.de.pseudobulk.* fields.
    log : logger
        Logger instance.
    """
    pb = cfg.de.pseudobulk

    # ── (a) Input validation ──────────────────────────────────────────

    for col_key in ("celltype_col", "sample_col", "contrast_column"):
        col = getattr(pb, col_key)
        if col not in adata.obs:
            log.error(
                "Column '%s' (de.pseudobulk.%s) not found in adata.obs",
                col,
                col_key,
            )
            return

    if not pb.contrast_treatment or not pb.contrast_baseline:
        log.error(
            "de.pseudobulk.contrast_treatment and contrast_baseline must be non-empty",
        )
        return

    valid_levels = set(adata.obs[pb.contrast_column].unique())
    for val_name, val in [
        ("contrast_treatment", pb.contrast_treatment),
        ("contrast_baseline", pb.contrast_baseline),
    ]:
        if val not in valid_levels:
            log.error(
                "de.pseudobulk.%s='%s' not found in adata.obs['%s']. Valid values: %s",
                val_name,
                val,
                pb.contrast_column,
                sorted(valid_levels),
            )
            return

    # ── (c) Extract raw counts ────────────────────────────────────────

    if adata.raw is not None:
        count_adata = adata.raw.to_adata()
        log.info(
            "Using adata.raw for counts (%d cells, %d genes)",
            count_adata.n_obs,
            count_adata.n_vars,
        )
    else:
        count_adata = adata
        log.warning(
            "adata.raw not found — using .X which may be "
            "log-transformed; DESeq2 requires integer counts",
        )

    # ── (d) Pseudobulk aggregation ────────────────────────────────────

    import decoupler as dc

    log.info(
        "Pseudobulk aggregation: sample_col=%s, groups_col=%s",
        pb.sample_col,
        pb.celltype_col,
    )
    pdata = dc.pp.pseudobulk(
        count_adata,
        sample_col=pb.sample_col,
        groups_col=pb.celltype_col,
        layer=None,
        mode="sum",
    )
    log.info(
        "Pseudobulk output: %d samples × %d genes",
        pdata.n_obs,
        pdata.n_vars,
    )

    # ── (e) Filter by expression ──────────────────────────────────────

    log.info("Filter by expression: group=%s", pb.contrast_column)
    dc.pp.filter_by_expr(
        pdata,
        group=pb.contrast_column,
        min_count=10,
        min_total_count=15,
    )
    log.info(
        "After filtering: %d samples, %d genes",
        pdata.n_obs,
        pdata.n_vars,
    )

    # ── (f) Per cell type QC ──────────────────────────────────────────

    if pb.celltype_col in pdata.obs:
        cell_types = pdata.obs[pb.celltype_col].cat.categories
    else:
        log.error(
            "celltype_col '%s' not in pseudobulk output",
            pb.celltype_col,
        )
        return

    output_subdir = pb.output_dir if pb.output_dir else "pseudobulk_de"
    out_dir = os.path.join(cfg.table_dir, output_subdir)
    os.makedirs(out_dir, exist_ok=True)

    fig_dir = os.path.join(cfg.figure_dir, "07_markers")
    os.makedirs(fig_dir, exist_ok=True)

    all_results = []

    for ct in cell_types:
        ct_mask = pdata.obs[pb.celltype_col] == ct
        ct_condition_counts = pdata.obs.loc[ct_mask, pb.contrast_column].value_counts()

        if any(ct_condition_counts < pb.min_cells_per_group):
            log.warning(
                "Cell type '%s': insufficient samples per condition (%s), skipping",
                ct,
                dict(ct_condition_counts),
            )
            continue

        ct_pdata = pdata[ct_mask].copy()
        if ct_pdata.n_vars == 0:
            log.warning(
                "Cell type '%s': no genes remaining after filtering, skipping",
                ct,
            )
            continue

        # ── (g) DESeq2 per cell type ──────────────────────────────────

        log.info(
            "DESeq2 for cell type '%s' (%d samples)",
            ct,
            ct_pdata.n_obs,
        )
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        dds = DeseqDataSet(
            adata=ct_pdata,
            design=pb.design,
            n_cpus=pb.n_jobs if pb.n_jobs > 0 else None,
        )
        dds.deseq2()

        stat = DeseqStats(
            dds,
            contrast=[
                pb.contrast_column,
                pb.contrast_treatment,
                pb.contrast_baseline,
            ],
        )
        stat.summary()

        if pb.lfc_shrink:
            # Determine the LFC coefficient name for shrinkage
            lfc_coeffs = [
                c
                for c in stat.LFC.columns
                if c != "Intercept" and c.startswith(pb.contrast_column)
            ]
            coeff = lfc_coeffs[0] if lfc_coeffs else None
            if coeff:
                stat.lfc_shrink(coeff=coeff)
            else:
                log.warning("Could not determine LFC coefficient for lfc_shrink, skipping")

        results_df = stat.results_df.copy()
        results_df = results_df.reset_index().rename(
            columns={"index": "gene"},
        )

        out_path = os.path.join(out_dir, f"{ct}_de.csv")
        results_df.to_csv(out_path, index=False)
        log.info("Exported: %s (%d genes)", out_path, len(results_df))
        all_results.append(results_df)

        # ── (h) Volcano plot ──────────────────────────────────────────

        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            plot_df = results_df.dropna(subset=["padj"]).copy()
            plot_df = plot_df[plot_df["padj"] > 0]

            if len(plot_df) == 0:
                log.warning(
                    "No valid padj values for cell type '%s', skipping volcano",
                    ct,
                )
                plt.close()
                continue

            plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"])

            fig, ax = plt.subplots(figsize=(8, 6))

            sig = plot_df["padj"] < pb.alpha
            ax.scatter(
                plot_df.loc[~sig, "log2FoldChange"],
                plot_df.loc[~sig, "neg_log10_padj"],
                c="gray",
                alpha=0.5,
                s=10,
                label="NS",
            )
            ax.scatter(
                plot_df.loc[sig, "log2FoldChange"],
                plot_df.loc[sig, "neg_log10_padj"],
                c="red",
                alpha=0.5,
                s=10,
                label=f"padj<{pb.alpha}",
            )

            ax.axhline(
                -np.log10(pb.alpha),
                color="blue",
                linestyle="--",
                linewidth=0.8,
            )
            ax.axvline(-1, color="gray", linestyle="--", linewidth=0.5)
            ax.axvline(1, color="gray", linestyle="--", linewidth=0.5)

            ax.set_xlabel("log2 Fold Change")
            ax.set_ylabel("-log10(padj)")
            ax.set_title(
                f"{ct} ({pb.contrast_treatment} vs {pb.contrast_baseline})",
            )
            ax.legend(loc="upper right")

            fig_path = os.path.join(
                fig_dir,
                f"pseudobulk_volcano_{ct}.png",
            )
            fig.savefig(fig_path, bbox_inches="tight", dpi=150)
            plt.close(fig)
            log.info("Volcano saved: %s", fig_path)
        except Exception as e:
            log.warning("Volcano plot failed for '%s': %s", ct, e)

    # ── Export combined results ───────────────────────────────────────

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_path = os.path.join(out_dir, "combined_de_results.csv")
        combined.to_csv(combined_path, index=False)
        log.info(
            "Combined results exported: %s (%d rows)",
            combined_path,
            len(combined),
        )
    else:
        log.warning(
            "No pseudobulk DE results produced — check sample/contrast configuration",
        )
