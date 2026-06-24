#!/usr/bin/env python3
"""
Step 09: GO/KEGG 富集分析
=============================
输入: Step 07 输出的 marker_genes_per_group.csv
      （每类细胞 vs 其他所有细胞的 Wilcoxon 标记基因）

方法 (通过 GSEApy + Enrichr API):
  ORA: 取每类上调基因 top N → 过表达分析
  Pre-ranked GSEA: 使用全部基因的 score 排序 → 无需 cutoff

输出:
  {table_dir}/enrichment/
    ora_{gene_set}_{cluster}.csv          — ORA 结果表
    prerank_{gene_set}_{cluster}.csv      — GSEA 结果表
    ora_{gene_set}_summary.csv            — 汇总（所有聚类合并）
  {figure_dir}/enrichment/
    ora_{gene_set}_bubble.pdf             — 气泡图
    ora_{gene_set}_dotplot.pdf            — 点图

依赖: pip install gseapy (需要 Rust 编译器)
      curl https://sh.rustup.rs -sSf | sh  # 如需要
"""
import json
import logging
import sys, os, time, argparse, warnings
from typing import Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config
import numpy as np
import pandas as pd
import scanpy as sc  # for loading annotated h5ad (quality check)

warnings.filterwarnings("ignore", category=FutureWarning)


def enrichr_with_retry(gene_list, gene_sets_library, max_retries=3, log=None):
    """Call gseapy.enrichr() with exponential backoff on HTTP 429.

    The Enrichr API returns 429 (Too Many Requests) when the pipeline
    floods it with concurrent calls (e.g. ThreadPoolExecutor ORA mode).
    On 429 we sleep 5s, 10s, 20s before retrying; any other exception
    is re-raised immediately so genuine API/parse errors aren't masked.

    Args:
        gene_list: Iterable of gene symbols to submit to Enrichr.
        gene_sets_library: Enrichr library name (e.g. 'KEGG_2021_Human').
        max_retries: Maximum number of retries after the initial call.
        log: Optional logger for WARNING-level retry messages; falls back
            to a module-level logger when not provided.

    Returns:
        The Enrichr result object returned by a successful gp.enrichr()
        call (carries .results, .res2d, etc. as appropriate).
    """
    import gseapy as gp

    if log is None:
        log = logging.getLogger(__name__)

    backoff_schedule = [5, 10, 20]
    last_err = None
    for attempt in range(max_retries + 1):
        try:
            return gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_sets_library,
                outdir=None,
                no_plot=True,
                verbose=False,
            )
        except Exception as e:
            last_err = e
            err_text = str(e)
            is_429 = (
                "429" in err_text
                or "Too Many Requests" in err_text
                or (hasattr(e, "response") and getattr(e.response, "status_code", None) == 429)
            )
            if not is_429:
                # Non-rate-limit errors are not retried — surface immediately
                # so genuine API/parse failures aren't masked as timeouts.
                raise
            if attempt >= max_retries:
                log.error(
                    "Enrichr 429 persistent failure (after %d retries): %s", max_retries, err_text
                )
                raise
            wait_s = backoff_schedule[min(attempt, len(backoff_schedule) - 1)]
            log.warning(
                "Enrichr returned 429 (attempt %d/%d), retrying in %.0fs: %s",
                attempt + 1, max_retries + 1, wait_s, err_text,
            )
            time.sleep(wait_s)
    raise last_err  # pragma: no cover


def _extract_gene_symbol(name_val):
    """Extract pure gene symbol from mixed 'ENSG...;ENSG...;\"SYMBOL\"' format.
    
    Handles standard formats:
    - "SYMBOL" → "SYMBOL" (direct)
    - "ENSG00001;ENSG00002;SYMBOL" → "SYMBOL" (semicolon-separated, take last)
    - "ENSG00001;ENSG00002;\"SYMBOL\"" → "SYMBOL" (quoted last element)
    """
    parts = str(name_val).replace('"', '').split(';')
    return parts[-1].strip() if len(parts) >= 3 else parts[0].strip()


def _normalize_gene_symbols(
    df: pd.DataFrame,
    col: str,
    log=None,
    context: str = "",
) -> pd.DataFrame:
    """Uppercase the gene column and drop case-converted duplicates.

    Enrichr's KEGG_2021_Human (and most other libraries) store gene
    symbols in all-uppercase, so mixed-case inputs from
    rank_genes_groups must be normalized before submission; otherwise
    hits are silently dropped and the analysis returns zero significant
    terms. Duplicates that collapse to the same symbol after
    uppercasing (e.g. "TP53" and "tp53") are removed — the first
    occurrence is preserved, matching the deterministic order of
    prerank/ORA inputs.

    Args:
        df: DataFrame containing a gene symbol column. Not modified.
        col: Name of the gene symbol column in `df`.
        log: Optional logger; emits a DEBUG message with the dedup
            count when duplicates are dropped.
        context: Short tag (e.g. cluster name) included in the log.

    Returns:
        New DataFrame with `col` uppercased and case-converted
        duplicates removed. NaN entries in `col` are also dropped.
    """
    out = df.copy()
    n_before = len(out)
    # Drop NaN/empty entries BEFORE str cast so None doesn't become "NONE"
    # and "" doesn't survive as an upper-cased empty string.
    mask_valid = out[col].notna() & out[col].astype(str).str.len().gt(0)
    out = out.loc[mask_valid].copy()
    out[col] = out[col].astype(str).str.upper()
    out = out.drop_duplicates(subset=col)
    n_removed = n_before - len(out)
    if n_removed > 0 and log is not None:
        log.debug(
            "  %s: removed %d duplicate genes after case conversion (original %d → unique %d)",
            context, n_removed, n_before, len(out),
        )
    return out


def read_marker_csv(table_dir: str, log) -> pd.DataFrame:
    """读取 Step 07 产出的标记基因 CSV"""
    path = os.path.join(table_dir, "marker_genes_per_group.csv")
    if not os.path.exists(path):
        log.error("Marker gene file not found: %s", path)
        log.error("Please run Step 07 (07_markers_de.py) first to generate this file.")
        sys.exit(1)
    df = pd.read_csv(path)
    log.info("Loaded marker genes: %d rows, %d groups",
             len(df), df['group'].nunique())
    return df


def _ora_one_group(
    grp,
    grp_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
):
    """Run ORA for a single group via Enrichr API (used by ThreadPoolExecutor)."""
    grp_up_df = (
        grp_df[grp_df['logfoldchanges'] > 0]
        .nsmallest(CFG.enrichment_n_top_genes, 'pvals_adj')
        [['names']]
        .copy()
    )
    grp_up_df['gene'] = grp_up_df['names'].apply(_extract_gene_symbol)
    grp_up_df = _normalize_gene_symbols(
        grp_up_df, 'gene', log=log, context=f"{grp} ORA"
    )
    grp_up = grp_up_df['gene'].tolist()
    if len(grp_up) < CFG.enrichment_min_size:
        log.info("  %s: insufficient up-regulated genes (%d < %d), skipping ORA",
                 grp, len(grp_up), CFG.enrichment_min_size)
        return (grp, None)

    try:
        enr = enrichr_with_retry(
            gene_list=grp_up,
            gene_sets_library=gene_set,
            log=log,
        )
    except Exception as e:
        log.warning("  %s ORA failed (%s): %s", grp, gene_set, e)
        return (grp, None)

    res = enr.results
    if res is None or len(res) == 0:
        return (grp, None)
    res['cluster'] = grp
    res['n_genes_input'] = len(grp_up)
    n_sig = (res['Adjusted P-value'] < CFG.enrichment_pval_cutoff).sum()
    log.info("  %s: %d/%d significant pathways (ORA, %s)",
             grp, n_sig, len(res), gene_set.split('_')[0])
    return (grp, res)


def run_ora(
    marker_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
) -> pd.DataFrame:
    """
    ORA (Over-Representation Analysis) via GSEApy Enrichr.

    对每个分组的 top N 上调基因（按 pvals_adj 排序），
    查询 Enrichr API 计算通路富集。
    """
    import gseapy as gp

    groups = marker_df['group'].unique()
    # ORA calls Enrichr API — conservative cap to avoid HTTP 429
    max_workers = min(5, getattr(CFG, 'n_jobs', 4))
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_grp = {
            executor.submit(
                _ora_one_group, grp,
                marker_df[marker_df['group'] == grp],
                gene_set, CFG, log,
            ): grp
            for grp in groups
        }
        for future in as_completed(future_to_grp):
            grp, res = future.result()
            results[grp] = res

    all_rows = [results[grp] for grp in groups if results.get(grp) is not None]

    if not all_rows:
        log.warning("  ORA no results (gene_set=%s)", gene_set)
        return pd.DataFrame()

    combined = pd.concat(all_rows, ignore_index=True)
    return combined


def _prerank_one_group(
    grp,
    grp_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
):
    """Run pre-ranked GSEA for a single group (used by ThreadPoolExecutor)."""
    import gseapy as gp

    df = grp_df.dropna(subset=['scores', 'names']).copy()
    df['gene'] = df['names'].apply(_extract_gene_symbol)
    if len(df) < CFG.enrichment_min_size:
        log.info("  %s: insufficient genes (%d < %d), skipping GSEA",
                 grp, len(df), CFG.enrichment_min_size)
        return (grp, None)

    df = _normalize_gene_symbols(df, 'gene', log=log, context=f"{grp} GSEA")
    if len(df) < CFG.enrichment_min_size:
        log.info("  %s: insufficient genes (%d < %d), skipping GSEA",
                 grp, len(df), CFG.enrichment_min_size)
        return (grp, None)

    rnk = df.set_index('gene')['scores']

    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_set,
            min_size=CFG.enrichment_min_size,
            max_size=CFG.enrichment_max_size,
            permutation_num=CFG.enrichment_permutations,
            threads=1,
            outdir=None,
            seed=CFG.random_seed,
            verbose=False,
            no_plot=True,
        )
    except Exception as e:
        log.warning("  %s GSEA failed (%s): %s", grp, gene_set, e)
        return (grp, None)

    res = pre_res.res2d
    if res is None or len(res) == 0:
        return (grp, None)
    res['cluster'] = grp
    res['n_genes_input'] = len(rnk)
    n_sig = (res['FDR q-val'] < CFG.enrichment_pval_cutoff).sum()
    log.info("  %s: %d/%d significant pathways (GSEA, %s)",
             grp, n_sig, len(res), gene_set.split('_')[0])
    return (grp, res)


def run_prerank(
    marker_df: pd.DataFrame,
    gene_set: str,
    CFG,
    log,
) -> pd.DataFrame:
    """
    Pre-ranked GSEA via GSEApy.

    对每个分组使用全部基因的 scores 作为排序指标，
    无需 cutoff，捕获微弱的协同变化。
    """
    import gseapy as gp

    groups = marker_df['group'].unique()
    # Prerank is local computation — no API rate limit, but GIL limits actual gain
    max_workers = getattr(CFG, 'n_jobs', 4)
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_grp = {
            executor.submit(
                _prerank_one_group, grp,
                marker_df[marker_df['group'] == grp],
                gene_set, CFG, log,
            ): grp
            for grp in groups
        }
        for future in as_completed(future_to_grp):
            grp, res = future.result()
            results[grp] = res

    all_rows = [results[grp] for grp in groups if results.get(grp) is not None]

    if not all_rows:
        log.warning(
            "  GSEA no results for '%s': no gene sets passed filtering. "
            "Consider adjusting enrichment_min_size=%d / enrichment_max_size=%d "
            "or using a different gene set library.",
            gene_set, CFG.enrichment_min_size, CFG.enrichment_max_size,
        )
        return pd.DataFrame()

    combined = pd.concat(all_rows, ignore_index=True)
    return combined


def save_results(
    ora_results: dict,
    prerank_results: dict,
    CFG,
    log,
) -> None:
    """保存富集结果 CSV + 绘图"""
    table_dir = os.path.join(CFG.table_dir, "enrichment")
    fig_dir = os.path.join(CFG.figure_dir, "enrichment")
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # ── 保存 ORA CSV ──
    for gs_name, df in ora_results.items():
        if df.empty:
            continue
        path = os.path.join(table_dir, f"ora_{gs_name}_summary.csv")
        df.to_csv(path, index=False)
        log.info("  ORA results exported: %s (%d rows)", path, len(df))

    # ── 保存 GSEA CSV ──
    for gs_name, df in prerank_results.items():
        if df.empty:
            continue
        path = os.path.join(table_dir, f"prerank_{gs_name}_summary.csv")
        df.to_csv(path, index=False)
        log.info("  GSEA results exported: %s (%d rows)", path, len(df))

    # ── 气泡图（每组 top 20 通路） ──
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        for gs_name, df in ora_results.items():
            if df.empty:
                continue
            plot_ora_bubble(df, gs_name, fig_dir, CFG, log)

        for gs_name, df in prerank_results.items():
            if df.empty:
                continue
            plot_prerank_bubble(df, gs_name, fig_dir, CFG, log)
    except Exception as e:
        log.warning("Plotting failed: %s", e)


def plot_ora_bubble(
    df: pd.DataFrame,
    gs_name: str,
    fig_dir: str,
    CFG,
    log,
) -> None:
    """ORA 气泡图: x=cluster, y=Term, size=Overlap, color=Adjusted P-value"""
    import matplotlib.pyplot as plt

    sig = df[df['Adjusted P-value'] < CFG.enrichment_pval_cutoff].copy()
    if sig.empty:
        sig = df.head(5)
    top_per_cluster = (
        sig.sort_values('Adjusted P-value')
        .groupby('cluster', observed=True)
        .head(CFG.enrichment_n_top_genes // max(1, sig['cluster'].nunique()))
    )
    if len(top_per_cluster) < 3:
        log.info("  Skipping bubble plot (%s): insufficient significant pathways", gs_name)
        return

    top_per_cluster['-log10_padj'] = -np.log10(
        top_per_cluster['Adjusted P-value'].clip(lower=1e-300)
    )
    # 简短 Term 名称
    top_per_cluster['Term_short'] = top_per_cluster['Term'].str.replace(
        r'\s*\(GO:\d+\)$', '', regex=True
    ).str[:60]

    fig, ax = plt.subplots(figsize=(
        max(8, 0.5 * top_per_cluster['cluster'].nunique()),
        max(6, 0.3 * top_per_cluster['Term_short'].nunique()),
    ))
    overlap_numeric = (
        top_per_cluster['Overlap']
        .astype(str).str.split('/').str[0].astype(float)
    )
    sc = ax.scatter(
        top_per_cluster['cluster'],
        top_per_cluster['Term_short'],
        s=overlap_numeric * 30,
        c=top_per_cluster['-log10_padj'],
        cmap='YlOrRd',
        edgecolors='grey', linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label='-log10(Adjusted P-value)')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('')
    ax.set_title(f'Enrichment: {gs_name}')
    fig.tight_layout()
    path = os.path.join(fig_dir, f"ora_{gs_name}_bubble.pdf")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  ORA bubble plot: %s", path)


def plot_prerank_bubble(
    df: pd.DataFrame,
    gs_name: str,
    fig_dir: str,
    CFG,
    log,
) -> None:
    """Pre-ranked GSEA 气泡图: color=NES, size=-log10(FDR)"""
    import matplotlib.pyplot as plt

    sig = df[df['FDR q-val'] < CFG.enrichment_pval_cutoff].copy()
    if sig.empty:
        sig = df.head(10)
    top_per_cluster = (
        sig.sort_values('FDR q-val')
        .groupby('cluster', observed=True)
        .head(CFG.enrichment_n_top_genes // max(1, sig['cluster'].nunique()))
    )
    if len(top_per_cluster) < 3:
        return

    top_per_cluster['Term_short'] = top_per_cluster['Term'].str[:60]
    top_per_cluster['-log10_fdr'] = -np.log10(
        top_per_cluster['FDR q-val'].clip(lower=1e-300)
    )
    # NES 颜色: 红色=上调, 蓝色=下调
    vmax = max(abs(top_per_cluster['NES'].min()),
               abs(top_per_cluster['NES'].max()))
    norm = plt.Normalize(-vmax, vmax)

    fig, ax = plt.subplots(figsize=(
        max(8, 0.5 * top_per_cluster['cluster'].nunique()),
        max(6, 0.3 * top_per_cluster['Term_short'].nunique()),
    ))
    sc = ax.scatter(
        top_per_cluster['cluster'],
        top_per_cluster['Term_short'],
        s=top_per_cluster['-log10_fdr'] * 20,
        c=top_per_cluster['NES'],
        cmap='RdBu_r', norm=norm,
        edgecolors='grey', linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label='NES')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('')
    ax.set_title(f'GSEA: {gs_name}')
    fig.tight_layout()
    path = os.path.join(fig_dir, f"prerank_{gs_name}_bubble.pdf")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info("  GSEA bubble plot: %s", path)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()
    CFG = resolve_config(args.config)
    log = setup_logger("09_enrichment",
                        os.path.join(CFG.log_dir, "09_enrichment.log"))
    log.info("Step 09: Enrichment analysis (GO/KEGG)")

    if not CFG.run_enrichment:
        log.info("Enrichment analysis disabled (run_enrichment=False)")
        return

    marker_df = read_marker_csv(CFG.table_dir, log)
    log.info("Gene set libraries: %s", CFG.enrichment_gene_sets)
    log.info("Method: %s", CFG.enrichment_method)

    # Quality awareness (v3.1.0+): check marker_validation from annotated h5ad
    try:
        _annotated_path = os.path.join(CFG.h5ad_dir, '05_annotated.h5ad')
        _quality_path = os.path.join(CFG.table_dir, '05_annotation_quality.json')
        _pass_rate = None
        if os.path.exists(_quality_path):
            import json as _json
            with open(_quality_path, 'r') as _f:
                _q = _json.load(_f)
            _pass_rate = _q.get('pass_rate')
        if _pass_rate is None and os.path.exists(_annotated_path):
            _a = sc.read(_annotated_path)
            if 'marker_validation' in _a.obs and _a.n_obs > 0:
                _pass_cells = (_a.obs['marker_validation'] == 'PASS').sum()
                _pass_rate = _pass_cells / _a.n_obs
        if _pass_rate is not None:
            _pass_rate_min = getattr(CFG, 'marker_validation_pass_rate_min', 0.1)
            if _pass_rate < _pass_rate_min:
                log.warning(
                    "⚠  marker_validation PASS rate %.1f%% (<%.0f%%) — "
                    "enrichment analysis is based on potentially unreliable "
                    "cell_type labels. Results should be interpreted with caution.",
                    _pass_rate * 100, _pass_rate_min * 100,
                )
    except Exception:
        pass  # Non-critical — enrichment proceeds regardless

    # ── 按基因集库循环 ──
    ora_results = {}
    prerank_results = {}

    for gs in CFG.enrichment_gene_sets:
        gs_name = gs.replace(' ', '_').replace('/', '_')

        if CFG.enrichment_method in ('ora', 'both'):
            log.info("[ORA] Gene set: %s", gs)
            ora_df = run_ora(marker_df, gs, CFG, log)
            if ora_df is not None and not ora_df.empty:
                ora_results[gs_name] = ora_df

        if CFG.enrichment_method in ('prerank', 'both'):
            log.info("[GSEA] Gene set: %s", gs)
            prerank_df = run_prerank(marker_df, gs, CFG, log)
            if prerank_df is not None and not prerank_df.empty:
                prerank_results[gs_name] = prerank_df
            else:
                log.info("GSEA mode yielded no results, ORA results available separately")

    total_ora = sum(len(df) for df in ora_results.values())
    total_gsea = sum(len(df) for df in prerank_results.values())
    log.info("Enrichment results summary: ORA %d rows, GSEA %d rows",
             total_ora, total_gsea)

    save_results(ora_results, prerank_results, CFG, log)

    # ── AI Biological Interpretation (optional) ──
    if CFG.ai.enabled and CFG.ai.ai_interpretation:
        log.info("AI: Generating biological interpretation...")
        try:
            summary_data = []
            for gs_name, df in ora_results.items():
                if df.empty:
                    continue
                sig = df[df['Adjusted P-value'] < CFG.enrichment_pval_cutoff]
                for cluster in sig['cluster'].unique():
                    cluster_sig = sig[sig['cluster'] == cluster].head(5)
                    summary_data.append({
                        "gene_set": gs_name,
                        "cluster": str(cluster),
                        "top_terms": cluster_sig[['Term', 'Adjusted P-value']].to_dict('records')
                    })

            if summary_data:
                system_prompt = "You are an expert computational biologist interpreting scRNA-seq enrichment results."
                user_prompt = f"Enrichment results summary:\n{json.dumps(summary_data, indent=2)}\n\nProvide biological interpretation: key pathways, cross-cell-type patterns, and testable hypotheses."

                from core.ai_caller import ai_query
                interpretation = ai_query(system_prompt, user_prompt, cfg=CFG.ai)

                interp_path = os.path.join(CFG.table_dir, "enrichment", "ai_interpretation.txt")
                os.makedirs(os.path.dirname(interp_path), exist_ok=True)
                with open(interp_path, "w") as f:
                    f.write(interpretation)
                log.info("AI interpretation saved to %s", interp_path)

                summary_lines = [f"Biological Interpretation — {'scRNA-seq enrichment'}"]
                summary_lines.append("=" * 60)
                summary_lines.append(interpretation[:2000])
                summary_path = os.path.join(CFG.table_dir, "enrichment", "ai_interpretation_summary.txt")
                with open(summary_path, "w") as f:
                    f.write("\n".join(summary_lines))
        except Exception as e:
            log.warning("AI interpretation skipped: %s", e)

    log.info("Step 09 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
