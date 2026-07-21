#!/usr/bin/env python3
"""Classify adata.obs columns as batch/biology via PCA variance decomposition,
Gini coefficient, cluster purity, Cramer's V collinearity, and permutation tests.

Exports: ColumnDiagnosis, BatchDiagnosisReport, diagnose_batch_candidates,
validate_harmony_preservation, plot_diagnosis_report,
_compute_anova_r2_per_pc, _compute_gini_criterion, _compute_purity_one_shot,
_compute_cramer_v, _permutation_test.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats as stats

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════
#  Data classes
# ═══════════════════════════════════════════════════════════


@dataclass
class ColumnDiagnosis:
    column: str
    gini_criterion: float
    purity_score: float
    n_unique: int
    cramer_v: dict[str, float]
    judgment: str  # "batch" | "biology" | "ambiguous" | "skip"
    permutation_pval: float | None = None
    recommendation: str = ""


@dataclass
class BatchDiagnosisReport:
    column_diagnoses: list[ColumnDiagnosis] = field(default_factory=list)
    batch_cols: list[str] = field(default_factory=list)
    biology_cols: list[str] = field(default_factory=list)
    ambiguous_cols: list[str] = field(default_factory=list)
    suggested_batch_key: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════


def _categorize(
    gini: float,
    perm_pval: float | None,
    gini_batch_threshold: float,
    gini_biology_threshold: float,
) -> str:
    if gini <= gini_batch_threshold:
        return "batch"
    if gini >= gini_biology_threshold:
        return "biology"
    if perm_pval is not None and perm_pval < 0.05:
        return "biology"
    return "ambiguous"


# ═══════════════════════════════════════════════════════════
#  Core diagnostics
# ═══════════════════════════════════════════════════════════


def _compute_anova_r2_per_pc(pca_matrix: np.ndarray, col_values) -> list[float]:
    """ANOVA R² per PC via f_oneway."""
    unique_vals = np.unique(col_values)
    groups = [pca_matrix[col_values == v] for v in unique_vals]
    n_pcs = pca_matrix.shape[1]
    r2s = []
    for pc in range(n_pcs):
        pc_groups = [g[:, pc] for g in groups]
        try:
            f_stat, _ = stats.f_oneway(*pc_groups)
        except Exception:
            r2s.append(0.0)
            continue
        df_between = len(groups) - 1
        df_error = len(col_values) - len(groups)
        if df_error <= 0 or f_stat < 0:
            r2s.append(0.0)
        else:
            r2 = (f_stat * df_between) / (f_stat * df_between + df_error)
            r2s.append(float(r2))
    return r2s


def _compute_gini_criterion(r2_array: np.ndarray) -> float:
    """Gini = 2*sum(i*ri)/(n*sum(r2)) - (n+1)/n."""
    r2_sorted = np.sort(r2_array)
    n = len(r2_sorted)
    total = float(np.sum(r2_sorted))
    if total == 0.0:
        return 0.0
    gini = (2.0 * np.sum((np.arange(1, n + 1)) * r2_sorted)) / (n * total) - (n + 1) / n
    return float(gini)


def _compute_purity_one_shot(adata, col: str, use_rep: str = "X_pca") -> float:
    """Cluster purity: copy internally, kNN + Leiden, crosstab max sum / total."""
    if use_rep not in adata.obsm:
        raise ValueError(f"'{use_rep}' not found in adata.obsm. Run PCA first.")
    if adata.n_obs < 3:
        return 1.0
    adata = adata.copy()
    try:
        sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep)
        sc.tl.leiden(
            adata,
            resolution=1.0,
            key_added="_diag_leiden",
            flavor="igraph",
            directed=False,
            n_iterations=2,
        )
    except Exception:
        return 1.0
    ct = pd.crosstab(adata.obs["_diag_leiden"], adata.obs[col])
    return float(ct.max(axis=1).sum() / ct.values.sum())


def _compute_cramer_v(col_a: pd.Series, col_b: pd.Series) -> float:
    """Cramer's V = sqrt(chi2 / (n * min(k-1, r-1)))."""
    ct = pd.crosstab(col_a, col_b)
    chi2 = float(stats.chi2_contingency(ct, correction=False)[0])
    n = int(ct.values.sum())
    k, r = ct.shape
    min_dim = min(k - 1, r - 1)
    if min_dim == 0 or n == 0:
        return 0.0
    return float(math.sqrt(chi2 / (n * min_dim)))


def _permutation_test(
    pca_matrix: np.ndarray,
    col_values,
    real_r2s: np.ndarray,
    n_perm: int = 100,
    random_state: int = 42,
) -> float:
    """p-value = fraction of permutations where mean R² exceeds real mean."""
    rng = np.random.RandomState(random_state)
    real_mean = float(np.mean(real_r2s))
    count = 0
    col_arr = np.array(col_values).copy()
    for _ in range(n_perm):
        rng.shuffle(col_arr)
        perm_r2s = _compute_anova_r2_per_pc(pca_matrix, col_arr)
        if float(np.mean(perm_r2s)) > real_mean:
            count += 1
    return count / n_perm


# ═══════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════


def _build_rec(judgment: str, col: str, gini: float, purity: float) -> str:
    if judgment == "batch":
        return f"Use '{col}' as batch key for Harmony (gini={gini:.3f}, purity={purity:.3f})"
    elif judgment == "biology":
        return f"'{col}' appears biological \u2014 do NOT use as batch key (gini={gini:.3f}, purity={purity:.3f})"
    elif judgment == "ambiguous":
        return f"'{col}' ambiguous \u2014 manual review (gini={gini:.3f}, purity={purity:.3f})"
    return f"'{col}' skipped (single value)."


def diagnose_batch_candidates(
    adata,
    n_pcs: int = 50,
    random_state: int = 42,
    gini_batch_threshold: float = 0.3,
    gini_biology_threshold: float = 0.6,
    permute_n: int = 100,
) -> BatchDiagnosisReport:
    """Classify every categorical obs column as batch / biology / ambiguous."""
    if "X_pca" not in adata.obsm:
        raise ValueError("'X_pca' not found in adata.obsm. Run sc.pp.pca first.")

    pca_mat = adata.obsm["X_pca"][:, : min(n_pcs, adata.obsm["X_pca"].shape[1])]

    cat_cols = [
        col for col in adata.obs.columns if isinstance(adata.obs[col].dtype, pd.CategoricalDtype)
    ]

    diagnoses: list[ColumnDiagnosis] = []
    batch_cols: list[str] = []
    biology_cols: list[str] = []
    ambiguous_cols: list[str] = []
    warnings: list[str] = []

    if not cat_cols:
        return BatchDiagnosisReport(suggested_batch_key=[], warnings=warnings)

    for col in cat_cols:
        n_unique = len(adata.obs[col].dropna().unique())
        if n_unique < 2:
            continue

        col_vals = adata.obs[col].values
        r2s = _compute_anova_r2_per_pc(pca_mat, col_vals)
        r2_arr = np.array(r2s)
        gini = _compute_gini_criterion(r2_arr)
        purity = _compute_purity_one_shot(adata, col, use_rep="X_pca")

        cramer_v: dict[str, float] = {}
        for oc in cat_cols:
            if oc != col:
                v = _compute_cramer_v(adata.obs[col], adata.obs[oc])
                cramer_v[oc] = v
                if v >= 1.0:
                    warnings.append(
                        f"Column '{col}' is perfectly collinear with '{oc}' "
                        f"(V=1.0) \u2014 redundant column."
                    )

        judgment: str = "ambiguous"
        perm_pval: float | None = None

        if gini <= gini_batch_threshold:
            judgment = "batch"
            perm_pval = None
        else:
            perm_pval = _permutation_test(
                pca_mat, col_vals, r2_arr, n_perm=permute_n, random_state=random_state
            )
            judgment = _categorize(gini, perm_pval, gini_batch_threshold, gini_biology_threshold)

        diagnoses.append(
            ColumnDiagnosis(
                column=col,
                gini_criterion=gini,
                purity_score=purity,
                n_unique=n_unique,
                cramer_v=cramer_v,
                judgment=judgment,
                permutation_pval=perm_pval,
                recommendation=_build_rec(judgment, col, gini, purity),
            )
        )

        if judgment == "batch":
            batch_cols.append(col)
        elif judgment == "biology":
            biology_cols.append(col)
        elif judgment == "ambiguous":
            ambiguous_cols.append(col)

    return BatchDiagnosisReport(
        column_diagnoses=diagnoses,
        batch_cols=batch_cols,
        biology_cols=biology_cols,
        ambiguous_cols=ambiguous_cols,
        suggested_batch_key=list(batch_cols),
        warnings=warnings,
    )


def validate_harmony_preservation(
    adata_before, adata_after, biology_cols: list[str]
) -> dict[str, float]:
    """purity_after/purity_before per biology col using X_pca vs X_integrated."""
    results: dict[str, float] = {}
    for col in biology_cols:
        pb = _compute_purity_one_shot(adata_before, col, use_rep="X_pca")
        pa = _compute_purity_one_shot(adata_after, col, use_rep="X_integrated")
        results[col] = pa / pb if pb > 0 else 1.0
    return results


def plot_diagnosis_report(report: BatchDiagnosisReport, save_path: str) -> None:
    """Cramer's V heatmap + Gini barplot + recommendation panel."""
    import matplotlib.pyplot as plt

    if not report.column_diagnoses:
        logger.info("No diagnoses to plot.")
        return

    cols = [d.column for d in report.column_diagnoses]
    if not cols:
        return

    n = len(cols)
    fig, axes = plt.subplots(1, 3, figsize=(18, max(5, n * 0.5 + 2)))

    # 1 — Cramer's V heatmap
    vmat = np.zeros((n, n))
    dmap = {d.column: d for d in report.column_diagnoses}
    for i, ci in enumerate(cols):
        for j, cj in enumerate(cols):
            vmat[i, j] = 1.0 if i == j else dmap[ci].cramer_v.get(cj, 0.0)

    ax0 = axes[0]
    im = ax0.imshow(vmat, cmap="Reds", vmin=0, vmax=1, aspect="auto")
    ax0.set_xticks(range(n))
    ax0.set_yticks(range(n))
    ax0.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
    ax0.set_yticklabels(cols, fontsize=8)
    ax0.set_title("Cramer's V", fontsize=12)
    fig.colorbar(im, ax=ax0, shrink=0.8)

    # 2 — Gini barplot
    cmap = {"batch": "#e74c3c", "biology": "#2ecc71", "ambiguous": "#f39c12"}
    ginis = [d.gini_criterion for d in report.column_diagnoses]
    colors = [cmap.get(d.judgment, "#95a5a6") for d in report.column_diagnoses]

    ax1 = axes[1]
    ax1.barh(range(len(ginis)), ginis, color=colors, alpha=0.8)
    ax1.set_yticks(range(len(ginis)))
    ax1.set_yticklabels(cols, fontsize=8)
    ax1.set_xlabel("Gini coefficient")
    ax1.set_title("Signal concentration", fontsize=12)
    ax1.axvline(0.3, color="red", ls="--", alpha=0.5, label="batch threshold")
    ax1.axvline(0.6, color="green", ls="--", alpha=0.5, label="biology threshold")
    ax1.legend(fontsize=8)
    ax1.set_xlim(0, 1)

    # 3 — Text panel
    ax2 = axes[2]
    ax2.axis("off")
    lines: list[str] = ["Recommendations:"]
    for d in report.column_diagnoses:
        lines.append(f"  \u2022 {d.column}: {d.judgment}")
        if d.recommendation:
            lines.append(f"    {d.recommendation}")
    if report.warnings:
        lines.extend(["", "Warnings:"] + [f"  \u26a0 {w}" for w in report.warnings])

    ax2.text(
        0,
        0.95,
        "\n".join(lines),
        transform=ax2.transAxes,
        fontsize=9,
        verticalalignment="top",
        family="monospace",
    )

    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Diagnosis report saved to %s", save_path)
