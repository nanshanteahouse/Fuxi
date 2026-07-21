"""Automatic biological sex detection from scRNA-seq data.

Scans sex-linked gene expression (Xist/XIST for female, Y-linked
genes for male) to predict per-cell biological sex.
"""

import os

import numpy as np
import pandas as pd

FEMALE_GENES = ["Xist", "XIST"]
MALE_GENES = [
    "Eif2s3y",
    "Ddx3y",
    "Uty",
    "Kdm5d",  # mouse
    "RPS4Y1",
    "DDX3Y",
    "UTY",
    "KDM5D",  # human
    "EIF2S3Y",
    "SRY",
    "ZFY",
    "RPS4Y2",  # additional human
]


def detect_sex(adata, cfg, log):
    """Auto-detect biological sex by scanning sex-linked gene expression.

    Adds a ``predicted_sex`` column to ``adata.obs`` with values
    ``'Female'``, ``'Male'``, ``'Ambiguous'``, or ``'Unknown'``,
    and saves a detailed ``sex_report.csv`` to ``CFG.table_dir``.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with optional ``.raw`` attribute.
    CFG : Config
        Pipeline config (uses ``CFG.table_dir``).
    log : logging.Logger
        Logger for progress and warnings.
    """
    # ── 1. Already annotated? ──────────────────────────────────────────
    for col in ("sex", "gender"):
        if col in adata.obs:
            log.info("Sex column '%s' already present in adata.obs", col)
            counts = adata.obs[col].value_counts()
            for val, cnt in counts.items():
                log.info("  %s: %d (%.1f%%)", val, cnt, cnt / adata.n_obs * 100)
            return

    # ── 2. Need raw data ───────────────────────────────────────────────
    if adata.raw is None:
        log.warning("No raw data available — skipping sex detection")
        return

    # ── 3. Locate sex-linked genes ─────────────────────────────────────
    female_found = [g for g in FEMALE_GENES if g in adata.raw.var_names]
    male_found = [g for g in MALE_GENES if g in adata.raw.var_names]

    if not female_found and not male_found:
        log.warning("No sex-linked genes found in gene set")
        return

    # Species hint
    if "Xist" in female_found:
        log.info("Sex detection using mouse panel (Xist + Y-linked)")
    elif "XIST" in female_found:
        log.info("Sex detection using human panel (XIST + Y-linked)")
    else:
        log.info("Sex detection using mixed/cross-species panel")

    # ── 4. Per-cell classification ─────────────────────────────────────
    def _expr_positive(gene):
        """Return bool mask of cells with UMI > 0 for *gene*."""
        gene_idx = list(adata.raw.var_names).index(gene)
        expr = adata.raw.X[:, gene_idx]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = np.asarray(expr).flatten()
        return expr > 0

    female_mask = np.zeros(adata.n_obs, dtype=bool)
    for g in female_found:
        female_mask |= _expr_positive(g)

    male_mask = np.zeros(adata.n_obs, dtype=bool)
    for g in male_found:
        male_mask |= _expr_positive(g)

    predicted = np.full(adata.n_obs, "Unknown", dtype=object)
    predicted[female_mask & ~male_mask] = "Female"
    predicted[male_mask & ~female_mask] = "Male"
    predicted[female_mask & male_mask] = "Ambiguous"

    adata.obs["predicted_sex"] = predicted

    # ── 5. Logging ─────────────────────────────────────────────────────
    n_female = int((predicted == "Female").sum())
    n_male = int((predicted == "Male").sum())
    n_ambig = int((predicted == "Ambiguous").sum())
    n_total = adata.n_obs

    log.info(
        "Sex detection: Female=%.1f%% (%d/%d), Male=%.1f%% (%d/%d), Ambiguous=%.1f%% (%d/%d)",
        n_female / n_total * 100,
        n_female,
        n_total,
        n_male / n_total * 100,
        n_male,
        n_total,
        n_ambig / n_total * 100,
        n_ambig,
        n_total,
    )

    female_ratio = n_female / n_total
    male_ratio = n_male / n_total
    if female_ratio > 0.05 and male_ratio > 0.05:
        _genes = sorted(set(female_found + male_found))
        _gene_hint = ", ".join(repr(g) for g in _genes[:6])
        if len(_genes) > 6:
            _gene_hint += ", ..."
        log.warning("Mixed-sex dataset detected — sex may act as batch effect.")
        log.warning("  Option 1 (mild):  CFG.integration.batch_key = 'predicted_sex'")
        log.warning("  Option 2 (strong): CFG.normalization.regress_out_genes = [%s]", _gene_hint)

    if "sample" in adata.obs:
        for sample_name, group in adata.obs.groupby("sample"):
            sex_counts = group["predicted_sex"].value_counts()
            total = len(group)
            parts = [f"{v}: {c} ({c / total * 100:.1f}%)" for v, c in sex_counts.items()]
            log.info("  Sample %s: %s", sample_name, ", ".join(parts))

    # ── 6. Save CSV report ─────────────────────────────────────────────
    report = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "predicted_sex": predicted,
            "female_marker_positive": female_mask,
            "male_marker_positive": male_mask,
        }
    )
    os.makedirs(cfg.table_dir, exist_ok=True)
    report_path = os.path.join(cfg.table_dir, "sex_report.csv")
    report.to_csv(report_path, index=False)
    log.info("Sex report saved: %s", report_path)
