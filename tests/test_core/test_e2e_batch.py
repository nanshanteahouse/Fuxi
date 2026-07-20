"""End-to-end batch diagnostics with synthetic AnnData — no real GEO data required."""

import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy  # harmony integration (no torch)
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from rna.utils.batch_diagnostics import (
    diagnose_batch_candidates,
    validate_harmony_preservation,
    plot_diagnosis_report,
    BatchDiagnosisReport,
)

sc.settings.verbosity = 0


def scenario_1() -> None:
    """sample (random) -> batch.  tissue (structured) -> biology."""
    rng = np.random.RandomState(42)
    n_cells, n_genes = 350, 20
    X_raw = rng.randn(n_cells, n_genes).astype(np.float32)
    sample = rng.choice(["S1", "S2", "S3", "S4"], n_cells)
    tissue = rng.choice(["TissueA", "TissueB"], n_cells)
    X_raw[tissue == "TissueA", :5] += 2.0
    X_raw[tissue == "TissueB", :5] -= 1.0
    adata = ad.AnnData(
        X=csr_matrix(X_raw),
        obs=pd.DataFrame({"sample": pd.Categorical(sample), "tissue": pd.Categorical(tissue)}),
    )
    sc.pp.pca(adata, n_comps=10, random_state=42)
    report = diagnose_batch_candidates(adata)
    assert "sample" in report.batch_cols, f"FAIL 1: 'sample' not in batch_cols {report.batch_cols}"
    print("PASS 1: 'sample' classified as batch")
    assert "tissue" in report.biology_cols, f"FAIL 1: 'tissue' not in biology_cols {report.biology_cols}"
    print("PASS 1: 'tissue' classified as biology")


def scenario_2() -> None:
    """Identical batch_x / disease columns -> Cramer's V = 1.0 warning."""
    rng = np.random.RandomState(42)
    n_cells, n_genes = 100, 10
    labels = rng.choice(["X", "Y"], n_cells)
    X_raw = rng.randn(n_cells, n_genes).astype(np.float32)
    X_raw[labels == "X", :3] += 1.5
    X_raw[labels == "Y", :3] -= 1.0
    adata = ad.AnnData(
        X=csr_matrix(X_raw),
        obs=pd.DataFrame({"batch_x": pd.Categorical(labels), "disease": pd.Categorical(labels.copy())}),
    )
    sc.pp.pca(adata, n_comps=5, random_state=42)
    report = diagnose_batch_candidates(adata)
    has_collinear = any("collinear" in w.lower() for w in report.warnings)
    assert has_collinear, f"FAIL 2: no collinearity warning: {report.warnings}"
    print("PASS 2: collinearity warning emitted")


def scenario_3() -> None:
    """No categorical obs -> empty BatchDiagnosisReport, no crash."""
    rng = np.random.RandomState(42)
    adata = ad.AnnData(
        X=csr_matrix(rng.randn(30, 8).astype(np.float32)),
        obs=pd.DataFrame({"n_counts": rng.poisson(1000, 30).astype(float)}),
    )
    sc.pp.pca(adata, n_comps=5, random_state=42)
    report = diagnose_batch_candidates(adata)
    assert isinstance(report, BatchDiagnosisReport), "FAIL 3: wrong type"
    assert len(report.column_diagnoses) == 0, f"FAIL 3: expected 0 diagnoses, got {len(report.column_diagnoses)}"
    print("PASS 3: numeric-only obs returns empty report gracefully")


def scenario_4() -> None:
    """Harmonize with wrong batch_key degrades purity; correct key preserves it."""
    rng = np.random.RandomState(42)
    n_cells, n_genes = 300, 15
    X = rng.randn(n_cells, n_genes).astype(np.float32)
    sample = rng.choice(["S1", "S2", "S3", "S4"], n_cells)
    for i, s in enumerate(["S1", "S2", "S3", "S4"]):
        X[sample == s, :4] += (i - 1.5) * 2.0
    tissue = rng.choice(["TissueA", "TissueB"], n_cells)
    X[tissue == "TissueA", 4:7] += 2.0
    X[tissue == "TissueB", 4:7] -= 1.0
    obs = pd.DataFrame({"sample": pd.Categorical(sample), "tissue": pd.Categorical(tissue)})
    adata = ad.AnnData(X=csr_matrix(X), obs=obs)
    sc.pp.pca(adata, n_comps=8, random_state=42)
    before = adata.copy()
    # Wrong batch_key: tissue (overcorrects biology)
    ad_wrong = adata.copy()
    ho_w = harmonypy.run_harmony(before.obsm["X_pca"], obs, "tissue", verbose=False, random_state=42)
    ad_wrong.obsm["X_pca_harmony"] = ho_w.Z_corr
    rw = validate_harmony_preservation(before, ad_wrong, ["tissue"]).get("tissue", 0.0)
    assert rw < 0.9, f"FAIL 4: wrong batch should degrade purity, got {rw:.4f}"
    print(f"PASS 4: wrong batch_key degrades purity ({rw:.4f} < 0.9)")
    # Correct batch_key: sample (preserves biology)
    ad_correct = adata.copy()
    ho_c = harmonypy.run_harmony(before.obsm["X_pca"], obs, "sample", verbose=False, random_state=42)
    ad_correct.obsm["X_pca_harmony"] = ho_c.Z_corr
    rc = validate_harmony_preservation(before, ad_correct, ["tissue"]).get("tissue", 0.0)
    assert rc >= 0.85, f"FAIL 4: correct batch should preserve purity, got {rc:.4f}"
    print(f"PASS 4: correct batch_key preserves purity ({rc:.4f} >= 0.85)")


def scenario_5() -> None:
    """plot_diagnosis_report writes a non-empty PDF."""
    rng = np.random.RandomState(42)
    adata = ad.AnnData(
        X=csr_matrix(rng.randn(100, 12).astype(np.float32)),
        obs=pd.DataFrame({"batch": pd.Categorical(rng.choice(["A", "B"], 100)),
                          "bio": pd.Categorical(rng.choice(["X", "Y"], 100))}),
    )
    sc.pp.pca(adata, n_comps=5, random_state=42)
    report = diagnose_batch_candidates(adata)
    save_path = "/tmp/e2e_test_diag.pdf"
    plot_diagnosis_report(report, save_path)
    assert os.path.exists(save_path), "FAIL 5: PDF not created"
    assert os.path.getsize(save_path) > 0, "FAIL 5: PDF is empty"
    os.remove(save_path)
    print("PASS 5: plot_diagnosis_report saves valid PDF")


if __name__ == "__main__":
    print("=" * 52)
    print("  E2E Batch Diagnostics - Synthetic AnnData")
    print("=" * 52)
    scenario_1(); scenario_2(); scenario_3()
    scenario_4(); scenario_5()
    print("\n" + "=" * 52)
    print("  ALL 5 SCENARIOS PASSED")
    print("=" * 52)
