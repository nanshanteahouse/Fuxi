"""Tests for load_meta persistence in perf_report.json + model-based memory estimate."""

import json

import pytest

from core.utils._perf import PerformanceReport, PerformanceSummary


def _summary_with_load_meta(load_meta):
    s = PerformanceSummary()
    s.pipeline_info = {"modality": "rna", "load_meta": load_meta}
    s.steps.append(
        PerformanceReport(
            step="00 Load raw data",
            wall_sec=10.0,
            cpu_sec=8.0,
            peak_rss_mib=5000.0,
            avg_cpu_pct=80.0,
            gpu_mem_mb=-1.0,
            n_cells=70_477,
            n_genes=36_601,
            checkpoint_mib=700.0,
            exit_status="completed",
        )
    )
    return s


def test_to_dict_roundtrip_preserves_load_meta(tmp_path):
    lm = {"n_cells": 70_477, "n_genes": 36_601, "nnz": 426_241_129, "format": "10X_h5"}
    s = _summary_with_load_meta(lm)
    p = tmp_path / "perf_report.json"
    s.save_json(str(p))
    loaded = PerformanceSummary.load_existing(str(p))
    assert loaded is not None
    assert loaded.pipeline_info["load_meta"] == lm


def test_load_existing_without_load_meta_is_empty(tmp_path):
    p = tmp_path / "perf_report.json"
    p.write_text(json.dumps({"pipeline": {"modality": "rna"}, "steps": [], "summary": {}}))
    loaded = PerformanceSummary.load_existing(str(p))
    assert loaded is not None
    assert "load_meta" not in loaded.pipeline_info


def test_estimate_memory_model_branch_uses_estimate_step_peak():
    lm = {"n_cells": 70_477, "n_genes": 36_601, "nnz": 426_241_129, "format": "10X_h5"}
    est = PerformanceSummary._estimate_memory(7000.0, 70_477, 36_601, load_meta=lm)
    # Multiome: cf=1.0 estimate 7.13 GiB at 70k cells; scales with cells.
    assert set(est) == {"50k", "100k", "200k", "500k"}
    # model-based values are monotonic in cell count
    vals = [est[k] for k in ("50k", "100k", "200k", "500k")]
    assert all(b >= a for a, b in zip(vals, vals[1:]))


def test_estimate_memory_fallback_matches_legacy_linear(tmp_path):
    # No load_meta -> legacy per-cell-gene linear extrapolation (byte-stable).
    est = PerformanceSummary._estimate_memory(5000.0, 50_000, 36_601)
    assert set(est) == {"50k", "100k", "200k", "500k"}
    # hand-computed legacy value: 5000 / (50000*36601) * 50000 * 36601 / 1024 = 4.9
    assert est["50k"] == pytest.approx(5000.0 / 1024, rel=0.01)
    assert est["100k"] == pytest.approx(5000.0 * 2 / 1024, rel=0.01)


def test_estimate_memory_zero_cells_returns_empty():
    assert PerformanceSummary._estimate_memory(0.0, 0, 0) == {}
    assert PerformanceSummary._estimate_memory(5000.0, 0, 36_601) == {}
    assert PerformanceSummary._estimate_memory(5000.0, 50_000, 0) == {}


def test_load_meta_nnz_scales_by_density():
    lm = {"n_cells": 50_000, "n_genes": 10_000, "nnz": 50_000_000, "format": "10X_h5"}
    est_low = PerformanceSummary._estimate_memory(3000.0, 50_000, 10_000, load_meta=lm)
    est_high = PerformanceSummary._estimate_memory(
        3000.0, 50_000, 10_000, load_meta={**lm, "nnz": 100_000_000}
    )
    # Same cell count but double density (nnz) -> strictly larger estimate.
    assert est_high["50k"] > est_low["50k"]
