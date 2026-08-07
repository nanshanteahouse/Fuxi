"""Modality-aware memory estimator registry tests.

Covers the ``modality=`` dispatch introduced for ``estimate_step_peak``:
RNA keeps its calibrated 0-12 formulas, atac/spatial/bulk share the
modality-agnostic step-0 estimator, and unknown (modality, step) pairs
return 0.0 with a warning instead of silently reusing the RNA formula
(step numbers collide across modalities).
"""

from __future__ import annotations

import logging

import pytest

from core.utils._memory import estimate_step_peak

# Real measured anchor (Li2026_Multiome, step-0 run): 70,477 cells,
# 36,601 genes, nnz 199,155,812. Step-0 formula is nnz*12*cf/1e9 +
# 1.5 + n_cells*536/1e9 + 1.0 (GB) -> ~4.93 GiB at concat_factor 1.0.
_ANCHOR_CELLS = 70_477
_ANCHOR_GENES = 36_601
_ANCHOR_NNZ = 199_155_812


def _step0_expected() -> float:
    raw_csr = _ANCHOR_NNZ * 12 * 1.0 / 1e9
    return raw_csr + 1.5 + _ANCHOR_CELLS * 536 / 1e9 + 1.0


class TestModalityStep0Shared:
    """atac/spatial/bulk step 0 == rna step 0 (shared formula)."""

    @pytest.mark.parametrize("modality", ["atac", "spatial", "bulk"])
    def test_step0_equals_rna(self, modality: str) -> None:
        rna = estimate_step_peak(0, _ANCHOR_CELLS, _ANCHOR_GENES, _ANCHOR_NNZ, concat_factor=1.0)
        other = estimate_step_peak(
            0,
            _ANCHOR_CELLS,
            _ANCHOR_GENES,
            _ANCHOR_NNZ,
            concat_factor=1.0,
            modality=modality,
        )
        assert other == pytest.approx(rna, abs=1e-9)
        assert other == pytest.approx(_step0_expected(), abs=1e-9)

    def test_step0_positive(self) -> None:
        v = estimate_step_peak(0, _ANCHOR_CELLS, _ANCHOR_GENES, _ANCHOR_NNZ, modality="atac")
        assert v > 0.0


class TestAtacStepEstimators:
    """ATAC-specific estimators registered in batch 3 (steps 1/2/4)."""

    @pytest.mark.parametrize("step", [1, 2, 4])
    def test_atac_step_positive(self, step: int) -> None:
        v = estimate_step_peak(step, 100_000, 50_000, modality="atac")
        assert v > 0.0

    def test_atac_step1_grows_with_cells(self) -> None:
        small = estimate_step_peak(1, 10_000, 50_000, modality="atac")
        large = estimate_step_peak(1, 1_000_000, 50_000, modality="atac")
        assert large > small

    def test_atac_step2_dense_float64_dominates(self) -> None:
        # 100k x 50k float64 X = 40 GB + obsm + KNN
        v = estimate_step_peak(2, 100_000, 50_000, modality="atac")
        assert v > 40.0

    def test_atac_step4_larger_than_step1(self) -> None:
        # per-cluster peak union ~1.5x pooled -> step 4 > step 1 at same size
        v1 = estimate_step_peak(1, 100_000, 0, modality="atac")
        v4 = estimate_step_peak(4, 100_000, 0, modality="atac")
        assert v4 >= v1


class TestModalityUnknownStep:
    """Unknown (modality, step) pairs -> 0.0 + warning, never RNA reuse."""

    @pytest.mark.parametrize(
        ("modality", "step"),
        [("atac", 3), ("atac", 13), ("spatial", 5), ("bulk", 3), ("metis", 0)],
    )
    def test_unknown_returns_zero_with_warning(self, modality: str, step: int, caplog) -> None:
        with caplog.at_level(logging.WARNING, logger="core.utils._memory"):
            v = (
                estimate_step_peak(0, 100, 100, nnz=0, modality=modality)
                if step == 0
                else (estimate_step_peak(step, 100, 100, nnz=0, modality=modality))
            )
        assert v == 0.0
        assert any(
            f"modality={modality}" in r.message and f"step={step}" in r.message
            for r in caplog.records
        )


class TestRnaUnchanged:
    """Default modality="rna" keeps calibrated behavior."""

    def test_rna_step4_anchor(self) -> None:
        v = estimate_step_peak(4, 110_000, 4000)
        assert v == pytest.approx(12.0, abs=1e-9)

    def test_rna_step0_anchor(self) -> None:
        v = estimate_step_peak(0, _ANCHOR_CELLS, _ANCHOR_GENES, _ANCHOR_NNZ, concat_factor=1.0)
        assert v == pytest.approx(_step0_expected(), abs=1e-9)

    @pytest.mark.parametrize("step", list(range(13)))
    def test_all_rna_steps_nonnegative(self, step: int) -> None:
        v = estimate_step_peak(step, 70_477, 36_601, 199_155_812)
        assert v >= 0.0

    def test_rna_unknown_step_zero(self) -> None:
        assert estimate_step_peak(99, 100, 100) == 0.0
