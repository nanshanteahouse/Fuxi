"""Tests for validate_pipeline_state — cross-step data flow validation.

Checks that required obs/obsm/obsp columns are enforced per step boundary.
"""

import numpy as np
import pytest
from anndata import AnnData

from core.utils import _STEP_REQUIREMENTS, validate_pipeline_state

# ── Helpers ────────────────────────────────────────────────────────────


def _make_adata(obs_cols=None, obsm_keys=None, obsp_keys=None):
    """Create a minimal AnnData with optional columns."""
    adata = AnnData(
        X=np.zeros((10, 5)),
        obs={k: ["a"] * 10 for k in (obs_cols or [])},
    )
    for k in obsm_keys or []:
        adata.obsm[k] = np.zeros((10, 2))
    for k in obsp_keys or []:
        adata.obsp[k] = np.eye(10)
    return adata


# ── RNA tests ──────────────────────────────────────────────────────────


class TestRNAStepValidation:
    def test_step01_missing_doublet_scores(self):
        """Step 01 requires doublet_scores and predicted_doublet."""
        adata = _make_adata()
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="01", modality="rna")
        assert "obs['doublet_scores']" in str(exc.value)
        assert "obs['predicted_doublet']" in str(exc.value)

    def test_step01_passes(self):
        """Step 01 passes with both required columns."""
        adata = _make_adata(
            obs_cols=["doublet_scores", "predicted_doublet"],
        )
        validate_pipeline_state(adata, step="01", modality="rna")  # no raise

    def test_step02_missing_predicted_doublet(self):
        """Step 02 requires doublet_scores and predicted_doublet."""
        adata = _make_adata(obs_cols=["doublet_scores"])
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="02", modality="rna")
        assert "obs['predicted_doublet']" in str(exc.value)

    def test_step04_missing_X_pca(self):
        """Step 04 requires obsm['X_pca']."""
        adata = _make_adata()
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="04", modality="rna")
        assert "obsm['X_pca']" in str(exc.value)

    def test_step04_passes(self):
        """Step 04 passes with X_pca in obsm."""
        adata = _make_adata(obsm_keys=["X_pca"])
        validate_pipeline_state(adata, step="04", modality="rna")  # no raise

    def test_step05_missing_leiden(self):
        """Step 05 requires leiden, X_umap, and X_pca."""
        adata = _make_adata(obsm_keys=["X_umap", "X_pca"])
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="05", modality="rna")
        assert "obs['leiden']" in str(exc.value)

    def test_step05_passes(self):
        """Step 05 passes with all required columns."""
        adata = _make_adata(
            obs_cols=["leiden"],
            obsm_keys=["X_umap", "X_pca"],
        )
        validate_pipeline_state(adata, step="05", modality="rna")  # no raise

    def test_step06_missing_cell_type(self):
        """Step 06 requires cell_type and X_umap."""
        adata = _make_adata(obsm_keys=["X_umap"])
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="06", modality="rna")
        assert "obs['cell_type']" in str(exc.value)

    def test_step06_passes(self):
        """Step 06 passes with cell_type and X_umap."""
        adata = _make_adata(
            obs_cols=["cell_type"],
            obsm_keys=["X_umap"],
        )
        validate_pipeline_state(adata, step="06", modality="rna")  # no raise

    def test_step07_requires_cell_type_and_leiden(self):
        """Step 07 needs both cell_type and leiden."""
        adata = _make_adata(obs_cols=["cell_type"])
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="07", modality="rna")
        assert "obs['leiden']" in str(exc.value)

    def test_step08_requires_cell_type_xpca_xumap(self):
        """Step 08 needs cell_type and PCA/UMAP embeddings."""
        adata = _make_adata(obs_cols=["cell_type"])
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="08", modality="rna")
        assert "obsm['X_pca']" in str(exc.value)
        assert "obsm['X_umap']" in str(exc.value)


# ── ATAC tests ─────────────────────────────────────────────────────────


class TestATACStepValidation:
    def test_atac_step02_requires_predicted_doublet(self):
        """ATAC step 02 requires predicted_doublet."""
        adata = _make_adata()
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="02", modality="atac")
        assert "obs['predicted_doublet']" in str(exc.value)

    def test_atac_step04_requires_leiden_and_umap(self):
        """ATAC step 04 requires leiden and X_umap."""
        adata = _make_adata()
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="04", modality="atac")
        assert "obs['leiden']" in str(exc.value)
        assert "obsm['X_umap']" in str(exc.value)

    def test_atac_step05_missing_cell_type(self):
        """ATAC step 05 requires cell_type and X_umap."""
        adata = _make_adata()
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="05", modality="atac")
        assert "obs['cell_type']" in str(exc.value)

    def test_atac_step09_passes(self):
        """ATAC step 09 passes with cell_type."""
        adata = _make_adata(obs_cols=["cell_type"])
        validate_pipeline_state(adata, step="09", modality="atac")  # no raise


# ── Spatial tests ──────────────────────────────────────────────────────


class TestSpatialStepValidation:
    def test_spatial_step04_missing_xpca(self):
        """Spatial step 04 requires X_pca."""
        adata = _make_adata()
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="04", modality="spatial")
        assert "obsm['X_pca']" in str(exc.value)

    def test_spatial_step06_requires_spatial_connectivities(self):
        """Spatial step 06 requires cell_type and spatial_connectivities."""
        adata = _make_adata(obs_cols=["cell_type"])
        with pytest.raises(AssertionError) as exc:
            validate_pipeline_state(adata, step="06", modality="spatial")
        assert "obsp['spatial_connectivities']" in str(exc.value)

    def test_spatial_step06_passes(self):
        """Spatial step 06 passes with cell_type and spatial_connectivities."""
        adata = _make_adata(
            obs_cols=["cell_type"],
            obsp_keys=["spatial_connectivities"],
        )
        validate_pipeline_state(adata, step="06", modality="spatial")  # no raise


# ── Edge cases ─────────────────────────────────────────────────────────


class TestEdgeCases:
    def test_step00_skips_validation(self):
        """Step 00 (no requirements) should not raise."""
        adata = _make_adata()
        validate_pipeline_state(adata, step="00", modality="rna")  # no raise

    def test_step00_atac_skips_validation(self):
        """ATAC step 00 (no requirements) should not raise."""
        adata = _make_adata()
        validate_pipeline_state(adata, step="00", modality="atac")  # no raise

    def test_unknown_modality_skips_gracefully(self):
        """Unknown modality with no requirements logs debug, does not raise."""
        adata = _make_adata()
        validate_pipeline_state(adata, step="05", modality="unknown")  # no raise

    def test_step_03_has_no_requirements(self):
        """RNA step 03 has no defined requirements."""
        adata = _make_adata()
        validate_pipeline_state(adata, step="03", modality="rna")  # no raise

    def test_all_modalities_have_requirements(self):
        """Every modality and step references valid keys."""
        for modality, steps in _STEP_REQUIREMENTS.items():
            for step_key, req in steps.items():
                assert isinstance(req, dict), f"{modality}/{step_key} not a dict"
                for key in ("obs", "obsm", "obsp"):
                    if key in req:
                        assert isinstance(req[key], list), (
                            f"{modality}/{step_key}.{key} not a list"
                        )
