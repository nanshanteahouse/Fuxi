"""Spatial pipeline smoke test — synthetic Visium-like h5ad through steps 00→07.

Runs the real pipeline steps as subprocesses (``core/run_pipeline.py
--modality spatial --step N --config ...``) on a synthetic ~150-spot × ~800-gene
h5ad carrying ``obsm['spatial']`` + a 2-level ``sample`` column, then asserts
checkpoint contents:

- 00_raw.h5ad carries spatial coords + sample column
- 01_qc.h5ad survives lenient QC (no spot loss)
- 03_processed.h5ad carries the spatial neighbor graph + PCA
- 04_clustered.h5ad carries ``leiden`` + ``spatial_domain`` (leiden_spatial) + X_umap
- 05_annotated.h5ad carries ``cell_type`` (score_genes fallback → leiden labels,
  no tissue_kb / no AI configured)
- 06_deconvolved.h5ad records ``uns["deconvolution"]["status"] == "skipped"``
  (deconv_method='none' must still produce the checkpoint)
- 06_svg.h5ad carries full-gene Moran's I results (``uns["moranI"]`` + ``uns["svg"]``)

Skip with ``SKIP_SLOW_TESTS=1``; opt-in collection via ``-m slow``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from core.utils import resolve_config  # noqa: E402

N_SPOTS = 150
N_GENES = 800
N_DOMAINS = 3
_STEPS = 8  # 00 → 07


def _skip_slow() -> bool:
    """Return True when SKIP_SLOW_TESTS env var is set to skip expensive tests."""
    return os.environ.get("SKIP_SLOW_TESTS", "").lower() in ("1", "true")


# ═════════════════════════════════════════════════════════════════════════
#  Synthetic data builder
# ═════════════════════════════════════════════════════════════════════════


def _make_synthetic_h5ad(path: Path, *, seed: int = 42) -> Path:
    """Synthetic Visium-like sparse-count h5ad with 3 marker-gene domains.

    Each spot draws from one of ``N_DOMAINS`` marker-gene blocks (higher Poisson
    rate) plus a random background; ``obs`` carries a 2-level ``sample`` column
    (simulates a merged multi-slide object) and ``in_tissue=1``; coordinates are
    stored in ``obsm['spatial']``.  Deterministic via ``np.random.RandomState``.
    """
    rng = np.random.RandomState(seed)
    markers_per = N_GENES // N_DOMAINS  # first block of genes = domain markers
    marker_sets = [set(range(d * markers_per, (d + 1) * markers_per)) for d in range(N_DOMAINS)]

    x = sp.lil_matrix((N_SPOTS, N_GENES), dtype=np.float32)
    for i in range(N_SPOTS):
        markers = sorted(marker_sets[rng.randint(0, N_DOMAINS)])
        x[i, markers] = rng.poisson(5.0, size=len(markers))
        bg = rng.choice(N_GENES, size=rng.randint(20, 60), replace=False)
        x[i, bg] = rng.poisson(0.5, size=len(bg))
    x = x.tocsr()

    adata = sc.AnnData(X=x)
    adata.var_names = [f"GENE_{i:04d}" for i in range(N_GENES)]
    adata.obs_names = [f"spot_{i:04d}" for i in range(N_SPOTS)]
    adata.obs["sample"] = [f"sample_{i % 2}" for i in range(N_SPOTS)]
    adata.obs["in_tissue"] = 1
    adata.obsm["spatial"] = np.column_stack(
        [rng.randint(0, 15, size=N_SPOTS), rng.randint(0, 10, size=N_SPOTS)]
    ).astype(np.float64)
    adata.write(path, compression="gzip")
    return path


# ═════════════════════════════════════════════════════════════════════════
#  Config / subprocess helpers
# ═════════════════════════════════════════════════════════════════════════


def _write_config(tmp: Path, h5ad_path: Path) -> Path:
    """Minimal spatial config for the synthetic run (lenient QC, off switches)."""
    cfg_path = tmp / "config_spatial.yaml"
    cfg_path.write_text(
        f"""modality: spatial
data_format: h5ad
project_dir: "{tmp}"
data_dir: "{tmp}"
tissue: test
species: human

data_input:
  input_h5ad: "{h5ad_path}"

spatial:
  platform: generic
  crop_image: false
  run_segmentation: false
  decontamination: none
  integration_method: none
  domain_method: leiden_spatial
  deconv_method: none
  svg_method: moran
  run_autocorr: true

qc:
  min_genes: 10
  max_genes: 10000
  max_pct_mito: 50.0
  min_genes_per_umi: 0.0
  min_cells_per_gene: 2
  use_adaptive_thresholds: false

hvg:
  n_top_genes: 200
  flavor: seurat

pca:
  n_pcs_use: 20

normalization:
  normalize_target_sum: 10000

clustering:
  cluster_selection_method: silhouette
  param_grid_n_neighbors: [10]
  param_grid_resolutions: [0.5]
  param_grid_min_dist: [0.3]
  param_grid_spread: [1.0]
  umap_min_dist: 0.3
  umap_spread: 1.0
  leiden_flavor: igraph
  leiden_n_iterations: 2
  stability_n_seeds: 3
  multi_metric_stability_top_k: 2

de:
  method: wilcoxon
  n_genes: 50

marker:
  marker_dict: {{}}

ai:
  enabled: false

execution:
  n_jobs: 2
  random_seed: 42
  device: cpu
  memory:
    guard: "off"
""",
        encoding="utf-8",
    )
    return cfg_path


def _run_step(step: int, cfg_path: Path, data_root: Path) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    env["FUXI_DATA_ROOT"] = str(data_root)
    return subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "core" / "run_pipeline.py"),
            "--modality",
            "spatial",
            "--step",
            str(step),
            "--config",
            str(cfg_path),
        ],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
    )


@pytest.fixture(scope="module")
def spatial_pipeline(tmp_path_factory):
    """Run steps 00→07 once per module; return (cfg_path, tmp)."""
    tmp = tmp_path_factory.mktemp("spatial_smoke")
    h5ad_path = _make_synthetic_h5ad(tmp / "synthetic_spatial.h5ad")
    cfg_path = _write_config(tmp, h5ad_path)
    for step in range(_STEPS):  # 00 → 07
        r = _run_step(step, cfg_path, tmp)
        assert r.returncode == 0, f"step {step} failed\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    return cfg_path, tmp


# ═════════════════════════════════════════════════════════════════════════
#  Tests
# ═════════════════════════════════════════════════════════════════════════


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_pipeline_chain_completes(spatial_pipeline):
    """Every step 00→07 produced its checkpoint file."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    for name in (
        "00_raw.h5ad",
        "01_qc.h5ad",
        "02_image.h5ad",
        "03_processed.h5ad",
        "04_clustered.h5ad",
        "05_annotated.h5ad",
        "06_deconvolved.h5ad",
        "06_svg.h5ad",
    ):
        assert os.path.exists(os.path.join(cfg.h5ad_dir, name)), f"missing checkpoint {name}"


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_00_raw_carries_spatial_coords_and_sample(spatial_pipeline):
    """00_load (h5ad branch) preserves coords + the 2-level sample column."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "00_raw.h5ad"))
    assert ad.shape == (N_SPOTS, N_GENES)
    assert "spatial" in ad.obsm and ad.obsm["spatial"].shape == (N_SPOTS, 2)
    assert "sample" in ad.obs
    assert ad.obs["sample"].nunique() == 2
    assert sp.isspmatrix_csr(ad.X)


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_01_qc_keeps_spots(spatial_pipeline):
    """Lenient QC thresholds must not empty the object."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "01_qc.h5ad"))
    assert ad.n_obs == N_SPOTS, "lenient QC should keep all synthetic spots"
    assert "spatial" in ad.obsm


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_03_processed_has_graph_and_pca(spatial_pipeline):
    """03_normalize builds the spatial neighbor graph + PCA on HVGs."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "03_processed.h5ad"))
    assert "spatial_connectivities" in ad.obsp, "missing spatial neighbor graph"
    assert "X_pca" in ad.obsm and ad.obsm["X_pca"].shape[1] == 20
    assert "highly_variable" in ad.var and ad.var["highly_variable"].sum() > 0


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_04_clustered_carries_leiden_and_spatial_domain(spatial_pipeline):
    """04_cluster auto-locks leiden + X_umap and adds the spatial-domain label."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "04_clustered.h5ad"))
    assert "leiden" in ad.obs and ad.obs["leiden"].nunique() >= 2
    assert "X_umap" in ad.obsm and ad.obsm["X_umap"].shape == (N_SPOTS, 2)
    assert "spatial_domain" in ad.obs, "leiden_spatial must add the domain column"
    assert ad.uns["spatial_domain"]["method"] == "leiden_spatial"
    assert "cell_type" not in ad.obs, "annotation must be deferred to step 05"


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_05_annotated_has_cell_type(spatial_pipeline):
    """score_genes fallback (no tissue_kb / AI) yields non-null cell_type."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "05_annotated.h5ad"))
    assert "cell_type" in ad.obs, "05_annotated missing cell_type"
    vals = ad.obs["cell_type"].astype(str).values
    assert all(v != "nan" for v in vals)
    assert ad.obs["cell_type"].nunique() >= 2


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_06_deconvolved_graceful_skip(spatial_pipeline):
    """deconv_method='none' → skipped status recorded, checkpoint still written."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "06_deconvolved.h5ad"))
    assert ad.uns["deconvolution"]["status"] == "skipped", ad.uns.get("deconvolution")
    assert ad.obs["cell_type"].nunique() >= 2  # upstream state preserved


@pytest.mark.slow
@pytest.mark.skipif(_skip_slow(), reason="Slow pipeline smoke test skipped via SKIP_SLOW_TESTS")
def test_07_svg_full_gene_autocorr(spatial_pipeline):
    """07_spatial_stats runs Moran's I on the full gene set and persists 06_svg.h5ad."""
    cfg_path, _ = spatial_pipeline
    cfg = resolve_config(str(cfg_path))
    ad = sc.read(os.path.join(cfg.h5ad_dir, "06_svg.h5ad"))
    moran = ad.uns["moranI"]
    n_full = ad.raw.n_vars if ad.raw is not None else ad.n_vars
    assert moran.shape[0] == n_full, "Moran's I must screen the FULL gene set"
    assert ad.uns["svg"]["method"] == "moran"
    assert ad.uns["svg"]["n_genes_tested"] == n_full
    # top SVGs marked on the HVG-subset var
    assert "spatially_variable" in ad.var and ad.var["spatially_variable"].sum() > 0
