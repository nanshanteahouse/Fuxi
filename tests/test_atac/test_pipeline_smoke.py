"""ATAC pipeline smoke test — synthetic fragments through steps 00→05 (+09).

Runs the real pipeline steps as subprocesses (``core/run_pipeline.py
--modality atac --step N --config ...``) on a synthetic fragments.tsv.gz,
then asserts checkpoint contents:

- 04_peaks.h5ad carries X_umap / X_spectral + merged ``leiden`` obs (batch 1
  regression: the per-cluster peak matrix is the downstream superset).
- 05_annotated.h5ad carries ``cell_type`` (AI fallback Cluster_N when no AI).
- 09_trajectory.h5ad with no terminal_cell_types marks
  ``uns["trajectory"]["status"] == "skipped"`` and does NOT fabricate a
  zero pseudotime column (anti-fake-data regression).
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from core.utils import resolve_config  # noqa: E402

# Signal regions per cluster on chr1 (synthetic "open chromatin").
_CLUSTER_REGIONS = [
    [(2_000_000, 2_001_000), (5_000_000, 5_001_000)],
    [(8_000_000, 8_001_000), (12_000_000, 12_001_000)],
    [(15_000_000, 15_001_000), (20_000_000, 20_001_000)],
]


def _make_fragments(path: Path, n_cells: int = 200, seed: int = 7) -> Path:
    """Synthetic fragments.tsv.gz: 3 clusters, 70% signal in cluster regions."""
    per = n_cells // 3
    recs: list[tuple[str, int, int]] = []
    for cl in range(3):
        for i in range(per):
            bc = f"cell_{cl}_{i}"
            rng = np.random.RandomState(seed + cl * 100 + i)
            fr: list[tuple[int, int]] = []
            for _ in range(rng.randint(60, 220)):
                if rng.rand() < 0.70:
                    s, e = _CLUSTER_REGIONS[cl][rng.randint(0, 2)]
                else:
                    s = rng.randint(0, 25_000_000)
                    e = min(s + rng.randint(100, 500), 25_000_000)
                fr.append((s, e))
            fr.sort()
            for s, e in fr:
                recs.append((bc, s, e))
    recs.sort(key=lambda r: (r[0], r[1]))
    with gzip.open(path, "wt") as f:
        for bc, s, e in recs:
            f.write(f"chr1\t{s}\t{e}\t{bc}\t1\n")
    return path


def _write_config(tmp: Path, frag: Path) -> Path:
    """Minimal ATAC config for the synthetic run (small grids, lenient QC)."""
    # chrom_sizes schema field is a str path — write a size table file.
    cs_path = tmp / "chrom_sizes.tsv"
    cs_path.write_text("chr1\t25000000\n", encoding="utf-8")
    cfg_path = tmp / "config_atac.yaml"
    cfg_path.write_text(
        f"""modality: atac
data_format: 10x_fragments
project_dir: "{tmp}"
data_dir: "{tmp}"

data_input:
  fragment_file: "{frag}"

atac:
  genome: hg38
  chrom_sizes: "{cs_path}"
  min_fragments: 10
  max_fragments: 100000
  min_tsse: 0.0
  peak_qval: 0.05
  n_features: 2000
  n_spectral: 10

clustering:
  n_neighbors: 10
  param_grid_n_neighbors: [10]
  param_grid_resolutions: [0.5]

execution:
  n_jobs: 2
  random_seed: 42
  memory:
    policy: speed
    budget: 8GB
    guard: "off"

tissue: test
species: human
""",
        encoding="utf-8",
    )
    return cfg_path


def _run_step(step: int, cfg_path: Path) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    return subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "core" / "run_pipeline.py"),
            "--modality",
            "atac",
            "--step",
            str(step),
            "--config",
            str(cfg_path),
        ],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )


@pytest.fixture(scope="module")
def atac_pipeline(tmp_path_factory):
    """Run 00→05 once per module; return (cfg, tmp)."""
    tmp = tmp_path_factory.mktemp("atac_smoke")
    frag = _make_fragments(tmp / "fragments.tsv.gz")
    cfg_path = _write_config(tmp, frag)
    for step in range(6):  # 00 → 05
        r = _run_step(step, cfg_path)
        assert r.returncode == 0, f"step {step} failed\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    return cfg_path, tmp


@pytest.mark.slow
def test_pipeline_chain_completes(atac_pipeline):
    cfg_path, tmp = atac_pipeline
    cfg = resolve_config(str(cfg_path))
    assert os.path.exists(cfg.raw_h5ad)
    assert os.path.exists(cfg.filtered_h5ad)
    assert os.path.exists(cfg.processed_h5ad)
    assert os.path.exists(cfg.clustered_h5ad)
    assert os.path.exists(cfg.peak_h5ad)
    assert os.path.exists(cfg.annotated_h5ad)


@pytest.mark.slow
def test_04_peaks_carries_metadata_superset(atac_pipeline):
    """Batch-1 regression: 04_peaks merges clustered obs + obsm."""
    cfg_path, _ = atac_pipeline
    cfg = resolve_config(str(cfg_path))
    import snapatac2 as snap

    pd = snap.read(cfg.peak_h5ad)
    assert "X_spectral" in pd.obsm, "04_peaks missing X_spectral"
    assert "X_umap" in pd.obsm, "04_peaks missing X_umap"
    assert "leiden" in pd.obs, "04_peaks missing merged leiden obs"
    n_unassigned = sum(1 for v in pd.obs["leiden"] if str(v) == "unassigned")
    assert pd.n_obs > 0
    # at least some cells matched upstream clusters
    assert n_unassigned < pd.n_obs


@pytest.mark.slow
def test_05_annotated_has_cell_type(atac_pipeline):
    cfg_path, _ = atac_pipeline
    cfg = resolve_config(str(cfg_path))
    import snapatac2 as snap

    ad = snap.read(cfg.annotated_h5ad)
    if hasattr(ad.obs, "to_pandas"):
        obs = ad.obs.to_pandas()
    else:
        obs = ad.obs
    assert "cell_type" in obs, "05_annotated missing cell_type"
    vals = list(obs["cell_type"])
    assert all(v is not None and str(v) != "nan" for v in vals)


@pytest.mark.slow
def test_09_trajectory_skips_without_terminal_cell_types(atac_pipeline):
    """Anti-fake-data: no terminal_cell_types -> skipped, no zero pseudotime."""
    cfg_path, tmp = atac_pipeline
    cfg = resolve_config(str(cfg_path))
    cfg.atac.terminal_cell_types = []  # ensure explicit
    r = _run_step(9, cfg_path)
    assert r.returncode == 0, r.stderr
    import snapatac2 as snap

    tr = snap.read(cfg.trajectory_h5ad)
    assert tr.uns["trajectory"]["status"] == "skipped", tr.uns.get("trajectory")
    if hasattr(tr.obs, "to_pandas"):
        obs = tr.obs.to_pandas()
    else:
        obs = tr.obs
    assert "pseudotime" not in obs or obs["pseudotime"].isna().all(), (
        "trajectory fabricated pseudotime on skip path"
    )
