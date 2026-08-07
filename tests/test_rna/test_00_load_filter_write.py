"""Filter/downsample-on-write tests for rna/steps/00_load.py (metis G4).

The step-00 final write filters and downsamples cells ON WRITE via a boolean
cell mask (no full-matrix copies). These tests pin the mask computation to
the UNMODIFIED ``core.downsample`` reference oracle: for every strategy /
filter combo the selected cell order must be RNG-identical to the old
``filter_by_config`` + ``downsample_by_config`` chain, and the written h5ad
must round-trip to the old-path output.

Coverage
--------
- reference oracle: stratified (seeded / non-default seed), random,
  max_per_sample, sample_keep, obs_filter, filter+downsample combos, no-op
  passthrough, all-cells-filtered edge
- write path: masked ``write_csr_chunked`` output re-read equals the old-path
  ``safe_write`` output (obs_names order + X content + var_names)
- raw-write compression (folded-in T7 wiring): ``_raw_h5ad_compression``
  resolves ``"zstd"`` by schema default; the masked writer stores the zstd
  filter (32015) when hdf5plugin is present, gzip otherwise
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata
import h5py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy.sparse import csr_matrix

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_filter_write_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

from core.downsample import downsample_by_config, filter_by_config  # noqa: E402
from core.utils import safe_write, write_csr_chunked  # noqa: E402

# The reference oracle (core/downsample.py L66) triggers a pandas FutureWarning
# on `per_sample_targets.iloc[i] += 1`; this repo turns warnings into errors.
pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning")


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════


def _make_cfg(
    target=None,
    strategy="stratified",
    random_seed=42,
    max_per_sample=5000,
    sample_keep=None,
    obs_filter="",
    per_step=None,
):
    """Duck-typed config carrying the same attributes 00_load reads.

    ``filter_by_config`` reads nested ``cfg.downsample.sample_keep/obs_filter``;
    ``downsample_by_config`` reads the FLAT ``cfg.downsample_target`` etc. —
    exactly the shape a resolved ``core.config.schema.Config`` exposes to the
    step-00 call site.
    """
    return SimpleNamespace(
        downsample=SimpleNamespace(sample_keep=sample_keep or [], obs_filter=obs_filter),
        downsample_target=target,
        downsample_strategy=strategy,
        downsample_random_seed=random_seed,
        downsample_max_per_sample=max_per_sample,
        per_step_h5ad_compression={"raw": "zstd"} if per_step is None else dict(per_step),
        h5ad_compression="gzip",
        h5ad_compression_opts=None,
    )


def _build_adata(n=100, seed=3) -> anndata.AnnData:
    rng = np.random.RandomState(seed)
    x = csr_matrix(rng.rand(n, 10) * 5).astype(np.float32)
    samples = np.random.RandomState(seed + 1).choice(["S1", "S2", "S3"], n)
    obs = pd.DataFrame({"sample": samples}, index=[f"cell_{i:04d}" for i in range(n)])
    var = pd.DataFrame(index=[f"GENE_{i}" for i in range(10)])
    return anndata.AnnData(X=x, obs=obs, var=var)


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_load_filter_write")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


def _oracle(adata, cfg) -> anndata.AnnData:
    """Reference oracle: the UNMODIFIED core.downsample chain (old path)."""
    return downsample_by_config(filter_by_config(adata, cfg, _make_logger()), cfg, _make_logger())


# ═════════════════════════════════════════════════════════════════════
#  Reference-oracle: mask selection == old chain (RNG-exact)
# ═════════════════════════════════════════════════════════════════════

_ORACLE_CASES = [
    ("noop", {}),
    ("sample_keep", {"sample_keep": ["S1", "S2"]}),
    ("obs_filter", {"obs_filter": "sample != 'S3'"}),
    (
        "obs_filter_eval_fallback",
        {"obs_filter": "True"},
    ),  # obs.eval returns scalar → per-row eval fallback
    ("stratified", {"target": 40}),
    ("stratified_seed7", {"target": 40, "random_seed": 7}),
    ("stratified_filter_combo", {"target": 25, "sample_keep": ["S1", "S2"]}),
    ("stratified_obs_filter", {"target": 30, "obs_filter": "sample != 'S3'"}),
    ("random", {"target": 50, "strategy": "random"}),
    ("max_per_sample", {"target": 10_000, "strategy": "max_per_sample", "max_per_sample": 20}),
    ("target_ge_n", {"target": 500}),
    ("obs_filter_false_all", {"obs_filter": "sample == 'NOPE'"}),
]


class TestMaskOracle:
    @pytest.mark.parametrize(
        "params", [c[1] for c in _ORACLE_CASES], ids=[c[0] for c in _ORACLE_CASES]
    )
    def test_selection_matches_oracle(self, params: dict) -> None:
        adata = _build_adata()
        cfg = _make_cfg(**params)
        mask = _mod._filter_downsample_mask(adata.obs, cfg, _make_logger())
        oracle_adata = _oracle(adata, cfg)
        # SAME cells in SAME order (RNG-exact) — obs_names order + set.
        assert list(adata.obs_names[mask]) == list(oracle_adata.obs_names)
        assert int(mask.sum()) == oracle_adata.n_obs

    def test_x_content_on_subset_matches_oracle(self) -> None:
        """The X rows selected by the mask equal the oracle subset's X."""
        adata = _build_adata()
        cfg = _make_cfg(target=40, sample_keep=["S1", "S2"])
        mask = _mod._filter_downsample_mask(adata.obs, cfg, _make_logger())
        oracle_adata = _oracle(adata, cfg)
        assert np.array_equal(adata.X[mask].toarray(), oracle_adata.X.toarray(), equal_nan=True)


# ═════════════════════════════════════════════════════════════════════
#  Write path: masked chunked write == old safe_write output
# ═════════════════════════════════════════════════════════════════════


class TestWritePath:
    @pytest.mark.parametrize(
        "params",
        [
            {"target": 40},
            {"sample_keep": ["S1", "S2"]},
            {"target": 10_000, "strategy": "max_per_sample", "max_per_sample": 25},
        ],
        ids=["stratified", "sample_keep", "max_per_sample"],
    )
    def test_masked_write_roundtrips_to_old_path(self, tmp_path: Path, params: dict) -> None:
        adata = _build_adata()
        cfg = _make_cfg(**params)
        mask = _mod._filter_downsample_mask(adata.obs, cfg, _make_logger())

        # New path: filter ON WRITE via the cell mask (same code shape as main()).
        compression, compression_opts = _mod._raw_h5ad_compression(cfg)
        new_path = str(tmp_path / "new.h5ad")
        write_csr_chunked(
            new_path,
            adata.X,
            cell_mask=mask,
            n_obs=int(mask.sum()),
            obs=adata.obs[mask].copy(),
            var=adata.var.copy(),
            compression=compression,
            compression_opts=compression_opts,
        )

        # Old path: filter + downsample into AnnData copies, then safe_write.
        old_path = str(tmp_path / "old.h5ad")
        safe_write(_oracle(adata, cfg), old_path)

        new = sc.read(new_path)
        old = sc.read(old_path)
        assert list(new.obs_names) == list(old.obs_names)
        assert list(new.var_names) == list(old.var_names)
        assert np.array_equal(new.X.toarray(), old.X.toarray(), equal_nan=True)

    def test_noop_write_identical_to_safe_write(self, tmp_path: Path) -> None:
        """No filter + no downsample → mask all-True; content identical to old."""
        adata = _build_adata()
        cfg = _make_cfg()
        mask = _mod._filter_downsample_mask(adata.obs, cfg, _make_logger())
        assert mask.all()

        compression, compression_opts = _mod._raw_h5ad_compression(cfg)
        new_path = str(tmp_path / "new.h5ad")
        write_csr_chunked(
            new_path,
            adata.X,
            cell_mask=mask,
            n_obs=int(mask.sum()),
            obs=adata.obs.copy(),
            var=adata.var.copy(),
            compression=compression,
            compression_opts=compression_opts,
        )
        old_path = str(tmp_path / "old.h5ad")
        safe_write(adata, old_path)

        new = sc.read(new_path)
        old = sc.read(old_path)
        assert list(new.obs_names) == list(old.obs_names)
        assert list(new.var_names) == list(old.var_names)
        assert np.array_equal(new.X.toarray(), old.X.toarray(), equal_nan=True)


# ═════════════════════════════════════════════════════════════════════
#  Raw-write compression (folded-in T7 wiring: step_alias="raw" → zstd)
# ═════════════════════════════════════════════════════════════════════


class TestRawCompression:
    def test_raw_compression_resolves_zstd_by_default(self) -> None:
        """_raw_h5ad_compression picks per_step["raw"]="zstd" over gzip."""
        compression, compression_opts = _mod._raw_h5ad_compression(_make_cfg())
        assert compression == "zstd"
        assert compression_opts is None

    def test_raw_compression_falls_back_to_global(self) -> None:
        """per_step dict without "raw" key → cfg.h5ad_compression."""
        cfg = _make_cfg(per_step={"integrated": "gzip"})
        compression, compression_opts = _mod._raw_h5ad_compression(cfg)
        assert compression == "gzip"
        assert compression_opts is None

    def test_masked_write_stores_zstd_filter(self, tmp_path: Path) -> None:
        """Masked writer with zstd → X/data filter 32015 (hdf5plugin) / gzip."""
        adata = _build_adata()
        cfg = _make_cfg(target=40)
        mask = _mod._filter_downsample_mask(adata.obs, cfg, _make_logger())
        compression, compression_opts = _mod._raw_h5ad_compression(cfg)
        path = str(tmp_path / "raw.h5ad")
        write_csr_chunked(
            path,
            adata.X,
            cell_mask=mask,
            n_obs=int(mask.sum()),
            obs=adata.obs[mask].copy(),
            var=adata.var.copy(),
            compression=compression,
            compression_opts=compression_opts,
        )
        with h5py.File(path, "r") as f:
            dset = f["X/data"]
            try:
                has_zstd = "32015" in dset._filters
            except AttributeError:
                has_zstd = False
        if compression == "zstd":
            assert has_zstd, f"expected zstd filter 32015, got {dset.compression}"
        else:  # hdf5plugin absent → gzip fallback
            assert not has_zstd
        # Round-trip still reads back.
        reread = sc.read(path)
        assert reread.n_obs == int(mask.sum())
