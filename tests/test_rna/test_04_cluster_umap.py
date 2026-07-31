"""Tests for rna/steps/04_cluster_umap.py — copy+append incremental h5ad write.

Item 1.7 (h5ad-incremental-io): the two ``safe_write`` sites in step 04 are
replaced by a copy+append helper ``_write_cluster_h5ad``:

- ``incremental_io=True``  → ``copy2(03_integrated.h5ad → 04_clustered.h5ad)``
  then ``write_h5ad_incremental`` with ONLY the keys step 04 owns
  (mode A: on failure the corrupt copy is deleted; 03 stays pristine).
- ``incremental_io=False`` → full ``safe_write`` fallback (escape hatch).

The OWNED-key contract is the core of the plan: everything the file inherited
from 03 (X, layers, raw, ``X_pca``/``X_integrated`` obsm, 03-era ``uns`` such
as ``pca`` variance_ratio / integration results, 03-era obs columns) must be
left untouched — only ``obsm[X_umap|umap_*]``, ``obsp[connectivities|distances]``,
step-04 ``uns`` keys and step-04 ``leiden*``/``funnel*`` obs columns are
appended/overwritten.

Run with::

    pytest tests/test_rna/test_04_cluster_umap.py -v --tb=short
"""

from __future__ import annotations

import importlib.util
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

# ── Ensure repo root is on sys.path (conftest.py also does this) ──────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 04_cluster_umap module via file path (name starts with digit) ──
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "04_cluster_umap.py")
_spec = importlib.util.spec_from_file_location("rna.steps._04_cluster_umap_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ── Fixtures / helpers ─────────────────────────────────────────────────


def _make_step04_state(n: int = 60, n_genes: int = 40, seed: int = 7) -> tuple[AnnData, AnnData]:
    """Return (step03, step04_adata) pair.

    ``step03`` is the pristine 03_integrated.h5ad content: X + X_pca +
    pca variance_ratio uns + a 03-era obs column + a 03-era uns key.
    ``step04_adata`` is the in-memory object after step 04 ran: it contains
    everything step03 has PLUS the keys step 04 owns.
    """
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=2.0, size=(n, n_genes)).astype(np.float32)

    step03 = AnnData(x)
    step03.var_names = [f"GENE_{i}" for i in range(n_genes)]
    step03.obs_names = [f"cell_{i}" for i in range(n)]
    step03.obs["sample"] = rng.choice(["A", "B"], n)
    step03.obsm["X_pca"] = rng.standard_normal((n, 8))
    step03.uns["pca"] = {"variance_ratio": rng.rand(8), "params": {"n_pcs": 8}}
    step03.uns["integration"] = {"method": "harmony"}

    adata = step03.copy()
    # ── step-04 owned obsm ──
    adata.obsm["X_umap"] = rng.standard_normal((n, 2))
    adata.obsm["umap_5_0.3"] = adata.obsm["X_umap"].copy()
    adata.obsm["umap_10_0.3"] = adata.obsm["X_umap"] * 1.5
    # ── step-04 owned obsp (sparse, the risk area) ──
    idx = rng.randint(0, n, size=2 * n)
    data = rng.rand(2 * n).astype(np.float32)
    conn = sp.csr_matrix((data, (np.arange(2 * n) % n, idx)), shape=(n, n))
    adata.obsp["connectivities"] = conn
    adata.obsp["distances"] = sp.csr_matrix(conn + conn.T)
    # ── step-04 owned uns ──
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": 5, "use_rep": "X_pca", "n_pcs": 8},
    }
    adata.uns["umap"] = {"params": {"min_dist": 0.3, "spread": 1.0}}
    adata.uns["leiden_5_0.3"] = {"params": {"resolution": 0.3, "n_clusters": 3}}
    adata.uns["leiden_10_0.3"] = {"params": {"resolution": 0.3, "n_clusters": 5}}
    adata.uns["grid_scan_mode"] = "full"
    adata.uns["best_resolution"] = 0.3
    adata.uns["best_n_neighbors"] = 5
    adata.uns["cluster_selection_method"] = "pareto_elbow"
    # ── step-04 owned obs columns ──
    adata.obs["leiden"] = pd.Categorical(rng.choice(["0", "1", "2"], n))
    adata.obs["leiden_5_0.3"] = adata.obs["leiden"].copy()
    adata.obs["leiden_10_0.3"] = pd.Categorical(rng.choice(["0", "1", "2", "3", "4"], n))
    return step03, adata


def _make_config(tmp_path, integrated_name="03_integrated.h5ad") -> MagicMock:
    cfg = MagicMock()
    cfg.integrated_h5ad = str(tmp_path / integrated_name)
    cfg.cluster_h5ad = str(tmp_path / "04_clustered.h5ad")
    return cfg


def _write_step03(cfg, step03: AnnData) -> None:
    os.makedirs(os.path.dirname(cfg.integrated_h5ad), exist_ok=True)
    step03.write_h5ad(cfg.integrated_h5ad)


# ── Tests ──────────────────────────────────────────────────────────────


class TestIncrementalWrite:
    """copy+append round-trip through the real write_h5ad_incremental engine."""

    def test_roundtrip_owned_keys_copied_and_appended(self, tmp_path) -> None:
        """Full round-trip: copy 03 → append → read back with scanpy.

        Given: 03_integrated.h5ad on disk (X + pca + 03-era obs/uns) and an
              in-memory adata holding every step-04 key.
        When:  _write_cluster_h5ad runs with incremental_io=True.
        Then:  the file is 03 + step-04 keys only — obsm values equal,
              obsp sparse structure (indices/indptr/data) identical, uns
              params equal, obs column values/dtypes equal, and no 03 key
              was rewritten (X/layers untouched).
        """
        step03, adata = _make_step04_state()
        cfg = _make_config(tmp_path)
        cfg.incremental_io = True
        _write_step03(cfg, step03)

        _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        out = sc.read(cfg.cluster_h5ad)

        # 03 keys preserved verbatim (copy2 semantics: nothing rewritten)
        assert list(out.obs.columns) == ["sample", "leiden", "leiden_5_0.3", "leiden_10_0.3"]
        np.testing.assert_array_equal(out.obsm["X_pca"], step03.obsm["X_pca"])
        np.testing.assert_allclose(
            out.uns["pca"]["variance_ratio"], step03.uns["pca"]["variance_ratio"]
        )
        assert out.uns["integration"] == {"method": "harmony"}
        assert out.obs["sample"].tolist() == step03.obs["sample"].tolist()
        # X / layers never touched by copy+append
        np.testing.assert_array_equal(np.asarray(out.X), np.asarray(step03.X))
        assert set(out.layers.keys()) == set(step03.layers.keys())

        # step-04 obsm appended with exact values
        for k in ("X_umap", "umap_5_0.3", "umap_10_0.3"):
            assert k in out.obsm, f"obsm key {k!r} missing after append"
            np.testing.assert_allclose(out.obsm[k], adata.obsm[k])

        # step-04 obsp: sparse structure must round-trip exactly
        for k in ("connectivities", "distances"):
            got = out.obsp[k]
            want = adata.obsp[k]
            assert sp.issparse(got), f"obsp {k!r} lost sparsity"
            got_csr = got.tocsr()
            want_csr = want.tocsr()
            assert got_csr.shape == want_csr.shape
            np.testing.assert_array_equal(got_csr.indices, want_csr.indices)
            np.testing.assert_array_equal(got_csr.indptr, want_csr.indptr)
            np.testing.assert_allclose(got_csr.data, want_csr.data)

        # step-04 uns appended
        assert out.uns["neighbors"]["params"]["n_neighbors"] == 5
        assert out.uns["umap"]["params"]["min_dist"] == 0.3
        assert out.uns["leiden_5_0.3"]["params"]["resolution"] == 0.3
        assert out.uns["grid_scan_mode"] == "full"
        assert out.uns["best_resolution"] == 0.3
        assert out.uns["best_n_neighbors"] == 5
        assert out.uns["cluster_selection_method"] == "pareto_elbow"

        # step-04 obs columns appended with values and dtypes
        assert isinstance(out.obs["leiden"].dtype, pd.CategoricalDtype)
        assert out.obs["leiden"].tolist() == adata.obs["leiden"].tolist()
        assert out.obs["leiden_5_0.3"].tolist() == adata.obs["leiden_5_0.3"].tolist()
        assert out.obs["leiden_10_0.3"].tolist() == adata.obs["leiden_10_0.3"].tolist()

    def test_only_step04_owned_keys_passed_to_engine(self, tmp_path) -> None:
        """The incremental call must receive ONLY step-04 keys.

        Given: adata holding 03-era keys (X_pca, pca, integration, sample)
              AND step-04 keys.
        When:  _write_cluster_h5ad runs with incremental_io=True.
        Then:  write_h5ad_incremental receives exactly {X_umap, umap_*},
              {connectivities, distances}, the step-04 uns keys and the
              step-04 obs columns — never X_pca / pca / integration / sample.
        """
        _, adata = _make_step04_state()
        cfg = _make_config(tmp_path)
        cfg.incremental_io = True
        _write_step03(cfg, _make_step04_state()[0])

        with patch("core.utils.write_h5ad_incremental") as mock_write:
            _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        assert mock_write.call_count == 1
        kwargs = mock_write.call_args.kwargs

        # obsm: step-04 only
        assert set(kwargs["obsm"]) == {"X_umap", "umap_5_0.3", "umap_10_0.3"}
        assert "X_pca" not in kwargs["obsm"]
        # obsp
        assert set(kwargs["obsp"]) == {"connectivities", "distances"}
        # uns: step-04 only — never pca / integration
        assert "pca" not in kwargs["uns"]
        assert "integration" not in kwargs["uns"]
        for k in (
            "neighbors",
            "umap",
            "leiden_5_0.3",
            "leiden_10_0.3",
            "grid_scan_mode",
            "best_resolution",
            "best_n_neighbors",
            "cluster_selection_method",
        ):
            assert k in kwargs["uns"]
        # obs: only the step-04 columns
        assert list(kwargs["obs"].columns) == ["leiden", "leiden_5_0.3", "leiden_10_0.3"]

    def test_rerun_overwrites_stale_values(self, tmp_path) -> None:
        """Re-running the step must replace stale results (overwrite semantics).

        Given: an existing 04_clustered.h5ad written from run #1.
        When:  _write_cluster_h5ad runs again with changed X_umap + leiden.
        Then:  the file reflects the NEW values, not the old ones.
        """
        step03, adata = _make_step04_state(seed=11)
        cfg = _make_config(tmp_path)
        cfg.incremental_io = True
        _write_step03(cfg, step03)
        _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        # "second run": mutate the step-04 keys and write again
        adata.obsm["X_umap"] = adata.obsm["X_umap"] * 100.0
        adata.obs["leiden"] = pd.Categorical(["7"] * adata.n_obs)
        _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        out = sc.read(cfg.cluster_h5ad)
        np.testing.assert_allclose(out.obsm["X_umap"], adata.obsm["X_umap"] * 1.0)
        assert out.obs["leiden"].tolist() == ["7"] * adata.n_obs
        # untouched step-04 keys from run #1 still present
        assert out.uns["best_n_neighbors"] == 5
        # 03 keys still intact (never blown away by either append)
        np.testing.assert_array_equal(out.obsm["X_pca"], step03.obsm["X_pca"])

    def test_failure_deletes_corrupt_copy_mode_a(self, tmp_path) -> None:
        """Mode A: a failing append deletes the corrupt 04 copy; 03 survives.

        Given: 03_integrated.h5ad on disk.
        When:  write_h5ad_incremental raises (simulated crash mid-append).
        Then:  04_clustered.h5ad is deleted and the exception propagates;
              03_integrated.h5ad is untouched.
        """
        step03, adata = _make_step04_state()
        cfg = _make_config(tmp_path)
        cfg.incremental_io = True
        _write_step03(cfg, step03)

        with patch(
            "core.utils.write_h5ad_incremental",
            side_effect=RuntimeError("simulated append crash"),
        ):
            with pytest.raises(RuntimeError, match="simulated append crash"):
                _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        assert not os.path.exists(cfg.cluster_h5ad), (
            "corrupt 04 copy must be deleted on failure (mode A)"
        )
        assert os.path.exists(cfg.integrated_h5ad), "03 source must survive"

    def test_incremental_io_false_uses_full_safe_write(self, tmp_path) -> None:
        """incremental_io=False → full safe_write fallback (escape hatch)."""
        _, adata = _make_step04_state()
        cfg = _make_config(tmp_path)
        cfg.incremental_io = False

        with patch.object(_mod, "safe_write") as mock_safe:
            _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        mock_safe.assert_called_once_with(adata, cfg.cluster_h5ad, cfg=cfg)

    def test_incremental_io_missing_from_old_config_defaults_true(self, tmp_path) -> None:
        """Old configs without incremental_io → default True (copy+append)."""
        import types

        _, adata = _make_step04_state()
        # MagicMock auto-creates attributes; use a plain namespace that genuinely
        # LACKS incremental_io to simulate a pre-T2 config.
        cfg = types.SimpleNamespace(
            integrated_h5ad=str(tmp_path / "03_integrated.h5ad"),
            cluster_h5ad=str(tmp_path / "04_clustered.h5ad"),
        )
        _write_step03(cfg, _make_step04_state()[0])

        with patch("core.utils.write_h5ad_incremental") as mock_write:
            _mod._write_cluster_h5ad(adata, cfg, MagicMock())

        assert mock_write.call_count == 1
