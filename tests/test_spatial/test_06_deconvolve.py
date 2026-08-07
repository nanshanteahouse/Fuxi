"""Tests for spatial/steps/06_deconvolve.py — cell2location deconvolution.

TDD RED phase (spatial pipeline rewrite phase 2 plan, todo 8 — new step 06,
original 06-12 slide to 07-13 in Wave 3; this todo only creates the file
plus its tests, no runner/schema wiring).

Behaviour under test (all with cell2location + mygene mocked — no real
training / network queries):

- ``deconv_method == "none"`` → graceful skip (uns['deconvolution']
  status='skipped'), no model calls, no fabricated proportions.
- cell2location not importable (``_get_models`` → None) → log.warning +
  graceful skip.
- No scRNA reference (rna_ref empty + auto-discovery empty) → graceful skip.
- Normal path: RegressionModel reference estimation → export_posterior →
  inf_aver (varm slice) → Cell2location → export_posterior → proportions
  (spot × cell-type, rows sum to 1) written to obsm['deconv_proportions']
  and table_dir/proportions.csv. ``train`` must be called with
  ``accelerator=`` and never with ``use_gpu=`` (scvi-tools 1.5 removed it).
- ENSG→symbol alignment via mygene: spatial var (ENSG) mapped to symbols
  present in the reference; unmapped genes dropped with a log; overlap
  < 50% → log.error + safe_write (state preserved) + sys.exit(1).
- NMF tissue zones when ``deconv_n_factors > 0`` → obs['tissue_zone'] +
  tissue_zone_composition.csv.
- Boundary cases: single cell type (proportions all 1.0), single spot
  (row sums to 1).
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

# ── Ensure repo root is on sys.path ──────────────────────────────────────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 06_deconvolve module via file path ──────────────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "spatial", "steps", "06_deconvolve.py")
_spec = importlib.util.spec_from_file_location(
    "spatial.steps._06_deconvolve_test",
    _STEP_PATH,
)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

# ── Shared test constants ────────────────────────────────────────────────
CELL_TYPES = ["T Cell", "B Cell", "Monocyte", "Epithelial"]
N_SPOTS = 12
N_GENES = 40


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════


def _make_spot_adata(
    n_spots: int = N_SPOTS,
    n_genes: int = N_GENES,
    seed: int = 7,
    ensembl: bool = False,
) -> AnnData:
    """Minimal spot-level AnnData (05_annotated shape: spots × genes)."""
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=2.0, size=(n_spots, n_genes)).astype(np.float32)
    adata = AnnData(x)
    if ensembl:
        adata.var_names = [f"ENSG0000000{i:04d}" for i in range(n_genes)]
    else:
        adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.obs["sample"] = ["S1", "S2"] * (n_spots // 2) if n_spots >= 2 else ["S1"]
    return adata


def _make_ref_adata(n_cells: int = 50, n_genes: int = N_GENES, seed: int = 11) -> AnnData:
    """Minimal scRNA reference (05_annotated shape: cells × genes, symbols)."""
    rng = np.random.RandomState(seed)
    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = rng.choice(CELL_TYPES, n_cells)
    return adata


def _make_cfg(tmp_path) -> MagicMock:
    """Config mock. All deconv fields set explicitly (getattr on a MagicMock
    would return another MagicMock instead of the default)."""
    cfg = MagicMock()

    cfg.spatial = MagicMock()
    cfg.spatial.deconv_method = "cell2location"
    cfg.spatial.deconv_n_factors = 0
    cfg.spatial.deconv_max_epochs = 30000
    cfg.spatial.deconv_batch_size = 256
    cfg.spatial.deconv_n_cells_per_location = 8
    cfg.spatial.deconv_detection_alpha = 20.0
    cfg.spatial.deconv_ref_max_epochs = 500
    cfg.spatial.deconv_ref_n_samples = 1000
    cfg.spatial.deconv_ref_batch_size = 2500

    cfg.rna_ref = ""  # must be explicitly "" (empty → auto-discover)

    cfg.execution = MagicMock()
    cfg.execution.device = "cpu"
    cfg.execution.random_seed = 42

    cfg.h5ad_dir = str(tmp_path)
    cfg.log_dir = str(tmp_path)
    cfg.table_dir = str(tmp_path / "tables")
    cfg.figure_dir = str(tmp_path)
    return cfg


def _make_fake_models(
    cell_types: list[str] | None = None,
    seed: int = 3,
) -> tuple[MagicMock, MagicMock, MagicMock]:
    """MagicMock 'cell2location.models' that mimics the standard flow.

    RegressionModel.export_posterior populates a real ``varm`` +
    ``uns['mod']`` so the step's inf_aver extraction runs on real data.
    Cell2location.export_posterior populates a real ``obsm`` abundance
    DataFrame so proportions are computed from real values.
    """
    cell_types = cell_types or list(CELL_TYPES)
    fake_models = MagicMock()

    def _fake_ref_export(adata, sample_kwargs=None, **kw):
        rng = np.random.RandomState(seed)
        varm = pd.DataFrame(
            rng.rand(adata.n_vars, len(cell_types)),
            index=adata.var_names,
            columns=[f"means_per_cluster_mu_fg_{ct}" for ct in cell_types],
        )
        adata.varm["means_per_cluster_mu_fg"] = varm
        adata.uns["mod"] = {"factor_names": list(cell_types)}
        return adata

    fake_sc_model = MagicMock()
    fake_sc_model.export_posterior.side_effect = _fake_ref_export
    fake_models.RegressionModel = MagicMock(return_value=fake_sc_model)

    def _fake_spot_export(adata, sample_kwargs=None, **kw):
        rng = np.random.RandomState(seed + 1)
        abund = pd.DataFrame(
            rng.rand(adata.n_obs, len(cell_types)),
            index=adata.obs_names,
            columns=list(cell_types),
        )
        adata.obsm["q05_cell_abundance_w_sf"] = abund
        return adata

    fake_c2l_model = MagicMock()
    fake_c2l_model.export_posterior.side_effect = _fake_spot_export
    fake_models.Cell2location = MagicMock(return_value=fake_c2l_model)
    return fake_models, fake_sc_model, fake_c2l_model


def _common_patches(spot_adata: AnnData, ref_adata: AnnData, cfg: MagicMock) -> list:
    """Shared patch chain (I/O + find_rna_h5ad auto-discovery)."""

    def _fake_read(path, *a, **k):
        if "rna_ref" in str(path):
            return ref_adata
        return spot_adata

    return [
        patch.object(
            _mod.argparse.ArgumentParser,
            "parse_args",
            return_value=argparse.Namespace(config="/tmp/test.yaml"),
        ),
        patch.object(_mod, "resolve_config", return_value=cfg),
        patch.object(_mod, "setup_logger", return_value=MagicMock()),
        patch.object(_mod.sc, "read", side_effect=_fake_read),
        patch.object(_mod, "find_rna_h5ad", return_value="/tmp/rna_ref/05_annotated.h5ad"),
    ]


def _run_main(
    spot_adata: AnnData,
    ref_adata: AnnData,
    cfg: MagicMock,
    extra_patches: list,
) -> None:
    """Enter all patch contexts via ExitStack, then run main()."""
    from contextlib import ExitStack

    with ExitStack() as stack:
        for p in _common_patches(spot_adata, ref_adata, cfg) + list(extra_patches):
            stack.enter_context(p)
        _mod.main()


def _ensure_input(tmp_path) -> None:
    """Create the Step 05 input file so the existence gate passes."""
    (tmp_path / "05_annotated.h5ad").touch()


def _capture(captured: list) -> callable:
    """safe_write side-effect that captures the final adata."""

    def _side_effect(adata, *args, **kwargs):
        captured.append(adata)

    return _side_effect


# ═══════════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════════


def test_deconv_method_none_skips(tmp_path) -> None:
    """deconv_method == 'none' → graceful skip with uns metadata.

    Given:  method='none'.
    When:   main() runs.
    Then:   06_deconvolved.h5ad written with uns['deconvolution']
            status='skipped'; no model training, no fabricated proportions.
    """
    adata = _make_spot_adata()
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    cfg.spatial.deconv_method = "none"
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    fake_models, fake_sc, fake_c2l = _make_fake_models()
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "_get_models", return_value=fake_models),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1, "safe_write should be called once (skip record)"
    result = captured[0]
    assert result.uns["deconvolution"]["status"] == "skipped"
    assert "deconv_method" in result.uns["deconvolution"]["reason"]
    assert "deconv_proportions" not in result.obsm, "no fabricated proportions on skip"
    fake_sc.train.assert_not_called()
    fake_c2l.train.assert_not_called()


def test_cell2location_missing_skips(tmp_path) -> None:
    """cell2location not importable → log.warning + graceful skip.

    Given:  _get_models() returns None (ImportError path).
    When:   main() runs.
    Then:   status='skipped' with a cell2location reason; no model calls.
    """
    adata = _make_spot_adata()
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "_get_models", return_value=None),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1
    result = captured[0]
    assert result.uns["deconvolution"]["status"] == "skipped"
    assert "cell2location" in result.uns["deconvolution"]["reason"]


def test_no_rna_ref_skips(tmp_path) -> None:
    """No scRNA reference → graceful skip with uns metadata.

    Given:  rna_ref empty AND auto-discovery returns None.
    When:   main() runs.
    Then:   status='skipped' with a reference reason; no training.
    """
    adata = _make_spot_adata()
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    fake_models, fake_sc, _ = _make_fake_models()
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "find_rna_h5ad", return_value=None),
            patch.object(_mod, "_get_models", return_value=fake_models),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1
    result = captured[0]
    assert result.uns["deconvolution"]["status"] == "skipped"
    assert "reference" in result.uns["deconvolution"]["reason"]
    fake_sc.train.assert_not_called()


def test_normal_path_produces_proportions(tmp_path) -> None:
    """Full mocked flow → spot×cell-type proportions, accelerator not use_gpu.

    Given:  symbol-var spot data + symbol reference, device='cpu'.
    When:   main() runs.
    Then:   train called with accelerator= (never use_gpu); Cell2location
            ctor receives inf_aver + N_cells_per_location + detection_alpha;
            obsm['deconv_proportions'] rows sum to 1; proportions.csv written;
            uns['deconvolution'] status='completed'.
    """
    adata = _make_spot_adata()
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    cfg.execution.device = "cpu"
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    fake_models, fake_sc, fake_c2l = _make_fake_models()
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "_get_models", return_value=fake_models),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1
    result = captured[0]

    # ── train uses accelerator, never use_gpu (scvi-tools 1.5 removed it) ──
    assert "accelerator" in fake_sc.train.call_args.kwargs
    assert fake_sc.train.call_args.kwargs["accelerator"] == "cpu"
    assert "use_gpu" not in fake_sc.train.call_args.kwargs
    assert "accelerator" in fake_c2l.train.call_args.kwargs
    assert fake_c2l.train.call_args.kwargs["accelerator"] == "cpu"
    assert "use_gpu" not in fake_c2l.train.call_args.kwargs

    # ── Cell2location ctor kwargs ──
    ctor_kwargs = fake_models.Cell2location.call_args.kwargs
    assert ctor_kwargs["N_cells_per_location"] == 8
    assert ctor_kwargs["detection_alpha"] == 20.0
    inf_aver = ctor_kwargs["cell_state_df"]
    assert isinstance(inf_aver, pd.DataFrame)
    assert inf_aver.shape == (N_GENES, len(CELL_TYPES))
    assert inf_aver.index.tolist() == [f"GENE_{i}" for i in range(N_GENES)]
    assert inf_aver.columns.tolist() == CELL_TYPES

    # ── Output h5ad ──
    assert "deconv_proportions" in result.obsm
    props = pd.DataFrame(
        result.obsm["deconv_proportions"],
        index=adata.obs_names,
        columns=CELL_TYPES,
    )
    assert props.shape == (N_SPOTS, len(CELL_TYPES))
    np.testing.assert_allclose(
        props.sum(axis=1).values,
        np.ones(N_SPOTS),
        atol=1e-6,
        err_msg="each spot's proportion row must sum to 1",
    )

    # ── proportions.csv ──
    csv_path = os.path.join(cfg.table_dir, "proportions.csv")
    assert os.path.exists(csv_path)
    csv_df = pd.read_csv(csv_path, index_col=0)
    assert csv_df.shape == (N_SPOTS, len(CELL_TYPES))
    assert list(csv_df.columns) == CELL_TYPES

    # ── progress record ──
    rec = result.uns["deconvolution"]
    assert rec["status"] == "completed"
    assert rec["method"] == "cell2location"
    assert rec["n_cell_types"] == len(CELL_TYPES)
    assert rec["n_spots"] == N_SPOTS
    assert rec["n_genes_mapped"] == N_GENES
    assert rec["accelerator"] == "cpu"
    assert "wall_sec" in rec


def test_ensg_to_symbol_mapping_aligns(tmp_path) -> None:
    """Spatial var is ENSG, reference is symbol → mygene mapping aligns them.

    Given:  spot var_names are ENSG ids; ref var_names are symbols.
    When:   main() runs with _mygene_query mocked to ENSG→symbol.
    Then:   mygene queried once; inf_aver index is symbols; n_genes_mapped
            recorded; proportions produced.
    """
    adata = _make_spot_adata(ensembl=True)
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    mapping = {g: f"GENE_{i}" for i, g in enumerate(adata.var_names)}
    captured: list[AnnData] = []
    fake_models, _, _ = _make_fake_models()
    with patch.object(_mod, "_mygene_query", return_value=mapping) as mg_mock:
        _run_main(
            adata,
            ref,
            cfg,
            [
                patch.object(_mod, "_get_models", return_value=fake_models),
                patch.object(_mod, "safe_write", side_effect=_capture(captured)),
            ],
        )

    assert len(captured) == 1
    result = captured[0]
    assert result.uns["deconvolution"]["status"] == "completed"
    assert result.uns["deconvolution"]["n_genes_mapped"] == N_GENES

    # mygene received the ENSG ids
    mg_mock.assert_called_once()
    queried = mg_mock.call_args.args[0]
    assert queried[:3] == ["ENSG00000000000", "ENSG00000000001", "ENSG00000000002"]

    # inf_aver index must be symbols (post-rename), matching ref genes
    inf_aver = fake_models.Cell2location.call_args.kwargs["cell_state_df"]
    assert inf_aver.shape == (N_GENES, len(CELL_TYPES))
    assert inf_aver.index.tolist()[:3] == ["GENE_0", "GENE_1", "GENE_2"]


def test_ensg_mapping_drops_unmapped_genes(tmp_path) -> None:
    """ENSG ids mapping to symbols absent from the reference are dropped.

    Given:  spot has 40 ENSG ids; ref only has symbols GENE_0..19 (50%
            overlap — boundary passes).
    When:   main() runs.
    Then:   only 20 genes survive; inf_aver is 20×n_types; proportions.csv
            uses the surviving gene set.
    """
    n_survive = 20
    adata = _make_spot_adata(n_genes=N_GENES, ensembl=True)
    ref = _make_ref_adata(n_genes=n_survive)
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    # only the first 20 ENSG ids map into the reference's symbol space
    mapping = {g: f"GENE_{i}" for i, g in enumerate(adata.var_names) if i < n_survive}
    captured: list[AnnData] = []
    fake_models, _, _ = _make_fake_models()
    with patch.object(_mod, "_mygene_query", return_value=mapping):
        _run_main(
            adata,
            ref,
            cfg,
            [
                patch.object(_mod, "_get_models", return_value=fake_models),
                patch.object(_mod, "safe_write", side_effect=_capture(captured)),
            ],
        )

    assert len(captured) == 1
    result = captured[0]
    assert result.uns["deconvolution"]["status"] == "completed"
    assert result.uns["deconvolution"]["n_genes_mapped"] == n_survive

    inf_aver = fake_models.Cell2location.call_args.kwargs["cell_state_df"]
    assert inf_aver.shape == (n_survive, len(CELL_TYPES))
    assert inf_aver.index.tolist() == [f"GENE_{i}" for i in range(n_survive)]


def test_low_gene_overlap_exits(tmp_path) -> None:
    """<50% gene overlap → log.error + safe_write (preserve) + sys.exit(1).

    Given:  spot has 40 genes, ref only 10 (25% overlap).
    When:   main() runs.
    Then:   SystemExit(1); state preserved via safe_write with
            uns['deconvolution'] status='failed'; no model training.
    """
    adata = _make_spot_adata(n_genes=N_GENES)
    ref = _make_ref_adata(n_genes=10)
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    fake_models, fake_sc, _ = _make_fake_models()
    with pytest.raises(SystemExit) as excinfo:
        _run_main(
            adata,
            ref,
            cfg,
            [
                patch.object(_mod, "_get_models", return_value=fake_models),
                patch.object(_mod, "safe_write", side_effect=_capture(captured)),
            ],
        )

    assert excinfo.value.code == 1
    assert len(captured) == 1, "state must be preserved before exit"
    assert captured[0].uns["deconvolution"]["status"] == "failed"
    assert "overlap" in captured[0].uns["deconvolution"]["reason"].lower()
    fake_sc.train.assert_not_called()


def test_nmf_tissue_zones(tmp_path) -> None:
    """deconv_n_factors > 0 → obs['tissue_zone'] + composition CSV.

    Given:  deconv_n_factors=3, mocked CoLocatedGroupsSklearnNMF.
    When:   main() runs.
    Then:   tissue_zone assigned per spot; tissue_zone_composition.csv
            written (cell types × factors).
    """
    adata = _make_spot_adata()
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    cfg.spatial.deconv_n_factors = 3
    _ensure_input(tmp_path)

    fake_models, _, _ = _make_fake_models()
    fake_nmf = MagicMock()
    rng = np.random.RandomState(9)

    def _fake_sample2df(node_name=None, ct_node_name=None, **kw):
        fake_nmf.location_factors_df = pd.DataFrame(
            rng.rand(N_SPOTS, 3),
            index=adata.obs_names,
            columns=[f"mean_location_factors_factor{i}" for i in range(3)],
        )
        fake_nmf.cell_type_fractions = pd.DataFrame(
            rng.rand(len(CELL_TYPES), 3),
            index=list(CELL_TYPES),
            columns=[f"mean_cell_type_factors_factor{i}" for i in range(3)],
        )

    fake_nmf.sample2df.side_effect = _fake_sample2df
    fake_models.CoLocatedGroupsSklearnNMF = MagicMock(return_value=fake_nmf)

    captured: list[AnnData] = []
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "_get_models", return_value=fake_models),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1
    result = captured[0]
    assert "tissue_zone" in result.obs, "NMF zones must be recorded in obs"
    assert result.obs["tissue_zone"].nunique() <= 3
    assert result.uns["deconvolution"]["n_factors"] == 3

    comp_csv = os.path.join(cfg.table_dir, "tissue_zone_composition.csv")
    assert os.path.exists(comp_csv)
    comp = pd.read_csv(comp_csv, index_col=0)
    assert comp.shape == (len(CELL_TYPES), 3)


def test_single_cell_type(tmp_path) -> None:
    """Reference with one cell type → proportions are all 1.0.

    Given:  ref obs['cell_type'] has a single unique value.
    When:   main() runs.
    Then:   proportions shape (n_spots, 1), every value 1.0.
    """
    adata = _make_spot_adata()
    ref = _make_ref_adata()
    ref.obs["cell_type"] = ["OnlyType"] * ref.n_obs
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    fake_models, _, _ = _make_fake_models(cell_types=["OnlyType"])
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "_get_models", return_value=fake_models),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1
    result = captured[0]
    props = pd.DataFrame(
        result.obsm["deconv_proportions"],
        index=adata.obs_names,
        columns=["OnlyType"],
    )
    assert props.shape == (N_SPOTS, 1)
    np.testing.assert_allclose(props.values, np.ones((N_SPOTS, 1)), atol=1e-6)


def test_single_spot(tmp_path) -> None:
    """Single spot → single proportion row summing to 1.

    Given:  spot data with n_obs == 1.
    When:   main() runs.
    Then:   proportions shape (1, n_types), row sums to 1.
    """
    adata = _make_spot_adata(n_spots=1)
    ref = _make_ref_adata()
    cfg = _make_cfg(tmp_path)
    _ensure_input(tmp_path)

    captured: list[AnnData] = []
    fake_models, _, _ = _make_fake_models()
    _run_main(
        adata,
        ref,
        cfg,
        [
            patch.object(_mod, "_get_models", return_value=fake_models),
            patch.object(_mod, "safe_write", side_effect=_capture(captured)),
        ],
    )

    assert len(captured) == 1
    result = captured[0]
    props = pd.DataFrame(
        result.obsm["deconv_proportions"],
        index=adata.obs_names,
        columns=CELL_TYPES,
    )
    assert props.shape == (1, len(CELL_TYPES))
    np.testing.assert_allclose(props.sum(axis=1).values, np.array([1.0]), atol=1e-6)
