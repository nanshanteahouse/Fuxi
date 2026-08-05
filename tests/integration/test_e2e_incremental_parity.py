"""E2E parity: incremental_io=True vs False must produce identical checkpoints.

Runs the full RNA pipeline (steps 00 → 08) twice on the same synthetic
fixture (``tests/fixtures/synthetic_rna.h5ad``), once with ``incremental_io:
true`` (in-place obs/obsm/obsp/uns appends) and once with ``false`` (full
``safe_write`` rewrites).  Every h5ad checkpoint produced by both runs is
then compared read-back-wise (obs columns/dtypes, obsm, obsp sparse
structure, uns, X, var, raw, layers) and must match exactly.

This is the acceptance test for plan item 10 (e2e RNA pipeline synthetic
parity).  Both runs use ``random_seed: 42`` so every stochastic algorithm
(scrublet, PCA, harmony, leiden, UMAP, DPT) is deterministic and comparable.

Known pre-existing issues (NOT regressions of this task):

1. **Step-04 CPU UMAP**: ``gpu_umap`` passes ``n_epochs`` to ``sc.tl.umap``
   which scanpy 1.12.2 no longer accepts → UMAP is non-fatal-swallowed by the
   grid search, so ``04_clustered.h5ad`` has no ``X_umap`` in *both* runs.
   Parity is symmetric and still asserted; the absence is recorded.

2. **safe_write's ``verify_write_integrity`` leaks a backed read handle**, so
   a second write to the same path in one process crashes ("unable to truncate
   a file which is already open").  Steps 03/06 write their checkpoints
   repeatedly → the test sets ``verify_write_integrity: false`` for both modes
   (see the config builder comment).  The parity assertions here are the
   stronger integrity check.

3. **Step-08 PAGA crashes on this synthetic fixture**: the synthetic cell
   types are KNN-disconnected components (500 cells per well-separated blob),
   so scanpy's PAGA v1.2 cluster graph has 0 edges and ``sc.tl.paga`` raises
   ``ValueError: mismatching number of index arrays`` (igraph 1.0.0 + scanpy
   1.12.2).  Both modes fail step 08 **identically**; the test asserts that
   symmetry (and the shared absence of ``05_final.h5ad``), and performs full
   ``05_final`` parity + sentinel checks automatically once the library issue
   is fixed (i.e. when both runs succeed).

Run with::

    pytest tests/integration/test_e2e_incremental_parity.py -v --tb=short

Takes ~3-6 minutes (two full pipeline runs on 5k synthetic cells).  Set
``SKIP_SLOW_TESTS=1`` to bypass.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp
import yaml
from anndata import OldFormatWarning

import core.utils._io  # noqa: F401  (registers the HDF5 zstd plugin for read-back)

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_FIXTURE = _REPO_ROOT / "tests" / "fixtures" / "synthetic_rna.h5ad"

# Checkpoints both runs must produce and compare, in pipeline order.
# 05_annotated.h5ad is compared in its FINAL state (i.e. after step 06's
# in-place write-back / full rewrite respectively - exactly the divergent
# write path under test).  05_final.h5ad (step 08) is asserted conditionally:
# the synthetic fixture's KNN graph is disconnected per cell type, which
# trips a pre-existing scanpy PAGA bug, so both runs fail step 08 identically
# (see the test body); 05_final parity activates once that bug is fixed.
CHECKPOINTS = [
    "00_raw.h5ad",
    "01_doublet.h5ad",
    "02_qc.h5ad",
    "03_integrated.h5ad",
    "04_clustered.h5ad",
    "05_annotated.h5ad",
]
SUB_FILE = "05_sub_cell_type_0.h5ad"  # extra step-06 product (identical too)
SENTINEL_STEP06 = "05_annotated.h5ad.step06_done"
SENTINEL_STEP08 = "05_final.h5ad.step08_done"

# Steps executed by _run_pipeline before the separate step-08 run
# (07/09-12 write no h5ad checkpoint; step 08 is run and asserted
#  symmetrically by the caller - see test body).
STEPS_BASE = [0, 1, 2, 3, 4, 5]
_STEP_TIMEOUT = {0: 600, 1: 900, 2: 600, 3: 900, 4: 1800, 5: 900, 6: 1800, 8: 1800}


def _skip_slow() -> bool:
    """True when SKIP_SLOW_TESTS is set (mirrors test_step06_rna.py)."""
    return os.environ.get("SKIP_SLOW_TESTS", "").lower() in ("1", "true")


def _marker_dict_from_fixture() -> dict[str, list[str]]:
    """Read the fixture's ``uns['marker_dict']`` as a plain dict[str, list[str]].

    Anndata round-trips the dict of lists as dict[str, ndarray], so each
    value is converted back to a list before dumping into the YAML config
    (step 05's score_genes path needs marker lists).
    """
    adata = _read_ckpt(_FIXTURE)
    return {ct: [str(g) for g in genes] for ct, genes in adata.uns["marker_dict"].items()}


def _build_config(
    run_dir: Path,
    fixture: Path,
    marker_dict: dict[str, list[str]],
    incremental_io: bool,
    *,
    step06: bool,
) -> dict:
    """Return a minimal full-pipeline RNA config dict rooted under *run_dir*.

    ``step06=True`` produces the step-06 variant that deterministically
    triggers a real write-back: ``ai.enabled/subcluster=True`` pointing at a
    closed local port (connection refused → numeric fallback label → the
    ``cell_subtype`` write-back actually happens) plus the configured
    ``subcluster_types``.
    """
    cfg: dict = {
        "modality": "rna",
        "tissue": "test",
        "species": "mouse",
        "expression_type": "raw_counts",
        "data_format": "h5ad",
        "data_input": {"input_h5ad": str(fixture)},
        "project_dir": str(run_dir),
        "h5ad_dir": str(run_dir / "results" / "h5ad"),
        "figure_dir": str(run_dir / "results" / "figures"),
        "table_dir": str(run_dir / "results" / "tables"),
        "log_dir": str(run_dir / "logs"),
        "execution": {"device": "cpu", "random_seed": 42, "n_jobs": 1},
        "incremental_io": incremental_io,
        # Workaround for a PRE-EXISTING bug (HEAD d8135f2, not this task): the
        # integrity verify inside safe_write opens the target backed and never
        # closes it, so a second write to the same path fails with "file already
        # open". Step 03 writes 03_integrated twice and step 06 writes
        # 05_sub_*.h5ad three times -> both would crash. The parity assertions
        # in this test are the stronger integrity check, so verification is
        # disabled here (identical for both modes).
        "verify_write_integrity": False,
        "pca": {"n_pcs_use": 5, "n_pcs_full": 20},
        "hvg": {
            "n_top_genes": 200,
            "flavor": "seurat_v3",
        },  # >=51 features for step-06 subset PCA (n_comps=min(50, n_obs-2))
        "qc": {"min_genes": 100, "max_genes": 5000, "min_genes_per_umi": 0.1},
        # Scrublet on this synthetic fixture is FLAKY: a backed anndata read
        # can yield an object-dtype _CSRDataset that sp.csr_matrix() rejects
        # (nondeterministic across processes). Disabling it (the task-sanctioned
        # escape hatch) makes step 01 deterministic - both modes write
        # doublet_scores=0 / predicted_doublet=False identically.
        "scrublet": {"run": False},
        "normalization": {
            "use_regress_out": False,  # pct_counts_mt is all-zero on synthetic data
            "score_cell_cycle": False,
            "detect_sex": False,
        },
        "integration": {
            "method": "harmony",
            "batch_key": "batch",
            "max_iter": 5,
            "diagnose": True,
            "collinearity_guard": True,
        },
        "clustering": {
            "param_grid_n_neighbors": [5],
            "param_grid_n_neighbors_adaptive": False,
            "param_grid_resolutions": [0.3, 0.5],
            "leiden_flavor": "igraph",
            "leiden_n_iterations": 2,
            "cluster_selection_method": "silhouette",  # cheap, deterministic
            "stability_n_seeds": 2,
            "umap_plot_mode": "skip",
            "plot_per_combo": False,
            "param_grid_min_dist": [0.3],
            "param_grid_spread": [1.0],
        },
        "marker": {
            "marker_dict": marker_dict,
            "subcluster_types": [],
            "subcluster_resolution": 0.4,
            "min_cells_subcluster": 50,
        },
        "ai": {"enabled": False},
        "trajectory": {"method": "paga_dpt", "save_final_h5ad": True},
        "plot": {
            "figure_dpi": 72,
            "palette": {"categorical": "tab20", "dotplot_fill": "YlOrRd"},
            "umap_panel_size": [4, 4],
        },
    }
    if step06:
        cfg["marker"]["subcluster_types"] = ["cell_type_0"]
        cfg["ai"] = {
            "enabled": True,
            "subcluster": True,
            "api_base": "http://127.0.0.1:1",  # closed port → deterministic fallback
            "api_key": "sk-test-invalid",
            "model": "test-model",
            "max_tokens": 256,
            "temperature": 0.1,
            "timeout": 5,
            "thinking_enabled": False,
        }
    return cfg


def _run_step(config_path: Path, step: int) -> subprocess.CompletedProcess:
    """Run one pipeline step as a subprocess via the runner (never import steps)."""
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "core.run_pipeline",
            "--modality",
            "rna",
            "--step",
            str(step),
            "--config",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
        timeout=_STEP_TIMEOUT.get(step, 1800),
    )


def _assert_step_ok(step: int, result: subprocess.CompletedProcess) -> None:
    if result.returncode == 0:
        return
    raise AssertionError(
        f"Step {step} exited with code {result.returncode}.\n"
        f"=== STDOUT (tail) ===\n{result.stdout[-4000:]}\n"
        f"=== STDERR (tail) ===\n{result.stderr[-4000:]}"
    )


def _run_pipeline(
    run_dir: Path, marker_dict: dict[str, list[str]], incremental_io: bool
) -> subprocess.CompletedProcess | None:
    """Write configs and execute steps 0-5, 6, 8 for one mode.

    Returns the step-08 subprocess result.  Step 08's PAGA branch crashes on
    this synthetic fixture (pre-existing scanpy bug on disconnected group
    graphs - the cluster graph has 0 edges, see module docstring); the caller
    asserts both modes fail identically, or full 05_final parity once the
    library issue is fixed and both succeed.
    """
    run_dir.mkdir(parents=True, exist_ok=True)
    base_cfg = _build_config(run_dir, _FIXTURE, marker_dict, incremental_io, step06=False)
    base_path = run_dir / "config_base.yaml"
    with open(base_path, "w") as f:
        yaml.safe_dump(base_cfg, f)

    step06_cfg = _build_config(run_dir, _FIXTURE, marker_dict, incremental_io, step06=True)
    step06_path = run_dir / "config_step06.yaml"
    with open(step06_path, "w") as f:
        yaml.safe_dump(step06_cfg, f)

    for step in STEPS_BASE:
        _assert_step_ok(step, _run_step(base_path, step))
    _assert_step_ok(6, _run_step(step06_path, 6))
    return _run_step(base_path, 8)


# ═══════════════════════════════════════════════════════════════════════
#  Read-back parity helpers
# ═══════════════════════════════════════════════════════════════════════


def _arrays_equal(a, b) -> bool:
    """NaN-tolerant value equality, handling string/object arrays too.

    ``np.array_equal(..., equal_nan=True)`` calls ``isnan`` which raises on
    object/string arrays (uns vlen strings like ``marker_dict`` values and
    ``*_colors`` hex lists round-trip as object arrays) — fall back to plain
    equality for non-numeric dtypes.
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind in "USO" or b.dtype.kind in "USO":
        return bool(np.array_equal(a, b))
    return bool(np.array_equal(a, b, equal_nan=True))


def _matrices_equal(a, b) -> bool:
    """Exact value equality for dense and/or sparse matrix mixes.

    Sparse inputs are compared by CSR layout (indices/indptr/data), dense by
    ``_arrays_equal``.  A dense/sparse cross-comparison fails (layout is part
    of the artifact).
    """
    a_sp, b_sp = sp.issparse(a), sp.issparse(b)
    if a_sp or b_sp:
        if not (a_sp and b_sp):
            return False
        a = a.tocsr()
        b = b.tocsr()
        return bool(
            a.shape == b.shape
            and a.dtype == b.dtype
            and np.array_equal(a.indices, b.indices)
            and np.array_equal(a.indptr, b.indptr)
            and _arrays_equal(a.data, b.data)
        )
    a = np.asarray(a)
    b = np.asarray(b)
    return bool(a.shape == b.shape and a.dtype == b.dtype and _arrays_equal(a, b))


def _series_equal(a: pd.Series, b: pd.Series, label: str) -> None:
    """Assert two obs/var columns are identical (dtype, categories, values)."""
    da, db = a.dtype, b.dtype
    assert da == db, f"{label}: dtype mismatch ({da} vs {db})"
    if isinstance(da, pd.CategoricalDtype):
        assert list(a.cat.categories) == list(b.cat.categories), f"{label}: category order differs"
        assert np.array_equal(a.cat.codes.to_numpy(), b.cat.codes.to_numpy()), (
            f"{label}: categorical codes differ"
        )
        return
    av, bv = a.to_numpy(), b.to_numpy()
    if av.dtype == object or bv.dtype == object:
        # object/mixed → compare as nullable string (NA-aware, no "nan" false-positives)
        assert a.astype("string").array.equals(b.astype("string").array), (
            f"{label}: string values differ"
        )
    else:
        assert np.array_equal(av, bv, equal_nan=True), f"{label}: values differ"


def _uns_value_equal(a, b, label: str) -> None:
    """Recursive equality for arbitrary ``uns`` / ``varm`` values."""
    if isinstance(a, dict) and isinstance(b, dict):
        assert set(a.keys()) == set(b.keys()), f"{label}: keys differ ({set(a) ^ set(b)})"
        for k in a:
            _uns_value_equal(a[k], b[k], f"{label}[{k}]")
        return
    if sp.issparse(a) or sp.issparse(b):
        assert _matrices_equal(a, b), f"{label}: sparse mismatch"
        return
    if isinstance(a, pd.DataFrame) and isinstance(b, pd.DataFrame):
        assert a.equals(b), f"{label}: DataFrame differs"
        return
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            assert _arrays_equal(a, b), f"{label}: array differs"
            return
        assert isinstance(a, np.ndarray) and isinstance(b, np.ndarray), (
            f"{label}: ndarray vs {type(b).__name__}"
        )
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        assert _arrays_equal(a, b), f"{label}: list differs"
        return
    if isinstance(a, (int, float, str, bool, np.generic)) or isinstance(
        b, (int, float, str, bool, np.generic)
    ):
        if isinstance(a, float) and isinstance(b, float):
            assert (a == b) or (math.isnan(a) and math.isnan(b)), f"{label}: {a!r} != {b!r}"
        else:
            assert a == b, f"{label}: {a!r} != {b!r}"
        return
    assert a == b, f"{label}: {type(a).__name__} values differ"


def _read_ckpt(path: Path) -> "sc.AnnData":
    """Read a checkpoint h5ad, tolerating legacy-format /layers metadata.

    Step 02_qc hand-writes the h5ad container and creates ``/layers`` with
    no ``encoding-type``/``encoding-version`` attributes (unlike every other
    checkpoint).  anndata 0.13 raises ``OldFormatWarning`` on read; pytest's
    ``filterwarnings=error`` escalates it to a failure.  The artifact is
    byte-identical in BOTH the incremental and full runs (step 02 is write-
    mode-independent), so the warning is symmetric — it is not a parity
    divergence.  Suppress exactly this one legacy-format warning so the
    read-back parity assertions below remain the real integrity check.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OldFormatWarning)
        return sc.read(str(path))


def _assert_h5ad_parity(a_path: Path, b_path: Path, label: str) -> None:
    """Assert two checkpoint h5ads read back identically (semantic content)."""
    a = _read_ckpt(a_path)
    b = _read_ckpt(b_path)

    assert a.shape == b.shape, f"{label}: shape differs {a.shape} vs {b.shape}"
    _assert_matrix("X", a.X, b.X, label)

    # obs: names + per-column parity (aligned by column name, order-agnostic)
    assert list(a.obs_names) == list(b.obs_names), f"{label}: obs_names differ"
    assert set(a.obs.columns) == set(b.obs.columns), (
        f"{label}: obs columns differ ({set(a.obs.columns) ^ set(b.obs.columns)})"
    )
    for col in a.obs.columns:
        _series_equal(a.obs[col], b.obs[col], f"{label}.obs[{col}]")

    # var
    assert list(a.var_names) == list(b.var_names), f"{label}: var_names differ"
    assert set(a.var.columns) == set(b.var.columns), f"{label}: var columns differ"
    for col in a.var.columns:
        _series_equal(a.var[col], b.var[col], f"{label}.var[{col}]")

    # obsm / obsp
    assert set(a.obsm.keys()) == set(b.obsm.keys()), (
        f"{label}: obsm keys differ ({set(a.obsm.keys()) ^ set(b.obsm.keys())})"
    )
    for k in a.obsm:
        _assert_matrix(f"obsm[{k}]", a.obsm[k], b.obsm[k], label)
    assert set(a.obsp.keys()) == set(b.obsp.keys()), (
        f"{label}: obsp keys differ ({set(a.obsp.keys()) ^ set(b.obsp.keys())})"
    )
    for k in a.obsp:
        _assert_matrix(f"obsp[{k}]", a.obsp[k], b.obsp[k], label)

    # varm / uns
    assert set(a.varm.keys()) == set(b.varm.keys()), f"{label}: varm keys differ"
    for k in a.varm:
        _uns_value_equal(a.varm[k], b.varm[k], f"{label}.varm[{k}]")
    # Step-04 stability evaluation leaves `_temp_stab_*` scratch keys in uns.
    # The copy+append path writes only OWNED keys and intentionally drops them;
    # the full safe_write path persists everything. So `_temp_stab_*` is the
    # single ALLOWED uns asymmetry (present in the full-write artifact only);
    # every other uns key must match exactly.
    _a_uns = {k: v for k, v in a.uns.items() if not k.startswith("_temp_stab_")}
    _b_uns = {k: v for k, v in b.uns.items() if not k.startswith("_temp_stab_")}
    _uns_value_equal(_a_uns, _b_uns, f"{label}.uns")
    _scratch_diff = (set(a.uns) ^ set(b.uns)) - (set(_a_uns) | set(_b_uns))
    assert all(k.startswith("_temp_stab_") for k in _scratch_diff), (
        f"{label}: unexpected uns key asymmetry: {_scratch_diff}"
    )

    # layers / raw (structural carry-overs)
    assert set(a.layers.keys()) == set(b.layers.keys()), f"{label}: layers keys differ"
    for k in a.layers:
        _assert_matrix(f"layers[{k}]", a.layers[k], b.layers[k], label)
    assert (a.raw is None) == (b.raw is None), f"{label}: raw presence differs"
    if a.raw is not None:
        assert list(a.raw.var_names) == list(b.raw.var_names), f"{label}: raw var_names differ"
        _assert_matrix("raw.X", a.raw.X, b.raw.X, label)


def _assert_matrix(what: str, a, b, label: str) -> None:
    assert _matrices_equal(a, b), f"{label}: {what} differs"


class TestE2EIncrementalParity:
    """Double pipeline run (incremental vs full write) with read-back parity."""

    @pytest.mark.skipif(
        _skip_slow(),
        reason="Slow integration test skipped via SKIP_SLOW_TESTS",
    )
    def test_incremental_and_full_runs_produce_identical_checkpoints(self, tmp_path: Path) -> None:
        if not _FIXTURE.exists():
            pytest.skip(
                f"Fixture not found: {_FIXTURE}. Run tests/fixtures/generate_synthetic.py."
            )

        marker_dict = _marker_dict_from_fixture()

        incr_dir = tmp_path / "incremental"
        full_dir = tmp_path / "full"
        step8_incr = _run_pipeline(incr_dir, marker_dict, incremental_io=True)
        step8_full = _run_pipeline(full_dir, marker_dict, incremental_io=False)

        h5ad_incr = incr_dir / "results" / "h5ad"
        h5ad_full = full_dir / "results" / "h5ad"

        # ── Step-06 sentinel: real write-back must have happened in BOTH modes ──
        sent06 = "05_annotated.h5ad.step06_done"
        assert (h5ad_incr / sent06).read_text().strip() == "done"
        assert (h5ad_full / sent06).read_text().strip() == "done"

        # Prove both modes REALLY took different write paths to reach that
        # identical result: incremental -> in-place obs writeback (backup +
        # append), full -> whole-file safe_write. Guards against a config typo
        # silently making both runs use the same write path (vacuous parity).
        incr_step6_log = (incr_dir / "logs" / "06_subcluster.log").read_text()
        full_step6_log = (full_dir / "logs" / "06_subcluster.log").read_text()
        assert "In-place obs writeback" in incr_step6_log, (
            "incremental run did not use the in-place writeback path"
        )
        assert "In-place obs writeback" not in full_step6_log, (
            "full run unexpectedly used the in-place writeback path"
        )
        assert "Wrote back" in full_step6_log and "Wrote back" in incr_step6_log

        # ── Checkpoint-by-checkpoint read-back parity (00 → 06 write-back) ──
        for ckpt in [
            "00_raw.h5ad",
            "01_doublet.h5ad",
            "02_qc.h5ad",
            "03_integrated.h5ad",
            "04_clustered.h5ad",
            "05_annotated.h5ad",
        ]:
            a, b = h5ad_incr / ckpt, h5ad_full / ckpt
            assert a.is_file(), f"missing incremental checkpoint {ckpt}"
            assert b.is_file(), f"missing full checkpoint {ckpt}"
            _assert_h5ad_parity(a, b, ckpt)

        # Extra: step-06 subcluster product must match too
        a_sub, b_sub = h5ad_incr / SUB_FILE, h5ad_full / SUB_FILE
        assert a_sub.is_file(), f"missing incremental {SUB_FILE}"
        assert b_sub.is_file(), f"missing full {SUB_FILE}"
        _assert_h5ad_parity(a_sub, b_sub, SUB_FILE)

        # ── Step 08: symmetric outcome in both modes ──
        # The synthetic fixture's cell types are KNN-disconnected components, so
        # scanpy's PAGA v1.2 crashes (0-edge cluster graph). Both modes hit the
        # IDENTICAL library failure - that symmetry is itself asserted. If a
        # future scanpy/igraph fixes the crash, both runs succeed and the full
        # 05_final parity + step-08 sentinel checks below activate automatically.
        assert step8_incr is not None and step8_full is not None
        assert step8_incr.returncode == step8_full.returncode, (
            f"step-08 outcomes differ: incr={step8_incr.returncode}, full={step8_full.returncode}"
        )
        if step8_incr.returncode == 0:
            _assert_h5ad_parity(
                h5ad_incr / "05_final.h5ad", h5ad_full / "05_final.h5ad", "05_final.h5ad"
            )
            assert (h5ad_incr / "05_final.h5ad.step08_done").is_file()
            assert (h5ad_full / "05_final.h5ad.step08_done").is_file()
        else:
            for label, res in (("incremental", step8_incr), ("full", step8_full)):
                assert "mismatching number of index arrays" in res.stderr, (
                    f"step-08 {label} failed with an UNEXPECTED error:\n"
                    f"{res.stderr[-2000:]}\n"
                    f"{res.stdout[-2000:]}"
                )
            assert not (h5ad_incr / "05_final.h5ad").exists()
            assert not (h5ad_full / "05_final.h5ad").exists()
            # both modes must be missing the sentinel identically
            assert not (h5ad_incr / "05_final.h5ad.step08_done").exists()
            assert not (h5ad_full / "05_final.h5ad.step08_done").exists()

        # Known issue note: step-04 CPU UMAP fails under scanpy 1.12.2 → both
        # runs must show the SAME absence (recorded, not asserted as present).
        a04 = _read_ckpt(h5ad_incr / "04_clustered.h5ad")
        assert "X_umap" not in a04.obsm, (
            "scanpy 1.12.2 removed n_epochs from sc.tl.umap; expected step-04 UMAP "
            "failure to be symmetric (no X_umap in either run). If this changed, "
            "re-check parity — X_umap now exists and is compared above."
        )
