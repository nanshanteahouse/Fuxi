"""Step-00 memory preflight guard tests (plan T9, metis G7/G8).

``rna/steps/00_load.py::_preflight_step00_meta`` is the metadata-only pre-scan
that feeds ``estimate_step_peak(0, ...)`` + ``check_memory_guard`` BEFORE the
heavy load.  These tests pin the probe itself (per format) and the guard
behaviour (warn continues / block raises), mirroring ``01_doublet.py``.

The guard is a no-op by default: the schema default ``execution.memory.guard``
is "warn", so a run with a sane budget only logs ``[memory-guard]`` and never
blocks unless the config opts into ``guard=block``.
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from scipy.sparse import csr_matrix

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_guard_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

from core.utils._memory import (  # noqa: E402
    check_memory_guard,
    estimate_step_peak,
    resolve_memory_settings,
)


def _write_10x_h5(path: Path, n_cells: int, n_genes: int, seed: int) -> None:
    """Minimal 10x v3 h5 with a dense counts matrix (like the T1b fixtures)."""
    rng = np.random.RandomState(seed)
    x = rng.randint(0, 4, size=(n_cells, n_genes)).astype(np.float32)
    csr = csr_matrix(x)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int32))
        g.create_dataset("data", data=csr.data.astype(np.int32))
        g.create_dataset("indices", data=csr.indices.astype(np.int64))
        g.create_dataset("indptr", data=csr.indptr.astype(np.int64))
        g.create_dataset(
            "barcodes", data=np.array([f"b{i}".encode() for i in range(n_cells)], dtype="S20")
        )
        feats = g.create_group("features")
        feats.create_dataset(
            "name", data=np.array([f"g{i}".encode() for i in range(n_genes)], dtype="S20")
        )
        feats.create_dataset(
            "id", data=np.array([f"e{i}".encode() for i in range(n_genes)], dtype="S20")
        )
        feats.create_dataset(
            "feature_type", data=np.array([b"Gene Expression"] * n_genes, dtype="S20")
        )
        feats.create_dataset("genome", data=np.array([b"GRCh38"] * n_genes, dtype="S20"))


def _make_h5_cfg(h5_dir: str, pattern: str = "*.h5") -> SimpleNamespace:
    return SimpleNamespace(
        data_format="10X_h5",
        data_input=SimpleNamespace(h5_dir=h5_dir, h5_file_pattern=pattern),
        data_dir=h5_dir,
    )


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_load_guard")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


def _exec_cfg(guard: str = "warn", budget: str = "64GB") -> SimpleNamespace:
    return SimpleNamespace(
        execution=SimpleNamespace(
            memory=SimpleNamespace(policy="speed", budget=budget, guard=guard)
        )
    )


# ═════════════════════════════════════════════════════════════════════
#  Metadata pre-scan probes
# ═════════════════════════════════════════════════════════════════════


def test_preflight_10x_h5_single_file(tmp_path: Path) -> None:
    d = tmp_path / "h5"
    d.mkdir()
    _write_10x_h5(d / "s1.h5", n_cells=50, n_genes=12, seed=1)
    meta = _mod._preflight_step00_meta(_make_h5_cfg(str(d)), _make_logger())
    assert meta is not None
    n_cells, n_genes, nnz, concat = meta
    assert n_cells == 50
    assert n_genes == 12
    assert nnz > 0 and nnz < 50 * 12  # sparse counts, exact from h5 data shape
    assert concat == 1.0  # single-file → no union-var growth


def test_preflight_10x_h5_multi_file_concat(tmp_path: Path) -> None:
    d = tmp_path / "h5"
    d.mkdir()
    _write_10x_h5(d / "s1.h5", n_cells=50, n_genes=12, seed=1)
    _write_10x_h5(d / "s2.h5", n_cells=70, n_genes=12, seed=2)
    meta = _mod._preflight_step00_meta(_make_h5_cfg(str(d)), _make_logger())
    assert meta is not None
    n_cells, n_genes, nnz, concat = meta
    assert n_cells == 120
    assert n_genes == 12
    assert concat == 1.3  # multi-file merge → union-var growth bound
    assert nnz > 0
    assert n_cells == 120
    assert concat == 1.3  # multi-file merge → union-var growth bound
    assert nnz > 0


def test_preflight_csv_table_estimates(tmp_path: Path) -> None:
    p = tmp_path / "matrix.csv"
    # genes×cells table with header (first row = barcode columns);
    # comma files need explicit csv_sep (as in GSE173180's config).
    genes, cells = ["G1", "G2", "G3"], ["C1", "C2"]
    with open(p, "w") as f:
        f.write(",".join([""] + cells) + "\n")
        for g in genes:
            f.write(f"{g},1,0\n")
    cfg = SimpleNamespace(
        data_format="csv_matrix",
        data_input=SimpleNamespace(matrix_file=str(p), csv_sep=","),
    )
    meta = _mod._preflight_step00_meta(cfg, _make_logger())
    assert meta is not None
    n_cells, n_genes, nnz, concat = meta
    assert n_cells == 2
    assert n_genes == 3
    assert 0 < nnz <= 2 * 3
    assert concat == 1.0


def test_preflight_csv_table_sep_autodetect(tmp_path: Path) -> None:
    # Tab-separated default (the branch's auto-detect path): header on one
    # line -> sep=None peek splits >1 column -> tab wins.
    p = tmp_path / "matrix.tsv"
    with open(p, "w") as f:
        f.write("\tC1\tC2\nG1\t1\t0\nG2\t0\t1\n")
    cfg = SimpleNamespace(
        data_format="csv_matrix",
        data_input=SimpleNamespace(matrix_file=str(p), csv_sep=None),
    )
    meta = _mod._preflight_step00_meta(cfg, _make_logger())
    assert meta is not None
    n_cells, n_genes, nnz, concat = meta
    assert n_cells == 2
    assert n_genes == 2
    assert 0 < nnz <= 4
    assert concat == 1.0
    assert n_cells == 2
    assert n_genes == 2
    assert 0 < nnz <= 4


def test_preflight_csv_mtx_header(tmp_path: Path) -> None:
    p = tmp_path / "matrix.mtx"
    p.write_text("%%MatrixMarket matrix coordinate integer general\n%\n4 3 5\n1 1 1\n")
    cfg = SimpleNamespace(
        data_format="csv_matrix",
        data_input=SimpleNamespace(matrix_file=str(p)),
    )
    meta = _mod._preflight_step00_meta(cfg, _make_logger())
    assert meta is not None
    n_cells, n_genes, nnz, concat = meta
    assert (n_cells, n_genes, nnz, concat) == (3, 4, 5, 1.0)


def test_preflight_preprocessed_shapes(tmp_path: Path) -> None:
    d = tmp_path / "pre"
    d.mkdir()
    f1 = d / "f1.tsv.gz"
    f2 = d / "f2.tsv.gz"
    import gzip

    for fp, n_rows in ((f1, 4), (f2, 3)):
        with gzip.open(fp, "wt") as f:
            f.write("barcode\tg1\tg2\n")
            for i in range(n_rows):
                f.write(f"b{i}\t1\t0\n")
    cfg = SimpleNamespace(
        data_format="preprocessed",
        data_input=SimpleNamespace(separator="\t"),
        data_dir=str(d),
    )
    meta = _mod._preflight_step00_meta(cfg, _make_logger())
    assert meta is not None
    n_cells, n_genes, nnz, concat = meta
    assert n_cells == 7
    assert n_genes == 2  # 3 header cols minus the barcode column
    assert 0 < nnz <= 7 * 2
    assert concat == 1.0


def test_preflight_unsupported_returns_none() -> None:
    cfg = SimpleNamespace(data_format="h5ad", data_input=SimpleNamespace())
    assert _mod._preflight_step00_meta(cfg, _make_logger()) is None


def test_mtx_header_info_gz(tmp_path: Path) -> None:
    import gzip

    p = tmp_path / "m.mtx.gz"
    with gzip.open(p, "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n2 9 7\n")
    assert _mod._mtx_header_info(str(p)) == (2, 9, 7)


# ═════════════════════════════════════════════════════════════════════
#  Guard behaviour (warn continues / block raises / off skips)
# ═════════════════════════════════════════════════════════════════════


def _estimate_from_preflight(meta) -> float:
    n_cells, n_genes, nnz, concat = meta
    return estimate_step_peak(0, n_cells, n_genes, nnz, concat_factor=concat)


def test_guard_warn_continues_over_budget(tmp_path: Path) -> None:
    d = tmp_path / "h5"
    d.mkdir()
    _write_10x_h5(d / "s1.h5", n_cells=5_000, n_genes=2_000, seed=3)
    meta = _mod._preflight_step00_meta(_make_h5_cfg(str(d)), _make_logger())
    assert meta is not None
    est = {0: _estimate_from_preflight(meta)}
    tiny_budget = int(2 * 2**30)  # 2 GiB — everything is over it
    # warn: logs a warning and returns True (run may proceed) — the DEFAULT.
    assert check_memory_guard(est, tiny_budget, "warn", logger_obj=_make_logger()) is True


def test_guard_block_raises_over_budget(tmp_path: Path) -> None:
    d = tmp_path / "h5"
    d.mkdir()
    _write_10x_h5(d / "s1.h5", n_cells=5_000, n_genes=2_000, seed=3)
    meta = _mod._preflight_step00_meta(_make_h5_cfg(str(d)), _make_logger())
    assert meta is not None
    est = {0: _estimate_from_preflight(meta)}
    tiny_budget = int(2 * 2**30)
    with pytest.raises(RuntimeError, match=r"\[memory-guard\]"):
        check_memory_guard(est, tiny_budget, "block", logger_obj=_make_logger())


def test_guard_block_within_budget_passes(tmp_path: Path) -> None:
    d = tmp_path / "h5"
    d.mkdir()
    _write_10x_h5(d / "s1.h5", n_cells=1_000, n_genes=500, seed=4)
    meta = _mod._preflight_step00_meta(_make_h5_cfg(str(d)), _make_logger())
    assert meta is not None
    est = {0: _estimate_from_preflight(meta)}
    assert check_memory_guard(est, int(64 * 2**30), "block", logger_obj=_make_logger()) is True


def test_guard_off_skips(tmp_path: Path) -> None:
    d = tmp_path / "h5"
    d.mkdir()
    _write_10x_h5(d / "s1.h5", n_cells=5_000, n_genes=2_000, seed=3)
    meta = _mod._preflight_step00_meta(_make_h5_cfg(str(d)), _make_logger())
    assert meta is not None
    est = {0: _estimate_from_preflight(meta)}
    assert check_memory_guard(est, int(2 * 2**30), "off", logger_obj=_make_logger()) is True


def test_estimate_step00_uses_cfg_resolution() -> None:
    # resolve_memory_settings default must be guard="warn" (schema default) so
    # existing step-00 runs are a no-op unless the config opts into "block".
    policy, budget, guard = resolve_memory_settings(_exec_cfg())
    assert guard == "warn"
    assert policy == "speed"
    assert budget > 0
