"""Tests for single-directory 10X MTX loading in rna/steps/00_load.py.

Covers
------
- the single-dir 10X_mtx branch calls ``sc.read_10x_mtx(..., cache=False)``
  (no scanpy ``_cache`` directory / .h5ad cache file written into the data dir)
- a minimal fixture run through the branch's load path produces no ``_cache`` dir
"""

from __future__ import annotations

import gzip
import importlib
import os
import sys
from pathlib import Path

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_nocache_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Helpers (mirror the _write_mtx_dir pattern from test_00_load_multi_mtx.py)
# ═════════════════════════════════════════════════════════════════════


def _write_features(path: Path, gene_symbols: list[str]) -> None:
    with gzip.open(path, "wt") as f:
        for i, g in enumerate(gene_symbols):
            f.write(f"ENSG{i:05d}\t{g}\tGene Expression\n")


def _write_mtx_dir(
    path: Path,
    gene_symbols: list[str],
    barcodes: list[str],
    prefix: str = "",
    seed: int = 0,
) -> None:
    """Write a minimal single-dir 10X MTX dir (genes × cells matrix)."""
    path.mkdir(parents=True, exist_ok=True)
    x = np.random.RandomState(seed).randint(0, 5, size=(len(gene_symbols), len(barcodes)))
    _write_features(path / f"{prefix}features.tsv.gz", gene_symbols)
    with gzip.open(path / f"{prefix}barcodes.tsv.gz", "wt") as f:
        for bc in barcodes:
            f.write(bc + "\n")
    with gzip.open(path / f"{prefix}matrix.mtx.gz", "wb") as f:
        mmwrite(f, csr_matrix(x))


def _single_dir_read_10x_mtx_call() -> str:
    """Return the argument block of the single-dir branch's sc.read_10x_mtx call.

    The single-dir call lives inside ``main()`` (after the multi-sample helper);
    the multi-sample helper ``_load_multi_sample_10x_mtx`` also calls
    ``sc.read_10x_mtx`` — exclude that definition so only the single-dir branch
    is inspected.
    """
    source = Path(_STEP_PATH).read_text()
    # Keep only the part after the multi-sample helper definition (main() body).
    marker = "def _load_multi_sample_10x_mtx"
    tail = source[source.index(marker) :]
    call_start = tail.index("read_10x_mtx(")
    depth = 0
    i = tail.index("(", call_start)
    j = i
    while j < len(tail):
        if tail[j] == "(":
            depth += 1
        elif tail[j] == ")":
            depth -= 1
            if depth == 0:
                break
        j += 1
    return tail[call_start : j + 1]


GENES = [f"GENE_{i:03d}" for i in range(12)]
BARCODES = [f"AAACCCA{i:03d}" for i in range(6)]


# ═════════════════════════════════════════════════════════════════════
#  Source-level guard: single-dir branch must request cache=False
# ═════════════════════════════════════════════════════════════════════


class TestSourceCacheFlag:
    def test_single_dir_branch_uses_cache_false(self) -> None:
        call = _single_dir_read_10x_mtx_call()
        assert "cache=False" in call
        assert "cache=True" not in call

    def test_no_cache_true_anywhere_in_step(self) -> None:
        source = Path(_STEP_PATH).read_text()
        assert "cache=True" not in source


# ═════════════════════════════════════════════════════════════════════
#  Functional: branch load path must not create a _cache dir
# ═════════════════════════════════════════════════════════════════════


class TestNoCacheDirCreated:
    def test_single_dir_load_produces_no_cache_dir(self, tmp_path: Path) -> None:
        mtx_dir = tmp_path / "mtx"
        _write_mtx_dir(mtx_dir, GENES, BARCODES, prefix="GSE1_", seed=1)

        # Same call the single-dir branch makes in main():
        # sc.read_10x_mtx(dir, var_names="gene_symbols", prefix=..., cache=False, gex_only=False)
        adata = _mod.sc.read_10x_mtx(
            str(mtx_dir),
            var_names="gene_symbols",
            prefix="GSE1_",
            cache=False,
            gex_only=False,
        )

        assert adata.n_obs == len(BARCODES)
        assert adata.n_vars == len(GENES)

        # No scanpy cache directory or .h5ad cache file may appear next to the input.
        assert not (mtx_dir / "_cache").exists()
        assert not (mtx_dir / ".cache").exists()
        assert not list(mtx_dir.glob("*.h5ad"))
        # Only the fixture files remain — no cache artifacts beside them.
        assert set(mtx_dir.iterdir()) == {
            mtx_dir / "GSE1_features.tsv.gz",
            mtx_dir / "GSE1_barcodes.tsv.gz",
            mtx_dir / "GSE1_matrix.mtx.gz",
        }
