"""Tests for multi-sample 10X MTX loading in rna/steps/00_load.py.

Covers
------
- fast path: identical gene sets across sample dirs → sparse hstack
- fallback path: differing gene sets → one-shot outer-join concat
- fail-fast: missing matrix file → SystemExit with clear message
- no matching dirs → SystemExit
- sample naming: dir basename + optional mtx_sample_regex
- sample_map remap of dir-name samples
- legacy genes.tsv.gz → features.tsv.gz conversion
- mtx_prefix detection
"""

from __future__ import annotations

import gzip
import importlib
import logging
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "00_load.py")
_spec = importlib.util.spec_from_file_location("rna.steps._00_load_test", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════


def _write_features(path: Path, gene_symbols: list[str]) -> None:
    with gzip.open(path, "wt") as f:
        for i, g in enumerate(gene_symbols):
            f.write(f"ENSG{i:05d}\t{g}\tGene Expression\n")


def _write_legacy_genes(path: Path, gene_symbols: list[str]) -> None:
    with gzip.open(path, "wt") as f:
        for i, g in enumerate(gene_symbols):
            f.write(f"ENSG{i:05d}\t{g}\n")


def _write_mtx_dir(
    path: Path,
    gene_symbols: list[str],
    barcodes: list[str],
    seed: int = 0,
    legacy: bool = False,
) -> None:
    """Write a minimal 10X MTX sample directory (genes × cells matrix)."""
    path.mkdir(parents=True, exist_ok=True)
    x = np.random.RandomState(seed).randint(0, 5, size=(len(gene_symbols), len(barcodes)))
    if legacy:
        _write_legacy_genes(path / "genes.tsv.gz", gene_symbols)
    else:
        _write_features(path / "features.tsv.gz", gene_symbols)
    with gzip.open(path / "barcodes.tsv.gz", "wt") as f:
        for bc in barcodes:
            f.write(bc + "\n")
    with gzip.open(path / "matrix.mtx.gz", "wb") as f:
        mmwrite(f, csr_matrix(x))


def _make_cfg(
    mtx_dir: str,
    mtx_dir_pattern: str = "*",
    mtx_sample_regex: str = "",
    mtx_concat_batch: int = 0,
    sample_map: dict | None = None,
    has_sample_mapping: bool = False,
) -> MagicMock:
    cfg = MagicMock()
    cfg.data_input.mtx_dir = mtx_dir
    cfg.data_input.mtx_dir_pattern = mtx_dir_pattern
    cfg.data_input.mtx_sample_regex = mtx_sample_regex
    cfg.data_input.gene_symbol_column = ""
    cfg.data_input.mtx_concat_batch = mtx_concat_batch
    cfg.sample_meta.sample_map = sample_map or {}
    cfg.sample_meta.barcode_parse_regex = ""
    cfg.sample_meta.barcode_parse_groups = {}
    cfg.has_sample_mapping.return_value = has_sample_mapping
    cfg.has_stage_mapping.return_value = False
    return cfg


def _make_logger() -> logging.Logger:
    log = logging.getLogger("test_00_multi_mtx")
    log.handlers = []
    log.addHandler(logging.NullHandler())
    return log


GENES_A = [f"GENE_{i:03d}" for i in range(20)]
GENES_B = GENES_A + ["GENE_EXTRA"]  # superset of A


# ═════════════════════════════════════════════════════════════════════
#  Fast path (identical gene sets)
# ═════════════════════════════════════════════════════════════════════


class TestFastPath:
    def test_hstack_identical_genes(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "sampleA", GENES_A, [f"AAACCCA{i:03d}" for i in range(5)], seed=1)
        _write_mtx_dir(parent / "sampleB", GENES_A, [f"AAACCCB{i:03d}" for i in range(7)], seed=2)

        adata = _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), _make_logger())

        assert adata.n_obs == 12
        assert adata.n_vars == len(GENES_A)
        assert list(adata.var_names) == GENES_A
        assert set(adata.obs["sample"]) == {"sampleA", "sampleB"}
        # barcode suffixes match the 10X_h5 multi-file convention (-0, -1, ...)
        assert all(str(bc).endswith("-0") for bc in adata.obs_names[:5])
        assert all(str(bc).endswith("-1") for bc in adata.obs_names[5:])
        assert adata.obs_names.is_unique

    def test_sample_map_remap(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "sampleA", GENES_A, [f"AAACCCA{i:03d}" for i in range(3)])
        _write_mtx_dir(parent / "sampleB", GENES_A, [f"AAACCCB{i:03d}" for i in range(4)])

        cfg = _make_cfg(
            str(parent),
            sample_map={"sampleA": "WT", "sampleB": "KO"},
            has_sample_mapping=True,
        )
        adata = _mod._load_multi_sample_10x_mtx(cfg, _make_logger())
        assert set(adata.obs["sample"]) == {"WT", "KO"}

    def test_sample_regex_extraction(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "GSM1_WT_10x", GENES_A, [f"AAACCCA{i:03d}" for i in range(3)])
        _write_mtx_dir(parent / "GSM2_KO_10x", GENES_A, [f"AAACCCB{i:03d}" for i in range(3)])

        cfg = _make_cfg(str(parent), mtx_sample_regex=r"GSM\d+_([A-Za-z]+)")
        adata = _mod._load_multi_sample_10x_mtx(cfg, _make_logger())
        assert set(adata.obs["sample"]) == {"WT", "KO"}


# ═════════════════════════════════════════════════════════════════════
#  Fallback path (differing gene sets)
# ═════════════════════════════════════════════════════════════════════


class TestFallbackPath:
    def test_outer_join_gene_union(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "sampleA", GENES_A, [f"AAACCCA{i:03d}" for i in range(5)], seed=1)
        _write_mtx_dir(parent / "sampleB", GENES_B, [f"AAACCCB{i:03d}" for i in range(7)], seed=2)

        adata = _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), _make_logger())

        assert adata.n_obs == 12
        assert adata.n_vars == len(GENES_B)  # union (B is superset)
        assert "GENE_EXTRA" in adata.var_names
        # GENE_EXTRA only expressed in sampleB cells (sampleA block is zero-filled)
        sample_a_mask = adata.obs["sample"].values == "sampleA"
        extra_col = list(adata.var_names).index("GENE_EXTRA")
        assert adata.X[sample_a_mask, extra_col].toarray().sum() == 0
        assert adata.X[~sample_a_mask, extra_col].toarray().sum() > 0

    def test_batched_concat(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        for i, genes in enumerate([GENES_A, GENES_B, GENES_A, GENES_B]):
            _write_mtx_dir(
                parent / f"s{i}", genes, [f"AAACCC{i:02d}{j:02d}" for j in range(4)], seed=i
            )

        cfg = _make_cfg(str(parent), mtx_concat_batch=2)
        adata = _mod._load_multi_sample_10x_mtx(cfg, _make_logger())
        assert adata.n_obs == 16
        assert adata.n_vars == len(GENES_B)


# ═════════════════════════════════════════════════════════════════════
#  Fail-fast
# ═════════════════════════════════════════════════════════════════════


class TestFailFast:
    def test_no_matching_dirs_exits(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        parent.mkdir()
        (parent / "not_mtx.txt").write_text("x")
        with pytest.raises(SystemExit):
            _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), _make_logger())

    def test_missing_matrix_exits_with_sample_name(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "ok", GENES_A, ["AAACCCA001"])
        broken = parent / "broken"
        broken.mkdir(parents=True)
        _write_features(broken / "features.tsv.gz", GENES_A)
        with gzip.open(broken / "barcodes.tsv.gz", "wt") as f:
            f.write("AAACCCA002\n")
        # no matrix.mtx.gz → read fails → SystemExit

        log = _make_logger()
        log.error = MagicMock()
        with pytest.raises(SystemExit):
            _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), log)
        assert log.error.called
        # error message should name the broken sample and its matrix file
        msg = " ".join(str(a) for a in log.error.call_args.args)
        assert "broken" in msg

    def test_missing_features_exits(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "sampleA", GENES_A, ["AAACCCA001"])
        broken = parent / "broken"
        broken.mkdir(parents=True)
        x = np.ones((5, 2))
        with gzip.open(broken / "matrix.mtx.gz", "wb") as f:
            mmwrite(f, csr_matrix(x))
        with pytest.raises(SystemExit):
            _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), _make_logger())


# ═════════════════════════════════════════════════════════════════════
#  Prefix detection & legacy conversion
# ═════════════════════════════════════════════════════════════════════


class TestHelpers:
    def test_detect_mtx_prefix(self, tmp_path: Path) -> None:
        d = tmp_path / "prefixed"
        d.mkdir()
        (d / "GSM1_matrix.mtx.gz").write_text("x")
        assert _mod._detect_mtx_prefix(str(d)) == "GSM1_"
        (d / "GSM1_matrix.mtx.gz").unlink()
        (d / "GSM1_matrix.mtx").write_text("x")
        assert _mod._detect_mtx_prefix(str(d)) == "GSM1_"
        d2 = tmp_path / "plain"
        d2.mkdir()
        (d2 / "matrix.mtx.gz").write_text("x")
        assert _mod._detect_mtx_prefix(str(d2)) == ""

    def test_legacy_genes_conversion(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        _write_mtx_dir(
            parent / "sampleA", GENES_A, [f"AAACCCA{i:03d}" for i in range(3)], legacy=True
        )

        log = _make_logger()
        adata = _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), log)

        assert adata.n_obs == 3
        assert adata.n_vars == len(GENES_A)
        # legacy genes.tsv.gz converted in place to 3-column features.tsv.gz
        features = parent / "sampleA" / "features.tsv.gz"
        assert features.exists()
        with gzip.open(features, "rt") as f:
            first = f.readline().rstrip("\n")
        assert first.count("\t") == 2  # id, symbol, Gene Expression

    def test_prefixed_sample_dirs(self, tmp_path: Path) -> None:
        parent = tmp_path / "mtx"
        for i in range(2):
            d = parent / f"GSM{i}_sample"
            d.mkdir(parents=True)
            _write_features(d / "GSM{}_sample_features.tsv.gz".format(i), GENES_A)
            with gzip.open(d / f"GSM{i}_sample_barcodes.tsv.gz", "wt") as f:
                f.write(f"AAACCCA{i}01\n")
            with gzip.open(d / f"GSM{i}_sample_matrix.mtx.gz", "wb") as f:
                mmwrite(f, csr_matrix(np.ones((20, 1))))

        adata = _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), _make_logger())
        assert adata.n_obs == 2
        assert set(adata.obs["sample"]) == {"GSM0_sample", "GSM1_sample"}
        adata = _mod._load_multi_sample_10x_mtx(_make_cfg(str(parent)), _make_logger())
        assert adata.n_obs == 2
        assert set(adata.obs["sample"]) == {"GSM0_sample", "GSM1_sample"}


# ═════════════════════════════════════════════════════════════════════
#  generate_config() multi-dir fix
# ═════════════════════════════════════════════════════════════════════


class TestGenerateConfigMultiMtx:
    """generate_config() must not silently drop multi-sample MTX dirs."""

    @staticmethod
    def _build_classification(dirs: list[Path]) -> dict:
        return {
            "tenx_mtx_dirs": {str(d): [str(p) for p in sorted(d.iterdir())] for d in dirs},
            "tenx_h5_dirs": {},
            "fragment_dirs": {},
            "tenx_peak_dirs": {},
            "h5ad_files": [],
            "csv_files": [],
            "metadata_files": [],
            "archives": [],
            "unmatched": [],
            "unsupported": [],
        }

    def test_multi_dir_emits_pattern(self, tmp_path: Path) -> None:
        import core.preprocess.matrix_loader as ml

        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "GSM1_sampleA", GENES_A, ["AAACCCA001"])
        _write_mtx_dir(parent / "GSM2_sampleB", GENES_A, ["AAACCCB001"])

        config_path = ml.generate_config(
            gse_id="GSE000999",
            modality="rna",
            classification=self._build_classification(
                [parent / "GSM1_sampleA", parent / "GSM2_sampleB"]
            ),
            file_list=[str(p) for p in (parent / "GSM1_sampleA").iterdir()],
            output_dir=str(tmp_path / "out"),
            data_root=str(tmp_path / "data"),
            force=True,
        )
        assert config_path is not None

        text = Path(config_path).read_text()
        assert "mtx_dir_pattern:" in text
        assert str(parent) in text  # absolute parent mtx_dir

    def test_multi_dir_includes_all_dirs(self, tmp_path: Path, capsys) -> None:
        import core.preprocess.matrix_loader as ml

        parent = tmp_path / "mtx"
        _write_mtx_dir(parent / "GSM1_sampleA", GENES_A, ["AAACCCA001"])
        _write_mtx_dir(parent / "GSM2_sampleB", GENES_A, ["AAACCCB001"])
        _write_mtx_dir(parent / "GSM3_sampleC", GENES_A, ["AAACCCC001"])

        dirs = [parent / "GSM1_sampleA", parent / "GSM2_sampleB", parent / "GSM3_sampleC"]
        ml.generate_config(
            gse_id="GSE000998",
            modality="rna",
            classification=self._build_classification(dirs),
            file_list=[str(p) for p in dirs[0].iterdir()],
            output_dir=str(tmp_path / "out"),
            data_root=str(tmp_path / "data"),
            force=True,
        )
        captured = capsys.readouterr()
        assert "GSM1_sampleA" in captured.out
        assert "GSM2_sampleB" in captured.out
        assert "GSM3_sampleC" in captured.out

    def test_single_dir_keeps_legacy(self, tmp_path: Path) -> None:
        import core.preprocess.matrix_loader as ml

        parent = tmp_path / "mtx"
        _write_mtx_dir(parent, GENES_A, ["AAACCCA001"])
        config_path = ml.generate_config(
            gse_id="GSE000997",
            modality="rna",
            classification=self._build_classification([parent]),
            file_list=[str(p) for p in parent.iterdir()],
            output_dir=str(tmp_path / "out"),
            data_root=str(tmp_path / "data"),
            force=True,
        )
        assert config_path is not None
        text = Path(config_path).read_text()
        assert "mtx_dir_pattern:" in text
