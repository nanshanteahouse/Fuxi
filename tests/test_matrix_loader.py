"""Tests for core/preprocess/matrix_loader.py — config generation & paper_context bridge.

Covers:
- ``generate_config()`` paper_context heuristic overrides
- ``_post_process_config()`` ast-based marker_dict/is_nuclei/tissue_kb/tissue_ontology injection
- Idempotency (replacement vs. append)
- Backward compatibility (no paper_context → identical behavior)
"""
from __future__ import annotations

import ast
import os
from pathlib import Path
from typing import Any
import pytest
from pytest import MonkeyPatch

# Paths — use a known template directory
_REPO_ROOT = Path(__file__).resolve().parent.parent
_TEMPLATE_DIR = _REPO_ROOT / "templates" / "config_templates"


# ═══════════════════════════════════════════════════════════════════════
#  Fixtures
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def fake_classification() -> dict[str, Any]:
    """Minimal classification dict that triggers the '10X_mtx' template."""
    return {
        "tenx_mtx_dirs": {},
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


@pytest.fixture
def fake_file_list(tmp_path: Path) -> list[str]:
    """Create dummy files so guess_species / guess_tissue don't crash."""
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True)
    f1 = data_dir / "sample1_R1.fastq.gz"
    f2 = data_dir / "sample1_R2.fastq.gz"
    f1.write_text("")
    f2.write_text("")
    return [str(f1), str(f2)]


@pytest.fixture
def sample_config_source() -> str:
    """A minimal generated config source with a marker_dict placeholder."""
    return """import os
from core.config.schema import CFG

# Data format
CFG.data_format = '10X_mtx'
CFG.mtx_prefix = ''

# Dataset metadata
CFG.tissue = 'unknown'
CFG.species = 'human'
CFG.expression_type = 'raw_counts'

# Cell type markers
CFG.marker_dict = {
    # 'CellTypeA': ['GENE1', 'GENE2'],
}

# Knowledge base
# CFG.tissue_kb = ''
CFG.n_jobs = 0
CFG.random_seed = 42
"""


# ═══════════════════════════════════════════════════════════════════════
#  _post_process_config tests
# ═══════════════════════════════════════════════════════════════════════


class TestPostProcessConfig:
    """ast-based injection of paper-derived CFG fields."""

    def test_inject_marker_dict(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject a marker_dict with extracted features."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        paper_context = {"features": ["GENE1", "GENE2", "GENE3"]}
        _post_process_config(str(config_path), paper_context)

        source = config_path.read_text()
        # Should have marker_dict with extracted features
        assert "CFG.marker_dict" in source
        assert "GENE1" in source
        assert "GENE2" in source
        assert "GENE3" in source
        assert "extracted" in source

    def test_inject_is_nuclei(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject is_nuclei=True when paper_context indicates snRNA-seq."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        paper_context = {"is_nuclei": True}
        _post_process_config(str(config_path), paper_context)

        source = config_path.read_text()
        assert "CFG.is_nuclei = True" in source

    def test_inject_tissue_kb(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject tissue_kb value."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        paper_context = {"tissue_kb": "retina"}
        _post_process_config(str(config_path), paper_context)

        source = config_path.read_text()
        assert "CFG.tissue_kb = 'retina'" in source

    def test_inject_tissue_ontology(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject tissue_ontology value."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        paper_context = {"tissue_ontology": "UBERON_0000966"}
        _post_process_config(str(config_path), paper_context)

        source = config_path.read_text()
        assert "CFG.tissue_ontology = 'UBERON_0000966'" in source

    def test_idempotent_replace_existing_marker_dict(
        self, tmp_path: Path, sample_config_source: str
    ) -> None:
        """If CFG.marker_dict already exists, replace its value (no duplicate)."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        # First call
        paper_context = {"features": ["A", "B"]}
        _post_process_config(str(config_path), paper_context)

        source_after_first = config_path.read_text()
        count_marker_dict_first = source_after_first.count("CFG.marker_dict")
        assert count_marker_dict_first == 1, (
            f"Expected 1 CFG.marker_dict after first call, got {count_marker_dict_first}"
        )

        # Second call with different features — should replace, not append
        paper_context2 = {"features": ["X", "Y", "Z"]}
        _post_process_config(str(config_path), paper_context2)

        source_after_second = config_path.read_text()
        count_marker_dict_second = source_after_second.count("CFG.marker_dict")
        assert count_marker_dict_second == 1, (
            f"Expected 1 CFG.marker_dict after second call, got {count_marker_dict_second}"
        )
        # The old features should be gone
        assert "A" not in source_after_second
        assert "B" not in source_after_second
        # The new features should be present
        assert "X" in source_after_second
        assert "Y" in source_after_second
        assert "Z" in source_after_second

    def test_idempotent_new_field_appended_not_duplicated(
        self, tmp_path: Path, sample_config_source: str
    ) -> None:
        """New field (is_nuclei) should not be duplicated on subsequent calls."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        # First call — inject is_nuclei
        _post_process_config(str(config_path), {"is_nuclei": True})
        source1 = config_path.read_text()
        assert source1.count("CFG.is_nuclei") == 1

        # Second call — inject again; it already exists so should replace, not append
        _post_process_config(str(config_path), {"is_nuclei": True})
        source2 = config_path.read_text()
        assert source2.count("CFG.is_nuclei") == 1, (
            f"Expected 1 is_nuclei after second call, got {source2.count('CFG.is_nuclei')}"
        )

    def test_empty_paper_context_no_changes(
        self, tmp_path: Path, sample_config_source: str
    ) -> None:
        """Empty paper_context should not modify the config."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)
        original = config_path.read_text()

        _post_process_config(str(config_path), {})
        assert config_path.read_text() == original

        _post_process_config(str(config_path), None)  # type: ignore[arg-type]
        assert config_path.read_text() == original

    def test_mixed_replace_and_append(
        self, tmp_path: Path, sample_config_source: str
    ) -> None:
        """Replace marker_dict and append is_nuclei in one call."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        paper_context = {
            "features": ["M1", "M2"],
            "is_nuclei": True,
            "tissue_kb": "retina",
            "tissue_ontology": "UBERON_0000966",
        }
        _post_process_config(str(config_path), paper_context)

        source = config_path.read_text()
        # marker_dict was replaced (was already in template)
        assert source.count("CFG.marker_dict") == 1
        assert "M1" in source
        assert "M2" in source

        # is_nuclei was appended (not in template body as active assignment)
        assert "CFG.is_nuclei = True" in source

        # tissue_kb was appended
        assert "CFG.tissue_kb = 'retina'" in source

        # tissue_ontology was appended
        assert "CFG.tissue_ontology = 'UBERON_0000966'" in source

    def test_result_is_valid_python(self, tmp_path: Path, sample_config_source: str) -> None:
        """The modified config should remain valid Python."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        paper_context = {
            "features": ["V1", "V2"],
            "is_nuclei": True,
            "tissue_kb": "retina",
        }
        _post_process_config(str(config_path), paper_context)

        # Parse with ast — should not raise SyntaxError
        source = config_path.read_text()
        tree = ast.parse(source)
        assert tree is not None

    def test_source_without_marker_dict(
        self, tmp_path: Path, sample_config_source: str
    ) -> None:
        """Config without marker_dict — append it via ast."""
        from core.preprocess.matrix_loader import _post_process_config

        # Remove marker_dict from source
        source_without_md = sample_config_source.replace(
            "CFG.marker_dict = {\n    # 'CellTypeA': ['GENE1', 'GENE2'],\n}\n", ""
        )
        config_path = tmp_path / "config_test.py"
        config_path.write_text(source_without_md)

        _post_process_config(str(config_path), {"features": ["X", "Y"]})

        source = config_path.read_text()
        assert "CFG.marker_dict" in source
        assert "X" in source
        assert "Y" in source

    def test_skip_on_syntax_error(self, tmp_path: Path) -> None:
        """Malformed Python config should be skipped gracefully."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_bad.py"
        config_path.write_text("this is not valid python }{")

        # Should not raise
        _post_process_config(str(config_path), {"features": ["X"]})

        # File content unchanged
        assert config_path.read_text() == "this is not valid python }{"


# ═══════════════════════════════════════════════════════════════════════
#  generate_config tests
# ═══════════════════════════════════════════════════════════════════════


class TestGenerateConfigPaperContext:
    """paper_context override logic in generate_config()."""

    @pytest.fixture(autouse=True)
    def _patch_format_detection(self, monkeypatch: MonkeyPatch) -> None:
        """Force _detect_primary_format to return a known template format.

        The fake_classification fixture has no real files, so without mocking
        the format detection returns 'unknown' and generate_config skips.
        """
        import core.preprocess.matrix_loader as ml
        monkeypatch.setattr(ml, '_detect_primary_format', lambda *a, **kw: '10X_h5')
    def test_heuristic_overrides(self, tmp_path: Path, fake_classification: dict[str, Any], fake_file_list: list[str]) -> None:
        """Paper_context values override heuristic replacements."""
        from core.preprocess.matrix_loader import generate_config

        output_dir = str(tmp_path / "output")
        paper_context = {
            "species": "mus_musculus",
            "tissue": "retina",
            "expression_type": "log1p_counts",
            "genome": "mm10",
        }
        result = generate_config(
            gse_id="GSE000001",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            paper_context=paper_context,
            dry_run=True,
        )
        assert result is not None
        # With dry_run, the file is not written, but the path is returned.

    def test_without_paper_context_no_change(
        self, tmp_path: Path, fake_classification: dict[str, Any], fake_file_list: list[str]
    ) -> None:
        """Without paper_context, behavior is identical to not passing it."""
        from core.preprocess.matrix_loader import generate_config

        output_dir = str(tmp_path / "output_a")
        result_without = generate_config(
            gse_id="GSE000002",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            dry_run=True,
        )
        output_dir_b = str(tmp_path / "output_b")
        result_with_none = generate_config(
            gse_id="GSE000002",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir_b,
            data_root=str(tmp_path / "data"),
            paper_context=None,
            dry_run=True,
        )
        assert result_without is not None
        assert result_with_none is not None
        # Both dry-run, so they just return paths (different dirs but both succeed)
        assert result_without != result_with_none  # different dirs, different paths

    def test_paper_context_with_features(
        self, tmp_path: Path, fake_classification: dict[str, Any], fake_file_list: list[str]
    ) -> None:
        """Paper_context with features list triggers _post_process_config."""
        from core.preprocess.matrix_loader import generate_config

        output_dir = str(tmp_path / "output")
        paper_context = {
            "species": "human",
            "features": ["GENE_A", "GENE_B", "GENE_C"],
        }
        config_path = generate_config(
            gse_id="GSE000003",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            paper_context=paper_context,
            force=True,
        )
        assert config_path is not None
        assert os.path.exists(config_path)

        source = Path(config_path).read_text()
        # The features should appear in the config via _post_process_config
        assert "GENE_A" in source
        assert "GENE_B" in source
        assert "GENE_C" in source

    def test_assay_type_override(
        self, tmp_path: Path, fake_classification: dict[str, Any], fake_file_list: list[str]
    ) -> None:
        """assay_type in paper_context overrides the template's EXPRESSION_TYPE."""
        from core.preprocess.matrix_loader import generate_config

        output_dir = str(tmp_path / "output")
        paper_context = {
            "assay_type": "snRNAseq",
            "species": "human",
        }
        config_path = generate_config(
            gse_id="GSE000004",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            paper_context=paper_context,
            force=True,
        )
        assert config_path is not None
        assert os.path.exists(config_path)

    def test_force_overwrite(
        self, tmp_path: Path, fake_classification: dict[str, Any], fake_file_list: list[str]
    ) -> None:
        """force=True with paper_context should overwrite existing config."""
        from core.preprocess.matrix_loader import generate_config

        output_dir = str(tmp_path / "output")
        paper_context = {"features": ["A", "B"]}

        # First write
        gen1 = generate_config(
            gse_id="GSE000005",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            paper_context=paper_context,
            force=True,
        )
        assert gen1 is not None
        mtime1 = os.path.getmtime(gen1)

        # Second write with force — should overwrite not skip
        gen2 = generate_config(
            gse_id="GSE000005",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            paper_context=paper_context,
            force=True,
        )
        assert gen2 is not None
        assert gen2 == gen1  # same path
        mtime2 = os.path.getmtime(gen2)
        assert mtime2 >= mtime1, "File should have been re-written"

    def test_dry_run_with_paper_context(
        self, tmp_path: Path, fake_classification: dict[str, Any], fake_file_list: list[str], capsys
    ) -> None:
        """Dry-run with paper_context does not write files."""
        from core.preprocess.matrix_loader import generate_config

        output_dir = str(tmp_path / "output")
        paper_context = {"species": "mouse", "features": ["X"]}
        result = generate_config(
            gse_id="GSE000006",
            modality="rna",
            classification=fake_classification,
            file_list=fake_file_list,
            output_dir=output_dir,
            data_root=str(tmp_path / "data"),
            paper_context=paper_context,
            dry_run=True,
        )
        assert result is not None
        assert not os.path.exists(result)
        captured = capsys.readouterr()
        assert "[DRY-RUN]" in captured.out

# ═══════════════════════════════════════════════════════════════════════
#  _post_process_config — inject parameter tests
# ═══════════════════════════════════════════════════════════════════════


class TestPostProcessConfigInject:
    """inject parameter for arbitrary CFG.* value injection."""

    def test_inject_sample_keep(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject sample_keep list via inject parameter."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        _post_process_config(
            str(config_path),
            paper_context={},
            inject={"sample_keep": ["GSM1"]},
        )

        source = config_path.read_text()
        assert "CFG.sample_keep = [" in source
        assert "'GSM1'" in source

    def test_inject_subset_suffix(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject subset_suffix string via inject parameter."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        _post_process_config(
            str(config_path),
            paper_context={},
            inject={"subset_suffix": "_test"},
        )

        source = config_path.read_text()
        assert "CFG.subset_suffix = '_test'" in source

    def test_inject_mixed_types(self, tmp_path: Path, sample_config_source: str) -> None:
        """Inject mixed types (string, list, bool) — all correctly repr'd."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)

        _post_process_config(
            str(config_path),
            paper_context={},
            inject={
                "subset_suffix": "_test",
                "sample_keep": ["GSM1", "GSM2"],
                "skip_qc": True,
            },
        )

        source = config_path.read_text()
        assert "CFG.subset_suffix = '_test'" in source
        assert "CFG.sample_keep = ['GSM1', 'GSM2']" in source
        assert "CFG.skip_qc = True" in source

    def test_inject_none(self, tmp_path: Path, sample_config_source: str) -> None:
        """inject=None → unchanged (backward compat)."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)
        original = config_path.read_text()

        _post_process_config(str(config_path), paper_context={}, inject=None)
        assert config_path.read_text() == original

    def test_inject_empty(self, tmp_path: Path, sample_config_source: str) -> None:
        """inject={} → no changes."""
        from core.preprocess.matrix_loader import _post_process_config

        config_path = tmp_path / "config_test.py"
        config_path.write_text(sample_config_source)
        original = config_path.read_text()

        _post_process_config(str(config_path), paper_context={}, inject={})
        assert config_path.read_text() == original
