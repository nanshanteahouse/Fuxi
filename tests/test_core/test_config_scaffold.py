"""Regression tests for the schema-driven config scaffold.

Guards the "spec drifted from schema" class of bugs:

1. Spec ↔ schema consistency (``validate_specs``) — schema renames or
   deleted fields fail immediately; free-form keys below dict-typed
   fields (e.g. marker_dict entries) stay allowed.
2. Every rendered starter config is valid YAML (no duplicate keys,
   which used to silently corrupt nested blocks) and resolves through
   ``resolve_config`` — rendering happens on demand, so the rendered
   text itself is the artifact under test.
3. ``generate_config`` end-to-end across modalities/formats, including the
   previously unsupported visium / preprocessed / unknown / multiome paths.
"""

from __future__ import annotations

import io
import pathlib

import pytest
import yaml

from core.preprocess.config_specs import (
    _path_is_valid,
    materialized_specs,
    validate_specs,
)
from core.utils._config import resolve_config

# ═══════════════════════════════════════════════════════════════════
# 1 — Spec ↔ schema consistency
# ═══════════════════════════════════════════════════════════════════


class TestValidateSpecs:
    def test_specs_consistent_with_schema(self) -> None:
        assert validate_specs() == []

    @pytest.mark.parametrize(
        "path, expected",
        [
            ("qc.min_genes", True),
            ("integration.scvi.batch_size", True),
            ("integration.scvi.datasplitter_kwargs.num_workers", True),  # dict → free-form
            ("marker.marker_dict.RPCs", True),  # dict → free-form
            ("qc.not_a_field", False),
            ("clustering.leiden_resolutions", False),  # ghost field removed in 2026
        ],
    )
    def test_path_validity(self, path: str, expected: bool) -> None:
        assert _path_is_valid(path) is expected


# ═══════════════════════════════════════════════════════════════════
# 2 — Rendered starter configs are valid and resolvable
# ═══════════════════════════════════════════════════════════════════


def _duplicate_keys(text: str) -> list[str]:
    """Return mapping keys that appear twice in the same YAML mapping."""
    dups: list[str] = []
    stack: list[list] = []  # [known_keys: set[str], expect_key: bool]
    seq_depth = 0
    for ev in yaml.parse(io.StringIO(text)):
        if isinstance(ev, yaml.events.MappingStartEvent):
            stack.append([set(), True])
        elif isinstance(ev, yaml.events.MappingEndEvent):
            stack.pop()
        elif isinstance(ev, yaml.events.SequenceStartEvent):
            if seq_depth == 0 and stack:
                stack[-1][1] = True  # the sequence consumed the value slot
            seq_depth += 1
        elif isinstance(ev, yaml.events.SequenceEndEvent):
            seq_depth -= 1
        elif isinstance(ev, yaml.events.ScalarEvent) and seq_depth == 0:
            if stack and stack[-1][1]:
                if ev.value in stack[-1][0]:
                    dups.append(ev.value)
                stack[-1][0].add(ev.value)
            if stack:
                stack[-1][1] = not stack[-1][1]
    return dups


class TestRenderedTemplates:
    @pytest.mark.parametrize("name", [s.key for s in materialized_specs()])
    def test_spec_renders_to_valid_yaml(self, name: str) -> None:
        from core.config.scaffold import render_template_text

        spec = next(s for s in materialized_specs() if s.key == name)
        data = yaml.safe_load(render_template_text(spec))
        assert isinstance(data, dict)
        assert data["data_format"] == spec.data_format

    @pytest.mark.parametrize("name", [s.key for s in materialized_specs()])
    def test_spec_render_has_no_duplicate_keys(self, name: str) -> None:
        from core.config.scaffold import render_template_text

        spec = next(s for s in materialized_specs() if s.key == name)
        assert _duplicate_keys(render_template_text(spec)) == []

    @pytest.mark.parametrize("name", [s.key for s in materialized_specs()])
    def test_rendered_spec_resolves(self, tmp_path: pathlib.Path, name: str) -> None:
        from core.config.scaffold import render_template_text

        spec = next(s for s in materialized_specs() if s.key == name)
        dst = tmp_path / f"{name}.yaml"
        dst.write_text(render_template_text(spec), encoding="utf-8")
        cfg = resolve_config(str(dst))
        assert cfg.data_format
        assert cfg.modality


# ═══════════════════════════════════════════════════════════════════
# 3 — scaffold CLI
# ═══════════════════════════════════════════════════════════════════


class TestScaffoldCli:
    def test_list_lists_all_formats(self, capsys: pytest.CaptureFixture) -> None:
        from core.config.scaffold import main

        assert main(["--list"]) == 0
        out = capsys.readouterr().out
        for spec in materialized_specs():
            assert spec.key in out

    def test_unknown_format_fails(self, capsys: pytest.CaptureFixture) -> None:
        from core.config.scaffold import main

        assert main(["--format", "not_a_format"]) == 1
        assert "not_a_format" in capsys.readouterr().err

    def test_out_requires_format(self) -> None:
        from core.config.scaffold import main

        assert main(["--out", "config.yaml"]) == 1

    def test_format_out_writes_file(self, tmp_path: pathlib.Path) -> None:
        from core.config.scaffold import main

        out = tmp_path / "config.yaml"
        assert main(["--format", "10X_h5", "--out", str(out)]) == 0
        assert yaml.safe_load(out.read_text(encoding="utf-8"))["data_format"] == "10X_h5"


# ═══════════════════════════════════════════════════════════════════
# 4 — generate_config end-to-end
# ═══════════════════════════════════════════════════════════════════


def _make_files(tmp_path: pathlib.Path, paths: list[str]) -> list[str]:
    out: list[str] = []
    for p in paths:
        full = tmp_path / p
        full.parent.mkdir(parents=True, exist_ok=True)
        full.touch()
        out.append(str(full))
    return out


class TestGenerateConfigEndToEnd:
    @pytest.mark.parametrize(
        "modality, classification_files, expect_format",
        [
            ("rna", {"tenx_h5_dirs": {"s1": ["filtered_feature_bc_matrix.h5"]}}, "10X_h5"),
            (
                "rna",
                {"tenx_mtx_dirs": {"s1": ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"]}},
                "10X_mtx",
            ),
            (
                "rna",
                {"csv_files": ["counts.csv"], "metadata_files": ["barcodes.tsv", "features.tsv"]},
                "csv_matrix",
            ),
            ("rna", {"h5ad_files": ["data.h5ad"]}, "10X_h5"),  # reuse
            ("rna", {"preprocessed_dirs": {"s1": ["preprocessed.tsv.gz"]}}, "preprocessed"),
            (
                "atac",
                {"fragment_dirs": {"s1": ["fragments.tsv.gz", "barcodes.tsv.gz"]}},
                "10x_fragments",
            ),
            ("atac", {"tenx_peak_dirs": {"s1": ["peaks.bed"]}}, "10x_fragments"),  # reuse
            ("spatial", {"visium_dirs": {"s1": ["filtered_feature_bc_matrix.h5"]}}, "visium"),
            ("bulk", {"csv_files": ["counts.csv"]}, "csv_matrix"),  # cross-modal fallback
            ("bulk", {"h5ad_files": ["bulk.h5ad"]}, "10X_h5"),  # cross-modal fallback
            ("rna", {}, "10X_h5"),  # unknown → rna default
            ("spatial", {}, "visium"),  # unknown → spatial default
            ("multiome", {"tenx_h5_dirs": {"s1": ["a.h5"]}}, "10X_h5"),  # inherits rna
            ("multiome", {"fragment_dirs": {"s1": ["fragments.tsv.gz"]}}, "10x_fragments"),
            ("bulk", {}, "count_matrix"),  # unknown → bulk default
        ],
    )
    def test_generate_config(
        self,
        tmp_path: pathlib.Path,
        modality: str,
        classification_files: dict,
        expect_format: str,
    ) -> None:
        from core.preprocess.matrix_loader import generate_config

        classification: dict = {}
        file_list: list[str] = []
        for key, value in classification_files.items():
            if isinstance(value, dict):
                resolved: dict = {}
                for subdir, files in value.items():
                    real = _make_files(tmp_path, [f"{subdir}/{f}" for f in files])
                    resolved[str(tmp_path / subdir)] = real
                    file_list.extend(real)
                classification[key] = resolved
            else:
                real = _make_files(tmp_path, value)
                classification[key] = real
                file_list.extend(real)

        out_dir = tmp_path / "out"
        result = generate_config(
            "GSE99999",
            modality,
            classification,
            file_list,
            output_dir=str(out_dir),
            data_root=str(tmp_path),
        )
        assert result is not None, (
            f"generate_config returned None for {modality}/{classification_files}"
        )
        assert pathlib.Path(result).is_file()

        with open(result, encoding="utf-8") as f:
            data = yaml.safe_load(f)
        assert data["data_format"] == expect_format

        cfg = resolve_config(result)
        assert cfg.data_format == expect_format

    def test_paper_context_overrides(self, tmp_path: pathlib.Path) -> None:
        from core.preprocess.matrix_loader import generate_config

        files = _make_files(
            tmp_path, ["s1/matrix.mtx.gz", "s1/barcodes.tsv.gz", "s1/features.tsv.gz"]
        )
        classification = {"tenx_mtx_dirs": {str(tmp_path / "s1"): files}}
        result = generate_config(
            "GSE88888",
            "rna",
            classification,
            files,
            output_dir=str(tmp_path / "out"),
            data_root=str(tmp_path),
            paper_context={
                "tissue": "retina",
                "species": "mus_musculus",
                "is_nuclei": True,
                "features": ["RHO", "GNAT1"],
            },
        )
        assert result is not None
        with open(result, encoding="utf-8") as f:
            data = yaml.safe_load(f)
        assert data["tissue"] == "retina"
        assert data["species"] == "mus_musculus"
        assert data["qc"]["is_nuclei"] is True
        assert data["marker"]["marker_dict"] == {"extracted": ["RHO", "GNAT1"]}
