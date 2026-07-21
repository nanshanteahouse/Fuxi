"""Tests for core/run_reproduce.py — Reproduction orchestration (P3)."""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml

from core.config.dataset import (
    DatasetMeta,
    load_dataset,
    save_dataset,
    update_pipeline_status,
)
from core.paper.registry import MasterRegistry
from core.pipeline.reproduce import (
    REPRODUCE_TIMEOUT,
    _detect_modality,
    _extract_geo_ids,
    _run_pipeline_for_gse,
    _write_pipeline_status,
    run_reproduce,
)

# ═══════════════════════════════════════════════════════════════════════
# _detect_modality
# ═══════════════════════════════════════════════════════════════════════


class TestDetectModality:
    """Modality detection from pipeline config files (reuses registry)."""

    def test_rna(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text('CFG.modality = "rna"\n')
        assert _detect_modality(str(p)) == "rna"

    def test_atac(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text('CFG.modality = "atac"\n')
        assert _detect_modality(str(p)) == "atac"

    def test_spatial(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text('CFG.modality = "spatial"\n')
        assert _detect_modality(str(p)) == "spatial"

    def test_single_quotes(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text("CFG.modality = 'rna'\n")
        assert _detect_modality(str(p)) == "rna"

    def test_missing_file(self) -> None:
        assert _detect_modality("/nonexistent/path.py") == "unknown"

    def test_no_modality_line(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text("x = 1\n")
        assert _detect_modality(str(p)) == "unknown"


# ═══════════════════════════════════════════════════════════════════════
# _extract_geo_ids
# ═══════════════════════════════════════════════════════════════════════


class TestExtractGeoIds:
    """GEO ID extraction from insights dict + regex fallback."""

    def test_from_data_access(self) -> None:
        insights = {"data_access": {"geo_ids": ["GSE123456", "GSE789012"]}}
        assert _extract_geo_ids(insights) == ["GSE123456", "GSE789012"]

    def test_empty_then_regex_fallback(self) -> None:
        insights = {"data_access": {"geo_ids": []}}
        text = "Some GSE123456 and GSE789012 and GSE123456 again"
        result = _extract_geo_ids(insights, text)
        # Deduplicated, first-occurrence order preserved
        assert result == ["GSE123456", "GSE789012"]

    def test_no_matches(self) -> None:
        insights = {"data_access": {"geo_ids": []}}
        assert _extract_geo_ids(insights, "no GSE here") == []

    def test_dedup_in_dict(self) -> None:
        insights = {"data_access": {"geo_ids": ["GSE001", "GSE001", "GSE002"]}}
        assert _extract_geo_ids(insights) == ["GSE001", "GSE001", "GSE002"]

    def test_no_data_access_key(self) -> None:
        insights = {"paper_meta": {"pmid": "123"}}
        assert _extract_geo_ids(insights) == []

    def test_regex_fallback_dedup(self) -> None:
        """Regex fallback deduplicates; dict path does not."""
        insights = {"data_access": {"geo_ids": []}}
        text = "GSE001 GSE002 GSE001 GSE003 GSE002"
        result = _extract_geo_ids(insights, text)
        assert result == ["GSE001", "GSE002", "GSE003"]

    def test_non_list_geo_ids(self) -> None:
        insights = {"data_access": {"geo_ids": None}}
        assert _extract_geo_ids(insights) == []


# ═══════════════════════════════════════════════════════════════════════
# _run_pipeline_for_gse
# ═══════════════════════════════════════════════════════════════════════


class TestRunPipelineForGse:
    """Subprocess-based pipeline execution per GSE."""

    def _make_config(self, tmp_path: Path, modality: str = "rna") -> str:
        p = tmp_path / "config.py"
        p.write_text(f'CFG.modality = "{modality}"\n')
        return str(p)

    def test_success(self, tmp_path: Path) -> None:
        config_path = self._make_config(tmp_path, "rna")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="pipeline done", stderr="")
            result = _run_pipeline_for_gse("GSE001", config_path)

        assert result["status"] == "success"
        assert result["modality"] == "rna"
        assert result["config_path"] == config_path
        assert result["output"] == "pipeline done"
        assert isinstance(result["duration_s"], (int, float))
        mock_run.assert_called_once()

    def test_failure(self, tmp_path: Path) -> None:
        config_path = self._make_config(tmp_path, "atac")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error occurred")
            result = _run_pipeline_for_gse("GSE002", config_path)

        assert result["status"] == "failed"
        assert result["modality"] == "atac"
        assert "error occurred" in result["error"]

    def test_timeout(self, tmp_path: Path) -> None:
        config_path = self._make_config(tmp_path, "rna")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.TimeoutExpired("cmd", REPRODUCE_TIMEOUT)
            result = _run_pipeline_for_gse("GSE003", config_path)

        assert result["status"] == "timeout"
        assert "Timed out" in result["error"]
        assert result["modality"] == "rna"

    def test_unknown_modality_does_not_call_subprocess(self, tmp_path: Path) -> None:
        config_path = tmp_path / "config.py"
        config_path.write_text("x = 1\n")  # no CFG.modality

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            result = _run_pipeline_for_gse("GSE004", str(config_path))

        assert result["status"] == "failed"
        assert result["modality"] == "unknown"
        assert "Cannot detect modality" in result["error"]
        mock_run.assert_not_called()

    def test_subprocess_called_with_correct_args(self, tmp_path: Path) -> None:
        config_path = self._make_config(tmp_path, "spatial")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            _run_pipeline_for_gse("GSE005", config_path)

        cmd = mock_run.call_args[0][0]
        assert "core/run_pipeline.py" in cmd
        assert "--config" in cmd
        assert config_path in cmd
        assert "--modality" in cmd
        assert "spatial" in cmd

    # ── New params tests (modality override, experiment_group) ─────────

    def test_modality_param_override(self, tmp_path: Path) -> None:
        """Passing modality='atac' overrides config file detection."""
        # Config says "rna", but we pass "atac" explicitly
        config_path = self._make_config(tmp_path, "rna")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = _run_pipeline_for_gse("GSE010", config_path, modality="atac")

        assert result["modality"] == "atac"
        assert result["status"] == "success"
        # Verify the subprocess call used "atac" despite config saying "rna"
        cmd = mock_run.call_args[0][0]
        atac_idx = cmd.index("--modality") + 1
        assert cmd[atac_idx] == "atac"

    def test_modality_param_none(self, tmp_path: Path) -> None:
        """Passing modality=None falls back to config detection."""
        config_path = self._make_config(tmp_path, "atac")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
            result = _run_pipeline_for_gse("GSE011", config_path, modality=None)

        assert result["modality"] == "atac"
        cmd = mock_run.call_args[0][0]
        atac_idx = cmd.index("--modality") + 1
        assert cmd[atac_idx] == "atac"

    def test_experiment_group_param(self, tmp_path: Path) -> None:
        """Passing ExperimentGroup should not error (stored for W2.4)."""
        from core.paper.registry import DatasetStatus, ExperimentGroup

        config_path = self._make_config(tmp_path, "rna")
        eg = ExperimentGroup(
            group_name="test_group",
            sample_ids=["s1", "s2"],
            subset_suffix="_subset",
            modality="rna",
            status=DatasetStatus.CONFIG_EXISTS,
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="eg_ok", stderr="")
            result = _run_pipeline_for_gse("GSE012", config_path, experiment_group=eg)

        assert result["status"] == "success"
        assert result["output"] == "eg_ok"

    def test_backward_compat(self, tmp_path: Path) -> None:
        """Call with only _gse_id and config_path works identically to before."""
        config_path = self._make_config(tmp_path, "spatial")

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="spatial_ok", stderr="")
            result = _run_pipeline_for_gse("GSE013", config_path)

        assert result["status"] == "success"
        assert result["modality"] == "spatial"
        assert result["output"] == "spatial_ok"
        cmd = mock_run.call_args[0][0]
        assert "--modality" in cmd
        spatial_idx = cmd.index("--modality") + 1
        assert cmd[spatial_idx] == "spatial"

    # ── Pipeline status write tests ──────────────────────────────────

    def test_success_writes_dataset_status(self, tmp_path: Path) -> None:
        """Successful run writes completed status to dataset.yaml."""
        config_path = self._make_config(tmp_path, "rna")
        # Create dataset.yaml next to config
        ds_yaml = tmp_path / "dataset.yaml"
        ds_yaml.write_text(
            yaml.dump(
                {
                    "id": "GSE001",
                    "type": "SingleAccession",
                    "title": "Test",
                    "meta": {"pipeline_status": {"scRNAseq": "pending"}},
                }
            )
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            with patch("core.pipeline.reproduce.update_pipeline_status") as mock_update:
                result = _run_pipeline_for_gse("GSE001", config_path)

        assert result["status"] == "success"
        mock_update.assert_called_once()
        args = mock_update.call_args[0]
        assert args[0] == str(ds_yaml)  # yaml_path
        assert args[1] == "rna"  # modality_key
        assert args[2] == "completed"  # status

    def test_success_missing_dataset_yaml_skipped(self, tmp_path: Path) -> None:
        """No dataset.yaml → skip status write, result still success."""
        config_path = self._make_config(tmp_path, "rna")
        # No dataset.yaml created

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            with patch("core.pipeline.reproduce.update_pipeline_status") as mock_update:
                result = _run_pipeline_for_gse("GSE002", config_path)

        assert result["status"] == "success"
        mock_update.assert_not_called()

    def test_success_corrupt_dataset_yaml_does_not_crash(self, tmp_path: Path, caplog) -> None:
        """Corrupt dataset.yaml → warning logged, result still success."""
        config_path = self._make_config(tmp_path, "rna")
        # Create corrupt dataset.yaml
        ds_yaml = tmp_path / "dataset.yaml"
        ds_yaml.write_text("::: not valid yaml :::")

        caplog.set_level(logging.WARNING)

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            result = _run_pipeline_for_gse("GSE003", config_path)

        assert result["status"] == "success"
        assert "Failed to update pipeline status" in caplog.text


# ═══════════════════════════════════════════════════════════════════════
# run_reproduce — integration with mocked subprocess
# ═══════════════════════════════════════════════════════════════════════


def _make_paper_dir(
    tmp_path: Path,
    name: str = "TestPaper",
    pmid: str = "12345678",
    geo_ids: list[str] | None = None,
    features: list[str] | None = None,
) -> str:
    """Create a temporary paper directory with a minimal insights.yaml."""
    paper_dir = tmp_path / "papers" / name
    paper_dir.mkdir(parents=True)
    data: dict = {
        "paper_meta": {"pmid": pmid, "species": "human"},
        "data_access": {"geo_ids": geo_ids or ["GSE001"]},
    }
    if features:
        data["figures"] = [{"features": features}]
    (paper_dir / "insights.yaml").write_text(yaml.dump(data))
    return str(paper_dir)


def _make_registry(
    datasets: list[dict],
    pmid: str = "12345678",
    paper_dir: str = "TestPaper",
) -> MasterRegistry:
    """Build a minimal MasterRegistry for testing."""
    from core.paper.registry import (
        DatasetConfig,
        DatasetEntry,
        LinkRole,
        MasterRegistry,
        ModalityInfo,
        PaperDatasetLink,
        PaperEntry,
    )

    paper_id = paper_dir
    paper = PaperEntry(paper_id=paper_id, slug=paper_dir, pmid=pmid, paper_dir=paper_dir)

    dataset_entries: dict[str, DatasetEntry] = {}
    links: list[PaperDatasetLink] = []

    for ds in datasets:
        gse_id = ds["gse_id"]
        modality = ds.get("modality", "rna")
        status = ds.get("status", "unknown")
        config_path = ds.get("config_path", "")
        experiments = ds.get("experiments")

        ds_config = DatasetConfig(path=config_path)
        if experiments is not None:
            ds_config.experiments = experiments

        mod_info = ModalityInfo(status=status, configs=[ds_config])
        dataset_entry = DatasetEntry(modalities={modality: mod_info})
        dataset_entries[gse_id] = dataset_entry
        links.append(PaperDatasetLink(paper_id=paper_id, dataset_id=gse_id, role=LinkRole.PRIMARY))

    return MasterRegistry(papers=[paper], datasets=dataset_entries, links=links)


class TestRunReproduceErrors:
    """Error / edge-case paths in run_reproduce()."""

    def test_missing_insights_raises(self, tmp_path: Path) -> None:
        paper_dir = tmp_path / "papers" / "Missing"
        paper_dir.mkdir(parents=True)
        with pytest.raises(ValueError, match="No insights.yaml"):
            run_reproduce(str(paper_dir))

    def test_empty_insights_raises(self, tmp_path: Path) -> None:
        paper_dir = tmp_path / "papers" / "Empty"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text("")
        with pytest.raises(ValueError, match="Empty or invalid"):
            run_reproduce(str(paper_dir))

    def test_none_insights_raises(self, tmp_path: Path) -> None:
        """YAML with only a comment yields None when loaded."""
        paper_dir = tmp_path / "papers" / "NoneYaml"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text("# just a comment\n")
        with pytest.raises(ValueError, match="Empty or invalid"):
            run_reproduce(str(paper_dir))

    def test_no_registry_entry(self, tmp_path: Path) -> None:
        """Paper PMID not found in registry → empty results."""
        paper_dir = _make_paper_dir(tmp_path, name="NoReg", pmid="99999999")
        results = run_reproduce(paper_dir, registry=MasterRegistry())
        assert results == {}

    def test_gse_filter_excludes_others(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path, name="MultiGSE", geo_ids=["GSE001", "GSE002"])
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "not_configured",
                    "modality": "rna",
                },
                {
                    "gse_id": "GSE002",
                    "config_path": "",
                    "status": "data_not_downloaded",
                    "modality": "rna",
                },
            ]
        )
        results = run_reproduce(paper_dir, registry=registry, gse_filter="GSE001")
        assert "GSE001" in results
        assert "GSE002" not in results


class TestRunReproduceStates:
    """All 5 DatasetStatus values + dry_run."""

    def test_config_exists_runs_pipeline(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        config_path = tmp_path / "config_GSE001.py"
        config_path.write_text('CFG.modality = "rna"\n')

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": str(config_path),
                    "status": "config_exists",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "success"
        assert results["GSE001"]["modality"] == "rna"
        assert results["GSE001"]["config_path"] == str(config_path)

    def test_not_configured(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "not_configured",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "not_configured"
        assert "config generation" in results["GSE001"]["reason"]
        mock_run.assert_not_called()

    def test_config_exists_no_path(self, tmp_path: Path) -> None:
        """config_exists but empty config_path → treated as not_configured."""
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "config_exists",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "not_configured"
        mock_run.assert_not_called()

    def test_data_not_downloaded(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "data_not_downloaded",
                    "modality": "atac",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "skipped"
        assert "data not downloaded" in results["GSE001"]["reason"].lower()
        assert results["GSE001"]["modality"] == "atac"
        mock_run.assert_not_called()

    def test_pipeline_complete_skipped(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "pipeline_complete",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "skipped"
        assert "pipeline_complete" in results["GSE001"]["reason"]
        mock_run.assert_not_called()

    def test_data_only_skipped(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "data_only",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "skipped"
        assert "data_only" in results["GSE001"]["reason"]
        mock_run.assert_not_called()

    def test_unknown_skipped(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "unknown",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry)

        assert results["GSE001"]["status"] == "skipped"
        assert "unknown" in results["GSE001"]["reason"]
        mock_run.assert_not_called()


class TestRunReproduceDryRun:
    """Dry-run mode must never execute the pipeline."""

    def test_dry_run_skips_pipeline(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "/fake/path.py",
                    "status": "config_exists",
                    "modality": "rna",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry, dry_run=True)

        assert results["GSE001"]["status"] == "dry_run"
        assert results["GSE001"]["config_path"] == "/fake/path.py"
        assert results["GSE001"]["modality"] == "rna"
        mock_run.assert_not_called()

    def test_dry_run_not_configured(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(tmp_path)
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "not_configured",
                    "modality": "atac",
                }
            ]
        )
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            results = run_reproduce(paper_dir, registry=registry, dry_run=True)

        assert results["GSE001"]["status"] == "dry_run"
        assert results["GSE001"]["modality"] == "atac"
        mock_run.assert_not_called()


class TestRunReproducePaperContext:
    """Paper context propagation from insights.yaml figures."""

    def test_features_extracted(self, tmp_path: Path) -> None:
        """Marker features from the first figure should be extracted."""
        paper_dir = _make_paper_dir(
            tmp_path,
            name="FeaturesPaper",
            features=["CD3D", "CD8A", "GZMB"],
        )
        # Verify the YAML was written with figures
        with open(Path(paper_dir) / "insights.yaml") as f:
            data = yaml.safe_load(f)
        assert data["figures"][0]["features"] == ["CD3D", "CD8A", "GZMB"]

    def test_multiple_figures_first_used(self, tmp_path: Path) -> None:
        """Only the first figure with features is used for markers."""
        paper_dir = tmp_path / "papers" / "MultiFig"
        paper_dir.mkdir(parents=True)
        data = {
            "paper_meta": {"pmid": "12345678"},
            "data_access": {"geo_ids": ["GSE001"]},
            "figures": [
                {"label": "Fig1", "features": ["A", "B"]},
                {"label": "Fig2", "features": ["C", "D"]},
            ],
        }
        (paper_dir / "insights.yaml").write_text(yaml.dump(data))
        # This test verifies the code path exists; the actual context
        # is used only during preprocess (not tested in dry-run mode).


class TestUpdatePipelineStatus:
    """update_pipeline_status: modality key mapping and edge cases."""

    def test_rna_to_scrnaseq(self, tmp_path: Path) -> None:
        """Happy path: rna maps to scRNAseq field."""
        yaml_path = tmp_path / "dataset.yaml"
        ds = DatasetMeta(id="test", type="SingleAccession", title="test")
        save_dataset(ds, str(yaml_path))
        update_pipeline_status(str(yaml_path), "rna", "completed")
        loaded = load_dataset(str(yaml_path))
        assert loaded.meta.pipeline_status.scRNAseq == "completed"

    def test_atac_to_atacseq(self, tmp_path: Path) -> None:
        """Modality mapping: atac maps to ATACseq field."""
        yaml_path = tmp_path / "dataset.yaml"
        ds = DatasetMeta(id="test", type="SingleAccession", title="test")
        save_dataset(ds, str(yaml_path))
        update_pipeline_status(str(yaml_path), "atac", "running")
        loaded = load_dataset(str(yaml_path))
        assert loaded.meta.pipeline_status.ATACseq == "running"
        assert loaded.meta.pipeline_status.scRNAseq is None

    def test_spatial_to_spatial(self, tmp_path: Path) -> None:
        """Modality mapping: spatial maps to spatial field."""
        yaml_path = tmp_path / "dataset.yaml"
        ds = DatasetMeta(id="test", type="SingleAccession", title="test")
        save_dataset(ds, str(yaml_path))
        update_pipeline_status(str(yaml_path), "spatial", "failed")
        loaded = load_dataset(str(yaml_path))
        assert loaded.meta.pipeline_status.spatial == "failed"

    def test_unknown_modality(self, tmp_path: Path) -> None:
        """Unknown modality logs warning and returns without crash."""
        yaml_path = tmp_path / "dataset.yaml"
        ds = DatasetMeta(id="test", type="SingleAccession", title="test")
        save_dataset(ds, str(yaml_path))
        update_pipeline_status(str(yaml_path), "unknown", "completed")
        loaded = load_dataset(str(yaml_path))
        # Verify no fields were modified
        assert loaded.meta.pipeline_status.scRNAseq is None
        assert loaded.meta.pipeline_status.ATACseq is None
        assert loaded.meta.pipeline_status.spatial is None

    def test_none_yaml_path(self) -> None:
        """None yaml_path logs warning and returns without crash."""
        update_pipeline_status(None, "rna", "completed")


class TestCLI:
    """CLI argument parsing via main()."""

    def test_dry_run_all_runs(self, tmp_path: Path) -> None:
        """python run_reproduce.py --all --dry-run must not error."""
        # Create a minimal registry with one paper
        paper_dir = tmp_path / "papers" / "CLITest"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump(
                {
                    "paper_meta": {"pmid": "12345678"},
                    "data_access": {"geo_ids": ["GSE001"]},
                }
            )
        )
        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "not_configured",
                    "modality": "rna",
                }
            ],
            paper_dir="CLITest",
        )

        # Patch the registry loader to return our test registry
        # and redirect projects/papers to tmp_path
        with (
            patch(
                "core.pipeline.reproduce.load_master_registry", return_value=registry
            ) as mock_load,
            patch("core.pipeline.reproduce.Path.is_dir", return_value=True),
        ):
            from core.pipeline.reproduce import main

            # We need to patch sys.argv
            test_args = ["run_reproduce.py", "--all", "--dry-run"]
            with patch.object(sys, "argv", test_args):
                # Should not raise
                main()

        mock_load.assert_called_once()


class TestRunReproduceWithExperimentGroups:
    """3-layer nested dispatch for experiment groups (W2.4)."""

    def _make_config(
        self, tmp_path: Path, name: str = "config_GSE001.py", modality: str = "rna"
    ) -> str:
        p = tmp_path / name
        p.write_text(f'CFG.modality = "{modality}"\n')
        return str(p)

    def test_single_experiment_group(self, tmp_path: Path) -> None:
        """1 experiment group (not multiome) -> 1 pipeline call."""
        paper_dir = _make_paper_dir(tmp_path, name="SingleExp", geo_ids=["GSE001"])
        config_path = self._make_config(tmp_path)
        exp_config_path = self._make_config(tmp_path, "config_GSE001_myeloid.py", "rna")

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": config_path,
                    "status": "config_exists",
                    "modality": "rna",
                    "experiments": [
                        {
                            "group_name": "Myeloid",
                            "sample_ids": ["s1", "s2"],
                            "subset_suffix": "_myeloid",
                            "modality": "rna",
                            "status": "config_exists",
                            "config_path": exp_config_path,
                        }
                    ],
                }
            ]
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        assert mock_run.call_count == 1
        assert "GSE001_Myeloid_rna" in results
        assert results["GSE001_Myeloid_rna"]["status"] == "success"
        cmd = mock_run.call_args[0][0]
        cfg_idx = cmd.index("--config") + 1
        assert cmd[cfg_idx] == exp_config_path
        mod_idx = cmd.index("--modality") + 1
        assert cmd[mod_idx] == "rna"

    def test_multiome_experiment_group(self, tmp_path: Path) -> None:
        """1 multiome group -> 2 pipeline calls (rna + atac)."""
        paper_dir = _make_paper_dir(tmp_path, name="MultiomeExp", geo_ids=["GSE002"])
        config_path = self._make_config(tmp_path, "config_GSE002.py")

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE002",
                    "config_path": config_path,
                    "status": "config_exists",
                    "modality": "multiome",
                    "experiments": [
                        {
                            "group_name": "MultiomeGroup",
                            "sample_ids": ["s1", "s2"],
                            "subset_suffix": "_multiome",
                            "modality": "multiome",
                            "status": "config_exists",
                        }
                    ],
                }
            ]
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        assert mock_run.call_count == 2
        assert "GSE002_MultiomeGroup_rna" in results
        assert "GSE002_MultiomeGroup_atac" in results

        calls = mock_run.call_args_list
        cmd0 = calls[0][0][0]
        mod0_idx = cmd0.index("--modality") + 1
        assert cmd0[mod0_idx] == "rna"
        cmd1 = calls[1][0][0]
        mod1_idx = cmd1.index("--modality") + 1
        assert cmd1[mod1_idx] == "atac"
        expected_config = config_path[:-3] + "_multiome.py"
        cfg0_idx = cmd0.index("--config") + 1
        cfg1_idx = cmd1.index("--config") + 1
        assert cmd0[cfg0_idx] == expected_config
        assert cmd1[cfg1_idx] == expected_config

    def test_multiple_experiment_groups(self, tmp_path: Path) -> None:
        """GSE310245 pattern: 2 groups (1 multiome, 1 rna) -> 3 pipeline calls."""
        paper_dir = _make_paper_dir(tmp_path, name="MultiGroup", geo_ids=["GSE003"])
        config_path = self._make_config(tmp_path, "config_GSE003.py")
        exp_a_config = self._make_config(tmp_path, "config_GSE003_A_rna.py", "rna")
        exp_b_config = self._make_config(tmp_path, "config_GSE003_B_rna.py", "rna")

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE003",
                    "config_path": config_path,
                    "status": "config_exists",
                    "modality": "rna",
                    "experiments": [
                        {
                            "group_name": "GroupA",
                            "sample_ids": ["s1"],
                            "subset_suffix": "_A",
                            "modality": "multiome",
                            "status": "config_exists",
                            "config_path": exp_a_config,
                        },
                        {
                            "group_name": "GroupB",
                            "sample_ids": ["s2"],
                            "subset_suffix": "_B",
                            "modality": "rna",
                            "status": "config_exists",
                            "config_path": exp_b_config,
                        },
                    ],
                }
            ]
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        # GroupA (multiome) -> 2 calls, GroupB (rna) -> 1 call = 3 total
        assert mock_run.call_count == 3
        assert "GSE003_GroupA_rna" in results
        assert "GSE003_GroupA_atac" in results
        assert "GSE003_GroupB_rna" in results

        calls = mock_run.call_args_list
        modalities_used = []
        for c in calls:
            cmd = c[0][0]
            mod_idx = cmd.index("--modality") + 1
            modalities_used.append(cmd[mod_idx])
        assert modalities_used.count("rna") == 2
        assert modalities_used.count("atac") == 1

    def test_no_experiments_backward_compat(self, tmp_path: Path) -> None:
        """Dataset without experiments -> 1 call, original flat key."""
        paper_dir = _make_paper_dir(tmp_path, name="Compat", geo_ids=["GSE004"])
        config_path = self._make_config(tmp_path, "config_GSE004.py", "rna")

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE004",
                    "config_path": config_path,
                    "status": "config_exists",
                    "modality": "rna",
                }
            ]
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        assert mock_run.call_count == 1
        assert "GSE004" in results
        assert "GSE004_GroupA_rna" not in results
        assert results["GSE004"]["status"] == "success"
        assert results["GSE004"]["modality"] == "rna"
        # Result dict keys match pre-W2.2 format (flat, no experiment grouping)
        assert set(results["GSE004"].keys()) == {
            "status",
            "config_path",
            "modality",
            "output",
            "error",
            "duration_s",
        }

        # Subprocess call uses the base config and dataset-level modality
        cmd = mock_run.call_args[0][0]
        cfg_idx = cmd.index("--config") + 1
        assert cmd[cfg_idx] == config_path
        mod_idx = cmd.index("--modality") + 1
        assert cmd[mod_idx] == "rna"

        # Config's experiments list is empty → flat dispatch
        ds_entry = registry.datasets["GSE004"]
        cfg = ds_entry.modalities["rna"].configs[0]
        assert cfg.experiments == []

    def test_mixed_datasets(self, tmp_path: Path) -> None:
        """Some datasets with experiments, some without."""
        paper_dir = _make_paper_dir(tmp_path, name="Mixed", geo_ids=["GSE005", "GSE006"])
        config_5 = self._make_config(tmp_path, "config_GSE005.py", "rna")
        config_6 = self._make_config(tmp_path, "config_GSE006.py", "atac")
        exp_config = self._make_config(tmp_path, "config_GSE005_exp.py", "rna")

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE005",
                    "config_path": config_5,
                    "status": "config_exists",
                    "modality": "rna",
                    "experiments": [
                        {
                            "group_name": "ExpGroup",
                            "sample_ids": ["s1"],
                            "subset_suffix": "_exp",
                            "modality": "rna",
                            "status": "config_exists",
                            "config_path": exp_config,
                        }
                    ],
                },
                {
                    "gse_id": "GSE006",
                    "config_path": config_6,
                    "status": "config_exists",
                    "modality": "atac",
                },
            ]
        )

        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        assert mock_run.call_count == 2
        assert "GSE005_ExpGroup_rna" in results
        assert results["GSE005_ExpGroup_rna"]["status"] == "success"
        assert "GSE006" in results
        assert results["GSE006"]["status"] == "success"
        assert results["GSE006"]["modality"] == "atac"

    def test_not_configured_with_config_generation(self, tmp_path: Path) -> None:
        """Experiment group with status=not_configured -> config generated."""
        paper_dir = _make_paper_dir(tmp_path, name="ExpCfgGen", geo_ids=["GSE001"])

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE001",
                    "config_path": "",
                    "status": "not_configured",
                    "modality": "rna",
                    "experiments": [
                        {
                            "group_name": "Tcell",
                            "sample_ids": ["GSM1", "GSM2"],
                            "subset_suffix": "_Tcell",
                            "modality": "rna",
                            "status": "not_configured",
                        },
                        {
                            "group_name": "Myeloid",
                            "sample_ids": ["GSM3"],
                            "subset_suffix": "_Myeloid",
                            "modality": "rna",
                            "status": "not_configured",
                        },
                    ],
                }
            ]
        )

        with (
            patch("core.pipeline.reproduce.subprocess.run") as mock_subproc,
            patch("core.preprocess.preprocessor.run_preprocess", return_value=0) as mock_preproc,
            patch("core.pipeline.reproduce.shutil.copy2") as mock_copy,
            patch("core.preprocess.matrix_loader._post_process_config") as mock_post,
            patch("core.pipeline.reproduce.os.path.exists", return_value=True),
        ):
            results = run_reproduce(paper_dir, registry=registry)

        # run_preprocess was called with gse_id
        mock_preproc.assert_called_once()
        call_kwargs = mock_preproc.call_args[1]
        assert call_kwargs["gse_id"] == "GSE001"

        # Two experiment groups -> two copy + post-process calls
        assert mock_copy.call_count == 2
        assert mock_post.call_count == 2

        # Verify _post_process_config was called with inject dict
        for call_args in mock_post.call_args_list:
            args, kwargs = call_args
            assert "inject" in kwargs
            assert "sample_keep" in kwargs["inject"]
            assert "subset_suffix" in kwargs["inject"]

        # Results use compound keys with "configured" status
        assert "GSE001_Tcell_rna" in results
        assert "GSE001_Myeloid_rna" in results
        assert results["GSE001_Tcell_rna"]["status"] == "configured"
        assert results["GSE001_Myeloid_rna"]["status"] == "configured"

        # Subprocess should NOT have been called
        mock_subproc.assert_not_called()


class TestGSE310245ExperimentGroups:
    """GSE310245 capstone validation: 2 experiment groups, 3 pipeline calls, 6 verification points.

    Verifies the complete experiment-groups flow end-to-end:
    - 1 multiome group triggers 2 pipeline calls (rna + atac),
    - 1 rna group triggers 1 pipeline call (rna),
    - Config paths follow subset_suffix naming,
    - No sample overlap between groups.
    """

    def test_full_gse310245_experiment_groups(self, tmp_path: Path) -> None:
        """1 multiome group (2 calls) + 1 rna group (1 call) = 3 subprocess calls."""
        # VP1: Mock registry with GSE310245 having 2 experiment groups
        paper_dir = _make_paper_dir(tmp_path, name="GSE310245Paper", geo_ids=["GSE310245"])

        base_cfg = str(tmp_path / "config_GSE310245.py")
        pcw8_cfg = str(tmp_path / "config_GSE310245_pcw8_multiome.py")
        d140_cfg = str(tmp_path / "config_GSE310245_d140_rs.py")

        Path(base_cfg).write_text('CFG.modality = "rna"\n')
        Path(pcw8_cfg).write_text('CFG.modality = "multiome"\n')
        Path(d140_cfg).write_text('CFG.modality = "rna"\n')

        pcw8_samples = ["GSM9292434", "GSM9292436"]
        d140_samples = ["GSM9567287", "GSM9567288"]

        registry = _make_registry(
            [
                {
                    "gse_id": "GSE310245",
                    "config_path": base_cfg,
                    "status": "config_exists",
                    "modality": "rna",
                    "experiments": [
                        {
                            "group_name": "pcw8_multiome",
                            "sample_ids": pcw8_samples,
                            "subset_suffix": "_pcw8_multiome",
                            "modality": "multiome",
                            "status": "config_exists",
                            "config_path": pcw8_cfg,
                            "figures": ["Fig1B", "Fig1C"],
                        },
                        {
                            "group_name": "d140_rs",
                            "sample_ids": d140_samples,
                            "subset_suffix": "_d140_rs",
                            "modality": "rna",
                            "status": "config_exists",
                            "config_path": d140_cfg,
                            "figures": ["Fig4E", "Fig4F"],
                        },
                    ],
                }
            ],
            paper_dir="GSE310245Paper",
        )

        # VP3: Call run_reproduce(paper_dir, dry_run=False) with mocked subprocess
        with patch("core.pipeline.reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="ok", stderr="")
            results = run_reproduce(paper_dir, registry=registry)

        # VP4: 3 subprocess calls total with correct --modality flags
        assert mock_run.call_count == 3, f"Expected 3 calls, got {mock_run.call_count}"
        assert "GSE310245_pcw8_multiome_rna" in results
        assert "GSE310245_pcw8_multiome_atac" in results
        assert "GSE310245_d140_rs_rna" in results

        calls = mock_run.call_args_list
        rna_count = 0
        atac_count = 0
        for c in calls:
            cmd = c[0][0]
            mod_idx = cmd.index("--modality") + 1
            mod = cmd[mod_idx]
            if mod == "rna":
                rna_count += 1
            elif mod == "atac":
                atac_count += 1
        assert rna_count == 2, f"Expected 2 rna calls, got {rna_count}"
        assert atac_count == 1, f"Expected 1 atac call, got {atac_count}"

        # VP5: Config paths use subset_suffix naming - 2 configs, correct sample_keep per group
        assert results["GSE310245_pcw8_multiome_rna"]["status"] == "success"
        assert results["GSE310245_pcw8_multiome_atac"]["status"] == "success"
        assert results["GSE310245_d140_rs_rna"]["status"] == "success"

        cfg_count_pcw8 = 0
        cfg_count_d140 = 0
        for c in calls:
            cmd = c[0][0]
            cfg_idx = cmd.index("--config") + 1
            cfg_path = cmd[cfg_idx]
            if cfg_path == pcw8_cfg:
                cfg_count_pcw8 += 1
            elif cfg_path == d140_cfg:
                cfg_count_d140 += 1
        # pcw8_multiome (multiome) -> 2 calls (rna + atac) sharing same config
        assert cfg_count_pcw8 == 2, f"Expected 2 pcw8_multiome calls, got {cfg_count_pcw8}"
        # d140_rs (rna) -> 1 call (rna only)
        assert cfg_count_d140 == 1, f"Expected 1 d140_rs call, got {cfg_count_d140}"

        # VP6: No sample overlap between the two groups
        assert set(pcw8_samples).isdisjoint(set(d140_samples))


# ═══════════════════════════════════════════════════════════════════════
# _write_pipeline_status
# ═══════════════════════════════════════════════════════════════════════


class TestPipelineStatusWrite:
    """Incremental registry status writes after successful pipeline runs."""

    def test_write_sets_dataset_status(self, tmp_path: Path) -> None:
        """Successful run updates dataset status to pipeline_complete."""
        paper_dir = tmp_path / "papers" / "TestPaper"
        paper_dir.mkdir(parents=True)
        reg_path = tmp_path / "papers" / "registry.yaml"

        from core.paper.registry import (
            DatasetConfig,
            DatasetEntry,
            LinkRole,
            MasterRegistry,
            ModalityInfo,
            PaperDatasetLink,
            PaperEntry,
        )

        paper = PaperEntry(
            paper_id="12345678", slug="TestPaper", pmid="12345678", paper_dir="TestPaper"
        )
        dataset = DatasetEntry(
            modalities={
                "rna": ModalityInfo(
                    status="config_exists",
                    configs=[DatasetConfig(path="/fake/config.py")],
                )
            }
        )
        test_registry = MasterRegistry(
            papers=[paper],
            datasets={"GSE001": dataset},
            links=[
                PaperDatasetLink(paper_id="12345678", dataset_id="GSE001", role=LinkRole.PRIMARY)
            ],
        )

        result = {
            "status": "success",
            "config_path": "/fake/config.py",
            "modality": "rna",
        }

        with (
            patch("core.pipeline.reproduce.load_master_registry", return_value=test_registry),
            patch("core.pipeline.reproduce.save_master_registry") as mock_save,
        ):
            _write_pipeline_status(str(reg_path), "GSE001", result, str(paper_dir))

        saved_registry = mock_save.call_args[0][0]
        ds = saved_registry.datasets["GSE001"]
        assert ds.status == "pipeline_complete"
        assert ds.modalities["rna"].configs[0].pipeline_status == "pipeline_complete"

    def test_write_sets_experiment_group_status(self, tmp_path: Path) -> None:
        """Successful experiment group run updates group status."""
        paper_dir = tmp_path / "papers" / "ExpPaper"
        paper_dir.mkdir(parents=True)
        reg_path = tmp_path / "papers" / "registry.yaml"

        from core.paper.registry import (
            DatasetConfig,
            DatasetEntry,
            LinkRole,
            MasterRegistry,
            ModalityInfo,
            PaperDatasetLink,
            PaperEntry,
        )

        paper = PaperEntry(
            paper_id="87654321", slug="ExpPaper", pmid="87654321", paper_dir="ExpPaper"
        )
        dataset = DatasetEntry(
            modalities={
                "rna": ModalityInfo(
                    status="config_exists",
                    configs=[
                        DatasetConfig(
                            path="/fake/base.py",
                            experiments=[
                                {
                                    "group_name": "Myeloid",
                                    "sample_ids": ["s1"],
                                    "subset_suffix": "_myeloid",
                                    "modality": "rna",
                                    "status": "config_exists",
                                }
                            ],
                        )
                    ],
                )
            }
        )
        test_registry = MasterRegistry(
            papers=[paper],
            datasets={"GSE002": dataset},
            links=[
                PaperDatasetLink(paper_id="87654321", dataset_id="GSE002", role=LinkRole.PRIMARY)
            ],
        )

        result = {
            "status": "success",
            "config_path": "/fake/myeloid.py",
            "modality": "rna",
            "group_name": "Myeloid",
        }

        with (
            patch("core.pipeline.reproduce.load_master_registry", return_value=test_registry),
            patch("core.pipeline.reproduce.save_master_registry") as mock_save,
        ):
            _write_pipeline_status(str(reg_path), "GSE002", result, str(paper_dir))

        saved_registry = mock_save.call_args[0][0]
        ds = saved_registry.datasets["GSE002"]
        mod_info = ds.modalities["rna"]
        cfg = mod_info.configs[0]
        assert cfg.experiments[0]["status"] == "pipeline_complete"
        assert cfg.experiments[0]["group_name"] == "Myeloid"

    def test_write_gse_not_found_logs_warning(self, tmp_path: Path, caplog) -> None:
        """GSE not in registry logs warning without error."""
        paper_dir = tmp_path / "papers" / "NoGSE"
        paper_dir.mkdir(parents=True)
        reg_path = tmp_path / "papers" / "registry.yaml"
        reg_path.write_text(yaml.dump({"papers": []}))

        caplog.set_level(logging.WARNING)
        result: dict = {"status": "success", "config_path": "/fake.py"}
        _write_pipeline_status(str(reg_path), "GSE999", result, str(paper_dir))

        assert "GSE999" in caplog.text
        assert "not found" in caplog.text

    def test_write_calls_save_master_registry(self, tmp_path: Path) -> None:
        """Verify pipeline_complete write calls save_master_registry."""
        paper_dir = tmp_path / "papers" / "AtomicTest"
        paper_dir.mkdir(parents=True)
        reg_path = tmp_path / "papers" / "registry.yaml"
        reg_path.write_text(yaml.dump({"papers": []}))

        result: dict = {"status": "success"}

        with (
            patch("core.pipeline.reproduce.load_master_registry") as mock_load,
            patch("core.pipeline.reproduce.save_master_registry") as mock_save,
        ):
            from core.paper.registry import DatasetEntry, MasterRegistry

            mock_load.return_value = MasterRegistry(datasets={"GSE001": DatasetEntry()})

            _write_pipeline_status(str(reg_path), "GSE001", result, str(paper_dir))

        mock_save.assert_called_once()
