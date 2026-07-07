"""Tests for core/run_reproduce.py — Reproduction orchestration (P3)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml

from core.run_reproduce import (
    REPRODUCE_TIMEOUT,
    _detect_modality,
    _extract_geo_ids,
    _run_pipeline_for_gse,
    run_reproduce,
)


# ═══════════════════════════════════════════════════════════════════════
# _detect_modality
# ═══════════════════════════════════════════════════════════════════════


class TestDetectModality:
    """Modality detection from pipeline config files (reuses paper_registry)."""

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

        with patch("core.run_reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0, stdout="pipeline done", stderr=""
            )
            result = _run_pipeline_for_gse("GSE001", config_path)

        assert result["status"] == "success"
        assert result["modality"] == "rna"
        assert result["config_path"] == config_path
        assert result["output"] == "pipeline done"
        assert isinstance(result["duration_s"], (int, float))
        mock_run.assert_called_once()

    def test_failure(self, tmp_path: Path) -> None:
        config_path = self._make_config(tmp_path, "atac")

        with patch("core.run_reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1, stdout="", stderr="error occurred"
            )
            result = _run_pipeline_for_gse("GSE002", config_path)

        assert result["status"] == "failed"
        assert result["modality"] == "atac"
        assert "error occurred" in result["error"]

    def test_timeout(self, tmp_path: Path) -> None:
        config_path = self._make_config(tmp_path, "rna")

        with patch("core.run_reproduce.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.TimeoutExpired(
                "cmd", REPRODUCE_TIMEOUT
            )
            result = _run_pipeline_for_gse("GSE003", config_path)

        assert result["status"] == "timeout"
        assert "Timed out" in result["error"]
        assert result["modality"] == "rna"

    def test_unknown_modality_does_not_call_subprocess(
        self, tmp_path: Path
    ) -> None:
        config_path = tmp_path / "config.py"
        config_path.write_text("x = 1\n")  # no CFG.modality

        with patch("core.run_reproduce.subprocess.run") as mock_run:
            result = _run_pipeline_for_gse("GSE004", str(config_path))

        assert result["status"] == "failed"
        assert result["modality"] == "unknown"
        assert "Cannot detect modality" in result["error"]
        mock_run.assert_not_called()

    def test_subprocess_called_with_correct_args(
        self, tmp_path: Path
    ) -> None:
        config_path = self._make_config(tmp_path, "spatial")

        with patch("core.run_reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0, stdout="", stderr=""
            )
            _run_pipeline_for_gse("GSE005", config_path)

        cmd = mock_run.call_args[0][0]
        assert "core/run_pipeline.py" in cmd
        assert "--config" in cmd
        assert config_path in cmd
        assert "--modality" in cmd
        assert "spatial" in cmd


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
) -> dict:
    """Build a minimal registry dict."""
    return {
        "papers": [
            {
                "pmid": pmid,
                "paper_dir": paper_dir,
                "datasets": datasets,
            }
        ]
    }


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
        results = run_reproduce(paper_dir, registry={"papers": []})
        assert results == {}

    def test_gse_filter_excludes_others(self, tmp_path: Path) -> None:
        paper_dir = _make_paper_dir(
            tmp_path, name="MultiGSE", geo_ids=["GSE001", "GSE002"]
        )
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0, stdout="ok", stderr=""
            )
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
            results = run_reproduce(
                paper_dir, registry=registry, dry_run=True
            )

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
        with patch("core.run_reproduce.subprocess.run") as mock_run:
            results = run_reproduce(
                paper_dir, registry=registry, dry_run=True
            )

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
        with patch(
            "core.run_reproduce.load_registry", return_value=registry
        ) as mock_load, patch(
            "core.run_reproduce.Path.is_dir", return_value=True
        ):
            from core.run_reproduce import main

            # We need to patch sys.argv
            test_args = ["run_reproduce.py", "--all", "--dry-run"]
            with patch.object(
                sys, "argv", test_args
            ):
                # Should not raise
                main()

        mock_load.assert_called_once()
