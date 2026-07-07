"""Tests for core/paper_registry.py — PaperRegistry YAML module."""
from __future__ import annotations

import os
import yaml
from pathlib import Path



from core.paper_registry import (
    DatasetStatus,
    DatasetEntry,
    PaperEntry,
    load_registry,
    save_registry,
    build_registry,
    detect_modality,
    _scan_insights_yamls,
    _scan_project_dirs,
    _find_data_only_entries,
    _paper_to_dict,
    _dict_to_paper,
)

from core.paper_registry_models import ExperimentGroup
from dataclasses import asdict
# ──────────────────────────────────────────────
# DatasetStatus enum
# ──────────────────────────────────────────────


class TestDatasetStatus:
    """Verify DatasetStatus enum values and members."""

    def test_values(self) -> None:
        assert DatasetStatus.NOT_CONFIGURED.value == "not_configured"
        assert DatasetStatus.CONFIG_EXISTS.value == "config_exists"
        assert DatasetStatus.PIPELINE_COMPLETE.value == "pipeline_complete"
        assert DatasetStatus.DATA_NOT_DOWNLOADED.value == "data_not_downloaded"
        assert DatasetStatus.DATA_ONLY.value == "data_only"
        assert DatasetStatus.UNKNOWN.value == "unknown"

    def test_member_count(self) -> None:
        assert len(DatasetStatus) == 6

    def test_is_str_enum(self) -> None:
        """DatasetStatus must be a valid str Enum for YAML compat."""
        assert isinstance(DatasetStatus.CONFIG_EXISTS, str)


# ──────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────


class TestDataModel:
    """PaperEntry / DatasetEntry construction and round-trip."""

    def test_dataset_entry_defaults(self) -> None:
        ds = DatasetEntry(gse_id="GSE123456")
        assert ds.gse_id == "GSE123456"
        assert ds.config_path == ""
        assert ds.status is DatasetStatus.NOT_CONFIGURED
        assert ds.modality == "rna"
        assert ds.notes == ""

    def test_paper_entry_defaults(self) -> None:
        paper = PaperEntry(pmid="12345678", paper_dir="test_dir")
        assert paper.pmid == "12345678"
        assert paper.paper_dir == "test_dir"
        assert paper.title == ""
        assert paper.journal == ""
        assert paper.year == ""
        assert paper.first_author == ""
        assert paper.doi == ""
        assert paper.datasets == []
        assert paper.insights_status == "generated"

    def test_paper_round_trip(self) -> None:
        """Serialize PaperEntry -> dict -> PaperEntry preserves data."""
        ds = DatasetEntry(
            gse_id="GSE123456",
            config_path="projects/rna/GSE123456/config_GSE123456.py",
            status=DatasetStatus.CONFIG_EXISTS,
            modality="rna",
            notes="linked via Menon 2019",
        )
        paper = PaperEntry(
            pmid="31653841",
            paper_dir="2019_Menon_Nature_Com_...",
            title="Single-cell transcriptomic atlas of the human retina",
            journal="Nature Communications",
            year="2019",
            first_author="Menon",
            doi="10.1038/s41467-019-12780-8",
            datasets=[ds],
            insights_status="generated",
        )
        d = _paper_to_dict(paper)
        paper2 = _dict_to_paper(d)

        assert paper2.pmid == paper.pmid
        assert paper2.paper_dir == paper.paper_dir
        assert paper2.title == paper.title
        assert paper2.journal == paper.journal
        assert paper2.year == paper.year
        assert paper2.first_author == paper.first_author
        assert paper2.doi == paper.doi
        assert paper2.insights_status == paper.insights_status
        assert len(paper2.datasets) == 1
        assert paper2.datasets[0].gse_id == "GSE123456"
        assert paper2.datasets[0].status is DatasetStatus.CONFIG_EXISTS
        assert paper2.datasets[0].modality == "rna"
        assert paper2.datasets[0].notes == "linked via Menon 2019"

    def test_paper_no_datasets_round_trip(self) -> None:
        paper = PaperEntry(pmid="12345678", paper_dir="empty", insights_status="no_geo")
        d = _paper_to_dict(paper)
        paper2 = _dict_to_paper(d)
        assert paper2.pmid == "12345678"
        assert paper2.datasets == []



class TestExperimentGroup:
    """Verify ExperimentGroup dataclass construction and round-trip."""

    def test_create_all_fields(self) -> None:
        eg = ExperimentGroup(
            group_name="Retina_Amacrine",
            sample_ids=["GSM001", "GSM002"],
            subset_suffix="_amacrine",
            modality="rna",
            status=DatasetStatus.CONFIG_EXISTS,
            config_path="projects/rna/GSE123456/config_GSE123456.py",
            figures=["fig1.png", "fig2.png"],
        )
        assert eg.group_name == "Retina_Amacrine"
        assert eg.sample_ids == ["GSM001", "GSM002"]
        assert eg.subset_suffix == "_amacrine"
        assert eg.modality == "rna"
        assert eg.status is DatasetStatus.CONFIG_EXISTS
        assert eg.config_path == "projects/rna/GSE123456/config_GSE123456.py"
        assert eg.figures == ["fig1.png", "fig2.png"]

    def test_default_config_path(self) -> None:
        eg = ExperimentGroup(
            group_name="Test",
            sample_ids=["GSM001"],
            subset_suffix="_test",
            modality="rna",
            status=DatasetStatus.NOT_CONFIGURED,
        )
        assert eg.config_path is None

    def test_default_figures(self) -> None:
        eg = ExperimentGroup(
            group_name="Test",
            sample_ids=["GSM001"],
            subset_suffix="_test",
            modality="rna",
            status=DatasetStatus.NOT_CONFIGURED,
        )
        assert eg.figures == []

    def test_round_trip_dict(self) -> None:
        eg = ExperimentGroup(
            group_name="Retina_Amacrine",
            sample_ids=["GSM001", "GSM002"],
            subset_suffix="_amacrine",
            modality="rna",
            status=DatasetStatus.CONFIG_EXISTS,
            config_path="projects/rna/GSE123456/config_GSE123456.py",
            figures=["fig1.png", "fig2.png"],
        )
        d = asdict(eg)
        eg2 = ExperimentGroup(**d)
        assert eg2 == eg

    def test_dataset_entry_with_experiments(self) -> None:
        eg = ExperimentGroup(
            group_name="Retina_Amacrine",
            sample_ids=["GSM001"],
            subset_suffix="_amacrine",
            modality="rna",
            status=DatasetStatus.CONFIG_EXISTS,
        )
        ds = DatasetEntry(
            gse_id="GSE123456",
            experiments=[eg],
        )
        assert ds.gse_id == "GSE123456"
        assert ds.experiments is not None
        assert len(ds.experiments) == 1
        assert ds.experiments[0].group_name == "Retina_Amacrine"
        assert ds.experiments[0].status is DatasetStatus.CONFIG_EXISTS

    def test_dataset_entry_no_experiments(self) -> None:
        ds = DatasetEntry(gse_id="GSE123456")
        assert ds.experiments is None

# ──────────────────────────────────────────────
# detect_modality
# ──────────────────────────────────────────────


class TestDetectModality:
    """Regex-based modality detection from config files."""

    def test_detect_rna(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text('CFG.modality = "rna"\n')
        assert detect_modality(str(p)) == "rna"

    def test_detect_atac(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text('CFG.modality = "atac"\n')
        assert detect_modality(str(p)) == "atac"

    def test_detect_spatial(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text('CFG.modality = "spatial"\n')
        assert detect_modality(str(p)) == "spatial"

    def test_detect_missing_file(self) -> None:
        assert detect_modality("/nonexistent/path.py") == "unknown"

    def test_detect_no_modality_line(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text("x = 1\n")
        assert detect_modality(str(p)) == "unknown"

    def test_detect_single_quotes(self, tmp_path: Path) -> None:
        p = tmp_path / "config_test.py"
        p.write_text("CFG.modality = 'rna'\n")
        assert detect_modality(str(p)) == "rna"

    def test_detect_on_real_config(self) -> None:
        """Use a real config file in the repo to verify detection."""
        real = "projects/rna/GSE107618/config_GSE107618.py"
        if os.path.exists(real):
            assert detect_modality(real) == "rna"


# ──────────────────────────────────────────────
# _scan_insights_yamls
# ──────────────────────────────────────────────


class TestScanInsightsYamls:
    """Scanning paper directories for insights.yaml metadata."""

    def test_basic(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "2024_Test_Paper"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {
                    "pmid": "12345678",
                    "title": "Test Paper",
                    "first_author": "Test",
                    "journal": "Test Journal",
                    "year": "2024",
                    "doi": "10.1234/test",
                },
                "data_access": {"geo_ids": ["GSE123456"]},
            })
        )
        result = _scan_insights_yamls(str(papers_dir))
        assert len(result) == 1
        r = result[0]
        assert r["pmid"] == "12345678"
        assert r["geo_ids"] == ["GSE123456"]
        assert r["paper_dir"] == "2024_Test_Paper"
        assert r["title"] == "Test Paper"

    def test_empty_geo_ids(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "Test"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "87654321", "title": "No GEO"},
                "data_access": {"geo_ids": []},
            })
        )
        result = _scan_insights_yamls(str(papers_dir))
        assert len(result) == 1
        assert result[0]["geo_ids"] == []

    def test_no_pmid_skipped(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "NoPmid"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"title": "No PMID"},
                "data_access": {"geo_ids": []},
            })
        )
        result = _scan_insights_yamls(str(papers_dir))
        assert len(result) == 0

    def test_corrupt_yaml_skipped(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "Corrupt"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text("{corrupt: yaml: [unclosed\n")
        result = _scan_insights_yamls(str(papers_dir))
        assert len(result) == 0

    def test_dedup_geo_ids(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "Dup"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "11111111"},
                "data_access": {"geo_ids": ["GSE001", "GSE001", "GSE002", "GSE001"]},
            })
        )
        result = _scan_insights_yamls(str(papers_dir))
        assert result[0]["geo_ids"] == ["GSE001", "GSE002"]

    def test_no_insights_file(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "NoInsights"
        paper_dir.mkdir(parents=True)
        # No insights.yaml in the directory
        result = _scan_insights_yamls(str(papers_dir))
        assert len(result) == 0

    def test_non_dir_skipped(self, tmp_path: Path) -> None:
        papers_dir = tmp_path / "papers"
        papers_dir.mkdir(parents=True)
        (papers_dir / "not_a_dir.txt").write_text("hello")
        result = _scan_insights_yamls(str(papers_dir))
        assert len(result) == 0


# ──────────────────────────────────────────────
# _scan_project_dirs
# ──────────────────────────────────────────────


class TestScanProjectDirs:
    """Scanning projects/{rna,atac,spatial} for GSE directories."""

    def test_basic(self, tmp_path: Path) -> None:
        projects_dir = tmp_path / "projects"
        rna_dir = projects_dir / "rna"
        gse_dir = rna_dir / "GSE123456"
        gse_dir.mkdir(parents=True)
        (gse_dir / "config_GSE123456.py").write_text("")
        result = _scan_project_dirs(str(projects_dir))
        assert len(result) == 1
        assert result[0]["gse_id"] == "GSE123456"
        assert result[0]["modality"] == "rna"
        assert result[0]["config_path"] != ""

    def test_empty_atac(self, tmp_path: Path) -> None:
        """Empty projects/atac/ should not cause errors."""
        projects_dir = tmp_path / "projects"
        (projects_dir / "rna" / "GSE001").mkdir(parents=True)
        (projects_dir / "rna" / "GSE001" / "config_GSE001.py").write_text("")
        (projects_dir / "atac").mkdir()
        result = _scan_project_dirs(str(projects_dir))
        assert len(result) == 1

    def test_gse_without_config(self, tmp_path: Path) -> None:
        projects_dir = tmp_path / "projects"
        gse_dir = projects_dir / "rna" / "GSE999999"
        gse_dir.mkdir(parents=True)
        result = _scan_project_dirs(str(projects_dir))
        assert len(result) == 1
        assert result[0]["config_path"] == ""

    def test_multiple_modalities(self, tmp_path: Path) -> None:
        projects_dir = tmp_path / "projects"
        (projects_dir / "rna" / "GSE001").mkdir(parents=True)
        (projects_dir / "rna" / "GSE001" / "config_GSE001.py").write_text("")
        (projects_dir / "atac" / "GSE002").mkdir(parents=True)
        (projects_dir / "atac" / "GSE002" / "config_GSE002.py").write_text("")
        (projects_dir / "spatial" / "GSE003").mkdir(parents=True)
        (projects_dir / "spatial" / "GSE003" / "config_GSE003.py").write_text("")
        result = _scan_project_dirs(str(projects_dir))
        assert len(result) == 3
        modalities = {e["gse_id"]: e["modality"] for e in result}
        assert modalities["GSE001"] == "rna"
        assert modalities["GSE002"] == "atac"
        assert modalities["GSE003"] == "spatial"


# ──────────────────────────────────────────────
# _find_data_only_entries
# ──────────────────────────────────────────────


class TestFindDataOnly:
    """Identifying GSE configs not linked to any paper."""

    def test_data_only_detected(self) -> None:
        gses = [
            {"gse_id": "GSE001", "config_path": "/some/path/config.py", "modality": "rna"},
            {"gse_id": "GSE002", "config_path": "/other/path/config.py", "modality": "rna"},
        ]
        all_paper_ids = {"GSE001"}  # GSE001 is in a paper, GSE002 is not
        result = _find_data_only_entries(gses, all_paper_ids)
        assert len(result) == 1
        assert result[0].gse_id == "GSE002"
        assert result[0].status is DatasetStatus.DATA_ONLY

    def test_none_when_all_linked(self) -> None:
        gses = [
            {"gse_id": "GSE001", "config_path": "config.py", "modality": "rna"},
        ]
        result = _find_data_only_entries(gses, {"GSE001"})
        assert len(result) == 0

    def test_only_config_gses(self) -> None:
        """Only GSEs with config_path set get data_only status."""
        gses = [
            {"gse_id": "GSE001", "config_path": "config.py", "modality": "rna"},
            {"gse_id": "GSE002", "config_path": "", "modality": "rna"},
        ]
        result = _find_data_only_entries(gses, set())
        assert len(result) == 1
        assert result[0].gse_id == "GSE001"


# ──────────────────────────────────────────────
# build_registry (integration)
# ──────────────────────────────────────────────


class TestBuildRegistry:
    """Full build_registry() with temp directory fixtures."""

    def test_basic_build(self, tmp_path: Path) -> None:
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"

        # Paper with geo_ids
        paper_dir = papers_dir / "2024_Test_Paper"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {
                    "pmid": "12345678",
                    "title": "Test Paper",
                    "first_author": "Test",
                    "journal": "Test J",
                    "year": "2024",
                    "doi": "10.1234/test",
                },
                "data_access": {"geo_ids": ["GSE999999"]},
            })
        )

        # Matching GSE directory with config
        gse_dir = projects_dir / "rna" / "GSE999999"
        gse_dir.mkdir(parents=True)
        (gse_dir / "config_GSE999999.py").write_text('CFG.modality = "rna"\n')

        registry = build_registry(str(papers_dir), str(projects_dir))
        assert len(registry["papers"]) == 1
        paper = registry["papers"][0]
        assert paper["pmid"] == "12345678"
        assert len(paper["datasets"]) == 1
        assert paper["datasets"][0]["gse_id"] == "GSE999999"
        assert paper["datasets"][0]["status"] == "config_exists"
        assert paper["datasets"][0]["modality"] == "rna"
        assert "data_only_datasets" not in registry

    def test_data_not_downloaded(self, tmp_path: Path) -> None:
        """Paper references GSE not found in projects/."""
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "Test"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "12345678"},
                "data_access": {"geo_ids": ["GSE999999"]},
            })
        )
        registry = build_registry(str(papers_dir), str(tmp_path / "projects"))
        paper = registry["papers"][0]
        assert len(paper["datasets"]) == 1
        assert paper["datasets"][0]["status"] == "data_not_downloaded"

    def test_data_only_entries(self, tmp_path: Path) -> None:
        """GSE config exists but no paper links to it."""
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"
        papers_dir.mkdir(parents=True)

        gse_dir = projects_dir / "rna" / "GSE999999"
        gse_dir.mkdir(parents=True)
        (gse_dir / "config_GSE999999.py").write_text('CFG.modality = "rna"\n')

        registry = build_registry(str(papers_dir), str(projects_dir))
        assert registry["papers"] == []
        assert "data_only_datasets" in registry
        assert len(registry["data_only_datasets"]) == 1
        assert registry["data_only_datasets"][0]["gse_id"] == "GSE999999"
        assert registry["data_only_datasets"][0]["status"] == "data_only"

    def test_no_geo_ids_paper(self, tmp_path: Path) -> None:
        """Paper with empty geo_ids gets insights_status = 'no_geo'."""
        papers_dir = tmp_path / "papers"
        paper_dir = papers_dir / "NoGEO"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "12345678", "title": "No GEO Paper"},
                "data_access": {"geo_ids": []},
            })
        )
        registry = build_registry(str(papers_dir), str(tmp_path / "projects"))
        assert len(registry["papers"]) == 1
        assert registry["papers"][0]["insights_status"] == "no_geo"
        assert registry["papers"][0].get("datasets", []) == []

    def test_empty_registry(self, tmp_path: Path) -> None:
        """No papers and no projects yield empty registry."""
        registry = build_registry(
            str(tmp_path / "empty_papers"), str(tmp_path / "empty_projects")
        )
        assert registry["papers"] == []

    def test_multiple_gses_per_paper(self, tmp_path: Path) -> None:
        """Paper with multiple geo_ids links to multiple GSEs."""
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"

        paper_dir = papers_dir / "MultiGSE"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "11111111"},
                "data_access": {"geo_ids": ["GSE001", "GSE002"]},
            })
        )

        (projects_dir / "rna" / "GSE001").mkdir(parents=True)
        (projects_dir / "rna" / "GSE001" / "config_GSE001.py").write_text("")
        (projects_dir / "rna" / "GSE002").mkdir(parents=True)
        (projects_dir / "rna" / "GSE002" / "config_GSE002.py").write_text("")

        registry = build_registry(str(papers_dir), str(projects_dir))
        assert len(registry["papers"]) == 1
        assert len(registry["papers"][0]["datasets"]) == 2

    def test_multiple_modalities_same_gse(self, tmp_path: Path) -> None:
        """Same GSE ID in rna and atac creates two dataset entries."""
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"

        paper_dir = papers_dir / "MultiMod"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "22222222"},
                "data_access": {"geo_ids": ["GSE001"]},
            })
        )

        # Same GSE in both rna and atac
        (projects_dir / "rna" / "GSE001").mkdir(parents=True)
        (projects_dir / "rna" / "GSE001" / "config_GSE001.py").write_text("")
        (projects_dir / "atac" / "GSE001").mkdir(parents=True)
        (projects_dir / "atac" / "GSE001" / "config_GSE001.py").write_text("")

        registry = build_registry(str(papers_dir), str(projects_dir))
        assert len(registry["papers"]) == 1
        assert len(registry["papers"][0]["datasets"]) == 2

    def test_corrupt_paper_skipped(self, tmp_path: Path) -> None:
        """Corrupt insights.yaml skips that paper without crashing."""
        papers_dir = tmp_path / "papers"
        projects_dir = tmp_path / "projects"

        # Good paper
        good = papers_dir / "Good"
        good.mkdir(parents=True)
        (good / "insights.yaml").write_text(
            yaml.dump({
                "paper_meta": {"pmid": "11111111"},
                "data_access": {"geo_ids": ["GSE001"]},
            })
        )
        (projects_dir / "rna" / "GSE001").mkdir(parents=True)
        (projects_dir / "rna" / "GSE001" / "config_GSE001.py").write_text("")

        # Corrupt paper
        bad = papers_dir / "Bad"
        bad.mkdir()
        (bad / "insights.yaml").write_text("{bad yaml: [}\n")

        registry = build_registry(str(papers_dir), str(projects_dir))
        assert len(registry["papers"]) == 1
        assert registry["papers"][0]["pmid"] == "11111111"


# ──────────────────────────────────────────────
# File I/O
# ──────────────────────────────────────────────


class TestFileIO:
    """save_registry / load_registry round-trips."""

    def test_save_load_round_trip(self, tmp_path: Path) -> None:
        reg_path = tmp_path / "registry.yaml"
        registry = {
            "papers": [
                {
                    "pmid": "31653841",
                    "paper_dir": "2019_Menon_Nature_Com_...",
                    "title": "Test",
                    "journal": "Nat Com",
                    "year": "2019",
                    "first_author": "Menon",
                    "doi": "10.1038/s41467-019-12780-8",
                    "insights_status": "generated",
                    "datasets": [
                        {
                            "gse_id": "GSE107618",
                            "config_path": "projects/rna/GSE107618/config_GSE107618.py",
                            "status": "config_exists",
                            "modality": "rna",
                            "notes": "",
                        }
                    ],
                }
            ]
        }
        save_registry(registry, str(reg_path))
        assert reg_path.exists()
        loaded = load_registry(str(reg_path))
        assert loaded == registry

    def test_load_missing_file(self, tmp_path: Path) -> None:
        result = load_registry(str(tmp_path / "nonexistent.yaml"))
        assert result == {"papers": []}

    def test_load_empty_file(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.yaml"
        path.write_text("")
        result = load_registry(str(path))
        assert result == {"papers": []}

    def test_yaml_format(self, tmp_path: Path) -> None:
        """Saved YAML should use block style (default_flow_style=False)."""
        reg_path = tmp_path / "registry.yaml"
        registry = {
            "papers": [
                {
                    "pmid": "12345678",
                    "paper_dir": "test",
                    "title": "T",
                    "journal": "J",
                    "year": "2024",
                    "first_author": "A",
                    "doi": "10.1234/test",
                    "insights_status": "generated",
                }
            ]
        }
        save_registry(registry, str(reg_path))
        content = reg_path.read_text()
        assert "pmid:" in content  # block-style keys
        assert "paper_dir:" in content
