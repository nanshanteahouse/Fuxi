"""Tests for core/paper_registry.py — PaperRegistry YAML module."""
from __future__ import annotations

import os
import sys
from unittest.mock import patch
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
    _dataset_to_dict,
    _dict_to_dataset,
    _reset_pipeline_status,
    main,
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

    def test_dataset_entry_empty_experiments(self) -> None:
        ds = DatasetEntry(gse_id='GSE123456', experiments=[])
        assert ds.experiments == []

    def test_yaml_round_trip_with_experiments(self) -> None:
        """DatasetEntry with experiments survives _dataset_to_dict -> _dict_to_dataset."""
        eg1 = ExperimentGroup(
            group_name='Retina_Amacrine',
            sample_ids=['GSM001', 'GSM002'],
            subset_suffix='_amacrine',
            modality='rna',
            status=DatasetStatus.CONFIG_EXISTS,
            config_path='projects/rna/GSE123456/config_GSE123456.py',
            figures=['fig1.png', 'fig2.png'],
        )
        eg2 = ExperimentGroup(
            group_name='Retina_RGC',
            sample_ids=['GSM003'],
            subset_suffix='_rgc',
            modality='rna',
            status=DatasetStatus.NOT_CONFIGURED,
        )
        ds = DatasetEntry(
            gse_id='GSE123456',
            experiments=[eg1, eg2],
        )
        d = _dataset_to_dict(ds)
        assert 'experiments' in d
        assert len(d['experiments']) == 2
        assert d['experiments'][0]['group_name'] == 'Retina_Amacrine'
        assert d['experiments'][0]['status'] == 'config_exists'
        assert d['experiments'][0]['config_path'] == 'projects/rna/GSE123456/config_GSE123456.py'
        assert d['experiments'][0]['figures'] == ['fig1.png', 'fig2.png']
        assert 'config_path' not in d['experiments'][1]  # skipped because None
        assert 'figures' not in d['experiments'][1]       # skipped because empty
        assert d['experiments'][1]['status'] == 'not_configured'
        ds2 = _dict_to_dataset(d)
        assert ds2.experiments is not None
        assert len(ds2.experiments) == 2
        assert ds2.experiments[0].group_name == 'Retina_Amacrine'
        assert ds2.experiments[0].status is DatasetStatus.CONFIG_EXISTS
        assert ds2.experiments[0].config_path == 'projects/rna/GSE123456/config_GSE123456.py'
        assert ds2.experiments[0].figures == ['fig1.png', 'fig2.png']
        assert ds2.experiments[1].config_path is None
        assert ds2.experiments[1].figures == []
        assert ds2.experiments[1].status is DatasetStatus.NOT_CONFIGURED

    def test_yaml_round_trip_no_experiments(self) -> None:
        """Dict without experiments key produces experiments is None on load."""
        d = {
            'gse_id': 'GSE123456',
            'config_path': '',
            'status': 'not_configured',
            'modality': 'rna',
            'notes': '',
        }
        ds = _dict_to_dataset(d)
        assert ds.experiments is None

    def test_yaml_round_trip_empty_experiments(self) -> None:
        """Dict with empty experiments list produces [] on load."""
        d = {
            'gse_id': 'GSE123456',
            'config_path': '',
            'status': 'not_configured',
            'modality': 'rna',
            'notes': '',
            'experiments': [],
        }
        ds = _dict_to_dataset(d)
        assert ds.experiments == []

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

# ──────────────────────────────────────────────
# build_registry — preserve experiment groups
# ──────────────────────────────────────────────


class TestBuildRegistryPreservesExperiments:
    """Verify that hand-declared experiment groups survive --build."""

    def test_experiments_preserved_when_old_entry_matches(self, tmp_path: Path) -> None:
        """Old registry experiments are copied to matching new entries."""
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"
        registry_path = tmp_path / "registry.yaml"

        # Old registry with a dataset that has experiments
        save_registry({
            "papers": [{
                "pmid": "12345678",
                "paper_dir": "2024_Test_Paper",
                "title": "Test",
                "journal": "J",
                "year": "2024",
                "first_author": "A",
                "doi": "10.1234/test",
                "insights_status": "generated",
                "datasets": [{
                    "gse_id": "GSE999999",
                    "config_path": "",
                    "status": "config_exists",
                    "modality": "rna",
                    "notes": "",
                    "experiments": [{
                        "group_name": "Test_Group",
                        "sample_ids": ["GSM001"],
                        "subset_suffix": "_test",
                        "modality": "rna",
                        "status": "config_exists",
                    }],
                }],
            }],
        }, str(registry_path))

        # Matching paper dir + project dir
        paper_dir = papers_dir / "2024_Test_Paper"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(yaml.dump({
            "paper_meta": {"pmid": "12345678", "title": "Test", "first_author": "A"},
            "data_access": {"geo_ids": ["GSE999999"]},
        }))
        gse_dir = projects_dir / "rna" / "GSE999999"
        gse_dir.mkdir(parents=True)
        (gse_dir / "config_GSE999999.py").write_text('CFG.modality = "rna"\n')

        registry = build_registry(
            papers_dir=str(papers_dir),
            projects_dir=str(projects_dir),
            registry_path=str(registry_path),
        )

        assert len(registry["papers"]) == 1
        ds = registry["papers"][0]["datasets"][0]
        assert "experiments" in ds
        assert len(ds["experiments"]) == 1
        assert ds["experiments"][0]["group_name"] == "Test_Group"
        assert ds["experiments"][0]["status"] == "config_exists"

    def test_new_datasets_without_old_entry_get_no_experiments(self, tmp_path: Path) -> None:
        """Datasets not present in old registry get no experiments field."""
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"
        registry_path = tmp_path / "registry.yaml"

        # Old registry — empty datasets for the paper
        save_registry({
            "papers": [{
                "pmid": "12345678",
                "paper_dir": "2024_Test_Paper",
                "title": "Test",
                "journal": "J",
                "year": "2024",
                "first_author": "A",
                "doi": "10.1234/test",
                "insights_status": "generated",
                "datasets": [],
            }],
        }, str(registry_path))

        # Paper referencing a dataset that WASN'T in old registry
        paper_dir = papers_dir / "2024_Test_Paper"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(yaml.dump({
            "paper_meta": {"pmid": "12345678"},
            "data_access": {"geo_ids": ["GSE999999"]},
        }))
        gse_dir = projects_dir / "rna" / "GSE999999"
        gse_dir.mkdir(parents=True)
        (gse_dir / "config_GSE999999.py").write_text('CFG.modality = "rna"\n')

        registry = build_registry(
            papers_dir=str(papers_dir),
            projects_dir=str(projects_dir),
            registry_path=str(registry_path),
        )

        assert len(registry["papers"]) == 1
        ds = registry["papers"][0]["datasets"][0]
        assert "experiments" not in ds

    def test_no_registry_yaml_does_not_crash(self, tmp_path: Path) -> None:
        """First build without an existing registry.yaml runs cleanly."""
        projects_dir = tmp_path / "projects"
        papers_dir = tmp_path / "papers"
        registry_path = tmp_path / "nonexistent.yaml"

        paper_dir = papers_dir / "Test"
        paper_dir.mkdir(parents=True)
        (paper_dir / "insights.yaml").write_text(yaml.dump({
            "paper_meta": {"pmid": "12345678"},
            "data_access": {"geo_ids": ["GSE999999"]},
        }))
        gse_dir = projects_dir / "rna" / "GSE999999"
        gse_dir.mkdir(parents=True)
        (gse_dir / "config_GSE999999.py").write_text('CFG.modality = "rna"\n')

        registry = build_registry(
            papers_dir=str(papers_dir),
            projects_dir=str(projects_dir),
            registry_path=str(registry_path),
        )

        assert len(registry["papers"]) == 1
        assert "experiments" not in registry["papers"][0]["datasets"][0]


# ──────────────────────────────────────────────
# _reset_pipeline_status
# ──────────────────────────────────────────────


class TestResetPipelineStatus:
    """Verify _reset_pipeline_status reverts pipeline_complete entries.

    Tests cover full reset, filtering by paper/GSE, dry-run mode, and
    edge cases with no pipeline_complete entries.
    """

    @staticmethod
    def _make_registry(tmp_path, papers_data=None):
        """Create a temporary registry.yaml from paper data dicts."""
        if papers_data is None:
            papers_data = [
                {
                    "pmid": "11111111",
                    "paper_dir": "Paper_A",
                    "title": "Paper A",
                    "journal": "J",
                    "year": "2024",
                    "first_author": "A",
                    "doi": "10.1/test",
                    "insights_status": "generated",
                    "datasets": [
                        {
                            "gse_id": "GSE001",
                            "config_path": "/cfg/GSE001.py",
                            "status": "pipeline_complete",
                            "modality": "rna",
                            "notes": "",
                        }
                    ],
                }
            ]
        registry = {"papers": papers_data}
        reg_path = tmp_path / "registry.yaml"
        save_registry(registry, str(reg_path))
        return str(reg_path)

    # ── tests ─────────────────────────────────

    def test_full_reset(self, tmp_path):
        """Multiple pipeline_complete entries are all reverted.

        Two papers with three datasets total, one with config_path and one
        without — verifies both CONFIG_EXISTS and NOT_CONFIGURED outcomes.
        """
        papers_data = [
            {
                "pmid": "11111111",
                "paper_dir": "Paper_A",
                "title": "A", "journal": "J", "year": "2024",
                "first_author": "A", "doi": "10.1/test",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE001",
                        "config_path": "/cfg/GSE001.py",
                        "status": "pipeline_complete",
                        "modality": "rna",
                        "notes": "",
                    },
                    {
                        "gse_id": "GSE002",
                        "config_path": "",
                        "status": "pipeline_complete",
                        "modality": "atac",
                        "notes": "",
                    },
                ],
            },
            {
                "pmid": "22222222",
                "paper_dir": "Paper_B",
                "title": "B", "journal": "J", "year": "2024",
                "first_author": "B", "doi": "10.2/test",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE003",
                        "config_path": "/cfg/GSE003.py",
                        "status": "pipeline_complete",
                        "modality": "spatial",
                        "notes": "",
                    }
                ],
            },
        ]
        reg_path = self._make_registry(tmp_path, papers_data)

        with patch("core.paper_registry._reset_single_dataset_yaml"):
            result = _reset_pipeline_status(reg_path)

        assert result["papers_affected"] == 2
        assert result["datasets_reset"] == 3
        assert result["dataset_yamls_updated"] == []

        registry = load_registry(reg_path)
        ds_a0 = registry["papers"][0]["datasets"][0]
        ds_a1 = registry["papers"][0]["datasets"][1]
        ds_b = registry["papers"][1]["datasets"][0]
        assert ds_a0["status"] == "config_exists"   # has config_path
        assert ds_a1["status"] == "not_configured"   # no config_path
        assert ds_b["status"] == "config_exists"    # has config_path

    def test_filter_paper(self, tmp_path):
        """Only the paper matching filter_paper is affected."""
        papers_data = [
            {
                "pmid": "11111111",
                "paper_dir": "KeepPaper",
                "title": "Keep", "journal": "J", "year": "2024",
                "first_author": "K", "doi": "10.1/keep",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE001",
                        "config_path": "/cfg/GSE001.py",
                        "status": "pipeline_complete",
                        "modality": "rna",
                        "notes": "",
                    }
                ],
            },
            {
                "pmid": "22222222",
                "paper_dir": "OtherPaper",
                "title": "Other", "journal": "J", "year": "2024",
                "first_author": "O", "doi": "10.2/other",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE002",
                        "config_path": "/cfg/GSE002.py",
                        "status": "pipeline_complete",
                        "modality": "rna",
                        "notes": "",
                    }
                ],
            },
        ]
        reg_path = self._make_registry(tmp_path, papers_data)

        # Filter for a paper that does NOT exist
        with patch("core.paper_registry._reset_single_dataset_yaml"):
            result = _reset_pipeline_status(reg_path, filter_paper="NonExistent")
        assert result["papers_affected"] == 0
        assert result["datasets_reset"] == 0

        # Filter for KeepPaper — should reset that entry
        result = _reset_pipeline_status(reg_path, filter_paper="KeepPaper")
        assert result["papers_affected"] == 1
        assert result["datasets_reset"] == 1
        registry = load_registry(reg_path)
        assert registry["papers"][0]["datasets"][0]["status"] == "config_exists"
        # OtherPaper unchanged
        assert registry["papers"][1]["datasets"][0]["status"] == "pipeline_complete"

    def test_filter_gse(self, tmp_path):
        """Only the dataset matching filter_gse is reset within a paper."""
        papers_data = [
            {
                "pmid": "11111111",
                "paper_dir": "Paper_A",
                "title": "A", "journal": "J", "year": "2024",
                "first_author": "A", "doi": "10.1/test",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE001",
                        "config_path": "/cfg/GSE001.py",
                        "status": "pipeline_complete",
                        "modality": "rna",
                        "notes": "",
                    },
                    {
                        "gse_id": "GSE002",
                        "config_path": "/cfg/GSE002.py",
                        "status": "pipeline_complete",
                        "modality": "atac",
                        "notes": "",
                    },
                ],
            },
        ]
        reg_path = self._make_registry(tmp_path, papers_data)

        with patch("core.paper_registry._reset_single_dataset_yaml"):
            result = _reset_pipeline_status(reg_path, filter_gse="GSE001")

        assert result["papers_affected"] == 1
        assert result["datasets_reset"] == 1

        registry = load_registry(reg_path)
        assert registry["papers"][0]["datasets"][0]["status"] == "config_exists"  # reset
        assert registry["papers"][0]["datasets"][1]["status"] == "pipeline_complete"  # intact

    def test_dry_run(self, tmp_path):
        """Dry-run mode inspects without modifying any files."""
        papers_data = [
            {
                "pmid": "11111111",
                "paper_dir": "Paper_A",
                "title": "A", "journal": "J", "year": "2024",
                "first_author": "A", "doi": "10.1/test",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE001",
                        "config_path": "/cfg/GSE001.py",
                        "status": "pipeline_complete",
                        "modality": "rna",
                        "notes": "",
                    }
                ],
            },
        ]
        reg_path = self._make_registry(tmp_path, papers_data)
        import pathlib
        original = pathlib.Path(reg_path).read_text()

        with patch("core.paper_registry._reset_single_dataset_yaml") as mock_reset:
            result = _reset_pipeline_status(reg_path, dry_run=True)

        assert result["papers_affected"] == 1
        assert result["datasets_reset"] == 1
        assert result["dataset_yamls_updated"] == []

        # Registry file unchanged
        assert pathlib.Path(reg_path).read_text() == original
        # _reset_single_dataset_yaml still called (dry_run applies inside it)
        mock_reset.assert_called_once()

    def test_no_pipeline_complete(self, tmp_path):
        """Registry with no pipeline_complete entries returns empty summary."""
        papers_data = [
            {
                "pmid": "11111111",
                "paper_dir": "Paper_A",
                "title": "A", "journal": "J", "year": "2024",
                "first_author": "A", "doi": "10.1/test",
                "insights_status": "generated",
                "datasets": [
                    {
                        "gse_id": "GSE001",
                        "config_path": "/cfg/GSE001.py",
                        "status": "config_exists",
                        "modality": "rna",
                        "notes": "",
                    }
                ],
            },
        ]
        reg_path = self._make_registry(tmp_path, papers_data)

        with patch("core.paper_registry._reset_single_dataset_yaml") as mock_reset:
            result = _reset_pipeline_status(reg_path)

        assert result["papers_affected"] == 0
        assert result["datasets_reset"] == 0
        assert result["dataset_yamls_updated"] == []
        mock_reset.assert_not_called()


# ──────────────────────────────────────────────
# CLI — --reset flag
# ──────────────────────────────────────────────


class TestCLIReset:
    """Verify --reset CLI flag routes correctly to _reset_pipeline_status."""

    def test_reset_default(self, tmp_path):
        """--reset calls _reset_pipeline_status with registry_path default."""
        reg_path = tmp_path / "registry.yaml"
        reg_path.write_text("papers: []\n")

        test_argv = ["paper_registry.py", "--reset", "--output", str(reg_path)]
        with patch.object(sys, "argv", test_argv):
            with patch("core.paper_registry._reset_pipeline_status") as mock_fn:
                mock_fn.return_value = {
                    "papers_affected": 0, "datasets_reset": 0,
                    "dataset_yamls_updated": [],
                }
                main()

        mock_fn.assert_called_once_with(
            registry_path=str(reg_path),
            filter_paper=None,
            filter_gse=None,
            dry_run=False,
        )

    def test_reset_with_paper(self, tmp_path):
        """--reset --paper <dir> passes filter_paper."""
        reg_path = tmp_path / "registry.yaml"
        reg_path.write_text("papers: []\n")

        test_argv = [
            "paper_registry.py", "--reset", "--paper", "2024_Test_Paper",
            "--output", str(reg_path),
        ]
        with patch.object(sys, "argv", test_argv):
            with patch("core.paper_registry._reset_pipeline_status") as mock_fn:
                mock_fn.return_value = {
                    "papers_affected": 1, "datasets_reset": 0,
                    "dataset_yamls_updated": [],
                }
                main()

        mock_fn.assert_called_once_with(
            registry_path=str(reg_path),
            filter_paper="2024_Test_Paper",
            filter_gse=None,
            dry_run=False,
        )

    def test_reset_with_gse(self, tmp_path):
        """--reset --gse <GSE_ID> passes filter_gse."""
        reg_path = tmp_path / "registry.yaml"
        reg_path.write_text("papers: []\n")

        test_argv = [
            "paper_registry.py", "--reset", "--gse", "GSE107618",
            "--output", str(reg_path),
        ]
        with patch.object(sys, "argv", test_argv):
            with patch("core.paper_registry._reset_pipeline_status") as mock_fn:
                mock_fn.return_value = {
                    "papers_affected": 1, "datasets_reset": 1,
                    "dataset_yamls_updated": [],
                }
                main()

        mock_fn.assert_called_once_with(
            registry_path=str(reg_path),
            filter_paper=None,
            filter_gse="GSE107618",
            dry_run=False,
        )

    def test_reset_dry_run(self, tmp_path):
        """--reset --dry-run passes dry_run=True."""
        reg_path = tmp_path / "registry.yaml"
        reg_path.write_text("papers: []\n")

        test_argv = [
            "paper_registry.py", "--reset", "--dry-run",
            "--output", str(reg_path),
        ]
        with patch.object(sys, "argv", test_argv):
            with patch("core.paper_registry._reset_pipeline_status") as mock_fn:
                mock_fn.return_value = {
                    "papers_affected": 0, "datasets_reset": 0,
                    "dataset_yamls_updated": [],
                }
                main()

        mock_fn.assert_called_once_with(
            registry_path=str(reg_path),
            filter_paper=None,
            filter_gse=None,
            dry_run=True,
        )

    def test_paper_and_gse_mutually_exclusive(self, tmp_path):
        """--paper and --gse together raise SystemExit."""
        import pytest
        test_argv = [
            "paper_registry.py", "--reset",
            "--paper", "TestPaper",
            "--gse", "GSE001",
            "--output", str(tmp_path / "registry.yaml"),
        ]
        with patch.object(sys, "argv", test_argv):
            with pytest.raises(SystemExit):
                main()
