#!/usr/bin/env python3
"""
core/registry.py — 五域统一论文登记表 (Master Registry)

数据模型 & YAML I/O & 查询 API
================================

统一管理五个区域的论文引用：
  1. projects/papers/           — 论文 XML / insights 解读
  2. projects/{rna,atac,spatial}/ — 管线运行产物
  3. notes/supplements/         — 论文附表
  4. $FUXI_DATA_ROOT/           — GSE 原始数据
  5. rna/tissue_ontologies/     — 专家注释知识库

用法:
    from core.registry import load_master_registry, MasterRegistry

    reg = load_master_registry()
    for ds_id, role in reg.get_dataset_links("41578023"):
        print(ds_id, role)
    for pmid, role in reg.get_paper_links("GSE118614"):
        print(pmid, role)
    orphans = reg.find_orphans()
"""
from __future__ import annotations

import logging
import os
import re
from enum import Enum
from pathlib import Path
from typing import Any, Optional
from collections.abc import Callable

import yaml
from pydantic import BaseModel, ConfigDict, Field

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════
# Enums
# ═══════════════════════════════════════════════════════


class InsightStatus(str, Enum):
    """AI 解读状态."""
    GENERATED = "generated"
    PENDING = "pending"
    FAILED = "failed"
    NO_GEO = "no_geo"
    PENDING_REVIEW = "pending_review"   # 脏数据/需人工审核
    PREPRINT = "preprint"               # 预印本，无 PMID
    PDF_ONLY = "pdf_only"               # 只有原始 PDF，未处理


class DatasetStatus(str, Enum):
    """GSE 数据集状态."""
    DATA_DOWNLOADED = "data_downloaded"
    CONFIG_EXISTS = "config_exists"
    PIPELINE_COMPLETE = "pipeline_complete"
    DATA_NOT_DOWNLOADED = "data_not_downloaded"
    ORPHAN = "orphan"                   # 无关联论文的数据集
    UNKNOWN = "unknown"


class ModalityStatus(str, Enum):
    """Per-modality 运行状态."""
    DATA_DOWNLOADED = "data_downloaded"       # 数据已下载但未配 config
    NOT_CONFIGURED = "not_configured"         # config 缺失
    CONFIG_EXISTS = "config_exists"           # config 文件已存在
    PIPELINE_COMPLETE = "pipeline_complete"     # pipeline 已跑通
    N_A = "n_a"                               # 非管线模态 (bulk/STARR/SuperSeries 容器)
    UNKNOWN = "unknown"


class LinkRole(str, Enum):
    """论文↔数据集关联角色."""
    PRIMARY = "primary"
    VALIDATION = "validation"
    SUPERSERIES = "superseries"
    PART_OF = "part_of"
    RELATED = "related"


class RepositoryType(str, Enum):
    """数据仓库类型."""
    GEO = "geo"
    ARRAYEXPRESS = "arrayexpress"
    SRA = "sra"
    HCA = "hca"
    LOCAL = "local"
    UNKNOWN = "unknown"


class RelationshipType(str, Enum):
    """数据集间关系类型."""
    SUPERSERIES_OF = "superseries_of"
    PART_OF = "part_of"
    RELATED = "related"


# ═══════════════════════════════════════════════════════
# Legacy types (used by run_reproduce.py)
# ═══════════════════════════════════════════════════════


class ExperimentGroup:
    """An experimental sub-grouping within a dataset."""

    def __init__(self, *, group_name, sample_ids, subset_suffix, modality, status,
                 config_path=None, figures=None):
        self.group_name: str = group_name
        self.sample_ids: list[str] = sample_ids
        self.subset_suffix: str = subset_suffix
        self.modality: str = modality
        self.status: str = status
        self.config_path: Optional[str] = config_path
        self.figures: list[str] = figures or []


def _dict_to_exp_group(data: dict[str, Any]) -> ExperimentGroup:
    return ExperimentGroup(
        group_name=data["group_name"],
        sample_ids=data["sample_ids"],
        subset_suffix=data["subset_suffix"],
        modality=data["modality"],
        status=data["status"],
        config_path=data.get("config_path"),
        figures=data.get("figures", []),
    )


def detect_modality(config_path: str) -> str:
    """Detect modality from a pipeline config file.

    Uses a regex scan for ``CFG.modality = "..."`` rather than importing
    the file (avoids side-effects). Returns the modality string or
    ``"unknown"`` if detection fails.
    """
    try:
        with open(config_path, encoding="utf-8") as f:
            content = f.read()
    except Exception:
        return "unknown"

    match = re.search(r"""CFG\.modality\s*=\s*["'](\w+)["']""", content)
    return match.group(1) if match else "unknown"

# ═══════════════════════════════════════════════════════
# Sub-models
# ═══════════════════════════════════════════════════════


class InsightEntry(BaseModel):
    """论文 AI 解读入口."""
    model_config = ConfigDict(extra="forbid")
    status: InsightStatus = InsightStatus.PDF_ONLY
    insights_path: Optional[str] = None
    pdf_raw: Optional[str] = None


class SupplementFile(BaseModel):
    """单个附表文件."""
    model_config = ConfigDict(extra="forbid")
    path: str
    description: str = ""


class SupplementEntry(BaseModel):
    """论文附表."""
    model_config = ConfigDict(extra="forbid")
    source: str = ""
    dir: str = ""
    files: list[SupplementFile] = []
    build_script: Optional[str] = None


class KbSourceEntry(BaseModel):
    """知识库来源."""
    model_config = ConfigDict(extra="forbid")
    kb_id: str
    path: str
    n_cell_types: int = 0
    n_markers: int = 0
    last_audited: Optional[str] = None
    flagged: bool = False


class DatasetConfig(BaseModel):
    """single config for one modality run."""
    model_config = ConfigDict(extra="forbid")
    path: str
    pipeline_status: str = "not_configured"
    exists: bool = True
    experiments: list[dict[str, Any]] = []


class ModalityInfo(BaseModel):
    """单个模态的信息.

    status 取值见 ModalityStatus enum (data_downloaded / not_configured /
    config_exists / pipeline_complete / n_a / unknown).
    """
    model_config = ConfigDict(extra="forbid")
    status: str = "unknown"
    configs: list[DatasetConfig] = []

class DatasetRelationship(BaseModel):
    """数据集间关系（SuperSeries / part_of）。"""
    model_config = ConfigDict(extra="forbid")
    type: RelationshipType
    dataset_id: str


class DatasetEntry(BaseModel):
    """数据集（GSE / E-MTAB / 本地数据）。"""
    model_config = ConfigDict(extra="forbid")
    repository: RepositoryType = RepositoryType.UNKNOWN
    modalities: dict[str, ModalityInfo] = Field(default_factory=dict)
    status: str = "unknown"          # ∈ {data_downloaded, not_downloaded, orphan, unknown}
    type: str = "SingleAccession"     # "SingleAccession" | "SuperSeries"
    non_pipeline: bool = False       # bulk / STARR / SuperSeries containers — pipeline never runs on these
    data_root: str = ""
    dataset_yaml: Optional[str] = None
    relationships: list[DatasetRelationship] = []
    subseries: list[dict[str, str]] = []   # for SuperSeries: [{id: GSE..., note: ""}, ...]
    notes: str = ""
    # ── 数据集补充字段（手工维护在 datasets.yaml 中）──
    species: str = ""
    tissue: str = ""
    data_format: str = ""
    size_desc: str = ""
    parent_series: str = ""
    n_samples: Optional[int] = None
    n_cells: Optional[int] = None
    sample_info: str = ""
    paper_pmids: list[str] = []


class PaperEntry(BaseModel):
    """论文条目。"""
    model_config = ConfigDict(extra="forbid")
    paper_id: str
    slug: str
    pmid: Optional[str] = None
    title: str = ""
    journal: str = ""
    year: str = ""
    first_author: str = ""
    doi: str = ""
    paper_dir: str = ""
    insights: InsightEntry = Field(default_factory=InsightEntry)
    supplements: list[SupplementEntry] = []
    kb_sources: list[KbSourceEntry] = []
    cross_references: dict = Field(
        default_factory=lambda: {"also_cited_by": [], "notes": ""}
    )


class PaperDatasetLink(BaseModel):
    """论文↔数据集关联（M:N 连接）。"""
    model_config = ConfigDict(extra="forbid")
    paper_id: str
    dataset_id: str
    role: LinkRole = LinkRole.PRIMARY


# ═══════════════════════════════════════════════════════
# Master Registry
# ═══════════════════════════════════════════════════════


class MasterRegistry(BaseModel):
    """完整的五域统一登记表。"""
    model_config = ConfigDict(extra="forbid")
    papers: list[PaperEntry] = []
    datasets: dict[str, DatasetEntry] = {}
    links: list[PaperDatasetLink] = []

    # ── 查询 API ────────────────────────────────────

    def get_dataset_links(
        self, paper_id: str,
    ) -> list[tuple[str, LinkRole]]:
        return [
            (ln.dataset_id, ln.role)
            for ln in self.links
            if ln.paper_id == paper_id
        ]

    def get_paper_links(
        self, dataset_id: str,
    ) -> list[tuple[str, LinkRole]]:
        return [
            (ln.paper_id, ln.role)
            for ln in self.links
            if ln.dataset_id == dataset_id
        ]

    def get_paper(self, paper_id: str) -> Optional[PaperEntry]:
        for p in self.papers:
            if p.paper_id == paper_id:
                return p
        return None

    def get_dataset(self, dataset_id: str) -> Optional[DatasetEntry]:
        return self.datasets.get(dataset_id)

    def get_paper_by_pmid(self, pmid: str) -> Optional[PaperEntry]:
        for p in self.papers:
            if p.pmid == pmid:
                return p
        return None

    def get_paper_by_slug(self, slug: str) -> Optional[PaperEntry]:
        for p in self.papers:
            if p.slug == slug:
                return p
        return None

    def find_orphans(self) -> list[tuple[str, DatasetEntry]]:
        linked_ids = {ln.dataset_id for ln in self.links}
        return [
            (ds_id, ds)
            for ds_id, ds in self.datasets.items()
            if ds_id not in linked_ids
        ]

    def find_orphan_supplements(
        self, supplements_root: str = "notes/supplements",
    ) -> list[str]:
        if not os.path.isdir(supplements_root):
            return []
        known_pmids = {p.pmid for p in self.papers if p.pmid}
        orphans: list[str] = []
        for entry in sorted(os.listdir(supplements_root)):
            entry_path = os.path.join(supplements_root, entry)
            if not os.path.isdir(entry_path):
                continue
            if entry not in known_pmids:
                orphans.append(entry)
        return orphans

    def verify(self) -> list[dict[str, Any]]:
        findings: list[dict[str, Any]] = []

        known_papers = {p.paper_id for p in self.papers}
        for ln in self.links:
            if ln.paper_id not in known_papers:
                findings.append({
                    "level": "error",
                    "message": f"Link 引用不存在的 paper_id: {ln.paper_id}",
                    "source": ln.dataset_id,
                })

        known_datasets = set(self.datasets.keys())
        for ln in self.links:
            if ln.dataset_id not in known_datasets:
                findings.append({
                    "level": "error",
                    "message": f"Link 引用不存在的 dataset_id: {ln.dataset_id}",
                    "source": ln.paper_id,
                })

        orphan_ids = {ln.dataset_id for ln in self.links}
        for ds_id, ds in self.datasets.items():
            if ds_id not in orphan_ids and ds.status != DatasetStatus.ORPHAN:
                findings.append({
                    "level": "warn",
                    "message": (
                        f"数据集 {ds_id} 无关联论文但 status={ds.status!r}, "
                        f"应为 'orphan'"
                    ),
                    "source": ds_id,
                })

        for ds_id, ds in self.datasets.items():
            for mod_key, mod_info in ds.modalities.items():
                for cfg in mod_info.configs:
                    abs_path = resolve_path(cfg.path)
                    if not os.path.exists(abs_path):
                        findings.append({
                            "level": "warn",
                            "message": (
                                f"config 路径不存在: {cfg.path} "
                                f"(dataset={ds_id}, modality={mod_key})"
                            ),
                            "source": ds_id,
                        })

        papers_root = "projects/papers"
        for p in self.papers:
            if not p.paper_dir:
                continue
            paper_path = os.path.join(papers_root, p.paper_dir)
            if not os.path.isdir(paper_path):
                findings.append({
                    "level": "warn",
                    "message": (
                        f"paper_dir 不存在: {p.paper_dir} "
                        f"(paper_id={p.paper_id})"
                    ),
                    "source": p.paper_id,
                })

        return findings


# ═══════════════════════════════════════════════════════
# Path resolution
# ═══════════════════════════════════════════════════════

_FUXI_DATA_ROOT_PATTERN = re.compile(r"\{FUXI_DATA_ROOT\}")


def resolve_path(path: str) -> str:
    data_root = os.environ.get("FUXI_DATA_ROOT", "")
    if data_root and _FUXI_DATA_ROOT_PATTERN.search(path):
        return _FUXI_DATA_ROOT_PATTERN.sub(data_root.rstrip("/"), path)
    return path


# ═══════════════════════════════════════════════════════
# YAML I/O — 两种模式:
#   目录模式: path 指向目录, 加载其下 papers.yaml / datasets.yaml / links.yaml
#   文件模式: path 指向单个 .yaml 文件 (向后兼容)
# ═══════════════════════════════════════════════════════


def _serialize_sections(reg: MasterRegistry) -> tuple[list, dict, list]:
    """序列化 registry 的三个独立块。"""
    papers = [
        p.model_dump(exclude_none=True, exclude_defaults=True, mode="json")
        for p in reg.papers
    ]
    datasets = {
        ds_id: ds.model_dump(exclude_none=True, exclude_defaults=True, mode="json")
        for ds_id, ds in reg.datasets.items()
    }
    links = [
        ln.model_dump(exclude_none=True, exclude_defaults=True, mode="json")
        for ln in reg.links
    ]
    return papers, datasets, links


def _dump_yaml(data: Any, path: str) -> None:
    """YAML 写入辅助, block 风格。"""
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(
            data, f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )


def _load_yaml_file(path: str) -> Any:
    """加载单个 YAML 文件，文件不存在返回 None。"""
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except Exception:
        logger.exception("Failed to load %s", path)
        return None


_REGISTRY_DIR = "projects/registry"


def _resolve_registry_path(path: str | None) -> str:
    """解析注册表路径，始终返回目录。"""
    if path is not None:
        return path
    return _REGISTRY_DIR


def load_master_registry(
    path: str | None = None,
) -> MasterRegistry:
    """加载注册表目录下的 papers.yaml / datasets.yaml / links.yaml。"""
    dir_path = _resolve_registry_path(path)
    if not os.path.isdir(dir_path):
        logger.info("Registry dir not found at %s, returning empty", dir_path)
        return MasterRegistry()

    papers = _load_yaml_file(os.path.join(dir_path, "papers.yaml")) or []
    datasets = _load_yaml_file(os.path.join(dir_path, "datasets.yaml")) or {}
    links = _load_yaml_file(os.path.join(dir_path, "links.yaml")) or []

    papers_parsed = [PaperEntry(**p) for p in papers]
    datasets_parsed = {
        ds_id: DatasetEntry(**ds)
        for ds_id, ds in datasets.items()
    }
    links_parsed = [PaperDatasetLink(**ln) for ln in links]

    return MasterRegistry(
        papers=papers_parsed,
        datasets=datasets_parsed,
        links=links_parsed,
    )


def save_master_registry(
    registry: MasterRegistry,
    path: str | None = None,
) -> None:
    """将注册表保存为 3 个子文件到指定目录。"""
    dir_path = _resolve_registry_path(path)
    os.makedirs(dir_path, exist_ok=True)
    papers_data, datasets_data, links_data = _serialize_sections(registry)
    _dump_yaml(papers_data, os.path.join(dir_path, "papers.yaml"))
    _dump_yaml(datasets_data, os.path.join(dir_path, "datasets.yaml"))
    _dump_yaml(links_data, os.path.join(dir_path, "links.yaml"))

# ═══════════════════════════════════════════════════════
# 命令行入口
# ═══════════════════════════════════════════════════════


def _print_report(registry: MasterRegistry, verbose: bool = False) -> None:
    n_papers = len(registry.papers)
    n_datasets = len(registry.datasets)
    n_links = len(registry.links)
    orphans = registry.find_orphans()

    print(f"论文: {n_papers}")
    print(f"数据集: {n_datasets}")
    print(f"关联: {n_links}")
    print(f"孤儿数据集: {len(orphans)}")

    if orphans:
        print("\n\u26a0\ufe0f  孤儿数据集（无关联论文）:")
        for ds_id, ds in orphans:
            status = ds.status
            modalities = ", ".join(ds.modalities.keys()) if ds.modalities else "?"
            print(f"  - {ds_id}  ({modalities}, status={status})")

    supp_orphans = registry.find_orphan_supplements()
    if supp_orphans:
        print(f"\n\u26a0\ufe0f  孤儿附表（目录存在但无对应论文 PMID）:")
        for pmid_dir in supp_orphans:
            print(f"  - notes/supplements/{pmid_dir}/")

    if verbose:
        findings = registry.verify()
        if findings:
            print(f"\n一致性校验发现 {len(findings)} 项:")
            for f in findings:
                print(f"  [{f['level'].upper()}] {f['source']}: {f['message']}")


def _cmd_reset_gse(
    registry: MasterRegistry, dataset_id: str,
    path: str,
) -> MasterRegistry:
    ds = registry.datasets.get(dataset_id)
    if ds is None:
        print(f"\u274c dataset_id '{dataset_id}' 不在 registry 中")
        return registry

    data_root_abs = resolve_path(ds.data_root) if ds.data_root else ""
    orphan_ids = {ln.dataset_id for ln in registry.links}
    is_orphan = dataset_id not in orphan_ids

    if is_orphan:
        ds.status = DatasetStatus.ORPHAN
        print(f"  \u2192 孤儿数据集，保持 status = orphan")
    elif not data_root_abs or not os.path.isdir(data_root_abs):
        ds.status = DatasetStatus.DATA_NOT_DOWNLOADED
        print(f"  \u2192 data_root 不存在，status = data_not_downloaded")
    else:
        ds.status = DatasetStatus.DATA_DOWNLOADED
        print(f"  \u2192 data_root 存在，status = data_downloaded")

    for mod_key, mod_info in ds.modalities.items():
        for cfg in mod_info.configs:
            if os.path.exists(resolve_path(cfg.path)):
                cfg.pipeline_status = "config_exists"
            else:
                cfg.pipeline_status = "not_configured"

    registry.datasets[dataset_id] = ds
    print(f"\u2705 {dataset_id} 已重置")
    return registry


_JOURNAL_ABBREVS = {
    "cell": "cell", "neuron": "neuron", "nature": "nature",
    "nature communications": "natcomms", "nature genetics": "natgenet",
    "cell genomics": "cellgenom", "cell reports": "cellrep",
    "scientific reports": "scirep", "scientific data": "scidata",
    "plos biology": "plosbiol", "plos genetics": "plosgenet",
    "genome biology": "genomebiol", "developmental cell": "devcell",
    "elife": "elife", "iscience": "iscience",
    "frontiers in immunology": "frontimmunol",
    "frontiers in genetics": "frontgenet",
    "proceedings of the national academy of sciences": "pnas",
    "protein & cell": "proteincell", "biorxiv": "biorxiv",
    "research square": "researchsq", "stem cell reports": "stemcellrep",
}


def _build_slug_local(first_author: str, year: str, journal: str, pmid: Optional[str]) -> str:
    author = re.sub(r"[^a-z]", "", (first_author or "unknown").lower())[:20] or "unknown"
    yr = (year or "0000")[:4]
    j = (journal or "").strip().lower()
    j = re.sub(r"\s+", " ", j)
    ab = _JOURNAL_ABBREVS.get(j, "unknown")
    return f"{author}{yr}_{ab}"

def _register_from_insights(
    registry: MasterRegistry,
    insights: dict[str, Any],
    subdir: str,
    dry_run: bool = False,
) -> MasterRegistry:
    """共享登记逻辑：从 insights.yaml 读 meta -> 写入 registry"""
    meta = insights.get("paper_meta", {}) or {}
    pmid = str(meta.get("pmid", "") or "")
    if not pmid:
        # no PMID --- 用 author-year-journal 构造 paper_id
        first_auth = str(meta.get("first_author", "") or "")
        year = str(meta.get("year", "") or "")
        journal = str(meta.get("journal", "") or "")
        pmid = f"no-pmid-{_build_slug_local(first_auth, year, journal, None)}"

    slug = _build_slug_local(
        str(meta.get("first_author", "")),
        str(meta.get("year", "")),
        str(meta.get("journal", "")),
        pmid if not pmid.startswith("no-pmid-") else None,
    )

    if registry.get_paper(pmid):
        print(f"\u26a0\ufe0f Paper {pmid} exists, skip")
        return registry

    paper_entry = PaperEntry(
        paper_id=pmid, slug=slug, pmid=pmid if not pmid.startswith("no-pmid-") else None,
        title=str(meta.get("title", "") or ""),
        journal=str(meta.get("journal", "") or ""),
        year=str(meta.get("year", "") or ""),
        first_author=str(meta.get("first_author", "") or ""),
        doi=str(meta.get("doi", "") or ""),
        paper_dir=subdir,
        insights=InsightEntry(status="generated", insights_path="insights.yaml"),
    )
    geo_ids = (insights.get("data_access", {}) or {}).get("geo_ids", []) or []
    new_links: list[PaperDatasetLink] = []
    for gse_id in geo_ids:
        if gse_id not in registry.datasets:
            registry.datasets[gse_id] = DatasetEntry(
                repository=RepositoryType.GEO,
                status="data_not_downloaded",
                data_root=f"{{FUXI_DATA_ROOT}}/{gse_id}",
            )
        new_links.append(PaperDatasetLink(
            paper_id=pmid, dataset_id=gse_id, role=LinkRole.PRIMARY,
        ))
    registry.papers.append(paper_entry)
    registry.links.extend(new_links)
    print(f"\u2705 Added: {slug}")
    print(f"   {str(meta.get('title', ''))[:100]}")
    print(f"   {meta.get('journal', '')} ({meta.get('year', '')})")
    print(f"   GSEs: {', '.join(geo_ids) or 'none'}")
    if dry_run:
        print("\n\U0001f4a1 --dry-run, not saved")
    return registry


def _cmd_add_paper(
    registry: MasterRegistry,
    pmid: str = "",
    xml: str = "",
    pdf: str = "",
    paper_dir: str = "",
    dry_run: bool = False,
    download: bool = False,
) -> MasterRegistry:
    import subprocess
    import sys

    papers_root = "projects/papers"
    insights = None
    subdir = ""

    # ── PMID / XML / PDF 模式：调 paper_insights ──
    if pmid or xml or pdf:
        cmd = [sys.executable, "core/paper_insights.py"]
        if pmid:
            cmd.extend(["--pmid", pmid])
            print(f"\U0001f50d pmid={pmid}: calling paper_insights ...")
        if xml:
            cmd.extend(["--xml", xml])
            print(f"\U0001f50d xml={xml}: calling paper_insights ...")
        elif pdf:
            cmd.extend(["--pdf", pdf])
            print(f"\U0001f50d pdf={pdf}: calling paper_insights ...")
        result = subprocess.run(cmd, text=True)
        if result.returncode != 0:
            print(f"\u274c paper_insights failed (exit={result.returncode})")
            if result.stderr:
                print(result.stderr[-500:])
            return registry
        # 扫描新生成的 insights.yaml
        for d in sorted(os.listdir(papers_root)):
            ipath = os.path.join(papers_root, d, "insights.yaml")
            if not os.path.isfile(ipath):
                continue
            try:
                with open(ipath) as f:
                    data = yaml.safe_load(f)
                m = data.get("paper_meta", {}) or {}
                # PMID 模式：按 PMID 匹配；XML 模式：匹配第一个
                if pmid and str(m.get("pmid", "")) == pmid:
                    insights, subdir = data, d
                    break
                if xml and insights is None:
                    insights, subdir = data, d
                    break
            except Exception:
                continue
        if insights is None:
            print(f"\u274c No insights.yaml found for PMID={pmid or xml}")
            return registry

    # ── paper-dir 模式：直接从已有目录读取 ──
    elif paper_dir:
        ipath = os.path.join(paper_dir, "insights.yaml")
        if not os.path.isfile(ipath):
            print(f"\u274c {ipath} not found")
            return registry
        try:
            with open(ipath) as f:
                insights = yaml.safe_load(f)
        except Exception as e:
            print(f"\u274c Failed to read {ipath}: {e}")
            return registry
        subdir = os.path.basename(paper_dir.rstrip("/"))

    else:
        print("\u274c Must specify --pmid, --xml, or --paper-dir")
        return registry

    return _register_from_insights(registry, insights, subdir, dry_run)


def _cmd_register_gse(
    registry: MasterRegistry,
    gse_id: str,
    dry_run: bool = False,
) -> MasterRegistry:
    """Register a GSE dataset by SOFT metadata \u2192 PMID link.

    1. Fetch SOFT metadata from NCBI to get PMID(s)
    2. Create DatasetEntry if not already present
    3. Link to each PMID found (if PMID exists in registry)
    """
    from core.geo_downloader import fetch_soft_metadata

    gse_id = gse_id.upper()

    try:
        meta = fetch_soft_metadata(gse_id)
    except Exception as e:
        print(f"\u274c Failed to fetch metadata for {gse_id}: {e}")
        return registry

    title = meta.get("title", "") or ""
    pmid_list = meta.get("pmid", []) or []
    print(f"\n\U0001f50d {gse_id}: {title[:100]}")
    print(f"   PMIDs: {', '.join(pmid_list) or 'none'}")

    if not pmid_list:
        print(f"\u26a0\ufe0f  No PMID found in GEO metadata for {gse_id}")
        print("   Registering dataset without paper link")

    if gse_id not in registry.datasets:
        data_root = f"{{FUXI_DATA_ROOT}}/{gse_id}"
        resolved = resolve_path(data_root)
        status = "data_downloaded" if os.path.isdir(resolved) else "data_not_downloaded"
        registry.datasets[gse_id] = DatasetEntry(
            repository=RepositoryType.GEO,
            status=status,
            data_root=data_root,
        )
        print(f"\u2705  Created dataset entry: {gse_id} ({status})")
    else:
        ds = registry.datasets[gse_id]
        print(f"\u2139\ufe0f  Dataset {gse_id} already exists (status={ds.status})")

    linked = 0
    for pmid in pmid_list:
        paper = registry.get_paper(pmid)
        if paper is None:
            print(f"\u26a0\ufe0f  PMID {pmid} not found in registry - skipping link")
            print(f"       Run: python -m core.registry register --pmid {pmid}")
            continue
        existing = any(
            ln.paper_id == pmid and ln.dataset_id == gse_id
            for ln in registry.links
        )
        if not existing:
            registry.links.append(PaperDatasetLink(
                paper_id=pmid, dataset_id=gse_id, role=LinkRole.RELATED,
            ))
            print(f"\u2705  Linked {gse_id} \u2192 {pmid} ({paper.slug})")
            linked += 1
        else:
            print(f"\u2139\ufe0f  Link {gse_id} \u2192 {pmid} already exists")

    if linked:
        print(f"\n   Total new links: {linked}")
    if dry_run:
        print("\n\U0001f4a1 --dry-run, not saved")
    return registry

def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="五域统一论文登记表 (Master Registry) 工具",
    )
    parser.add_argument(
        "--registry", "-r", default=None,
        help="路径 (目录或 .yaml 文件, 默认自动检测)",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    p_report = sub.add_parser("report", help="打印汇总报告")
    p_report.add_argument("--verbose", "-v", action="store_true")

    sub.add_parser("verify", help="一致性校验")

    p_reset = sub.add_parser("reset-gse", help="重置数据集状态")
    p_reset.add_argument("dataset_id", help="数据集 ID (如 GSE12345)")

    sub.add_parser("find-orphans", help="查找孤儿数据集")

    p_register = sub.add_parser("register", help="注册论文/数据集（--pmid | --gse | --xml | --pdf）")
    p_register.add_argument("--pmid", default=None, help="PubMed ID (NCBI 自动下载)")
    p_register.add_argument("--gse", default=None, help="GEO 数据集 ID (如 GSE164044)")
    p_register.add_argument("--paper-dir", default=None, help="已有 paper 目录 (含 insights.yaml)")
    p_register.add_argument("--xml", dest="xml_path", default=None, help="本地 PMC XML 文件路径")
    p_register.add_argument("--pdf", default=None, help="PDF 文件路径 (pymupdf4llm → md → LLM)")
    p_register.add_argument("--dry-run", action="store_true", help="预览不写入")
    p_register.add_argument("--download", action="store_true",
                            help="Auto-download GSE datasets from NCBI GEO after paper import.")

    p_add = sub.add_parser("add-paper", help="[DEPRECATED] 请改用 register --pmid")
    p_add.add_argument("--pmid", default=None, help="PubMed ID (NCBI 自动下载)")
    p_add.add_argument("--paper-dir", default=None, help="已有 paper 目录 (含 insights.yaml)")
    p_add.add_argument("--xml", dest="xml_path", default=None, help="本地 PMC XML 文件路径")
    p_add.add_argument("--pdf", default=None, help="PDF 文件路径 (pymupdf4llm → md → LLM)")
    p_add.add_argument("--dry-run", action="store_true", help="预览不写入")
    p_add.add_argument("--download", action="store_true",
                        help="Auto-download GSE datasets from NCBI GEO after paper import.")
    args = parser.parse_args()
    reg_path = args.registry if hasattr(args, "registry") else None

    if args.command in ("report", "verify", "reset-gse", "find-orphans", "add-paper", "register"):
        registry = load_master_registry(reg_path)

    if args.command == "report":
        _print_report(registry, verbose=args.verbose)

    elif args.command == "verify":
        findings = registry.verify()
        if not findings:
            print("\u2705 一致性校验通过，无发现")
        else:
            for f in findings:
                icon = {"error": "\u274c", "warn": "\u26a0\ufe0f", "info": "\u2139\ufe0f"}.get(
                    f["level"], "?"
                )
                print(f"  {icon} [{f['level'].upper()}] {f['source']}: {f['message']}")

    elif args.command == "reset-gse":
        registry = _cmd_reset_gse(registry, args.dataset_id, str(reg_path))
        save_master_registry(registry, reg_path)

    elif args.command == "find-orphans":
        orphans = registry.find_orphans()
        if not orphans:
            print("\u2705 无孤儿数据集")
        else:
            print(f"\u26a0\ufe0f  发现 {len(orphans)} 个孤儿数据集:")
            for ds_id, ds in orphans:
                mod_str = ", ".join(ds.modalities.keys()) if ds.modalities else "?"
                print(f"  - {ds_id}  ({mod_str}, status={ds.status})")
        supp_orphans = registry.find_orphan_supplements()
    elif args.command == "register":
        if args.gse:
            registry = _cmd_register_gse(registry, args.gse, dry_run=args.dry_run)
        elif args.pmid or args.xml_path or args.pdf or args.paper_dir:
            registry = _cmd_add_paper(registry, pmid=args.pmid or "", xml=args.xml_path or "",
                                      pdf=args.pdf or "", paper_dir=args.paper_dir or "",
                                      dry_run=args.dry_run, download=args.download)
        else:
            print("\u274c Must specify --pmid, --gse, --xml, --pdf, or --paper-dir")
            return registry
        if not args.dry_run:
            save_master_registry(registry, reg_path)
    elif args.command == "add-paper":
        print("\u26a0\ufe0f  [DEPRECATED] add-paper \u5df2\u5f03\u7528\uff0c\u8bf7\u6539\u7528: register --pmid <PMId>")
        registry = _cmd_add_paper(registry, pmid=args.pmid or "", xml=args.xml_path or "",
                                  pdf=args.pdf or "", paper_dir=args.paper_dir or "",
                                  dry_run=args.dry_run, download=args.download)
        if not args.dry_run:
            save_master_registry(registry, reg_path)


if __name__ == "__main__":
    main()
