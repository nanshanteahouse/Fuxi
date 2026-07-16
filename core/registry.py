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
    """单个数据集配置."""
    model_config = ConfigDict(extra="forbid")
    path: str
    pipeline_status: str = "not_configured"
    exists: bool = True


class ModalityInfo(BaseModel):
    """单个模态的信息."""
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
    status: str = "unknown"
    data_root: str = ""
    dataset_yaml: Optional[str] = None
    relationships: list[DatasetRelationship] = []
    notes: str = ""


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


_REGISTRY_DIR = "projects/papers/registry"


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

    args = parser.parse_args()
    reg_path = args.registry if hasattr(args, "registry") else None

    if args.command in ("report", "verify", "reset-gse", "find-orphans"):
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
        if supp_orphans:
            print(f"\n\u26a0\ufe0f  发现 {len(supp_orphans)} 个孤儿附表目录:")
            for pmid_dir in supp_orphans:
                print(f"  - notes/supplements/{pmid_dir}/")


if __name__ == "__main__":
    main()
