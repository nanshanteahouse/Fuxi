#!/usr/bin/env python3
"""
core/registry_migrate.py — 从 v1 迁移到 v2 统一登记表

从以下源构建 master_registry.yaml:
  1. projects/papers/registry.yaml  (v1 注册表)
  2. projects/{rna,atac,spatial}/  (管线产物目录)
  3. notes/supplements/{PMID}/      (附表目录)
  4. rna/tissue_ontologies/retina/sources/*.yaml  (KB 来源)

用法:
    python core/registry_migrate.py                          # 迁移动
    python core/registry_migrate.py --report                 # 只报告不写文件
    python core/registry_migrate.py --verify                 # 写入后校验
    python core/registry_migrate.py --dry-run                # 模拟运行
"""
from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any, Optional

import yaml

from core.registry import (
    DatasetConfig,
    DatasetEntry,
    DatasetRelationship,
    DatasetStatus,
    InsightEntry,
    InsightStatus,
    KbSourceEntry,
    LinkRole,
    MasterRegistry,
    ModalityInfo,
    PaperDatasetLink,
    PaperEntry,
    RelationshipType,
    RepositoryType,
    SupplementEntry,
    SupplementFile,
    load_master_registry,
    resolve_path,
    save_master_registry,
)

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════
# 路径常量
# ═══════════════════════════════════════════════════════

PAPERS_ROOT = "projects/papers"
SUPPLEMENTS_ROOT = "notes/supplements"
KB_SOURCES_ROOT = "rna/tissue_ontologies/retina/sources"
MODALITY_ROOTS = {
    "rna": "projects/rna",
    "atac": "projects/atac",
    "spatial": "projects/spatial",
}

# Journal slug 缩写映射（用于生成 slug）
_JOURNAL_SLUGS: dict[str, str] = {
    "cell": "cell",
    "neuron": "neuron",
    "nature": "nature",
    "nature communications": "natcomms",
    "nature genetics": "natgenet",
    "cell genomics": "cellgenom",
    "cell reports": "cellrep",
    "scientific reports": "scirep",
    "scientific data": "scidata",
    "plos biology": "plosbiol",
    "plos genetics": "plosgenet",
    "genome biology": "genomebiol",
    "developmental cell": "devcell",
    "nucleic acids research": "nar",
    "stem cell reports": "stemcellrep",
    "stem cells translational medicine": "stemcellstransl",
    "elife": "elife",
    "experimental eye research": "expeyeres",
    "human molecular genetics": "hummolgenet",
    "frontiers in immunology": "frontimmunol",
    "frontiers in genetics": "frontgenet",
    "frontiers in molecular neuroscience": "frontmolneurosci",
    "iscience": "iscience",
    "protein & cell": "proteincell",
    "journal of cellular and molecular medicine": "jcellmolmed",
    "development (cambridge, england)": "development",
    "biorxiv": "biorxiv",
    "research square": "researchsq",
    "proceedings of the national academy of sciences of the united states of america": "pnas",
}


def _normalize_journal(raw: Optional[str]) -> str:
    """标准化期刊名称为缩写 slug。"""
    if not raw or raw == "null" or raw.strip() == "":
        return "unknown"
    key = raw.strip().lower()
    # 去除多余空格
    key = re.sub(r"\s+", " ", key)
    return _JOURNAL_SLUGS.get(key, "unknown")


def _build_slug(
    first_author: str,
    year: str,
    journal_raw: Optional[str],
    pmid: Optional[str],
) -> str:
    """构建人类可读 slug。

    格式: {first_author_lowercase}{year}_{journal_abbrev}
    示例: li2026_natgenet
    """
    author = re.sub(r"[^a-z]", "", (first_author or "unknown").lower())[:20]
    if not author:
        author = "unknown"
    year = year[:4] if year and year != "null" else "0000"
    journal = _normalize_journal(journal_raw)
    return f"{author}{year}_{journal}"


def _build_paper_id(
    pmid: Optional[str],
    first_author: str,
    year: str,
    journal_raw: Optional[str],
) -> str:
    """构建统一主键 paper_id。

    有 PMID → 直接用 PMID
    无 PMID → 用 slug 格式（{author}{year}_{journal}）
    """
    if pmid and pmid not in ("null", ""):
        return pmid
    return _build_slug(first_author, year, journal_raw, pmid)


def _parse_insights_status(raw: Optional[str]) -> InsightStatus:
    """从旧格式解析 InsightStatus。"""
    if not raw or raw == "null":
        return InsightStatus.PDF_ONLY
    mapping = {
        "generated": InsightStatus.GENERATED,
        "pending": InsightStatus.PENDING,
        "failed": InsightStatus.FAILED,
        "no_geo": InsightStatus.NO_GEO,
    }
    return mapping.get(raw.lower(), InsightStatus.PDF_ONLY)


def _modality_root_to_repository(mod_root: str) -> RepositoryType:
    """根据模态根目录推断仓库类型。"""
    return RepositoryType.GEO


def _detect_repository(dataset_id: str) -> RepositoryType:
    """根据 dataset_id 格式推断数据仓库。"""
    if dataset_id.startswith("GSE"):
        return RepositoryType.GEO
    if dataset_id.startswith("E-MTAB-"):
        return RepositoryType.ARRAYEXPRESS
    if dataset_id.startswith("PRJNA"):
        return RepositoryType.SRA
    if dataset_id.startswith("local_"):
        return RepositoryType.LOCAL
    return RepositoryType.UNKNOWN


# ═══════════════════════════════════════════════════════
# Scanner: supplements
# ═══════════════════════════════════════════════════════


def _scan_supplements(pmid: str) -> SupplementEntry:
    """扫描 notes/supplements/{pmid}/ 目录。"""
    supp_dir = os.path.join(SUPPLEMENTS_ROOT, pmid)
    entry = SupplementEntry(
        source="NIHMS" if pmid.startswith("NIHMS") else "Publisher",
        dir=supp_dir,
    )
    if not os.path.isdir(supp_dir):
        return entry

    # 检测 build script
    for fname in sorted(os.listdir(supp_dir)):
        if fname.startswith("build_") and fname.endswith(".py"):
            entry.build_script = fname
            break

    # 检测附件
    for fname in sorted(os.listdir(supp_dir)):
        fpath = os.path.join(supp_dir, fname)
        if not os.path.isfile(fpath):
            continue
        if fname == "ref.txt":
            continue
        entry.files.append(SupplementFile(
            path=fname,
        ))
    return entry


# ═══════════════════════════════════════════════════════
# Scanner: KB sources
# ═══════════════════════════════════════════════════════


def _scan_kb_sources() -> dict[str, dict[str, Any]]:
    """扫描 KB sources 目录，按 pmid 索引。

    Returns:
        {pmid: {kb_id, path, n_cell_types, last_audited, flagged}}
    """
    result: dict[str, dict[str, Any]] = {}
    sources_dir = Path(KB_SOURCES_ROOT)
    if not sources_dir.is_dir():
        return result

    for fpath in sorted(sources_dir.glob("*.yaml")):
        if fpath.name == "schema.yaml":
            continue
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
            meta = data.get("source_meta", {}) if data else {}
            kb_id = meta.get("id", fpath.stem)
            pmid = str(meta.get("pmid", "")) if meta.get("pmid") else ""
            if not pmid:
                continue

            # 统计 cell types
            markers = data.get("markers", {}) if data else {}
            n_types = len(markers) if markers else 0

            # 取得最新的 audit 日期
            last_audited = ""
            for ct_name, ct_data in markers.items():
                audit = ct_data.get("audit", {}) if isinstance(ct_data, dict) else {}
                audited_raw = audit.get("last_audited", "") or ""
                audited_str = str(audited_raw) if not isinstance(audited_raw, str) else audited_raw
                if audited_str > last_audited:
                    last_audited = audited_str

            result[pmid] = {
                "kb_id": kb_id,
                "path": str(fpath.relative_to(".").as_posix()),
                "n_cell_types": n_types,
                "n_markers": _count_markers(markers),
                "last_audited": last_audited or None,
                "flagged": _any_flagged(markers),
            }
        except Exception as e:
            logger.warning("Failed to parse KB source %s: %s", fpath, e)

    return result


def _count_markers(markers: dict) -> int:
    """粗略统计 markers 数量。"""
    total = 0
    for ct_data in markers.values():
        if not isinstance(ct_data, dict):
            continue
        for key in ("confirm", "add"):
            genes = ct_data.get(key, [])
            if isinstance(genes, list):
                total += len(genes)
    return total


def _any_flagged(markers: dict) -> bool:
    """检查是否有任何 cell type 被标记为 flagged。"""
    for ct_data in markers.values():
        if not isinstance(ct_data, dict):
            continue
        audit = ct_data.get("audit", {}) if isinstance(ct_data, dict) else {}
        if isinstance(audit, dict) and audit.get("flagged"):
            return True
    return False


# ═══════════════════════════════════════════════════════
# Scanner: 管线产物目录
# ═══════════════════════════════════════════════════════


def _scan_project_gses() -> dict[str, dict[str, Any]]:
    """扫描 projects/{rna,atac,spatial}/{GSE_ID}/ 目录。

    Returns:
        {gse_id: {modality, configs: [{path, exists}], data_root_exists}}
    """
    result: dict[str, dict[str, Any]] = {}
    for modality, root in MODALITY_ROOTS.items():
        if not os.path.isdir(root):
            continue
        for entry in sorted(os.listdir(root)):
            entry_path = os.path.join(root, entry)
            if not os.path.isdir(entry_path):
                continue
            # 跳过非 GSE 命名的目录
            if not re.match(r"^(GSE|E-MTAB|PRJNA)\d+", entry):
                continue

            if entry not in result:
                result[entry] = {"modalities": {}, "configs": []}

            mod_info: dict[str, Any] = {"configs": []}
            # 找到所有 config 文件
            for fname in os.listdir(entry_path):
                if fname.endswith(".yaml") and "config" in fname.lower():
                    cfg_path = f"projects/{modality}/{entry}/{fname}"
                    mod_info["configs"].append({
                        "path": cfg_path,
                        "exists": os.path.isfile(os.path.join(entry_path, fname)),
                    })
                elif fname == "dataset.yaml":
                    pass  # 这个是在 FUXI_DATA_ROOT 下的

            # 检测 data_root 存在性
            data_root = os.environ.get("FUXI_DATA_ROOT", "")
            data_root_exists = False
            if data_root:
                dr_path = os.path.join(data_root, entry)
                data_root_exists = os.path.isdir(dr_path)

            result[entry]["modalities"][modality] = mod_info
            result[entry]["data_root_exists"] = data_root_exists

    return result


# ═══════════════════════════════════════════════════════
# Scanner: 待整理 PDF
# ═══════════════════════════════════════════════════════


def _scan_pending_pdfs() -> dict[str, str]:
    """扫描 projects/papers/待整理/ 中的 PDF，建立 PMID→filename 映射。

    Returns:
        {pmid: filename}
    """
    pending_dir = os.path.join(PAPERS_ROOT, "待整理")
    result: dict[str, str] = {}
    if not os.path.isdir(pending_dir):
        return result

    for fname in os.listdir(pending_dir):
        if not fname.endswith(".pdf"):
            continue
        # 尝试从文件名提取 PMID（纯数字）
        stem = os.path.splitext(fname)[0]
        if stem.isdigit():
            result[stem] = fname
    return result


def _parse_audit_table(audit_path: str) -> dict[str, dict[str, Any]]:
    if not os.path.isfile(audit_path):
        return {}
    with open(audit_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    result: dict[str, dict[str, Any]] = {}
    in_table = False
    for line in lines:
        s = line.strip()
        if not s.startswith("|"):
            continue
        if s.startswith("| **GSE"):
            in_table = True
        if not in_table or "---" in s:
            continue
        cells = [c.strip() for c in s.split("|")[1:-1]]
        if len(cells) < 11:
            continue
        gse_raw = cells[0].replace("**", "").strip()
        if not gse_raw.startswith("GSE"):
            continue
        species = cells[1].replace("*", "").strip()
        tissue = cells[2].replace("*", "").strip()
        data_format = cells[4].strip()
        size_desc = cells[5].strip() if cells[5] != "—" else ""
        parent_raw = cells[9].strip()
        notes_raw = cells[10].strip()
        paper_raw = cells[11].strip()
        pmids = re.findall(r"\[(\d+)\]\(https://pubmed[^)]+\)", paper_raw)
        n_samples, n_cells, sinfo = _parse_sample_info(notes_raw)
        parent_series = ""
        if parent_raw and parent_raw != "—":
            m = re.search(r"(GSE\d+)", parent_raw)
            parent_series = m.group(1) if m else ""
        result[gse_raw] = {
            "species": species, "tissue": tissue,
            "data_format": data_format, "size_desc": size_desc,
            "parent_series": parent_series,
            "n_samples": n_samples, "n_cells": n_cells,
            "sample_info": sinfo or notes_raw,
            "paper_pmids": list(dict.fromkeys(pmids)),
        }
    return result


def _parse_sample_info(text: str) -> tuple[Optional[int], Optional[int], str]:
    ns, nc = None, None
    m = re.search(r"(\d+)\s*(?:供体|donors?|samples?|样本|eyes?|时间点)", text, re.I)
    if m:
        ns = int(m.group(1))
    m = re.search(r"[>~]?(\d{1,3}(?:,\d{3})*(?:\.\d+)?)\s*[Kk]\b", text)
    if m:
        nc = int(float(m.group(1).replace(",", "")) * 1000)
    else:
        m = re.search(r"(\d{1,3}(?:,\d{3})*)\s*(?:cells?|核|nuclei)", text, re.I)
        if m:
            nc = int(m.group(1).replace(",", ""))
    return ns, nc, text


# ═══════════════════════════════════════════════════════


def _load_v1_registry(path: str = "projects/papers/registry.yaml") -> list[dict[str, Any]]:
    """加载旧版 registry.yaml（v1 格式）。"""
    from core.paper_registry_models import load_registry
    data = load_registry(path)
    return data.get("papers", [])


def _sanitize_str(val: Any) -> str:
    """清理字符串值，处理 null/None 和垃圾数据。"""
    if val is None:
        return ""
    s = str(val).strip()
    # 检测 pymupdf4llm 污染
    if s.startswith("pymupdf4llm_"):
        return ""
    if s.lower() in ("null", "none", ""):
        return ""
    return s


def _v1_to_paper_entry(
    v1_paper: dict[str, Any],
    kb_index: dict[str, dict[str, Any]],
    pending_pdfs: dict[str, str],
) -> Optional[PaperEntry]:
    """将 v1 registry 中的单条 paper dict 转换为 PaperEntry。"""
    pmid = _sanitize_str(v1_paper.get("pmid"))
    first_author = v1_paper.get("first_author", "") or ""
    year = str(v1_paper.get("year", "") or "")
    journal_raw = v1_paper.get("journal")
    title = _sanitize_str(v1_paper.get("title"))
    paper_dir = v1_paper.get("paper_dir", "") or ""
    doi = _sanitize_str(v1_paper.get("doi"))

    # 跳过完全空白的条目
    if not pmid and not title and not paper_dir:
        return None

    # paper_id: 有 PMID 就用 PMID
    paper_id = pmid if pmid else _build_paper_id(None, first_author, year, journal_raw)
    slug = _build_slug(first_author, year, journal_raw, pmid)

    # insights
    raw_status = _sanitize_str(v1_paper.get("insights_status"))
    insights_status = _parse_insights_status(raw_status)

    # 如果 title 是空（污染数据），标记
    if not title and paper_dir:
        if insights_status in (InsightStatus.GENERATED, InsightStatus.PENDING):
            insights_status = InsightStatus.PENDING_REVIEW

    insights = InsightEntry(
        status=insights_status,
        insights_path="insights.yaml" if paper_dir else None,
        pdf_raw=pending_pdfs.get(pmid) if pmid else None,
    )

    # KB sources
    kb_list: list[KbSourceEntry] = []
    if pmid in kb_index:
        kb_info = kb_index[pmid]
        kb_list.append(KbSourceEntry(
            kb_id=kb_info["kb_id"],
            path=kb_info["path"],
            n_cell_types=kb_info["n_cell_types"],
            n_markers=kb_info["n_markers"],
            last_audited=kb_info["last_audited"],
            flagged=kb_info["flagged"],
        ))

    # supplements
    supp_list: list[SupplementEntry] = []
    if pmid:
        supp = _scan_supplements(pmid)
        if supp.dir and (supp.files or supp.build_script):
            supp_list.append(supp)

    return PaperEntry(
        paper_id=paper_id,
        slug=slug,
        pmid=pmid if pmid else None,
        title=title,
        journal=_sanitize_str(journal_raw),
        year=year,
        first_author=first_author,
        doi=doi,
        paper_dir=paper_dir,
        insights=insights,
        supplements=supp_list,
        kb_sources=kb_list,
    )


def _build_datasets_and_links(
    v1_papers: list[dict[str, Any]],
    scanned_gses: dict[str, dict[str, Any]],
):
    """从 v1 registry + 扫描结果构建 datasets 和 links。

    处理:
      - M:N 关系（多个论文共享 GSE）
      - 多模态（一个 GSE 可能 rna+atac 都有）
      - 孤儿数据集（scanned_gses 中有但不在 registry 中）

    Returns:
        (datasets_dict, links_list, orphan_gse_ids)
    """
    datasets: dict[str, DatasetEntry] = {}
    links: list[PaperDatasetLink] = []
    all_linked_gses: set[str] = set()

    for v1_paper in v1_papers:
        pmid = _sanitize_str(v1_paper.get("pmid"))
        first_author = v1_paper.get("first_author", "") or ""
        year = str(v1_paper.get("year", "") or "")
        journal_raw = v1_paper.get("journal")

        paper_id = pmid if pmid else _build_paper_id(None, first_author, year, journal_raw)

        v1_datasets = v1_paper.get("datasets", [])
        if not v1_datasets:
            continue

        for v1_ds in v1_datasets:
            gse_id = v1_ds.get("gse_id", "")
            if not gse_id:
                continue

            all_linked_gses.add(gse_id)

            # 构建 link
            links.append(PaperDatasetLink(
                paper_id=paper_id,
                dataset_id=gse_id,
                role=LinkRole.PRIMARY,
            ))

            # 如果 datasets 已存在（被另一篇论文共享），跳过重复定义
            if gse_id in datasets:
                continue

            # 推断仓库类型
            repo = _detect_repository(gse_id)

            # 从 v1 数据提取
            v1_status = v1_ds.get("status", "unknown") or "unknown"
            v1_notes = v1_ds.get("notes", "") or ""
            v1_modality = v1_ds.get("modality", "rna") or "rna"
            v1_config_path = v1_ds.get("config_path", "") or ""

            # 多模态构建
            modalities: dict[str, ModalityInfo] = {}
            # 从 v1 获取一个 modality
            mod_configs: list[DatasetConfig] = []
            if v1_config_path:
                mod_configs.append(DatasetConfig(
                    path=v1_config_path,
                    pipeline_status=v1_status,
                    exists=os.path.exists(v1_config_path),
                ))

            modalities[v1_modality] = ModalityInfo(
                status=v1_status,
                configs=mod_configs,
            )

            # 合并 scanned 数据（补充 configs 可能遗漏的 config 文件）
            scanned = scanned_gses.get(gse_id, {})
            scanned_mods = scanned.get("modalities", {})
            for smod, sinfo in scanned_mods.items():
                if smod not in modalities:
                    modalities[smod] = ModalityInfo(status="data_downloaded")
                for scfg in sinfo.get("configs", []):
                    # 去重
                    existing_paths = {c.path for c in modalities[smod].configs}
                    if scfg["path"] not in existing_paths:
                        modalities[smod].configs.append(DatasetConfig(
                            path=scfg["path"],
                            pipeline_status="not_configured",
                            exists=scfg["exists"],
                        ))

            # 检测 data_root
            data_root = f"{{FUXI_DATA_ROOT}}/{gse_id}"
            data_root_exists = scanned.get("data_root_exists", False)
            # 对所有数据集检测 FUXI_DATA_ROOT
            if not data_root_exists:
                fdr = os.environ.get("FUXI_DATA_ROOT", "")
                if fdr:
                    data_root_exists = os.path.isdir(os.path.join(fdr, gse_id))
            ds_status = DatasetStatus.DATA_DOWNLOADED if data_root_exists else DatasetStatus.DATA_NOT_DOWNLOADED

            # 如果有 config 存在，升级状态
            has_valid_config = any(
                cfg.exists
                for mod_info in modalities.values()
                for cfg in mod_info.configs
            )
            if has_valid_config:
                ds_status = DatasetStatus.CONFIG_EXISTS

            datasets[gse_id] = DatasetEntry(
                repository=repo,
                modalities=modalities,
                status=ds_status.value,
                data_root=data_root,
                notes=v1_notes if v1_notes else ""
            )

    # 孤儿数据集检测
    orphan_gses = []
    for gse_id in scanned_gses:
        if gse_id not in all_linked_gses:
            orphan_gses.append(gse_id)
            scanned = scanned_gses[gse_id]
            # 构建 modaliies
            mods: dict[str, ModalityInfo] = {}
            for smod, sinfo in scanned.get("modalities", {}).items():
                mods[smod] = ModalityInfo(
                    status="data_downloaded",
                    configs=[
                        DatasetConfig(path=cfg["path"], exists=cfg["exists"])
                        for cfg in sinfo.get("configs", [])
                    ],
                )
            data_root = f"{{FUXI_DATA_ROOT}}/{gse_id}"
            datasets[gse_id] = DatasetEntry(
                repository=RepositoryType.GEO,
                modalities=mods,
                status=DatasetStatus.ORPHAN.value,
                data_root=data_root,
                notes="由迁移脚本自动发现的孤儿数据集"
            )
    return datasets, links, orphan_gses


def build_master_registry(
    v1_registry_path: str = "projects/papers/registry.yaml",
    dry_run: bool = False,
    report_only: bool = False,
) -> MasterRegistry:
    """构建完整的主注册表。

    Args:
        v1_registry_path: v1 registry.yaml 路径
        dry_run: 如果为 True，不写文件
        report_only: 如果为 True，只打印报告不写文件

    Returns:
        构建好的 MasterRegistry
    """
    print("🔍 扫描源数据...")

    # 1. 加载 v1 registry
    v1_papers = _load_v1_registry(v1_registry_path)
    print(f"  ✓ v1 registry: {len(v1_papers)} 篇论文")

    # 2. 扫描 KB sources
    kb_index = _scan_kb_sources()
    print(f"  ✓ KB sources: {len(kb_index)} 个来源")

    # 3. 扫描待整理 PDF
    pending_pdfs = _scan_pending_pdfs()
    print(f"  ✓ 待整理 PDF: {len(pending_pdfs)} 个")

    # 4. 扫描管线产物目录
    scanned_gses = _scan_project_gses()
    print(f"  ✓ 管线目录: {len(scanned_gses)} 个数据集目录")

    # 5. 转换 papers
    print("\n📄 转换论文条目...")
    papers: list[PaperEntry] = []
    v1_skipped = 0
    for v1_paper in v1_papers:
        entry = _v1_to_paper_entry(v1_paper, kb_index, pending_pdfs)
        if entry is None:
            v1_skipped += 1
            continue
        papers.append(entry)
    print(f"  ✓ {len(papers)} 篇转换成功" +
          (f", {v1_skipped} 条跳过" if v1_skipped else ""))

    # 6. 构建 datasets + links
    print("\n🔗 构建数据集与关联...")
    datasets, links, orphan_gses = _build_datasets_and_links(v1_papers, scanned_gses)
    print(f"  ✓ {len(datasets)} 个数据集 (其中 {len(orphan_gses)} 个孤儿)")
    print(f"  ✓ {len(links)} 条关联")

    # 7. 组装
    registry = MasterRegistry(
        papers=papers,
        datasets=datasets,
        links=links,
    )

    # 8. 校验
    print("\n🔎 一致性校验...")
    findings = registry.verify()
    has_errors = any(f["level"] == "error" for f in findings)
    for f in findings:
        level_icon = {"error": "❌", "warn": "⚠️", "info": "ℹ️"}.get(f["level"], "?")
        print(f"  {level_icon} [{f['level'].upper()}] {f['source']}: {f['message']}")

    # 9. 检查附表孤儿
    supp_orphans = registry.find_orphan_supplements()
    if supp_orphans:
        print(f"\n⚠️  孤儿附表目录: {len(supp_orphans)} 个")

    if not findings and not supp_orphans:
        print("\n✅ 校验通过，无发现问题")

    # 10. 报告
    print(f"\n{'='*50}")
    print(f"📊 汇总")
    print(f"    论文: {len(papers)}")
    print(f"    数据集: {len(datasets)}  (其中 {len(orphan_gses)} 个孤儿)")
    print(f"    关联: {len(links)}")
    print(f"    KB 来源: {len(kb_index)}")
    print(f"    附表孤儿: {len(supp_orphans)}")
    if has_errors:
        print(f"\n❌ 存在 ERROR 级别的问题，请检查后重试")
    print(f"{'='*50}")

    return registry


# ═══════════════════════════════════════════════════════
# 命令行入口
# ═══════════════════════════════════════════════════════


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="迁移 v1 registry.yaml → v2 master_registry.yaml",
    )
    parser.add_argument(
        "--v1-path", default="projects/papers/registry.yaml",
        help="v1 registry.yaml 路径",
    )
    parser.add_argument(
        "--output", default="projects/papers/registry",
        help="输出路径 (默认: %(default)s)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="模拟运行，不写文件",
    )
    parser.add_argument(
        "--report", action="store_true",
        help="只打印报告，不写文件",
    )
    parser.add_argument(
        "--verify", action="store_true",
        help="迁移后校验并报告",
    )

    args = parser.parse_args()

    reg = build_master_registry(
        v1_registry_path=args.v1_path,
        dry_run=args.dry_run,
        report_only=args.report,
    )

    if args.report or args.dry_run:
        print("\n💡 未写入文件 (dry-run / report 模式)")

    if args.report:
        return

    if not args.dry_run:
        save_master_registry(reg, args.output)
        print(f"\n💾 已写入: {args.output}/"
              f"{{papers,datasets,links}}.yaml")
        if args.verify:
            loaded = load_master_registry(args.output)
            findings = loaded.verify()
            if findings:
                print(f"\n⚠️  二次校验发现 {len(findings)} 项:")
                for f in findings:
                    print(f"  [{f['level'].upper()}] {f['source']}: {f['message']}")
            else:
                print("\n✅ 二次校验通过")


if __name__ == "__main__":
    main()
