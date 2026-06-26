#!/usr/bin/env python3
"""
dataset_schema.py — dataset.yaml Python 数据模型
==================================================

定义 dataset.yaml 文件对应的 Python 数据类，用于:
  - 类型安全的读写操作
  - 与 dataset_detector.py 配合自动生成
  - 与 dataset_validator.py 配合验证完整性

实际格式基于 GEO 数据集的 dataset.yaml 文件。

用法:
    from core.dataset_schema import DatasetMeta, load_dataset, save_dataset
    import os
    data_root = os.environ['FUXI_DATA_ROOT']
    ds = load_dataset(os.path.join(data_root, "your_dataset", "dataset.yaml"))
    print(ds.modalities[0].name)  # "scRNA-seq"
"""

import os
import yaml
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class FileEntry:
    """单个数据文件描述"""
    file: str
    format: str


@dataclass
class SampleEntry:
    """样本描述 — 包含 RNA/ATAC/Spatial 文件列表"""
    id: str
    label: str
    group: Optional[str] = None
    rna: list = field(default_factory=list)      # list[FileEntry]
    atac: list = field(default_factory=list)     # list[FileEntry]
    spatial: list = field(default_factory=list)  # list[FileEntry]
    spots: list = field(default_factory=list)    # spatial-specific
    species: Optional[str] = None                # cross-species datasets
    note: Optional[str] = None


@dataclass
class ModalityEntry:
    """组学类型声明"""
    name: str            # scRNA-seq, scATAC-seq, spatial_transcriptomics, sc_multiome
    status: str          # downloaded, partial, not_downloaded
    format: str
    file_count: int = 0
    total_size_gb: float = 0.0
    subseries: Optional[str] = None
    note: Optional[str] = None


@dataclass
class Comparison:
    """实验比较设计"""
    name: str
    type: str            # condition, time_series, perturbation
    groups: list = field(default_factory=list)


@dataclass
class Resources:
    """外部资源引用"""
    genome: Optional[str] = None
    ortholog_map: Optional[str] = None
    technology: Optional[str] = None


@dataclass
class PipelineStatus:
    """管线运行状态"""
    scRNAseq: Optional[str] = None   # pending, running, completed, failed
    ATACseq: Optional[str] = None
    spatial: Optional[str] = None


@dataclass
class Meta:
    """元数据的元数据"""
    created: Optional[str] = None
    updated: Optional[str] = None
    generated_by: Optional[str] = None
    pipeline_status: PipelineStatus = field(default_factory=PipelineStatus)


@dataclass
class DatasetMeta:
    """完整的 dataset.yaml 数据模型"""
    id: str
    type: str           # SingleAccession, SuperSeries
    title: str
    species: Optional[str] = None
    species_key: Optional[str] = None   # normalised pipeline key (e.g. 'human', 'mouse')
    tissue: Optional[str] = None
    note: Optional[str] = None
    description: Optional[str] = None
    pubmed_id: Optional[str] = None
    parent_superseries: Optional[str] = None

    modalities: list = field(default_factory=list)      # list[ModalityEntry]
    samples: list = field(default_factory=list)          # list[SampleEntry]
    subseries: list = field(default_factory=list)        # SuperSeries only
    comparisons: list = field(default_factory=list)      # list[Comparison]
    resources: Optional[Resources] = None
    meta: Meta = field(default_factory=Meta)


def load_dataset(yaml_path: str) -> DatasetMeta:
    """从 YAML 文件加载数据集元数据"""
    with open(yaml_path, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    # Parse modalities
    modalities = []
    for m in data.get('modalities', []):
        modalities.append(ModalityEntry(
            name=m.get('name', ''),
            status=m.get('status', ''),
            format=m.get('format', ''),
            file_count=m.get('file_count', 0),
            total_size_gb=m.get('total_size_gb', 0.0),
            subseries=m.get('subseries'),
            note=m.get('note'),
        ))

    # Parse samples
    samples = []
    for s in data.get('samples', []):
        samples.append(SampleEntry(
            id=s.get('id', ''),
            label=s.get('label', ''),
            group=s.get('group'),
            rna=s.get('rna', []),
            atac=s.get('atac', []),
            spatial=s.get('spatial', []),
            spots=s.get('spots'),
            species=s.get('species'),
            note=s.get('note'),
        ))

    # Parse comparisons
    comparisons = []
    for c in data.get('comparisons', []):
        comparisons.append(Comparison(
            name=c.get('name', ''),
            type=c.get('type', ''),
            groups=c.get('groups', []),
        ))

    # Parse resources
    res = data.get('resources')
    resources = Resources(
        genome=res.get('genome') if res else None,
        ortholog_map=res.get('ortholog_map') if res else None,
        technology=res.get('technology') if res else None,
    ) if res else None

    # Parse meta
    m = data.get('meta', {})
    ps = m.get('pipeline_status', {})
    meta = Meta(
        created=m.get('created'),
        updated=m.get('updated'),
        generated_by=m.get('generated_by'),
        pipeline_status=PipelineStatus(
            scRNAseq=ps.get('scRNAseq'),
            ATACseq=ps.get('ATACseq'),
            spatial=ps.get('spatial'),
        ),
    )

    return DatasetMeta(
        id=data.get('id', ''),
        type=data.get('type', 'SingleAccession'),
        title=data.get('title', ''),
        species=data.get('species'),
        species_key=data.get('species_key'),
        tissue=data.get('tissue'),
        note=data.get('note'),
        description=data.get('description'),
        pubmed_id=data.get('pubmed_id'),
        parent_superseries=data.get('parent_superseries'),
        modalities=modalities,
        samples=samples,
        subseries=data.get('subseries', []),
        comparisons=comparisons,
        resources=resources,
        meta=meta,
    )


def _sample_to_dict(s):
    """Serialize a SampleEntry to a plain dict for YAML emission."""
    d = {"id": s.id, "label": s.label}
    if s.group is not None:
        d["group"] = s.group
    if s.species is not None:
        d["species"] = s.species
    if s.note is not None:
        d["note"] = s.note
    if s.rna:
        d["rna"] = [{"file": f.file, "format": f.format} if isinstance(f, FileEntry) else f
                     for f in s.rna]
    if s.atac:
        d["atac"] = [{"file": f.file, "format": f.format} if isinstance(f, FileEntry) else f
                      for f in s.atac]
    if s.spatial:
        d["spatial"] = [{"file": f.file, "format": f.format} if isinstance(f, FileEntry) else f
                         for f in s.spatial]
    if s.spots:
        d["spots"] = [{"file": f.file, "format": f.format} if isinstance(f, FileEntry) else f
                       for f in s.spots]
    return d


def save_dataset(ds: DatasetMeta, yaml_path: str) -> None:
    """将 DatasetMeta 保存为 YAML 文件"""
    import yaml
    from datetime import datetime

    data = {
        "id": ds.id,
        "type": ds.type,
        "title": ds.title,
    }
    if ds.species is not None:
        data["species"] = ds.species
    if ds.species_key is not None:
        data["species_key"] = ds.species_key
    if ds.tissue is not None:
        data["tissue"] = ds.tissue
    if ds.note is not None:
        data["note"] = ds.note
    if ds.description is not None:
        data["description"] = ds.description
    if ds.pubmed_id is not None:
        data["pubmed_id"] = ds.pubmed_id
    if ds.parent_superseries is not None:
        data["parent_superseries"] = ds.parent_superseries

    data["modalities"] = [
        {
            "name": m.name,
            "status": m.status,
            "format": m.format,
            "file_count": m.file_count,
            "total_size_gb": m.total_size_gb,
        }
        for m in ds.modalities
    ]

    data["samples"] = [_sample_to_dict(s) for s in ds.samples]

    if ds.subseries:
        data["subseries"] = ds.subseries

    data["comparisons"] = [
        {"name": c.name, "type": c.type, "groups": c.groups}
        for c in ds.comparisons
    ]

    res = ds.resources
    if res is not None:
        data["resources"] = {}
        if res.genome:
            data["resources"]["genome"] = res.genome
        if res.ortholog_map:
            data["resources"]["ortholog_map"] = res.ortholog_map
        if res.technology:
            data["resources"]["technology"] = res.technology

    ps = ds.meta.pipeline_status
    data["meta"] = {
        "created": ds.meta.created or datetime.now().isoformat(),
        "updated": ds.meta.updated or datetime.now().isoformat(),
        "generated_by": ds.meta.generated_by or "fuxi_preprocess",
        "pipeline_status": {
            "scRNAseq": ps.scRNAseq,
            "ATACseq": ps.ATACseq,
            "spatial": ps.spatial,
        },
    }

    os.makedirs(os.path.dirname(yaml_path) or ".", exist_ok=True)
    with open(yaml_path, 'w', encoding='utf-8') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
