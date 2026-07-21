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

import logging
import os
from typing import Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field


class FileEntry(BaseModel):
    """单个数据文件描述"""

    model_config = ConfigDict(extra="forbid")
    file: str
    format: str


class SampleEntry(BaseModel):
    """样本描述 — 包含 RNA/ATAC/Spatial 文件列表"""

    model_config = ConfigDict(extra="forbid")
    id: str
    label: str
    group: Optional[str] = None
    rna: list[FileEntry] = Field(default_factory=list)
    atac: list[FileEntry] = Field(default_factory=list)
    spatial: list[FileEntry] = Field(default_factory=list)
    spots: list[FileEntry] = Field(default_factory=list)
    species: Optional[str] = None
    note: Optional[str] = None


class ModalityEntry(BaseModel):
    """组学类型声明"""

    model_config = ConfigDict(extra="forbid")
    name: str  # scRNA-seq, scATAC-seq, spatial_transcriptomics, sc_multiome
    status: str  # downloaded, partial, not_downloaded
    format: str
    file_count: int = 0
    total_size_gb: float = 0.0
    assay_type: Optional[str] = None  # "scRNAseq" | "snRNAseq" | None
    subseries: Optional[str] = None
    note: Optional[str] = None


class Comparison(BaseModel):
    """实验比较设计"""

    model_config = ConfigDict(extra="forbid")
    name: str
    type: str  # condition, time_series, perturbation
    groups: list[str] = Field(default_factory=list)


class Resources(BaseModel):
    """外部资源引用"""

    model_config = ConfigDict(extra="forbid")
    genome: Optional[str] = None
    ortholog_map: Optional[str] = None
    technology: Optional[str] = None


class PipelineStatus(BaseModel):
    """管线运行状态"""

    model_config = ConfigDict(extra="forbid")
    scRNAseq: Optional[str] = None  # noqa: N815
    ATACseq: Optional[str] = None
    spatial: Optional[str] = None


class Meta(BaseModel):
    """元数据的元数据"""

    model_config = ConfigDict(extra="forbid")
    created: Optional[str] = None
    updated: Optional[str] = None
    generated_by: Optional[str] = None
    pipeline_status: PipelineStatus = Field(default_factory=PipelineStatus)


class DatasetMeta(BaseModel):
    """完整的 dataset.yaml 数据模型"""

    model_config = ConfigDict(extra="forbid")
    id: str
    type: str  # SingleAccession, SuperSeries
    title: str
    species: Optional[str] = None
    species_key: Optional[str] = None  # normalised pipeline key (e.g. 'human', 'mouse')
    tissue: Optional[str] = None
    note: Optional[str] = None
    description: Optional[str] = None
    pubmed_id: Optional[str] = None
    parent_superseries: Optional[str] = None
    assay_type: Optional[str] = None

    modalities: list[ModalityEntry] = Field(default_factory=list)
    samples: list[SampleEntry] = Field(default_factory=list)
    subseries: list[dict[str, str]] = Field(default_factory=list)
    comparisons: list[Comparison] = Field(default_factory=list)
    resources: Optional[Resources] = None
    meta: Meta = Field(default_factory=Meta)


def load_dataset(yaml_path: str) -> DatasetMeta:
    """从 YAML 文件加载数据集元数据"""
    with open(yaml_path, "r", encoding="utf-8") as f:
        return DatasetMeta.model_validate(yaml.safe_load(f))


def save_dataset(ds: DatasetMeta, yaml_path: str) -> None:
    """将 DatasetMeta 保存为 YAML 文件"""
    os.makedirs(os.path.dirname(yaml_path) or ".", exist_ok=True)
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.dump(
            ds.model_dump(exclude_none=True, exclude_defaults=True, by_alias=True),
            f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
            indent=2,
        )


# Modality key mapping for pipeline status updates
_MODALITY_TO_STATUS_KEY = {
    "rna": "scRNAseq",
    "atac": "ATACseq",
    "spatial": "spatial",
}


def update_pipeline_status(yaml_path: Optional[str], modality_key: str, status: str) -> None:
    """Update pipeline status for a given modality in a dataset.yaml file.

    Handles unknown modality gracefully (log warning, return without crash).
    Handles None yaml_path gracefully (log warning, return without crash).
    """
    if yaml_path is None:
        logging.warning("yaml_path is None, skipping pipeline status update")
        return

    if modality_key not in _MODALITY_TO_STATUS_KEY:
        logging.warning(
            f"Unknown modality '{modality_key}'. "
            f"Valid modalities: {list(_MODALITY_TO_STATUS_KEY.keys())}"
        )
        return

    ds = load_dataset(yaml_path)
    status_key = _MODALITY_TO_STATUS_KEY[modality_key]
    # Immutable-style update via model_copy instead of setattr
    new_pipeline_status = ds.meta.pipeline_status.model_copy(update={status_key: status})
    new_meta = ds.meta.model_copy(update={"pipeline_status": new_pipeline_status})
    ds = ds.model_copy(update={"meta": new_meta})
    save_dataset(ds, yaml_path)
