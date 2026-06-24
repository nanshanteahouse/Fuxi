#!/usr/bin/env python3
"""
dataset_detector.py — 自动检测数据组学类型
============================================

根据文件名模式自动扫描目录，推断组学类型和样本结构。
用于辅助生成 dataset.yaml。

核心映射表来自 architecture_review 第 5.4 节。

用法:
    python -m core.dataset_detector $FUXI_DATA_ROOT/your_dataset
    python -m core.dataset_detector $FUXI_DATA_ROOT/ --all
"""

import os
import sys
import json
import glob as glob_mod
from collections import defaultdict
from typing import Optional

# ── 核心映射表 ──────────────────────────────────────────────────────────
# (pattern, (modality_name, format_name))
MODALITY_PATTERNS = [
    # scRNA-seq
    ('filtered_feature_bc_matrix.h5',    ('scRNA-seq', '10X_h5')),
    ('features.tsv',                      ('scRNA-seq', 'features')),
    ('matrix.mtx',                        ('scRNA-seq', '10X_mtx')),
    ('count.*mat.*csv',                   ('scRNA-seq', 'csv_matrix')),
    ('_rna_',                             ('scRNA-seq', '10X_h5')),
    ('filtered_feature_bc_matrix',        ('scRNA-seq', '10X_mtx')),
    ('gene_names.txt',                    ('scRNA-seq', 'genes')),
    ('sample_annotations',                ('scRNA-seq', 'metadata')),

    # scATAC-seq
    ('fragments.tsv',                     ('scATAC-seq', 'fragments')),
    ('filtered_peak_bc_matrix',           ('scATAC-seq', '10X_peak_h5')),
    ('_atac_',                            ('scATAC-seq', 'fragments')),
    ('per_barcode_metrics',               ('scATAC-seq', 'metrics')),
    ('peak_annotation',                   ('scATAC-seq', 'annotation')),

    # Spatial
    ('spatial',                           ('spatial_transcriptomics', 'visium')),
    ('visium',                            ('spatial_transcriptomics', 'visium')),
    ('tissue_hires_image',                ('spatial_transcriptomics', 'visium')),
    ('tissue_lowres_image',               ('spatial_transcriptomics', 'visium')),
    ('scalefactors_json',                 ('spatial_transcriptomics', 'visium')),

    # Bulk RNA-seq
    ('_bulk_',                            ('bulk_RNA_seq', 'tsv_counts')),
    ('counts.txt',                        ('bulk_RNA_seq', 'tsv_counts')),
    ('tpm.csv',                           ('bulk_RNA_seq', 'tpm_matrix')),
]


def detect_modality_from_files(file_list: list[str]) -> dict:
    """根据文件名列表推断组学类型。

    返回:
        {modality_name: [matching_files]}
    """
    results = defaultdict(list)
    for fpath in file_list:
        fname = os.path.basename(fpath).lower()
        for pattern, (modality, fmt) in MODALITY_PATTERNS:
            if pattern.lower() in fname:
                results[modality].append(fpath)
                break
    return dict(results)


def scan_directory(directory: str) -> dict:
    """扫描目录，按推断的样本分组文件。

    返回:
        {
            "modalities": [{"name": ..., "format": ..., "file_count": ..., "total_size_gb": ...}],
            "samples": [{"id": ..., "files": ...}],
            "unmatched_files": [...]
        }
    """
    all_files = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.startswith('.') or f == 'dataset.yaml':
                continue
            full_path = os.path.join(root, f)
            rel_path = os.path.relpath(full_path, directory)
            all_files.append(rel_path)

    # Detect modalities
    modality_files = detect_modality_from_files(all_files)

    # Group by sample (naive: group by common prefix in filename)
    samples = defaultdict(list)
    for f in all_files:
        # Extract sample prefix (everything before the second underscore or before _atac/_rna_)
        fname = os.path.basename(f)
        parts = fname.split('_')
        sample_id = '_'.join(parts[:2]) if len(parts) >= 2 else parts[0]
        samples[sample_id].append(f)

    # Calculate sizes
    modality_summary = []
    matched_files = set()
    for modality, files in modality_files.items():
        matched_files.update(files)
        total_size = 0
        for f in files:
            full = os.path.join(directory, f)
            if os.path.exists(full):
                total_size += os.path.getsize(full)
        # Infer format from first matching pattern
        fmt = "unknown"
        for f in files:
            fname = os.path.basename(f).lower()
            for pattern, (mod, f) in MODALITY_PATTERNS:
                if pattern.lower() in fname and mod == modality:
                    fmt = f
                    break
        modality_summary.append({
            "name": modality,
            "format": fmt,
            "file_count": len(files),
            "total_size_gb": round(total_size / 1e9, 2),
        })

    unmatched = [f for f in all_files if f not in matched_files]

    return {
        "modalities": modality_summary,
        "samples": [{"id": sid, "files": flist} for sid, flist in sorted(samples.items())],
        "unmatched_files": unmatched,
    }


def generate_skeleton(directory: str, geo_id: Optional[str] = None) -> str:
    """生成 dataset.yaml 骨架文本。

    返回可写入文件的 YAML 字符串。
    """
    if geo_id is None:
        geo_id = os.path.basename(os.path.abspath(directory))

    result = scan_directory(directory)

    lines = [
        f"id: {geo_id}",
        f"title: \"\"",
        f"species: homo_sapiens",
        f"tissue: unknown",
        "",
        "modalities:",
    ]

    for m in result["modalities"]:
        lines.append(f"  - name: {m['name']}")
        lines.append(f"    status: downloaded")
        lines.append(f"    format: {m['format']}")
        lines.append(f"    file_count: {m['file_count']}")
        lines.append(f"    total_size_gb: {m['total_size_gb']}")

    lines.append("")
    lines.append("samples:")
    for s in result["samples"]:
        lines.append(f"  - id: {s['id']}")
        lines.append(f"    label: \"\"")
        for f in s["files"]:
            lines.append(f"    rna:")
            lines.append(f"      - file: {f}")
            lines.append(f"        format: auto")

    lines.append("")
    lines.append("meta:")
    lines.append("  created: \"\"")
    lines.append("  generated_by: auto_detect")
    lines.append("  pipeline_status: {}")

    return "\n".join(lines)


def main():
    if len(sys.argv) < 2:
        print(f"Usage: python {__file__} <directory> [--all] [--generate-skeleton]")
        print(f"       python {__file__} $FUXI_DATA_ROOT/GSE246169")
        print(f"       python {__file__} $FUXI_DATA_ROOT/ --all --generate-skeleton")
        sys.exit(1)

    directory = sys.argv[1]
    generate = "--generate-skeleton" in sys.argv
    scan_all = "--all" in sys.argv

    if scan_all:
        for entry in sorted(os.listdir(directory)):
            full = os.path.join(directory, entry)
            if os.path.isdir(full):
                result = scan_directory(full)
                print(f"\n{'=' * 60}")
                print(f"  {entry}")
                print(f"{'=' * 60}")
                print(f"  Modalities: {[m['name'] for m in result['modalities']]}")
                print(f"  Samples: {len(result['samples'])}")
                if result['unmatched_files']:
                    print(f"  Unmatched: {len(result['unmatched_files'])} files")
                if generate:
                    skeleton = generate_skeleton(full, entry)
                    print(skeleton)
    else:
        result = scan_directory(directory)
        print(json.dumps(result, indent=2, default=str))
        if generate:
            print("\n--- dataset.yaml skeleton ---\n")
            print(generate_skeleton(directory))


if __name__ == "__main__":
    main()
