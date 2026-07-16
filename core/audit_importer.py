#!/usr/bin/env python3
"""
core/audit_importer.py — 一次性将 FUXI_DATA_ROOT/dataset_audit.md 导入 datasets.yaml

用法:
    source .env && python core/audit_importer.py

将 audit.md 表格中的物种、组织、数据格式、大小、样本/细胞数、父系列、论文 PMID
写入 projects/papers/registry/datasets.yaml 的对应 DatasetEntry 字段。
"""
from __future__ import annotations

import os
import re
import sys
from typing import Any, Optional

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


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


def _dataset_audit_fields(ad: dict) -> dict[str, Any]:
    return {
        k: v
        for k in (
            "species", "tissue", "data_format",
            "size_desc", "parent_series",
            "n_samples", "n_cells", "sample_info",
            "paper_pmids",
        )
        if (v := ad.get(k))
    }


def main() -> None:
    data_root = os.environ.get("FUXI_DATA_ROOT", "")
    if not data_root:
        print("FUXI_DATA_ROOT not set. Source .env first.")
        sys.exit(1)

    audit_path = os.path.join(data_root, "dataset_audit.md")
    if not os.path.isfile(audit_path):
        print(f"audit.md not found at {audit_path}")
        sys.exit(1)

    datasets_path = "projects/papers/registry/datasets.yaml"
    if not os.path.isfile(datasets_path):
        print(f"datasets.yaml not found at {datasets_path}. Run migration first.")
        sys.exit(1)

    print(f"Parsing {audit_path} ...")
    audit_data = _parse_audit_table(audit_path)
    print(f"  {len(audit_data)} entries parsed")

    with open(datasets_path, "r", encoding="utf-8") as f:
        datasets = yaml.safe_load(f) or {}

    updated = 0
    for gse_id, ad in audit_data.items():
        if gse_id in datasets:
            datasets[gse_id].update(_dataset_audit_fields(ad))
            updated += 1

    with open(datasets_path, "w", encoding="utf-8") as f:
        yaml.dump(datasets, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"Updated {updated}/{len(datasets)} datasets in {datasets_path}")


if __name__ == "__main__":
    main()
