#!/usr/bin/env python3
"""
core/label_transfer.py — 可复用的 Label Transfer 验证工具
===========================================================

对标 sc.tl.ingest，封装 Adult Reference → Fetal Query 的标注验证流程，
输出结构化报告（CSV + JSON），供后续跨验证分析使用。

用法
----
    from core.label_transfer import run_label_transfer

    report = run_label_transfer(
        ref_h5ad="path/to/reference.h5ad",
        query_h5ad="path/to/04_clustered.h5ad",
        ref_label_col="majorclass",
        query_cluster_col="leiden",
        kb_annotations="path/to/cell_type_annotations.csv",
        label_map={"Cone": "Cone_Photoreceptor", ...},
        output_dir="results/tables",
    )
    print(report.summary())          # 核心指标
    print(report.per_cluster())      # per-cluster DataFrame
    report.confusion_matrix()        # 打印混淆矩阵
    report.save_all()                # 写 CSV + JSON
"""

import json
import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc


# ── Public data structures ────────────────────────────────────────────


@dataclass
class LabelTransferReport:
    """Label Transfer 的结构化报告。

    调用 :func:`run_label_transfer` 返回此对象，提供便捷方法查看/保存结果。
    """

    query_adata: "sc.AnnData"  # 已标注转移标签的 query（含原始 obs）
    ref_adata: "sc.AnnData"    # 训练用的 reference 子集（含 PCA/UMAP）
    label_map: dict            # ref 标签名 → KB 细胞类型名
    kb_map: dict               # cluster → KB 细胞类型名
    cluster_col: str           # query 的 cluster 列名
    _ref_label_col: str        # reference 的标签列名（summary 用）

    per_cluster_df: pd.DataFrame = field(init=False)
    cell_agreement: int = field(init=False)
    cell_total: int = field(init=False)
    cluster_agreement: int = field(init=False)
    cluster_total: int = field(init=False)
    confusion: pd.DataFrame = field(init=False)
    confusion_pct: pd.DataFrame = field(init=False)

    def __post_init__(self):
        self._ref_label_col = self._ref_label_col  # ensure stored
        self._compute_agreement()

    # ── internal ──────────────────────────────────────────────────────

    def _compute_agreement(self):
        q = self.query_adata
        cc = self.cluster_col
        kb_to_ref = {v: k for k, v in self.label_map.items()}

        q.obs["kb_cell_type"] = q.obs[cc].map(self.kb_map)
        q.obs["kb_ref_label"] = (
            q.obs["kb_cell_type"].map(kb_to_ref).fillna("Unknown")
        )
        q.obs["transfer_label"] = q.obs["transfer_label"]  # already set
        q.obs["kb_x_transfer_match"] = (
            q.obs["kb_ref_label"] == q.obs["transfer_label"]
        )

        # Per-cluster
        rows = []
        for cluster in sorted(q.obs[cc].unique()):
            mask = q.obs[cc] == cluster
            n = mask.sum()
            kb_label = self.kb_map.get(cluster, "Unknown")
            tc = q.obs.loc[mask, "transfer_label"].value_counts()
            top = tc.index[0]
            top_pct = tc.iloc[0] / n * 100
            second = (
                f"{tc.index[1]}({tc.iloc[1]/n*100:.0f}%)"
                if len(tc) > 1
                else ""
            )
            kb_r = kb_to_ref.get(kb_label, "N/A")
            match_icon = "✓" if top == kb_r else "✗"
            rows.append(
                {
                    "cluster": cluster,
                    "n": n,
                    "kb_label": kb_label,
                    "kb_ref_label": kb_r,
                    "transfer": top,
                    "transfer%": top_pct,
                    "2nd_best": second,
                    "match": match_icon,
                }
            )

        self.per_cluster_df = pd.DataFrame(rows).sort_values("cluster")
        self.cluster_agreement = int((self.per_cluster_df["match"] == "✓").sum())
        self.cluster_total = len(rows)
        self.cell_agreement = int(q.obs["kb_x_transfer_match"].sum())
        self.cell_total = q.n_obs

        # Confusion
        self.confusion = pd.crosstab(
            q.obs["kb_ref_label"],
            q.obs["transfer_label"],
            margins=True,
            margins_name="Total",
            normalize=False,
        )
        self.confusion_pct = (
            pd.crosstab(
                q.obs["kb_ref_label"],
                q.obs["transfer_label"],
                margins=True,
                margins_name="Total",
                normalize="index",
            )
            * 100
        )

    # ── public accessors ──────────────────────────────────────────────

    def summary(self) -> dict:
        """核心指标，适合快速查看。"""
        return {
            "cluster_agreement": {
                "n": self.cluster_agreement,
                "total": self.cluster_total,
                "pct": round(self.cluster_agreement / self.cluster_total * 100, 1),
            },
            "cell_agreement": {
                "n": self.cell_agreement,
                "total": self.cell_total,
                "pct": round(self.cell_agreement / self.cell_total * 100, 1),
            },
            "reference_cell_types": sorted(
                self.ref_adata.obs[self._ref_label_col].unique().tolist()
            ),
            "common_genes": self.ref_adata.n_vars,
        }

    def per_cluster(self) -> pd.DataFrame:
        return self.per_cluster_df

    def print_cluster_report(self):
        """在终端打印 per-cluster 结果。"""
        print(f"\n{'=' * 60}")
        print("PER-CLUSTER LABEL TRANSFER RESULTS")
        print(f"{'=' * 60}")
        for _, r in self.per_cluster_df.iterrows():
            print(
                f"  C{r['cluster']:>4s}  ({r['n']:>5d})  "
                f"KB={r['kb_label']:22s}  →  "
                f"TL={r['transfer']:15s} ({r['transfer%']:5.1f}%) "
                f"{r['match']}  {r['2nd_best']}"
            )
        print(
            f"\nCluster-level agreement: "
            f"{self.cluster_agreement}/{self.cluster_total} "
            f"({self.cluster_agreement/self.cluster_total*100:.1f}%)"
        )
        print(
            f"Cell-level agreement:   "
            f"{self.cell_agreement}/{self.cell_total} "
            f"({self.cell_agreement/self.cell_total*100:.1f}%)"
        )

    def print_mismatches(self):
        """只打印不匹配的 cluster（若有）。"""
        mism = self.per_cluster_df[self.per_cluster_df["match"] == "✗"]
        if len(mism):
            print(f"\nMISMATCHED CLUSTERS ({len(mism)}):")
            for _, r in mism.iterrows():
                print(
                    f"  C{r['cluster']}: KB says {r['kb_label']:22s} → "
                    f"TL says {r['transfer']:15s} "
                    f"({r['transfer%']:.1f}%, 2nd={r['2nd_best']})"
                )

    def print_confusion(self):
        """终端打印混淆矩阵（count + pct）。"""
        print("\nCONFUSION MATRIX (KB rows × TL columns)")
        print(self.confusion.to_string())
        print("\nAs % of row total:")
        print(self.confusion_pct.round(1).to_string())

    def print_type_agreement(self):
        """打印 per-type 一致率。"""
        print("\nPER-TYPE AGREEMENT")
        for ct_name in sorted(self.label_map.keys()):
            mask = self.query_adata.obs["kb_ref_label"] == ct_name
            n = mask.sum()
            if n == 0:
                agree = 0
                pct = 0
            else:
                agree = (
                    self.query_adata.obs.loc[mask, "transfer_label"] == ct_name
                ).sum()
                pct = agree / n * 100
            print(f"  {ct_name:15s}  {agree:>6d}/{n:<6d} ({pct:5.1f}%)")

    # ── file output ────────────────────────────────────────────────────

    def to_dict(self, **extra_meta) -> dict:
        """导出为可序列化的 dict（含 JSON 友好值）。"""
        mismatches = self.per_cluster_df[self.per_cluster_df["match"] == "✗"]
        s = self.summary()
        return {
            "reference": {
                "n_cells": self.ref_adata.n_obs,
                "n_genes": self.ref_adata.n_vars,
                "cell_types": s["reference_cell_types"],
            },
            "query": {
                "n_cells": self.cell_total,
                "common_genes": s["common_genes"],
            },
            "transfer": s,
            "mismatched_clusters": [
                {
                    "cluster": r["cluster"],
                    "kb_label": r["kb_label"],
                    "transfer_label": r["transfer"],
                    "transfer_pct": r["transfer%"],
                }
                for _, r in mismatches.iterrows()
            ],
            "method": "sc.tl.ingest (PCA+UMAP embedding mapping)",
            **extra_meta,
        }

    def write_tables(self, output_dir: str, prefix: str = "label_transfer"):
        """写所有 CSV + JSON 到 output_dir。"""
        os.makedirs(output_dir, exist_ok=True)
        self.per_cluster_df.to_csv(
            os.path.join(output_dir, f"{prefix}_results.csv"), index=False
        )
        self.confusion.to_csv(
            os.path.join(output_dir, f"{prefix}_confusion.csv")
        )
        self.confusion_pct.to_csv(
            os.path.join(output_dir, f"{prefix}_confusion_pct.csv")
        )
        # per-cell mapping
        self.query_adata.obs[
            [
                self.cluster_col,
                "kb_cell_type",
                "kb_ref_label",
                "transfer_label",
                "kb_x_transfer_match",
            ]
        ].to_csv(os.path.join(output_dir, f"{prefix}_per_cell.csv"))

        with open(
            os.path.join(output_dir, f"{prefix}_summary.json"), "w"
        ) as f:
            json.dump(self.to_dict(), f, indent=2)

    def save_annotated_query(self, output_h5ad: str):
        """保存已标注转移标签的 query h5ad。"""
        self.query_adata.write(output_h5ad)


# ── Main entry point ──────────────────────────────────────────────────


def run_label_transfer(
    ref_h5ad: str,
    query_h5ad: str,
    ref_label_col: str,
    *,
    query_cluster_col: str = "leiden",
    kb_annotations: Optional[str] = None,
    kb_cluster_col: str = "cluster",
    kb_label_col: str = "cell_type",
    label_map: Optional[dict] = None,
    pca_comps: int = 100,
    n_neighbors: int = 30,
    n_pcs: int = 50,
    umap_min_dist: float = 0.3,
    embedding_method=("pca", "umap"),
    output_dir: Optional[str] = None,
    output_prefix: str = "label_transfer",
    save_h5ad: Optional[str] = None,
    verbosity: int = 2,
) -> LabelTransferReport:
    """运行 Label Transfer 并返回结构化报告。

    Parameters
    ----------
    ref_h5ad : str
        Reference h5ad 路径（要求含 log-normalised X）。
    query_h5ad : str
        Query h5ad 路径（通常是 04_clustered.h5ad）。
    ref_label_col : str
        Reference obs 中的标签列名（e.g. ``'majorclass'``）。
    query_cluster_col : str
        Query obs 中的 cluster 列名（默认 ``'leiden'``）。
    kb_annotations : str or None
        KB 标注 CSV 路径（需包含 cluster→cell_type 映射）。若为 None，
        报告中的 kb_label 将全部为 ``'N/A'``。
    kb_cluster_col, kb_label_col : str
        KB 标注 CSV 的列名。
    label_map : dict or None
        标签映射 dict: ``{ref_label: kb_type_name}``。
        若为 None，视为 ref_label 与 kb_type 同名。
    pca_comps, n_neighbors, n_pcs, umap_min_dist, embedding_method
        传给 ``sc.pp.pca`` / ``sc.pp.neighbors`` / ``sc.tl.umap`` /
        ``sc.tl.ingest`` 的参数，一般无需修改。
    output_dir : str or None
        若提供，自动写 CSV + JSON 到该目录。
    output_prefix : str
        输出文件名前缀。
    save_h5ad : str or None
        若提供，保存含转移标签的 query h5ad。
    verbosity : int
        scanpy 日志等级（0=silent, 1=warn, 2=info）。

    Returns
    -------
    LabelTransferReport
    """
    sc.settings.verbosity = verbosity

    # ── 1. Load ───────────────────────────────────────────────────────
    ref = sc.read(ref_h5ad)
    query = sc.read(query_h5ad)
    print(f"Reference: {ref.n_obs}c × {ref.n_vars}g")
    print(f"Query:     {query.n_obs}c × {query.n_vars}g (raw={query.raw.n_vars}g)")

    # ── 2. 共有基因子集 ──────────────────────────────────────────────
    common = np.intersect1d(ref.var_names, query.raw.var_names)
    print(f"Common genes: {len(common)}")

    query_raw = query.raw.to_adata()
    query_raw = query_raw[:, common].copy()
    ref_sub = ref[:, common].copy()

    # ── 3. Reference PCA / UMAP ──────────────────────────────────────
    print("Computing reference embedding ...", end=" ", flush=True)
    sc.pp.pca(ref_sub, n_comps=pca_comps, svd_solver="arpack", mask_var=None)
    sc.pp.neighbors(ref_sub, n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.umap(ref_sub, min_dist=umap_min_dist)
    print("Done.")

    # ── 4. Ingest ─────────────────────────────────────────────────────
    print("Running sc.tl.ingest ...", end=" ", flush=True)
    sc.tl.ingest(
        query_raw, ref_sub, obs=ref_label_col, embedding_method=embedding_method
    )
    print("Done.")

    # ══ 整理 obs 列名 ══
    # ingest 添加的列名 = ref_label_col（对于分类 obs，直接添加）
    query_raw.obs["transfer_label"] = query_raw.obs[ref_label_col].astype(str)

    # ── 5. KB 映射 ──────────────────────────────────────────────────
    kb_map: dict = {}
    if kb_annotations:
        ann_df = pd.read_csv(kb_annotations)
        kb_map = dict(zip(ann_df[kb_cluster_col].astype(str), ann_df[kb_label_col]))

    if label_map is None:
        # 默认：ref 标签名 = KB 细胞类型名
        label_map = {
            ct: ct
            for ct in sorted(ref.obs[ref_label_col].unique())
        }

    # ── 6. 构建报告 ──────────────────────────────────────────────────
    report = LabelTransferReport(
        query_adata=query_raw,
        ref_adata=ref_sub,
        label_map=label_map,
        kb_map=kb_map,
        cluster_col=query_cluster_col,
        _ref_label_col=ref_label_col,
    )

    report.print_cluster_report()
    report.print_mismatches()
    report.print_confusion()
    report.print_type_agreement()

    # ── 7. 写出 ──────────────────────────────────────────────────────
    if output_dir:
        report.write_tables(output_dir, prefix=output_prefix)
        print(f"\nTables -> {output_dir}/")
    if save_h5ad:
        report.save_annotated_query(save_h5ad)
        print(f"Query -> {save_h5ad}")

    return report
