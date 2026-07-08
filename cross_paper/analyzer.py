#!/usr/bin/env python3
"""
Generic Cross-Paper Comparison Analyzer
========================================
Loads annotated h5ad files from multiple pipeline runs, computes pathway
gene-set activity per cell type, supports per-condition splitting
(sample/stage/treatment), and generates structured comparison reports.

Not restricted to RA signaling — any pathway with defined gene sets works.

Usage:
    >>> from cross_paper.analyzer import CrossPaperAnalyzer, DatasetEntry
    >>> analyzer = CrossPaperAnalyzer(
    ...     gene_sets={
    ...         "synthesis": {"label": "Synthesis", "genes": ["ALDH1A1", ...]},
    ...         "receptors": {"label": "Receptors", "genes": ["RARA", ...]},
    ...     },
    ...     tf_list=["RARA", "RARB", ...],
    ...     datasets=[DatasetEntry(...), ...],
    ...     individual_genes=["OPN1SW", "OPN1MW", ...],
    ...     comparison_metrics={"ML_ratio": ...},
    ... )
    >>> analyzer.run()
    >>> analyzer.compare_conditions()
    >>> analyzer.to_csv("comparison.csv")

CLI usage via project presets:
    >>> analyzer = CrossPaperAnalyzer.from_yaml(
    ...     "cross_paper/pathway_config.yaml",
    ...     "cross_paper/dataset_registry.yaml",
    ... )
"""

from __future__ import annotations

import json
import os
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────────────
#  Data types
# ──────────────────────────────────────────────────────────────────────


@dataclass
class DatasetEntry:
    """A single dataset (or dataset subset) to include in comparison.

    When ``split_col`` is set, analysis is performed per-split within
    the same dataset (e.g., per sample, per treatment group).  When
    ``condition_map`` is provided, raw values are translated to
    human-readable condition labels.
    """

    label: str  # display name, e.g. "Wohl2026_D140"
    paper: str  # paper label, e.g. "Wohlschlegel 2026"
    h5ad_path: str  # path to 05_annotated.h5ad (absolute or project-relative)
    grn_table: Optional[str] = None  # path to GRN tf_activity_per_cell_type.csv
    cell_type_col: str = "cell_type"  # obs column for cell type labels
    condition_label: str = ""  # fallback label when not splitting
    split_col: Optional[str] = None  # obs column to split by (e.g. "sample")
    condition_map: dict[str, str] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    _adata: object = field(default=None, repr=False, init=False)


@dataclass
class CellTypeResult:
    """Per-cell-type pathway metrics for one dataset × condition."""

    dataset: str
    paper: str
    cell_type: str
    condition: str
    n_cells: int
    cell_pct: float
    # Gene set scores (mean expression of all genes in set)
    gene_set_scores: dict[str, float] = field(default_factory=dict)
    # Individual gene expression
    individual_expr: dict[str, float] = field(default_factory=dict)
    # Derived comparison metrics
    derived_metrics: dict[str, float] = field(default_factory=dict)
    # GRN TF activity
    grn: dict[str, float] = field(default_factory=dict)


# ──────────────────────────────────────────────────────────────────────
#  Core Analyzer
# ──────────────────────────────────────────────────────────────────────


class CrossPaperAnalyzer:
    """Load multiple datasets and compare pathway activity.

    Parameters
    ----------
    gene_sets : dict[str, dict]
        Nested dict::

            {
                "synthesis": {"label": "RA Synthesis", "genes": ["ALDH1A1", ...]},
                "degradation": {"label": "RA Degradation", "genes": ["CYP26A1", ...]},
            }

    tf_list : list[str]
        Transcription factor names to extract from GRN tables.
    datasets : list[DatasetEntry]
        One entry per dataset group.
    individual_genes : list[str], optional
        Additional individual gene expressions to track (not grouped).
    comparison_metrics : dict[str, dict], optional
        Derived metrics computed from individual genes::

            {
                "ML_opsin_ratio": {
                    "label": "M/L opsin ratio",
                    "formula": "(OPN1MW + OPN1LW) / (OPN1SW + OPN1MW + OPN1LW + RHO)"
                }
            }

        Formulas support +, -, *, /, ( ) and reference individual gene names
        as variables.
    split_col : str or None
        If set, auto-split all datasets by this obs column (e.g. "sample").
        Overridden by per-dataset ``DatasetEntry.split_col`` if present.
    """

    def __init__(
        self,
        gene_sets: dict[str, dict],
        tf_list: list[str],
        datasets: Sequence[DatasetEntry],
        individual_genes: Optional[list[str]] = None,
        comparison_metrics: Optional[dict[str, dict]] = None,
        split_col: Optional[str] = None,
    ):
        self.gene_sets = gene_sets
        self.tf_list = list(tf_list)
        self.datasets = list(datasets)
        self.individual_genes = individual_genes or []
        self.comparison_metrics = comparison_metrics or {}
        self.split_col = split_col
        self._results: list[CellTypeResult] = []
        self._gene_expr_cache: dict[str, pd.DataFrame] = {}

        # All genes we care about (for expression matrix extraction)
        self._all_genes: set[str] = set()
        for gs in self.gene_sets.values():
            self._all_genes.update(gs.get("genes", []))
        self._all_genes.update(self.individual_genes)

        # Validate comparison metrics formulas
        self._validate_metrics()

    def _validate_metrics(self):
        """Simple validation: metric formulas only use known gene names + basic ops."""
        allowed = self._all_genes | {"+", "-", "*", "/", "(", ")", " ", "nan"}
        for name, cfg in self.comparison_metrics.items():
            formula = cfg.get("formula", "")
            tokens = set()
            for tok in formula.replace("(", " ").replace(")", " ").replace(
                "+", " "
            ).replace("-", " ").replace("*", " ").replace("/", " ").split():
                tok = tok.strip()
                if tok:
                    tokens.add(tok)
            unknown = tokens - self._all_genes
            if unknown:
                raise ValueError(
                    f"Metric '{name}' formula uses unknown genes: {unknown}"
                )

    # ── YAML factory ─────────────────────────────────────────────────

    @classmethod
    def from_yaml(
        cls,
        pathway_yaml: str,
        registry_yaml: str,
    ) -> CrossPaperAnalyzer:
        """Create analyzer from YAML config files.

        Parameters
        ----------
        pathway_yaml : str
            Path to pathway gene set YAML (project-relative or absolute).
        registry_yaml : str
            Path to dataset registry YAML.

        Returns
        -------
        CrossPaperAnalyzer
        """
        from cross_paper import load_yaml

        repo_root = Path(__file__).resolve().parent.parent

        def _resolve(p: str) -> str:
            path = Path(p)
            return str(path if path.is_absolute() else repo_root / p)

        pw = load_yaml(pathway_yaml)
        reg = load_yaml(registry_yaml)

        # Select first pathway (or allow specifying pathway name later)
        pathway_name = list(pw.keys())[0]
        pathway = pw[pathway_name]

        gene_sets = pathway.get("gene_sets", {})
        tf_list = pathway.get("tf_list", [])
        individual_genes = pathway.get("individual_genes", [])
        comparison_metrics = pathway.get("comparison_metrics", {})

        # Build DatasetEntry list from registry
        datasets = []
        for entry in reg.get("datasets", []):
            datasets.append(
                DatasetEntry(
                    label=entry["label"],
                    paper=entry["paper"],
                    h5ad_path=_resolve(entry["h5ad"]),
                    grn_table=_resolve(entry["grn_table"])
                    if entry.get("grn_table")
                    else None,
                    cell_type_col=entry.get("cell_type_col", "cell_type"),
                    condition_label=entry.get("condition_label", ""),
                    split_col=entry.get("split_col"),
                    condition_map=entry.get("condition_map", {}),
                    metadata=entry.get("metadata", {}),
                )
            )

        return cls(
            gene_sets=gene_sets,
            tf_list=tf_list,
            datasets=datasets,
            individual_genes=individual_genes,
            comparison_metrics=comparison_metrics,
        )

    # ── Loading ──────────────────────────────────────────────────────

    def _load(self, ds: DatasetEntry):
        """Load h5ad (lazy, cached per DatasetEntry)."""
        if ds._adata is not None:
            return ds._adata

        import scanpy as sc

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata = sc.read_h5ad(ds.h5ad_path)
        ds._adata = adata
        return adata

    # ── Gene expression matrix ───────────────────────────────────────

    def _get_expression_matrix(
        self, ds: DatasetEntry, cache_key: str, cell_mask: np.ndarray,
    ) -> pd.DataFrame:
        """Return cell-type x gene mean expression for cells matching mask."""
        if cache_key in self._gene_expr_cache:
            return self._gene_expr_cache[cache_key]

        adata = self._load(ds)

        use_raw = adata.raw is not None
        if use_raw:
            gene_list = adata.raw.var_names
            mat = adata.raw.X
        else:
            gene_list = adata.var_names
            mat = adata.X

        from scipy.sparse import issparse

        if issparse(mat):
            mat = mat.toarray()

        available = [g for g in self._all_genes if g in gene_list]
        gene_to_idx = {g: i for i, g in enumerate(gene_list)}
        idx_list = [gene_to_idx[g] for g in available]

        # Subset to cells of interest
        mat = mat[cell_mask]
        obs = adata.obs.iloc[cell_mask]

        rows = []
        for ct in sorted(obs[ds.cell_type_col].dropna().unique()):
            ct_mask = (obs[ds.cell_type_col] == ct).values
            ct_expr = mat[ct_mask][:, idx_list]
            row = {"cell_type": ct}
            for i, g in enumerate(available):
                row[g] = float(np.mean(ct_expr[:, i]))
            rows.append(row)

        df = pd.DataFrame(rows).set_index("cell_type")
        self._gene_expr_cache[cache_key] = df
        return df

    # ── GRN loading ──────────────────────────────────────────────────

    def _load_grn(self, ds: DatasetEntry) -> Optional[pd.DataFrame]:
        if ds.grn_table is None or not os.path.exists(ds.grn_table):
            return None
        df = pd.read_csv(ds.grn_table)
        col0 = df.columns[0]
        if col0 in ("Unnamed: 0", ""):
            df = df.rename(columns={col0: "cell_type"})
        elif "cell_type" not in df.columns:
            df = df.rename(columns={df.columns[0]: "cell_type"})
        df = df.set_index("cell_type")
        return df

    # ── Core analysis ────────────────────────────────────────────────

    def run(self) -> list[CellTypeResult]:
        """Run cross-paper analysis, producing per-cell-type results."""
        self._results = []

        for ds in self.datasets:
            adata = self._load(ds)
            split_by = ds.split_col or self.split_col

            if split_by and split_by in adata.obs.columns:
                condition_groups = sorted(adata.obs[split_by].unique().tolist())
                cond_map = ds.condition_map if ds.condition_map else {}
                for raw_cond in condition_groups:
                    cond_label = cond_map.get(raw_cond, str(raw_cond))
                    cond_mask = (adata.obs[split_by] == raw_cond).values
                    self._analyze_group(
                        ds=ds,
                        adata=adata,
                        mask=cond_mask,
                        cache_key=f"{ds.label}__{cond_label}",
                        condition=cond_label,
                    )
            else:
                self._analyze_group(
                    ds=ds,
                    adata=adata,
                    mask=np.ones(adata.n_obs, dtype=bool),
                    cache_key=ds.label,
                    condition=ds.condition_label or "",
                )

        return self._results

    def _analyze_group(
        self,
        ds: DatasetEntry,
        adata,
        mask: np.ndarray,
        cache_key: str,
        condition: str,
    ):
        total_cells = int(mask.sum())
        if total_cells == 0:
            return

        expr_df = self._get_expression_matrix(ds, cache_key, mask)
        grn_df = self._load_grn(ds)
        ct_col = ds.cell_type_col
        cell_types = sorted(adata.obs[ct_col].loc[mask].dropna().unique())

        all_genes_in_expr = set(expr_df.columns)

        def gene_expr(g: str) -> float:
            return float(expr_df.loc[ct, g]) if g in all_genes_in_expr else float("nan")

        for ct in cell_types:
            n_cells = int((adata.obs[ct_col] == ct).values.sum())
            # Use per-condition count
            ct_cond_mask = (adata.obs[ct_col] == ct).values & mask
            n_cells_cond = int(ct_cond_mask.sum())
            if n_cells_cond < 10:
                continue

            # Gene set scores
            gs_scores: dict[str, float] = {}
            for gs_key, gs_cfg in self.gene_sets.items():
                genes = [g for g in gs_cfg.get("genes", []) if g in all_genes_in_expr]
                gs_scores[gs_key] = float(expr_df.loc[ct, genes].mean()) if genes else float("nan")

            # Individual genes
            ind_expr: dict[str, float] = {}
            for g in self.individual_genes:
                ind_expr[g] = gene_expr(g)
            # Also include any genes from gene sets that are individually useful
            for gs_cfg in self.gene_sets.values():
                for g in gs_cfg.get("genes", []):
                    if g not in ind_expr:
                        ind_expr[g] = gene_expr(g)

            # Derived metrics
            derived: dict[str, float] = {}
            for m_name, m_cfg in self.comparison_metrics.items():
                formula = m_cfg.get("formula", "")
                try:
                    val = self._eval_formula(formula, ind_expr)
                    derived[m_name] = val
                except Exception:
                    derived[m_name] = float("nan")

            # GRN TF activity
            grn_vals: dict[str, float] = {}
            if grn_df is not None:
                for tf in self.tf_list:
                    if tf in grn_df.columns and ct in grn_df.index:
                        grn_vals[tf] = float(grn_df.loc[ct, tf])

            self._results.append(
                CellTypeResult(
                    dataset=ds.label,
                    paper=ds.paper,
                    cell_type=ct,
                    condition=condition,
                    n_cells=n_cells_cond,
                    cell_pct=round(n_cells_cond / total_cells * 100, 2),
                    gene_set_scores=gs_scores,
                    individual_expr=ind_expr,
                    derived_metrics=derived,
                    grn=grn_vals,
                )
            )

    @staticmethod
    def _eval_formula(formula: str, values: dict[str, float]) -> float:
        """Evaluate a derived metric formula safely.

        Formula supports: +, -, *, /, (, ) and variable names matching gene names.
        """
        # Build a safe evaluation dict
        safe: dict[str, float] = {}
        for k, v in values.items():
            if not np.isnan(v):
                safe[k] = v
        if not safe:
            return float("nan")

        try:
            result = eval(formula, {"__builtins__": {}}, safe)
            return float(result) if np.isfinite(result) else float("nan")
        except Exception:
            return float("nan")

    # ── Condition comparison ─────────────────────────────────────────

    def compare_conditions(
        self,
        metrics: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """Compute within-paper condition-to-condition fold changes.

        For each paper × cell_type, compute ratio of each metric
        between conditions (e.g. WIN / CTRL).
        """
        if not self._results:
            raise RuntimeError("Call run() first.")

        df = self.to_dataframe()
        df = df[df["condition"] != ""].copy()

        if metrics is None:
            # Auto-detect: gene set scores + individual genes + derived metrics
            metrics = self._available_metric_names(df)

        comparisons = []
        for (paper, dataset), group in df.groupby(["paper", "dataset"]):
            conditions = sorted(group["condition"].unique())
            for i, ca in enumerate(conditions):
                for cb in conditions[i + 1 :]:
                    ga = group[group["condition"] == ca].set_index("cell_type")
                    gb = group[group["condition"] == cb].set_index("cell_type")
                    common_ct = ga.index.intersection(gb.index)
                    for ct in common_ct:
                        for m in metrics:
                            va = ga.loc[ct, m]
                            vb = gb.loc[ct, m]
                            if pd.isna(va) or pd.isna(vb):
                                continue
                            fc = round(va / vb, 4) if vb != 0 else float("inf")
                            comparisons.append(
                                {
                                    "paper": paper,
                                    "dataset": dataset,
                                    "cell_type": ct,
                                    "metric": m,
                                    "condition_a": ca,
                                    "condition_b": cb,
                                    f"value_{ca}": round(va, 6),
                                    f"value_{cb}": round(vb, 6),
                                    "fold_change": fc,
                                    "log2_fc": round(np.log2(fc), 4)
                                    if fc > 0
                                    else float("-inf"),
                                }
                            )

        return pd.DataFrame(comparisons)

    @staticmethod
    def _available_metric_names(df: pd.DataFrame) -> list[str]:
        """Return metric column names from dataframe, excluding metadata cols."""
        meta = {"dataset", "paper", "condition", "cell_type", "n_cells", "cell_pct"}
        return [c for c in df.columns if c not in meta and not c.startswith("GRN_")]

    # ── Output ───────────────────────────────────────────────────────

    def to_dataframe(self) -> pd.DataFrame:
        records = []
        for r in self._results:
            rec = {
                "dataset": r.dataset,
                "paper": r.paper,
                "condition": r.condition,
                "cell_type": r.cell_type,
                "n_cells": r.n_cells,
                "cell_pct": r.cell_pct,
            }
            # Gene set scores
            for gs_key, score in r.gene_set_scores.items():
                gs_label = self.gene_sets.get(gs_key, {}).get("label", gs_key)
                rec[gs_label] = round(score, 4) if not np.isnan(score) else float("nan")
            # Individual genes
            for g_name, expr in r.individual_expr.items():
                rec[g_name] = round(expr, 4) if not np.isnan(expr) else float("nan")
            # Derived metrics
            for m_name, val in r.derived_metrics.items():
                rec[m_name] = round(val, 4) if not np.isnan(val) else float("nan")
            # GRN
            for tf in self.tf_list:
                val = r.grn.get(tf, float("nan"))
                rec[f"GRN_{tf}"] = round(val, 4) if not np.isnan(val) else float("nan")
            records.append(rec)
        return pd.DataFrame(records)

    def to_pivot(self, metric: str = "RA Synthesis") -> pd.DataFrame:
        df = self.to_dataframe()
        if metric not in df.columns:
            raise ValueError(
                f"Unknown metric: {metric}. Available: {list(df.columns)}"
            )
        has_condition = "condition" in df.columns and df["condition"].astype(bool).any()
        index_col = (
            ["cell_type", "condition"] if has_condition else "cell_type"
        )
        return df.pivot_table(
            index=index_col, columns="dataset", values=metric, aggfunc="mean"
        )

    def to_csv(self, output_path: str):
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        self.to_dataframe().to_csv(output_path, index=False)
        print(f"[CrossPaperAnalyzer] Comparison written: {output_path}")

    def to_json(self, output_path: str):
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        report = {
            "datasets": [ds.label for ds in self.datasets],
            "gene_sets": {
                k: v.get("genes", []) for k, v in self.gene_sets.items()
            },
            "results": [
                {
                    "dataset": r.dataset,
                    "paper": r.paper,
                    "condition": r.condition,
                    "cell_type": r.cell_type,
                    "n_cells": r.n_cells,
                    "cell_pct": r.cell_pct,
                    "gene_set_scores": {
                        k: v for k, v in r.gene_set_scores.items()
                    },
                    "individual_expr": r.individual_expr,
                    "derived_metrics": r.derived_metrics,
                    "grn": {k: v for k, v in r.grn.items()},
                }
                for r in self._results
            ],
        }
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2, default=str)
        print(f"[CrossPaperAnalyzer] JSON report written: {output_path}")

    def print_summary(self):
        df = self.to_dataframe()
        # Auto-select columns: metadata + gene set scores + individual_genes
        meta_cols = [
            "dataset",
            "condition",
            "paper",
            "cell_type",
            "n_cells",
            "cell_pct",
        ]
        metric_cols = [
            c
            for c in df.columns
            if c not in set(meta_cols) and not c.startswith("GRN_")
        ][:12]  # limit to 12 metric columns for readability
        cols = [c for c in meta_cols if c in df.columns] + metric_cols
        print(df[cols].to_string(index=False))
