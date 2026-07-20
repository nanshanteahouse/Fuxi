#!/usr/bin/env python3
"""
core/cross_dataset_meta.py — Cross-dataset meta-clustering analysis engine.

Discovers completed pipeline output directories under ``projects/{modality}/*/results/tables/``,
loads cluster-level marker genes, and performs unsupervised meta-clustering (IDF-weighted
cosine similarity + community detection) to identify groups of clusters that appear
consistently across multiple datasets.

Key outputs
-----------
1. Cross-dataset communities — clusters from different datasets that share similar marker
   gene profiles, independent of their KB-sourced cell type labels.
2. Novelty classification — each community is flagged as ``NOVEL_CANDIDATE`` (contains
   unknown clusters from ≥3 datasets), ``KNOWN`` (consistent with KB), or ``LOW_PURITY``.
3. KB consistency audit — for each KB cell type, measures average cosine similarity among
   all clusters that received that label across datasets.

CLI usage
---------
::

    # Run with default discovery under projects/rna/
    python core/cross_dataset_meta.py

    # Specify a different modality or project root
    python core/cross_dataset_meta.py --modality atac
    python core/cross_dataset_meta.py --project-dir /custom/path --modality rna

    # Adjust analysis parameters
    python core/cross_dataset_meta.py --top-markers 20 --min-datasets 4 --cos-threshold 0.30
"""

from __future__ import annotations

import argparse
import sys
import os
from pathlib import Path
from collections import Counter, defaultdict
from typing import Iterable, Optional

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import cosine_similarity
import networkx as nx

# ── Repo root discovery ─────────────────────────────────────────────
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO))

from core.kb import load_kb


# ── Defaults ────────────────────────────────────────────────────────
DEFAULT_MODALITY = "rna"
DEFAULT_PROJECT_DIR = _REPO / "projects"
N_TOP_MARKERS = 15
MIN_DATASETS = 3
COSINE_THRESHOLD = 0.25


def _load_kb_markers() -> dict[str, set]:
    """Load the retina KB and return ``{cell_type: set(marker_genes)}``."""
    kb = load_kb("retina")
    kb_types = {
        k for k in kb
        if not k.startswith("_") and k != "expert_rules" and isinstance(kb[k], dict)
    }
    kb_markers: dict[str, set] = {}
    for ct in kb_types:
        m: set = set()
        for tier in ("confirm", "add"):
            m.update(kb[ct].get("markers", {}).get(tier, {}).keys())
        if m:
            kb_markers[ct] = m
    return kb_markers


# ── Dataset discovery ───────────────────────────────────────────────

def discover_datasets(
    project_dir: Path = DEFAULT_PROJECT_DIR,
    modality: str = DEFAULT_MODALITY,
) -> list[tuple[str, Path, str]]:
    """Auto-discover completed pipeline output directories.

    Scans ``{project_dir}/{modality}/*/results/`` for ``tables/`` or
    ``tables_*/`` directories containing ``cell_type_annotations.csv``.

    Returns a list of ``(dataset_id, tables_path, label)`` tuples.
    """
    datasets: list[tuple[str, Path, str]] = []
    modality_dir = project_dir / modality
    if not modality_dir.is_dir():
        print(f"  WARNING: {modality_dir} not found")
        return datasets

    for gse_dir in sorted(modality_dir.iterdir()):
        if not gse_dir.is_dir() or gse_dir.name.startswith("."):
            continue
        results_dir = gse_dir / "results"
        if not results_dir.is_dir():
            continue

        # Main tables dir
        main_tables = results_dir / "tables"
        if (main_tables / "cell_type_annotations.csv").exists():
            datasets.append(
                (gse_dir.name, main_tables, f"{gse_dir.name} {modality}")
            )

        # Subset tables dirs (e.g. tables_pcw8_multiome)
        for sub_dir in sorted(results_dir.iterdir()):
            if not sub_dir.is_dir() or sub_dir.name == "tables":
                continue
            if sub_dir.name.startswith("tables"):
                if (sub_dir / "cell_type_annotations.csv").exists():
                    subset_id = f"{gse_dir.name}_{sub_dir.name.replace('tables_', '')}"
                    datasets.append(
                        (subset_id, sub_dir, f"{gse_dir.name}/{sub_dir.name}")
                    )
    return datasets


# ── Cluster loading ─────────────────────────────────────────────────

def load_clusters(
    datasets: list[tuple[str, Path, str]],
    n_top: int = N_TOP_MARKERS,
) -> dict[str, dict]:
    """Load cluster-level marker data from discovered datasets.

    Returns a dict ``{cluster_name: {...}}`` where each value has keys:
    ``markers``, ``cell_type``, ``confidence``, ``diagnostic``,
    ``dataset``, ``cluster_id``, ``label``.
    """
    all_clusters: dict[str, dict] = {}

    for ds_id, tables_dir, label in datasets:
        ann_f = tables_dir / "cell_type_annotations.csv"
        m_f = tables_dir / "marker_genes_unified.csv"
        g_f = tables_dir / "marker_genes_per_group_cell_type_filtered.csv"

        if not ann_f.exists():
            continue

        ann = pd.read_csv(ann_f)
        mdf = pd.read_csv(m_f) if m_f.exists() else (
            pd.read_csv(g_f) if g_f.exists() else pd.DataFrame()
        )
        if mdf.empty:
            print(f"  SKIP {ds_id}: no marker file")
            continue

        has_logfc = "logfoldchanges" in mdf.columns
        has_grp = "group" in mdf.columns

        for _, row in ann.iterrows():
            cid = str(int(row.get("cluster", -1)))
            ct = str(row.get("cell_type", "unknown"))
            conf = str(row.get("confidence", "unknown"))
            diag = str(row.get("diagnostic_category", ""))

            if has_grp:
                cm = mdf[mdf["group"] == ct]
                if cm.empty:
                    cm = mdf[mdf["group"].astype(str) == cid]
            else:
                cm = mdf

            if cm.empty:
                continue

            if has_logfc:
                cm = cm.dropna(subset=["logfoldchanges"])
                top = cm.nlargest(n_top, "logfoldchanges")
            elif "scores" in cm.columns:
                top = cm.nlargest(n_top, "scores")
            else:
                top = cm.head(n_top)

            gene_col = "names" if "names" in top.columns else top.columns[0]
            markers = {str(g).upper().strip() for g in top[gene_col]}
            if len(markers) < 3:
                continue

            name = f"{ds_id}__{cid}"
            all_clusters[name] = {
                "markers": markers,
                "cell_type": ct,
                "confidence": conf,
                "diagnostic": diag,
                "dataset": ds_id,
                "cluster_id": cid,
                "label": label,
            }

    return all_clusters


# ── Vectorization ───────────────────────────────────────────────────

def compute_idf_vectors(
    all_clusters: dict[str, dict],
) -> tuple[np.ndarray, list[str], dict[str, float]]:
    """Build IDF-weighted TF vectors from cluster marker sets.

    Returns ``(tfidf_norm, cluster_names, idf)``.
    """
    names = sorted(all_clusters)
    n = len(names)

    all_genes: set = set()
    for nm in names:
        all_genes.update(all_clusters[nm]["markers"])
    genes = sorted(all_genes)
    g2i = {g: i for i, g in enumerate(genes)}
    m = len(genes)

    df = Counter()
    for nm in names:
        for g in all_clusters[nm]["markers"]:
            df[g] += 1

    idf = {g: np.log((n + 1.0) / (df[g] + 1.0)) for g in genes}

    data, rows, cols = [], [], []
    for i, nm in enumerate(names):
        for g in all_clusters[nm]["markers"]:
            data.append(idf[g])
            rows.append(i)
            cols.append(g2i[g])

    X = csr_matrix((data, (rows, cols)), shape=(n, m))
    X_norm = normalize(X, norm="l2", axis=1)
    return X_norm, names, idf


# ── Community detection ──────────────────────────────────────────────

def detect_communities(
    sim: np.ndarray,
    names: list[str],
    threshold: float = COSINE_THRESHOLD,
) -> list[set[int]]:
    """Build similarity graph and detect communities (Louvain)."""
    G = nx.Graph()
    for i in range(len(names)):
        G.add_node(i)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if sim[i, j] >= threshold:
                G.add_edge(i, j, weight=float(sim[i, j]))

    if G.number_of_edges() == 0:
        return []

    try:
        from networkx.algorithms.community import louvain_communities
        return list(louvain_communities(G, weight="weight", seed=42))
    except Exception:
        return list(nx.connected_components(G))


# ── KB matching ──────────────────────────────────────────────────────

def match_kb_markers(
    markers: set,
    kb_markers: dict[str, set],
    min_overlap: int = 3,
) -> tuple[Optional[str], int]:
    """Find the KB cell type with the most marker overlap."""
    best, best_count = None, 0
    for ct, kg in kb_markers.items():
        ov = len(markers & kg)
        if ov > best_count:
            best_count, best = ov, ct
    return (best, best_count) if best_count >= min_overlap else (None, best_count)


# ── Main analysis ────────────────────────────────────────────────────

def analyze(
    all_clusters: dict[str, dict],
    kb_markers: dict[str, set],
    min_datasets: int = MIN_DATASETS,
    cos_threshold: float = COSINE_THRESHOLD,
) -> dict:
    """Run the full meta-clustering analysis pipeline.

    Parameters
    ----------
    all_clusters : dict
        Cluster data from :func:`load_clusters`.
    kb_markers : dict
        KB marker sets from :func:`_load_kb_markers`.
    min_datasets : int
        Minimum datasets for ``NOVEL_CANDIDATE`` classification.
    cos_threshold : float
        Cosine similarity threshold for graph edges.

    Returns
    -------
    dict with keys:
        ``communities`` — list of community dicts
        ``kb_consistency`` — per-cell-type cross-dataset similarity
        ``summary`` — classification counts
        ``similarity_matrix`` — (n x n) array
        ``cluster_names`` — list matching matrix order
    """
    names = sorted(all_clusters)
    if len(names) < 2:
        return {
            "communities": [],
            "kb_consistency": {},
            "summary": {},
            "similarity_matrix": np.array([]),
            "cluster_names": names,
        }

    X_norm, names, idf = compute_idf_vectors(all_clusters)

    # Cosine similarity
    S = cosine_similarity(X_norm)

    # Communities
    comms = detect_communities(S, names, cos_threshold)

    # Characterize
    results: list[dict] = []
    for ci, comm in enumerate(comms):
        if len(comm) < 2:
            continue
        cn = [names[i] for i in comm]
        cc = [all_clusters[nn] for nn in cn]
        ds = set(c["dataset"] for c in cc)
        nd = len(ds)
        cts = [c["cell_type"] for c in cc]
        ctc = Counter(cts)
        dom = ctc.most_common(1)[0][0]
        pur = ctc[dom] / len(cts)
        unks = [
            c for c in cc
            if c["cell_type"].lower() in ("unknown", "uncertain")
            or c["diagnostic"] in ("no_kb_match", "true_unknown")
        ]
        nu = len(unks)
        am: set = set()
        for c in cc:
            am.update(c["markers"])
        kh, ko = match_kb_markers(am, kb_markers)

        if nu > 0 and nd >= min_datasets:
            nov = "NOVEL_CANDIDATE"
        elif pur >= 0.7 and dom.lower() != "unknown":
            nov = "KNOWN"
        elif pur >= 0.5:
            nov = "MIXED"
        else:
            nov = "LOW_PURITY"

        results.append({
            "id": f"C{ci:02d}",
            "size": len(comm),
            "n_datasets": nd,
            "datasets": sorted(ds),
            "dominant": dom,
            "purity": pur,
            "ct_distribution": dict(ctc),
            "n_unknown": nu,
            "kb_best": kh,
            "kb_overlap": ko,
            "novelty": nov,
            "top_markers": sorted(
                am, key=lambda g: idf.get(g, 0), reverse=True
            )[:15],
            "members": [
                (c["dataset"], c["cluster_id"], c["cell_type"], c["confidence"])
                for c in cc
            ],
        })
    results.sort(key=lambda r: (-r["n_datasets"], -r["n_unknown"], -r["purity"], -r["size"]))

    # KB consistency
    kb_consistency: dict = {}
    for nm in names:
        ct = all_clusters[nm]["cell_type"]
        if ct.lower() not in ("unknown", "uncertain"):
            kb_consistency.setdefault(ct, []).append(nm)
    for ct in list(kb_consistency):
        mems = kb_consistency[ct]
        if len(mems) < 3:
            del kb_consistency[ct]
            continue
        idxs = [names.index(m) for m in mems]
        sims = [
            S[idxs[i], idxs[j]]
            for i in range(len(idxs))
            for j in range(i + 1, len(idxs))
        ]
        ds_list = sorted(set(all_clusters[m]["dataset"] for m in mems))
        kb_consistency[ct] = {
            "n_clusters": len(mems),
            "n_datasets": len(ds_list),
            "avg_cosine": float(np.mean(sims)),
            "datasets": ds_list,
        }

    # Summary
    summary = Counter(r["novelty"] for r in results)

    return {
        "communities": results,
        "kb_consistency": kb_consistency,
        "summary": dict(summary),
        "similarity_matrix": S,
        "cluster_names": names,
    }


# ── CLI ──────────────────────────────────────────────────────────────

def print_report(result: dict) -> None:
    """Print a human-readable analysis report."""
    communities = result["communities"]
    kb_consistency = result["kb_consistency"]
    summary = result["summary"]

    print("\n" + "=" * 80)
    print(f"COMPLETED: {result['similarity_matrix'].shape[0]} clusters, "
          f"{len(communities)} communities")
    print("=" * 80)

    for r in communities:
        flag = {
            "NOVEL_CANDIDATE": "🆕",
            "KNOWN": "✅",
            "MIXED": "⚠️",
        }.get(r["novelty"], "❓")
        print(
            f"\n{flag} {r['id']} ds={r['n_datasets']} n={r['size']} "
            f"pur={r['purity']:.2f} unk={r['n_unknown']} {r['novelty']}"
        )
        print(f"   ds={r['datasets']}  ct={r['dominant']}  "
              f"kb={r['kb_best']}(ov={r['kb_overlap']})")
        print(f"   top: {', '.join(r['top_markers'][:10])}")
        ms = "; ".join(f"{d}/{c}={t}" for d, c, t, _ in r["members"][:6])
        print(f"   members: {ms}")

    # KB consistency
    print("\n=== KB consistency ===")
    for ct in sorted(kb_consistency):
        kc = kb_consistency[ct]
        flag = "⚠️" if kc["avg_cosine"] < 0.15 else "✅"
        print(
            f"  {flag} {ct}: {kc['n_clusters']} clusters × {kc['n_datasets']} ds, "
            f"avg_cos={kc['avg_cosine']:.3f} | {kc['datasets']}"
        )

    # Summary
    print("\n=== SUMMARY ===")
    for cat in ("KNOWN", "NOVEL_CANDIDATE", "MIXED", "LOW_PURITY"):
        count = summary.get(cat, 0)
        if count:
            print(f"  {cat}: {count}")

    novel = [r for r in communities if r["novelty"] == "NOVEL_CANDIDATE"]
    if novel:
        print("\n--- NOVEL CANDIDATES ---")
        for r in novel:
            print(
                f"\n  🆕 {r['id']}: {r['n_datasets']} ds, {r['n_unknown']} unknowns"
            )
            print(f"     KB: {r['kb_best']}(ov={r['kb_overlap']})  "
                  f"markers: {', '.join(r['top_markers'][:12])}")
            for d, c, t, conf in r["members"]:
                print(f"       {d}/c{c}: {t} ({conf})")
    else:
        print("\n  NOVEL_CANDIDATE: 0 — no cross-dataset unknown communities found")
        print("  → Try more datasets or enable unconstrained_annotation mode")


# ── CLI entry point ────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cross-dataset meta-clustering analysis engine"
    )
    parser.add_argument(
        "--project-dir",
        type=str,
        default=str(DEFAULT_PROJECT_DIR),
        help=f"Project root directory (default: {DEFAULT_PROJECT_DIR})",
    )
    parser.add_argument(
        "--modality",
        type=str,
        default=DEFAULT_MODALITY,
        help=f"Modality subdirectory (default: {DEFAULT_MODALITY})",
    )
    parser.add_argument(
        "--top-markers",
        type=int,
        default=N_TOP_MARKERS,
        help=f"Top N markers per cluster (default: {N_TOP_MARKERS})",
    )
    parser.add_argument(
        "--min-datasets",
        type=int,
        default=MIN_DATASETS,
        help=f"Min datasets for NOVEL classification (default: {MIN_DATASETS})",
    )
    parser.add_argument(
        "--cos-threshold",
        type=float,
        default=COSINE_THRESHOLD,
        help=f"Cosine similarity threshold (default: {COSINE_THRESHOLD})",
    )
    args = parser.parse_args()

    project_dir = Path(args.project_dir)
    modality = args.modality

    print(f"Project: {project_dir}")
    print(f"Modality: {modality}")
    print(f"Parameters: top_markers={args.top_markers}, "
          f"min_datasets={args.min_datasets}, cos_threshold={args.cos_threshold}")

    # Discover
    print("\n=== Discovering datasets ===")
    datasets = discover_datasets(project_dir, modality)
    if not datasets:
        print("  No completed datasets found.")
        return
    print(f"  Found {len(datasets)} dataset(s)")
    for ds_id, _, label in datasets:
        print(f"  - {ds_id}: {label}")

    # Load
    print("\n=== Loading clusters ===")
    kb_markers = _load_kb_markers()
    print(f"KB: {len(kb_markers)} cell types with markers")
    all_clusters = load_clusters(datasets, args.top_markers)
    print(f"Total: {len(all_clusters)} clusters")

    if not all_clusters:
        print("No clusters loaded. Check your data.")
        return

    # Analyze
    print("\n=== Analyzing ===")
    result = analyze(all_clusters, kb_markers, args.min_datasets, args.cos_threshold)

    # Report
    print_report(result)


if __name__ == "__main__":
    main()
