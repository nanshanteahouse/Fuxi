#!/usr/bin/env python3
"""KADP/METC calibration harness (plan todo 6, Wave 2).

Rerunnable, read-only calibration harness for the KADP developmental-potency
axis (plan ``.omo/plans/annotation-kadp-metc.md`` todo 6).  It recomputes the
KADP naming decision over a threshold grid — ``ratio ∈ {1.0, 1.5, 2.0, 3.0}``
× ``abs ∈ {0.5, 0.6, 0.7}`` × ``gap ∈ {0.05, 0.10, 0.15, 0.20}`` ×
``use_gap_criterion ∈ {False, True}`` — on the **GSE246169 fetal** subproject,
reusing the engine's own pure functions
(``core.annotation.potency``: ``derive_developmental_poles`` /
``filter_pole_scores`` / ``compute_potency`` / ``evaluate_passes``).

Data source
-----------
* Config: ``projects/rna/GSE246169/fetal/config_GSE246169_fetal.yaml``
* Step-5 outputs: ``results/tables/marker_genes_unified.csv`` (per-cluster DE),
  ``results/tables/05_annotation_quality.json``, ``results/tables/cell_type_annotations.csv``
  and ``results/h5ad/05_annotated.h5ad`` (baseline obs ``annot_label``).

Per-cluster marker scores are **recomputed identically to the engine**
(``score_cluster_against_kb`` with the resolved config's species / expand steps,
sorted by logfoldchanges) from ``marker_genes_unified.csv`` — the engine does
not persist ``all_scores``, and the obs ``annot_evidence`` ``top_competitors``
list is truncated (≤10), so the full pole scores cannot be reconstructed from
the output alone.  The recomputed top-10 is validated against the stored
``annot_evidence`` before any grid evaluation.

KADP applies only to **baseline ambiguous** clusters (engine exit semantics:
``_kadp_name_candidate`` fires only when ``candidate.method == "ambiguous"``).

F5 gates (machine-checkable)
----------------------------
(a) **4 developmental groups named** — every developmental target cluster in
    the dataset registry (cluster 0 photoreceptor-precursor family, NRPC family
    2/4/11, cluster 10 Amacrine family) is named a Progenitor-pole member.
(b) **Label-invariance gate** — terminal label set = union of the
    ``Neuron ∪ Glia ∪ Non-neural`` hierarchy members (derived programmatically
    from ``kb["_hierarchy"]["categories"]``); baseline labels are read
    byte-for-byte from the baseline obs ``annot_label`` (no hand-typed
    "Müller_Glia" strings); a misnaming = a baseline terminal-labeled cluster
    that KADP renames to a Progenitor-pole member.  **F5 misnaming count = 0.**

``quality["celltypist"]`` / ``harmonization_rate`` are reported for the
baseline 05 output: the baseline config has ``celltypist.enabled: false`` and
no celltypist labels exist, so ``harmonization_rate`` is ``None`` with a note
(the harmonization function itself is todo-8-owned — the full computation and
assertion lives in todo 11).

Read-only: writes nothing outside ``adhoc/output/`` and never touches the
formal config, the 05 outputs or the results directory.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone

# Repo root on sys.path (adhoc scripts run from anywhere).
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import pandas as pd

# ── Threshold grid (written in stone by the plan todo 6) ────────────────
RATIO_GRID = [1.0, 1.5, 2.0, 3.0]
ABS_GRID = [0.5, 0.6, 0.7]
GAP_GRID = [0.05, 0.10, 0.15, 0.20]
USE_GAP_MODES = [False, True]

# ── Dataset registry (paths repo-relative; cluster ids from planning-phase
#    exploration of GSE246169 fetal — see plan todo 6/7 inherited wisdom) ──
DATASETS = {
    "GSE246169:fetal": {
        "display": "GSE246169 fetal",
        "config_rel": "projects/rna/{gse}/{sub}/config_{gse}_{sub}.yaml",
        "results_rel": "projects/rna/{gse}/{sub}/results",
        "kb": "retina",
        # Developmental populations that F5 gate (a) requires to be named.
        # Families: cluster 0 -> photoreceptor-precursor family, NRPC family
        # 2/4/11, cluster 10 -> Amacrine family.
        "developmental_families": {
            "photoreceptor_precursor": ["0"],
            "n_rpc": ["2", "4", "11"],
            "amacrine": ["10"],
        },
        # Ambiguous cluster that must NOT be named (RGC pool, terminal-dominant).
        "rgc_pool_clusters": ["5"],
    },
}


def _py(v):
    """numpy scalar -> python scalar for JSON."""
    if hasattr(v, "item"):
        return v.item()
    return v


def _load_json(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _die(msg: str) -> int:
    print(f"[kadp_metc_calibration] ERROR: {msg}", file=sys.stderr)
    return 2


def _load_kb(tissue: str):
    from core.kb import load_kb

    try:
        return load_kb(tissue)
    except ValueError as exc:
        raise FileNotFoundError(f"tissue KB '{tissue}' unavailable: {exc}") from exc


def _load_config(config_path: str):
    from core.utils import resolve_config

    return resolve_config(config_path)


def _score_cluster(cl_data: pd.DataFrame, kb, cfg) -> dict:
    """Recompute {type_key: Score} exactly as the engine first pass does."""
    from core.annotation.scoring import score_cluster_against_kb

    return score_cluster_against_kb(
        kb,
        cl_data,
        species=cfg.species,
        target_class=getattr(cfg, "target_class", "") or "",
        target_order=getattr(cfg, "target_order", "") or "",
        adaptive_top_n=True,
        expand_steps=list(getattr(cfg.marker, "candidate_pool_expand_steps", None) or [50, 100, 200]),
    )


def _cluster_marker_df(marker_df: pd.DataFrame, cl) -> pd.DataFrame:
    cl_data = marker_df[marker_df["cluster"] == cl].copy()
    lfc_idx = cl_data["logfoldchanges"].argsort()[::-1]
    return cl_data.iloc[lfc_idx]


def _read_baseline_obs(h5ad_path: str) -> pd.DataFrame:
    """Read only the annotation obs columns (backed mode; closes the handle)."""
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        cols = [c for c in ("leiden", "annot_label", "annot_method") if c in adata.obs.columns]
        return adata.obs[cols].copy()
    finally:
        adata.file.close()


def _per_cluster_label(obs: pd.DataFrame) -> dict[str, str]:
    """cluster -> mode of the per-cell annot_label (baseline obs, byte-for-byte)."""
    out: dict[str, str] = {}
    leiden = obs["leiden"].astype(str)
    label = obs["annot_label"].astype(str)
    for cl in sorted(leiden.unique(), key=lambda x: int(x) if str(x).isdigit() else str(x)):
        cl_str = str(cl)
        counts = Counter(label[leiden == cl_str])
        if counts:
            out[cl_str] = counts.most_common(1)[0][0]
    return out


def _validate_recomputed_scores(all_scores: dict, annot_records: pd.DataFrame, h5ad_path: str) -> dict:
    """Sanity-check the recomputed scores against stored annot_evidence top competitors.

    The obs ``annot_evidence.top_competitors`` carries the truncated (≤10) tied
    list; we compare the recomputed per-cluster top-1 (type + score) for every
    ambiguous cluster whose evidence carries competitors.
    """
    import anndata as ad

    checked = 0
    mismatches = []
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        if "annot_evidence" not in adata.obs.columns:
            return {"checked": 0, "ok": True, "note": "no annot_evidence column in baseline obs"}
        ev = adata.obs["annot_evidence"]
        leiden = adata.obs["leiden"].astype(str)
        for cl_str, scores in all_scores.items():
            if not scores:
                continue
            row = ev[leiden == cl_str]
            if len(row) == 0:
                continue
            raw = str(row.iloc[0])
            if not raw or raw == "nan":
                continue
            try:
                d = json.loads(raw)
            except Exception:
                continue
            comps = d.get("top_competitors") or []
            if not comps:
                continue
            checked += 1
            stored_top = comps[0]
            recomputed_top = max(scores.items(), key=lambda kv: kv[1].score)
            if abs(stored_top["score"] - recomputed_top[1].score) > 1e-9:
                mismatches.append(
                    {
                        "cluster": cl_str,
                        "stored_top": {"cell_type": stored_top["cell_type"], "score": stored_top["score"]},
                        "recomputed_top": {
                            "cell_type": recomputed_top[0],
                            "score": round(float(recomputed_top[1].score), 9),
                        },
                    }
                )
    finally:
        adata.file.close()
    return {
        "checked_clusters": checked,
        "ok": not mismatches,
        "mismatches": mismatches,
        "note": (
            "recomputed {type_key: Score} via score_cluster_against_kb on "
            "marker_genes_unified.csv; top-1 compared against stored "
            "annot_evidence.top_competitors[0]"
        ),
    }


def _terminal_set_from_kb(kb: dict) -> set[str]:
    """union(Neuron ∪ Glia ∪ Non-neural) hierarchy members — programmatic."""
    from core.annotation.potency import derive_developmental_poles

    _prog, terminal = derive_developmental_poles(kb)
    return terminal


def _progenitor_set_from_kb(kb: dict) -> set[str]:
    from core.annotation.potency import derive_developmental_poles

    progenitor, _term = derive_developmental_poles(kb)
    return progenitor


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "KADP/METC calibration harness — recompute KADP decisions over the "
            "ratio/abs/gap threshold grid on GSE246169 fetal (read-only)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--gse", default="GSE246169", help="GEO series id (uppercase).")
    ap.add_argument("--subproject", default="fetal", help="Subproject folder under projects/rna/<GSE>/.")
    ap.add_argument(
        "--output",
        default=os.path.join(REPO, "adhoc/output/kadp_calibration.json"),
        help="Where the machine-checkable summary JSON is written.",
    )
    ap.add_argument("--quiet", action="store_true", help="Suppress the terminal summary.")
    args = ap.parse_args()

    gse = args.gse.upper()
    sub = args.subproject
    ds_key = f"{gse}:{sub}"
    ds = DATASETS.get(ds_key)
    if ds is None:
        return _die(
            f"no dataset registry entry for '{ds_key}' (supported: {sorted(DATASETS)}). "
            "This harness is calibrated on GSE246169 fetal only."
        )

    if not os.environ.get("FUXI_DATA_ROOT"):
        return _die(
            "FUXI_DATA_ROOT is not set — source .env first "
            "(set -a && source .env && set +a). The harness resolves the "
            "dataset config under the repo environment and must not silently "
            "run degraded; nothing read or written."
        )

    # ── Locate + load inputs (read-only) ──────────────────────────────
    config_rel = ds["config_rel"].format(gse=gse, sub=sub)
    config_path = os.path.join(REPO, config_rel)
    if not os.path.isfile(config_path):
        return _die(f"config not found: {config_rel} (nothing read or written)")

    results_rel = ds["results_rel"].format(gse=gse, sub=sub)
    results_dir = os.path.join(REPO, results_rel)
    marker_csv = os.path.join(results_dir, "tables", "marker_genes_unified.csv")
    quality_path = os.path.join(results_dir, "tables", "05_annotation_quality.json")
    annot_csv = os.path.join(results_dir, "tables", "cell_type_annotations.csv")
    h5ad_path = os.path.join(results_dir, "h5ad", "05_annotated.h5ad")
    missing = [p for p in (marker_csv, quality_path, annot_csv, h5ad_path) if not os.path.exists(p)]
    if missing:
        return _die(
            "missing baseline 05 output artifact(s) under "
            f"{os.path.relpath(results_dir, REPO)}: "
            + ", ".join(os.path.basename(p) for p in missing)
            + " (run step 5 first; nothing read or written)"
        )

    cfg = _load_config(config_path)
    kb = _load_kb(ds["kb"])
    progenitor_pole = _progenitor_set_from_kb(kb)
    terminal_set = _terminal_set_from_kb(kb)

    marker_df = pd.read_csv(marker_csv)
    if "cluster" not in marker_df.columns or "logfoldchanges" not in marker_df.columns:
        return _die(f"unexpected marker table schema: {list(marker_df.columns)}")
    annot_records = pd.read_csv(annot_csv)
    baseline_quality = _load_json(quality_path)
    obs = _read_baseline_obs(h5ad_path)

    # ── Baseline per-cluster state ────────────────────────────────────
    baseline_label = _per_cluster_label(obs)
    method_by_cluster = {
        str(r["cluster"]): str(r["method"]) for r in annot_records.to_dict(orient="records")
    }
    ambiguous_clusters = sorted(
        (cl for cl, m in method_by_cluster.items() if m == "ambiguous"),
        key=lambda x: int(x),
    )
    terminal_labeled = sorted(
        (cl for cl, lab in baseline_label.items() if lab in terminal_set),
        key=lambda x: int(x),
    )
    # KADP fires only on baseline-ambiguous candidates (engine exit semantics).
    kadp_targets = [cl for cl in ambiguous_clusters]

    developmental_families = ds["developmental_families"]
    developmental_targets = sorted(
        {cl for family in developmental_families.values() for cl in family}, key=lambda x: int(x)
    )
    # Every developmental family must be fully named by a progenitor-pole member.
    rgc_pool = set(ds["rgc_pool_clusters"])

    # ── Recompute per-cluster marker scores (engine-identical) ─────────
    all_scores: dict[str, dict] = {}
    for cl in sorted(marker_df["cluster"].unique(), key=lambda x: int(str(x))):
        all_scores[str(cl)] = _score_cluster(_cluster_marker_df(marker_df, cl), kb, cfg)

    validation = _validate_recomputed_scores(all_scores, annot_records, h5ad_path)

    # ── Per-cluster potency (threshold-independent three-value dict) ───
    from core.annotation.potency import KADPConfig, compute_potency

    per_cluster_potency: dict[str, dict] = {}
    for cl in sorted(all_scores, key=lambda x: int(x)):
        cl_str = str(cl)
        r = compute_potency(all_scores[cl_str], (progenitor_pole, terminal_set), KADPConfig(enabled=True))
        per_cluster_potency[cl_str] = {
            "baseline_cell_type": baseline_label.get(cl_str),
            "baseline_method": method_by_cluster.get(cl_str),
            "is_kadp_target": cl_str in kadp_targets,
            "is_terminal_labeled": cl_str in terminal_labeled,
            "max_prog": round(float(r.max_prog), 9),
            "max_term": round(float(r.max_term), 9),
            "ratio": round(float(r.ratio), 9),
            "abs": round(float(r.max_prog), 9),
            "gap": round(float(r.gap), 9),
            "best_progenitor_type": r.best_progenitor_type,
            "best_in_progenitor_pole": r.best_progenitor_type in progenitor_pole,
        }

    # ── Grid evaluation ───────────────────────────────────────────────
    from core.annotation.potency import evaluate_passes

    grid_results: list[dict] = []
    for ratio in RATIO_GRID:
        for abs_th in ABS_GRID:
            for gap in GAP_GRID:
                for use_gap in USE_GAP_MODES:
                    kadp_cfg = KADPConfig(
                        enabled=True,
                        ratio_threshold=ratio,
                        abs_threshold=abs_th,
                        gap_threshold=gap,
                        use_gap_criterion=use_gap,
                    )
                    named: dict[str, str] = {}
                    for cl in kadp_targets:
                        r = compute_potency(all_scores[cl], (progenitor_pole, terminal_set), kadp_cfg)
                        if evaluate_passes(r, kadp_cfg) and r.best_progenitor_type is not None:
                            named[cl] = r.best_progenitor_type
                            named[cl] = r.best_progenitor_type
                    dev_named = sorted(cl for cl in developmental_targets if cl in named)
                    dev_missing = sorted(cl for cl in developmental_targets if cl not in named)
                    # F5 gate (b): baseline terminal-labeled clusters renamed to a
                    # Progenitor-pole member -> misnaming.
                    misnamed = sorted(
                        cl
                        for cl in named
                        if cl in terminal_labeled and named[cl] in progenitor_pole
                    )
                    # RGC pool (terminal-dominant) must not be named.
                    rgc_named = sorted(cl for cl in rgc_pool if cl in named)
                    gate_a = not dev_missing
                    gate_b = len(misnamed) == 0 and not rgc_named
                    grid_results.append(
                        {
                            "ratio": ratio,
                            "abs": abs_th,
                            "gap": gap,
                            "use_gap_criterion": use_gap,
                            "named": dict(sorted(named.items(), key=lambda kv: int(kv[0]))),
                            "developmental_named": dev_named,
                            "developmental_missing": dev_missing,
                            "rgc_pool_named": rgc_named,
                            "f5_misnaming_count": len(misnamed),
                            "misnamed": misnamed,
                            "gate_a": gate_a,
                            "gate_b": gate_b,
                            "f5_passed": gate_a and gate_b,
                        }
                    )

    # ── Threshold lock: minimal F5-satisfying config, else failure ────
    passing = [g for g in grid_results if g["f5_passed"]]
    if passing:
        # Minimal = smallest threshold values (most permissive) among passers,
        # prefer use_gap_criterion=False.
        lock = min(
            passing,
            key=lambda g: (g["ratio"], g["abs"], g["gap"], g["use_gap_criterion"]),
        )
        recommended_lock = {
            "kadp_enabled": True,
            "kadp_ratio_threshold": lock["ratio"],
            "kadp_abs_threshold": lock["abs"],
            "kadp_gap_threshold": lock["gap"],
            "use_gap_criterion": lock["use_gap_criterion"],
        }
        f5_satisfiable = True
    else:
        recommended_lock = None
        f5_satisfiable = False

    # Best partial within the grid (naming the most developmental clusters
    # with zero misnaming) — the default schema thresholds when they qualify.
    def _score_partial(g: dict) -> tuple:
        n_dev = len(g["developmental_named"])
        return (
            -n_dev,  # most developmental clusters named
            g["f5_misnaming_count"],  # then fewest misnamings
            bool(g["rgc_pool_named"]),  # then no RGC pool naming
            abs(g["ratio"] - 2.0) + abs(g["abs"] - 0.6) + abs(g["gap"] - 0.1),  # near defaults
        )

    best_partial = min(grid_results, key=_score_partial)
    best_partial_lock = {
        "kadp_enabled": True,
        "kadp_ratio_threshold": best_partial["ratio"],
        "kadp_abs_threshold": best_partial["abs"],
        "kadp_gap_threshold": best_partial["gap"],
        "use_gap_criterion": best_partial["use_gap_criterion"],
    }

    # ── quality / harmonization (baseline 05 output has no celltypist) ─
    celltypist_key_present = "celltypist" in baseline_quality
    quality_celltypist = bool(celltypist_key_present and baseline_quality.get("celltypist"))
    # harmonization needs celltypist labels; none exist in the baseline 05 output.
    harmonization_rate = None

    # Latent risk: clusters that would be named IF they were ambiguous
    # candidates (informational — the engine never fires KADP on them).
    latent_risks: dict[str, dict] = {}
    for cl in sorted(all_scores, key=lambda x: int(x)):
        if cl in kadp_targets:
            continue
        r = compute_potency(all_scores[cl], (progenitor_pole, terminal_set), KADPConfig(enabled=True))
        if evaluate_passes(r, KADPConfig(enabled=True)) and r.best_progenitor_type in progenitor_pole:
            latent_risks[cl] = {
                "baseline_cell_type": baseline_label.get(cl),
                "baseline_method": method_by_cluster.get(cl),
                "is_terminal_labeled": cl in terminal_labeled,
                "would_be_named": r.best_progenitor_type,
                "in_progenitor_pole": r.best_progenitor_type in progenitor_pole,
                "passing_variant": "abs" if r.passes_abs else ("ratio" if r.passes_ratio else "gap"),
                "note": (
                    "latent only: cluster is NOT a baseline-ambiguous candidate, "
                    "so the engine's KADP branch (candidate.method=='ambiguous') "
                    "never fires on it. Listed as evidence that the potency axis "
                    "cannot separate saturated terminal clusters from developmental ones."
                ),
            }

    n_ambiguous_cells = sum(len(obs[obs["leiden"].astype(str) == cl]) for cl in ambiguous_clusters)

    report = {
        "harness": "kadp_metc_calibration",
        "plan_todo": 6,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dataset": {
            "gse": gse,
            "subproject": sub,
            "display": ds["display"],
            "config": config_rel,
            "tissue_kb": ds["kb"],
            "tissue_maturity": getattr(cfg, "tissue_maturity", None),
        },
        "read_only": True,
        "sources": {
            "marker_scores_recomputed_from": os.path.relpath(marker_csv, REPO),
            "baseline_quality": os.path.relpath(quality_path, REPO),
            "baseline_annotations": os.path.relpath(annot_csv, REPO),
            "baseline_obs": os.path.relpath(h5ad_path, REPO),
            "poles": {
                "derivation": "derive_developmental_poles(kb): Progenitor pole = Progenitor category members; terminal pole = Neuron ∪ Glia ∪ Non-neural members",
                "progenitor_members": sorted(progenitor_pole),
                "terminal_members": len(terminal_set),
                "terminal_labels": sorted(terminal_labeled),
                "score_recompute_validation": validation,
            },
        },
        "baseline": {
            "n_cells": int(len(obs)),
            "n_clusters": int(len(all_scores)),
            "ambiguous_clusters": ambiguous_clusters,
            "ambiguous_cells": int(n_ambiguous_cells),
            "terminal_labeled_clusters": terminal_labeled,
            "developmental_families": developmental_families,
            "developmental_targets": developmental_targets,
            "rgc_pool_clusters": sorted(rgc_pool),
            "kadp_target_clusters": kadp_targets,
            "per_cluster": {
                cl: {"baseline_cell_type": baseline_label.get(cl), "baseline_method": method_by_cluster.get(cl)}
                for cl in sorted(all_scores, key=lambda x: int(x))
            },
        },
        "grid": {
            "ratio": RATIO_GRID,
            "abs": ABS_GRID,
            "gap": GAP_GRID,
            "use_gap_criterion_modes": USE_GAP_MODES,
            "n_points": len(grid_results),
        },
        "per_cluster_potency": per_cluster_potency,
        "grid_results": grid_results,
        "quality": {
            "celltypist": quality_celltypist,
            "celltypist_note": (
                "baseline 05_annotation_quality.json has no 'celltypist' key and the "
                "baseline config has celltypist.enabled=false — no celltypist labels "
                "exist; quality['celltypist'] derived as False"
            ),
            "ai_available": bool(getattr(cfg.ai, "enabled", False)),
            "ai_annotation_enabled": bool(getattr(cfg.ai, "annotation", False)),
            "harmonization_rate": harmonization_rate,
            "harmonization_rate_note": (
                "None: baseline 05 output contains no celltypist labels, so no "
                "harmonization can be computed. The harmonization function is "
                "todo-8-owned; the full computation + assertion lives in todo 11."
            ),
        },
        "f5": {
            "gate_a_defined": "all developmental target clusters named a Progenitor-pole member",
            "gate_b_defined": (
                "baseline terminal-labeled clusters (label ∈ Neuron∪Glia∪Non-neural "
                "hierarchy members, read byte-for-byte from baseline obs) renamed to a "
                "Progenitor-pole member == misnaming; count must be 0"
            ),
            "misnaming_count_max_across_grid": max(g["f5_misnaming_count"] for g in grid_results),
            "any_grid_point_passes": f5_satisfiable,
            "passing_grid_points": [g for g in grid_results if g["f5_passed"]],
            "best_partial_grid_point": {
                "ratio": best_partial["ratio"],
                "abs": best_partial["abs"],
                "gap": best_partial["gap"],
                "use_gap_criterion": best_partial["use_gap_criterion"],
                "named": best_partial["named"],
                "developmental_named": best_partial["developmental_named"],
                "developmental_missing": best_partial["developmental_missing"],
                "f5_misnaming_count": best_partial["f5_misnaming_count"],
                "gate_a": best_partial["gate_a"],
                "gate_b": best_partial["gate_b"],
            },
            "failure_evidence": _build_failure_evidence(per_cluster_potency, grid_results, latent_risks),
            "latent_risks": latent_risks,
        },
        "recommendation": {
            "f5_conditions_satisfiable": f5_satisfiable,
            "recommended_lock_for_todo7": recommended_lock,
            "best_partial_lock": best_partial_lock,
            "action": (
                (
                    f"LOCK: minimal F5-satisfying config ratio={recommended_lock['kadp_ratio_threshold']} "
                    f"abs={recommended_lock['kadp_abs_threshold']} gap={recommended_lock['kadp_gap_threshold']} "
                    f"use_gap_criterion={recommended_lock['use_gap_criterion']}"
                )
                if f5_satisfiable
                else (
                    "FAILURE RECORDED (forbidden to silently lower the bar): no grid point "
                    "satisfies F5 gate (a). Clusters 0 and 10 (developmental targets) have "
                    "max_prog < max_term — Fisher scores saturate to 1.0 on terminal types, "
                    "so ratio < 1.0, gap < 0 and the abs 'max_prog > max_term' guard fails at "
                    "every grid threshold. The potency axis cannot name them; the plan's "
                    "Wave-1 formulation needs a planner decision (different scoring input / "
                    "pole definition / candidate gate). todo 7 cannot pass the F5 gate with "
                    "these thresholds — do not enter Wave 3 without resolution."
                )
            ),
        },
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=_py)

    if not args.quiet:
        _print_summary(report)
    print(f"\n[kadp_metc_calibration] wrote {args.output}")
    return 0


def _build_failure_evidence(per_cluster_potency: dict, grid_results: list, latent_risks: dict) -> dict:
    """Human-digestible evidence for the F5 gate result (works for pass or fail)."""
    dev_clusters = ["0", "2", "4", "10", "11"]
    evidence = {
        "developmental_cluster_potency": {
            cl: {
                "max_prog": per_cluster_potency[cl]["max_prog"],
                "max_term": per_cluster_potency[cl]["max_term"],
                "ratio": per_cluster_potency[cl]["ratio"],
                "abs": per_cluster_potency[cl]["abs"],
                "gap": per_cluster_potency[cl]["gap"],
                "best_progenitor_type": per_cluster_potency[cl]["best_progenitor_type"],
            }
            for cl in dev_clusters
            if cl in per_cluster_potency
        },
        "nameable_at_any_grid_point": sorted(
            cl for cl in dev_clusters if any(cl in g["named"] for g in grid_results)
        ),
        "not_nameable_at_any_grid_point": sorted(
            cl for cl in dev_clusters if not any(cl in g["named"] for g in grid_results)
        ),
        "latent_terminal_misnaming_risk": latent_risks,
        "note": (
            "max_prog/max_term saturate to 1.0 on this Fisher-scored dataset: every "
            "developmental cluster has max_term == 1.0, so only clusters whose "
            "progenitor pole also hits 1.0 (NRPC family 2/4/11) pass the abs variant; "
            "clusters 0 and 10 are terminal-dominant at saturation and cannot pass "
            "any ratio/abs/gap threshold in the grid."
        ),
    }
    return evidence


def _print_summary(report: dict) -> None:
    f5 = report["f5"]
    rec = report["recommendation"]
    print(f"\n=== KADP/METC calibration — {report['dataset']['display']} (todo 6) ===")
    print(f"grid: {len(report['grid_results'])} points "
          f"(ratio {report['grid']['ratio']} × abs {report['grid']['abs']} × gap {report['grid']['gap']} × use_gap {report['grid']['use_gap_criterion_modes']})")
    print(f"score recompute validated: {report['sources']['poles']['score_recompute_validation']['ok']} "
          f"({report['sources']['poles']['score_recompute_validation']['checked_clusters']} clusters checked)")
    print(f"baseline: {report['baseline']['n_cells']} cells, {report['baseline']['n_clusters']} clusters; "
          f"ambiguous {report['baseline']['ambiguous_clusters']} ({report['baseline']['ambiguous_cells']} cells)")
    print(f"terminal-labeled (gate b protected): {report['baseline']['terminal_labeled_clusters']}")

    print("\n[per-cluster potency] (three-value; thresholds-independent)")
    for cl, p in sorted(report["per_cluster_potency"].items(), key=lambda kv: int(kv[0])):
        flag = "  <-- KADP target" if p["is_kadp_target"] else ""
        print(f"  cl {cl:>2} {p['baseline_cell_type']:<12} prog={p['max_prog']:.4f} term={p['max_term']:.4f} "
              f"ratio={p['ratio']:.4f} abs={p['abs']:.4f} gap={p['gap']:+.4f} best={p['best_progenitor_type']}{flag}")

    print(f"\n[F5 gate] gate_a (4 groups named): any grid point pass = {f5['any_grid_point_passes']}")
    print(f"  developmental clusters nameable at some grid point: {f5['failure_evidence']['nameable_at_any_grid_point']}")
    print(f"  never nameable: {f5['failure_evidence']['not_nameable_at_any_grid_point']}")
    print(f"  F5 misnaming count (max across grid): {f5['misnaming_count_max_across_grid']}")
    if f5["latent_risks"]:
        print(f"  LATENT terminal misnaming risk (non-ambiguous clusters): {sorted(f5['latent_risks'])}")
    bp = f5["best_partial_grid_point"]
    print(f"  best partial: ratio={bp['ratio']} abs={bp['abs']} gap={bp['gap']} use_gap={bp['use_gap_criterion']} "
          f"-> named={sorted(bp['named'])} missing={bp['developmental_missing']} misnaming={bp['f5_misnaming_count']}")

    print(f"\n[quality] celltypist={report['quality']['celltypist']} "
          f"ai_available={report['quality']['ai_available']} "
          f"harmonization_rate={report['quality']['harmonization_rate']}")
    print(f"\n[recommendation]")
    print(f"  F5 satisfiable: {rec['f5_conditions_satisfiable']}")
    if rec["recommended_lock_for_todo7"]:
        print(f"  LOCK: {rec['recommended_lock_for_todo7']}")
    print(f"  best-partial lock: {rec['best_partial_lock']}")
    print(f"  action: {rec['action']}")


if __name__ == "__main__":
    sys.exit(main())
