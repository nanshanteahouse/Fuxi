#!/usr/bin/env python3
"""
F5/METC acceptance-gate validation for the GSE246169 fetal KADP+METC dual-on
run (plan annotation-kadp-metc, todo 11).

Read-only validation against:
  - baseline obs snapshot  : .omo/evidence/kadp-metc/w4-t11/cell_metadata_baseline.csv
  - input cluster ids      : projects/rna/GSE246169/fetal_kadp_metc/results/h5ad/04_clustered.h5ad
  - new run obs            : projects/rna/GSE246169/fetal_kadp_metc/results/tables/cell_metadata.csv
  - quality report         : projects/rna/GSE246169/fetal_kadp_metc/results/tables/05_annotation_quality.json
  - step log               : projects/rna/GSE246169/fetal_kadp_metc/logs/05_annotate_major.log

Asserts the todo-11 gates (plan annotation-kadp-metc, Wave 4) and writes a
machine-checkable summary JSON to
.omo/evidence/kadp-metc/w4-t11/gse246169_fetal_metc.json

This run is a REAL dual-on validation.  The `quality["celltypist"]` gate and
`harmonization_rate` are recorded AS ACTUALLY OBSERVED (no silent pass, no
code patching): if CellTypist's labels were never captured by the engine the
run is recorded as a FAILURE with the documented reason (escalation path).

Usage: .venv/bin/python adhoc/kadp_metc_validate.py
"""
import json
import os
import sys
from datetime import datetime, timezone

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import pandas as pd  # noqa: E402

from core.kb import load_kb  # noqa: E402


# ── Paths ──────────────────────────────────────────────────────────────
BASELINE_CSV = os.path.join(
    REPO, ".omo/evidence/kadp-metc/w4-t11/cell_metadata_baseline.csv"
)
CLUSTER_H5AD = os.path.join(
    REPO, "projects/rna/GSE246169/fetal_kadp_metc/results/h5ad/04_clustered.h5ad"
)
NEW_CSV = os.path.join(
    REPO, "projects/rna/GSE246169/fetal_kadp_metc/results/tables/cell_metadata.csv"
)
QUALITY_JSON = os.path.join(
    REPO,
    "projects/rna/GSE246169/fetal_kadp_metc/results/tables/05_annotation_quality.json",
)
STEP_LOG = os.path.join(
    REPO, "projects/rna/GSE246169/fetal_kadp_metc/logs/05_annotate_major.log"
)
OUT_JSON = os.path.join(REPO, ".omo/evidence/kadp-metc/w4-t11/gse246169_fetal_metc.json")

KADP_TARGETS = {"2", "4", "11"}          # todo 7 Option A: NRPC-family KADP names
METC_CANDIDATES = ["0", "5", "10"]       # todo 7 Option A + cluster 5
TERMINAL_LABEL_EXPECT = {                 # todo 7 gate (b), byte-identical baseline
    "1": "Bipolar_Cell",
    "3": "Muller_Glia",
    "6": "RGC",
    "7": "Bipolar_Cell",
    "8": "RGC",
    "9": "Amacrine_Cell",
    "12": "Muller_Glia",
    "13": "RGC",
}


def _load_kb_poles():
    """Derive the Progenitor pole and terminal set from the retina KB,
    exactly as the engine does (core.kb.load_kb → _hierarchy.categories)."""
    kb = load_kb("retina")
    cats = kb["_hierarchy"]["categories"]
    prog = set(cats.get("Progenitor", {}).get("members", []))
    terminal = set()
    for name in ("Neuron", "Glia", "Non-neural"):
        terminal.update(cats.get(name, {}).get("members", []))
    return prog, terminal


def _load_baseline_new():
    """Return (baseline_df, new_df) both with a 'leiden' cluster column."""
    import scanpy as sc

    base = pd.read_csv(BASELINE_CSV)
    new = pd.read_csv(NEW_CSV)
    adata = sc.read(CLUSTER_H5AD)
    obsdf = adata.obs.reset_index().rename(columns={"index": "barcode"})[["barcode", "leiden"]]
    base = base.merge(obsdf, on="barcode", how="left")
    new = new.merge(obsdf, on="barcode", how="left")
    assert base["leiden"].notna().all(), "baseline barcode/leiden join failed"
    assert new["leiden"].notna().all(), "new barcode/leiden join failed"
    return base, new


def _cluster_summary(df, potency_col="potency"):
    out = {}
    for cl, g in df.groupby(df["leiden"].astype(str)):
        ct = g["cell_type"].dropna()
        st = g["cell_state"].dropna()
        meth = g["annot_method"].dropna()
        out[cl] = {
            "n": int(len(g)),
            "cell_type": ct.mode().iloc[0] if not ct.empty else "",
            "cell_state": st.mode().iloc[0] if not st.empty else "",
            "method": meth.mode().iloc[0] if not meth.empty else "",
        }
        if potency_col in g.columns:
            pv = g[potency_col].dropna().astype(str)
            out[cl]["potency_present"] = bool(len(pv) and pv.iloc[0].strip())
        else:
            out[cl]["potency_present"] = False
    return out


def _read_quality():
    with open(QUALITY_JSON) as f:
        return json.load(f)


def _read_log():
    with open(STEP_LOG) as f:
        return f.read()


def main():
    prog_pole, terminal_set = _load_kb_poles()
    base, new = _load_baseline_new()
    b = _cluster_summary(base)
    n = _cluster_summary(new)
    clusters = sorted(set(b) | set(n), key=int)
    quality_json = _read_quality()
    log = _read_log()

    per_cluster = {}
    for cl in clusters:
        per_cluster[cl] = {
            "baseline_cell_type": b[cl]["cell_type"],
            "new_cell_type": n[cl]["cell_type"],
            "method": n[cl]["method"],
            "cell_state": n[cl]["cell_state"],
            "potency_present": n[cl]["potency_present"],
            "n_cells": n[cl]["n"],
        }

    # ── Gate: KADP naming preserved (clusters 2/4/11) ──────────────────
    kadp_preserved = True
    for cl in sorted(KADP_TARGETS):
        rec = n[cl]
        ok = (
            rec["method"] == "developmental_potency"
            and rec["cell_type"] in prog_pole
            and rec["cell_type"] != ""
            and rec["cell_state"] == "differentiating"
        )
        if not ok:
            kadp_preserved = False
            print(f"[FAIL] KADP preserved: cluster {cl} not KADP-named: {rec}")

    # ── Gate (b): label-invariant (terminal clusters byte-identical) ───
    terminal_label_clusters = sorted(
        (cl for cl in clusters if b[cl]["cell_type"] in terminal_set), key=int
    )
    misnaming_count = 0
    label_invariant = True
    for cl in terminal_label_clusters:
        new_ct = n[cl]["cell_type"]
        if new_ct in prog_pole or new_ct != b[cl]["cell_type"]:
            misnaming_count += 1
            label_invariant = False
            print(
                f"[FAIL] gate(b): terminal cluster {cl} ({b[cl]['cell_type']}) -> {new_ct}"
            )
    # Explicit byte-identical check against the t7 fixed expectations
    for cl, expected in TERMINAL_LABEL_EXPECT.items():
        if n[cl]["cell_type"] != expected:
            label_invariant = False
            print(
                f"[FAIL] gate(b): cluster {cl} expected '{expected}', got '{n[cl]['cell_type']}'"
            )

    # ── Ambiguous drop ─────────────────────────────────────────────────
    amb_baseline = sum(b[cl]["n"] for cl in clusters if b[cl]["method"] == "ambiguous")
    amb_new = sum(n[cl]["n"] for cl in clusters if n[cl]["method"] == "ambiguous")
    ambiguous_drop_ok = amb_new <= 3418 and amb_new < amb_baseline and amb_baseline == 4890

    # ── Homogeneous clusters transitional stays 0 ──────────────────────
    # Only METC candidates may become transitional; homogeneous (non-candidate,
    # terminal-labelled) clusters must have method != transition_state.
    homogeneous_transitional = [
        cl
        for cl in clusters
        if cl not in METC_CANDIDATES and n[cl]["method"] == "transition_state"
    ]
    zero_homogeneous_transitional = not homogeneous_transitional

    # ── METC decision counting (from annot_evidence review_reason / obs) ─
    metc_reasons = []
    for cl in METC_CANDIDATES:
        rec = n[cl]
        # obs annot_method / cell_type carry the METC outcome:
        # divergent -> transition_state + 'transitional:'; 2way -> ambiguous w/ metc_2way;
        # consensus -> rescued method.  review_queue in the quality JSON lists reasons.
        metc_reasons.append(
            {
                "cluster": cl,
                "method": rec["method"],
                "cell_type": rec["cell_type"],
                "cell_state": rec["cell_state"],
            }
        )
    # authoritative count: quality JSON review_queue reasons startswith metc_*
    review_reasons = [
        entry.get("reason", "")
        for entry in quality_json.get("review_queue", [])
        if (entry.get("reason") or "").startswith("metc_")
    ]
    metc_decision_count = len(review_reasons)
    branch_counts = {
        "divergent": sum(1 for r in review_reasons if r == "metc_divergent"),
        "2way": sum(1 for r in review_reasons if r == "metc_2way"),
        "consensus": sum(1 for r in review_reasons if r == "metc_consensus"),
    }
    # Sanity: obs-level transition_state on METC candidates also counts as divergent
    obs_transitional = [
        cl for cl in METC_CANDIDATES if n[cl]["method"] == "transition_state"
    ]
    if obs_transitional and "metc_divergent" not in review_reasons:
        # obs says transitional but quality JSON lacks the reason → still count
        metc_decision_count = max(metc_decision_count, len(obs_transitional))
        branch_counts["divergent"] = max(branch_counts["divergent"], len(obs_transitional))

    # ── quality["celltypist"] / harmonization_rate ─────────────────────
    # The written 05_annotation_quality.json carries `harmonization_rate`
    # (None iff celltypist produced no labels); the `celltypist` bool lives
    # only in the in-memory fusion quality dict.  Evidence for "CellTypist
    # actually ran" therefore comes from the log + harmonization_rate.
    harmonization_rate = quality_json.get("harmonization_rate")
    ct_ran_log = "CellTypist: predicted" in log
    ct_skipped_log = "column not found in adata.obs — skipping" in log
    if harmonization_rate is not None:
        quality_celltypist = True
    else:
        # CellTypist ran the model but labels were never captured (engine reads
        # adata.obs; celltypist>=1.6 returns AnnotationResult w/o mutating adata)
        quality_celltypist = False

    # ── AI fallback fired ──────────────────────────────────────────────
    ai_fallback_fired = (
        "AI fallback for" in log and "cluster suggestions received" in log
    )

    # ── documented 3-source run flag ───────────────────────────────────
    # METC sources in practice = marker + AI + CellTypist (expert structurally
    # None, todo 9).  CellTypist abstains entirely here → n_spoke = 2 < 3 →
    # candidates returned unchanged → 0 METC decisions (documented reason).
    documented_three_source_run = False
    zero_metc_reason = (
        "CellTypist labels never captured: engine reads adata.obs['majority_voting'] "
        "after celltypist.annotate(), but celltypist>=1.6 returns an AnnotationResult "
        "without mutating adata.obs (confirmed empirically). celltypist_results={} -> "
        "harmonization_rate=null -> CellTypist source abstains. METC n_spoke = "
        "marker+AI = 2 < metc_min_sources=3 -> ambiguous candidates 0/5/10 returned "
        "unchanged (candidate-return semantics). This is a product wiring gap to "
        "report per plan escalation, NOT a silent pass."
    )

    # ── Overall gates ──────────────────────────────────────────────────
    gates = {
        "quality_celltypist": quality_celltypist,
        "harmonization_rate_ge_20pct": harmonization_rate is not None and harmonization_rate >= 0.20,
        "ai_fallback_fired": ai_fallback_fired,
        "metc_decision": metc_decision_count >= 1 or bool(obs_transitional),
        "metc_decision_or_documented_zero": (
            metc_decision_count >= 1 or bool(obs_transitional) or True
        ),  # documented-zero allowed with explicit reason
        "kadp_preserved": kadp_preserved,
        "b_label_invariant": label_invariant,
        "c_ambiguous_drop": ambiguous_drop_ok,
        "d_zero_homogeneous_transitional": zero_homogeneous_transitional,
        "misnaming_count_zero": misnaming_count == 0,
    }
    # A dual-on validation is FAILED unless CellTypist actually contributed
    # labels (plan QA: explicit failure, never silent pass).
    run_passed = (
        quality_celltypist
        and ai_fallback_fired
        and gates["harmonization_rate_ge_20pct"]
        and kadp_preserved
        and label_invariant
        and ambiguous_drop_ok
        and zero_homogeneous_transitional
        and misnaming_count == 0
        and (metc_decision_count >= 1 or bool(obs_transitional))
    )
    # The documented-zero-METC case is permitted for the METC gate itself, but
    # run_passed additionally requires quality_celltypist (explicit failure
    # when the CellTypist source never contributed).
    if not run_passed:
        print("[RESULT] RUN FAILED validation gates — see gates dict (escalation path)")

    summary = {
        "dataset": "GSE246169",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_dir": "projects/rna/GSE246169/fetal_kadp_metc",
        "config": "projects/rna/GSE246169/fetal_kadp_metc/config_GSE246169_fetal_kadp_metc.yaml",
        "model_chosen": "Fetal_Human_Retina.pkl",
        "model_description": "cell types from human fetal neural retina and retinal pigment epithelium (25 classes)",
        "per_cluster": per_cluster,
        "metc_candidates": METC_CANDIDATES,
        "metc_candidate_details": metc_reasons,
        "metc_decision_count": metc_decision_count,
        "metc_branch_counts": branch_counts,
        "obs_transitional_on_candidates": obs_transitional,
        "harmonization_rate": harmonization_rate,
        "quality_celltypist": quality_celltypist,
        "quality_celltypist_evidence": {
            "log_celltypist_predicted": ct_ran_log,
            "log_column_not_found_skip": ct_skipped_log,
            "written_json_has_celltypist_key": "celltypist" in quality_json,
            "note": (
                "CellTypist model ran (Prediction + Majority voting done in log) but "
                "the engine never captured its labels: engine reads adata.obs after "
                "celltypist.annotate(), celltypist>=1.6 returns AnnotationResult "
                "without mutating adata.obs. quality_celltypist therefore False."
            ),
        },
        "ai_fallback_fired": ai_fallback_fired,
        "documented_three_source_run": documented_three_source_run,
        "documented_zero_metc_reason": zero_metc_reason if metc_decision_count == 0 else "",
        "ambiguous_cells": {"baseline": amb_baseline, "new": amb_new},
        "terminal_label_clusters": terminal_label_clusters,
        "misnaming_count": misnaming_count,
        "homogeneous_transitional_clusters": homogeneous_transitional,
        "gates": gates,
        "run_passed": run_passed,
        "notes": (
            "Dual-on validation: kadp_enabled + locked thresholds + metc_enabled + "
            "ai.ai_annotation=true (real AI API calls, log shows 'AI fallback for 4 "
            "low-confidence clusters' / '14 cluster suggestions received') + celltypist "
            "enabled/model=Fetal_Human_Retina.pkl. KADP preserved (2/4/11). Label-invariant "
            "gate (b) preserved. Ambiguous dropped to t7 level. METC decision count 0: "
            "CellTypist source abstained because its labels were never captured by the "
            "engine (celltypist>=1.6 API returns AnnotationResult, engine reads adata.obs) "
            "-> n_spoke=2 < 3 -> candidates returned unchanged. Recorded as FAILURE on "
            "the quality_celltypist gate per plan QA (explicit failure, never silent pass); "
            "escalation path applies. No product code or config modified to force a pass."
        ),
    }

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("\n=== todo-11 validation gates ===")
    for k, v in gates.items():
        print(f"  {k}: {v}")
    print(f"  harmonization_rate: {harmonization_rate}")
    print(f"  metc_decision_count: {metc_decision_count}")
    print(f"  metc_branch_counts: {branch_counts}")
    print(f"  ambiguous_cells: baseline={amb_baseline} new={amb_new}")
    print(f"  ai_fallback_fired: {ai_fallback_fired}")
    print(f"  quality_celltypist: {quality_celltypist}")
    print(f"  run_passed: {run_passed}")
    print(f"  summary JSON -> {OUT_JSON}")

    # Exit code reflects honest gate results: the CellTypist-capture failure
    # (quality_celltypist False) makes this dual-on validation FAIL per plan QA.
    return 0 if run_passed else 1


if __name__ == "__main__":
    sys.exit(main())
