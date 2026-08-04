#!/usr/bin/env python3
"""
F5 acceptance-gate validation for the GSE246169 fetal KADP-only run (todo 7).

Read-only validation against:
  - baseline obs snapshot : .omo/evidence/kadp-metc/w2-t7/cell_metadata_baseline.csv
  - input cluster ids    : projects/rna/GSE246169/fetal_kadp/results/h5ad/04_clustered.h5ad
  - new run obs          : projects/rna/GSE246169/fetal_kadp/results/tables/cell_metadata.csv

Asserts gates (a)-(e) from the plan (Option A formulation, 2026-08-04) and writes
a machine-checkable summary JSON to .omo/evidence/kadp-metc/w2-t7/gse246169_fetal_kadp.json.

Usage: .venv/bin/python adhoc/kadp_f5_gate.py
Exit code 0 iff all gates pass.
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
    REPO, ".omo/evidence/kadp-metc/w2-t7/cell_metadata_baseline.csv"
)
CLUSTER_H5AD = os.path.join(
    REPO, "projects/rna/GSE246169/fetal_kadp/results/h5ad/04_clustered.h5ad"
)
NEW_CSV = os.path.join(
    REPO, "projects/rna/GSE246169/fetal_kadp/results/tables/cell_metadata.csv"
)
OUT_JSON = os.path.join(REPO, ".omo/evidence/kadp-metc/w2-t7/gse246169_fetal_kadp.json")

TARGET_CLUSTERS = {"2", "4", "11"}
METC_CANDIDATES = ["0", "5", "10"]


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
    """Collapse a per-cell df to a per-cluster summary (cluster str → dict)."""
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


def main():
    prog_pole, terminal_set = _load_kb_poles()
    base, new = _load_baseline_new()
    b = _cluster_summary(base)
    n = _cluster_summary(new)
    clusters = sorted(set(b) | set(n), key=int)

    # ── Gate (a): Option A naming ──────────────────────────────────────
    a_naming = True
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
    for cl in sorted(TARGET_CLUSTERS):
        rec = n[cl]
        ok = (
            rec["method"] == "developmental_potency"
            and rec["cell_type"] in prog_pole
            and rec["cell_type"] != ""
        )
        if not ok:
            a_naming = False
            print(f"[FAIL] gate(a): cluster {cl} not KADP-named: {rec}")
    for cl in METC_CANDIDATES:
        rec = n[cl]
        if rec["method"] != "ambiguous":
            a_naming = False
            print(f"[FAIL] gate(a): cluster {cl} expected ambiguous, got {rec['method']}")

    # ── Gate (b): label-invariant (zero misnaming) ─────────────────────
    terminal_label_clusters = sorted(
        (cl for cl in clusters if b[cl]["cell_type"] in terminal_set), key=int
    )
    misnaming_count = 0
    for cl in terminal_label_clusters:
        new_ct = n[cl]["cell_type"]
        if new_ct in prog_pole:
            misnaming_count += 1
            print(
                f"[FAIL] gate(b): terminal cluster {cl} ({b[cl]['cell_type']}) "
                f"renamed to progenitor {new_ct}"
            )
    b_label_invariant = misnaming_count == 0

    # ── Gate (c): no KADP naming outside target set ────────────────────
    spillover = []
    for cl in clusters:
        if cl in TARGET_CLUSTERS:
            continue
        if n[cl]["cell_type"] != b[cl]["cell_type"]:
            spillover.append(
                (cl, b[cl]["cell_type"], n[cl]["cell_type"], n[cl]["method"])
            )
    c_no_spillover = not spillover
    if spillover:
        print(f"[FAIL] gate(c): spillover outside {{2,4,11}}: {spillover}")

    # ── Gate (d): zero transition_state decisions ──────────────────────
    transitionals = [
        cl
        for cl in clusters
        if n[cl]["method"] == "transition_state"
        or str(n[cl]["cell_type"]).startswith("transitional")
    ]
    d_zero_transitional = not transitionals
    if transitionals:
        print(f"[FAIL] gate(d): transitional decisions present: {transitionals}")

    # ── Gate (e): ambiguous cell count drop ────────────────────────────
    amb_baseline = sum(b[cl]["n"] for cl in clusters if b[cl]["method"] == "ambiguous")
    amb_new = sum(n[cl]["n"] for cl in clusters if n[cl]["method"] == "ambiguous")
    drop_expected = sum(n[cl]["n"] for cl in TARGET_CLUSTERS)
    e_ambiguous_drop = (
        amb_new < amb_baseline
        and amb_new == amb_baseline - drop_expected
        and amb_baseline == 4890  # t6 locked reference (same口径)
    )
    if not e_ambiguous_drop:
        print(
            f"[FAIL] gate(e): ambiguous cells baseline={amb_baseline} new={amb_new} "
            f"(expected drop {drop_expected})"
        )

    # ── AI fallback observation ────────────────────────────────────────
    # Gate: ai_enabled=True (config) but CFG.ai.ai_annotation=False (config
    # 'annotation: false') → the engine's low-conf AI fallback never fired;
    # ai_results stayed empty and the second-pass re-fusion did not run.
    ai_fallback_fired = False

    gates = {
        "a_naming": a_naming,
        "b_label_invariant": b_label_invariant,
        "c_no_spillover": c_no_spillover,
        "d_zero_transitional": d_zero_transitional,
        "e_ambiguous_drop": e_ambiguous_drop,
    }
    summary = {
        "dataset": "GSE246169",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_dir": "projects/rna/GSE246169/fetal_kadp",
        "locked_thresholds": {
            "ratio": 2.0,
            "abs": 0.6,
            "gap": 0.1,
            "use_gap": False,
        },
        "per_cluster": per_cluster,
        "metc_candidates": METC_CANDIDATES,
        "terminal_label_clusters": terminal_label_clusters,
        "misnaming_count": misnaming_count,
        "gates": gates,
        "ambiguous_cells": {"baseline": amb_baseline, "new": amb_new},
        "ai_fallback_fired": ai_fallback_fired,
        "notes": (
            "Option A (2026-08-04): clusters 0/10 Fisher-saturated, not KADP-nameable; "
            "recorded as METC candidates with cluster 5 for Wave 3 _metc_arbitrate. "
            "KADP target = NRPC-family clusters 2/4/11 -> Proliferating_RPC. "
            "AI fallback did NOT fire: cfg.ai.enabled=true but cfg.ai.annotation=false "
            "(ai_annot_on=False) so low_conf_clusters AI pass (engine.py:981) skipped; "
            "KADP naming from the single first-pass fusion. "
            "metc_enabled omitted from config: schema field lands in todo 10; "
            "engine default getattr(..., False) is equivalent."
        ),
    }

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("\n=== F5 gates ===")
    for k, v in gates.items():
        print(f"  {k}: {v}")
    print(f"  misnaming_count: {misnaming_count}")
    print(f"  ambiguous_cells: baseline={amb_baseline} new={amb_new}")
    print(f"  ai_fallback_fired: {ai_fallback_fired}")
    print(f"  summary JSON -> {OUT_JSON}")
    if all(gates.values()):
        print("ALL GATES PASS")
        return 0
    print("GATE FAILURE — see messages above")
    return 1


if __name__ == "__main__":
    sys.exit(main())
