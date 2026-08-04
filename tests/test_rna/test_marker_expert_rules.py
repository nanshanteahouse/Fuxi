"""Tests for rna/utils/marker_expert_rules.py — apply_expert_rules corroboration.

Covers the D4 ``corroborators`` (any-of) schema on top of
``markers_present``:

1. Corroboration hit — a corroborator in the same top-N + pval ``de_subset``
   yields ``corroborated=True`` with the hit recorded.
2. Uncorroborated — markers present but no corroborator in the subset: the
   match is *returned* (never silently rejected) with ``corroborated=False``;
   the arbitration/downgrade is the caller's job (fuse_evidence, plan task 7).
3. Legacy rules without a ``corroborators`` field stay ``corroborated=True``.
4. The corroborator check shares the *identical* de_subset (top-N + pval
   filter) used for the primary marker check — genes outside it never
   corroborate, whether excluded by rank or by pvals_adj.
5. Hits are copies of the KB rule dict — no cross-cluster pollution when the
   same kb is reused.
6. ``pvals_adj`` column absent → pval filter silently skipped; corroboration
   still works.

The downgrade-to-weak arbitration lives in fuse_evidence (plan task 7) and
is deliberately out of scope here — these tests pin the flag semantics of
``corroborated`` / ``corroborators_hit`` at the apply_expert_rules layer.
"""

import pandas as pd

from rna.utils.marker_expert_rules import apply_expert_rules

# ── Helpers ────────────────────────────────────────────────────────────


def _make_rule(
    action: str = "Type_A",
    markers_present: dict[str, float] | None = None,
    corroborators: list[str] | None = None,
    markers_absent: list[str] | None = None,
    priority: int = 1,
) -> dict:
    """Build an expert-rule dict in the KB schema (priority/condition/action)."""
    condition: dict = {"markers_present": markers_present or {"A": 1.0}}
    if markers_absent:
        condition["markers_absent"] = markers_absent
    if corroborators:
        condition["corroborators"] = corroborators
    return {"priority": priority, "condition": condition, "action": action}


def _make_kb(rules: list[dict]) -> dict:
    return {"expert_rules": rules}


def _de_df(
    gene_order: list[str],
    lfc_map: dict[str, float] | None = None,
    pval_map: dict[str, float] | None = None,
    with_pvals: bool = True,
    default_pval: float = 0.001,
) -> pd.DataFrame:
    """Build a cluster DE DataFrame in rank order (rank 1 first).

    ``names`` / ``logfoldchanges`` are always present; ``pvals_adj`` is only
    added when *with_pvals* is True (mirrors a DE table that may lack it).
    """
    rows: list[dict] = []
    for i, gene in enumerate(gene_order, start=1):
        row = {"names": gene, "logfoldchanges": (lfc_map or {}).get(gene, 1.0 + 0.1 * i)}
        if with_pvals:
            row["pvals_adj"] = (pval_map or {}).get(gene, default_pval)
        rows.append(row)
    return pd.DataFrame(rows)


# ── Import / API ───────────────────────────────────────────────────────


class TestApplyExpertRulesImport:
    def test_import_apply_expert_rules(self) -> None:
        assert apply_expert_rules is not None


# ── Case 1: corroboration hit ──────────────────────────────────────────


class TestCorroborationHit:
    def test_corroborated_when_corroborator_in_subset(self) -> None:
        """Given a rule markers_present={A:1.0} corroborators=[B, C] and a
        cluster whose top-50 DE contains both A and B, When apply_expert_rules
        runs, Then the rule matches with corroborated=True and
        corroborators_hit == [B]."""
        kb = _make_kb([_make_rule(markers_present={"A": 1.0}, corroborators=["B", "C"])])
        de = _de_df(["A", "B", "F1", "F2", "F3"], lfc_map={"A": 2.0, "B": 1.5})

        action, all_matched = apply_expert_rules(kb, de)

        assert action == "Type_A"
        assert len(all_matched) == 1
        hit = all_matched[0]
        assert hit["corroborated"] is True
        assert hit["corroborators_hit"] == ["B"]


# ── Case 2: uncorroborated ─────────────────────────────────────────────


class TestUncorroborated:
    def test_uncorroborated_match_is_returned_not_rejected(self) -> None:
        """Given the same rule but a top-50 containing only A, When
        apply_expert_rules runs, Then the match is still returned (the caller
        must downgrade, not be denied) with corroborated=False and
        corroborators_hit == []."""
        kb = _make_kb([_make_rule(markers_present={"A": 1.0}, corroborators=["B", "C"])])
        de = _de_df(["A", "F1", "F2", "F3"], lfc_map={"A": 2.0})

        action, all_matched = apply_expert_rules(kb, de)

        # Not a silent rejection: the winner is still reported.
        assert action == "Type_A"
        assert len(all_matched) == 1
        hit = all_matched[0]
        assert hit["corroborated"] is False
        assert hit["corroborators_hit"] == []


# ── Case 3: legacy rule without corroborators ──────────────────────────


class TestLegacyRuleNoCorroborators:
    def test_rule_without_corroborators_is_always_corroborated(self) -> None:
        """Given a traditional rule with no ``corroborators`` field, When it
        matches, Then corroborated=True (legacy behavior unchanged) and
        corroborators_hit == []."""
        kb = _make_kb([_make_rule(markers_present={"A": 1.0})])
        de = _de_df(["A", "F1", "F2"], lfc_map={"A": 2.0})

        action, all_matched = apply_expert_rules(kb, de)

        assert action == "Type_A"
        assert len(all_matched) == 1
        hit = all_matched[0]
        assert hit["corroborated"] is True
        assert hit["corroborators_hit"] == []
        # The copy must not carry a spurious corroborators key.
        assert "corroborators" not in hit["condition"]


# ── Case 4: shared de_subset constraint ────────────────────────────────


class TestSharedDeSubsetConstraint:
    def test_corroborator_beyond_top_n_never_corroborates(self) -> None:
        """Given a corroborator at rank > top_n (here rank 60) with a high
        lfc, When apply_expert_rules runs, Then it is excluded from de_subset
        by the top-N gate and does not corroborate — proving the corroborator
        check uses the same top-N subset as the primary marker check."""
        kb = _make_kb([_make_rule(markers_present={"A": 1.0}, corroborators=["B"])])
        # 60 genes: A at rank 1, B at rank 60 (index 59).
        gene_order = ["A"] + [f"F{i}" for i in range(1, 59)] + ["B"]
        de = _de_df(gene_order, lfc_map={"A": 2.0, "B": 8.0})

        # B genuinely exists in the full DE table — exclusion is a subset gate.
        assert "B" in set(de["names"])

        action, all_matched = apply_expert_rules(kb, de)

        assert action == "Type_A"
        hit = all_matched[0]
        assert hit["corroborated"] is False
        assert hit["corroborators_hit"] == []

    def test_corroborator_within_top_n_but_bad_pval_never_corroborates(self) -> None:
        """Given a corroborator inside the top-N window but failing the pval
        filter (pvals_adj 0.9 >= 0.05), When apply_expert_rules runs, Then it
        is dropped from de_subset by the pval gate and does not corroborate —
        the corroborator check shares the identical pval-filtered subset."""
        kb = _make_kb([_make_rule(markers_present={"A": 1.0}, corroborators=["B"])])
        # B at rank 21 — inside head(50) — but pvals_adj 0.9 fails the cutoff.
        gene_order = ["A"] + [f"F{i}" for i in range(1, 20)] + ["B"]
        de = _de_df(
            gene_order,
            lfc_map={"A": 2.0, "B": 5.0},
            pval_map={"A": 0.001, "B": 0.9},
        )

        action, all_matched = apply_expert_rules(kb, de)

        assert action == "Type_A"
        hit = all_matched[0]
        assert hit["corroborated"] is False
        assert hit["corroborators_hit"] == []


# ── Case 5: no KB pollution across clusters ────────────────────────────


class TestKbIsolationAcrossClusters:
    def test_first_cluster_corroborators_hit_not_leaked_to_second(self) -> None:
        """Given the same kb dict reused for two cluster calls in sequence,
        When cluster 1 corroborates via B, Then cluster 2's hit carries its
        own empty corroborators_hit — no state leaks through the shared kb."""
        rule = _make_rule(markers_present={"A": 1.0}, corroborators=["B"])
        kb = _make_kb([rule])

        de_first = _de_df(["A", "B", "F1"], lfc_map={"A": 2.0, "B": 1.5})
        action1, matched1 = apply_expert_rules(kb, de_first)

        assert action1 == "Type_A"
        assert matched1[0]["corroborated"] is True
        assert matched1[0]["corroborators_hit"] == ["B"]

        de_second = _de_df(["A", "F1"], lfc_map={"A": 2.0})
        action2, matched2 = apply_expert_rules(kb, de_second)

        assert action2 == "Type_A"
        assert matched2[0]["corroborated"] is False
        assert matched2[0]["corroborators_hit"] == []

    def test_returned_hit_is_a_copy_not_the_kb_rule(self) -> None:
        """Given a matched rule, When the caller mutates the returned dict,
        Then the KB rule stays pristine (the hit is a shallow copy), so reuse
        across clusters cannot be contaminated."""
        rule = _make_rule(markers_present={"A": 1.0}, corroborators=["B"])
        kb = _make_kb([rule])
        de = _de_df(["A", "B"], lfc_map={"A": 2.0, "B": 1.5})

        _, all_matched = apply_expert_rules(kb, de)
        hit = all_matched[0]

        assert hit is not kb["expert_rules"][0]
        hit["corroborated"] = False
        hit["corroborators_hit"] = ["INJECTED"]
        hit["extra"] = "mutated"

        assert "corroborated" not in kb["expert_rules"][0]
        assert "corroborators_hit" not in kb["expert_rules"][0]
        assert "extra" not in kb["expert_rules"][0]


# ── Case 6: pvals_adj absent ───────────────────────────────────────────


class TestMissingPvalsAdj:
    def test_pval_filter_silently_skipped_corroboration_still_works(self) -> None:
        """Given a DE table without a ``pvals_adj`` column, When
        apply_expert_rules runs, Then the pval filter is silently skipped and
        corroboration works normally on the top-N subset."""
        kb = _make_kb([_make_rule(markers_present={"A": 1.0}, corroborators=["B"])])
        de = _de_df(["A", "B", "F1"], lfc_map={"A": 2.0, "B": 1.5}, with_pvals=False)

        assert "pvals_adj" not in de.columns

        action, all_matched = apply_expert_rules(kb, de)

        assert action == "Type_A"
        hit = all_matched[0]
        assert hit["corroborated"] is True
        assert hit["corroborators_hit"] == ["B"]
