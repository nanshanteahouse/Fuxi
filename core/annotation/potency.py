"""KADP developmental-potency pure functions (layer 3).

A multi-peak / ambiguous cluster is a plausible *transitional* population when
its marker scores are dominated by **progenitor** KB types rather than terminal
(Neuron / Glia / Non-neural) types.  This module computes that dominance as
three independent variants over the poles derived from ``kb["_hierarchy"]``:

- ``ratio`` — ``max_prog / max(max_term, epsilon)``
- ``abs``   — ``max_prog`` (with a ``max_prog > max_term`` guard against
  saturation false-positives)
- ``gap``   — ``max_prog - max_term``

All functions are **pure** (no IO, no rna-layer imports); ``core`` must not
depend on ``rna``.  ``FusionDecision.potency`` and ``annot_evidence`` consume
the three-value dict produced by :func:`to_potency_dict` / ``PotencyResult.to_potency_dict``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

__all__ = [
    "KADPConfig",
    "PotencyResult",
    "derive_developmental_poles",
    "filter_pole_scores",
    "compute_potency",
    "evaluate_passes",
    "to_potency_dict",
]


# ═══════════════════════════════════════════════════════════════════════
#  Configuration
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class KADPConfig:
    """Configuration for the KADP potency axis.

    Thresholds are locked by calibration (todo 6/7); defaults keep KADP
    disabled so existing pipelines are bit-for-bit unchanged.
    """

    enabled: bool = False
    ratio_threshold: float = 2.0
    abs_threshold: float = 0.6
    gap_threshold: float = 0.1
    use_gap_criterion: bool = False
    epsilon: float = 1e-9


# ═══════════════════════════════════════════════════════════════════════
#  Result
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class PotencyResult:
    """Potency of one cluster along the progenitor → terminal axis.

    Attributes
    ----------
    max_prog, max_term : float
        Highest positive score among progenitor / terminal pole members
        (0.0 when a pole has no positive scores).
    ratio, gap : float
        ``max_prog / max(max_term, epsilon)`` and ``max_prog - max_term``.
    best_progenitor_type : str or None
        Argmax progenitor type (deterministic tie-break); ``None`` when
        ``max_prog <= 0`` — naming requires ``max_prog > 0``.
    passes_ratio, passes_abs, passes_gap : bool
        Per-variant pass decisions (hard-coded; ``passes_abs`` includes the
        ``max_prog > max_term`` guard).
    """

    max_prog: float
    max_term: float
    ratio: float
    gap: float
    best_progenitor_type: Optional[str]
    passes_ratio: bool
    passes_abs: bool
    passes_gap: bool

    def to_potency_dict(self) -> dict[str, float]:
        """Three-value dict for ``FusionDecision.potency`` / ``annot_evidence``.

        All three variants are preserved — never a single float.
        """
        return {"ratio": self.ratio, "abs": self.max_prog, "gap": self.gap}


# ═══════════════════════════════════════════════════════════════════════
#  Pole derivation
# ═══════════════════════════════════════════════════════════════════════

_PROGENITOR_CATEGORY = "Progenitor"
_TERMINAL_CATEGORIES = ("Neuron", "Glia", "Non-neural")


def derive_developmental_poles(kb: Optional[dict[str, Any]]) -> tuple[set[str], set[str]]:
    """Derive the progenitor and terminal poles from the KB hierarchy.

    ``Progenitor`` category members form the progenitor pole; members of the
    ``Neuron``, ``Glia`` and ``Non-neural`` categories form the terminal pole.
    Returns ``(set(), set())`` when ``kb`` carries no ``_hierarchy``.

    Parameters
    ----------
    kb : dict or None
        The knowledge-base dict (raw or hierarchy-attached) — only
        ``kb["_hierarchy"]["categories"]`` is read.
    """
    categories = (kb or {}).get("_hierarchy", {}).get("categories")
    if not categories:
        return set(), set()

    progenitor = set(categories.get(_PROGENITOR_CATEGORY, {}).get("members", []))
    terminal: set[str] = set()
    for cat in _TERMINAL_CATEGORIES:
        terminal.update(categories.get(cat, {}).get("members", []))
    return progenitor, terminal


# ═══════════════════════════════════════════════════════════════════════
#  Pole scoring
# ═══════════════════════════════════════════════════════════════════════


def filter_pole_scores(pole_members, marker_scores: dict[str, Any]) -> dict[str, float]:
    """Keep pole members with a strictly positive score in ``marker_scores``.

    Ghost members (present in the pole but absent from ``marker_scores``) are
    dropped, and so are zero-scored entries — ``score_cluster_against_kb``
    emits ``Score(0.0, ...)`` objects for types without markers, which must
    not inflate a pole.

    Parameters
    ----------
    pole_members : iterable of str
        Pole members from :func:`derive_developmental_poles`.
    marker_scores : dict
        ``{type_key: score}`` — values may be raw floats or
        ``core.annotation.scoring.Score`` objects (duck-typed via ``.score``).

    Returns
    -------
    dict
        ``{type_key: float_score}`` for members with ``score > 0``.
    """
    filtered: dict[str, float] = {}
    for member in pole_members:
        if member not in marker_scores:
            continue
        raw = marker_scores[member]
        score = float(getattr(raw, "score", raw))
        if score > 0:
            filtered[member] = score
    return filtered


# ═══════════════════════════════════════════════════════════════════════
#  Potency computation
# ═══════════════════════════════════════════════════════════════════════


def compute_potency(
    marker_scores: dict[str, Any],
    poles: tuple[set[str], set[str]],
    cfg: KADPConfig,
) -> PotencyResult:
    """Compute the three potency variants for one cluster.

    Parameters
    ----------
    marker_scores : dict
        ``{type_key: score}`` — the per-cluster KB scoring output (Score
        objects or raw floats), e.g. from ``score_cluster_against_kb``.
    poles : tuple of 2 sets
        ``(progenitor_members, terminal_members)`` from
        :func:`derive_developmental_poles`.
    cfg : KADPConfig
        Thresholds and epsilon.

    Returns
    -------
    PotencyResult
        All three variants, the argmax progenitor type (``None`` when
        ``max_prog <= 0``) and the three hard-coded pass decisions.
    """
    prog_members, term_members = poles
    prog_scores = filter_pole_scores(prog_members, marker_scores)
    term_scores = filter_pole_scores(term_members, marker_scores)

    max_prog = max(prog_scores.values(), default=0.0)
    max_term = max(term_scores.values(), default=0.0)

    ratio = max_prog / max(max_term, cfg.epsilon)
    gap = max_prog - max_term

    # Naming precondition: max_prog > 0 (F9).  filter_pole_scores already
    # guarantees positivity, kept explicit per the plan spec.
    best_progenitor_type: Optional[str] = None
    if max_prog > 0 and prog_scores:
        # Deterministic argmax: highest score, tie → lexicographically
        # largest type name (set iteration order is otherwise arbitrary).
        best_progenitor_type = max(prog_scores.items(), key=lambda kv: (kv[1], kv[0]))[0]

    passes_ratio = ratio >= cfg.ratio_threshold
    passes_gap = gap >= cfg.gap_threshold
    # Oracle r2 MINOR 1: a one-sided abs check false-passes on saturated data —
    # the progenitor score must be strictly higher than the terminal score.
    passes_abs = max_prog >= cfg.abs_threshold and max_prog > max_term

    return PotencyResult(
        max_prog=max_prog,
        max_term=max_term,
        ratio=ratio,
        gap=gap,
        best_progenitor_type=best_progenitor_type,
        passes_ratio=passes_ratio,
        passes_abs=passes_abs,
        passes_gap=passes_gap,
    )


def evaluate_passes(result: PotencyResult, cfg: KADPConfig) -> bool:
    """Hard-coded combination of the three variant decisions.

    ``passes_ratio or (use_gap_criterion and passes_gap) or passes_abs``.
    Extracted as a pure function so fusion (todo 4) and this module share one
    source of truth for the written-in-stone rule.
    """
    return (
        result.passes_ratio or (cfg.use_gap_criterion and result.passes_gap) or result.passes_abs
    )


def to_potency_dict(result: PotencyResult) -> dict[str, float]:
    """Module-level convenience: ``PotencyResult → {"ratio", "abs", "gap"}``."""
    return result.to_potency_dict()
