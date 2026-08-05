"""Species-aware forced-gene selection from the merged KB.

The HVG selection in ``rna/steps/03_integrate.py`` force-keeps genes listed in
``cfg.hvg.forced_genes`` (and ``cfg.marker.marker_dict``) so that low-variance
but biologically critical genes (transcription factors, sparse lineage markers)
survive the top-N variance cut.  This module turns the *corrected* KB consensus
(``core/kb/merge.py`` — effective source counts with shared-root discount,
negative-evidence penalty, consensus levels) into that gene list:

* markers are taken from the ``confirm`` tier only,
* filtered by the target species (a gene counts only when at least one of its
  supporting sources is annotated for that species),
* filtered by a minimum consensus level — the *corrected* level from
  ``consensus_effective_counts``, which already bakes in shared-root discount
  and the negative-evidence penalty (see ``merge.compute_effective_source_count``).

Negative evidence is therefore NOT re-applied here: ``RPE65`` stays ``high``
(effective 4.0) even though several sources recorded ``cross_species_validated:
False`` for it, because all of those checks collapsed onto a single shared
validation dataset (GSE118614) and the penalty is applied once.  Re-filtering
here would double-penalize.
"""

from __future__ import annotations

import logging
import os

_log = logging.getLogger(__name__)

_LEVEL_ORDER = {"low": 1, "medium": 2, "high": 3, "gold": 4}

# Ubiquitous genes excluded from forced selection: they pass the consensus
# gates (species + threshold) but are broadly expressed, so forcing them into
# the HVG set wastes a slot without adding cluster-separation signal.
# B2M (MHC-I light chain) was observed expressed in 38/38 validation datasets
# (see notes/research/2026-08-03_marker_claim_audit.md).  This list is a
# conservative denylist for the *forced-HVG* use case only; the genes remain
# valid markers for annotation/scoring.
_UBIQUITOUS_DENYLIST: frozenset[str] = frozenset({"B2M"})

# Consensus levels accepted at each threshold.
_THRESHOLDS = {
    "gold": {"gold"},
    "high": {"gold", "high"},
    "medium": {"gold", "high", "medium"},
    "any": {"gold", "high", "medium", "low"},
}


def _source_species(sources: list[dict]) -> dict[str, set[str]]:
    """Map source_id -> set of species (normalised to pipeline slugs)."""
    from core.preprocess.format_detector import _normalise_species

    out: dict[str, set[str]] = {}
    for src in sources:
        sp = {_normalise_species(s) for s in (src.get("meta", {}).get("species") or []) if s}
        out[src["source_id"]] = sp
    return out


def build_forced_genes(
    tissue: str,
    target_species: str,
    threshold: str = "high",
) -> list[str]:
    """Return a sorted, deduplicated forced-gene list for *target_species*.

    Parameters
    ----------
    tissue : str
        Tissue directory name under ``core/kb/`` (e.g. ``"retina"``).
    target_species : str
        Species filter; a gene must be supported by at least one source
        annotated for this species.  Accepts either a pipeline slug
        ("human", "mouse", "zebrafish") or a Latin binomial
        ("Homo sapiens") — normalisation happens inside.
    threshold : str
        One of ``"gold"``, ``"high"``, ``"medium"``, ``"any"``.  The level
        compared is the corrected ``consensus_levels`` (effective source
        count), not the raw source count.
    """
    from core.kb.merge import (
        build_final_kb,
        load_all_sources,
        load_source_independence,
        merge_markers,
    )

    kb_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), tissue)
    sources_dir = os.path.join(kb_dir, "sources")
    if not os.path.isdir(sources_dir):
        raise ValueError(f"No KB sources dir for tissue '{tissue}': {sources_dir}")

    sources = load_all_sources(sources_dir)
    indep = load_source_independence(os.path.join(sources_dir, "_source_independence.yaml"))
    merged = merge_markers(sources)
    kb = build_final_kb(merged, [], sources, indep)
    sp_of = _source_species(sources)

    # Normalise the requested species once to a pipeline slug; matches the
    # slugs _source_species emits, so slug and Latin binomial both work.
    from core.preprocess.format_detector import _normalise_species

    target_key = _normalise_species(target_species)

    accepted = _THRESHOLDS.get(threshold)
    if accepted is None:
        raise ValueError(f"Unknown threshold '{threshold}' (use gold/high/medium/any)")

    picked: dict[str, str] = {}  # gene -> best level (highest wins)
    for type_key, entry in kb.items():
        if not isinstance(entry, dict) or "markers" not in entry:
            continue
        conf = entry.get("markers", {}).get("confirm", {})
        levels = entry.get("consensus_levels", {})
        for gene, src_ids in conf.items():
            if not isinstance(src_ids, (list, tuple)):
                continue
            src_set = set(src_ids)
            # Species gate: at least one supporting source for target species.
            if not any(sp_of.get(s, set()) and target_key in sp_of.get(s, set()) for s in src_set):
                continue
            level = levels.get(gene, "low")
            if level not in accepted:
                continue
            if gene in _UBIQUITOUS_DENYLIST:
                continue
            # Keep the best (highest) level seen across cell types.
            cur = picked.get(gene)
            if cur is None or _LEVEL_ORDER[level] > _LEVEL_ORDER[cur]:
                picked[gene] = level

    ordered = sorted(picked)
    _log.info(
        "forced_genes[%s/%s]: %d genes (threshold=%s)",
        tissue,
        target_species,
        len(ordered),
        threshold,
    )
    return ordered


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Generate species-aware forced genes from KB")
    ap.add_argument("--tissue", default="retina")
    ap.add_argument("--species", default="Homo sapiens")
    ap.add_argument("--threshold", default="high", choices=["gold", "high", "medium", "any"])
    ap.add_argument("--yaml", action="store_true", help="Print as YAML list for config")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO)
    genes = build_forced_genes(args.tissue, args.species, threshold=args.threshold)
    if args.yaml:
        for g in genes:
            print(f"    - {g}")
    else:
        print(" ".join(genes))
