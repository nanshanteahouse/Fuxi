"""
rna/annotation_standardizer.py — 6-tier name-standardization system
for cell type annotations from any source (AI, Unified KB, Score_genes).

Provides:
  - StandardOntology: Load KB + synonyms → 6-tier matching → standardize names
  - map_annotations: Batch-standardize the annotations dict from ai_annotate()
  - validate: Marker cross-validation per cluster using KB markers

Usage:
    >>> from rna.annotation_standardizer import StandardOntology
    >>> std = StandardOntology("retina")
    >>> std.standardize("Muller glial cell")
    ('Muller_Glia', 'Muller Glia', 'high')
"""

from __future__ import annotations

import sys
import os
from typing import Any

# Add project root so imports (tissue_ontologies) resolve correctly
# in both pip-installed and standalone usage.
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(_script_dir)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)
if _script_dir not in sys.path:
    sys.path.append(_script_dir)

import scanpy as sc


class StandardOntology:
    """6-tier cell-type name standardizer backed by a tissue-specific
    knowledge base and comprehensive synonym dictionary.

    The class loads a KB and synonym dictionary for the requested tissue,
    then provides methods to standardize cell type names through a 6-tier
    matching algorithm, batch-process annotation dictionaries from
    :func:`ai_annotate`, and cross-validate assignments against known
    marker genes.

    Parameters
    ----------
    tissue : str
        Tissue identifier (e.g. ``"retina"``).  Must be a supported tissue
        with a KB and synonyms module available.

    Raises
    ------
    NotImplementedError
        If the requested tissue is not yet supported.

    Examples
    --------
    >>> std = StandardOntology("retina")
    >>> std.standardize("Muller glial cell")
    ('Muller_Glia', 'Muller Glia', 'high')
    >>> std.get_candidates()
    ['Amacrine Cell', 'Amacrine Precursor', ...]
    """

    SUPPORTED_TISSUES: list[str] = ["retina"]

    def __init__(self, tissue: str) -> None:
        if tissue not in self.SUPPORTED_TISSUES:
            raise NotImplementedError(
                f"Tissue '{tissue}' is not supported. "
                f"Supported tissues: {self.SUPPORTED_TISSUES}"
            )
        self._tissue = tissue
        self._synonyms = self._load_synonyms(tissue)
        self._kb = self._load_kb(tissue)
        self._syn_map = self._build_syn_map()
        self._canonical = self._build_canonical()

    # ── Private loading helpers ──────────────────────────────────────────

    @staticmethod
    def _load_synonyms(tissue: str) -> dict[str, Any]:
        """Load the synonym dictionary for **tissue**.

        Parameters
        ----------
        tissue : str
            Tissue identifier.

        Returns
        -------
        dict[str, Any]
            Synonym dict mapping canonical keys to
            ``{"display_name": str, "synonyms": list[str]}``.

        Raises
        ------
        NotImplementedError
            If **tissue** has no synonyms module.
        """
        if tissue == "retina":
            from tissue_ontologies.retina.synonyms import RETINA_SYNONYMS
            return RETINA_SYNONYMS
        raise NotImplementedError(f"No synonyms for tissue: {tissue}")

    @staticmethod
    def _load_kb(tissue: str) -> dict[str, Any]:
        """Load the knowledge base for **tissue**.

        Parameters
        ----------
        tissue : str
            Tissue identifier.

        Returns
        -------
        dict[str, Any]
            KB dict consumable by marker scoring and evidence fusion.

        Raises
        ------
        NotImplementedError
            If the tissue KB is not available.
        """
        from tissue_ontologies import load_kb
        try:
            return load_kb(tissue)
        except ValueError as exc:
            raise NotImplementedError(
                f"Tissue '{tissue}' KB not available: {exc}"
            ) from exc

    # ── Index builders ───────────────────────────────────────────────────

    def _build_syn_map(self) -> dict[str, list[str]]:
        """Build a case-insensitive synonym → [canonical_key] lookup.

        All synonym strings are lowercased and stripped for O(1) matching
        during the 6-tier standardize algorithm.

        Returns
        -------
        dict[str, list[str]]
            Lowercased synonym → list of canonical keys that include it.
        """
        syn_map: dict[str, list[str]] = {}
        for key, syn_info in self._synonyms.items():
            synonyms_list = (
                syn_info.get("synonyms", [])
                if isinstance(syn_info, dict)
                else []
            )
            for syn in synonyms_list:
                syn_lower = syn.lower().strip()
                if syn_lower not in syn_map:
                    syn_map[syn_lower] = []
                syn_map[syn_lower].append(key)
        return syn_map

    def _build_canonical(self) -> dict[str, dict[str, Any]]:
        """Build a canonical-key → metadata index from KB and synonyms.

        Each entry contains:
          - ``display_name``: human-readable name from synonyms
          - ``markers``: list of known marker genes from the KB
            (union of ``confirm``, ``add``, and ``refine`` sub-dict keys)

        Returns
        -------
        dict[str, dict[str, Any]]
            Canonical key → ``{"display_name": str, "markers": list[str]}``.
        """
        canonical: dict[str, dict[str, Any]] = {}
        for key, syn_info in self._synonyms.items():
            display_name = (
                syn_info.get("display_name", key)
                if isinstance(syn_info, dict)
                else key
            )

            # Extract marker genes from the KB for this type
            kb_entry = self._kb.get(key, {})
            markers_raw = kb_entry.get("markers", {})
            marker_genes: list[str] = []
            for sub_key in ("confirm", "add", "refine"):
                marker_genes.extend(markers_raw.get(sub_key, {}).keys())

            canonical[key] = {
                "display_name": display_name,
                "markers": marker_genes,
            }
        return canonical

    # ── Levenshtein similarity (no external deps) ────────────────────────

    @staticmethod
    def _lev_similarity(a: str, b: str) -> float:
        """Levenshtein ratio = 1 - (edit_distance / max(len(a), len(b))).

        Pure-Python implementation, no external dependencies.

        Parameters
        ----------
        a : str
            First string.
        b : str
            Second string.

        Returns
        -------
        float
            Similarity between 0.0 and 1.0.
        """
        if not a and not b:
            return 1.0
        if not a or not b:
            return 0.0

        # Ensure `a` is the shorter string for O(min(m,n)) memory
        if len(a) > len(b):
            a, b = b, a

        prev = list(range(len(b) + 1))
        curr = [0] * (len(b) + 1)

        for i, ca in enumerate(a):
            curr[0] = i + 1
            for j, cb in enumerate(b):
                cost = 0 if ca == cb else 1
                curr[j + 1] = min(
                    curr[j] + 1,       # insertion
                    prev[j + 1] + 1,   # deletion
                    prev[j] + cost,    # substitution
                )
            prev, curr = curr, prev

        distance = prev[len(b)]
        max_len = max(len(a), len(b))
        return 1.0 - (distance / max_len)

    # ── Public API ───────────────────────────────────────────────────────

    def get_candidates(self) -> list[str]:
        """Return sorted display names for LLM prompt constraint.

        Returns
        -------
        list[str]
            Sorted human-readable cell type names (length = 32 for retina).
        """
        return sorted(
            info["display_name"]
            for info in self._canonical.values()
        )

    # ── 6-tier standardize ───────────────────────────────────────────

    def standardize(self, name: str) -> tuple[str, str, str]:
        """6-tier matching → (canonical_key, display_name, confidence).

        +------+----------------------------------------+------------+
        | Tier | Strategy                               | Confidence |
        +------+----------------------------------------+------------+
        | 1    | Exact match against canonical key (CS)  | high       |
        | 2    | Exact match against display_name (CI)   | high       |
        | 3    | Exact match against synonym (CI)        | high       |
        | 4    | Contains match (substring, CI)          | medium     |
        | 5    | Levenshtein ratio > 0.85 (CI)           | medium     |
        | 6    | Fallback — return original name as-is   | low        |
        +------+----------------------------------------+------------+

        Parameters
        ----------
        name : str
            Cell type name to standardize.

        Returns
        -------
        tuple[str, str, str]
            ``(canonical_key, display_name, confidence)`` where confidence
            is one of ``"high"``, ``"medium"``, ``"low"``.
        """
        if not name or not isinstance(name, str):
            return (name, name, "low")

        name_lower = name.lower().strip()

        # Tier 1: exact match against canonical key (case-sensitive)
        if name in self._canonical:
            return (name, self._canonical[name]["display_name"], "high")

        # Tier 2: exact match against display_name (case-insensitive)
        for key, info in self._canonical.items():
            if info["display_name"].lower() == name_lower:
                return (key, info["display_name"], "high")

        # Tier 3: exact match against any synonym (case-insensitive)
        if name_lower in self._syn_map:
            key = self._syn_map[name_lower][0]
            return (key, self._canonical[key]["display_name"], "high")

        # Tier 4: contains match (synonym in name OR name in synonym, CI)
        for key, info in self._canonical.items():
            display_lower = info["display_name"].lower()
            if name_lower in display_lower or display_lower in name_lower:
                return (key, info["display_name"], "medium")
        for synonym, keys in self._syn_map.items():
            if name_lower in synonym or synonym in name_lower:
                return (keys[0], self._canonical[keys[0]]["display_name"], "medium")

        # Tier 5: Levenshtein ratio > 0.85
        best_score = 0.0
        best_match: str | None = None
        for key, info in self._canonical.items():
            score = self._lev_similarity(name_lower, key.lower())
            if score > best_score:
                best_score, best_match = score, key
            score = self._lev_similarity(name_lower, info["display_name"].lower())
            if score > best_score:
                best_score, best_match = score, key
        for synonym, keys in self._syn_map.items():
            score = self._lev_similarity(name_lower, synonym)
            if score > best_score:
                best_score, best_match = score, keys[0]
        if best_score > 0.85 and best_match is not None:
            return (best_match, self._canonical[best_match]["display_name"], "medium")

        # Tier 6: fallback — return original name as-is
        return (name, name, "low")

    # ── Batch annotations handler ────────────────────────────────────

    def map_annotations(
        self, annotations: dict[str, dict[str, Any]] | None
    ) -> dict[str, dict[str, Any]] | None:
        """Batch-standardize the annotations dict from ``ai_annotate()``.

        Adds the following keys to each cluster entry:
          - ``cell_type_std`` — standardized canonical key
          - ``cell_type_raw`` — original name from AI
          - ``marker_validation`` — initially ``"PENDING"``

        Parameters
        ----------
        annotations : dict[str, dict] | None
            Annotations dict from ``ai_annotate()`` in the shape::

                {cluster_id: {cell_type, state, subtype, confidence, reasoning}}

        Returns
        -------
        dict[str, dict] | None
            Modified annotations with ``cell_type_std``, ``cell_type_raw``,
            and ``marker_validation`` added, or ``None`` if input was ``None``.
        """
        if annotations is None:
            return None

        result: dict[str, dict[str, Any]] = {}
        for cid, ann in annotations.items():
            cell_type = ann.get("cell_type", "")
            std_key, _, _ = self.standardize(cell_type)
            result[cid] = {
                **ann,
                "cell_type_std": std_key,
                "cell_type_raw": cell_type,
                "marker_validation": "PENDING",
            }
        return result

    # ── Marker cross-validation ──────────────────────────────────────

    def validate(
        self, adata: Any, top_n: int | None = None,
        min_overlap: float | None = None,
        marginal_threshold: float | None = None,
    ) -> list[dict[str, Any]]:
        """Marker cross-validation per cluster using KB markers.

        For each cluster, the *top_n* differentially expressed genes are
        checked against the assigned cell type's KB markers.  If
        ``rank_genes_groups`` has not been computed yet it is run
        automatically using ``sc.tl.rank_genes_groups`` with
        ``method='wilcoxon'``.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.  Must have ``.obs['leiden']`` and either
            ``.obs['cell_type_std']`` or ``.obs['cell_type']``.
        top_n : int or None
            Number of top marker genes to consider per cluster.
            If None, reads from ``CFG.marker_validation_n_top_genes``
            (default: 15).
        min_overlap : float or None
            Minimum overlap ratio for PASS status.
            If None, reads from ``CFG.marker_validation_min_overlap``
            (default: 0.5).
        marginal_threshold : float or None
            Threshold for MARGINAL tier (PASS > MARGINAL > LOW > FAIL).
            If None, reads from ``CFG.marker_validation_marginal_threshold``
            (default: 0.25).  Set to 0 to disable MARGINAL tier.

        Returns
        -------
        list[dict]
            One dict per cluster with keys:
              - ``cluster`` — cluster ID
              - ``assigned_type`` — cell type assigned to this cluster
              - ``markers_found`` — number of KB markers found in top *top_n*
              - ``markers_total`` — total KB markers for this cell type
              - ``status`` — ``"PASS"`` (≥min_overlap), ``"MARGINAL"``
                (≥marginal_threshold), ``"LOW"`` (>0.0),
                ``"FAIL"`` (0.0), or ``"NO_ONTOLOGY"`` (no KB markers)
              - ``score`` — overlap ratio (0.0 – 1.0)
        """
        # Resolve thresholds: explicit args > CFG > built-in defaults
        try:
            from core.config import CFG
            _top_n = (
                top_n if top_n is not None
                else getattr(CFG, 'marker_validation_n_top_genes', 15)
            )
            _min_overlap = (
                min_overlap if min_overlap is not None
                else getattr(CFG, 'marker_validation_min_overlap', 0.5)
            )
            _marginal = (
                marginal_threshold if marginal_threshold is not None
                else getattr(CFG, 'marker_validation_marginal_threshold', 0.25)
            )
        except (ImportError, AttributeError):
            _top_n = top_n if top_n is not None else 15
            _min_overlap = min_overlap if min_overlap is not None else 0.5
            _marginal = marginal_threshold if marginal_threshold is not None else 0.25
        # Ensure rank_genes_groups is available
        if "rank_genes_groups" not in adata.uns:
            sc.tl.rank_genes_groups(
                adata, groupby="leiden", method="wilcoxon"
            )

        # Determine which obs column holds the standardised cell type
        type_col = (
            "cell_type_std" if "cell_type_std" in adata.obs else "cell_type"
        )

        # Build cluster → majority cell type mapping
        cluster_to_type: dict[str, str] = {}
        for cluster in sorted(
            adata.obs["leiden"].unique(), key=lambda x: int(x)
        ):
            mask = adata.obs["leiden"] == cluster
            type_vals = adata.obs.loc[mask, type_col]
            try:
                majority = type_vals.mode()
                cluster_to_type[str(cluster)] = (
                    majority.iloc[0] if len(majority) > 0 else "unknown"
                )
            except Exception:
                cluster_to_type[str(cluster)] = "unknown"

        results: list[dict[str, Any]] = []
        for cluster_str in sorted(
            cluster_to_type.keys(), key=lambda x: int(x)
        ):
            assigned_type = cluster_to_type[cluster_str]

            # Get top_n marker genes for this cluster
            try:
                marker_df = sc.get.rank_genes_groups_df(
                    adata, group=cluster_str
                )
                top_genes = [
                    g
                    for g in marker_df.head(_top_n)["names"].tolist()
                    if g
                ]
            except (KeyError, ValueError):
                top_genes = []

            # Look up KB markers for the assigned type
            std_key, _, _ = self.standardize(assigned_type)
            kb_markers: list[str] = []
            if std_key in self._canonical:
                kb_markers = self._canonical[std_key].get("markers", [])

            if not kb_markers:
                results.append({
                    "cluster": cluster_str,
                    "assigned_type": assigned_type,
                    "markers_found": 0,
                    "markers_total": 0,
                    "status": "NO_ONTOLOGY",
                    "score": 0.0,
                })
                continue

            # Calculate overlap between top marker genes and KB markers
            top_set = {g.upper() for g in top_genes}
            kb_set = {g.upper() for g in kb_markers if g}
            overlap = top_set & kb_set
            score = len(overlap) / len(kb_set) if kb_set else 0.0

            if score >= _min_overlap:
                status = "PASS"
            elif _marginal > 0 and score >= _marginal:
                status = "MARGINAL"
            elif score > 0.0:
                status = "LOW"
            else:
                status = "FAIL"

            results.append({
                "cluster": cluster_str,
                "assigned_type": assigned_type,
                "markers_found": len(overlap),
                "markers_total": len(kb_set),
                "status": status,
                "score": round(score, 4),
            })

        return results
