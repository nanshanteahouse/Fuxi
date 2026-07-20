"""
core/kb_validator.py — Empirical marker validation against the tissue KB.

Usage::

    from core.kb_validator import KbValidator    
    validator = KbValidator("retina")    
    df = validator.validate(adata, annotation_col="cell_type")
    df.to_csv("validation.csv", index=False)

CLI::

    python core/kb_validator.py --h5ad data.h5ad --annotation cell_type
    python core/kb_validator.py --h5ad data.h5ad --annotation cell_type --output result.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
import os
from typing import Optional

import pandas as pd
from anndata import AnnData

_log = logging.getLogger(__name__)

# Repo root so core imports resolve regardless of invocation directory
_script_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(_script_dir)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from core.kb import load_kb                   # noqa: E402
from core.annotation.scoring import _normalize_gene_name   # noqa: E402

# Threshold: at least 30% of cells must express a marker for it to validate
DEFAULT_PCT_THRESHOLD: float = 0.3


class KbValidator:
    """Cross-validate scRNA-seq cell-type annotations against KB marker genes.

    For each cell type in the annotated AnnData, looks up its markers in the
    tissue knowledge base and computes the fraction of cells expressing each
    marker.  Markers with ``pct_expressed >= pct_threshold`` are marked
    ``validated``; those below the threshold are ``not_validated``.

    Parameters
    ----------
    tissue : str
        Tissue identifier (e.g. ``"retina"``).
    pct_threshold : float
        Minimum fraction of cells expressing a gene for it to pass (default 0.3).
    use_ontology : bool
        If True (default), use ``StandardOntology`` to fuzzy-match annotation
        labels to KB cell-type keys.  If False, match labels directly.
    """

    def __init__(
        self,
        tissue: str = "retina",
        pct_threshold: float = DEFAULT_PCT_THRESHOLD,
        use_ontology: bool = True,
    ) -> None:
        self.tissue = tissue
        self.pct_threshold = pct_threshold
        self.kb = load_kb(tissue)
        self.ontology = None

        if use_ontology:
            try:
                from core.annotation.standardizer import StandardOntology
                self.ontology = StandardOntology(tissue)
            except (ImportError, NotImplementedError, ValueError):
                _log.info("StandardOntology unavailable for '%s' — direct matching only", tissue)

    # ── Public API ───────────────────────────────────────────────────────

    def validate(
        self,
        adata: AnnData,
        annotation_col: str = "cell_type",
    ) -> pd.DataFrame:
        """Validate KB markers against annotated cells.

        Parameters
        ----------
        adata : AnnData
            Annotated data with cell-type labels in ``.obs[annotation_col]``.
        annotation_col : str
            Column in ``adata.obs`` containing cell-type labels.

        Returns
        -------
        pd.DataFrame
            Columns: ``cell_type``, ``gene``, ``tier``, ``validated``,
            ``pct_expressed``, ``mean_expression``.
        """
        # 1. Map annotation labels → KB canonical keys
        kb_keys = self._map_labels_to_kb(adata.obs[annotation_col])

        # 2. Validate markers for each mapped cell type
        rows: list[dict] = []
        for original_label, kb_key in kb_keys.items():
            kb_entry = self.kb.get(kb_key)
            if kb_entry is None:
                continue

            cell_mask = (adata.obs[annotation_col] == original_label)
            cell_subset = adata[cell_mask]

            for tier in ("confirm", "add", "refine"):
                markers = kb_entry.get("markers", {}).get(tier, {})
                if isinstance(markers, str):
                    # Handle old-style string markers (shouldn't happen, but safe)
                    rows.extend(
                        self._validate_single_gene(
                            tier, kb_key, markers,
                            cell_subset,
                        )
                    )
                elif isinstance(markers, dict):
                    for gene in markers:
                        rows.extend(
                            self._validate_single_gene(
                                tier, kb_key, gene, cell_subset,
                            )
                        )

        if not rows:
            return pd.DataFrame(columns=[
                "cell_type", "gene", "tier", "validated",
                "pct_expressed", "mean_expression",
            ])

        return pd.DataFrame(rows)

    # ── Helpers ──────────────────────────────────────────────────────────

    def _map_labels_to_kb(self, labels: pd.Series) -> dict[str, str]:
        """Map each unique annotation label to a KB canonical key.

        When the ontology is available, uses ``standardize()``; otherwise
        does a direct membership check against the KB keys.
        """
        mapping: dict[str, str] = {}
        for label in labels.dropna().unique():
            label_str = str(label).strip()
            if not label_str:
                continue
            if self.ontology is not None:
                kb_key, _, _ = self.ontology.standardize(label_str)
                if kb_key in self.kb:
                    mapping[label_str] = kb_key
            else:
                # Direct key match
                if label_str in self.kb:
                    mapping[label_str] = label_str
        return mapping

    def _validate_single_gene(
        self,
        tier: str,
        cell_type: str,
        gene: str,
        cell_subset: AnnData,
    ) -> list[dict]:
        """Validate one marker gene against a cell subset.

        Returns a list of dicts (typically one per gene, but callers may
        append multiple — e.g. for tooltips).
        """
        matched_var = self._find_gene(gene, cell_subset.var_names)

        if matched_var is None:
            return [{
                "cell_type": cell_type,
                "gene": gene,
                "tier": tier,
                "validated": False,
                "pct_expressed": 0.0,
                "mean_expression": 0.0,
            }]

        # Extract expression vector
        X = cell_subset[:, matched_var].X
        if hasattr(X, "toarray"):
            expr = X.toarray().flatten()  # type: ignore[union-attr]
            expr = X.toarray().flatten()
        else:
            expr = X.flatten()

        pct = float((expr > 0).mean())
        mean = float(expr.mean())
        validated = pct >= self.pct_threshold

        return [{
            "cell_type": cell_type,
            "gene": gene,
            "tier": tier,
            "validated": validated,
            "pct_expressed": round(pct, 4),
            "mean_expression": round(mean, 4),
        }]

    @staticmethod
    def _find_gene(gene: str, var_names: pd.Index) -> Optional[str]:
        """Case-insensitive gene matching using ``_normalize_gene_name()``.

        Strips Macaca-specific Ensembl suffixes (``_p``, ``_n``) and compares
        uppercase-normalised names.
        """
        target = _normalize_gene_name(str(gene))
        for v in var_names:
            if _normalize_gene_name(str(v)) == target:
                return str(v)
        return None


# ═════════════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════════════


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Empirical KB marker validation for annotated scRNA-seq data.",
    )
    p.add_argument(
        "--h5ad", required=True, type=str,
        help="Path to annotated .h5ad file (e.g. projects/rna/GSE123456/results/h5ad/05_annotated.h5ad)",
    )
    p.add_argument(
        "--annotation", default="cell_type", type=str,
        help="Column in adata.obs with cell-type labels (default: cell_type)",
    )
    p.add_argument(
        "--tissue", default="retina", type=str,
        help="Tissue KB to validate against (default: retina)",
    )
    p.add_argument(
        "--output", "-o", default=None, type=str,
        help="CSV output path; if omitted, prints summary to stdout",
    )
    p.add_argument(
        "--no-ontology", action="store_true",
        help="Disable StandardOntology fuzzy-matching; match labels directly",
    )
    return p


def _print_summary(df: pd.DataFrame, tissue: str) -> None:
    """Print a concise validation summary to stdout."""
    total = len(df)
    validated = int(df["validated"].sum())
    overall_rate = validated / total * 100 if total else 0.0

    print(f"KB Validation: {tissue}")
    print(f"  Total markers tested: {total}")
    print(f"  Validated:  {validated}  ({overall_rate:.1f}%)")
    print(f"  Not validated:  {total - validated}")

    # Per-tier summary
    for tier in ("confirm", "add", "refine"):
        subset = df[df["tier"] == tier]
        if len(subset) == 0:
            continue
        t_validated = int(subset["validated"].sum())
        t_rate = t_validated / len(subset) * 100
        print(f"  {tier:>7s}:  {t_validated}/{len(subset)}  ({t_rate:.1f}%)")

    # Per-cell-type summary
    print("\n  Per cell type:")
    for ct in sorted(df["cell_type"].unique()):
        subset = df[df["cell_type"] == ct]
        ct_validated = int(subset["validated"].sum())
        ct_rate = ct_validated / len(subset) * 100 if len(subset) else 0.0
        print(f"    {ct}:  {ct_validated}/{len(subset)}  ({ct_rate:.1f}%)")



def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    # Avoid huge scanpy import for help/version
    import scanpy as sc

    if not os.path.exists(args.h5ad):
        parser.error(f"h5ad file not found: {args.h5ad}")

    adata = sc.read_h5ad(args.h5ad)
    if args.annotation not in adata.obs.columns:
        available = ", ".join(adata.obs.columns.tolist())
        parser.error(
            f"Annotation column '{args.annotation}' not found in adata.obs. "
            f"Available: {available}"
        )

    validator = KbValidator(
        tissue=args.tissue,
        use_ontology=not args.no_ontology,
    )
    df = validator.validate(adata, annotation_col=args.annotation)

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Validation results written to: {args.output}")

    _print_summary(df, args.tissue)


if __name__ == "__main__":
    main()

# ═════════════════════════════════════════════════════════════════════════════
# update_yaml_audit — Write validation results back into YAML source audit sections
# ═════════════════════════════════════════════════════════════════════════════



def update_yaml_audit(
    yaml_path: str,
    validation_df: pd.DataFrame,
    dataset_id: str,
) -> None:
    """Write validation results into a YAML source file's ``audit`` sections.

    For each cell type in the YAML that has matching validation data, updates::

    * ``expression_validated`` -- appends ``dataset_id`` if not already present.
    * ``pct_expressed`` -- per-gene ``{gene: pct_value}`` for validated (passing) genes.
    * ``cross_species_validated`` -- per-marker dict ``{gene: true/false}``.
    * ``last_audited`` -- today's date (ISO-8601).
    * ``flagged`` / ``flagged_reason`` -- set when any marker falls below threshold.

    Cell-type name matching: YAML keys are first checked directly, then resolved
    through the merged KB synonym table so that e.g. "Retinal_Ganglion_Cell" in
    the YAML matches validation results for "RGC".

    Parameters
    ----------
    yaml_path : str
        Path to the YAML source file
        (e.g. "rna/tissue_ontologies/retina/sources/hu2019.yaml").
    validation_df : pd.DataFrame
        DataFrame from ``KbValidator.validate()`` with columns
        ``cell_type, gene, tier, validated, pct_expressed, mean_expression``.
    dataset_id : str
        Dataset identifier (e.g. "GSE123456") to record in ``expression_validated``.
    """
    import yaml
    from datetime import date

    if not os.path.exists(yaml_path):
        _log.warning("YAML file not found: %s", yaml_path)
        return

    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    # Break YAML anchors (&id001 / *id001) by assigning new independent lists
    # to each cell type's audit section.  deepcopy() preserves shared refs so
    # we must do this explicitly.
    markers_section = data.get("markers", {})
    if not markers_section:
        return

    for ct_entry in markers_section.values():
        audit = ct_entry.get("audit")
        if isinstance(audit, dict):
            ev = audit.get("expression_validated")
            if isinstance(ev, list):
                audit["expression_validated"] = list(ev)
            sv = audit.get("supplement_verified")
            if isinstance(sv, list):
                audit["supplement_verified"] = list(sv)

    today = date.today().isoformat()

    # --- Build name->canonical-KB-key lookup from the merged KB -------------
    tissue = data.get("source_meta", {}).get("tissue", "retina")
    name_to_key: dict[str, str] = {}
    try:
        kb = load_kb(tissue)
        for key, entry in kb.items():
            if isinstance(entry, dict):
                name_to_key[key] = key
                for syn in entry.get("synonyms", []):
                    name_to_key[syn] = key
    except Exception:
        kb = {}

    # --- Iterate each cell type in the YAML ---------------------------------
    for yaml_ct, ct_entry in markers_section.items():
        # Resolve YAML cell-type name -> canonical KB key
        canonical_key = name_to_key.get(yaml_ct)
        if canonical_key is None:
            if yaml_ct in kb:
                canonical_key = yaml_ct
            elif yaml_ct in validation_df["cell_type"].values:
                canonical_key = yaml_ct
            else:
                continue

        # Get validation rows for this cell type
        ct_val = validation_df[validation_df["cell_type"] == canonical_key]
        if ct_val.empty:
            continue

        # --- Guarantee an audit dict ----------------------------------------
        audit = ct_entry.setdefault("audit", {})
        if not isinstance(audit, dict):
            audit = {}
            ct_entry["audit"] = audit

        # 1. expression_validated -- track which dataset validated this type
        ev = audit.get("expression_validated")
        if ev is None or not isinstance(ev, list):
            ev = []
            audit["expression_validated"] = ev
        if dataset_id not in ev:
            ev.append(dataset_id)

        # 2. pct_expressed -- {gene: pct} for passing markers (merge)
        validated = ct_val[ct_val["validated"]]
        if len(validated) > 0:
            existing_pct = audit.get("pct_expressed", {})
            if not isinstance(existing_pct, dict):
                existing_pct = {}
            existing_pct.update(
                dict(zip(validated["gene"], validated["pct_expressed"])),
            )
            audit["pct_expressed"] = existing_pct

        # 3. cross_species_validated -- per-marker dict (merge)
        existing_csv = audit.get("cross_species_validated", {})
        if not isinstance(existing_csv, dict):
            existing_csv = {}
        existing_csv.update(
            dict(zip(ct_val["gene"], ct_val["validated"])),
        )
        audit["cross_species_validated"] = existing_csv

        # 4. last_audited
        audit["last_audited"] = today

        # 5. Flag markers below threshold
        non_validated = ct_val[~ct_val["validated"]]
        if not non_validated.empty:
            audit["flagged"] = True
            audit["flagged_reason"] = "pct_expression_below_threshold"
        else:
            audit["flagged"] = False
            audit.pop("flagged_reason", None)

    # --- Write back ---------------------------------------------------------
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.dump(
            data,
            f,
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=True,
        )
    _log.info("Wrote audit updates to %s", yaml_path)
