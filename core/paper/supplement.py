"""
core/supplement_parser.py — Parse supplementary table Excel files into KB marker dicts.

Handles multiple publisher formats:
  - cluster_marker:      Hu s012 (PLoS Biology) — cluster + gene + avg_diff
  - gene_score_matrix:   Menon MOESM5 (Nat Commun) — Gene x cell type score matrix
  - per_type_sheet:      Zuo MOESM6 (Nat Commun) — one sheet per cell type
                         Li S7A (Nat Genetics) — group column as classifier
  - unknown:             Peng supplement-9 (Cell) — comma-separated signature genes
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

logger = logging.getLogger(__name__)

# Column normalisation regexes
_LOGFC_PATTERN = re.compile(r"(?i)(?:avg_|average_)?log[2]?[fF][cC](?:hanges?)?|avg_diff")
_PVAL_PATTERN = re.compile(r"(?i)p(?:val(?:ue)?|_val|vals_adj|_adjusted|adj)")
_GENE_PATTERN = re.compile(r"(?i)(?:gene|feature|names|symbol|name|Gene)\s*$")

_HEADER_KEYWORDS = {
    "gene",
    "names",
    "name",
    "symbol",
    "cluster",
    "group",
    "logfoldchanges",
    "logfc",
    "scores",
    "score",
    "pvals",
    "pval",
    "pvals_adj",
    "avg_diff",
    "cell",
    "marker",
    "signature",
}

_CELL_TYPE_ABBREVS = {
    "rod",
    "rods",
    "cone",
    "cones",
    "rgc",
    "rgcs",
    "bc",
    "bp",
    "bps",
    "bipolar",
    "ac",
    "acs",
    "amacrine",
    "hc",
    "hcs",
    "horizontal",
    "mg",
    "muller",
    "macroglia",
    "microglia",
    "rpc",
    "prpc",
    "nrpc",
    "rpe",
    "astrocyte",
    "vascular",
    "endothelial",
    "pericyte",
    "fibroblast",
    "oligodendrocyte",
}


# ── Helpers ──────────────────────────────────────────────────────────────────


def _normalize_col(col_name: str) -> str:
    """Normalize a column name to a canonical form."""
    s = str(col_name).strip()
    if _GENE_PATTERN.match(s):
        return "gene"
    if _LOGFC_PATTERN.match(s):
        return "logFC"
    if _PVAL_PATTERN.match(s):
        return "pval_adj"
    if re.match(r"(?i)scores?", s):
        return "score"
    if re.match(r"(?i)cluster", s):
        return "cluster"
    if re.match(r"(?i)group", s):
        return "group"
    if re.match(r"(?i)marker", s):
        return "marker"
    return s.lower()


def _normalize_columns(df: pd.DataFrame) -> dict[str, str]:
    """Return mapping {original_col_name: normalized_name}."""
    return {c: _normalize_col(c) for c in df.columns}


def _is_cell_type_sheet_name(name: str) -> bool:
    """Check if a sheet name represents a cell type."""
    nl = name.strip().lower()
    if nl in _CELL_TYPE_ABBREVS:
        return True
    for ab in _CELL_TYPE_ABBREVS:
        if len(ab) >= 2 and re.search(r"\b" + re.escape(ab) + r"\b", nl):
            return True
    return False


def _pre_normalize_cell_type(name: str) -> str:
    """Pre-normalize cell type name before ontology matching."""
    s = name.strip()
    if s.endswith("s") and len(s) > 3 and s[:-1].lower() in _CELL_TYPE_ABBREVS:
        s = s[:-1]
    abbrev_map = {
        "bc": "Bipolar Cell",
        "bps": "Bipolar Cell",
        "bipolar": "Bipolar Cell",
        "ac": "Amacrine Cell",
        "acs": "Amacrine Cell",
        "amacrine": "Amacrine Cell",
        "hc": "Horizontal Cell",
        "hcs": "Horizontal Cell",
        "horizontal": "Horizontal Cell",
        "mg": "Muller Glia",
        "muller": "Muller Glia",
        "macroglia": "Muller Glia",
        "microglia": "Microglia",
        "rgc": "RGC",
        "rgcs": "RGC",
        "ganglion": "RGC",
        "prpc": "Proliferating RPC",
        "nrpc": "Proliferating RPC",
        "rpc": "RPC",
        "retinal ganglion cell": "RGC",
        "endothelial cell": "Vascular Endothelial",
    }
    sl = s.lower()
    if sl in abbrev_map:
        return abbrev_map[sl]
    return s


def _standardize_cell_type(name: str) -> tuple[str, str, str]:
    """Standardize a cell type name via StandardOntology."""
    try:
        from core.annotation.standardizer import StandardOntology  # noqa: PLC0415

        onto = StandardOntology("retina")
        pre_norm = _pre_normalize_cell_type(name)
        return onto.standardize(pre_norm)
    except Exception:
        pre_norm = _pre_normalize_cell_type(name)
        return (pre_norm.replace(" ", "_"), pre_norm, "low")


def _read_sheet(xls: pd.ExcelFile, sheet_name: str) -> pd.DataFrame:
    """Read a sheet from an open ExcelFile, handling merged-cell headers.

    Some publishers (e.g. Nature Genetics) put a merged title cell in row 0
    and the real headers in row 1.  This detects that pattern and re-reads
    with the correct header row.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df_raw = pd.read_excel(xls, sheet_name=sheet_name)

        if len(df_raw) == 0:
            return df_raw

        # Check if a data row contains header-like keywords (merged-cell headers)
        row0 = [str(v).lower().strip() for v in df_raw.iloc[0].tolist()]
        nn0 = [v for v in row0 if v and v != "nan"]
        kw0 = any(v in _HEADER_KEYWORDS for v in nn0)

        # Merged-title case: row 0 is all-NaN, row 1 has real headers
        kw1 = False
        if not kw0 and len(nn0) == 0 and len(df_raw) >= 2:
            row1 = [str(v).lower().strip() for v in df_raw.iloc[1].tolist()]
            nn1 = [v for v in row1 if v and v != "nan"]
            kw1 = any(v in _HEADER_KEYWORDS for v in nn1)

        if kw0 and len(nn0) >= 3:
            return pd.read_excel(xls, sheet_name=sheet_name, header=1)

        if kw1:
            return pd.read_excel(xls, sheet_name=sheet_name, header=2)

        return df_raw


# ── SupplementTableParser ────────────────────────────────────────────────────


class SupplementTableParser:
    """Parse supplementary table Excel files into KB marker dicts."""

    def __init__(
        self,
        max_genes_per_type: int = 20,
        max_sheets: int = 20,
        confirm_logfc_threshold: float = 1.0,
        add_logfc_threshold: float = 0.5,
        confirm_score_threshold: float = 2.0,
        add_score_threshold: float = 1.0,
    ):
        self.max_genes_per_type = max_genes_per_type
        self.max_sheets = max_sheets
        self.confirm_logfc_threshold = confirm_logfc_threshold
        self.add_logfc_threshold = add_logfc_threshold
        self.confirm_score_threshold = confirm_score_threshold
        self.add_score_threshold = add_score_threshold

    # -- detect_format --------------------------------------------------------

    def detect_format(self, excel_path: str) -> str:
        """Detect the supplementary table format.

        Returns: "cluster_marker" | "gene_score_matrix" | "per_type_sheet" | "unknown"
        Raises FileNotFoundError or ValueError.
        """
        path = Path(excel_path)
        if not path.exists():
            raise FileNotFoundError(f"Excel file not found: {excel_path}")

        try:
            with pd.ExcelFile(excel_path) as xls:
                return self._detect_format_internal(xls)
        except Exception as exc:
            raise ValueError(f"Cannot open Excel file: {excel_path}") from exc

    def _detect_format_internal(self, xls: pd.ExcelFile) -> str:
        sheet_names = xls.sheet_names

        ct_sheets = [s for s in sheet_names if _is_cell_type_sheet_name(s)]
        if len(ct_sheets) >= 1:
            return "per_type_sheet"

        if not sheet_names:
            return "unknown"

        try:
            df = _read_sheet(xls, sheet_names[0])
        except Exception:
            return "unknown"

        norm_map = _normalize_columns(df)
        norm_cols = set(norm_map.values())

        if "cluster" in norm_cols and "gene" in norm_cols:
            return "cluster_marker"

        if "gene" in norm_cols:
            gene_orig = next(c for c, n in norm_map.items() if n == "gene")
            numeric_cols = [
                c for c in df.columns if c != gene_orig and pd.api.types.is_numeric_dtype(df[c])
            ]
            if len(numeric_cols) >= 3:
                return "gene_score_matrix"

        stats_cols = {"score", "logfc", "pval_adj"}
        if "group" in norm_cols and "gene" in norm_cols:
            if norm_cols & stats_cols:
                return "per_type_sheet"

        return "unknown"

    # -- parse_to_kb ----------------------------------------------------------

    def parse_to_kb(
        self, excel_path: str, source_meta: dict[str, Any]
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        """Parse supplementary table into KB marker dict and enhanced metadata."""
        path = Path(excel_path)
        if not path.exists():
            raise FileNotFoundError(f"Excel file not found: {excel_path}")

        try:
            with pd.ExcelFile(excel_path) as xls:
                fmt = self._detect_format_internal(xls)
                pmid = str(source_meta.get("pmid", ""))
                enhanced = dict(source_meta)

                if fmt == "cluster_marker":
                    markers = self._parse_cluster_marker(xls, pmid)
                elif fmt == "gene_score_matrix":
                    markers = self._parse_gene_score_matrix(xls, pmid)
                elif fmt == "per_type_sheet":
                    markers = self._parse_per_type_sheet(xls, pmid)
                else:
                    markers = self._parse_unknown(xls, pmid)
        except FileNotFoundError:
            raise
        except Exception as exc:
            raise ValueError(f"Cannot open Excel file: {excel_path}") from exc

        n_types = len(markers)
        if n_types > 0:
            enhanced["n_groups"] = max(enhanced.get("n_groups", 0), n_types)

        return markers, enhanced

    # -- to_yaml_source -------------------------------------------------------

    def to_yaml_source(
        self,
        markers_dict: dict[str, Any],
        source_meta: dict[str, Any],
        output_dir: str,
    ) -> str:
        """Write markers_dict + source_meta as a YAML source file."""
        os.makedirs(output_dir, exist_ok=True)

        source_id = source_meta.get("id", "unknown_source")
        filename = f"{source_id}.yaml"
        filepath = os.path.join(output_dir, filename)

        output = {
            "source_meta": _clean_meta(source_meta),
            "markers": _clean_markers(markers_dict),
        }

        with open(filepath, "w") as f:
            yaml.dump(
                output,
                f,
                default_flow_style=False,
                sort_keys=False,
                allow_unicode=True,
                indent=2,
            )

        logger.info("YAML source written: %s", filepath)
        return filepath

    # -- Format-specific parsers ---------------------------------------------

    def _parse_cluster_marker(
        self,
        xls: pd.ExcelFile,
        pmid: str,
    ) -> dict[str, Any]:
        """Parse cluster_marker format (Hu s012)."""
        sheet = xls.sheet_names[0]
        df = _read_sheet(xls, sheet)
        norm = _normalize_columns(df)

        cluster_col = next(c for c, n in norm.items() if n == "cluster")
        gene_col = next(c for c, n in norm.items() if n == "gene")
        logfc_col = next(
            (c for c, n in norm.items() if n == "logFC"),
            next((c for c, n in norm.items() if n == "score"), None),
        )

        markers: dict[str, dict[str, Any]] = {}

        for cluster_id in sorted(df[cluster_col].unique()):
            cluster_df = df[df[cluster_col] == cluster_id].copy()
            if logfc_col:
                cluster_df = cluster_df.sort_values(logfc_col, ascending=False)

            canonical, display, conf = _standardize_cell_type(f"Cluster_{cluster_id}")
            if conf == "low":
                display = str(cluster_id)
                canonical = f"Cluster_{display}"

            confirm: dict[str, list[str]] = {}
            add: dict[str, list[str]] = {}

            for _, row in cluster_df.iterrows():
                gene = str(row[gene_col]).strip().upper()
                if not gene or gene == "NAN":
                    continue
                logfc_val = (
                    float(row[logfc_col]) if logfc_col and pd.notna(row[logfc_col]) else None
                )
                if len(confirm) + len(add) >= self.max_genes_per_type:
                    break
                if logfc_val is not None and logfc_val > self.confirm_logfc_threshold:
                    confirm.setdefault(gene, []).append(pmid)
                elif logfc_val is not None and logfc_val > self.add_logfc_threshold:
                    add.setdefault(gene, []).append(pmid)

            if confirm or add:
                markers[canonical] = {"confirm": confirm, "add": add}

        return markers

    def _parse_gene_score_matrix(
        self,
        xls: pd.ExcelFile,
        pmid: str,
    ) -> dict[str, Any]:
        """Parse gene_score_matrix format (Menon MOESM5)."""
        sheet = xls.sheet_names[0]
        df = _read_sheet(xls, sheet)

        norm = _normalize_columns(df)
        gene_col = next(c for c, n in norm.items() if n == "gene")

        ct_cols = [c for c in df.columns if c != gene_col and pd.api.types.is_numeric_dtype(df[c])]

        markers: dict[str, dict[str, Any]] = {}

        for ct_col in ct_cols:
            canonical, display, _conf = _standardize_cell_type(ct_col)
            top_df = df.nlargest(self.max_genes_per_type, ct_col)

            confirm: dict[str, list[str]] = {}
            add: dict[str, list[str]] = {}

            for _, row in top_df.iterrows():
                gene = str(row[gene_col]).strip().upper()
                if not gene or gene == "NAN":
                    continue
                score = float(row[ct_col]) if pd.notna(row[ct_col]) else 0.0
                if score > self.confirm_score_threshold:
                    confirm.setdefault(gene, []).append(pmid)
                elif score > self.add_score_threshold:
                    add.setdefault(gene, []).append(pmid)

            if confirm or add:
                markers[canonical] = {"confirm": confirm, "add": add}

        return markers

    def _parse_per_type_sheet(
        self,
        xls: pd.ExcelFile,
        pmid: str,
    ) -> dict[str, Any]:
        """Parse per_type_sheet format (Zuo + Li variants)."""
        markers: dict[str, dict[str, Any]] = {}

        # Variant: sheets with group column (Li S7A)
        for sheet in xls.sheet_names:
            if "s7a" in sheet.lower():
                try:
                    self._parse_grouped_sheet(xls, sheet, pmid, markers)
                except Exception:
                    continue

        # Variant: cell-type-named sheets
        ct_sheets = [
            s
            for s in xls.sheet_names
            if _is_cell_type_sheet_name(s) and "s7a" not in s.lower() and "s7b" not in s.lower()
        ]
        for sheet in ct_sheets[: self.max_sheets]:
            try:
                self._parse_one_type_sheet(xls, sheet, pmid, markers)
            except Exception:
                continue

        return markers

    def _parse_one_type_sheet(
        self,
        xls: pd.ExcelFile,
        sheet: str,
        pmid: str,
        markers: dict[str, dict[str, Any]],
    ) -> None:
        """Parse a single cell-type sheet (Zuo or Li per-type format)."""
        try:
            df = _read_sheet(xls, sheet)
        except Exception:
            return

        norm = _normalize_columns(df)

        # Li per-type: has 'marker' column with comma-separated genes
        if "marker" in norm.values():
            self._parse_marker_column_sheet(df, norm, sheet, pmid, markers)
            return

        # Zuo-style: gene/names column + stats
        try:
            gene_col = next(c for c, n in norm.items() if n == "gene")
        except StopIteration:
            return

        logfc_col = next(
            (c for c, n in norm.items() if n == "logFC"),
            next((c for c, n in norm.items() if n == "score"), None),
        )
        score_col = next(
            (c for c, n in norm.items() if n == "score"),
            logfc_col,
        )

        canonical, display, _conf = _standardize_cell_type(sheet)

        sort_col = score_col or logfc_col
        df_sorted = (
            df.sort_values(sort_col, ascending=False)
            if sort_col and sort_col in df.columns
            else df
        )

        confirm: dict[str, list[str]] = {}
        add: dict[str, list[str]] = {}

        for _, row in df_sorted.iterrows():
            gene = str(row[gene_col]).strip().upper()
            if not gene or gene == "NAN":
                continue
            if len(confirm) + len(add) >= self.max_genes_per_type:
                break

            logfc = (
                float(row[logfc_col])
                if logfc_col and logfc_col in df.columns and pd.notna(row[logfc_col])
                else None
            )
            score = (
                float(row[score_col])
                if score_col and score_col in df.columns and pd.notna(row[score_col])
                else logfc
            )
            sc_val = score if score is not None else logfc
            if sc_val is None:
                continue

            if (logfc is not None and logfc > self.confirm_logfc_threshold) or (
                score is not None and score > self.confirm_score_threshold
            ):
                confirm.setdefault(gene, []).append(pmid)
            elif (logfc is not None and logfc > self.add_logfc_threshold) or (
                score is not None and score > self.add_score_threshold
            ):
                add.setdefault(gene, []).append(pmid)

        if confirm or add:
            markers[canonical] = {"confirm": confirm, "add": add}

    def _parse_marker_column_sheet(
        self,
        df: pd.DataFrame,
        norm: dict[str, str],
        sheet: str,
        pmid: str,
        markers: dict[str, dict[str, Any]],
    ) -> None:
        """Parse a sheet with comma-separated 'marker' column (Li per-type)."""
        marker_orig = next(c for c, n in norm.items() if n == "marker")
        ct_orig = next(
            (c for c, n in norm.items() if n == "celltype"),
            next((c for c, n in norm.items() if n == "gene"), None),
        )

        if ct_orig is None:
            canonical, display, _conf = _standardize_cell_type(sheet)
            confirm: dict[str, list[str]] = {}
            add: dict[str, list[str]] = {}
            for _, row in df.iterrows():
                genes_str = str(row[marker_orig]).strip()
                if not genes_str or genes_str == "nan":
                    continue
                genes = [g.strip().upper() for g in genes_str.split(",") if g.strip()]
                for g in genes[: self.max_genes_per_type]:
                    add.setdefault(g, []).append(pmid)
            if confirm or add:
                markers[canonical] = {"confirm": confirm, "add": add}
            return

        for ct_name in df[ct_orig].unique():
            ct_df = df[df[ct_orig] == ct_name]
            canonical, display, _conf = _standardize_cell_type(str(ct_name))
            confirm: dict[str, list[str]] = {}
            add: dict[str, list[str]] = {}
            for _, row in ct_df.iterrows():
                genes_str = str(row[marker_orig]).strip()
                if not genes_str or genes_str == "nan":
                    continue
                genes = [g.strip().upper() for g in genes_str.split(",") if g.strip()]
                for g in genes:
                    if len(confirm) + len(add) >= self.max_genes_per_type:
                        break
                    add.setdefault(g, []).append(pmid)
            if confirm or add:
                _merge_markers(markers, canonical, confirm, add)

    def _parse_grouped_sheet(
        self,
        xls: pd.ExcelFile,
        sheet: str,
        pmid: str,
        markers: dict[str, dict[str, Any]],
    ) -> None:
        """Parse a grouped sheet with 'group' column (Li S7A variant)."""
        try:
            df = _read_sheet(xls, sheet)
        except Exception:
            return

        norm = _normalize_columns(df)
        if "group" not in norm.values() or "gene" not in norm.values():
            return

        group_col = next(c for c, n in norm.items() if n == "group")
        gene_col = next(c for c, n in norm.items() if n == "gene")
        logfc_col = next((c for c, n in norm.items() if n == "logFC"), None)
        score_col = next((c for c, n in norm.items() if n == "score"), None)

        for group_name in df[group_col].unique():
            group_df = df[df[group_col] == group_name].copy()
            gs = str(group_name).lower().strip()
            if gs in {"group", "names", ""}:
                continue

            canonical, display, _conf = _standardize_cell_type(str(group_name))

            if score_col and score_col in df.columns:
                group_df = group_df.sort_values(score_col, ascending=False)
            elif logfc_col and logfc_col in df.columns:
                group_df = group_df.sort_values(logfc_col, ascending=False)

            confirm: dict[str, list[str]] = {}
            add: dict[str, list[str]] = {}

            for _, row in group_df.iterrows():
                gene = str(row[gene_col]).strip().upper()
                if not gene or gene == "NAN":
                    continue
                if len(confirm) + len(add) >= self.max_genes_per_type:
                    break

                logfc = (
                    float(row[logfc_col])
                    if logfc_col and logfc_col in df.columns and pd.notna(row[logfc_col])
                    else None
                )
                score = (
                    float(row[score_col])
                    if score_col and score_col in df.columns and pd.notna(row[score_col])
                    else None
                )
                if (logfc is not None and logfc > self.confirm_logfc_threshold) or (
                    score is not None and score > self.confirm_score_threshold
                ):
                    confirm.setdefault(gene, []).append(pmid)
                elif (logfc is not None and logfc > self.add_logfc_threshold) or (
                    score is not None and score > self.add_score_threshold
                ):
                    add.setdefault(gene, []).append(pmid)

            if confirm or add:
                _merge_markers(markers, canonical, confirm, add)

    def _parse_unknown(
        self,
        xls: pd.ExcelFile,
        pmid: str,
    ) -> dict[str, Any]:
        """Fallback parser for unknown formats (Peng MarkerGenes, etc.)."""
        markers: dict[str, dict[str, Any]] = {}

        for sheet in xls.sheet_names:
            try:
                df = _read_sheet(xls, sheet)
            except Exception:
                continue

            cols_lower = [str(c).lower().strip() for c in df.columns]
            has_sig = any("signature" in c for c in cols_lower)
            has_ct = any("cell class" in c for c in cols_lower)
            has_abbrev = any("abbreviation" in c for c in cols_lower)

            if has_sig and (has_ct or has_abbrev):
                sig_col = next(c for c in df.columns if "signature" in str(c).lower())
                ct_col = next(
                    (c for c in df.columns if "cell class" in str(c).lower()),
                    next((c for c in df.columns if "abbreviation" in str(c).lower()), None),
                )
                if ct_col is None:
                    continue

                for _, row in df.iterrows():
                    ct_name = str(row[ct_col]).strip()
                    genes_str = str(row[sig_col]).strip()
                    if not ct_name or not genes_str or genes_str == "nan":
                        continue

                    canonical, display, _conf = _standardize_cell_type(ct_name)
                    genes = [g.strip().upper() for g in genes_str.split(",") if g.strip()]
                    if not genes:
                        continue

                    conf = {g: [pmid] for g in genes[:5]}
                    ad = {g: [pmid] for g in genes[5 : self.max_genes_per_type]}

                    if canonical in markers:
                        existing = markers[canonical]
                        for g, ps in conf.items():
                            existing["confirm"].setdefault(g, []).extend(ps)
                        for g, ps in ad.items():
                            existing["add"].setdefault(g, []).extend(ps)
                    else:
                        markers[canonical] = {"confirm": conf, "add": ad}

        return markers


# ── YAML output helpers ──────────────────────────────────────────────────────


def _clean_meta(meta: dict[str, Any]) -> dict[str, Any]:
    """Remove empty values from source_meta for clean YAML output."""
    required = [
        "id",
        "short_name",
        "pmid",
        "journal",
        "year",
        "species",
        "tissue",
        "class",
        "order",
    ]
    result: dict[str, Any] = {}
    for key in required:
        if key in meta:
            result[key] = meta[key]
    for key in ["regions", "n_cells", "n_subtypes", "n_groups"]:
        if key in meta and meta[key]:
            result[key] = meta[key]
    return result


def _clean_markers(markers: dict[str, Any]) -> dict[str, Any]:
    """Remove empty sub-dicts from markers for clean YAML output."""
    result: dict[str, Any] = {}
    for ct, entry in markers.items():
        cleaned: dict[str, Any] = {}
        for key in ("confirm", "add", "refine"):
            if key in entry and entry[key]:
                cleaned[key] = entry[key]
        if cleaned:
            result[ct] = cleaned
    return result


def _merge_markers(
    markers: dict[str, dict[str, Any]],
    canonical: str,
    confirm: dict[str, list[str]],
    add: dict[str, list[str]],
) -> None:
    """Merge confirm/add markers into the markers dict."""
    if canonical in markers:
        existing = markers[canonical]
        for g, ps in confirm.items():
            existing["confirm"].setdefault(g, []).extend(ps)
        for g, ps in add.items():
            existing["add"].setdefault(g, []).extend(ps)
    else:
        markers[canonical] = {"confirm": confirm, "add": add}
