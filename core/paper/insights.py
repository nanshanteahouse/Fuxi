#!/usr/bin/env python3
"""
paper_insights.py — AI-assisted extraction of structured paper insights

Reads a paper's PDF-to-text markdown file and uses LLM (via core.ai_caller)
to generate a structured insights.yaml for reproduction tracking.

Usage:
    python core/paper_insights.py <paper.md> [--output OUTPUT] [--force] [--methodology]
"""

import argparse
import json
import logging
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# Ensure repo root is on sys.path for standalone CLI usage
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from core.ai.caller import ai_query  # noqa: E402
from core.ai.prompts import (  # noqa: E402
    PAPER_FIGURE_SYSTEM_PROMPT,
    PAPER_FIGURE_USER_TEMPLATE,
    PAPER_META_SYSTEM_PROMPT,
    PAPER_META_USER_TEMPLATE,
    PAPER_METHODOLOGY_SYSTEM_PROMPT,
    PAPER_METHODOLOGY_USER_TEMPLATE,
    PAPER_METHODS_SYSTEM_PROMPT,
    PAPER_METHODS_USER_TEMPLATE,
)
from core.paper.converter import (  # noqa: E402
    MarkdownSource,
    PaperSource,
    PmcXmlSource,
    PubmedSource,
    Pymupdf4llmSource,
)

logger = logging.getLogger(__name__)

# -- Section / figure regex patterns --------------------------------------------
_SECTION_RE = re.compile(
    r"^(?:#+\s*)?(SUMMARY|Abstract|Introduction|Results|Discussion|Methods|Materials\s*(?:and|&)\s*Methods|Experimental\s*Procedures)(?!(?-i:[a-z]))",
    re.MULTILINE | re.IGNORECASE,
)
_FIGURE_RE = re.compile(r"(?:Figure|Fig\.?)\s+\d+[a-z]?", re.IGNORECASE)

_METHOD_KEYWORDS = [
    "Seurat",
    "Scanpy",
    "Harmony",
    "UMAP",
    "t-SNE",
    "tsne",
    "Wilcoxon",
    "Mann-Whitney",
    "MAST",
    "DESeq2",
    "edgeR",
    "Monocle",
    "Slingshot",
    "Velocity",
    "scVelo",
    "CellChat",
    "NicheNet",
    "LIANA",
    "CellPhoneDB",
    "SCENIC",
    "pySCENIC",
    "AUCell",
    "ROGUE",
    "scran",
    "scran.js",
    "SoupX",
    "DoubletFinder",
    "Scrublet",
    "MAGIC",
    "SCTransform",
    "SCVI",
    "scGPT",
    "BayesSpace",
    "SPOTlight",
    "CARD",
    "RCTD",
]


@dataclass
class _LLMConfig:
    """Minimal LLM config for standalone CLI use."""

    model: str = "deepseek-v4-flash"
    api_base: str = "https://api.deepseek.com/v1"
    api_key: str = ""
    max_tokens: int = 16384
    temperature: float = 0.3
    ai_cache_responses: bool = True
    thinking_enabled: bool = False


def _parse_filename_meta(md_path: str) -> dict:
    """Parse year/author/journal from filename pattern {year}_{author}_{journal}_{title}.

    .. deprecated:: paper-module
       Use ``PaperSource.get_metadata()`` instead.
    """
    stem = Path(md_path).stem
    parts = stem.split("_")
    meta: dict = {
        "filename": Path(md_path).name,
        "stem": stem,
        "year": None,
        "first_author": None,
        "journal": None,
        "title": stem,
    }
    if len(parts) >= 3 and parts[0].isdigit() and len(parts[0]) in (2, 4):
        meta["year"] = parts[0]
        meta["first_author"] = parts[1]
        meta["journal"] = parts[2]
        meta["title"] = "_".join(parts[3:]) if len(parts) > 3 else parts[2]
    return meta


def _safe_json_parse(raw: str, fallback_label: str = "unknown") -> dict:
    """Parse JSON from LLM response, stripping markdown fences. Returns {} on failure."""
    text = raw.strip()
    if text.startswith("```"):
        idx = text.find("\n")
        if idx != -1:
            text = text[idx + 1 :]
        text = text.removesuffix("```").strip()
    if text.startswith("`") and text.endswith("`"):
        text = text[1:-1].strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        logger.warning("Failed to parse JSON for %s: %.200s", fallback_label, text)
        return {}


def _yaml_dump(data: Any, indent: int = 0) -> list[str]:
    """Serialize Python value to YAML lines (no PyYAML dependency)."""
    pfx = "  " * indent
    if data is None:
        return [f"{pfx}null"]
    if isinstance(data, bool):
        return [f"{pfx}{str(data).lower()}"]
    if isinstance(data, (int, float)):
        return [f"{pfx}{data}"]
    if isinstance(data, str):
        if "\n" in data:
            return [f"{pfx}|"] + [f"{pfx}  {line}" for line in data.split("\n")]
        escaped = data.replace("\\", "\\\\").replace('"', '\\"')
        return [f'{pfx}"{escaped}"']
    if isinstance(data, dict):
        if not data:
            return [f"{pfx}{{}}"]
        lines: list[str] = []
        for k, v in data.items():
            ks = str(k)
            sub = _yaml_dump(v, indent + 1)
            if isinstance(v, dict):
                lines.append(f"{pfx}{ks}:")
                lines.extend(sub)
            elif len(sub) == 1 and not sub[0].startswith("  " * (indent + 1) + "-"):
                lines.append(f"{pfx}{ks}: {sub[0].strip()}")
            else:
                lines.append(f"{pfx}{ks}:")
                lines.extend(sub)
        return lines
    if isinstance(data, (list, tuple)):
        if not data:
            return [f"{pfx}[]"]
        lines = []
        for item in data:
            sub = _yaml_dump(item, indent + 1)
            if len(sub) == 1:
                lines.append(f"{pfx}- {sub[0].strip()}")
            elif isinstance(item, dict) and len(sub) > 1:
                first = sub[0].lstrip()
                remaining = sub[1:]
                lines.append(f"{pfx}- {first}")
                lines.extend(remaining)
            else:
                lines.append(f"{pfx}-")
                lines.extend(sub)
        return lines
    return [f"{pfx}{data}"]


def _build_cfg_from_env() -> _LLMConfig:
    """Build LLM config from environment variables."""
    return _LLMConfig(
        api_key=os.environ.get("LLM_API_KEY", ""),
        api_base=os.environ.get("LLM_API_BASE", "https://api.deepseek.com/v1"),
        model=os.environ.get("LLM_MODEL", "deepseek-v4-flash"),
    )


def _extract_geo_ids(text: str) -> list[str]:
    """Extract GEO accession IDs (GSE\\d{4,8}) from text via regex."""
    return re.findall(r"GSE\d{4,8}", text)


def _extract_data_access(meta: dict | None, methods: dict | None, full_text: str = "") -> dict:
    """Extract geo_ids/sra_ids from meta.data_access or regex-fallback from data_notes.

    Falls back to regex scan of full_text if both primary path and data_notes yield empty.
    """
    result: dict[str, list] = {"geo_ids": [], "sra_ids": []}

    # Primary path: meta has data_access sub-object
    if meta:
        da = meta.get("data_access", {})
        if isinstance(da, dict):
            result["geo_ids"] = list(da.get("geo_ids", []))
            result["sra_ids"] = list(da.get("sra_ids", []))

    # Fallback: scan meta and methods data_notes
    if not result["geo_ids"] and not result["sra_ids"]:
        notes = []
        if meta:
            notes.extend(meta.get("data_notes", []))
        if methods:
            notes.extend(methods.get("data_notes", []))
        for note in notes:
            if isinstance(note, str):
                geo = re.findall(r"\bGSE\d{4,}\b", note)
                sra = re.findall(r"\bSRP\d{4,}\b", note)
                result["geo_ids"].extend(geo)
                result["sra_ids"].extend(sra)

    # Final fallback: regex scan of full text
    if not result["geo_ids"] and full_text:
        result["geo_ids"] = _extract_geo_ids(full_text)

    # Deduplicate while preserving order
    result["geo_ids"] = list(dict.fromkeys(result["geo_ids"]))
    result["sra_ids"] = list(dict.fromkeys(result["sra_ids"]))
    return result


def _extract_key_methods(text: str) -> list[str]:
    """Extract known method keywords from text via case-insensitive word-boundary matching.

    Uses _METHOD_KEYWORDS list. Returns deduplicated list preserving order of appearance.
    """
    found: list[str] = []
    seen: set[str] = set()
    for kw in _METHOD_KEYWORDS:
        if re.search(r"\b" + re.escape(kw) + r"\b", text, re.IGNORECASE):
            if kw not in seen:
                seen.add(kw)
                found.append(kw)
    return found


def _extract_methods_summary(methods_data: dict | None, full_text: str = "") -> dict:
    """Extract methods summary from methods LLM extraction output.

    Falls back to regex keyword scan of full_text when LLM yields empty key_methods.
    """
    result = {
        "key_methods": [],
        "software_versions": {},
        "reference_genome": None,
        "sequencing_platforms": [],
    }
    if methods_data:
        result["key_methods"] = list(methods_data.get("key_methods", []))
        result["software_versions"] = dict(methods_data.get("software_versions", {}))
        result["reference_genome"] = methods_data.get("reference_genome")
        result["sequencing_platforms"] = list(methods_data.get("sequencing_platforms", []))

    # Fallback: regex keyword scan when LLM yielded empty key_methods
    if not result["key_methods"] and full_text:
        result["key_methods"] = _extract_key_methods(full_text)

    return result


# ---------------------------------------------------------------------------
# Methodology-safe truncation: clip methods to a max token budget for LLM
# ---------------------------------------------------------------------------
_METHODOLOGY_METHODS_MAX_CHARS = 8000


def _resolve_doi_to_pmid(doi: str) -> str | None:
    """Resolve a DOI to PubMed ID via NCBI E-utilities esearch.

    Returns PMID string or None on failure. Includes rate limiting (0.35s)
    and exponential backoff on 429/503.
    """
    if not doi or not doi.strip():
        return None
    doi = doi.strip()
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=pubmed&term={doi}[doi]&retmode=json"
    )
    time.sleep(0.35)  # NCBI rate limit: ~3 req/s
    for attempt in range(3):
        try:
            req = Request(url)
            req.add_header("User-Agent", "Fuxi/1.0 (paper-insights; academic use)")
            with urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                idlist = data.get("esearchresult", {}).get("idlist", [])
                if idlist:
                    return str(idlist[0])
                return None
        except HTTPError as e:
            if e.code in (429, 503):
                time.sleep(0.5 * (2**attempt))
                continue
            return None
        except (URLError, OSError, json.JSONDecodeError):
            return None
    return None


class PaperInsights:
    """Extract structured insights from paper markdown using AI prompts."""

    @staticmethod
    def split_sections(md_text: str) -> dict[str, str]:
        """Split markdown into sections by header.

        Returns {abstract, introduction, results, discussion, methods} with
        lowercase keys. Missing sections get empty string.
        """
        sections: dict[str, str] = {
            k: "" for k in ("abstract", "introduction", "results", "discussion", "methods")
        }
        matches = list(_SECTION_RE.finditer(md_text))
        if not matches:
            sections["results"] = md_text.strip()
            return sections
        _section_key_map = {
            "summary": "abstract",
            "abstract": "abstract",
            "introduction": "introduction",
            "results": "results",
            "discussion": "discussion",
            "methods": "methods",
            "experimental procedures": "methods",
            "experimentalprocedures": "methods",
            "materials and methods": "methods",
            "materials & methods": "methods",
        }
        for i, m in enumerate(matches):
            raw_key = m.group(1).lower().strip()
            key = _section_key_map.get(raw_key, "unknown")
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(md_text)
            sections[key] = md_text[start:end].strip()
        # Capture pre-heading text as abstract if no Abstract heading was matched
        if not sections["abstract"] and matches and matches[0].start() > 0:
            pre_text = md_text[: matches[0].start()].strip()
            if pre_text:
                sections["abstract"] = pre_text
        return sections

    @staticmethod
    def extract_figure_blocks(results_text: str) -> list[str]:
        """Split Results text into figure blocks by Figure/Fig references."""
        if not results_text.strip():
            return []
        blocks = _FIGURE_RE.split(results_text)
        refs = _FIGURE_RE.findall(results_text)
        result: list[str] = []
        for i, block in enumerate(blocks[1:], start=1):
            header = refs[i - 1] if i - 1 < len(refs) else f"Figure_{i}"
            cleaned = block.strip().lstrip(".: ")
            combined = f"{header}. {cleaned}"
            if combined.strip():
                result.append(combined)
        if not result and results_text.strip():
            result.append(results_text.strip())
        return result

    @staticmethod
    def extract_metadata(abstract_text: str, cfg, full_text: str = "") -> dict:
        """Extract experimental design, key findings, data notes via LLM.

        Uses both the abstract and the first ~2000 chars of the full text
        (which usually contains title/authors/journal) for bibliographic fields.

        Returns parsed dict or {} on failure. Skips gracefully if abstract is empty.
        """
        if not abstract_text.strip():
            logger.info("Abstract empty -- skipping metadata extraction")
            return {}
        full_beginning = full_text[:2000].strip() if full_text else ""
        raw = ai_query(
            PAPER_META_SYSTEM_PROMPT,
            PAPER_META_USER_TEMPLATE.format(
                full_text_beginning=full_beginning or "(not available)",
                abstract_text=abstract_text,
            ),
            cfg,
            expect_json=True,
        )
        if raw is None:
            return {}
        return _safe_json_parse(raw, "abstract metadata")

    @staticmethod
    def extract_figure(figure_block: str, cfg) -> dict:
        """Extract structured info from figure legend via LLM. Returns {} on failure."""
        if not figure_block.strip():
            return {}
        raw = ai_query(
            PAPER_FIGURE_SYSTEM_PROMPT,
            PAPER_FIGURE_USER_TEMPLATE.format(figure_text=figure_block),
            cfg,
            expect_json=True,
        )
        if raw is None:
            return {}
        return _safe_json_parse(raw, "figure")

    @staticmethod
    def extract_methods(methods_text: str, cfg) -> dict:
        """Extract methods, software, versions via LLM. Returns {} on failure."""
        if not methods_text.strip():
            logger.info("Methods empty -- skipping methods extraction")
            return {}
        raw = ai_query(
            PAPER_METHODS_SYSTEM_PROMPT,
            PAPER_METHODS_USER_TEMPLATE.format(methods_text=methods_text),
            cfg,
            expect_json=True,
        )
        if raw is None:
            return {}
        return _safe_json_parse(raw, "methods")

    @staticmethod
    def extract_methodology(
        abstract: str, methods_text: str, key_findings: list[str], cfg
    ) -> dict:
        """Extract methodology patterns (5 dimensions) from paper sections via LLM.

        Returns dict with 'methodology_patterns' key, or {} on failure.
        Skips gracefully if both abstract and methods are empty.
        """
        if not abstract.strip() and not methods_text.strip():
            logger.info("Abstract + methods both empty -- skipping methodology extraction")
            return {}

        kf_text = (
            "\n".join(f"- {kf}" for kf in key_findings)
            if key_findings
            else "(no key findings available)"
        )

        user_prompt = PAPER_METHODOLOGY_USER_TEMPLATE.format(
            abstract=abstract.strip() or "(abstract not available)",
            methods=methods_text.strip()[:_METHODOLOGY_METHODS_MAX_CHARS]
            or "(methods not available)",
            key_findings=kf_text,
        )
        raw = ai_query(PAPER_METHODOLOGY_SYSTEM_PROMPT, user_prompt, cfg, expect_json=True)
        if raw is None:
            return {}
        return _safe_json_parse(raw, "methodology_patterns")

    @staticmethod
    def merge_to_insights(
        meta: dict,
        figures: list[dict],
        methods_data: dict,
        paper_meta: dict,
        full_text: str = "",
        methods_text: str = "",
        methodology: dict | None = None,
    ) -> dict:
        """Merge LLM-extracted metadata, figures, methods, and methodology into a single insights dict.

        Args:
            meta: Metadata dict from extract_metadata().
            figures: List of figure dicts from extract_figure().
            methods_data: Methods dict from extract_methods().
            paper_meta: Paper metadata from PaperSource.get_metadata().
            full_text: Full markdown text for regex fallback in data_access extraction.
            methods_text: Methods section text for regex fallback in methods extraction.
            methodology: Optional methodology_patterns dict from extract_methodology().
        """
        # Deduplicate figures by id while preserving order
        seen: set[str] = set()
        unique_figs: list[dict] = []
        for fig in figures:
            fid = fig.get("id", "")
            if fid not in seen:
                seen.add(fid)
                unique_figs.append(fig)
        unique_figs.sort(key=lambda x: x.get("id", ""))

        meta_notes = list(meta.get("data_notes", [])) if meta else []
        methods_notes = list(methods_data.get("data_notes", [])) if methods_data else []

        # Build reproduction_status with tracking + computed fields
        fig_count = len(unique_figs)
        repro_count = sum(1 for f in unique_figs if f.get("reproducible"))

        result = {
            "paper_meta": paper_meta,
            "experimental_design": meta.get("experimental_design", {}) if meta else {},
            "key_findings": list(meta.get("key_findings", [])) if meta else [],
            "data_access": _extract_data_access(meta, methods_data, full_text=full_text),
            "methods": _extract_methods_summary(methods_data, full_text=methods_text),
            "figures": unique_figs,
            "data_notes": meta_notes + [n for n in methods_notes if n not in meta_notes],
            "reproduction_status": {
                "pipeline_run": "not_started",
                "overall_match": None,
                "total_figures": fig_count,
                "reproducible_count": repro_count,
                "verified_figures": [],
                "notes": "",
            },
        }

        # ── Bibliographic backfill ──
        # The AI bibliographic extractor has access to the full text beginning
        # (title, author list, journal, etc.) and is generally more reliable
        # than filename-derived metadata. Override filename-derived values
        # when the AI provides a substantively better value.
        biblio = (meta or {}).get("bibliographic", {}) or {}
        if biblio:

            def _is_good(v):
                return (
                    v is not None
                    and str(v).strip() != ""
                    and not str(v).startswith("pymupdf4llm_")
                )

            def _is_bad(v):
                return v is None or str(v).strip() == "" or str(v).startswith("pymupdf4llm_")

            for _key, _bkey in [
                ("pmid", "pmid"),
                ("title", "title"),
                ("journal", "journal"),
                ("year", "year"),
                ("first_author", "first_author"),
                ("doi", "doi"),
            ]:
                ai_v = biblio.get(_bkey)
                pm_v = paper_meta.get(_key)
                if not _is_good(ai_v):
                    continue  # AI has nothing useful for this field
                if _is_bad(pm_v):
                    paper_meta[_key] = ai_v
                elif _key in ("title", "journal"):
                    if len(str(ai_v)) > len(str(pm_v)):
                        paper_meta[_key] = ai_v
                else:
                    paper_meta[_key] = ai_v

        # ── NCBI PMID resolution ──
        # When AI couldn't extract PMID from the text but we have a DOI,
        # resolve it via NCBI E-utilities esearch.
        if not paper_meta.get("pmid") and paper_meta.get("doi"):
            resolved = _resolve_doi_to_pmid(str(paper_meta["doi"]))
            if resolved:
                logger.info("Resolved DOI %s → PMID %s", paper_meta["doi"], resolved)
                paper_meta["pmid"] = resolved

        if methodology:
            result["methodology_patterns"] = methodology.get("methodology_patterns", methodology)

        return result

    def run(
        self,
        source: Union[PaperSource, str],
        cfg,
        output_path: str | None = None,
        force: bool = False,
        extract_methodology: bool = False,
    ) -> str | None:
        """Full pipeline: read, split, extract, merge, write.

        Returns path to written file, "SKIPPED" if output exists and not force,
        or None if the source could not be read.

        Auto-creates subdirectory for output if needed.

        Args:
            extract_methodology: If True, also run methodology pattern extraction (5 dimensions).
        """
        if isinstance(source, str):
            source = MarkdownSource(source)
        if output_path:
            out_path = Path(output_path)
        else:
            out_path = Path("projects/papers") / source.get_paper_name() / "insights.yaml"

        if out_path.exists() and not force:
            logger.info("Output exists at %s, skipping (use --force to overwrite)", out_path)
            return "SKIPPED"

        try:
            md_text = source.get_text()
        except (RuntimeError, HTTPError, URLError) as exc:
            logger.error("Failed to read source: %s", exc)
            return None
        sections = self.split_sections(md_text)
        logger.info("Found sections: %s", [k for k, v in sections.items() if v])

        meta: dict = {}
        if sections.get("abstract"):
            logger.info("Extracting metadata from abstract...")
            meta = self.extract_metadata(sections["abstract"], cfg, full_text=md_text)

        figures: list[dict] = []
        if sections.get("results"):
            blocks = self.extract_figure_blocks(sections["results"])
            logger.info("Found %d figure blocks", len(blocks))
            for i, block in enumerate(blocks):
                sys.stderr.write(f"\r  Figure {i + 1}/{len(blocks)}...")
                sys.stderr.flush()
                fig = self.extract_figure(block, cfg)
                if fig:
                    figures.append(fig)
            sys.stderr.write("\n")

        methods: dict = {}
        if sections.get("methods"):
            logger.info("Extracting methods...")
            methods = self.extract_methods(sections["methods"], cfg)

        # Optional methodology pattern extraction
        methodology: dict | None = None
        if extract_methodology:
            logger.info("Extracting methodology patterns...")
            key_findings = list(meta.get("key_findings", [])) if meta else []
            abstract_for_meth = sections.get("abstract", "") or md_text[:4000]
            methods_for_meth = sections.get("methods", "") or md_text
            methodology = self.extract_methodology(
                abstract_for_meth,
                methods_for_meth,
                key_findings,
                cfg,
            )

        insights = self.merge_to_insights(
            meta,
            figures,
            methods_data=methods,
            paper_meta=source.get_metadata(),
            full_text=md_text,
            methods_text=sections.get("methods", ""),
            methodology=methodology,
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(_yaml_dump(insights)) + "\n", encoding="utf-8")
        logger.info("Wrote insights to %s", out_path)
        return str(out_path)


def _resolve_source(args) -> PaperSource:
    """Determine PaperSource from CLI arguments."""
    if args.source == "pmc":
        if not (args.pmid or args.doi or args.xml):
            raise ValueError("--source pmc requires --pmid, --doi, or --xml")
        return PmcXmlSource(pmid=args.pmid, doi=args.doi, xml_path=args.xml)

    elif args.source == "pdf":
        if not args.pdf:
            raise ValueError("--source pdf requires --pdf <path>")
        return Pymupdf4llmSource(args.pdf, pmid=args.pmid)

    elif args.source == "md":
        if not args.positional:
            raise ValueError("--source md requires a positional .md file argument")
        return MarkdownSource(args.positional)

    # --source auto (default): PMC -> PubMed -> PDF -> MD fallback
    pubmed_pmid = args.pmid or (_resolve_doi_to_pmid(args.doi) if args.doi else None)

    if args.pmid or args.doi or args.xml:
        try:
            return PmcXmlSource(pmid=args.pmid, doi=args.doi, xml_path=args.xml)
        except (RuntimeError, ValueError, HTTPError, URLError) as e:
            logger.warning("PMC source failed: %s. Trying PubMed fallback...", e)
            if pubmed_pmid:
                try:
                    return PubmedSource(pubmed_pmid)
                except (RuntimeError, ValueError, HTTPError, URLError) as e2:
                    logger.warning("PubMed fallback also failed: %s", e2)

    if args.pdf:
        try:
            return Pymupdf4llmSource(args.pdf, pmid=pubmed_pmid)
        except ImportError:
            logger.warning("pymupdf4llm not installed, skipping PDF fallback")

    if args.positional:
        return MarkdownSource(args.positional)

    raise ValueError(
        "No valid source found. Provide --pmid/--doi/--xml/--pdf or a positional .md file"
    )


def main() -> None:
    """CLI entry point for paper_insights."""
    parser = argparse.ArgumentParser(description="Extract structured paper insights using AI.")
    parser.add_argument(
        "positional",
        nargs="?",
        type=str,
        help="Paper markdown file (optional if --pmid/--doi/--xml/--pdf given)",
    )
    parser.add_argument("--pmid", type=str, default=None, help="PubMed ID")
    parser.add_argument("--doi", type=str, default=None, help="DOI")
    parser.add_argument("--xml", type=str, default=None, help="Local PMC XML file path")
    parser.add_argument("--pdf", type=str, default=None, help="PDF file path")
    parser.add_argument(
        "--source",
        type=str,
        default="auto",
        choices=["auto", "pmc", "pdf", "md"],
        help="Source selection strategy",
    )
    parser.add_argument("--output", "-o", type=str, default=None, help="Output path")
    parser.add_argument(
        "--force", "-f", action="store_true", default=False, help="Overwrite existing output"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=False, help="Enable debug logging"
    )
    parser.add_argument(
        "--methodology",
        "-m",
        action="store_true",
        default=False,
        help="Also extract methodology patterns (5 dimensions: archetype, strategy, narrative, toolbox, contribution)",
    )
    args = parser.parse_args()

    if args.pmid and args.xml:
        parser.error("--pmid and --xml are mutually exclusive")

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s [%(name)s] %(message)s",
        stream=sys.stderr,
    )
    logging.getLogger("httpx").setLevel(logging.WARNING)

    cfg = _build_cfg_from_env()
    logger.info("Using env-based LLM config")

    source = _resolve_source(args)
    result = PaperInsights().run(
        source=source,
        cfg=cfg,
        output_path=args.output,
        force=args.force,
        extract_methodology=args.methodology,
    )
    if result is None:
        sys.exit(1)
    print("SKIPPED" if result == "SKIPPED" else f"Done: {result}")


if __name__ == "__main__":
    main()
