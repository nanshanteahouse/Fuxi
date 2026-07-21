"""
Paper source conversion utilities for Fuxi (伏羲).

Handles conversion from various paper sources (PMC XML, Markdown, PDF)
into a unified internal format for downstream processing.

Known limitation: ``clean_text`` applies regex-based fixes that may split
biology mixed-case terms like ``snRNA-seq`` → ``sn RNA-seq``. This is
acceptable for P0 — the benefit of fixing ~80% of concatenation artifacts
outweighs ~5% false positives for LLM processing context.
"""

import re


def clean_text(text: str) -> str:
    """Post-process text from PDF-to-Markdown conversion.

    Fixes word concatenation artifacts produced by ``markitdown`` and similar
    tools that omit spaces between some words (e.g. ``theelderly.However``
    becomes ``the elderly. However``).

    This post-processor is primarily designed for :class:`Pymupdf4llmSource` output but is
    safe to apply universally --- it is a no-op on already-clean text (PMC XML,
    existing ``.md`` files).

    Parameters
    ----------
    text : str
        Raw text to clean.

    Returns
    -------
    str
        Cleaned text with whitespace normalised.
    """
    # Step 1: Watermark/boilerplate line removal (line-context-aware)
    # Remove standalone lines matching "Author Manuscript", "HHS Public Access",
    # Only remove when they appear as standalone lines (^\s*pattern\s*$)
    text = re.sub(r"^ *Author Manuscript *$", "", text, flags=re.MULTILINE)
    text = re.sub(r"^ *HHS Public Access *$", "", text, flags=re.MULTILINE)
    text = re.sub(r"^ *Author *$", "", text, flags=re.MULTILINE)
    text = re.sub(r"^ *Manuscript *$", "", text, flags=re.MULTILINE)

    # Step 2: Character-spaced text merge (line-context-aware)
    # Lines consisting entirely of single-letter sequences separated by whitespace
    # e.g. "A u t h o r   M a n u s c r i p t" → characters merged, spaces removed
    text = re.sub(
        r"^(?:[A-Za-z]\s+){2,}[A-Za-z]\s*$",
        lambda m: m.group(0).replace(" ", ""),
        text,
        flags=re.MULTILINE,
    )

    # Step 3: Line-number stripping at line starts
    # Strip leading 1-3 digit numbers followed by space
    text = re.sub(r"^\d{1,3}\s+", "", text, flags=re.MULTILINE)

    # Step 4: Garbage character sequence suppression
    # Remove runs of 8+ identical repeated characters (bioRxiv header artifacts)
    text = re.sub(r"(.)\1{8,}", "", text)

    # P0: Insert space between lowercase letter and following uppercase letter
    # (catches concatenated words like ``ofAMD`` → ``of AMD``)
    text = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", text)

    # P0: Insert space after period followed by uppercase letter
    # (catches sentence-boundary concatenation like ``elderly.However``)
    text = re.sub(r"(?<=\.)(?=[A-Z])", " ", text)

    # Collapse multiple spaces to single space
    text = re.sub(r" +", " ", text)

    # Collapse 3+ consecutive newlines to exactly 2 (paragraph separator)
    text = re.sub(r"\n{3,}", "\n\n", text)

    return text.strip()


import json  # noqa: E402
import logging  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402
import xml.etree.ElementTree as ET  # noqa: E402
from abc import ABC, abstractmethod  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Optional  # noqa: E402
from urllib.error import HTTPError, URLError  # noqa: E402
from urllib.request import Request, urlopen  # noqa: E402

logger = logging.getLogger(__name__)


# ── Constants ───────────────────────────────────────────────────────────────────

JATS_NS = "http://www.ncbi.nlm.nih.gov/JATS1"

_SECTION_KEY_MAP: dict[str, str] = {
    "abstract": "abstract",
    "introduction": "introduction",
    "background": "introduction",
    "results": "results",
    "discussion": "discussion",
    "methods": "methods",
    "materials and methods": "methods",
    "materials & methods": "methods",
    "experimental procedures": "methods",
}

# Matches markdown or plain-text section headings (fallback split)
_FALLBACK_SECTION_RE = re.compile(
    r"^(?:#+\s*)?(SUMMARY|Abstract|Introduction|Results|Discussion|Methods|"
    r"Experimental\s*Procedures)\b",
    re.MULTILINE | re.IGNORECASE,
)


# ── Helpers ────────────────────────────────────────────────────────────────────────


def _strip_xml_namespaces(raw: str) -> str:
    """Remove xmlns declarations AND strip prefix from element/attribute names.
    This avoids "unbound prefix" errors when DTD is stripped but namespace
    declarations were the only source of prefix definitions."""
    # 1. Remove xmlns declarations from attributes
    raw = re.sub(r'\s+xmlns(?:\:\w+)?="[^"]*"', "", raw)
    # 2. Strip namespace prefixes from element tag names (opening AND closing)
    #    e.g. <ali:license_ref → <license_ref,  </ali:license_ref → </license_ref
    raw = re.sub(r"(</?)[\w-]+:(?=[\w-])", r"\1", raw)
    # 3. Strip namespace prefixes from attribute names
    #    e.g. xlink:href="..." → href="..."
    raw = re.sub(r"(?<=\s)[\w-]+:(?=[\w-])", "", raw)
    return raw


def _elem_text(el: Optional[ET.Element]) -> str:
    """Return all text within *el*, stripped, or '' if *el* is None."""
    if el is None:
        return ""
    return "".join(el.itertext()).strip()


def _slugify(value: str, max_len: int = 60) -> str:
    """Replace non-alphanumeric chars (except hyphen) with underscore, capped at max_len."""
    value = re.sub(r"[^\w\s-]", "_", value)
    value = re.sub(r"[-\s]+", "_", value)
    value = value.strip("_")
    return value[:max_len].rstrip("_")


# ═══════════════════════════════════════════════════════════════════════════════════
#  PaperSource — abstract base class
# ═══════════════════════════════════════════════════════════════════════════════════


class PaperSource(ABC):
    """Abstract base for all paper source converters.

    Subclasses must implement all five abstract methods to provide
    a uniform interface for downstream insight extraction.
    """

    @abstractmethod
    def get_sections(self) -> dict[str, str]:
        """Return section text keyed by canonical name.

        Keys include ``abstract``, ``introduction``, ``results``,
        ``discussion``, ``methods``, plus any additional named sections.
        Missing sections are omitted rather than returned as empty strings.
        """
        ...

    @abstractmethod
    def get_figure_blocks(self) -> list[str]:
        """Return list of figure captions.

        Each entry is like ``'Fig 1. Sample size and demographics...'``.
        """
        ...

    @abstractmethod
    def get_metadata(self) -> dict:
        """Return metadata dict with keys: pmid, doi, title, first_author, journal, year.

        Missing fields are ``None`` or empty string.
        """
        ...

    @abstractmethod
    def get_paper_name(self) -> str:
        """Return a filesystem-safe paper name.

        Format: ``{year}_{first_author}_{journal}_{title_slug}``.
        Used as the output directory name for processed paper data.
        """
        ...

    @abstractmethod
    def get_text(self) -> str:
        """Return the full paper text as a single concatenated string.

        Suitable for LLM context injection. Typically joins all section
        texts with double newlines.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════════════════
#  PmcXmlSource — JATS XML from NCBI PMC
# ═══════════════════════════════════════════════════════════════════════════════════


class PmcXmlSource(PaperSource):
    """Parse a JATS XML paper from NCBI PMC.

    Accepts either a ``pmid``, ``doi``, or a local ``xml_path``.
    If a PMID or DOI is given, the XML is fetched live from NCBI E-utilities
    with automatic rate limiting, retry with exponential backoff, and
    local disk caching.

    Parameters
    ----------
    pmid : str or None
        PubMed ID (e.g. ``'31467224'``).
    doi : str or None
        Digital Object Identifier (e.g. ``'10.1038/s41467-019-10874-1'``).
    xml_path : str or None
        Path to a local JATS XML file.

    Raises
    ------
    ValueError
        If none of *pmid*, *doi*, or *xml_path* is provided.
    RuntimeError
        If online resolution or XML parsing fails.
    """

    def __init__(
        self,
        pmid: Optional[str] = None,
        doi: Optional[str] = None,
        xml_path: Optional[str] = None,
    ) -> None:
        if pmid is None and doi is None and xml_path is None:
            raise ValueError("Provide at least one of: pmid, doi, xml_path")

        self._pmid: Optional[str] = pmid
        self._doi: Optional[str] = doi
        self._xml_path: Optional[str] = xml_path

        self._pmcid: Optional[str] = None
        self._raw_xml: Optional[str] = None
        self._root: Optional[ET.Element] = None

        # Lazy-loaded caches
        self._sections: Optional[dict[str, str]] = None
        self._figures: Optional[list[str]] = None
        self._metadata: Optional[dict] = None
        self._paper_name: Optional[str] = None

        self._load()

    # ── Loading ────────────────────────────────────────────────────────────────────

    def _load(self) -> None:
        """Load XML from disk or NCBI, then parse."""
        if self._xml_path:
            self._raw_xml = Path(self._xml_path).read_text(encoding="utf-8")
        else:
            self._fetch_from_ncbi()
        self._parse_xml()
        # Cache fetched XML for offline reuse
        if not self._xml_path and self._raw_xml:
            self._cache_xml()

    def _cache_xml(self) -> None:
        """Save fetched XML to ``projects/papers/{paper_name}/{pmcid}.xml``.

        Non-fatal on failure (disk issues, missing metadata).
        """
        pmcid = self._pmcid
        if not pmcid:
            return
        try:
            name = self.get_paper_name()
            cache_dir = Path("projects") / "papers" / name
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_path = cache_dir / f"{pmcid}.xml"
            if not cache_path.exists():
                cache_path.write_text(self._raw_xml or "", encoding="utf-8")
                logger.info("Cached XML to %s", cache_path)
        except Exception:
            logger.warning("Failed to cache XML", exc_info=True)

    # ── NCBI E-utilities ──────────────────────────────────────────────────────────

    def _fetch_from_ncbi(self) -> None:
        """Resolve PMID/DOI to PMCID, then fetch full XML."""
        pmcid: Optional[str] = None

        if self._pmid:
            pmcid = self._pmid_to_pmcid(self._pmid)
        elif self._doi:
            pmcid = self._doi_to_pmcid(self._doi)

        if not pmcid:
            raise RuntimeError(f"Could not resolve PMCID from pmid={self._pmid} doi={self._doi}")

        self._pmcid = pmcid

        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=pmc&id={pmcid}&retmode=xml"
        )
        self._raw_xml = self._ncbi_fetch(url)

    def _pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """Use NCBI elink to convert a PMID to a PMCID."""
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
            f"?dbfrom=pubmed&db=pmc&linkname=pubmed_pmc&id={pmid}&retmode=json"
        )
        data = json.loads(self._ncbi_fetch(url))
        try:
            for ls in data.get("linksets", []):
                for lsd in ls.get("linksetdbs", []):
                    links = lsd.get("links", [])
                    if links:
                        return str(links[0])
        except (KeyError, IndexError, ValueError):
            pass
        return None

    def _doi_to_pmcid(self, doi: str) -> Optional[str]:
        """Use NCBI esearch to convert a DOI to a PMCID."""
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=pmc&term={doi}[doi]&retmode=json"
        )
        data = json.loads(self._ncbi_fetch(url))
        try:
            idlist = data.get("esearchresult", {}).get("idlist", [])
            if idlist:
                return str(idlist[0])
        except (KeyError, IndexError, ValueError):
            pass
        return None

    def _ncbi_fetch(self, url: str) -> str:
        """Fetch *url* with User-Agent header, rate limiting, and backoff.

        Retries up to 3 times with exponential backoff on HTTP 429/503.
        """
        # Rate limit: NCBI allows ~3 requests/sec without API key
        time.sleep(0.35)

        last_error: Optional[Exception] = None
        for attempt in range(3):
            try:
                req = Request(url)
                req.add_header("User-Agent", "Fuxi/1.0 (paper-converter; academic use)")
                with urlopen(req, timeout=15) as resp:
                    return resp.read().decode("utf-8")
            except HTTPError as e:
                if e.code in (429, 503):
                    wait = 0.5 * (2**attempt)
                    logger.warning(
                        "NCBI HTTP %d for %s, retrying in %.1fs...",
                        e.code,
                        url,
                        wait,
                    )
                    time.sleep(wait)
                    last_error = e
                else:
                    raise
            except (URLError, OSError) as e:
                last_error = e
                if attempt < 2:
                    time.sleep(0.5 * (2**attempt))
                else:
                    raise

        raise RuntimeError(f"NCBI fetch failed after 3 attempts: {last_error}")

    # ── XML Parsing ────────────────────────────────────────────────────────────────

    def _parse_xml(self) -> None:
        """Strip DOCTYPE and namespaces, parse into ElementTree.
        Handles both bare <article> (fixture) and <pmc-articleset><article> (live NCBI)."""
        raw = self._raw_xml or ""
        # Prevent ElementTree from fetching external DTDs
        raw = re.sub(r"<!DOCTYPE[^>]+>", "", raw)
        # Strip namespace declarations for tag-name access
        raw = _strip_xml_namespaces(raw)
        self._root = ET.fromstring(raw)
        # Unwrap pmc-articleset wrapper from live NCBI responses
        if self._root.tag == "pmc-articleset":
            children = list(self._root)
            if children:
                self._root = children[0]

    # ── Section parsing ────────────────────────────────────────────────────────────

    def _parse_sections(self) -> dict[str, str]:
        """Extract sections from ``<body><sec>`` elements, plus abstract from ``<front>``."""
        root = self._root
        if root is None:
            return {}

        sections: dict[str, str] = {}

        # Extract abstract from <front><article-meta><abstract>
        front_abstract = root.find("./front/article-meta/abstract")
        if front_abstract is not None:
            abs_text = "".join(front_abstract.itertext()).strip()
            if abs_text:
                sections["abstract"] = abs_text

        body = root.find("body")
        if body is None:
            return sections

        secs = body.findall("sec")
        if not secs:
            body_text = "".join(body.itertext()).strip()
            if body_text:
                sections = self._fallback_split(body)
            return sections

        for sec in secs:
            title_el = sec.find("title")
            title = _elem_text(title_el) if title_el is not None else ""
            text = "".join(sec.itertext()).strip()
            # itertext() includes the title text at the start — separate it
            if title and text.startswith(title):
                body = text[len(title) :].lstrip()
                text = f"{title}\n{body}" if body else title

            key = self._section_key(title)
            if key in sections:
                sections[key] += "\n\n" + text
            else:
                sections[key] = text

        return sections

    @staticmethod
    def _section_key(title: str) -> str:
        """Map a section title to its canonical key name."""
        normalised = re.sub(r"\s+", " ", title.strip().lower())

        # Direct lookup
        if normalised in _SECTION_KEY_MAP:
            return _SECTION_KEY_MAP[normalised]

        # Pattern-based matching
        if normalised.startswith("introduction") or normalised.startswith("background"):
            return "introduction"
        if normalised.startswith("result"):
            return "results"
        if normalised.startswith("discussion"):
            return "discussion"
        if (
            normalised.startswith("method")
            or normalised.startswith("material")
            or normalised.startswith("experimental")
        ):
            return "methods"
        if normalised.startswith("abstract"):
            return "abstract"

        # Unknown section → slug the original title
        return normalised.replace(" ", "_")

    def _fallback_split(self, body: ET.Element) -> dict[str, str]:
        """Fallback: textually split ``<body>`` when no ``<sec>`` elements found.

        Uses the same regex pattern as the markdown source converter.
        """
        text = "".join(body.itertext()).strip()
        if not text:
            return {}
        sections: dict[str, str] = {}
        matches = list(_FALLBACK_SECTION_RE.finditer(text))
        if not matches:
            sections["results"] = text
            return sections

        for i, m in enumerate(matches):
            key = self._section_key(m.group(1))
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
            content = text[start:end].strip()
            if key in sections:
                sections[key] += "\n\n" + content
            else:
                sections[key] = content

        return sections

    # ── Figure parsing ─────────────────────────────────────────────────────────────

    def _parse_figures(self) -> list[str]:
        """Extract figure captions from ``<body>``."""
        root = self._root
        if root is None:
            return []

        body = root.find("body")
        if body is None:
            return []

        figs = body.findall(".//fig")
        result: list[str] = []
        for i, fig in enumerate(figs, start=1):
            label = _elem_text(fig.find("label"))
            caption = _elem_text(fig.find("caption"))
            if not label:
                label = f"Figure {i}"
            if caption:
                combined = f"{label}. {caption}"
            else:
                combined = label
            result.append(combined)

        return result

    # ── Metadata parsing ───────────────────────────────────────────────────────────

    def _parse_metadata(self) -> dict:
        """Extract metadata from ``<front><article-meta>``."""
        root = self._root
        if root is None:
            return {}

        meta_el = root.find("./front/article-meta")
        if meta_el is None:
            # Some JATS variants use front-stub
            meta_el = root.find("./front/front-stub")

        meta: dict = {
            "pmid": None,
            "doi": None,
            "title": "",
            "first_author": "",
            "journal": "",
            "year": None,
        }

        if meta_el is None:
            return meta

        # PMID
        pmid_el = meta_el.find('article-id[@pub-id-type="pmid"]')
        if pmid_el is not None and pmid_el.text:
            meta["pmid"] = pmid_el.text.strip()

        # DOI
        doi_el = meta_el.find('article-id[@pub-id-type="doi"]')
        if doi_el is not None and doi_el.text:
            meta["doi"] = doi_el.text.strip()

        # Title
        title_el = meta_el.find("title-group/article-title")
        if title_el is not None:
            meta["title"] = _elem_text(title_el)

        # First author (surname of first contrib-type='author')
        contribs = meta_el.findall('.//contrib[@contrib-type="author"]')
        if contribs:
            surname_el = contribs[0].find("./name/surname")
            if surname_el is not None and surname_el.text:
                meta["first_author"] = surname_el.text.strip()

        # Journal title
        journal_el = root.find("./front/journal-meta/journal-title-group/journal-title")
        if journal_el is not None and journal_el.text:
            meta["journal"] = journal_el.text.strip()

        # Year
        year_el = meta_el.find("pub-date/year")
        if year_el is not None and year_el.text:
            meta["year"] = year_el.text.strip()
        else:
            # Fallback: scan any pub-date
            for pd in meta_el.findall("pub-date"):
                y = pd.find("year")
                if y is not None and y.text:
                    meta["year"] = y.text.strip()
                    break

        return meta

    # ── Paper name ─────────────────────────────────────────────────────────────────

    def _compute_paper_name(self) -> str:
        """Build paper name in ``{pmid}_{year}_{first_author}_{journal}_{title_slug}`` format.

        PMID prefix ensures global uniqueness and easy PubMed lookup.
        Falls back to ``{year}_{first_author}_{journal}_{title_slug}`` if no PMID.
        """
        meta = self._metadata if self._metadata is not None else self._parse_metadata()
        pmid = str(meta.get("pmid", "")).strip()
        year = str(meta.get("year", "XXXX"))
        author = _slugify(str(meta.get("first_author", "Unknown")), max_len=20)
        journal = _slugify(str(meta.get("journal", "Journal")), max_len=10)
        title_slug = _slugify(str(meta.get("title", "")), max_len=40)

        if pmid:
            name = f"{pmid}_{year}_{author}_{journal}_{title_slug}".strip("_")
        else:
            name = f"{year}_{author}_{journal}_{title_slug}".strip("_")
        # Enforce 120-char limit (was 100, +9 for PMID prefix)
        name = name[:120].rstrip("_")
        # Replace any remaining special characters
        name = _slugify(name, max_len=120)
        return name or "Unknown_Paper"

    # ── Public API ─────────────────────────────────────────────────────────────────

    def get_sections(self) -> dict[str, str]:
        """Return parsed sections by canonical key."""
        if self._sections is None:
            self._sections = self._parse_sections()
        return dict(self._sections)

    def get_figure_blocks(self) -> list[str]:
        """Return parsed figure captions."""
        if self._figures is None:
            self._figures = self._parse_figures()
        return list(self._figures)

    def get_metadata(self) -> dict:
        """Return extracted metadata dict."""
        if self._metadata is None:
            self._metadata = self._parse_metadata()
        return dict(self._metadata)

    def get_paper_name(self) -> str:
        """Return computed paper name."""
        if self._paper_name is None:
            self._paper_name = self._compute_paper_name()
        return self._paper_name

    def get_text(self) -> str:
        """Return full paper text by joining all sections."""
        sections = self.get_sections()
        return "\n\n".join(sections.values())


# ═══════════════════════════════════════════════════════════════════════════════════
#  PubmedSource — PubMed metadata + abstract (PMC fallback)
# ═══════════════════════════════════════════════════════════════════════════════════


def _fetch_pubmed_metadata(pmid: str) -> dict:
    """Fetch paper metadata from NCBI PubMed esummary.

    Returns dict with keys: pmid, doi, title, first_author, journal, year.
    Used as reliable metadata source when PMC XML is unavailable.
    Mirrors the rate-limit + retry pattern of ``PmcXmlSource._ncbi_fetch``.
    """
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=pubmed&id={pmid}&retmode=json"
    )
    time.sleep(0.35)
    for attempt in range(3):
        try:
            req = Request(url)
            req.add_header("User-Agent", "Fuxi/1.0 (paper-converter; academic use)")
            with urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                article = list(data.get("result", {}).values())[1]
        except (HTTPError, URLError, OSError, json.JSONDecodeError, IndexError) as e:
            if isinstance(e, HTTPError) and e.code in (429, 503):
                time.sleep(0.5 * (2**attempt))
                continue
            return {}
        break
    else:
        return {}

    meta: dict = {
        "pmid": pmid,
        "doi": None,
        "title": article.get("title", "").rstrip("."),
        "first_author": None,
        "journal": article.get("source", ""),
        "year": None,
    }

    # Authors: first surname
    authors = article.get("authors", [])
    if authors:
        first = authors[0]
        meta["first_author"] = first.get("name", "").split()[0] if first.get("name") else None

    # Year from pubdate: "2020 Jun 19" → "2020"
    pubdate = article.get("pubdate", "")
    if pubdate:
        yr_match = re.match(r"(\d{4})", str(pubdate))
        if yr_match:
            meta["year"] = yr_match.group(1)

    # DOI from elocationid: "pii: dev185660. doi: 10.1242/dev.185660"
    eloc = article.get("elocationid", "")
    if eloc:
        doi_match = re.search(r"\b(10\.\d{4,}/[^\s]+)", str(eloc))
        if doi_match:
            meta["doi"] = doi_match.group(1).rstrip(".")

    return meta


def _fetch_pubmed_abstract(pmid: str) -> str:
    """Fetch abstract text from NCBI PubMed efetch."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={pmid}&retmode=text&rettype=abstract"
    )
    time.sleep(0.35)
    for attempt in range(3):
        try:
            req = Request(url)
            req.add_header("User-Agent", "Fuxi/1.0 (paper-converter; academic use)")
            with urlopen(req, timeout=15) as resp:
                return resp.read().decode("utf-8").strip()
        except HTTPError as e:
            if e.code in (429, 503):
                time.sleep(0.5 * (2**attempt))
                continue
            return ""
        except (URLError, OSError):
            return ""
    return ""


class PubmedSource(PaperSource):
    """Fetch paper metadata + abstract from PubMed (PMC-unavailable fallback).

    Uses NCBI E-utilities esummary for bibliographic metadata and efetch
    for the abstract text. Does NOT provide full body text or figures.

    Pair with ``Pymupdf4llmSource`` via ``--pdf`` for full-text analysis.

    Parameters
    ----------
    pmid : str
        PubMed ID.
    """

    def __init__(self, pmid: str) -> None:
        if not pmid:
            raise ValueError("PubmedSource requires a PMID")
        self._pmid = str(pmid)
        self._metadata: Optional[dict] = None
        self._abstract: Optional[str] = None
        self._paper_name: Optional[str] = None

    # ── Metadata ────────────────────────────────────────────────────────────────────

    def get_metadata(self) -> dict:
        if self._metadata is None:
            self._metadata = _fetch_pubmed_metadata(self._pmid)
        return dict(self._metadata)

    # ── Abstract ─────────────────────────────────────────────────────────────────────

    def _get_abstract(self) -> str:
        if self._abstract is None:
            self._abstract = _fetch_pubmed_abstract(self._pmid)
        return self._abstract

    def get_sections(self) -> dict[str, str]:
        abstract = self._get_abstract()
        if not abstract:
            return {}
        # Strip leading citation line ("1. Journal. Year... doi:...") from
        # efetch output so the AI only sees the abstract text.
        text_parts = abstract.split("\n\n", 1)
        body = text_parts[1] if len(text_parts) > 1 else abstract
        return {"abstract": body}

    def get_figure_blocks(self) -> list[str]:
        return []

    def get_text(self) -> str:
        sections = self.get_sections()
        return "\n\n".join(sections.values())

    # ── Paper name ─────────────────────────────────────────────────────────────────

    def _compute_paper_name(self) -> str:
        meta = self.get_metadata()
        pmid = str(meta.get("pmid", "")).strip()
        year = str(meta.get("year") or "XXXX")
        author = _slugify(str(meta.get("first_author") or "Unknown"), max_len=20)
        journal = _slugify(str(meta.get("journal") or "Journal"), max_len=10)
        title_slug = _slugify(str(meta.get("title") or ""), max_len=40)
        if pmid:
            name = f"{pmid}_{year}_{author}_{journal}_{title_slug}".strip("_")
        else:
            name = f"{year}_{author}_{journal}_{title_slug}".strip("_")
        name = name[:120].rstrip("_")
        name = _slugify(name, max_len=120)
        return name or "Unknown_Paper"

    def get_paper_name(self) -> str:
        if self._paper_name is None:
            self._paper_name = self._compute_paper_name()
        return self._paper_name


# ═══════════════════════════════════════════════════════════════════════════════════
#  MarkdownSource — Markdown paper file
# ═══════════════════════════════════════════════════════════════════════════════════


class MarkdownSource(PaperSource):
    """Parse a markdown paper file into sections, figures, and metadata.

    Accepts a path to a ``.md`` file produced by PDF-to-markdown conversion
    (e.g. via ``markitdown``).  Splits sections by heading regex and extracts
    figure references using the same patterns as :class:`PaperInsights`.

    Parameters
    ----------
    md_path : str
        Path to the markdown file.
    """

    # Section heading regex (identical to paper_insights._SECTION_RE)
    _SECTION_RE = re.compile(
        r"^(?:#+\s*)?(SUMMARY|Abstract|Introduction|Results|Discussion|Methods|"
        r"Experimental\s*Procedures)\b",
        re.MULTILINE | re.IGNORECASE,
    )
    # Figure reference regex (identical to paper_insights._FIGURE_RE)
    _FIGURE_RE = re.compile(r"(?:Figure|Fig\.?)\s+\d+[a-z]?", re.IGNORECASE)

    def __init__(self, md_path: str) -> None:
        self._md_path = md_path

        raw = Path(md_path).read_text(encoding="utf-8")
        self._raw_text = clean_text(raw)

        # Lazy-loaded caches
        self._sections: Optional[dict[str, str]] = None
        self._figures: Optional[list[str]] = None
        self._metadata: Optional[dict] = None
        self._paper_name: Optional[str] = None

    # ── Section parsing ────────────────────────────────────────────────────────────

    def _parse_sections(self) -> dict[str, str]:
        """Split markdown into sections by header.

        Returns a dict with lowercase keys (abstract, introduction, results,
        discussion, methods).  Missing sections are empty strings.
        No matching headers places all text under ``results``.
        """
        text = self._raw_text

        sections: dict[str, str] = {
            k: "" for k in ("abstract", "introduction", "results", "discussion", "methods")
        }

        matches = list(self._SECTION_RE.finditer(text))
        if not matches:
            sections["results"] = text.strip()
            return sections

        # Section key map (mirrors paper_insights.PaperInsights.split_sections)
        section_key_map: dict[str, str] = {
            "summary": "abstract",
            "abstract": "abstract",
            "introduction": "introduction",
            "results": "results",
            "discussion": "discussion",
            "methods": "methods",
            "experimental procedures": "methods",
            "experimentalprocedures": "methods",
        }

        for i, m in enumerate(matches):
            raw_key = m.group(1).lower().strip()
            key = section_key_map.get(raw_key, "unknown")
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
            sections[key] = text[start:end].strip()

        return sections

    # ── Figure parsing ─────────────────────────────────────────────────────────────

    def _parse_figures(self) -> list[str]:
        """Extract figure blocks from the full paper text.

        Splits on ``Figure`` / ``Fig`` references (same pattern as
        :meth:`PaperInsights.extract_figure_blocks`) but operates on the
        entire paper text rather than only the Results section.
        """
        text = self._raw_text
        if not text.strip():
            return []

        blocks = self._FIGURE_RE.split(text)
        refs = self._FIGURE_RE.findall(text)
        result: list[str] = []
        for i, block in enumerate(blocks[1:], start=1):
            header = refs[i - 1] if i - 1 < len(refs) else f"Figure_{i}"
            cleaned = block.strip().lstrip(".: ")
            combined = f"{header}. {cleaned}"
            if combined.strip():
                result.append(combined)
        if not result and text.strip():
            result.append(text.strip())
        return result

    # ── Metadata parsing ───────────────────────────────────────────────────────────

    def _parse_filename_meta(self) -> dict:
        """Parse year/author/journal from filename.

        Expects format: ``{year}_{author}_{journal}_{title}.md``
        (mirrors :func:`paper_insights._parse_filename_meta`).
        """
        stem = Path(self._md_path).stem
        parts = stem.split("_")
        meta: dict = {
            "filename": Path(self._md_path).name,
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

    # ── Paper name ─────────────────────────────────────────────────────────────────

    def _compute_paper_name(self) -> str:
        """Build paper name in ``{pmid}_{year}_{first_author}_{journal}_{title}`` format.

        For PDF sources, PMID is extracted from filename if named as '<PMID>.pdf'.
        """
        meta = self._metadata if self._metadata is not None else self._parse_filename_meta()
        pmid = str(meta.get("pmid", "")).strip()
        year = str(meta.get("year") or "XXXX")
        author = _slugify(str(meta.get("first_author") or "Unknown"), max_len=20)
        journal = _slugify(str(meta.get("journal") or "Journal"), max_len=10)
        title_slug = _slugify(str(meta.get("title") or ""), max_len=40)

        if pmid:
            name = f"{pmid}_{year}_{author}_{journal}_{title_slug}".strip("_")
        else:
            name = f"{year}_{author}_{journal}_{title_slug}".strip("_")
        name = name[:120].rstrip("_")
        name = _slugify(name, max_len=120)
        return name or "Unknown_Paper"

    # ── Public API ─────────────────────────────────────────────────────────────────

    def get_sections(self) -> dict[str, str]:
        """Return parsed sections by canonical key."""
        if self._sections is None:
            self._sections = self._parse_sections()
        return dict(self._sections)

    def get_figure_blocks(self) -> list[str]:
        """Return extracted figure blocks."""
        if self._figures is None:
            self._figures = self._parse_figures()
        return list(self._figures)

    def get_metadata(self) -> dict:
        """Return filename-derived metadata (year, first_author, journal, title)."""
        if self._metadata is None:
            self._metadata = self._parse_filename_meta()
        return dict(self._metadata)

    def get_paper_name(self) -> str:
        """Return computed paper name."""
        if self._paper_name is None:
            self._paper_name = self._compute_paper_name()
        return self._paper_name

    def get_text(self) -> str:
        """Return the full raw paper text."""
        return self._raw_text


# ═══════════════════════════════════════════════════════════════════════════════════
#  Pymupdf4llmSource — PDF via pymupdf4llm (optional)
# ═══════════════════════════════════════════════════════════════════════════════════


class Pymupdf4llmSource(PaperSource):
    """Convert a PDF paper via pymupdf4llm (optional), delegate to MarkdownSource."""

    def __init__(self, pdf_path: str, pmid: str | None = None) -> None:
        # Lazy import — raises ImportError if pymupdf4llm not installed
        try:
            import fitz  # pymupdf
            import pymupdf4llm
        except ImportError:
            raise ImportError(
                "pymupdf4llm not installed. Run: pip install -r requirements/paper.txt"
            ) from None

        # Convert PDF to markdown text
        doc = fitz.open(pdf_path)
        md_text = pymupdf4llm.to_markdown(doc)
        doc.close()

        # Apply clean_text post-processor
        md_text = clean_text(md_text)

        # Write to temp .md file and delegate to MarkdownSource
        import tempfile

        tmp_fd, tmp_path = tempfile.mkstemp(suffix=".md", prefix="pymupdf4llm_", text=True)
        with os.fdopen(tmp_fd, "w", encoding="utf-8") as f:
            f.write(md_text)

        self._md_source = MarkdownSource(tmp_path)
        self._pdf_path = pdf_path
        self._tmp_path = tmp_path  # store for cleanup in __del__
        self._pmid = pmid  # optional PubMed ID for authoritative metadata

    def __del__(self) -> None:
        """Clean up the temporary .md file on garbage collection."""
        if hasattr(self, "_tmp_path") and os.path.exists(self._tmp_path):
            try:
                os.unlink(self._tmp_path)
            except OSError:
                pass

    # ── Metadata from PDF filename (not temp .md file) ─────────────────

    def _parse_pdf_filename_meta(self) -> dict:
        """Parse year/author/journal from the original PDF filename.

        Supports two formats:
          {year}_{author}_{journal}_{title}.pdf
          {pmid}.pdf  (just a PubMed ID)
        """
        fname = Path(self._pdf_path).stem  # e.g. '2026_Wohlschlegel_CellRep_RA-Foveal-Development'
        parts = fname.split("_")
        meta: dict = {
            "filename": Path(self._pdf_path).name,
            "stem": fname,
            "year": None,
            "first_author": None,
            "journal": None,
            "title": fname,
        }
        # Pattern: {year}_{author}_{journal}_{title}
        if len(parts) >= 3 and parts[0].isdigit() and len(parts[0]) in (2, 4):
            meta["year"] = parts[0]
            meta["first_author"] = parts[1]
            meta["journal"] = parts[2]
            meta["title"] = "_".join(parts[3:]) if len(parts) > 3 else parts[2]
        else:
            # Fallback: maybe the PDF file name is just a PMID like "32467236.pdf"
            if parts[0].isdigit() and len(parts[0]) in (7, 8):
                meta["pmid"] = parts[0]
        return meta

    def _parse_text_title(self) -> str | None:
        """Extract the paper title from the markdown text.

        Strategy: find hash-marked headings, prefer the longest one that is not
        a known section heading. The paper title is typically the most prominent
        (longest) heading in the preamble before the main sections.
        """
        _section_keywords = {
            "abstract",
            "introduction",
            "results",
            "discussion",
            "methods",
            "materials",
            "acknowledgements",
            "acknowledgments",
            "references",
            "supplementary",
            "figures",
            "figure legends",
            "abbreviations",
            "conflict",
            "author contributions",
            "data availability",
            "keywords",
            "key words",
            "acknowledgment",
            "conclusions",
            "declarations",
            "consent",
            "ethics",
            "stem cells",
            "techniques",
            "resources",
        }
        candidates = []
        for line in self._md_source.get_text().split("\n"):
            stripped = line.strip()
            if not stripped or not stripped.startswith("#"):
                continue
            text = stripped.lstrip("#").strip().rstrip(".")
            if len(text) < 20:
                continue
            lower = text.lower()
            # Skip known section headings
            if lower in _section_keywords or any(lower.startswith(w) for w in _section_keywords):
                continue
            # Skip copyright / DOI / license lines that start with a hash
            if any(text.lower().startswith(w) for w in ["copyright", "©", "doi:", "the author"]):
                continue
            candidates.append(text)
        # Return the longest hash line (most likely the title)
        if candidates:
            candidates.sort(key=len, reverse=True)
            return candidates[0]
        # Fallback: first non-empty line that looks like a title
        for line in self._md_source.get_text().split("\n"):
            stripped = line.strip()
            if len(stripped) > 30 and not any(
                stripped.lower().startswith(w)
                for w in [
                    "figure",
                    "fig.",
                    "abstract",
                    "introduction",
                    "copyright",
                    "©",
                    "published",
                    "journal",
                    "the author",
                    "doi:",
                ]
            ):
                return stripped
        return None

    def get_metadata(self) -> dict:
        """Return paper metadata, preferring NCBI PubMed when PMID is available.

        When ``--pmid`` is provided, use authoritative NCBI esummary data.
        Falls back to PDF-filename-derived metadata otherwise, with text-title
        extraction when the filename doesn't follow the structured pattern.
        """
        if self._pmid:
            return _fetch_pubmed_metadata(self._pmid)  # authoritative

        meta = self._parse_pdf_filename_meta()
        if (
            meta.get("title")
            and not meta["title"].startswith("pymupdf4llm_")
            and not str(meta.get("year") or "").isdigit()
        ):
            # PDF filename gave us nothing useful — try text-derived title
            text_title = self._parse_text_title()
            if text_title:
                meta["title"] = text_title[:200]
        return meta

    def get_paper_name(self) -> str:
        """Build paper name from PDF-derived metadata instead of temp-file name."""
        meta = self.get_metadata()
        year = str(meta.get("year") or "XXXX")
        author = _slugify(str(meta.get("first_author") or "PDF"), max_len=20)
        journal = _slugify(str(meta.get("journal") or "Paper"), max_len=10)
        title_slug = _slugify(str(meta.get("title") or ""), max_len=40)

        name = f"{year}_{author}_{journal}_{title_slug}".strip("_")
        name = name[:100].rstrip("_")
        name = _slugify(name, max_len=100)
        return name or "Unknown_Paper"

    def get_sections(self) -> dict[str, str]:
        return self._md_source.get_sections()

    def get_figure_blocks(self) -> list[str]:
        return self._md_source.get_figure_blocks()

    def get_text(self) -> str:
        return self._md_source.get_text()
