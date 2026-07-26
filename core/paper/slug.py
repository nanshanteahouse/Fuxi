"""Slug generation for registered papers.

Format
------
    {author}{year}_{topic}

- author: NFKD-folded, lowercase, alnum-only, <=15 chars.
  Examples: "Macosko" -> macosko, "Kläsch" -> klasch,
            "Vidal-Vázquez" -> vidalvazquez, "Wohlschlegel" -> wohlschlegel.
- year: 4-digit year.
- topic: 2-3 substantive lowercase words derived from insights.yaml
  `experimental_design` (preferred) or `paper_meta.title`, with species names
  and generic filler filtered out.

PMID (`paper_id`) remains the canonical lookup key. Slug is a human-readable
alias for display and filesystem navigation.

Usage
-----
    from core.paper.slug import build_slug

    slug = build_slug(
        first_author=meta["first_author"],
        year=meta["year"],
        paper_meta=meta,
        insights=insights_dict,
    )

Design notes
------------
- The topic vocabulary is Fuxi-aware (retina-focused). Domain-redundant words
  like "retina", "cell", "single-cell" are stopped because every paper in the
  registry has them — they carry no discriminative signal.
- Insights.tissue_info is preferred over summary because the latter is
  verb-heavy and yields noisy topics.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Any

# ---------------------------------------------------------------------------
# Vocabulary filters
# ---------------------------------------------------------------------------

_SPECIES = {
    "mouse",
    "mice",
    "murine",
    "human",
    "humans",
    "zebrafish",
    "rat",
    "rats",
    "macaque",
    "macaques",
    "marmoset",
    "marmosets",
    "cynomolgus",
    "rhesus",
    "xenopus",
    "drosophila",
    "cichlid",
    "cichlids",
    "killifish",
    "shark",
    "sharks",
    "salmon",
    "bovine",
    "porcine",
    "chick",
    "chicken",
    "chicks",
    "fly",
    "flies",
    "frog",
    "frogs",
    "fish",
    "primate",
    "primates",
    "vertebrate",
    "vertebrates",
    "catshark",
}

_STOPWORDS = _SPECIES | {
    # articles / prepositions / conjunctions
    "the",
    "of",
    "and",
    "or",
    "in",
    "on",
    "at",
    "from",
    "by",
    "for",
    "to",
    "with",
    "without",
    "via",
    "based",
    "during",
    "into",
    "across",
    "between",
    "through",
    "their",
    "this",
    "that",
    "these",
    "those",
    "a",
    "an",
    "after",
    "before",
    "versus",
    "vs",
    "including",
    "include",
    "includes",
    # generic scRNA vocabulary
    "cell",
    "cells",
    "cellular",
    "tissue",
    "tissues",
    "type",
    "types",
    "subtype",
    "subtypes",
    "single",
    "single-cell",
    "singlecell",
    "transcriptomic",
    "transcriptomics",
    "rnaseq",
    "rna",
    "seq",
    "scrnaseq",
    "analysis",
    "atlas",
    "atlases",
    "profiling",
    "expression",
    "study",
    "characterization",
    "characterising",
    "characterization",
    "comprehensive",
    "molecular",
    "systematic",
    "spatiotemporal",
    "spatialtemporal",
    "highresolution",
    "high",
    "resolution",
    "genome-wide",
    "genomewide",
    "genome",
    "wide",
    "parallel",
    "highly",
    "novel",
    "new",
    "key",
    "regulatory",
    "regulation",
    "regulated",
    "identification",
    "identify",
    "identifies",
    "identifying",
    "identified",
    "reveals",
    "revealed",
    "revealing",
    "elucidates",
    "elucidate",
    "deciphering",
    "decipher",
    "dissecting",
    "dissect",
    "deconstructing",
    "deconstruct",
    "mapping",
    "mapped",
    "investigate",
    "investigates",
    "investigating",
    "investigated",
    "explore",
    "explores",
    "explored",
    "exploring",
    "compare",
    "comparing",
    "compared",
    "comparison",
    "differential",
    "differentially",
    "dynamics",
    "dynamic",
    "landscape",
    "landscapes",
    "diversity",
    "heterogeneity",
    "signature",
    "signatures",
    "profile",
    "profiles",
    "population",
    "populations",
    "association",
    "associated",
    "logic",
    "program",
    "programs",
    "programming",
    "reprogramming",
    # verb / process noise (especially from summaries)
    "applied",
    "apply",
    "uses",
    "used",
    "using",
    "use",
    "generated",
    "generates",
    "generating",
    "generate",
    "integrated",
    "integrates",
    "integrating",
    "integrate",
    "performed",
    "performs",
    "performing",
    "perform",
    "run",
    "running",
    "runs",
    "contains",
    "containing",
    "contain",
    "obtained",
    "obtains",
    "obtaining",
    "obtain",
    "derived",
    "derives",
    "deriving",
    "derive",
    "cultured",
    "cultures",
    "culturing",
    "culture",
    "induced",
    "induces",
    "inducing",
    "induce",
    "depleted",
    "depletes",
    "depleting",
    "deplete",
    "encompassing",
    "encompasses",
    "encompass",
    "covered",
    "covers",
    "cover",
    "lined",
    "lines",
    "line",
    "enriched",
    "enrich",
    "selected",
    "select",
    "sequencing",
    "sequence",
    "sequences",
    "full",
    "complete",
    "entire",
    "transcriptionally",
    "transcriptional",
    # developmental / generic descriptors
    "development",
    "developmental",
    "developing",
    "developed",
    "develop",
    "differentiation",
    "differentiated",
    "differentiating",
    "differentiate",
    "pluripotent",
    "stem",
    "invitro",
    "vivo",
    "vitro",
    "in vitro",
    "in vivo",
    "wildtype",
    "wild-type",
    "wild",
    "wt",
    "normal",
    "control",
    "controls",
    "adult",
    "postnatal",
    "prenatal",
    "fetal",
    "embryonic",
    "embryo",
    "embryos",
    "early",
    "late",
    "stage",
    "stages",
    "time",
    "timepoint",
    "timepoints",
    "aged",
    "ages",
    "age",
    "senescence",
    "senescent",
    "post",
    "pre",
    "peri",
    "neonates",
    "neonatal",
    "newborn",
    "birth",
    # number-as-words
    "one",
    "two",
    "three",
    "four",
    "five",
    "six",
    "seven",
    "eight",
    "nine",
    "ten",
    "eleven",
    "twelve",
    "thirteen",
    "fourteen",
    "fifteen",
    "sixteen",
    "first",
    "second",
    "third",
    "fourth",
    "fifth",
    # generic / metadata descriptors
    "approximately",
    "approx",
    "around",
    "whole",
    "multiple",
    "double",
    "triple",
    "independent",
    "independently",
    "published",
    "unpublished",
    "diverse",
    "various",
    "several",
    "some",
    "any",
    "all",
    "each",
    "every",
    "samples",
    "sample",
    "specimens",
    "specimen",
    "donors",
    "donor",
    "patients",
    "patient",
    "cases",
    "case",
    "weeks",
    "week",
    "days",
    "day",
    "months",
    "month",
    "years",
    "year",
    "hours",
    "hour",
    "points",
    "point",
    "non",
    "under",
    "within",
    "also",
    "than",
    "then",
    "female",
    "male",
    "sex",
    "sexes",
    # platform / library noise
    "platforms",
    "platform",
    "methods",
    "method",
    "data",
    "datasets",
    "dataset",
    # filler verbs / aux
    "is",
    "are",
    "was",
    "were",
    "be",
    "been",
    "being",
    "have",
    "has",
    "had",
    "do",
    "does",
    "did",
    "will",
    "would",
    "can",
    "could",
    "may",
    "might",
    "we",
    "our",
    "their",
    "his",
    "her",
    "its",
    "i",
    "not",
    "no",
    # domain-redundant (Fuxi is retina-focused)
    "retina",
    "retinal",
    "retinas",
    "eye",
    "eyes",
    "ocular",
    "organoid",
    "organoids",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ascii_fold(text: str) -> str:
    """NFKD fold to ASCII (Müller -> Muller, Kläsch -> Klasch)."""
    if not text:
        return ""
    nfkd = unicodedata.normalize("NFKD", text)
    return "".join(c for c in nfkd if not unicodedata.combining(c))


def normalize_author(name: str) -> str:
    """NFKD fold, lowercase, alnum-only, <=15 chars. 'unknown' on empty."""
    if not name:
        return "unknown"
    folded = _ascii_fold(name).lower()
    alnum = re.sub(r"[^a-z]", "", folded)
    return alnum[:15] or "unknown"


def _extract_topic_tokens(text: str, max_tokens: int = 3) -> list[str]:
    """ASCII fold, lowercase, strip parentheticals, drop stopwords, take first N."""
    if not text:
        return []
    text = _ascii_fold(text).lower()
    text = re.sub(r"\([^)]*\)", " ", text)  # drop (BCs), (WT), (hPSC) etc.
    text = re.sub(r"[^a-z\s\-]", " ", text)  # keep alpha + hyphen
    text = text.replace("-", " ")
    raw = re.findall(r"[a-z]+", text)
    out: list[str] = []
    seen: set[str] = set()
    for tok in raw:
        if tok in _STOPWORDS or len(tok) < 3:
            continue
        if tok in seen:
            continue
        seen.add(tok)
        out.append(tok)
        if len(out) >= max_tokens:
            break
    return out


def derive_topic(
    paper_meta: dict[str, Any], insights: dict[str, Any] | None = None
) -> tuple[str, str]:
    """Return (topic, source_field).

    Priority chain -- prefer the first source yielding >=2 substantive tokens:

      1. insights.experimental_design.tissue_info  (preferred: concise phrase)
      2. paper_meta.title                          (curated, usually clean)
      3. insights.experimental_design.summary      (verb-heavy, last resort)
      4. insights.experimental_design.tissue       (single word, e.g. 'retina')
      5. 'study' fallback

    Returns
    -------
    (topic, source) where source in
        {'tissue_info', 'title', 'summary', 'tissue', 'fallback'}.
    """
    ins = insights or {}

    tissue_info = ins.get("experimental_design", {}).get("tissue_info", "")
    if tissue_info:
        toks = _extract_topic_tokens(tissue_info)
        if len(toks) >= 2:
            return "_".join(toks), "tissue_info"

    title = paper_meta.get("title", "")
    if title:
        toks = _extract_topic_tokens(title)
        if len(toks) >= 2:
            return "_".join(toks), "title"

    # Last-resort chain: accept even 1-token topics
    if title:
        toks = _extract_topic_tokens(title)
        if toks:
            return "_".join(toks), "title"

    if tissue_info:
        toks = _extract_topic_tokens(tissue_info)
        if toks:
            return "_".join(toks), "tissue_info"

    summary = ins.get("experimental_design", {}).get("summary", "")
    if summary:
        toks = _extract_topic_tokens(summary)
        if toks:
            return "_".join(toks), "summary"

    tissue = ins.get("experimental_design", {}).get("tissue", "")
    if tissue:
        return tissue, "tissue"

    return "study", "fallback"


def build_slug(
    first_author: str,
    year: str,
    paper_meta: dict[str, Any],
    insights: dict[str, Any] | None = None,
) -> str:
    """Build a slug: '{author}{year}_{topic}'.

    Parameters
    ----------
    first_author : str
        Surname of first author (e.g. 'Shekhar', 'Kläsch', 'Vidal-Vázquez').
    year : str
        4-digit publication year.
    paper_meta : dict
        Must contain at least 'title'. Used for title-based topic derivation.
    insights : dict, optional
        Full insights.yaml content. Preferred source for topic derivation
        via experimental_design.tissue_info.

    Returns
    -------
    slug : str
        Example: 'shekhar2016_classification_bipolar_neurons'
    """
    author = normalize_author(first_author)
    yr = (year or "0000")[:4]
    topic, _ = derive_topic(paper_meta, insights)
    # Sanitize topic for filesystem safety
    clean_topic = re.sub(r"[^a-z0-9_]", "", topic.lower()) or "study"
    return f"{author}{yr}_{clean_topic}"
