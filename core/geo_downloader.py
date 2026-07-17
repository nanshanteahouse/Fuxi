#!/usr/bin/env python3
"""
geo_downloader.py — GEO dataset downloader for Fuxi
====================================================

Downloads GEO (Gene Expression Omnibus) datasets from NCBI servers
using wget/curl as the transfer backend. Fetches SOFT-format metadata
and supplementary data files to ``$FUXI_DATA_ROOT/GSE*``.

Integration points:
  1. Standalone CLI:  ``python core/geo_downloader.py --gse GSE123456``
  2. Preprocess:      ``preprocessor.py --gse GSE... --download``
  3. Registry:        ``python -m core.registry add-paper --pmid X --download``

Usage::

    python core/geo_downloader.py --gse GSE118614
    python core/geo_downloader.py --gse GSE118614 --dry-run
    python core/geo_downloader.py --gse GSE118614 --skip-soft
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional

# Add repo root to sys.path (consistent with all step scripts)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════

_NCBI_USER_AGENT = "Fuxi/1.0 (academic-use; geo-research)"
_NCBI_RATE_LIMIT = 0.35  # seconds between requests (3 req/s default limit)

# URL templates
_SOFT_GZ_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_nnn}/{gse}/soft/{gse}_family.soft.gz"
_SUPPL_LIST_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_nnn}/{gse}/suppl/"

# SOFT field mapping: !header → returned dict key
_SOFT_FIELD_MAP: dict[str, str] = {
    "Series_title":              "title",
    "Series_geo_accession":      "gse_id",
    "Series_status":             "status",
    "Series_submission_date":    "submission_date",
    "Series_last_update_date":   "last_update_date",
    "Series_summary":            "summary",
    "Series_overall_design":     "overall_design",
    "Series_type":               "series_type",
    "Platform_title":            "platform_title",
}

# These accumulate to a list
_SOFT_LIST_FIELDS = {"Series_pubmed_id", "Series_contributor"}

# ═══════════════════════════════════════════════════════════════════
# Downloader detection (wget / curl)
# ═══════════════════════════════════════════════════════════════════

def _detect_downloader() -> str:
    """Return the path to the best available downloader.

    Prefers ``wget`` (better resume UX) with ``curl`` as fallback.

    Returns:
        Full path string to the detected downloader binary.

    Raises:
        RuntimeError: Neither ``wget`` nor ``curl`` is available.
    """
    for cmd in ("wget", "curl"):
        path = shutil.which(cmd)
        if path:
            log.debug("Downloader detected: %s at %s", cmd, path)
            return cmd  # return the short name — used for command dispatch
    raise RuntimeError(
        "Neither 'wget' nor 'curl' found. "
        "Install one:  sudo apt install wget  (recommended) "
        "or:  sudo apt install curl"
    )


def _get_downloader_name() -> str:
    """Cached lazy detection of downloader name."""
    if not hasattr(_get_downloader_name, "_cached"):
        _get_downloader_name._cached = _detect_downloader()  # type: ignore[attr-defined]
    return _get_downloader_name._cached  # type: ignore[attr-defined]


# ═══════════════════════════════════════════════════════════════════
# NCBI HTTP helpers
# ═══════════════════════════════════════════════════════════════════

def _ncbi_api_key() -> Optional[str]:
    """Read NCBI_API_KEY from environment (empty string → None)."""
    key = os.environ.get("NCBI_API_KEY", "").strip()
    return key or None


def _ncbi_fetch(url: str, raw: bool = False) -> bytes | str:
    """Rate-limited HTTPS GET with NCBI-appropriate User-Agent.

    Args:
        url: Full HTTPS URL to fetch.
        raw: If True, return bytes; otherwise decode as UTF-8.

    Returns:
        Response body as bytes or str.

    Raises:
        urllib.error.HTTPError: On non-2xx responses.
        urllib.error.URLError: On network failures.
    """
    time.sleep(_NCBI_RATE_LIMIT)
    req = urllib.request.Request(url, headers={"User-Agent": _NCBI_USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = resp.read()
    return data if raw else data.decode("utf-8", errors="replace")


# ═══════════════════════════════════════════════════════════════════
# URL construction
# ═══════════════════════════════════════════════════════════════════

def _gse_nnn(gse_id: str) -> str:
    """Convert GSE accession to NCBI directory pattern.

    NCBI groups series by prefix: everything except the last 3 digits → ``nnn``.
    ``GSE107618`` → ``GSE107nnn``  (6-digit, prefix=107)
    ``GSE81905``  → ``GSE81nnn``   (5-digit, prefix=81)
    """
    m = re.match(r'GSE(\d+)', gse_id.upper())
    if not m:
        raise ValueError(f"Invalid GSE accession: {gse_id!r}")
    digits = m.group(1)
    prefix = digits[:-3] if len(digits) > 3 else digits
    return f"GSE{prefix}nnn"


def _build_soft_url(gse_id: str) -> str:
    """Build the SOFT.gz download URL for a GSE accession."""
    nnn = _gse_nnn(gse_id)
    return _SOFT_GZ_URL.format(gse_nnn=nnn, gse=gse_id)


def _build_suppl_list_url(gse_id: str) -> str:
    """Build the supplementary file listing URL."""
    nnn = _gse_nnn(gse_id)
    return _SUPPL_LIST_URL.format(gse_nnn=nnn, gse=gse_id)


def _build_file_url(gse_id: str, filename: str) -> str:
    """Build the download URL for a specific supplementary file."""
    nnn = _gse_nnn(gse_id)
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{gse_id}/suppl/{filename}"


# ═══════════════════════════════════════════════════════════════════
# HUMAN-READABLE SIZE
# ═══════════════════════════════════════════════════════════════════

def _human_size(n_bytes: int) -> str:
    """Convert byte count to human-readable string.

    >>> _human_size(264799)
    '258.6 KB'
    """
    if n_bytes < 1024:
        return f"{n_bytes} B"
    for unit in ("KB", "MB", "GB", "TB"):
        n_bytes /= 1024.0
        if n_bytes < 1024:
            return f"{n_bytes:.1f} {unit}"
    return f"{n_bytes:.1f} PB"


# ═══════════════════════════════════════════════════════════════════
# SOFT metadata parser
# ═══════════════════════════════════════════════════════════════════

def _parse_soft_metadata(text: str) -> dict:
    """Parse SOFT format text into a structured metadata dict.

    Extracts all ``!Series_*``, ``!Platform_*`` header fields,
    sample blocks (``^SAMPLE``), and organism info.

    Args:
        text: Raw SOFT format text (after gzip decompression).

    Returns:
        Dict with keys: gse_id, title, pmid (list), summary,
        overall_design, platform_title, contributors (list),
        organism, n_samples, sample_list, is_superseries, ...
    """
    meta: dict = {
        "gse_id": "",
        "title": "",
        "pmid": [],
        "summary": "",
        "overall_design": "",
        "status": "",
        "submission_date": "",
        "last_update_date": "",
        "series_type": "",
        "platform_title": "",
        "contributors": [],
        "organism": "",
        "n_samples": 0,
        "sample_list": [],
        "is_superseries": False,
    }

    lines = text.split("\n")

    # ── Detect SuperSeries ──
    for line in lines[:50]:
        if "SuperSeries" in line and ("composed of" in line or "This SuperSeries" in line):
            meta["is_superseries"] = True
            break

    # ── Parse header fields ──
    for line in lines:
        if not line.startswith("!"):
            continue
        line = line[1:]  # strip leading '!'
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip()

        if key in _SOFT_FIELD_MAP:
            meta[_SOFT_FIELD_MAP[key]] = value
        elif key in _SOFT_LIST_FIELDS:
            if _SOFT_FIELD_MAP.get(key):
                meta[_SOFT_FIELD_MAP[key]] = value  # type: ignore[index]
            elif key == "Series_pubmed_id":
                meta["pmid"].append(value)
            elif key == "Series_contributor":
                meta["contributors"].append(value)

    # ── Parse sample blocks ──
    sample_blocks = text.split("^SAMPLE")[1:]  # skip pre-first-sample text
    organisms: set[str] = set()
    for block in sample_blocks:
        sample: dict = {"accession": "", "title": "", "organism": ""}
        for line in block.split("\n"):
            line = line.strip()
            if line.startswith("!Sample_geo_accession"):
                sample["accession"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_title"):
                sample["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_organism_ch1"):
                org = line.split("=", 1)[1].strip()
                sample["organism"] = org
                organisms.add(org)
        if sample["accession"]:
            meta["sample_list"].append(sample)

    meta["n_samples"] = len(meta["sample_list"])
    meta["organism"] = _resolve_organism(organisms)

    return meta


def _resolve_organism(organisms: set[str]) -> str:
    """Pick the most specific common organism from a set.

    If there are multiple species, returns the first one (sorted).
    """
    if not organisms:
        return ""
    organisms.discard("")
    if len(organisms) == 1:
        return organisms.pop()
    # Multiple: return comma-joined or pick first
    return ", ".join(sorted(organisms))


# ═══════════════════════════════════════════════════════════════════
# SOFT metadata fetching (public)
# ═══════════════════════════════════════════════════════════════════

def fetch_soft_metadata(gse_id: str) -> dict:
    """Download and parse SOFT.gz metadata for a GEO series.

    Args:
        gse_id: GEO accession ID (e.g. ``GSE118614``).

    Returns:
        Structured metadata dict. See :func:`_parse_soft_metadata`
        for the full schema.

    Raises:
        urllib.error.HTTPError: If the SOFT file is not available.
        ValueError: If the GSE ID format is invalid.
    """
    url = _build_soft_url(gse_id)
    log.info("Fetching SOFT metadata: %s", url)

    raw_bytes = _ncbi_fetch(url, raw=True)
    if not isinstance(raw_bytes, bytes):
        raise RuntimeError(f"Expected bytes from SOFT fetch, got {type(raw_bytes)}")

    try:
        text = gzip.decompress(raw_bytes).decode("utf-8", errors="replace")
    except gzip.BadGzipFile:
        raise RuntimeError(f"Failed to decompress SOFT.gz for {gse_id} — corrupted file")

    meta = _parse_soft_metadata(text)
    meta["gse_id"] = gse_id.upper()
    return meta


# ═══════════════════════════════════════════════════════════════════
# Supplementary file listing
# ═══════════════════════════════════════════════════════════════════

def _parse_suppl_html(html: str) -> list[dict]:
    """Parse NCBI FTP directory listing into file metadata.

    NCBI serves an HTML table where each data file row looks like::

        <a href="FILENAME">FILENAME</a>   YYYY-MM-DD HH:MM  SIZE

    We parse line-by-line to avoid cross-line regex matches.

    Returns:
        List of dicts with keys: name, size_bytes, size_human, is_raw_tar.
        Directory entries and non-data files are excluded.
    """
    files: list[dict] = []
    # Match each file row: <a href="NAME">...</a>  DATE  TIME  SIZE
    file_row = re.compile(
        r'<a\s+href="([^"]+?)">[^<]+</a>'
        r'\s+(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2})'
        r'\s+([\d.]+[KMGT]?)',
    )

    for line in html.splitlines():
        m = file_row.search(line)
        if not m:
            continue
        name = m.group(1).strip()
        size_str = m.group(3).strip()

        # Skip directory navigation links and non-data files
        if name in ("/", "..", "Parent Directory"):
            continue
        if name.startswith("/") or name.startswith("/geo/"):
            continue
        if name.endswith(".html"):
            continue

        size_bytes = _parse_ftp_size(size_str)
        files.append({
            "name": name,
            "size_bytes": size_bytes,
            "size_human": _human_size(size_bytes) if size_bytes else "?",
            "is_raw_tar": bool(re.search(r'RAW.*\.tar', name, re.IGNORECASE)),
        })
    return files


def _parse_ftp_size(size_str: str) -> int:
    """Parse NCBI FTP size string ('373M', '1.1M', '631K') to bytes."""
    size_str = size_str.strip().upper()
    multipliers: dict[str, int] = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    for suffix, mult in multipliers.items():
        if size_str.endswith(suffix):
            try:
                return int(float(size_str[:-1]) * mult)
            except ValueError:
                return 0
    try:
        return int(float(size_str))
    except ValueError:
        return 0


def list_suppl_files(gse_id: str) -> list[dict]:
    """Fetch the supplementary file listing for a GEO series.

    Args:
        gse_id: GEO accession ID.

    Returns:
        List of file descriptors with download URLs attached.

    Raises:
        urllib.error.HTTPError: If the listing page is not available.
    """
    url = _build_suppl_list_url(gse_id)
    log.info("Fetching suppl listing: %s", url)

    html = _ncbi_fetch(url)
    if not isinstance(html, str):
        raise RuntimeError(f"Expected str from suppl listing, got {type(html)}")

    files = _parse_suppl_html(html)

    # Attach full download URLs
    for f in files:
        f["url"] = _build_file_url(gse_id, f["name"])

    # Sort: RAW.tar files first, then by name
    files.sort(key=lambda x: (not x["is_raw_tar"], x["name"].lower()))
    return files


# ═══════════════════════════════════════════════════════════════════
# File download via wget / curl
# ═══════════════════════════════════════════════════════════════════

def download_file(
    url: str,
    dest_path: str,
    resume: bool = True,
    show_progress: bool = True,
) -> bool:
    """Download a single file using wget or curl.

    Args:
        url: Full download URL.
        dest_path: Local path to write the file to.
        resume: Use ``--continue`` / ``-C -`` to resume partial downloads.
        show_progress: Show transfer progress bar.

    Returns:
        ``True`` if download succeeded (exit code 0).
    """
    downloader = _get_downloader_name()
    os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)

    if downloader == "wget":
        cmd = ["wget"]
        if resume:
            cmd.append("--continue")
        if show_progress:
            cmd.append("--show-progress")
        else:
            cmd.append("--quiet")
        cmd += ["-O", dest_path, url]
    else:  # curl
        cmd = ["curl", "-L", "--fail"]
        if resume:
            cmd += ["-C", "-"]
        if show_progress:
            cmd.append("-#")
        else:
            cmd.append("--silent")
        cmd += ["-o", dest_path, url]

    log.debug("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=False)
        return result.returncode == 0
    except Exception as exc:
        return False
    except Exception as exc:
        log.error("Download error: %s — %s", url, exc)
        return False


# ═══════════════════════════════════════════════════════════════════
# Main orchestrator
# ═══════════════════════════════════════════════════════════════════

def download_gse(
    gse_id: str,
    dest_dir: str,
    dry_run: bool = False,
    skip_existing: bool = True,
    fetch_meta: bool = True,
    quiet: bool = False,
) -> dict:
    """Download all supplementary data for a GEO series.

    Orchestrates the full download pipeline:
    1. Fetch SOFT.gz metadata
    2. List supplementary files from NCBI FTP
    3. Download each file via wget/curl (with resume)

    Args:
        gse_id: GEO accession ID.
        dest_dir: Directory to write files into (usually
            ``$FUXI_DATA_ROOT/GSEXXXXXX``).
        dry_run: Report what *would* happen without downloading.
        skip_existing: Skip files that already exist with matching size.
        fetch_meta: Fetch SOFT metadata (set False to skip).
        quiet: Suppress per-file progress output.

    Returns:
        Dict with keys: gse_id, dry_run, metadata, files (list of
        {name, status, size_human}), total_files, downloaded, skipped,
        failed.
    """
    gse_id = gse_id.upper()
    os.makedirs(dest_dir, exist_ok=True)

    result: dict = {
        "gse_id": gse_id,
        "dry_run": dry_run,
        "metadata": {},
        "files": [],
        "total_files": 0,
        "downloaded": 0,
        "skipped": 0,
        "failed": 0,
        "total_size_bytes": 0,
    }

    # ── Phase 1: Metadata ────────────────────────────────────────
    if fetch_meta:
        if not quiet:
            print(f"\n{'=' * 60}")
            print(f"  GEO Download: {gse_id}")
            print(f"{'=' * 60}")

        try:
            meta = fetch_soft_metadata(gse_id)
            result["metadata"] = meta

            if not quiet:
                print(f"\n  [METADATA] {gse_id}_family.soft.gz")
                print(f"    Title:    {meta.get('title', '?')[:120]}")
                print(f"    Organism: {meta.get('organism', '?')}")
                print(f"    Platform: {meta.get('platform_title', '?')}")
                print(f"    PMID:     {', '.join(meta.get('pmid', [])) or 'none'}")
                print(f"    Samples:  {meta.get('n_samples', 0)}")
                if meta.get("is_superseries"):
                    print(f"    ⚠  SuperSeries detected — individual sub-series "
                          f"should be downloaded separately")

            # Save metadata cache
            meta_path = os.path.join(dest_dir, ".geo_meta.json")
            if not dry_run:
                with open(meta_path, "w", encoding="utf-8") as f:
                    json.dump(meta, f, indent=2, ensure_ascii=False, default=str)
        except Exception as exc:
            if not quiet:
                print(f"  [WARNING] Failed to fetch SOFT metadata: {exc}")
                print(f"    Continuing with file listing only...")

    # ── Phase 2: File listing ─────────────────────────────────────
    if not quiet:
        print(f"\n  [SUPPL] Listing files...")

    try:
        files = list_suppl_files(gse_id)
    except Exception as exc:
        print(f"  [ERROR] Failed to list supplementary files: {exc}")
        result["failed"] = -1
        return result

    result["total_files"] = len(files)
    result["total_size_bytes"] = sum(f.get("size_bytes", 0) for f in files)
    total_human = _human_size(result["total_size_bytes"])

    if not quiet:
        print(f"    {len(files)} file(s) found ({total_human} total)")
        for f in files:
            tag = "★ RAW" if f["is_raw_tar"] else "  "
            print(f"    {tag} {f['name'][:55]:55s} {f['size_human']:>10s}")

    if not files:
        print(f"  No supplementary files found for {gse_id}.")
        return result

    # ── Phase 3: Download ─────────────────────────────────────────
    if dry_run:
        print(f"\n  [DRY-RUN] Would download to: {dest_dir}/")
        for f in files:
            print(f"    → {f['name']}  ({f['size_human']})")
        return result

    if not quiet:
        print(f"\n  [DOWNLOAD]")

    for f in files:
        dest_path = os.path.join(dest_dir, f["name"])

        # Check existing
        if skip_existing and os.path.exists(dest_path):
            existing_size = os.path.getsize(dest_path)
            expected_size = f.get("size_bytes", 0)
            if expected_size and existing_size == expected_size:
                if not quiet:
                    print(f"    [SKIP] {f['name']} (exists, size OK)")
                result["files"].append({**f, "status": "skipped"})
                result["skipped"] += 1
                continue
            elif expected_size and existing_size != expected_size:
                if not quiet:
                    print(f"    [RESUME] {f['name']} (partial: "
                          f"{_human_size(existing_size)} / {f['size_human']})")

        if not quiet:
            print(f"    [{_human_size(f['size_bytes']):>8s}] {f['name']}")

        success = download_file(
            f["url"], dest_path,
            resume=True,
            show_progress=not quiet,
        )

        if success:
            result["files"].append({**f, "status": "downloaded"})
            result["downloaded"] += 1
        else:
            if not quiet:
                print(f"    [FAILED] {f['name']}")
            result["files"].append({**f, "status": "failed"})
            result["failed"] += 1

    # ── Phase 4: Summary ──────────────────────────────────────────
    if not quiet:
        downloaded_size = sum(
            f["size_bytes"] for f in result["files"]
            if f.get("status") == "downloaded"
        )
        print(f"\n  {'=' * 60}")
        print(f"  [SUMMARY] {gse_id}")
        print(f"  {'=' * 60}")
        print(f"    Downloaded: {result['downloaded']}  "
              f"(~{_human_size(downloaded_size)})")
        print(f"    Skipped:    {result['skipped']}")
        print(f"    Failed:     {result['failed']}")
        if result["failed"]:
            failed_names = [
                f["name"] for f in result["files"]
                if f.get("status") == "failed"
            ]
            print(f"    Failed files: {', '.join(failed_names)}")
        print(f"    Output:     {dest_dir}/")
        print()
        print(f"  Next step:")
        print(f"    python core/preprocess/preprocessor.py "
              f"--gse {gse_id}")

    return result


# ── Species name normalisation (NCBI taxon → common slug) ─────────

_SPECIES_TO_SLUG: dict[str, str] = {
    "mus musculus":         "mouse",
    "homo sapiens":         "human",
    "gallus gallus":       "chick",
    "danio rerio":          "zebrafish",
    "macaca fascicularis":  "macaque",
    "rattus norvegicus":    "rat",
    "macaca mulatta":       "macaque",
    "xenopus laevis":       "xenopus",
    "xenopus tropicalis":   "xenopus",
    "drosophila melanogaster": "fruit_fly",
    "caenorhabditis elegans":  "c_elegans",
}


def _normalise_species(organism: str) -> str:
    """Convert NCBI organism name to a short slug."""
    return _SPECIES_TO_SLUG.get(organism.lower(), organism.lower().replace(" ", "_"))


def _load_cached_meta(gse_dir: str) -> dict | None:
    """Read cached SOFT metadata from .geo_meta.json if it exists."""
    meta_path = os.path.join(gse_dir, ".geo_meta.json")
    if not os.path.isfile(meta_path):
        return None
    try:
        with open(meta_path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def enrich_dataset_from_soft(
    gse_id: str,
    soft_meta: dict,
    dry_run: bool = False,
) -> bool:
    """Enrich a registry dataset entry with SOFT-derived metadata.

    Only fills fields that are currently empty/null — never overwrites
    manually curated data. Creates paper-dataset links if PMIDs match
    existing papers.

    Args:
        gse_id: GEO accession ID.
        soft_meta: Parsed SOFT metadata dict.
        dry_run: Report only, don't write.

    Returns:
        True if any enrichment was applied.
    """
    try:
        from core.registry import (
            load_master_registry,
            save_master_registry,
            LinkRole,
            PaperDatasetLink,
        )
    except ImportError:
        log.warning("Cannot import core.registry — skipping enrichment")
        return False

    registry = load_master_registry()
    ds = registry.datasets.get(gse_id)
    if ds is None:
        log.info("Dataset %s not in registry — nothing to enrich", gse_id)
        return False

    changed = False

    # ── Fill species (empty → from SOFT organism) ─────────────────
    if not ds.species and soft_meta.get("organism"):
        species_slug = _normalise_species(soft_meta["organism"])
        ds.species = species_slug
        log.info("  species: '' → %s", species_slug)
        changed = True

    # ── Fill n_samples (None → from SOFT) ────────────────────────
    n = soft_meta.get("n_samples", 0)
    if ds.n_samples is None and n > 0:
        ds.n_samples = n
        log.info("  n_samples: None → %d", n)
        changed = True

    # ── Append paper_pmids from SOFT (merge, no duplicate) ───────
    soft_pmids = soft_meta.get("pmid", [])
    if soft_pmids:
        for pmid in soft_pmids:
            if pmid not in ds.paper_pmids:
                ds.paper_pmids.append(pmid)
                log.info("  paper_pmids: + %s", pmid)
                changed = True

    # ── Create links to existing papers ──────────────────────────
    for pmid in soft_pmids:
        paper = registry.get_paper(pmid)
        if paper is None:
            continue
        existing_links = registry.get_dataset_links(paper.paper_id)
        if not any(ln[0] == gse_id for ln in existing_links):
            registry.links.append(PaperDatasetLink(
                paper_id=paper.paper_id,
                dataset_id=gse_id,
                role=LinkRole.PRIMARY,
            ))
            log.info("  link: %s ↔ %s", paper.paper_id, gse_id)
            changed = True

    # ── Append summary excerpt to notes (if notes empty) ────────
    if not ds.notes and soft_meta.get("summary"):
        excerpt = soft_meta["summary"][:200].strip()
        ds.notes = excerpt
        log.info("  notes: filled from SOFT summary")
        changed = True

    if changed and not dry_run:
        save_master_registry(registry)
        log.info("Registry enriched for %s", gse_id)

    return changed


def update_registry_after_download(
    gse_id: str,
    data_dir: str | None = None,
    dry_run: bool = False,
) -> bool:
    """Update Master Registry after successful download.

    1. Changes dataset status to ``data_downloaded``.
    2. Enriches dataset metadata from cached SOFT ``.geo_meta.json``
       (species, n_samples, paper_pmids, links).

    Args:
        gse_id: GEO accession ID.
        data_dir: Directory where data was downloaded (to find .geo_meta.json).
            Defaults to ``$FUXI_DATA_ROOT/{gse_id}``.
        dry_run: Report only, don't write.

    Returns:
        ``True`` if the registry was updated.
    """
    try:
        from core.registry import (
            load_master_registry,
            save_master_registry,
            DatasetStatus,
        )
    except ImportError:
        log.warning("Cannot import core.registry — skipping status update")
        return False

    gse_id = gse_id.upper()
    registry = load_master_registry()
    ds = registry.datasets.get(gse_id)

    if ds is None:
        log.info("Dataset %s not in registry — nothing to update", gse_id)
        return False

    updated = False

    # Step 1: Update status
    if ds.status != DatasetStatus.DATA_DOWNLOADED:
        ds.status = DatasetStatus.DATA_DOWNLOADED
        log.info("Updated registry: %s → data_downloaded", gse_id)
        updated = True

    if updated and not dry_run:
        save_master_registry(registry)

    # Step 2: Enrich from SOFT metadata
    if data_dir is None:
        data_root = os.environ.get("FUXI_DATA_ROOT", "")
        if data_root:
            data_dir = os.path.join(data_root, gse_id)

    if data_dir:
        soft_meta = _load_cached_meta(data_dir)
        if soft_meta:
            enriched = enrich_dataset_from_soft(gse_id, soft_meta, dry_run)
            updated = updated or enriched

    return updated


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fuxi GEO Downloader — Download datasets from NCBI GEO",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python core/geo_downloader.py --gse GSE118614
  python core/geo_downloader.py --gse GSE118614 --dry-run
  python core/geo_downloader.py --gse GSE118614 --skip-soft --force
  python core/geo_downloader.py --gse GSE118614 --data-root /mnt/e/data
""",
    )
    parser.add_argument(
        "--gse", type=str, required=True,
        help="GEO accession ID (e.g., GSE118614)",
    )
    parser.add_argument(
        "--data-root", type=str, default=None,
        help="Override FUXI_DATA_ROOT (default: from environment variable)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show what would be downloaded without actually downloading",
    )
    parser.add_argument(
        "--skip-soft", action="store_true",
        help="Skip SOFT metadata fetch (use when metadata already cached)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download files even if they already exist",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Minimal output",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    gse_id = args.gse.upper()

    # Resolve data root
    data_root = args.data_root
    if data_root is None:
        data_root = os.environ.get("FUXI_DATA_ROOT", "")
    if not data_root:
        print("[ERROR] FUXI_DATA_ROOT is not set.", file=sys.stderr)
        print("        Set it in .env or pass --data-root.",
              file=sys.stderr)
        sys.exit(1)

    dest_dir = os.path.join(data_root, gse_id)

    # Quick pre-check: does downloader exist?
    try:
        _get_downloader_name()
    except RuntimeError as exc:
        if not args.dry_run:
            print(f"[ERROR] {exc}", file=sys.stderr)
            sys.exit(1)

    result = download_gse(
        gse_id=gse_id,
        dest_dir=dest_dir,
        dry_run=args.dry_run,
        skip_existing=not args.force,
        fetch_meta=not args.skip_soft,
        quiet=args.quiet,
    )

    # Update registry on success
    if not args.dry_run and result["failed"] == 0:
        update_registry_after_download(gse_id, dry_run=False)

    if result["failed"] > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
