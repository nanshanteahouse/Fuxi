#!/usr/bin/env python3
"""
superseries_detector.py — SuperSeries detection
================================================

Detects whether a GEO dataset directory contains a SuperSeries
(collection of sub-series, potentially multi-modality or multi-species).

Three detection strategies, tried in order:

  1. **Directory structure**: Subdirectories named ``GSE\\d+`` after
     extraction indicate sub-series.

  2. **GEO Series Matrix parsing**: ``*_series_matrix.txt(.gz)`` files
     often contain ``!Series_relation = SuperSeries of: GSE...`` in the
     first few hundred lines.

  3. **NCBI E-utilities API** (opt-in via ``--query-ncbi``): Queries
     ``esummary.fcgi?db=gds&id={accession}`` for the entry type field.
     SuperSeries entries return ``entrytype: "GSE"`` with child links.

Strategy 3 is only executed if the caller explicitly enables it, because
it requires internet access and adds ~2s latency.  Results from the NCBI
API are cached for 24 hours in ``~/.fuxi/cache/ncbi_{accession}.json``.
"""

import os
import re
import gzip
import json
import time
from datetime import datetime, timezone
from typing import Optional
from urllib.request import Request, urlopen
from urllib.error import URLError


# ── Constants ─────────────────────────────────────────────────────────

_GSE_DIR_RE = re.compile(r'^GSE\d+$')
_SERIES_RELATION_RE = re.compile(
    r'!Series_relation\s*=\s*Super[Ss]eries\s+(?:of|is|:)\s*(.*)',
    re.IGNORECASE,
)
_SERIES_TYPE_RE = re.compile(
    r'!Series_type\s*=\s*(.*)',
    re.IGNORECASE,
)
_CACHE_DIR = os.path.expanduser('~/.fuxi/cache')
_CACHE_TTL_SECONDS = 86400  # 24 hours


# ── Strategy 1: Directory structure ───────────────────────────────────

def detect_subseries_dirs(root_dir: str) -> list[str]:
    """Scan *root_dir* for immediate subdirectories named like GEO accessions.

    Returns a list of directory names (not full paths).
    """
    subseries = []
    try:
        for entry in sorted(os.listdir(root_dir)):
            entry_path = os.path.join(root_dir, entry)
            if os.path.isdir(entry_path) and _GSE_DIR_RE.match(entry):
                subseries.append(entry)
    except OSError:
        pass
    return subseries


# ── Strategy 2: Series Matrix parsing ─────────────────────────────────

def _read_first_lines(filepath: str, n: int = 500) -> str:
    """Read the first *n* lines of a file (handles .gz transparently).

    Returns '' if the file cannot be read or is not valid gzip.
    """
    if filepath.lower().endswith('.gz'):
        try:
            with gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') as f:
                lines = []
                for i, line in enumerate(f):
                    if i >= n:
                        break
                    lines.append(line)
                return '\n'.join(lines)
        except (gzip.BadGzipFile, OSError, EOFError):
            return ''
        except Exception:
            return ''
    else:
        try:
            with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
                lines = []
                for i, line in enumerate(f):
                    if i >= n:
                        break
                    lines.append(line)
                return '\n'.join(lines)
        except (UnicodeDecodeError, OSError):
            return ''


def parse_series_matrix(filepath: str) -> Optional[dict]:
    """Parse a GEO Series Matrix file for SuperSeries metadata.

    Returns:
        {'is_superseries': bool,
         'series_type': str or None,
         'child_accessions': [str]}  or  None if parsing fails.
    """
    text = _read_first_lines(filepath, n=500)
    if not text:
        return None

    is_super = False
    child_accessions: list[str] = []
    series_type: Optional[str] = None

    for line in text.split('\n'):
        m = _SERIES_RELATION_RE.search(line)
        if m:
            is_super = True
            # Extract accession IDs from the rest of the line
            rest = m.group(1)
            child_accessions.extend(re.findall(r'GSE\d+', rest))

        m2 = _SERIES_TYPE_RE.search(line)
        if m2:
            series_type = m2.group(1).strip().strip('"')

    if is_super or series_type == 'SuperSeries':
        return {
            'is_superseries': True,
            'series_type': series_type,
            'child_accessions': child_accessions,
        }

    return None


def detect_from_series_matrix_files(file_list: list[str]) -> Optional[dict]:
    """Try to detect SuperSeries from any Series Matrix files in *file_list*.

    Returns the first successful parse result, or None.
    """
    from .format_detector import is_geo_series_matrix

    for f in file_list:
        if is_geo_series_matrix(f):
            result = parse_series_matrix(f)
            if result and result.get('is_superseries'):
                return result
    return None


# ── Strategy 3: NCBI E-utilities API ──────────────────────────────────

def _ncbi_cache_path(accession: str) -> str:
    """Return the cache file path for an NCBI API response."""
    os.makedirs(_CACHE_DIR, exist_ok=True)
    return os.path.join(_CACHE_DIR, f'ncbi_{accession}.json')


def _ncbi_cached(accession: str) -> Optional[dict]:
    """Load cached NCBI response if it exists and is still fresh."""
    cache_path = _ncbi_cache_path(accession)
    if not os.path.exists(cache_path):
        return None
    try:
        with open(cache_path, 'r') as f:
            data = json.load(f)
        timestamp = data.get('_cached_at', '')
        if timestamp:
            cached_time = datetime.fromisoformat(timestamp)
            age = (datetime.now(timezone.utc) - cached_time).total_seconds()
            if age < _CACHE_TTL_SECONDS:
                return data
    except (json.JSONDecodeError, KeyError, ValueError):
        pass
    return None


def _ncbi_save_cache(accession: str, data: dict) -> None:
    """Save an NCBI API response to disk cache."""
    cache_path = _ncbi_cache_path(accession)
    data['_cached_at'] = datetime.now(timezone.utc).isoformat()
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, 'w') as f:
        json.dump(data, f, indent=2)


def _ncbi_esummary(accession: str) -> Optional[dict]:
    """Query NCBI E-utilities esummary for a GEO accession.

    Returns the JSON result dict, or None on failure.
    """
    url = (
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        f'?db=gds&id={accession}&retmode=json'
    )
    req = Request(url)
    req.add_header('User-Agent', 'Fuxi/1.0 (single-cell pipeline; academic use)')
    try:
        with urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        return data
    except (URLError, json.JSONDecodeError, OSError) as e:
        return {'_error': str(e)}


def _ncbi_elink_superseries(accession: str) -> Optional[list[str]]:
    """Query NCBI E-utilities elink to get child accessions of a SuperSeries.

    Returns a list of child accession IDs, or None on failure.
    """
    url = (
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
        f'?dbfrom=gds&db=gds&linkname=gds_gds_superseries&id={accession}&retmode=json'
    )
    req = Request(url)
    req.add_header('User-Agent', 'Fuxi/1.0 (single-cell pipeline; academic use)')
    try:
        with urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        ids = []
        for linkset in data.get('linksets', []):
            for link in linkset.get('linksetdbs', []):
                ids.extend(link.get('links', []))
        return ids
    except (URLError, json.JSONDecodeError, OSError):
        return None


def query_ncbi(accession: str) -> Optional[dict]:
    """Query NCBI E-utilities for SuperSeries metadata.

    Caches results for 24 hours.

    Args:
        accession: GEO accession ID (e.g. 'GSE12345').

    Returns:
        {'is_superseries': bool,
         'entry_type': str,
         'child_accessions': [str],
         'title': str,
         'summary': str}  or  None on failure.
    """
    # Check cache first
    cached = _ncbi_cached(accession)
    if cached:
        return cached

    # Rate-limit: NCBI allows 3 req/sec without API key
    time.sleep(0.34)

    esummary = _ncbi_esummary(accession)
    if esummary is None or '_error' in esummary:
        return None

    result_data = esummary.get('result', {})
    entry = result_data.get(accession, {})

    entry_type = entry.get('entrytype', '')
    is_super = entry_type.upper() in ('GSE', 'SUPERSERIES')

    child_accessions: list[str] = []
    if is_super:
        time.sleep(0.34)
        children = _ncbi_elink_superseries(accession)
        if children:
            child_accessions.extend(children)

    result = {
        'is_superseries': is_super,
        'entry_type': entry_type,
        'child_accessions': child_accessions,
        'title': entry.get('title', ''),
        'summary': entry.get('summary', ''),
        'species': entry.get('taxon', ''),
    }

    _ncbi_save_cache(accession, result)
    return result


# ── Top-level detection orchestrator ──────────────────────────────────

def detect_superseries(root_dir: str,
                       file_list: list[str],
                       gse_id: str,
                       query_ncbi_flag: bool = False) -> dict:
    """Run all SuperSeries detection strategies and return a combined result.

    Args:
        root_dir:   Top-level dataset directory.
        file_list:  All file paths under root_dir (after extraction).
        gse_id:     GEO accession ID (e.g. 'GSE12345').
        query_ncbi_flag:  If True, also query NCBI E-utilities.

    Returns:
        {
            'is_superseries': bool,
            'detected_by': str or None,   # 'dirs' | 'series_matrix' | 'ncbi' | None
            'subseries_dirs': [str],      # directory names under root_dir
            'child_accessions': [str],    # accession IDs (from matrix or NCBI)
            'title': '',
            'summary': '',
        }
    """
    result: dict = {
        'is_superseries': False,
        'detected_by': None,
        'subseries_dirs': [],
        'child_accessions': [],
        'title': '',
        'summary': '',
    }

    # Strategy 1: Directory structure
    subseries_dirs = detect_subseries_dirs(root_dir)
    if subseries_dirs:
        result['is_superseries'] = True
        result['detected_by'] = 'dirs'
        result['subseries_dirs'] = subseries_dirs

    # Strategy 2: Series Matrix files
    matrix_result = detect_from_series_matrix_files(file_list)
    if matrix_result:
        result['is_superseries'] = True
        if not result['detected_by']:
            result['detected_by'] = 'series_matrix'
        existing = set(result['child_accessions'])
        for acc in matrix_result.get('child_accessions', []):
            if acc not in existing:
                result['child_accessions'].append(acc)

    # Strategy 3: NCBI API (opt-in)
    if query_ncbi_flag:
        ncbi_result = query_ncbi(gse_id)
        if ncbi_result:
            if ncbi_result.get('is_superseries'):
                result['is_superseries'] = True
                if not result['detected_by']:
                    result['detected_by'] = 'ncbi'
            result['title'] = result['title'] or ncbi_result.get('title', '')
            result['summary'] = result['summary'] or ncbi_result.get('summary', '')
            existing = set(result['child_accessions'])
            for acc in ncbi_result.get('child_accessions', []):
                if acc not in existing:
                    result['child_accessions'].append(acc)

    return result
