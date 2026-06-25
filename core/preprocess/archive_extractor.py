#!/usr/bin/env python3
"""
archive_extractor.py — Safe, idempotent archive extraction
===========================================================

Extracts the following archive formats using only stdlib:
  - .tar.gz  / .tgz       → tarfile (gzip)
  - .tar.bz2               → tarfile (bzip2)
  - .tar.xz                → tarfile (xz / lzma)
  - .tar                   → tarfile (uncompressed)
  - .zip                   → zipfile
  - .gz (single file)      → gzip (decompress to sibling without .gz)
  - .bz2 (single file)     → bzip2 (decompress to sibling without .bz2)

Idempotency via sentinel files:
  After extraction, an ``.extracted_{archive_name}.json`` file is written
  beside the archive containing:
    - archive path (relative)
    - SHA-256 checksum at extraction time
    - timestamp
    - file count
  On subsequent runs, if the sentinel exists and checksum matches, extraction
  is skipped.

Safety:
  - Extracts to a sibling directory named ``{archive_stem}_extracted/``.
  - Path traversal protection via ``tarfile.extractall(filter='data')``
    on Python 3.12+; manual path check for older Pythons.
  - Never overwrites existing files outside the extraction target.
  - Incremental SHA-256 (64KB chunks) — safe for large archives.
"""

import os
import sys
import json
import hashlib
import shutil
import tarfile
import zipfile
import gzip
import bz2
from datetime import datetime
from typing import Optional


# ── Checksum helpers ──────────────────────────────────────────────────

_CHUNK_SIZE = 64 * 1024  # 64 KB


def compute_sha256(filepath: str) -> str:
    """Compute SHA-256 hash of a file (streaming, safe for large files)."""
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while True:
            chunk = f.read(_CHUNK_SIZE)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _sentinel_path(archive_path: str) -> str:
    """Return the path of the sentinel file for an archive."""
    archive_name = os.path.basename(archive_path)
    parent = os.path.dirname(archive_path) or '.'
    return os.path.join(parent, f'.extracted_{archive_name}.json')


def _should_skip(archive_path: str, verbose: bool = False) -> bool:
    """Return True if the archive has already been extracted (checksum match)."""
    sentinel = _sentinel_path(archive_path)
    if not os.path.exists(sentinel):
        return False
    try:
        with open(sentinel, 'r') as f:
            data = json.load(f)
        current = compute_sha256(archive_path)
        if data.get('checksum') == current:
            return True
        if verbose:
            print(f"  Checksum changed for {os.path.basename(archive_path)} — re-extracting")
    except (json.JSONDecodeError, KeyError, FileNotFoundError):
        pass
    return False


def _write_sentinel(archive_path: str, extract_dir: str, file_count: int) -> None:
    """Write extraction sentinel file."""
    sentinel = _sentinel_path(archive_path)
    data = {
        'archive_path': os.path.basename(archive_path),
        'checksum': compute_sha256(archive_path),
        'extracted_at': datetime.now().isoformat(),
        'extract_dir': os.path.basename(extract_dir),
        'file_count': file_count,
    }
    with open(sentinel, 'w') as f:
        json.dump(data, f, indent=2)


# ── Extraction destination ────────────────────────────────────────────

def _extract_dest_dir(archive_path: str) -> str:
    """Compute the sibling extraction directory for an archive.

    Example: 'GSE12345_RAW.tar.gz' → 'GSE12345/RAW_extracted/'
    """
    parent = os.path.dirname(archive_path) or '.'
    stem = os.path.basename(archive_path)
    # Remove all known archive extensions
    for ext in ('.tar.gz', '.tgz', '.tar.bz2', '.tar.xz', '.tar', '.zip', '.gz', '.bz2'):
        if stem.lower().endswith(ext):
            stem = stem[:-len(ext)]
            break
    return os.path.join(parent, f'{stem}_extracted')


# ── Path traversal protection ─────────────────────────────────────────

def _is_safe_extraction(member_path: str, dest_dir: str) -> bool:
    """Return True if *member_path* resolves inside *dest_dir*."""
    dest_abs = os.path.abspath(dest_dir)
    member_abs = os.path.abspath(os.path.join(dest_dir, member_path))
    return os.path.commonpath([dest_abs, member_abs]) == dest_abs


# ── Tar extraction ────────────────────────────────────────────────────

def _extract_tar(archive_path: str, dest_dir: str, mode: str) -> int:
    """Extract a tar archive. Returns count of extracted files."""
    os.makedirs(dest_dir, exist_ok=True)
    count = 0

    # Python 3.12+ has filter='data'; older versions need manual check
    if sys.version_info >= (3, 12):
        with tarfile.open(archive_path, mode=mode) as tf:
            tf.extractall(path=dest_dir, filter='data')
            count = len(tf.getmembers())
    else:
        with tarfile.open(archive_path, mode=mode) as tf:
            for member in tf.getmembers():
                # Path traversal check
                if member.name.startswith('/') or '..' in member.name.split('/'):
                    print(f"  [WARNING] Skipping unsafe path: {member.name}")
                    continue
                try:
                    tf.extract(member, path=dest_dir)
                    count += 1
                except Exception as e:
                    print(f"  [WARNING] Failed to extract {member.name}: {e}")

    return count


# ── Zip extraction ────────────────────────────────────────────────────

def _extract_zip(archive_path: str, dest_dir: str) -> int:
    """Extract a zip archive. Returns count of extracted files."""
    os.makedirs(dest_dir, exist_ok=True)
    count = 0
    with zipfile.ZipFile(archive_path, 'r') as zf:
        for member in zf.infolist():
            # Path traversal check
            if member.filename.startswith('/') or '..' in member.filename.split('/'):
                print(f"  [WARNING] Skipping unsafe path: {member.filename}")
                continue
            try:
                zf.extract(member, path=dest_dir)
                count += 1
            except Exception as e:
                print(f"  [WARNING] Failed to extract {member.filename}: {e}")
    return count


# ── Single-file gzip / bzip2 ──────────────────────────────────────────

def _extract_gzip_single(archive_path: str, dest_dir: str) -> int:
    """Decompress a single .gz file into dest_dir (without .gz suffix)."""
    os.makedirs(dest_dir, exist_ok=True)
    basename = os.path.basename(archive_path)
    if basename.lower().endswith('.gz'):
        out_name = basename[:-3]
    else:
        out_name = basename + '.decompressed'
    out_path = os.path.join(dest_dir, out_name)
    try:
        with gzip.open(archive_path, 'rb') as f_in:
            with open(out_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except gzip.BadGzipFile as e:
        raise ValueError(f"Not a valid gzip file: {archive_path}") from e
    return 1


def _extract_bzip2_single(archive_path: str, dest_dir: str) -> int:
    """Decompress a single .bz2 file into dest_dir (without .bz2 suffix)."""
    os.makedirs(dest_dir, exist_ok=True)
    basename = os.path.basename(archive_path)
    if basename.lower().endswith('.bz2'):
        out_name = basename[:-4]
    else:
        out_name = basename + '.decompressed'
    out_path = os.path.join(dest_dir, out_name)
    try:
        with bz2.open(archive_path, 'rb') as f_in:
            with open(out_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        raise ValueError(f"Not a valid bzip2 file: {archive_path}") from e
    return 1


# ── Main extraction dispatcher ────────────────────────────────────────

EXTRACTION_MODES: dict[str, str] = {
    'tar.gz':  'r:gz',
    'tar.bz2': 'r:bz2',
    'tar.xz':  'r:xz',
    'tar':     'r',
}


def extract_archive(archive_path: str, fmt: str,
                    dest_dir: Optional[str] = None,
                    verbose: bool = False) -> dict:
    """Extract a single archive.

    Args:
        archive_path: Absolute path to the archive file.
        fmt: Archive format name (from format_detector.ARCHIVE_PATTERNS).
        dest_dir: Override extraction destination directory.
        verbose: Print progress messages.

    Returns:
        {'status': 'skipped'|'extracted'|'error',
         'dest_dir': str,
         'file_count': int,
         'error': str or None}
    """
    if dest_dir is None:
        dest_dir = _extract_dest_dir(archive_path)

    # Check sentinel
    if _should_skip(archive_path, verbose=verbose):
        if verbose:
            print(f"  Skipping (already extracted): {os.path.basename(archive_path)}")
        # Read file_count from sentinel
        sentinel = _sentinel_path(archive_path)
        with open(sentinel, 'r') as f:
            data = json.load(f)
        return {
            'status': 'skipped',
            'dest_dir': os.path.dirname(archive_path) or '.',
            'file_count': data.get('file_count', 0),
            'error': None,
        }

    if verbose:
        print(f"  Extracting: {os.path.basename(archive_path)} → {os.path.basename(dest_dir)}/")

    try:
        if fmt in EXTRACTION_MODES:
            mode = EXTRACTION_MODES[fmt]
            count = _extract_tar(archive_path, dest_dir, mode)
        elif fmt == 'zip':
            count = _extract_zip(archive_path, dest_dir)
        elif fmt == 'gzip_single':
            count = _extract_gzip_single(archive_path, dest_dir)
        elif fmt == 'bzip2_single':
            count = _extract_bzip2_single(archive_path, dest_dir)
        else:
            return {
                'status': 'error',
                'dest_dir': dest_dir,
                'file_count': 0,
                'error': f'Unknown format: {fmt}',
            }

        _write_sentinel(archive_path, dest_dir, count)

        result = {
            'status': 'extracted',
            'dest_dir': dest_dir,
            'file_count': count,
            'error': None,
        }
        if verbose:
            print(f"    Done: {count} files extracted")
        return result

    except Exception as e:
        return {
            'status': 'error',
            'dest_dir': dest_dir,
            'file_count': 0,
            'error': str(e),
        }


def extract_all_archives(archive_list: list[tuple[str, str]],
                         verbose: bool = False) -> list[dict]:
    """Extract all archives in *archive_list*.

    Args:
        archive_list: List of (path, format_name) tuples.
        verbose: Print progress messages.

    Returns:
        List of result dicts (see extract_archive).
    """
    results = []
    for path, fmt in archive_list:
        result = extract_archive(path, fmt, verbose=verbose)
        results.append(result)
    return results


def collect_file_tree(root_dir: str) -> list[str]:
    """Recursively collect all file paths under *root_dir*.

    Excludes:
      - Hidden files / directories (start with '.')
      - Sentinel files (start with '.extracted_')
    """
    all_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Exclude hidden directories and extraction directories
        dirnames[:] = [d for d in dirnames
                       if not d.startswith('.')]
        for fname in filenames:
            if fname.startswith('.') or fname.startswith('.extracted_'):
                continue
            all_files.append(os.path.join(dirpath, fname))
    return all_files
