#!/usr/bin/env python3
"""
format_detector.py — Extended file format detection for preprocessing
=====================================================================

Extends the modality-pattern table in core/dataset_detector.py with:
  - Archive format patterns (tar.gz, zip, gz, bz2, etc.)
  - Unsupported format patterns (7z, rar)
  - 10X directory structure co-occurrence detection
  - Improved sample grouping by stripping known suffixes

Design: all functions are pure — they take file lists / directory trees
and return classification results. No side effects. No I/O beyond what
the caller provides.
"""

import os
import re
from collections import defaultdict
from typing import Optional

# ── Archive format patterns ─────────────────────────────────────────
# Ordered so that compound extensions match before their components
# (e.g. .tar.gz before .gz alone).
ARCHIVE_PATTERNS = [
    (re.compile(r'\.tar\.gz$',  re.IGNORECASE), 'tar.gz'),
    (re.compile(r'\.tgz$',      re.IGNORECASE), 'tar.gz'),
    (re.compile(r'\.tar\.bz2$', re.IGNORECASE), 'tar.bz2'),
    (re.compile(r'\.tar\.xz$',  re.IGNORECASE), 'tar.xz'),
    (re.compile(r'\.tar$',      re.IGNORECASE), 'tar'),
    (re.compile(r'\.zip$',      re.IGNORECASE), 'zip'),
    (re.compile(r'\.gz$',       re.IGNORECASE), 'gzip_single'),
    (re.compile(r'\.bz2$',      re.IGNORECASE), 'bzip2_single'),
]

# ── Unsupported format patterns (cannot handle with stdlib) ─────────
UNSUPPORTED_PATTERNS = [
    re.compile(r'\.7z$',   re.IGNORECASE),
    re.compile(r'\.rar$',  re.IGNORECASE),
    re.compile(r'\.tar\.Z$', re.IGNORECASE),
]

# ── 10X MTX directory co-occurrence set ─────────────────────────────
# A directory containing >=N_MATCH of these files is classified as 10X MTX.
# We check whether basenames *end with* any of these patterns (handles
# prefixed filenames like "<GEO_ID>_retina_aggr_10_matrix.mtx.gz").
TENX_MTX_SUFFIXES = [
    '_matrix.mtx.gz', '_matrix.mtx',
    '_counts.mtx.gz', '_counts.mtx',
    '_barcodes.tsv.gz', '_barcodes.tsv', '_barcodes.csv.gz', '_barcodes.csv',
    '_features.tsv.gz', '_features.tsv', '_features.csv.gz', '_features.csv',
    '_genes.tsv.gz', '_genes.tsv', '_genes.csv.gz', '_genes.csv',
    '_gene_names.txt.gz', '_gene_names.txt',
    '_sample_annotations.tsv.gz', '_sample_annotations.txt.gz',
]
TENX_MTX_MIN_MATCH = 2

# ── TSV-based count matrix identifiers ───────────────────────────────
# Datasets where each GSM sample is a single cell-by-gene TSV/CSV file.
CSV_MATRIX_EXTENSIONS = {'.csv', '.csv.gz', '.tsv', '.tsv.gz', '.txt', '.txt.gz'}

# ── 10X H5 identifiers ──────────────────────────────────────────────
# Directory-level patterns: any file whose basename *ends with* one of
# these substrings signals 10X HDF5 format.
TENX_H5_SUFFIXES = ['_filtered_feature_bc_matrix.h5', '_raw_feature_bc_matrix.h5']

# ── 10X ATAC peak identifiers ───────────────────────────────────────
TENX_PEAK_SUFFIXES = ['_filtered_peak_bc_matrix.h5']
TENX_PEAK_FILES = {'peaks.bed', 'singlecell.csv'}
TENX_PEAK_MIN_MATCH = 2

# ── ATAC fragment identifiers ───────────────────────────────────────
FRAGMENT_SUFFIXES = ['_fragments.tsv.gz', '_fragments.tsv']

# ── Known suffixes for stripping when deriving sample names ─────────
KNOWN_SUFFIXES = [
    '_filtered_feature_bc_matrix.h5',
    '_raw_feature_bc_matrix.h5',
    '_filtered_peak_bc_matrix.h5',
    '_fragments.tsv.gz',
    '_fragments.tsv',
    '_counts.csv.gz',
    '_counts.csv',
    '_count_mat.csv.gz',
    '_TPM.csv.gz',
    '_atac_fragments.tsv.gz',
    '_barcodes.tsv.gz',
    '_barcodes.tsv',
    '_features.tsv.gz',
    '_features.tsv',
    '_genes.tsv.gz',
    '_genes.tsv',
    '_matrix.mtx.gz',
    '_matrix.mtx',
    '_metadata.csv',
    '_metadata.csv.gz',
    '_barcodes.csv',
    '_barcodes.csv.gz',
    '_series_matrix.txt.gz',
    '_series_matrix.txt',
]

# ── GEO-specific patterns ───────────────────────────────────────────
GEO_SERIES_MATRIX_RE = re.compile(r'GSE\d+_series_matrix', re.IGNORECASE)
GSE_DIR_RE = re.compile(r'^GSE\d+$')

# ── Species keyword hints in filenames ──────────────────────────────
SPECIES_KEYWORDS = [
    (['human', 'homo_sapiens', 'Hs_', 'GRCh38', 'hg38', 'hgnc'], 'homo_sapiens'),
    (['mouse', 'mus_musculus', 'Mm_', 'GRCm38', 'GRCm39', 'mm10', 'mgi'], 'mus_musculus'),
    (['rat', 'rattus'], 'rattus_norvegicus'),
    (['zebrafish', 'danio', 'DRerio'], 'danio_rerio'),
    (['cow', 'bos_taurus', 'cattle'], 'bos_taurus'),
    (['pig', 'sus_scrofa'], 'sus_scrofa'),
    (['macaque', 'rhesus', 'macaca'], 'macaca_mulatta'),
    (['marmoset', 'callithrix'], 'callithrix_jacchus'),
    (['ferret', 'mustela'], 'mustela_putorius_furo'),
    (['sheep', 'ovis_aries'], 'ovis_aries'),
    (['chicken', 'gallus'], 'gallus_gallus'),
    (['xenopus', 'tropicalis'], 'xenopus_tropicalis'),
    (['drosophila', 'melanogaster', 'dm6'], 'drosophila_melanogaster'),
    (['c_elegans', 'elegans', 'ce11', 'WBcel'], 'caenorhabditis_elegans'),
]


# ── Known compound suffixes for gzip data files that are NOT archives ──
# These are pipeline-readable compressed data files, not general-purpose
# archives to extract.  Exclude them from the 'gzip_single' archive category.
_GZ_DATA_SUFFIXES = {
    '.mtx.gz', '.tsv.gz', '.csv.gz', '.txt.gz', '.bed.gz', '.narrowPeak.gz',
    '.broadPeak.gz', '.gff.gz', '.gtf.gz', '.fastq.gz', '.fq.gz',
    '.tbi.gz', '.h5.gz',
    '.bam.gz',  # questionable but encountered
}


def is_archive(filepath: str) -> Optional[str]:
    """Return archive format name if *filepath* matches an archive pattern, else None.

    ``.gz`` files whose stem is a known data format (``.mtx``, ``.tsv``, ``.csv``,
    ``.txt``, ``.bed``, ``.gff``, ``.gtf``, ``.fastq``) are NOT classified as
    archives — the pipeline can read them natively.
    """
    basename = os.path.basename(filepath)
    for pattern, fmt in ARCHIVE_PATTERNS:
        if pattern.search(basename):
            # Exclude gzip-compressed data files from archive classification
            if fmt == 'gzip_single':
                basename_lower = basename.lower()
                for data_sfx in _GZ_DATA_SUFFIXES:
                    if basename_lower.endswith(data_sfx):
                        return None
            return fmt
    return None


def is_unsupported(filepath: str) -> Optional[str]:
    """Return unsupported format description if *filepath* matches, else None."""
    basename = os.path.basename(filepath)
    for pattern in UNSUPPORTED_PATTERNS:
        if pattern.search(basename):
            return pattern.pattern
    return None


def detect_10x_mtx_directory(dir_files: list[str]) -> bool:
    """Return True if *dir_files* contain enough 10X MTX markers."""
    basenames = [os.path.basename(f).lower() for f in dir_files]
    hits = 0
    for b in basenames:
        for sfx in TENX_MTX_SUFFIXES:
            if b.endswith(sfx.lower()):
                hits += 1
                break
    return hits >= TENX_MTX_MIN_MATCH


def _basename_ends_with_any(filename: str, suffixes: list[str]) -> bool:
    """Return True if *filename* ends with any of the given suffixes."""
    lower = filename.lower()
    return any(lower.endswith(sfx.lower()) for sfx in suffixes)


def detect_10x_h5_directory(dir_files: list[str]) -> bool:
    """Return True if *dir_files* contain a 10X HDF5 file."""
    basenames = [os.path.basename(f) for f in dir_files]
    for b in basenames:
        for sfx in TENX_H5_SUFFIXES:
            if b.lower().endswith(sfx.lower()):
                return True
    return False


def detect_10x_peak_directory(dir_files: list[str]) -> bool:
    """Return True if *dir_files* contain 10X ATAC peak markers."""
    basenames = {os.path.basename(f) for f in dir_files}
    # Check suffixed peak H5 files
    peak_h5_hits = sum(1 for b in basenames
                       for sfx in TENX_PEAK_SUFFIXES
                       if b.lower().endswith(sfx.lower()))
    other_hits = len(basenames & TENX_PEAK_FILES)
    return (peak_h5_hits + other_hits) >= TENX_PEAK_MIN_MATCH


def detect_fragment_file(dir_files: list[str]) -> bool:
    """Return True if *dir_files* contain ATAC fragment files."""
    basenames = [os.path.basename(f) for f in dir_files]
    for b in basenames:
        for sfx in FRAGMENT_SUFFIXES:
            if b.lower().endswith(sfx.lower()):
                return True
    return False


def is_geo_series_matrix(filepath: str) -> bool:
    """Return True if the file looks like a GEO Series Matrix file."""
    return bool(GEO_SERIES_MATRIX_RE.search(os.path.basename(filepath)))


def is_gse_directory(dirname: str) -> bool:
    """Return True if *dirname* matches the GSE accession pattern."""
    return bool(GSE_DIR_RE.match(dirname))


def classify_files_by_format(file_list: list[str]) -> dict:
    """Classify a list of file paths by detected format.

    Returns:
        {
            'archives':           [(path, format_name), ...],
            'unsupported':        [path, ...],
            'tenx_mtx_dirs':      {dirpath: [files], ...},
            'tenx_h5_dirs':       {dirpath: [files], ...},
            'tenx_peak_dirs':     {dirpath: [files], ...},
            'fragment_dirs':      {dirpath: [files], ...},
            'h5ad_files':         [path, ...],
            'csv_files':          [path, ...],
            'series_matrix_files': [path, ...],
            'metadata_files':     [path, ...],
            'unmatched':          [path, ...],
        }
    """
    result = {
        'archives': [],
        'unsupported': [],
        'tenx_mtx_dirs': defaultdict(list),
        'tenx_h5_dirs': defaultdict(list),
        'tenx_peak_dirs': defaultdict(list),
        'fragment_dirs': defaultdict(list),
        'h5ad_files': [],
        'csv_files': [],
        'series_matrix_files': [],
        'metadata_files': [],
        'unmatched': [],
    }

    # First pass: identify archives and unsupported files
    non_archive_files = []
    for f in file_list:
        arc_fmt = is_archive(f)
        if arc_fmt:
            result['archives'].append((f, arc_fmt))
            continue
        unsup = is_unsupported(f)
        if unsup:
            result['unsupported'].append(f)
            continue
        non_archive_files.append(f)

    # Second pass: group by directory for co-occurrence detection
    dir_groups = defaultdict(list)
    for f in non_archive_files:
        d = os.path.dirname(f) or '.'
        dir_groups[d].append(f)

    for dirpath, files in dir_groups.items():
        if detect_10x_mtx_directory(files):
            result['tenx_mtx_dirs'][dirpath].extend(files)
            continue
        if detect_10x_h5_directory(files):
            result['tenx_h5_dirs'][dirpath].extend(files)
            continue
        if detect_10x_peak_directory(files):
            result['tenx_peak_dirs'][dirpath].extend(files)
            continue
        if detect_fragment_file(files):
            result['fragment_dirs'][dirpath].extend(files)
            continue
    # Third pass: directory-level csv_matrix detection — if >=3 files in a dir
    # are unused .tsv.gz / .csv.gz / .txt.gz files that are NOT archives and NOT
    # already classified elsewhere, treat them as csv_matrix (one-per-sample).
    for dirpath, files in dir_groups.items():
        if not files:
            continue
        # Count files that look like count-matrix CSV/TSV
        candidates = []
        for f in files:
            b = os.path.basename(f).lower()
            ext_match = any(
                b.endswith(ext) for ext in ('.csv', '.csv.gz', '.tsv', '.tsv.gz', '.txt', '.txt.gz')
            )
            if not ext_match:
                continue
            # Exclude known metadata/files already classified
            if detect_10x_mtx_directory([f]) or detect_10x_h5_directory([f]) or detect_fragment_file([f]):
                continue
            if is_archive(f) or is_unsupported(f):
                continue
            candidates.append(f)
        if len(candidates) >= 3:
            result['csv_files'].extend(candidates)
            continue

        # Classify remaining files individually
        for f in files:
            _classify_single_file(f, result, dirpath)

    # Convert defaultdicts to plain dicts
    for key in ('tenx_mtx_dirs', 'tenx_h5_dirs', 'tenx_peak_dirs', 'fragment_dirs'):
        result[key] = dict(result[key])

    return result


def _classify_single_file(filepath: str, result: dict, dirpath: str = '') -> None:
    """Classify a single file that wasn't part of a directory-level pattern."""
    basename = os.path.basename(filepath).lower()

    if basename.endswith('.h5ad'):
        result['h5ad_files'].append(filepath)
    elif basename.endswith('.csv') or basename.endswith('.csv.gz') or basename.endswith('.tsv') or basename.endswith('.tsv.gz') or basename.endswith('.txt') or basename.endswith('.txt.gz'):
        # Special: if it's the ONLY file-like material in the dir + looks
        # like count data, classify as csv_matrix rather than metadata.
        if 'metadata' in basename or 'barcodes' in basename or 'meta' in basename or 'sample_annotation' in basename:
            result['metadata_files'].append(filepath)
        elif is_geo_series_matrix(filepath):
            result['series_matrix_files'].append(filepath)
        elif 'count' in basename or 'tpm' in basename or 'expr' in basename or 'matrix' in basename or 'gene' in basename:
            result['csv_files'].append(filepath)
        elif basename.endswith('.txt.gz') and ('barcode' in basename or 'raw' in basename):
            result['metadata_files'].append(filepath)
        else:
            result['metadata_files'].append(filepath)
    elif is_geo_series_matrix(filepath):
        result['series_matrix_files'].append(filepath)
    elif basename.endswith('.h5') and not basename.endswith('.h5ad'):
        # Check if it's a 10X H5 file with a sample prefix
        if _basename_ends_with_any(basename, TENX_H5_SUFFIXES):
            d = os.path.dirname(filepath) or '.'
            result['tenx_h5_dirs'].setdefault(d, []).append(filepath)
        else:
            result['unmatched'].append(filepath)
    elif basename.endswith('.bed') or basename.endswith('.bed.gz'):
        result['metadata_files'].append(filepath)
    elif basename.endswith('.gz') and not any(basename.endswith(s) for s in ('.csv.gz', '.tsv.gz', '.txt.gz', '.bed.gz', '.mtx.gz', '.tbi.gz', '.h5.gz')):
        # gzip file with unknown content — could be a fragment file backup etc.
        result['unmatched'].append(filepath)
    else:
        result['unmatched'].append(filepath)


def strip_known_suffix(filename: str) -> str:
    """Remove a known GEO/10X suffix from a filename to derive a sample prefix.

    Returns the original filename if no suffix matches.
    """
    for suffix in KNOWN_SUFFIXES:
        if filename.lower().endswith(suffix.lower()):
            return filename[:-len(suffix)]
    # Also try stripping .gz and re-checking
    if filename.lower().endswith('.gz'):
        inner = filename[:-3]
        for suffix in KNOWN_SUFFIXES:
            if inner.lower().endswith(suffix.lower()):
                return inner[:-len(suffix)]
    # Strip common extensions
    for ext in ('.csv', '.tsv', '.txt', '.gz', '.h5', '.h5ad', '.mtx', '.bed'):
        if filename.lower().endswith(ext):
            return filename[:-len(ext)]
    return filename


def group_files_by_sample(file_list: list[str]) -> dict[str, list[str]]:
    """Group files by inferred sample name.

    Strategy:
      1. First, identify files that are themselves data-bearing samples
         (10X H5, fragment files) — each of these becomes its own sample.
      2. Group remaining files by directory boundary.
      3. For files with recognizable suffixes, strip suffix to derive
         sample name.
      4. Merge groups whose sample IDs are identical.

    Returns:
        {sample_id: [file_paths]}
    """
    samples: dict[str, list[str]] = {}

    # Step 1: separate data files that each represent a distinct sample
    # from auxiliary files that belong to samples
    data_files = []
    aux_files = []
    for f in file_list:
        b = os.path.basename(f).lower()
        is_data = False
        for sfx in TENX_H5_SUFFIXES:
            if b.endswith(sfx.lower()):
                is_data = True
                break
        if not is_data:
            for sfx in FRAGMENT_SUFFIXES:
                if b.endswith(sfx.lower()):
                    is_data = True
                    break
        if is_data:
            data_files.append(f)
        else:
            aux_files.append(f)

    # Each data file becomes its own sample (strip the known suffix for the ID)
    for f in data_files:
        sample_id = strip_known_suffix(os.path.basename(f))
        samples.setdefault(sample_id, []).append(f)

    # Step 2: group auxiliary files by parent directory
    dir_groups: dict[str, list[str]] = defaultdict(list)
    for f in aux_files:
        d = os.path.dirname(f) or os.path.basename(f)
        dir_groups[d].append(f)

    # Step 3: within each directory group, try to find common prefix
    for dirname, files in dir_groups.items():
        if len(files) == 1:
            prefix = strip_known_suffix(os.path.basename(files[0]))
            # If prefix exists as an existing data sample, merge there
            existing_sample = _find_matching_sample(prefix, list(samples.keys()))
            if existing_sample:
                samples[existing_sample].extend(files)
            else:
                samples.setdefault(prefix, []).extend(files)
        else:
            basenames_stripped = [strip_known_suffix(os.path.basename(f)) for f in files]
            common = os.path.commonprefix(basenames_stripped)
            common = common.rstrip('_.- ')
            if common and len(common) >= 3:
                samples.setdefault(common, []).extend(files)
            else:
                # Use directory name as group id
                dir_key = os.path.basename(dirname.rstrip('/\\')) or dirname
                samples.setdefault(dir_key, []).extend(files)

    return dict(samples)


def _find_matching_sample(prefix: str, existing_ids: list[str]) -> str | None:
    """Return an existing sample id that matches *prefix*, or None."""
    prefix_lower = prefix.lower()
    for sid in existing_ids:
        if sid.lower() == prefix_lower:
            return sid
        if sid.lower().startswith(prefix_lower) or prefix_lower.startswith(sid.lower()):
            return sid
    return None


def _scan_file_head(filepath: str, n_lines: int = 200) -> str:
    """Read the first *n_lines* of a file (handles .gz transparently).

    Returns '' silently on any error (bad gzip, binary content, etc.).
    """
    import gzip as _gzip
    try:
        if filepath.lower().endswith('.gz'):
            with _gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') as fh:
                lines = []
                for i, line in enumerate(fh):
                    if i >= n_lines:
                        break
                    lines.append(line)
                return '\n'.join(lines)
        else:
            with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
                lines = []
                for i, line in enumerate(fh):
                    if i >= n_lines:
                        break
                    lines.append(line)
                return '\n'.join(lines)
    except Exception:
        return ''


def _sniff_species_from_text(text: str) -> Optional[str]:
    """Scan free text for known species names (common name / Latin binomial)."""
    import re as _re
    lower = text.lower()
    _BINOMIAL_MAP: list[tuple[_re.Pattern, str]] = [
        (_re.compile(r'homo\s+sapiens', _re.IGNORECASE), 'homo_sapiens'),
        (_re.compile(r'mus\s+musculus', _re.IGNORECASE), 'mus_musculus'),
        (_re.compile(r'rattus\s+norvegicus', _re.IGNORECASE), 'rattus_norvegicus'),
        (_re.compile(r'danio\s+rerio', _re.IGNORECASE), 'danio_rerio'),
        (_re.compile(r'bos\s+taurus', _re.IGNORECASE), 'bos_taurus'),
        (_re.compile(r'sus\s+scrofa', _re.IGNORECASE), 'sus_scrofa'),
        (_re.compile(r'macaca\s+mulatta', _re.IGNORECASE), 'macaca_mulatta'),
        (_re.compile(r'callithrix\s+jacchus', _re.IGNORECASE), 'callithrix_jacchus'),
        (_re.compile(r'gallus\s+gallus', _re.IGNORECASE), 'gallus_gallus'),
        (_re.compile(r'drosophila\s+melanogaster', _re.IGNORECASE), 'drosophila_melanogaster'),
        (_re.compile(r'caenorhabditis\s+elegans', _re.IGNORECASE), 'caenorhabditis_elegans'),
        (_re.compile(r'xenopus\s+tropicalis', _re.IGNORECASE), 'xenopus_tropicalis'),
    ]
    for pattern, species in _BINOMIAL_MAP:
        if pattern.search(lower):
            return species
    return None


def guess_species(file_list: list[str],
                  content_scan: bool = True) -> str:
    """Guess species from filenames AND optionally from file contents.

    1.  Filename keyword heuristics (fast, always tried).
    2.  Content scan: read first 200 lines of any CSV/TSV/txt/gz metadata
        or series-matrix file and look for Latin binomial strings.

    Returns a species name string in pipeline-key format (e.g. 'human',
    'mouse', 'zebrafish'), or 'unknown'.

    .. note::

       The return value is now normalised to a **pipeline key** (common name,
       no underscores) that can be used directly in ``CFG.species`` and will
       match the ``_SPECIES_SYNONYMS`` table in ``rna/utils/marker_scoring.py``.
       The older underscore-form keys (e.g. ``homo_sapiens``) are still
       returned as well for backward compatibility, but NEW code should
       prefer the normalised form.
    """
    # Layer 1: filename heuristics
    all_names = ' '.join(os.path.basename(f) for f in file_list).lower()
    for keywords, species_underscored in SPECIES_KEYWORDS:
        for kw in keywords:
            if kw.lower() in all_names:
                return _normalise_species(species_underscored)

    if not content_scan:
        return 'unknown'

    # Layer 2: content scan (metadata / series matrix files)
    for fpath in file_list:
        bn = os.path.basename(fpath).lower()
        if not (bn.endswith(('.csv', '.tsv', '.txt', '.csv.gz', '.tsv.gz', '.txt.gz'))):
            continue
        # Skip pure barcode-only files (very unlikely to carry species info)
        if bn.endswith(('.mtx', '.mtx.gz', '.h5', '.h5ad')):
            continue
        text = _scan_file_head(fpath, n_lines=200)
        if not text:
            continue
        species = _sniff_species_from_text(text)
        if species:
            return _normalise_species(species)

    return 'unknown'


# ── Species name normalisation ──────────────────────────────────────
# Maps both underscored names (homo_sapiens) and Latin binomials
# (Homo sapiens) from guess_species() / _sniff_species_from_text()
# to the canonical pipeline key used by CFG.species and
# _SPECIES_SYNONYMS in rna/utils/marker_scoring.py.
_SPECIES_NORMALISE: Dict[str, str] = {
    # Human
    "homo_sapiens":    "human",
    "Homo sapiens":    "human",
    # Mouse
    "mus_musculus":    "mouse",
    "Mus musculus":    "mouse",
    # Rat
    "rattus_norvegicus":     "rat",
    "Rattus norvegicus":     "rat",
    # Zebrafish
    "danio_rerio":     "zebrafish",
    "Danio rerio":     "zebrafish",
    # Cow
    "bos_taurus":      "cow",
    "Bos taurus":      "cow",
    # Pig
    "sus_scrofa":      "pig",
    "Sus scrofa":      "pig",
    # Macaque (rhesus)
    "macaca_mulatta":   "macaque",
    "Macaca mulatta":   "macaque",
    # Marmoset
    "callithrix_jacchus":   "marmoset",
    "Callithrix jacchus":   "marmoset",
    # Chicken
    "gallus_gallus":   "chicken",
    "Gallus gallus":   "chicken",
    # Fruit fly
    "drosophila_melanogaster":   "drosophila",
    "Drosophila melanogaster":   "drosophila",
    # Worm
    "caenorhabditis_elegans":    "c_elegans",
    "Caenorhabditis elegans":    "c_elegans",
    # Frog
    "xenopus_tropicalis":  "frog",
    "Xenopus tropicalis":  "frog",
    # Macaca fascicularis (cynomolgus) — also maps to macaque
    "Macaca fascicularis":  "macaque",
    # Ictidomys tridecemlineatus
    "Ictidomys tridecemlineatus":  "squirrel",
    # Anolis sagrei
    "Anolis sagrei":  "lizard",
    # Monodelphis domestica
    "Monodelphis domestica":  "opossum",
    # Didelphis marsupialis
    "Didelphis marsupialis":  "opossum",
    # Tupaia chinensis
    "Tupaia chinensis":  "tree_shrew",
    # Tupaia belangeri
    "Tupaia belangeri":  "tree_shrew",
    # Callithrix jacchus (already normalised above; just in case)
    "Callithrix jacchus":  "marmoset",
    # Mustela putorius furo
    "Mustela putorius furo":  "ferret",
    # Ovis aries
    "Ovis aries":  "sheep",
    # Peromyscus maniculatus
    "Peromyscus maniculatus":  "peromyscus",
    # Rhabdomys pumilio
    "Rhabdomys pumilio":  "rhabdomys",
    # Petromyzon marinus
    "Petromyzon marinus":  "lamprey",
}


def _normalise_species(raw: str) -> str:
    """Normalise a species name to its pipeline key.

    Accepts underscored names (``homo_sapiens``) and Latin binomials
    (``Homo sapiens``).  Falls back to *raw* when no normalisation
    entry is found.
    """
    return _SPECIES_NORMALISE.get(raw, _SPECIES_NORMALISE.get(raw.lower(), raw))


def guess_tissue(file_list: list[str]) -> str:
    """Guess tissue from filenames using keyword heuristics.

    Returns a tissue name string or 'unknown'.
    """
    known_tissues = [
        'retina', 'brain', 'liver', 'heart', 'kidney', 'lung',
        'muscle', 'skin', 'intestine', 'pancreas', 'spleen',
        'thymus', 'testis', 'ovary', 'bone_marrow', 'blood',
        'hypothalamus', 'cortex', 'cerebellum', 'hippocampus',
        'rpe', 'optic_nerve', 'spinal_cord',
    ]
    all_names = ' '.join(os.path.basename(f) for f in file_list).lower()
    for tissue in known_tissues:
        if tissue in all_names:
            return tissue
    return 'unknown'


def guess_genome(species: str) -> str:
    """Return a reference genome identifier for a given species name."""
    _map = {
        'homo_sapiens':           'hg38',
        'mus_musculus':           'mm10',
        'rattus_norvegicus':      'rn6',
        'danio_rerio':            'danRer11',
        'bos_taurus':             'bosTau9',
        'sus_scrofa':             'susScr11',
        'macaca_mulatta':         'rheMac10',
        'gallus_gallus':          'galGal6',
        'drosophila_melanogaster': 'dm6',
        'caenorhabditis_elegans': 'ce11',
        'xenopus_tropicalis':     'xenTro10',
    }
    return _map.get(species, '')
