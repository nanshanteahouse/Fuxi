"""Shared fragment-file chromosome-size auto-detection for ATAC steps.

SnapATAC2's ``import_fragments`` requires exact chromosome bounds.  The
observed max fragment end per chromosome (+ a small buffer) is species
agnostic — hg38 sizes are only used as an early-exit completion signal,
never as the returned size, so mm10 / other genomes get correct bounds.
"""

from __future__ import annotations

import gzip

# hg38 autosome + XY + M sizes (early-exit signal only).
HG38_CHROM_SIZES = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}
_N_STANDARD_CHROMS = len(HG38_CHROM_SIZES)


def auto_chrom_sizes(fragment_file: str) -> dict:
    """Auto-detect chromosome sizes from a fragments.tsv.gz file.

    Tracks the *observed* max end per chromosome (plus a 10kb buffer) and
    stops early once all hg38-standard chromosomes have been seen.
    """
    chrom_max = {}
    chroms_found = set()
    with gzip.open(fragment_file, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            c, end = parts[0], int(parts[2])
            if c not in chrom_max or end > chrom_max[c]:
                chrom_max[c] = end + 10000
            if c in HG38_CHROM_SIZES:
                chroms_found.add(c)
                if len(chroms_found) >= _N_STANDARD_CHROMS:
                    break
    return chrom_max
