#!/usr/bin/env python3
"""
rna/ortholog.py — Cross-species ortholog gene mapping for scRNA-seq pipelines
==================================================================================

Converts species-specific Ensembl gene IDs (e.g. ENSSSCG..., ENSMPUG...) to human
gene symbols via Ensembl 1:1 ortholog mapping.  This enables:

1. LLM annotation with recognisable gene names (human gene symbols)
2. Marker-dict cross-validation (KB uses human gene names)
3. Cross-species KB scoring (marker overlap checks)

Ortholog mappings are pre-downloaded via ``data/orthologs/fetch_all_orthologs.py``
and cached as JSON files in ``data/orthologs/``.

Usage::

    from rna.ortholog import OrthologMapper, convert_species_gene_names

    adata = convert_species_gene_names(adata, species="pig")

Reference
---------
- Hahn et al. (2023) Nature 624:415-424 — cross-species retina atlases
- Ensembl BioMart REST API: https://www.ensembl.org/biomart/
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# ── Ensembl gene ID pattern ────────────────────────────────────────────
# ENS + 0-4 species prefix chars + G + 11 digits
_ENSEMBL_PATTERN = re.compile(r'^ENS[A-Z]{0,4}G\d{11}$')

# ── Species name → BioMart dataset (verified 2026-06-23) ───────────────
SPECIES_TO_DATASET: Dict[str, str] = {
    "cow":          "btaurus_gene_ensembl",
    "ferret":       "mpfuro_gene_ensembl",
    "lizard":       "acarolinensis_gene_ensembl",
    "marmoset":     "cjacchus_gene_ensembl",
    "opossum":      "mdomestica_gene_ensembl",
    "peromyscus":   "pmbairdii_gene_ensembl",
    "pig":          "sscrofa_gene_ensembl",
    "rhabdomys":    None,  # not in Ensembl BioMart
    "sheep":        "oaries_gene_ensembl",
    "squirrel":     "itridecemlineatus_gene_ensembl",
    "tree_shrew":   "tbelangeri_gene_ensembl",
    "zebrafish":    "drerio_gene_ensembl",
    "human":        "hsapiens_gene_ensembl",
    "mouse":        "mmusculus_gene_ensembl",
    "macaque":      "mfascicularis_gene_ensembl",
}

# ── Species name → Ensembl REST API species name ───────────────────────
SPECIES_TO_REST_NAME: Dict[str, str] = {
    "cow":          "bos_taurus",
    "ferret":       "mustela_putorius_furo",
    "lizard":       "anolis_sagrei",
    "marmoset":     "callithrix_jacchus",
    "opossum":      "monodelphis_domestica",
    "peromyscus":   "peromyscus_maniculatus",
    "pig":          "sus_scrofa",
    "rhabdomys":    "rhabdomys_pumilio",
    "sheep":        "ovis_aries",
    "squirrel":     "ictidomys_tridecemlineatus",
    "tree_shrew":   "tupaia_chinensis",
    "zebrafish":    "danio_rerio",
    "human":        "homo_sapiens",
    "mouse":        "mus_musculus",
    "macaque":      "macaca_fascicularis",
}

# ── Species name → taxonomic class (纲) ──────────────────────────────
# Used by phylogenetic_weight() to penalise/reward cross-class matches.
# The KB stores class annotations per cell type (e.g. "Mammalia");
# score_cluster_against_kb() compares target_class against source classes
# to apply taxonomic-distance weighting.
SPECIES_TO_CLASS: Dict[str, str] = {
    "human":        "Mammalia",
    "mouse":        "Mammalia",
    "macaque":      "Mammalia",
    "marmoset":     "Mammalia",
    "tree_shrew":   "Mammalia",
    "cow":          "Mammalia",
    "pig":          "Mammalia",
    "sheep":        "Mammalia",
    "ferret":       "Mammalia",
    "squirrel":     "Mammalia",
    "opossum":      "Mammalia",
    "peromyscus":   "Mammalia",
    "rhabdomys":    "Mammalia",
    "lizard":       "Reptilia",
    "chicken":      "Aves",
    "frog":         "Amphibia",
    "zebrafish":    "Teleostei",
    "lamprey":      "Petromyzontida",
}


# ═══════════════════════════════════════════════════════════════════════
#  Gene ID type detection
# ═══════════════════════════════════════════════════════════════════════

def detect_gene_id_type(var_names) -> str:
    """Detect whether gene names are Ensembl IDs, standard symbols, or mixed.

    Parameters
    ----------
    var_names : array-like of str
        Gene identifiers (e.g. ``adata.var_names``).

    Returns
    -------
    str
        ``"ensembl_species"`` — >70% Ensembl IDs.
        ``"mixed"``          — some Ensembl IDs, some symbols.
        ``"gene_symbol"``    — no or very few Ensembl IDs.
    """
    sample = list(var_names[:200]) if hasattr(var_names, '__iter__') else []
    if not sample:
        return "gene_symbol"

    n_ensembl = sum(1 for g in sample if _ENSEMBL_PATTERN.match(str(g)))
    n_total = len(sample)
    if n_total == 0:
        return "gene_symbol"

    ratio = n_ensembl / n_total
    if ratio > 0.7:
        return "ensembl_species"
    elif ratio > 0.0:
        return "mixed"
    return "gene_symbol"


# ═══════════════════════════════════════════════════════════════════════
#  Orphan gene patterns (genes that are definitely not standard symbols)
# ═══════════════════════════════════════════════════════════════════════

# Genes matching these patterns won't be treated as recognisable symbols
_ORPHAN_PATTERNS = [
    re.compile(r'^FUN-\d+$'),        # Rhabdomys unannotated
    re.compile(r'^LOC\d+$'),         # NCBI uncharacterised loci
    re.compile(r'^LINC-.*$'),        # Long intergenic non-coding RNAs (provisional)
    re.compile(r'^LORF\d+.*$'),      # Large open reading frames
    re.compile(r'^RF\d+$'),          # RNA family
    re.compile(r'^C\d+orf\d+$', re.IGNORECASE),  # Chromosome ORF
    re.compile(r'^AC\d+\.\d+$'),     # Clone-based gene models
    re.compile(r'^AP\d+\.\d+$'),     # Clone-based gene models
    re.compile(r'^RP\d+[-.]'),       # Ribosomal protein (pseudogene-like)
]

def _is_orphan_name(gene: str) -> bool:
    """Return True if *gene* looks like a non-standard provisional ID."""
    for pat in _ORPHAN_PATTERNS:
        if pat.match(gene):
            return True
    return False


# ═══════════════════════════════════════════════════════════════════════
#  OrthologMapper
# ═══════════════════════════════════════════════════════════════════════

class OrthologMapper:
    """Apply cached 1:1 ortholog mappings to convert gene names across species.

    Parameters
    ----------
    ortholog_cache_dir : str or Path
        Directory for cached ortholog JSON files.
    """

    def __init__(self, ortholog_cache_dir: str = "data/orthologs"):
        self.cache_dir = Path(ortholog_cache_dir)
        self._mapping_cache: Dict[str, dict] = {}

    # ── Public API ────────────────────────────────────────────────────

    def get_1to1_orthologs(self, species: str, target: str = "human") -> dict:
        """Return ``{species_gene_id → human_gene_symbol}`` for 1:1 orthologs."""
        cache_key = f"{species}_to_{target}"
        if cache_key in self._mapping_cache:
            return self._mapping_cache[cache_key]

        cache_file = self.cache_dir / f"{cache_key}_orthologs.json"
        if cache_file.exists():
            mapping = json.loads(cache_file.read_text())
            self._mapping_cache[cache_key] = mapping
            logger.info("Loaded %d orthologs from cache: %s", len(mapping), cache_file)
            return mapping

        # No cache → species not supported; return empty
        logger.warning("No ortholog cache for %s→%s (file: %s)", species, target, cache_file)
        return {}

    def convert_var_names(self, adata, species: str, target: str = "human"):
        """Convert ``adata.var_names`` from species identifiers to human symbols.

        Modifies *adata* in-place:

        - Sets ``adata.var['original_gene']`` to the original names.
        - Replaces ``var_names`` with human orthologs where possible.
        - Prefixed ``UNMAPPED_`` for genes without a 1:1 ortholog.
        - Keeps genes that already match standard human symbols as-is.

        Parameters
        ----------
        adata : AnnData
        species : str
            Source species common name (e.g. ``"pig"``).
        target : str
            Target species (default ``"human"``).

        Returns
        -------
        AnnData
            Same object (in-place modification).
        """
        gene_type = detect_gene_id_type(adata.var_names)
        logger.info("Gene ID type for %s: %s", species, gene_type)

        original_names = adata.var_names.astype(str).tolist()
        adata.var['original_gene'] = original_names

        # Strategy depends on gene ID type:
        # - ensembl_species / mixed: use ortholog mapping
        # - gene_symbol: most genes already use standard symbols → keep as-is
        if gene_type == "gene_symbol":
            # Only replace truly orphan identifiers (FUN-*, LOC*, etc.)
            new_names = [
                f"UNMAPPED_{g}" if _is_orphan_name(g) else g
                for g in original_names
            ]
            n_mapped = len(original_names) - sum(1 for g in original_names if _is_orphan_name(g))
            n_ensembl_mapped = 0
        else:
            mapping = self.get_1to1_orthologs(species, target)
            new_names = []
            n_mapped = 0
            n_ensembl_mapped = 0
            for g in original_names:
                if g in mapping:
                    new_names.append(mapping[g])
                    n_ensembl_mapped += 1
                elif _ENSEMBL_PATTERN.match(str(g)):
                    # Ensembl ID without ortholog → keep as unmapped
                    new_names.append(f"UNMAPPED_{g}")
                else:
                    # Already a gene symbol (mixed dataset)
                    new_names.append(g)
                    n_mapped += 1
            n_mapped += n_ensembl_mapped

        adata.var_names = new_names
        pct = n_mapped / max(len(original_names), 1) * 100
        logger.info(
            "Ortholog conversion (%s→%s): %d/%d genes mapped (%.1f%%, "
            "ensembl_mapped=%d)",
            species, target, n_mapped, len(original_names), pct, n_ensembl_mapped,
        )
        return adata


# ═══════════════════════════════════════════════════════════════════════
#  Convenience function
# ═══════════════════════════════════════════════════════════════════════

def convert_species_gene_names(adata, species: str, target: str = "human",
                               cache_dir: str = "data/orthologs"):
    """Convert gene names in *adata* from *species* identifiers to *target*.

    Convenience wrapper around :class:`OrthologMapper`.

    Parameters
    ----------
    adata : AnnData
    species : str
        Source species common name.
    target : str
        Target species (default ``"human"``).
    cache_dir : str
        Directory containing the pre-downloaded ortholog JSON files.

    Returns
    -------
    AnnData
        Same object, modified in-place.
    """
    mapper = OrthologMapper(cache_dir)
    return mapper.convert_var_names(adata, species, target)
