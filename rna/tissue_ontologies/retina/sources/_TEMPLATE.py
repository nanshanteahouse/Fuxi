"""
_TEMPLATE.py — Schema template for retina cell type source data files.

This file documents the required structure for each source file in this directory.
Each source file defines pure Python dict/list constants only — NO business logic,
NO executable code beyond dict/list definitions.

=== SCHEMA ===

1. source_meta (dict):
    Metadata describing the publication/dataset this source file represents.
    Fields:
        id (str):           Unique source identifier, lowercase, e.g. "peng2019"
        short_name (str):   Readable label, e.g. "Peng 2019 Cell"
        pmid (str):         PubMed ID of the primary publication
        journal (str):      Journal abbreviation
        year (int):         Publication year
        species (list[str]): List of species covered, using scientific names
        tissue (str):       "retina"
        regions (list[str]): Anatomical regions, e.g. ["fovea", "periphery"]
        n_cells (int):      Approximate number of cells profiled
        n_subtypes (int):   Number of reported cell subtypes
        n_groups (int):     Number of major cell groups

2. markers (dict):
    Cell-type-level marker gene annotations reported by this source.
    Key: cell type name (e.g. "Rod_Photoreceptor")
    Value: dict with keys:
        confirm (dict): Markers this paper confirms from other sources.
            Key: gene symbol (str, uppercase HGNC)
            Value: list[str] of PMIDs supporting this marker assignment.
        add (dict): Markers this paper reports that other papers may not.
            Key: gene symbol (str, uppercase HGNC)
            Value: list[str] of PMIDs.
        refine (dict): Specialized annotations for individual markers.
            Key: gene symbol (str, uppercase HGNC)
            Value: dict with optional fields:
                note (str): Free-text annotation
                threshold (str): Expression threshold string
                pmid (str): PMID for this refinement.

3. novel_types (list[dict]):
    Cell types this paper reports that are NOT found in other sources.
    Each entry:
        name (str):         Cell type name
        parent (str):       Parent ontology type, e.g. "Amacrine_Cell"
        markers (list[str]): Marker gene symbols
        species (list[str]): Species in which this type was found
        source (str):       Source ID (same as source_meta.id)

4. expert_rules (list[dict]):
    Expert-defined decision rules for cell type annotation.
    Each entry:
        priority (int):     Rule priority (lower = higher priority)
        condition (dict):   Trigger condition with markers_present sub-dict:
            markers_present (dict): Key=gene symbol, value=float threshold
        action (str):       Cell type name to assign when condition is met
        source (str):       Source ID for this rule
        notes (str):        Explanation of the rule

5. conflicts (list[dict]):
    Known cross-source annotation conflicts.
    Each entry:
        type_a (dict):      First conflicting type with cell_type and marker keys
        type_b (dict):      Second conflicting type with cell_type and marker keys
        notes (str):        Description of the conflict
        source (dict):      Source IDs: a and b

=== REFERENCE EXPERT RULES (minimum 8 across all sources) ===
Rod:         RHO + NRL
Cone:        ARR3 + GNAT2
Muller_Glia: RLBP1 + GLUL
RGC:         POU4F2 + NEFM
RPC:         VSX2 + SOX2
RPE:         RPE65 + BEST1
Microglia:   AIF1 + P2RY12
Horizontal_Cell: ONECUT1 + PROX1

=== CORE CELL TYPES (each source must cover >=3 markers per type) ===
Rod_Photoreceptor, Cone_Photoreceptor, RGC, Bipolar_Cell,
Amacrine_Cell, Horizontal_Cell, Muller_Glia, RPC, RPE, Microglia,
Vascular_Endothelial, Pericyte, Astrocyte, Oligodendrocyte, Fibroblast
"""

# === METADATA ===
# source_meta = {
#     "id": "template",
#     "short_name": "Template Source",
#     "pmid": "00000000",
#     "journal": "Template Journal",
#     "year": 2024,
#     "species": ["Homo sapiens"],
#     "tissue": "retina",
#     "regions": ["fovea", "periphery"],
#     "n_cells": 0,
#     "n_subtypes": 0,
#     "n_groups": 0,
# }

# === MARKERS ===
# markers = {
#     "Rod_Photoreceptor": {
#         "confirm": {
#             "RHO": ["00000000"],
#             "NRL": ["00000000"],
#             "GNAT1": ["00000000"],
#         },
#         "add": {
#             "SAG": ["00000000"],
#         },
#         "refine": {
#             "CRX": {
#                 "note": "pan-photoreceptor marker",
#                 "threshold": "log2FC > 1.0",
#                 "pmid": "00000000",
#             },
#         },
#     },
# }

# === NOVEL TYPES ===
# novel_types = [
#     # {
#     #     "name": "Subtype_Name",
#     #     "parent": "Parent_Type",
#     #     "markers": ["GENE1", "GENE2"],
#     #     "species": ["Homo sapiens"],
#     #     "source": "source_id",
#     # },
# ]

# === EXPERT RULES ===
# expert_rules = [
#     # {
#     #     "priority": 1,
#     #     "condition": {
#     #         "markers_present": {
#     #             "RHO": 1.0,
#     #             "NRL": 1.0,
#     #         },
#     #     },
#     #     "action": "Rod_Photoreceptor",
#     #     "source": "source_id",
#     #     "notes": "RHO+NRL co-expression specifies rod photoreceptor fate",
#     # },
# ]

# === CONFLICTS ===
# conflicts = [
#     # {
#     #     "type_a": {
#     #         "cell_type": "...",
#     #         "marker": "...",
#     #     },
#     #     "type_b": {
#     #         "cell_type": "...",
#     #         "marker": "...",
#     #     },
#     #     "notes": "...",
#     #     "source": {"a": "...", "b": "..."},
#     # },
# ]
