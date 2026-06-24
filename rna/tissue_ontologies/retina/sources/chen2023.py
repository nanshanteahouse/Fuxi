"""
chen2023.py — Chen et al. 2023 (bioRxiv → published, PMID: 37388908)
A multi-omics atlas of the human retina at single-cell resolution.

Species: Homo sapiens (adult, 25 donors)
Tissue: retina (fovea/macula/periphery)
Platform: 10X Multiome (snRNA-seq + snATAC-seq)
Cells: >250K nuclei (snRNA-seq) + >150K nuclei (snATAC-seq)
Subtypes: >60
"""

source_meta = {
    "id": "chen2023",
    "short_name": "Chen 2023 Multiome",
    "pmid": "37388908",
    "journal": "bioRxiv / Nature Communications",
    "year": 2023,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["fovea", "macula", "periphery"],
    "n_cells": 250000,
    "n_subtypes": 60,
    "n_groups": 6,
    "class": "Mammalia",
    "order": "Primates",
}

markers = {
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["37388908"],
            "NRL": ["37388908"],
            "GNAT1": ["37388908"],
        },
        "add": {
            "PDE6B": ["37388908"],
            "SAG": ["37388908"],
            "NR2E3": ["37388908"],
        },
        "refine": {
            "CRX": {
                "note": "pan-photoreceptor; detected in both snRNA-seq and snATAC-seq",
                "threshold": "log2FC > 1.0",
                "pmid": "37388908",
            },
        },
    },
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["37388908"],
            "GNAT2": ["37388908"],
            "OPN1SW": ["37388908"],
        },
        "add": {
            "PDE6C": ["37388908"],
            "PDE6H": ["37388908"],
            "OPN1LW": ["37388908"],
        },
        "refine": {
            "OPN1LW": {
                "note": "M/L cone opsin; chromatin accessibility confirmed at OPN1LW locus",
                "threshold": "log2FC > 1.5",
                "pmid": "37388908",
            },
        },
    },
    "RGC": {
        "confirm": {
            "POU4F2": ["37388908"],
            "NEFM": ["37388908"],
            "POU4F1": ["37388908"],
        },
        "add": {
            "NEFL": ["37388908"],
            "ELAVL4": ["37388908"],
            "TUBB3": ["37388908"],
        },
        "refine": {},
    },
    "Bipolar_Cell": {
        "confirm": {
            "VSX1": ["37388908"],
            "PRKCA": ["37388908"],
            "TRPM1": ["37388908"],
        },
        "add": {
            "GRM6": ["37388908"],
            "CABP5": ["37388908"],
            "GRIK1": ["37388908"],
        },
        "refine": {},
    },
    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["37388908"],
            "GAD2": ["37388908"],
            "SLC6A9": ["37388908"],
        },
        "add": {
            "TFAP2A": ["37388908"],
            "TFAP2B": ["37388908"],
            "SLC6A1": ["37388908"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["37388908"],
            "ONECUT2": ["37388908"],
            "PROX1": ["37388908"],
        },
        "add": {
            "CALB1": ["37388908"],
            "LHX1": ["37388908"],
        },
        "refine": {},
    },
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["37388908"],
            "GFAP": ["37388908"],
            "GLUL": ["37388908"],
        },
        "add": {
            "VIM": ["37388908"],
            "SLC1A3": ["37388908"],
            "CRYAB": ["37388908"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["37388908"],
            "P2RY12": ["37388908"],
            "CSF1R": ["37388908"],
        },
        "add": {
            "ITGAM": ["37388908"],
            "CX3CR1": ["37388908"],
            "CD74": ["37388908"],
        },
        "refine": {},
    },
    "RPE": {
        "confirm": {
            "RPE65": ["37388908"],
            "BEST1": ["37388908"],
            "LRAT": ["37388908"],
        },
        "add": {
            "CRALBP": ["37388908"],
            "TIMP3": ["37388908"],
            "RDH10": ["37388908"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["37388908"],
            "VWF": ["37388908"],
            "CDH5": ["37388908"],
        },
        "add": {
            "CLDN5": ["37388908"],
            "EGFL7": ["37388908"],
            "FLT1": ["37388908"],
        },
        "refine": {},
    },
    "Pericyte": {
        "confirm": {
            "PDGFRB": ["37388908"],
            "CSPG4": ["37388908"],
            "ACTA2": ["37388908"],
        },
        "add": {
            "RGS5": ["37388908"],
            "ANPEP": ["37388908"],
        },
        "refine": {},
    },
    "RPC": {
        "confirm": {
            "VSX2": ["37388908"],
            "PAX6": ["37388908"],
            "SOX2": ["37388908"],
        },
        "add": {
            "HES1": ["37388908"],
            "NOTCH1": ["37388908"],
        },
        "refine": {},
    },
    "Astrocyte": {
        "confirm": {
            "GFAP": ["37388908"],
            "AQP4": ["37388908"],
            "ALDH1L1": ["37388908"],
        },
        "add": {
            "S100B": ["37388908"],
            "SOX9": ["37388908"],
        },
        "refine": {},
    },
    "Oligodendrocyte": {
        "confirm": {
            "MBP": ["37388908"],
            "PLP1": ["37388908"],
            "MOG": ["37388908"],
        },
        "add": {
            "OLIG2": ["37388908"],
            "SOX10": ["37388908"],
            "MOBP": ["37388908"],
        },
        "refine": {},
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["37388908"],
            "DCN": ["37388908"],
            "LUM": ["37388908"],
        },
        "add": {
            "COL3A1": ["37388908"],
            "COL1A2": ["37388908"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "Macula_Specific_Cone",
        "parent": "Cone_Photoreceptor",
        "markers": ["OPN1LW", "ARR3", "GNAT2", "PDE6H"],
        "species": ["Homo sapiens"],
        "source": "chen2023",
    },
    {
        "name": "Peripheral_Amacrine",
        "parent": "Amacrine_Cell",
        "markers": ["SLC6A1", "GAD1", "TFAP2A"],
        "species": ["Homo sapiens"],
        "source": "chen2023",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RHO": 1.0, "NRL": 1.0},
        },
        "action": "Rod_Photoreceptor",
        "source": "chen2023",
        "notes": "RHO+NRL co-expression defines rod photoreceptors in multi-omics atlas",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ARR3": 1.0, "GNAT2": 1.0},
        },
        "action": "Cone_Photoreceptor",
        "source": "chen2023",
        "notes": "ARR3+GNAT2 co-expression defines cone photoreceptors in multi-omics atlas",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RLBP1": 1.0, "GLUL": 1.0},
        },
        "action": "Muller_Glia",
        "source": "chen2023",
        "notes": "RLBP1+GLUL co-expression defines Muller glia in multi-omics atlas",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "RGC",
        "source": "chen2023",
        "notes": "POU4F2+NEFM co-expression defines RGCs in multi-omics atlas",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RPE65": 1.0, "BEST1": 1.0},
        },
        "action": "RPE",
        "source": "chen2023",
        "notes": "RPE65+BEST1 co-expression defines RPE cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"AIF1": 1.0, "P2RY12": 1.0},
        },
        "action": "Microglia",
        "source": "chen2023",
        "notes": "AIF1+P2RY12 co-expression defines retinal microglia in multi-omics atlas",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "chen2023",
        "notes": "ONECUT1+PROX1 co-expression defines horizontal cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"MBP": 1.0, "PLP1": 1.0},
        },
        "action": "Oligodendrocyte",
        "source": "chen2023",
        "notes": "MBP+PLP1 co-expression defines oligodendrocytes in the retina",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "Oligodendrocyte",
            "marker": "MBP",
        },
        "type_b": {
            "cell_type": "Glial_Cell",
            "marker": "MBP",
        },
        "notes": "MBP is a specific marker for oligodendrocytes; some sources use the more generic 'glial cell' label for all MBP+ cells",
        "source": {"a": "chen2023", "b": "hu2019"},
    },
]
