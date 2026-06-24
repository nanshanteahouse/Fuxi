"""
zhang2024.py — Zhang & Peng et al. 2024 (PNAS, PMID: 38598343)
Evolutionary and developmental specialization of foveal cell types in the marmoset.

Species: Callithrix jacchus (common marmoset)
Tissue: retina (neonate + adult, fovea + periphery)
Platform: 10X Genomics
Cells: 44,267 (29,169 foveal + 15,098 peripheral)
Subtypes: 68 (3 PR + 2 HC + 13 BC + 30 AC + 16 RGC + 4 NN)
Groups: 6 (PR, HC, BC, AC, RGC, Non-neuronal)
"""

source_meta = {
    "id": "zhang2024",
    "short_name": "Zhang 2024 PNAS",
    "pmid": "38598343",
    "journal": "PNAS",
    "year": 2024,
    "species": ["Callithrix jacchus"],
    "tissue": "retina",
    "regions": ["fovea", "periphery"],
    "n_cells": 44267,
    "n_subtypes": 68,
    "n_groups": 6,
    "class": "Mammalia",
    "order": "Primates",
}

markers = {
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["38598343"],
            "NRL": ["38598343"],
            "GNAT1": ["38598343"],
        },
        "add": {
            "PDE6B": ["38598343"],
            "SAG": ["38598343"],
            "NR2E3": ["38598343"],
        },
        "refine": {
            "CRX": {
                "note": "pan-photoreceptor marker in marmoset retina",
                "threshold": "log2FC > 1.0",
                "pmid": "38598343",
            },
        },
    },
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["38598343"],
            "GNAT2": ["38598343"],
            "OPN1SW": ["38598343"],
        },
        "add": {
            "PDE6C": ["38598343"],
            "PDE6H": ["38598343"],
            "OPN1LW": ["38598343"],
        },
        "refine": {
            "OPN1LW": {
                "note": "M/L cone opsin; marmoset fovea-enriched",
                "threshold": "log2FC > 2.0",
                "pmid": "38598343",
            },
            "OPN1SW": {
                "note": "S cone opsin; marmoset periphery-enriched",
                "threshold": "log2FC > 1.5",
                "pmid": "38598343",
            },
        },
    },
    "RGC": {
        "confirm": {
            "POU4F2": ["38598343"],
            "NEFM": ["38598343"],
            "POU4F1": ["38598343"],
        },
        "add": {
            "NEFL": ["38598343"],
            "ELAVL4": ["38598343"],
            "GAP43": ["38598343"],
        },
        "refine": {},
    },
    "Bipolar_Cell": {
        "confirm": {
            "VSX1": ["38598343"],
            "PRKCA": ["38598343"],
            "TRPM1": ["38598343"],
        },
        "add": {
            "GRM6": ["38598343"],
            "CABP5": ["38598343"],
            "GRIK1": ["38598343"],
        },
        "refine": {
            "PRKCA": {
                "note": "ON/rod bipolar marker in marmoset",
                "threshold": "log2FC > 1.0",
                "pmid": "38598343",
            },
        },
    },
    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["38598343"],
            "GAD2": ["38598343"],
            "SLC6A9": ["38598343"],
        },
        "add": {
            "TFAP2A": ["38598343"],
            "TFAP2B": ["38598343"],
            "CALB1": ["38598343"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["38598343"],
            "ONECUT2": ["38598343"],
            "PROX1": ["38598343"],
        },
        "add": {
            "CALB1": ["38598343"],
            "LHX1": ["38598343"],
        },
        "refine": {},
    },
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["38598343"],
            "GFAP": ["38598343"],
            "GLUL": ["38598343"],
        },
        "add": {
            "VIM": ["38598343"],
            "SLC1A3": ["38598343"],
            "CRYAB": ["38598343"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["38598343"],
            "P2RY12": ["38598343"],
            "CSF1R": ["38598343"],
        },
        "add": {
            "ITGAM": ["38598343"],
            "CX3CR1": ["38598343"],
            "CD74": ["38598343"],
        },
        "refine": {},
    },
    "RPC": {
        "confirm": {
            "VSX2": ["38598343"],
            "PAX6": ["38598343"],
            "SOX2": ["38598343"],
        },
        "add": {
            "HES1": ["38598343"],
            "NOTCH1": ["38598343"],
            "MKI67": ["38598343"],
        },
        "refine": {
            "MKI67": {
                "note": "proliferation marker; labels neonatal RPCs",
                "threshold": "log2FC > 1.5",
                "pmid": "38598343",
            },
        },
    },
    "RPE": {
        "confirm": {
            "RPE65": ["38598343"],
            "BEST1": ["38598343"],
            "LRAT": ["38598343"],
        },
        "add": {
            "TIMP3": ["38598343"],
            "RDH10": ["38598343"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["38598343"],
            "VWF": ["38598343"],
            "CDH5": ["38598343"],
        },
        "add": {
            "CLDN5": ["38598343"],
            "EGFL7": ["38598343"],
        },
        "refine": {},
    },
    "Pericyte": {
        "confirm": {
            "PDGFRB": ["38598343"],
            "CSPG4": ["38598343"],
            "ACTA2": ["38598343"],
        },
        "add": {
            "RGS5": ["38598343"],
            "ANPEP": ["38598343"],
        },
        "refine": {},
    },
    "Astrocyte": {
        "confirm": {
            "GFAP": ["38598343"],
            "AQP4": ["38598343"],
            "ALDH1L1": ["38598343"],
        },
        "add": {
            "S100B": ["38598343"],
            "SOX9": ["38598343"],
        },
        "refine": {},
    },
    "Oligodendrocyte": {
        "confirm": {
            "MBP": ["38598343"],
            "PLP1": ["38598343"],
            "MOG": ["38598343"],
        },
        "add": {
            "OLIG2": ["38598343"],
            "SOX10": ["38598343"],
        },
        "refine": {},
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["38598343"],
            "DCN": ["38598343"],
            "LUM": ["38598343"],
        },
        "add": {
            "COL3A1": ["38598343"],
            "COL1A2": ["38598343"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "Marmoset_Foveal_Cones",
        "parent": "Cone_Photoreceptor",
        "markers": ["OPN1LW", "OPN1SW", "GNAT2", "ARR3"],
        "species": ["Callithrix jacchus"],
        "source": "zhang2024",
    },
    {
        "name": "Neonatal_RPCs",
        "parent": "RPC",
        "markers": ["MKI67", "TOP2A", "HES1", "SOX2"],
        "species": ["Callithrix jacchus"],
        "source": "zhang2024",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RHO": 1.0, "NRL": 1.0},
        },
        "action": "Rod_Photoreceptor",
        "source": "zhang2024",
        "notes": "RHO+NRL co-expression defines marmoset rod photoreceptors",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ARR3": 1.0, "GNAT2": 1.0},
        },
        "action": "Cone_Photoreceptor",
        "source": "zhang2024",
        "notes": "ARR3+GNAT2 co-expression defines marmoset cone photoreceptors",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RLBP1": 1.0, "GLUL": 1.0},
        },
        "action": "Muller_Glia",
        "source": "zhang2024",
        "notes": "RLBP1+GLUL co-expression defines marmoset Muller glia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "RGC",
        "source": "zhang2024",
        "notes": "POU4F2+NEFM co-expression defines marmoset RGCs",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"VSX2": 1.0, "SOX2": 1.0},
        },
        "action": "RPC",
        "source": "zhang2024",
        "notes": "VSX2+SOX2 co-expression defines marmoset retinal progenitor cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"AIF1": 1.0, "P2RY12": 1.0},
        },
        "action": "Microglia",
        "source": "zhang2024",
        "notes": "AIF1+P2RY12 co-expression defines marmoset microglia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "zhang2024",
        "notes": "ONECUT1+PROX1 co-expression defines marmoset horizontal cells",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"GAD1": 1.0, "GAD2": 1.0},
        },
        "action": "Amacrine_Cell",
        "source": "zhang2024",
        "notes": "GAD1+GAD2 expression defines marmoset GABAergic amacrine cells",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "Cone_Photoreceptor",
            "marker": "OPN1LW",
        },
        "type_b": {
            "cell_type": "Cone_Photoreceptor",
            "marker": "OPN1SW",
        },
        "notes": "Marmoset cones show similar M/L vs S opsin mutual exclusivity as human and macaque; opsin expression patterns are conserved across primates",
        "source": {"a": "zhang2024", "b": "peng2019"},
    },
]
