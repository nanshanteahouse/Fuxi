"""
peng2020.py — Peng et al. 2020 (Scientific Reports, PMID: 32555229)
Cell Atlas of the Human Fovea and Peripheral Retina.

Species: Homo sapiens (adult)
Tissue: retina (fovea + periphery, 7 donors)
Platform: 10X Genomics droplet sequencing
Cells: 84,982 (fovea 55,736 + periphery 29,246)
Subtypes: 58 (3 PR + 2 HC + 12 BC + 15 AC + 12 RGC + 5 non-neuronal)
Groups: 6 (PR, HC, BC, AC, RGC, Non-neuronal)
"""

source_meta = {
    "id": "peng2020",
    "short_name": "Peng 2020 Sci Rep",
    "pmid": "32555229",
    "journal": "Scientific Reports",
    "year": 2020,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["fovea", "periphery"],
    "n_cells": 84982,
    "n_subtypes": 58,
    "n_groups": 6,
    "class": "Mammalia",
    "order": "Primates",
}

markers = {
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["32555229"],
            "NRL": ["32555229"],
            "GNAT1": ["32555229"],
        },
        "add": {
            "PDE6B": ["32555229"],
            "SAG": ["32555229"],
            "NR2E3": ["32555229"],
        },
        "refine": {
            "CRX": {
                "note": "pan-photoreceptor marker; also in cone precursors",
                "threshold": "log2FC > 1.0",
                "pmid": "32555229",
            },
        },
    },
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["32555229"],
            "GNAT2": ["32555229"],
            "OPN1SW": ["32555229"],
        },
        "add": {
            "PDE6C": ["32555229"],
            "PDE6H": ["32555229"],
            "OPN1LW": ["32555229"],
        },
        "refine": {
            "OPN1LW": {
                "note": "M/L cone opsin; human fovea-enriched",
                "threshold": "log2FC > 2.0",
                "pmid": "32555229",
            },
            "OPN1SW": {
                "note": "S cone opsin; enriched in peripheral retina",
                "threshold": "log2FC > 1.5",
                "pmid": "32555229",
            },
        },
    },
    "RGC": {
        "confirm": {
            "POU4F2": ["32555229"],
            "NEFM": ["32555229"],
            "POU4F1": ["32555229"],
        },
        "add": {
            "NEFL": ["32555229"],
            "ELAVL4": ["32555229"],
            "TUBB3": ["32555229"],
        },
        "refine": {},
    },
    "Bipolar_Cell": {
        "confirm": {
            "VSX1": ["32555229"],
            "PRKCA": ["32555229"],
            "TRPM1": ["32555229"],
        },
        "add": {
            "GRM6": ["32555229"],
            "CABP5": ["32555229"],
            "GRIK1": ["32555229"],
        },
        "refine": {
            "PRKCA": {
                "note": "ON/rod bipolar cell marker in human retina",
                "threshold": "log2FC > 1.0",
                "pmid": "32555229",
            },
        },
    },
    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["32555229"],
            "GAD2": ["32555229"],
            "SLC6A9": ["32555229"],
        },
        "add": {
            "TFAP2A": ["32555229"],
            "TFAP2B": ["32555229"],
            "SLC6A1": ["32555229"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["32555229"],
            "ONECUT2": ["32555229"],
            "PROX1": ["32555229"],
        },
        "add": {
            "CALB1": ["32555229"],
            "LHX1": ["32555229"],
        },
        "refine": {},
    },
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["32555229"],
            "GFAP": ["32555229"],
            "GLUL": ["32555229"],
        },
        "add": {
            "VIM": ["32555229"],
            "SLC1A3": ["32555229"],
            "CRYAB": ["32555229"],
        },
        "refine": {},
    },
    "RPE": {
        "confirm": {
            "RPE65": ["32555229"],
            "BEST1": ["32555229"],
            "LRAT": ["32555229"],
        },
        "add": {
            "CRALBP": ["32555229"],
            "TIMP3": ["32555229"],
            "RDH10": ["32555229"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["32555229"],
            "P2RY12": ["32555229"],
            "CSF1R": ["32555229"],
        },
        "add": {
            "ITGAM": ["32555229"],
            "CX3CR1": ["32555229"],
            "CD74": ["32555229"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["32555229"],
            "VWF": ["32555229"],
            "CDH5": ["32555229"],
        },
        "add": {
            "CLDN5": ["32555229"],
            "EGFL7": ["32555229"],
        },
        "refine": {},
    },
    "Pericyte": {
        "confirm": {
            "PDGFRB": ["32555229"],
            "CSPG4": ["32555229"],
            "ACTA2": ["32555229"],
        },
        "add": {
            "RGS5": ["32555229"],
            "ANPEP": ["32555229"],
        },
        "refine": {},
    },
    "RPC": {
        "confirm": {
            "VSX2": ["32555229"],
            "PAX6": ["32555229"],
            "SOX2": ["32555229"],
        },
        "add": {
            "HES1": ["32555229"],
            "NOTCH1": ["32555229"],
        },
        "refine": {},
    },
    "Astrocyte": {
        "confirm": {
            "GFAP": ["32555229"],
            "AQP4": ["32555229"],
            "ALDH1L1": ["32555229"],
        },
        "add": {
            "S100B": ["32555229"],
            "SOX9": ["32555229"],
        },
        "refine": {},
    },
    "Oligodendrocyte": {
        "confirm": {
            "MBP": ["32555229"],
            "PLP1": ["32555229"],
            "MOG": ["32555229"],
        },
        "add": {
            "OLIG2": ["32555229"],
            "SOX10": ["32555229"],
        },
        "refine": {},
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["32555229"],
            "DCN": ["32555229"],
            "LUM": ["32555229"],
        },
        "add": {
            "COL3A1": ["32555229"],
            "COL1A2": ["32555229"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "Foveal_S_Cones",
        "parent": "Cone_Photoreceptor",
        "markers": ["OPN1SW", "ARR3", "GNAT2"],
        "species": ["Homo sapiens"],
        "source": "peng2020",
    },
    {
        "name": "Foveal_ML_Cones",
        "parent": "Cone_Photoreceptor",
        "markers": ["OPN1LW", "ARR3", "GNAT2"],
        "species": ["Homo sapiens"],
        "source": "peng2020",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RHO": 1.0, "NRL": 1.0},
        },
        "action": "Rod_Photoreceptor",
        "source": "peng2020",
        "notes": "RHO+NRL co-expression defines rod photoreceptors in human",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ARR3": 1.0, "GNAT2": 1.0},
        },
        "action": "Cone_Photoreceptor",
        "source": "peng2020",
        "notes": "ARR3+GNAT2 co-expression defines cone photoreceptors in human",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RLBP1": 1.0, "GLUL": 1.0},
        },
        "action": "Muller_Glia",
        "source": "peng2020",
        "notes": "RLBP1+GLUL co-expression defines human Muller glia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "RGC",
        "source": "peng2020",
        "notes": "POU4F2+NEFM co-expression defines human RGCs",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RPE65": 1.0, "BEST1": 1.0},
        },
        "action": "RPE",
        "source": "peng2020",
        "notes": "RPE65+BEST1 co-expression defines retinal pigment epithelium",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"AIF1": 1.0, "P2RY12": 1.0},
        },
        "action": "Microglia",
        "source": "peng2020",
        "notes": "AIF1+P2RY12 co-expression defines human retinal microglia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "peng2020",
        "notes": "ONECUT1+PROX1 co-expression defines human horizontal cells",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"VSX1": 1.0, "TRPM1": 1.0},
        },
        "action": "Bipolar_Cell",
        "source": "peng2020",
        "notes": "VSX1+TRPM1 co-expression defines human bipolar cells",
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
        "notes": "M/L cones and S cones are both cone subtypes distinguished by opsin expression; OPN1LW and OPN1SW are mutually exclusive in mature cones",
        "source": {"a": "peng2020", "b": "peng2019"},
    },
]
