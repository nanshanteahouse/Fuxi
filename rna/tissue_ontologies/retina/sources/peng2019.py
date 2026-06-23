"""
peng2019.py — Peng et al. 2019 (Cell, PMID: 30712875)
Molecular specification of cell types underlying central and peripheral vision in primates.

Species: Macaca fascicularis (cynomolgus monkey)
Tissue: retina (fovea + periphery)
Platform: 10X Genomics
Cells: 165,679 (fovea 92,627 + periphery 73,052)
Subtypes: >60
Groups: 6 (Photoreceptors, Horizontal, Bipolar, Amacrine, RGC, Non-neuronal)
"""

source_meta = {
    "id": "peng2019",
    "short_name": "Peng 2019 Cell",
    "pmid": "30712875",
    "journal": "Cell",
    "year": 2019,
    "species": ["Macaca fascicularis"],
    "tissue": "retina",
    "regions": ["fovea", "periphery"],
    "n_cells": 165679,
    "n_subtypes": 60,
    "n_groups": 6,
}

markers = {
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["30712875"],
            "NRL": ["30712875"],
            "GNAT1": ["30712875"],
        },
        "add": {
            "PDE6B": ["30712875"],
            "SAG": ["30712875"],
            "NR2E3": ["30712875"],
        },
        "refine": {
            "CRX": {
                "note": "pan-photoreceptor; expressed in both rods and cones",
                "threshold": "log2FC > 1.5",
                "pmid": "30712875",
            },
        },
    },
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["30712875"],
            "GNAT2": ["30712875"],
            "OPN1SW": ["30712875"],
        },
        "add": {
            "PDE6C": ["30712875"],
            "PDE6H": ["30712875"],
            "OPN1LW": ["30712875"],
        },
        "refine": {
            "OPN1LW": {
                "note": "M/L cone opsin; fovea-enriched expression",
                "threshold": "log2FC > 2.0",
                "pmid": "30712875",
            },
            "OPN1SW": {
                "note": "S cone opsin; periphery-enriched relative to fovea",
                "threshold": "log2FC > 1.5",
                "pmid": "30712875",
            },
        },
    },
    "RGC": {
        "confirm": {
            "POU4F2": ["30712875"],
            "NEFM": ["30712875"],
            "POU4F1": ["30712875"],
        },
        "add": {
            "NEFL": ["30712875"],
            "ELAVL4": ["30712875"],
            "TUBB3": ["30712875"],
        },
        "refine": {},
    },
    "Bipolar_Cell": {
        "confirm": {
            "VSX1": ["30712875"],
            "PRKCA": ["30712875"],
            "TRPM1": ["30712875"],
        },
        "add": {
            "GRM6": ["30712875"],
            "CABP5": ["30712875"],
            "VSX2": ["30712875"],
        },
        "refine": {
            "PRKCA": {
                "note": "ON bipolar cell marker; labels rod bipolar cells",
                "threshold": "log2FC > 1.0",
                "pmid": "30712875",
            },
        },
    },
    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["30712875"],
            "GAD2": ["30712875"],
            "SLC6A9": ["30712875"],
        },
        "add": {
            "TFAP2A": ["30712875"],
            "TFAP2B": ["30712875"],
            "CALB1": ["30712875"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["30712875"],
            "ONECUT2": ["30712875"],
            "PROX1": ["30712875"],
        },
        "add": {
            "CALB1": ["30712875"],
            "LHX1": ["30712875"],
        },
        "refine": {},
    },
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["30712875"],
            "GFAP": ["30712875"],
            "GLUL": ["30712875"],
        },
        "add": {
            "VIM": ["30712875"],
            "SLC1A3": ["30712875"],
            "CRYAB": ["30712875"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["30712875"],
            "P2RY12": ["30712875"],
            "CSF1R": ["30712875"],
        },
        "add": {
            "ITGAM": ["30712875"],
            "CX3CR1": ["30712875"],
            "CD74": ["30712875"],
        },
        "refine": {},
    },
    "RPE": {
        "confirm": {
            "RPE65": ["30712875"],
            "BEST1": ["30712875"],
            "LRAT": ["30712875"],
        },
        "add": {
            "TIMP3": ["30712875"],
            "RDH10": ["30712875"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["30712875"],
            "VWF": ["30712875"],
            "CDH5": ["30712875"],
        },
        "add": {
            "CLDN5": ["30712875"],
            "EGFL7": ["30712875"],
        },
        "refine": {},
    },
    "Pericyte": {
        "confirm": {
            "PDGFRB": ["30712875"],
            "CSPG4": ["30712875"],
            "ACTA2": ["30712875"],
        },
        "add": {
            "RGS5": ["30712875"],
            "ANPEP": ["30712875"],
        },
        "refine": {},
    },
    "RPC": {
        "confirm": {
            "VSX2": ["30712875"],
            "SOX2": ["30712875"],
            "PAX6": ["30712875"],
        },
        "add": {
            "HES1": ["30712875"],
            "NOTCH1": ["30712875"],
        },
        "refine": {},
    },
    "Astrocyte": {
        "confirm": {
            "GFAP": ["30712875"],
            "AQP4": ["30712875"],
            "ALDH1L1": ["30712875"],
        },
        "add": {
            "S100B": ["30712875"],
            "SOX9": ["30712875"],
        },
        "refine": {},
    },
    "Oligodendrocyte": {
        "confirm": {
            "MBP": ["30712875"],
            "PLP1": ["30712875"],
            "MOG": ["30712875"],
        },
        "add": {
            "OLIG2": ["30712875"],
            "SOX10": ["30712875"],
        },
        "refine": {},
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["30712875"],
            "DCN": ["30712875"],
            "LUM": ["30712875"],
        },
        "add": {
            "COL3A1": ["30712875"],
            "COL1A2": ["30712875"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "Foveal_Cones",
        "parent": "Cone_Photoreceptor",
        "markers": ["OPN1LW", "OPN1SW", "ARR3", "GNAT2"],
        "species": ["Macaca fascicularis"],
        "source": "peng2019",
    },
    {
        "name": "Peripheral_Rods",
        "parent": "Rod_Photoreceptor",
        "markers": ["RHO", "NRL", "NR2E3", "GNAT1"],
        "species": ["Macaca fascicularis"],
        "source": "peng2019",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RHO": 1.0, "NRL": 1.0},
        },
        "action": "Rod_Photoreceptor",
        "source": "peng2019",
        "notes": "RHO+NRL co-expression specifies rod photoreceptor identity",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ARR3": 1.0, "GNAT2": 1.0},
        },
        "action": "Cone_Photoreceptor",
        "source": "peng2019",
        "notes": "ARR3+GNAT2 co-expression specifies cone photoreceptor identity",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RLBP1": 1.0, "GLUL": 1.0},
        },
        "action": "Muller_Glia",
        "source": "peng2019",
        "notes": "RLBP1+GLUL co-expression defines Muller glial cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "RGC",
        "source": "peng2019",
        "notes": "POU4F2+NEFM co-expression defines retinal ganglion cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"AIF1": 1.0, "P2RY12": 1.0},
        },
        "action": "Microglia",
        "source": "peng2019",
        "notes": "AIF1+P2RY12 co-expression defines retinal microglia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "peng2019",
        "notes": "ONECUT1+PROX1 co-expression defines horizontal cells",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"GAD1": 1.0, "GAD2": 1.0},
        },
        "action": "Amacrine_Cell",
        "source": "peng2019",
        "notes": "GAD1+GAD2 expression is characteristic of amacrine cells (GABAergic)",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"VSX1": 1.0, "PRKCA": 1.0},
        },
        "action": "Bipolar_Cell",
        "source": "peng2019",
        "notes": "VSX1+PRKCA co-expression defines bipolar cells",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "Rod_Photoreceptor",
            "marker": "NR2E3",
        },
        "type_b": {
            "cell_type": "Cone_Photoreceptor",
            "marker": "NR2E3",
        },
        "notes": "NR2E3 is rod-specific in this study but some sources report low expression in developing cones",
        "source": {"a": "peng2019", "b": "yan2020"},
    },
    {
        "type_a": {
            "cell_type": "Amacrine_Cell",
            "marker": "CALB1",
        },
        "type_b": {
            "cell_type": "Horizontal_Cell",
            "marker": "CALB1",
        },
        "notes": "CALB1 expressed in both amacrine and horizontal cells; not specific alone",
        "source": {"a": "peng2019", "b": "peng2020"},
    },
]
