"""
hoang2023.py — Hoang et al. 2023 (preprint, GSE246169)
ASCL1 induces neurogenesis in human Muller glia.

Species: Homo sapiens (fetal D59, D76 + ASCL1-overexpressing Muller glia culture)
Tissue: fetal retina + in vitro MG reprogramming system
Platform: 10X Multiome (snRNA-seq + snATAC-seq)
Cells: moderate scale
GEO: GSE246169
"""

"""

# NOTE (2026-06-24): This source has been DISABLED from KB merging.

source_meta = {
    "id": "hoang2023",
    "short_name": "Hoang 2023 ASCL1",
    "pmid": "00000000",
    "journal": "preprint / HRCA",
    "year": 2023,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["fetal retina", "MG culture"],
    "n_cells": 50000,
    "n_subtypes": 18,
    "n_groups": 9,
    "class": "Mammalia",
    "order": "Primates",
}

markers = {
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["GSE246169"],
            "GLUL": ["GSE246169"],
            "VIM": ["GSE246169"],
        },
        "add": {
            "SLC1A3": ["GSE246169"],
            "CRYAB": ["GSE246169"],
            "CLU": ["GSE246169"],
        },
        "refine": {
            "ASCL1": {
                "note": "ectopic ASCL1 expression drives MG reprogramming toward neuronal fate",
                "threshold": "log2FC > 2.0 (overexpression)",
                "pmid": "GSE246169",
            },
            "RLBP1": {
                "note": "canonical MG marker; downregulated upon ASCL1-induced reprogramming",
                "threshold": "log2FC > 1.0",
                "pmid": "GSE246169",
            },
        },
    },
    "RPC": {
        "confirm": {
            "VSX2": ["GSE246169"],
            "PAX6": ["GSE246169"],
            "SOX2": ["GSE246169"],
        },
        "add": {
            "HES1": ["GSE246169"],
            "NOTCH1": ["GSE246169"],
            "MKI67": ["GSE246169"],
        },
        "refine": {
            "MKI67": {
                "note": "proliferation marker; high in proliferating RPCs and ASCL1-reprogrammed MG",
                "threshold": "log2FC > 1.5",
                "pmid": "GSE246169",
            },
        },
    },
    "Proliferating_MG": {
        "confirm": {
            "MKI67": ["GSE246169"],
            "TOP2A": ["GSE246169"],
            "PCNA": ["GSE246169"],
        },
        "add": {
            "ASCL1": ["GSE246169"],
            "DLL1": ["GSE246169"],
            "HES6": ["GSE246169"],
        },
        "refine": {
            "ASCL1": {
                "note": "proneural factor; induces MG re-entry into cell cycle and neurogenesis",
                "threshold": "log2FC > 1.5",
                "pmid": "GSE246169",
            },
        },
    },
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["GSE246169"],
            "NRL": ["GSE246169"],
            "GNAT1": ["GSE246169"],
        },
        "add": {
            "PDE6B": ["GSE246169"],
            "NR2E3": ["GSE246169"],
        },
        "refine": {},
    },
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["GSE246169"],
            "GNAT2": ["GSE246169"],
            "RCVRN": ["GSE246169"],
        },
        "add": {
            "PDE6C": ["GSE246169"],
            "THRB": ["GSE246169"],
        },
        "refine": {},
    },
    "RGC": {
        "confirm": {
            "POU4F2": ["GSE246169"],
            "NEFM": ["GSE246169"],
            "ELAVL4": ["GSE246169"],
        },
        "add": {
            "GAP43": ["GSE246169"],
            "TUBB3": ["GSE246169"],
        },
        "refine": {},
    },
    "Bipolar_Cell": {
        "confirm": {
            "VSX1": ["GSE246169"],
            "PRKCA": ["GSE246169"],
            "TRPM1": ["GSE246169"],
        },
        "add": {
            "GRM6": ["GSE246169"],
            "OTX2": ["GSE246169"],
        },
        "refine": {},
    },
    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["GSE246169"],
            "GAD2": ["GSE246169"],
            "TFAP2A": ["GSE246169"],
        },
        "add": {
            "CALB2": ["GSE246169"],
            "SLC6A9": ["GSE246169"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["GSE246169"],
            "ONECUT2": ["GSE246169"],
            "PROX1": ["GSE246169"],
        },
        "add": {
            "CALB1": ["GSE246169"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["GSE246169"],
            "CSF1R": ["GSE246169"],
            "P2RY12": ["GSE246169"],
        },
        "add": {
            "ITGAM": ["GSE246169"],
            "CX3CR1": ["GSE246169"],
        },
        "refine": {},
    },
    "RPE": {
        "confirm": {
            "RPE65": ["GSE246169"],
            "BEST1": ["GSE246169"],
            "LRAT": ["GSE246169"],
        },
        "add": {
            "TIMP3": ["GSE246169"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["GSE246169"],
            "VWF": ["GSE246169"],
            "CDH5": ["GSE246169"],
        },
        "add": {
            "CLDN5": ["GSE246169"],
        },
        "refine": {},
    },
    "Pericyte": {
        "confirm": {
            "PDGFRB": ["GSE246169"],
            "CSPG4": ["GSE246169"],
            "ACTA2": ["GSE246169"],
        },
        "add": {
            "RGS5": ["GSE246169"],
        },
        "refine": {},
    },
    "Astrocyte": {
        "confirm": {
            "GFAP": ["GSE246169"],
            "AQP4": ["GSE246169"],
            "ALDH1L1": ["GSE246169"],
        },
        "add": {
            "S100B": ["GSE246169"],
        },
        "refine": {},
    },
    "Oligodendrocyte": {
        "confirm": {
            "MBP": ["GSE246169"],
            "PLP1": ["GSE246169"],
            "MOG": ["GSE246169"],
        },
        "add": {
            "OLIG2": ["GSE246169"],
        },
        "refine": {},
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["GSE246169"],
            "DCN": ["GSE246169"],
            "LUM": ["GSE246169"],
        },
        "add": {
            "COL3A1": ["GSE246169"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "ASCL1_Reprogrammed_MG",
        "parent": "Muller_Glia",
        "markers": ["ASCL1", "MKI67", "TOP2A", "DLL1", "HES6"],
        "species": ["Homo sapiens"],
        "source": "hoang2023",
    },
    {
        "name": "Proliferating_RPC",
        "parent": "RPC",
        "markers": ["MKI67", "TOP2A", "PCNA", "HES1"],
        "species": ["Homo sapiens"],
        "source": "hoang2023",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RLBP1": 1.0, "GLUL": 1.0},
        },
        "action": "Muller_Glia",
        "source": "hoang2023",
        "notes": "RLBP1+GLUL co-expression defines Muller glia in fetal retina and culture",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"VSX2": 1.0, "SOX2": 1.0},
        },
        "action": "RPC",
        "source": "hoang2023",
        "notes": "VSX2+SOX2 co-expression defines fetal retinal progenitor cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RHO": 1.0, "NRL": 1.0},
        },
        "action": "Rod_Photoreceptor",
        "source": "hoang2023",
        "notes": "RHO+NRL co-expression specifies rod photoreceptors",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ARR3": 1.0, "GNAT2": 1.0},
        },
        "action": "Cone_Photoreceptor",
        "source": "hoang2023",
        "notes": "ARR3+GNAT2 co-expression specifies cone photoreceptors",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "RGC",
        "source": "hoang2023",
        "notes": "POU4F2+NEFM co-expression defines RGCs",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RPE65": 1.0, "BEST1": 1.0},
        },
        "action": "RPE",
        "source": "hoang2023",
        "notes": "RPE65+BEST1 co-expression defines RPE",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"AIF1": 1.0, "P2RY12": 1.0},
        },
        "action": "Microglia",
        "source": "hoang2023",
        "notes": "AIF1+P2RY12 co-expression defines retinal microglia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"MKI67": 1.0, "TOP2A": 1.0, "ASCL1": 1.0},
        },
        "action": "Proliferating_MG",
        "source": "hoang2023",
        "notes": "MKI67+TOP2A+ASCL1 marks ASCL1-reprogrammed Muller glia re-entering cell cycle",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "hoang2023",
        "notes": "ONECUT1+PROX1 co-expression defines horizontal cells",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "Proliferating_MG",
            "marker": "MKI67",
        },
        "type_b": {
            "cell_type": "Muller_Glia",
            "marker": "RLBP1",
        },
        "notes": "ASCL1-reprogrammed MG co-express proliferation markers (MKI67) while retaining some MG markers (RLBP1). This transitional state is uniquely identified in Hoang 2023",
        "source": {"a": "hoang2023", "b": "peng2020"},
    },
]

"""
