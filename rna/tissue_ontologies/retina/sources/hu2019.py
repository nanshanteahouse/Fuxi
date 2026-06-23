"""
hu2019.py — Hu et al. 2019 (PLoS Biology, PMID: 31269016)
Dissecting the transcriptome landscape of the human fetal neural retina
and retinal pigment epithelium by single-cell RNA-seq analysis.

Species: Homo sapiens (fetal, 5-24 GW)
Tissue: fetal neural retina + RPE
Platform: Modified Smart-seq2 (384-well, mouth pipette single-cell picking)
Cells: 2,421
Groups: 9 major types
GEO: GSE107618
"""

source_meta = {
    "id": "hu2019",
    "short_name": "Hu 2019 PLoS Biol",
    "pmid": "31269016",
    "journal": "PLoS Biology",
    "year": 2019,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["fetal retina", "fetal RPE"],
    "n_cells": 2421,
    "n_subtypes": 9,
    "n_groups": 9,
}

markers = {
    "RPE": {
        "confirm": {
            "RPE65": ["31269016"],
            "BEST1": ["31269016"],
            "LRAT": ["31269016"],
        },
        "add": {
            "TIMP3": ["31269016"],
            "RDH10": ["31269016"],
            "CRALBP": ["31269016"],
        },
        "refine": {
            "RPE65": {
                "note": "RPE-specific marker in fetal tissue; high expression by GW8",
                "threshold": "log2FC > 2.0",
                "pmid": "31269016",
            },
        },
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["31269016"],
            "DCN": ["31269016"],
            "LUM": ["31269016"],
        },
        "add": {
            "COL3A1": ["31269016"],
            "COL1A2": ["31269016"],
            "PDGFRA": ["31269016"],
        },
        "refine": {},
    },
    "RPC": {
        "confirm": {
            "VSX2": ["31269016"],
            "PAX6": ["31269016"],
            "SOX2": ["31269016"],
        },
        "add": {
            "HES1": ["31269016"],
            "NOTCH1": ["31269016"],
            "MKI67": ["31269016"],
        },
        "refine": {
            "VSX2": {
                "note": "canonical RPC marker; expressed throughout fetal stages 5-24 GW",
                "threshold": "log2FC > 1.0",
                "pmid": "31269016",
            },
        },
    },
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["31269016"],
            "GFAP": ["31269016"],
            "GLUL": ["31269016"],
        },
        "add": {
            "VIM": ["31269016"],
            "SLC1A3": ["31269016"],
        },
        "refine": {},
    },
    "Retinal_Ganglion_Cell": {
        "confirm": {
            "POU4F2": ["31269016"],
            "NEFM": ["31269016"],
            "ELAVL4": ["31269016"],
        },
        "add": {
            "TUBB3": ["31269016"],
            "DCX": ["31269016"],
        },
        "refine": {},
    },
    "Photoreceptor_Precursor": {
        "confirm": {
            "CRX": ["31269016"],
            "RCVRN": ["31269016"],
            "OTX2": ["31269016"],
        },
        "add": {
            "NRL": ["31269016"],
            "RHO": ["31269016"],
        },
        "refine": {
            "CRX": {
                "note": "early photoreceptor precursor marker; precedes opsin expression",
                "threshold": "log2FC > 1.0",
                "pmid": "31269016",
            },
        },
    },
    "Bipolar_Precursor": {
        "confirm": {
            "VSX1": ["31269016"],
            "OTX2": ["31269016"],
            "VSX2": ["31269016"],
        },
        "add": {
            "PRKCA": ["31269016"],
            "GRM6": ["31269016"],
        },
        "refine": {},
    },
    "Amacrine_Precursor": {
        "confirm": {
            "GAD1": ["31269016"],
            "TFAP2A": ["31269016"],
            "TFAP2B": ["31269016"],
        },
        "add": {
            "CALB2": ["31269016"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["31269016"],
            "PROX1": ["31269016"],
            "CALB1": ["31269016"],
        },
        "add": {
            "ONECUT2": ["31269016"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["31269016"],
            "VWF": ["31269016"],
            "CDH5": ["31269016"],
        },
        "add": {
            "CLDN5": ["31269016"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["31269016"],
            "CSF1R": ["31269016"],
            "CX3CR1": ["31269016"],
        },
        "add": {
            "ITGAM": ["31269016"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "Fetal_RPE",
        "parent": "RPE",
        "markers": ["RPE65", "BEST1", "LRAT", "TIMP3"],
        "species": ["Homo sapiens"],
        "source": "hu2019",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RPE65": 1.0, "BEST1": 1.0},
        },
        "action": "RPE",
        "source": "hu2019",
        "notes": "RPE65+BEST1 co-expression defines fetal RPE cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"VSX2": 1.0, "SOX2": 1.0},
        },
        "action": "RPC",
        "source": "hu2019",
        "notes": "VSX2+SOX2 co-expression defines fetal retinal progenitor cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"CRX": 1.0, "RCVRN": 1.0},
        },
        "action": "Photoreceptor_Precursor",
        "source": "hu2019",
        "notes": "CRX+RCVRN expression defines photoreceptor precursors in fetal retina",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"COL1A1": 1.0, "DCN": 1.0},
        },
        "action": "Fibroblast",
        "source": "hu2019",
        "notes": "COL1A1+DCN co-expression defines fetal fibroblasts",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "Retinal_Ganglion_Cell",
        "source": "hu2019",
        "notes": "POU4F2+NEFM co-expression defines fetal RGCs",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "hu2019",
        "notes": "ONECUT1+PROX1 co-expression defines fetal horizontal cells",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"GAD1": 1.0, "TFAP2A": 1.0},
        },
        "action": "Amacrine_Precursor",
        "source": "hu2019",
        "notes": "GAD1+TFAP2A expression defines developing amacrine cells",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"VSX1": 1.0, "OTX2": 1.0},
        },
        "action": "Bipolar_Precursor",
        "source": "hu2019",
        "notes": "VSX1+OTX2 expression defines developing bipolar cells",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "Muller_Glia",
            "marker": "GFAP",
        },
        "type_b": {
            "cell_type": "Glial_Cell_Oligodendrocyte",
            "marker": "GFAP",
        },
        "notes": "GFAP is expressed by both Muller glia and astrocytes; this study labels all GFAP+ neural retinal cells as Muller glia",
        "source": {"a": "hu2019", "b": "yan2020"},
    },
]
