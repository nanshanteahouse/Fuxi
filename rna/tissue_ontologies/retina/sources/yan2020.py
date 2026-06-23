"""
yan2020.py — Yan/Clark et al. 2020 (Cell Reports, PMID: 32386599)
Single-cell analysis of human retina identifies evolutionarily conserved and
species-specific mechanisms controlling development.

Species: Homo sapiens (fetal, GW9-27, postnatal day 8, adult 86y)
Tissue: retina (developing + adult)
Platform: 10X Genomics v2 (GSE116106, GSE122970, GSE138002)
Cells: >45,000 (across 18 time points)
Groups: RPCs, RGCs, Cones, Rods, Bipolar, Horizontal, Amacrine,
        Muller glia, Astrocytes, Microglia + intermediate precursors
"""

source_meta = {
    "id": "yan2020",
    "short_name": "Yan 2020 Cell Rep",
    "pmid": "32386599",
    "journal": "Cell Reports",
    "year": 2020,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["fetal retina", "organoid"],
    "n_cells": 45000,
    "n_subtypes": 15,
    "n_groups": 10,
}

markers = {
    "RPC": {
        "confirm": {
            "VSX2": ["32386599"],
            "PAX6": ["32386599"],
            "SOX2": ["32386599"],
        },
        "add": {
            "HES1": ["32386599"],
            "NOTCH1": ["32386599"],
            "MKI67": ["32386599"],
        },
        "refine": {
            "MKI67": {
                "note": "proliferation marker; labels actively cycling RPCs in fetal retina",
                "threshold": "log2FC > 1.5",
                "pmid": "32386599",
            },
            "SOX2": {
                "note": "retinal progenitor marker; also in Muller glia later in development",
                "threshold": "log2FC > 1.0",
                "pmid": "32386599",
            },
        },
    },
    "Retinal_Ganglion_Cell": {
        "confirm": {
            "POU4F2": ["32386599"],
            "NEFM": ["32386599"],
            "ELAVL4": ["32386599"],
        },
        "add": {
            "DCX": ["32386599"],
            "TUBB3": ["32386599"],
            "GAP43": ["32386599"],
        },
        "refine": {
            "GAP43": {
                "note": "axon growth marker; high in developing RGCs during axon extension",
                "threshold": "log2FC > 1.5",
                "pmid": "32386599",
            },
        },
    },
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["32386599"],
            "GNAT2": ["32386599"],
            "OPN1SW": ["32386599"],
        },
        "add": {
            "PDE6C": ["32386599"],
            "RCVRN": ["32386599"],
            "THRB": ["32386599"],
        },
        "refine": {
            "THRB": {
                "note": "thyroid hormone receptor; specifies S-cone fate during development",
                "threshold": "log2FC > 1.0",
                "pmid": "32386599",
            },
        },
    },
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["32386599"],
            "NRL": ["32386599"],
            "GNAT1": ["32386599"],
        },
        "add": {
            "NR2E3": ["32386599"],
            "PDE6B": ["32386599"],
        },
        "refine": {
            "NRL": {
                "note": "rod-determining transcription factor; expressed before RHO in development",
                "threshold": "log2FC > 1.0",
                "pmid": "32386599",
            },
        },
    },
    "Bipolar_Cell": {
        "confirm": {
            "VSX1": ["32386599"],
            "PRKCA": ["32386599"],
            "TRPM1": ["32386599"],
        },
        "add": {
            "GRM6": ["32386599"],
            "CABP5": ["32386599"],
            "OTX2": ["32386599"],
        },
        "refine": {
            "OTX2": {
                "note": "photoreceptor+bipolar cell fate determinant in development",
                "threshold": "log2FC > 1.0",
                "pmid": "32386599",
            },
        },
    },
    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["32386599"],
            "GAD2": ["32386599"],
            "TFAP2A": ["32386599"],
        },
        "add": {
            "TFAP2B": ["32386599"],
            "SLC6A9": ["32386599"],
            "CALB2": ["32386599"],
        },
        "refine": {},
    },
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["32386599"],
            "ONECUT2": ["32386599"],
            "PROX1": ["32386599"],
        },
        "add": {
            "CALB1": ["32386599"],
            "LHX1": ["32386599"],
        },
        "refine": {},
    },
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["32386599"],
            "GFAP": ["32386599"],
            "GLUL": ["32386599"],
        },
        "add": {
            "VIM": ["32386599"],
            "SLC1A3": ["32386599"],
            "CRYAB": ["32386599"],
        },
        "refine": {},
    },
    "Microglia": {
        "confirm": {
            "AIF1": ["32386599"],
            "P2RY12": ["32386599"],
            "CSF1R": ["32386599"],
        },
        "add": {
            "ITGAM": ["32386599"],
            "CX3CR1": ["32386599"],
            "CD74": ["32386599"],
        },
        "refine": {},
    },
    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["32386599"],
            "VWF": ["32386599"],
            "CDH5": ["32386599"],
        },
        "add": {
            "CLDN5": ["32386599"],
            "EGFL7": ["32386599"],
        },
        "refine": {},
    },
    "Pericyte": {
        "confirm": {
            "PDGFRB": ["32386599"],
            "CSPG4": ["32386599"],
            "ACTA2": ["32386599"],
        },
        "add": {
            "RGS5": ["32386599"],
        },
        "refine": {},
    },
    "RPE": {
        "confirm": {
            "RPE65": ["32386599"],
            "BEST1": ["32386599"],
            "LRAT": ["32386599"],
        },
        "add": {
            "TIMP3": ["32386599"],
        },
        "refine": {},
    },
    "Fibroblast": {
        "confirm": {
            "COL1A1": ["32386599"],
            "DCN": ["32386599"],
            "LUM": ["32386599"],
        },
        "add": {
            "COL3A1": ["32386599"],
        },
        "refine": {},
    },
    "Astrocyte": {
        "confirm": {
            "GFAP": ["32386599"],
            "AQP4": ["32386599"],
            "ALDH1L1": ["32386599"],
        },
        "add": {
            "S100B": ["32386599"],
        },
        "refine": {},
    },
    "Oligodendrocyte": {
        "confirm": {
            "MBP": ["32386599"],
            "PLP1": ["32386599"],
            "MOG": ["32386599"],
        },
        "add": {
            "OLIG2": ["32386599"],
        },
        "refine": {},
    },
}

novel_types = [
    {
        "name": "Developing_AC_HC_Precursors",
        "parent": "Amacrine_Cell",
        "markers": ["TFAP2A", "TFAP2B", "ONECUT1", "PROX1"],
        "species": ["Homo sapiens"],
        "source": "yan2020",
    },
    {
        "name": "Developing_BC_Photo_Precursors",
        "parent": "Bipolar_Cell",
        "markers": ["VSX1", "OTX2", "CRX", "RCVRN"],
        "species": ["Homo sapiens"],
        "source": "yan2020",
    },
]

expert_rules = [
    {
        "priority": 1,
        "condition": {
            "markers_present": {"VSX2": 1.0, "SOX2": 1.0},
        },
        "action": "RPC",
        "source": "yan2020",
        "notes": "VSX2+SOX2 co-expression specifies retinal progenitor identity in fetal retina",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"POU4F2": 1.0, "NEFM": 1.0},
        },
        "action": "Retinal_Ganglion_Cell",
        "source": "yan2020",
        "notes": "POU4F2+NEFM co-expression defines developing RGCs",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RHO": 1.0, "NRL": 1.0},
        },
        "action": "Rod_Photoreceptor",
        "source": "yan2020",
        "notes": "RHO+NRL co-expression specifies developing rod photoreceptors",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ARR3": 1.0, "GNAT2": 1.0},
        },
        "action": "Cone_Photoreceptor",
        "source": "yan2020",
        "notes": "ARR3+GNAT2 co-expression specifies developing cone photoreceptors",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"RLBP1": 1.0, "GLUL": 1.0},
        },
        "action": "Muller_Glia",
        "source": "yan2020",
        "notes": "RLBP1+GLUL co-expression defines developing Muller glia",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"ONECUT1": 1.0, "PROX1": 1.0},
        },
        "action": "Horizontal_Cell",
        "source": "yan2020",
        "notes": "ONECUT1+PROX1 co-expression defines developing horizontal cells",
    },
    {
        "priority": 1,
        "condition": {
            "markers_present": {"AIF1": 1.0, "P2RY12": 1.0},
        },
        "action": "Microglia",
        "source": "yan2020",
        "notes": "AIF1+P2RY12 co-expression defines fetal retinal microglia",
    },
    {
        "priority": 2,
        "condition": {
            "markers_present": {"DCX": 1.0, "TUBB3": 1.0},
        },
        "action": "Retinal_Ganglion_Cell",
        "source": "yan2020",
        "notes": "DCX+TUBB3 marks differentiating neurons; most abundant in early RGCs in fetal retina",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "RPC",
            "marker": "MKI67",
        },
        "type_b": {
            "cell_type": "Proliferating_Cell",
            "marker": "MKI67",
        },
        "notes": "MKI67 is a general proliferation marker. In fetal retina, proliferating cells are predominantly RPCs, not a separate cell type",
        "source": {"a": "yan2020", "b": "hoang2023"},
    },
]
