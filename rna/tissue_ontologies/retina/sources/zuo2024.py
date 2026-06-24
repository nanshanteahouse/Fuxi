"""
zuo2024.py — Zuo et al. 2024 (Nature Communications, PMID: 39117640)
Single cell dual-omic atlas of the human developing retina.

Species: Homo sapiens (fetal, PCW 8-23)
Tissue: developing retina (foveal + peripheral, matched pairs)
Platform: 10X Multiome (snRNA-seq + snATAC-seq from same nuclei)
Nuclei: ~220,000 (14 donors)
Cell classes: 9 major classes (PRPC, NRPC, MG, RGC, HC, AC, BC, Cone, Rod)
GEO: GSE268630

Marker sources:
  - Supplementary Data 6 (top DE genes per major cell class, Wilcoxon rank-sum).
    Top 30 genes by score per class; genes in top 10 with p_adj=0 placed in
    `confirm`, the remainder in `add`.
  - Supplementary Data 11 (NRPC TF fate prediction table — 95 TFs with
    fate assignments validated against literature).
  - Supplementary Data 3 (Ground Truth: 55 TFs with experimental validation
    from loss-of-function animal models).
"""

source_meta = {
    "id": "zuo2024",
    "short_name": "Zuo 2024 Nat Commun",
    "pmid": "39117640",
    "journal": "Nature Communications",
    "year": 2024,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["fovea", "periphery", "developing"],
    "n_cells": 220000,
    "n_subtypes": 22,
    "n_groups": 9,
    "class": "Mammalia",
    "order": "Primates",
}

markers = {
    # ── PRPC: Proliferating Retinal Progenitor Cells ────────────────
    "PRPC": {
        "confirm": {
            "VSX2": ["39117640"],
            "SOX2": ["39117640"],
            "PAX6": ["39117640"],
        },
        "add": {
            "SORCS1": ["39117640"],
            "ZHX2": ["39117640"],
            "PARD3B": ["39117640"],
            "GLI3": ["39117640"],
            "PCDH11X": ["39117640"],
            "SHROOM3": ["39117640"],
            "PTPRK": ["39117640"],
            "ADAMTS18": ["39117640"],
            "PIP5K1B": ["39117640"],
            "MOB3B": ["39117640"],
            "ELL2": ["39117640"],
            "NFIA": ["39117640"],
            "CDK6": ["39117640"],
            "ASPM": ["39117640"],
            "MKI67": ["39117640"],
            "MXD3": ["39117640"],
            "NPAS3": ["39117640"],
            "ZNF367": ["39117640"],
            "MECOM": ["39117640"],
            "PROX1": ["39117640"],
            "FOXN4": ["39117640"],
            "ECT2": ["39117640"],
        },
        "refine": {
            "ARHGAP11B": {
                "note": "Human-specific gene expressed in transitioning PRPCs; "
                       "no mouse homolog. Identified in zuo2024 NRPC GRN.",
                "threshold": "log2FC > 0.5",
                "pmid": "39117640",
            },
        },
    },

    # ── NRPC: Neurogenic RPCs ───────────────────────────────────────
    "NRPC": {
        "confirm": {
            "ATOH7": ["39117640"],
            "PRDM13": ["39117640"],
            "OTX2": ["39117640"],
        },
        "add": {
            "NEUROD1": ["39117640"],
            "HES6": ["39117640"],
            "TFAP2A": ["39117640"],
            "TFAP2B": ["39117640"],
            "CRX": ["39117640"],
            "NRL": ["39117640"],
            "VSX1": ["39117640"],
            "ISL1": ["39117640"],
            "POU4F2": ["39117640"],
            "THRB": ["39117640"],
            "PRDM1": ["39117640"],
            "ONECUT1": ["39117640"],
            "ONECUT2": ["39117640"],
        },
        "refine": {
            "ATOH7": {
                "note": "Master regulator of RGC fate. Ground Truth validated "
                       "(mouse Atoh7 KO → complete loss of RGCs). PMID: 11723443",
                "threshold": "log2FC > 1.0",
                "pmid": "39117640",
            },
        },
    },

    # ── Rod Photoreceptor (developing + adult marker overlap) ──────
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["39117640"],
            "NRL": ["39117640"],
            "NR2E3": ["39117640"],
        },
        "add": {
            "EYS": ["39117640"],
            "ANO2": ["39117640"],
            "RP1": ["39117640"],
            "EPS8": ["39117640"],
            "TMEM244": ["39117640"],
            "RCVRN": ["39117640"],
            "USH2A": ["39117640"],
            "PDC": ["39117640"],
            "PDE6B": ["39117640"],
            "GNAT1": ["39117640"],
            "CRX": ["39117640"],
            "SAG": ["39117640"],
            "PDE6A": ["39117640"],
        },
        "refine": {
            "NR2E3": {
                "note": "Rod fate TF; validated in Ground Truth (ESCS patients "
                       "with NR2E3 mutations, PMID: 18547563)",
                "threshold": "log2FC > 1.0",
                "pmid": "39117640",
            },
        },
    },

    # ── Cone Photoreceptor ──────────────────────────────────────────
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["39117640"],
            "GNAT2": ["39117640"],
            "THRB": ["39117640"],
        },
        "add": {
            "EYS": ["39117640"],
            "PRDM1": ["39117640"],
            "PDE1C": ["39117640"],
            "EGFLAM": ["39117640"],
            "DCT": ["39117640"],
            "TMEM108": ["39117640"],
            "MCC": ["39117640"],
            "WWC1": ["39117640"],
            "ANKRD33B": ["39117640"],
            "OPN1SW": ["39117640"],
            "PDE6H": ["39117640"],
            "CNGB3": ["39117640"],
            "GUCY2E": ["39117640"],
        },
        "refine": {
            "THRB": {
                "note": "Cone fate TF; validated in Ground Truth (mouse TrB2 "
                       "KO → complete loss of cones, PMID: 11138006)",
                "threshold": "log2FC > 1.0",
                "pmid": "39117640",
            },
        },
    },

    # ── Bipolar Cell ────────────────────────────────────────────────
    "Bipolar_Cell": {
        "confirm": {
            "VSX2": ["39117640"],
            "GRM6": ["39117640"],
            "TRPM1": ["39117640"],
        },
        "add": {
            "CA10": ["39117640"],
            "FSTL5": ["39117640"],
            "IGSF21": ["39117640"],
            "GRIK1": ["39117640"],
            "VSX1": ["39117640"],
            "PRDM8": ["39117640"],
            "OTX2": ["39117640"],
            "PTPRR": ["39117640"],
            "GABRR3": ["39117640"],
            "TMEM215": ["39117640"],
            "ISL1": ["39117640"],
            "PRKCA": ["39117640"],
            "CABP5": ["39117640"],
            "PKC": ["39117640"],
            "SLC4A10": ["39117640"],
        },
    },

    # ── Amacrine Cell ───────────────────────────────────────────────
    "Amacrine_Cell": {
        "confirm": {
            "TFAP2A": ["39117640"],
            "GAD1": ["39117640"],
            "GAD2": ["39117640"],
        },
        "add": {
            "MYT1L": ["39117640"],
            "PDE4D": ["39117640"],
            "MDGA2": ["39117640"],
            "ASIC2": ["39117640"],
            "NRP1": ["39117640"],
            "CACNA1D": ["39117640"],
            "TFAP2B": ["39117640"],
            "DPP6": ["39117640"],
            "TENM2": ["39117640"],
            "CALB2": ["39117640"],
            "SLC32A1": ["39117640"],
            "SLC6A9": ["39117640"],
            "GLRA1": ["39117640"],
        },
    },

    # ── Horizontal Cell ─────────────────────────────────────────────
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["39117640"],
            "PROX1": ["39117640"],
            "CALB1": ["39117640"],
        },
        "add": {
            "CNTN4": ["39117640"],
            "ONECUT2": ["39117640"],
            "STMN2": ["39117640"],
            "TFAP2B": ["39117640"],
            "NDST3": ["39117640"],
            "ONECUT3": ["39117640"],
            "MEGF10": ["39117640"],
            "NTRK2": ["39117640"],
            "SYN3": ["39117640"],
            "LHX1": ["39117640"],
            "EBF1": ["39117640"],
        },
    },

    # ── Retinal Ganglion Cell ───────────────────────────────────────
    "RGC": {
        "confirm": {
            "POU4F2": ["39117640"],
            "RBPMS": ["39117640"],
            "NEFM": ["39117640"],
        },
        "add": {
            "EBF1": ["39117640"],
            "POU6F2": ["39117640"],
            "EBF3": ["39117640"],
            "GAP43": ["39117640"],
            "ELAVL4": ["39117640"],
            "ELAVL2": ["39117640"],
            "KLHL1": ["39117640"],
            "RBFOX3": ["39117640"],
            "ISL1": ["39117640"],
            "POU4F1": ["39117640"],
            "NEFL": ["39117640"],
            "SLC17A6": ["39117640"],
            "THY1": ["39117640"],
            "SNCG": ["39117640"],
            "ATOH7": ["39117640"],
        },
        "refine": {
            "POU4F2": {
                "note": "RGC master TF; Ground Truth validated (mouse Pou4f2 "
                       "KO → loss of most RGCs, PMID: 25775587)",
                "threshold": "log2FC > 1.0",
                "pmid": "39117640",
            },
        },
    },

    # ── Müller Glia ─────────────────────────────────────────────────
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["39117640"],
            "GLUL": ["39117640"],
            "SLC1A3": ["39117640"],
        },
        "add": {
            "PARD3B": ["39117640"],
            "RARB": ["39117640"],
            "SPON1": ["39117640"],
            "PCDH11X": ["39117640"],
            "SORCS1": ["39117640"],
            "SAT1": ["39117640"],
            "NFIA": ["39117640"],
            "ADAMTS6": ["39117640"],
            "GLI3": ["39117640"],
            "CRYM": ["39117640"],
            "CLU": ["39117640"],
            "VIM": ["39117640"],
            "AQP4": ["39117640"],
            "DKK3": ["39117640"],
            "SOX2": ["39117640"],
        },
    },

    # ── Microglia ───────────────────────────────────────────────────
    "Microglia": {
        "confirm": {
            "AIF1": ["39117640"],
            "P2RY12": ["39117640"],
            "CSF1R": ["39117640"],
        },
        "add": {
            "CX3CR1": ["39117640"],
            "ITGAM": ["39117640"],
            "TREM2": ["39117640"],
            "C1QA": ["39117640"],
            "C1QB": ["39117640"],
            "TYROBP": ["39117640"],
        },
    },

    # ── Astrocyte ───────────────────────────────────────────────────
    "Astrocyte": {
        "confirm": {
            "GFAP": ["39117640"],
            "S100B": ["39117640"],
            "ALDH1L1": ["39117640"],
        },
        "add": {
            "AQP4": ["39117640"],
            "GJA1": ["39117640"],
            "SLC1A2": ["39117640"],
        },
    },

    # ── RPE ─────────────────────────────────────────────────────────
    "RPE": {
        "confirm": {
            "RPE65": ["39117640"],
            "BEST1": ["39117640"],
            "LRAT": ["39117640"],
        },
        "add": {
            "TIMP3": ["39117640"],
            "RDH10": ["39117640"],
            "TYRP1": ["39117640"],
            "PMEL": ["39117640"],
            "DCT": ["39117640"],
        },
    },
}

novel_types = [
    {
        "name": "PRPC",
        "parent": "RPC",
        "markers": ["VSX2", "SOX2", "PAX6", "MKI67", "ASPM", "MXD3", "NPAS3",
                    "ZNF367", "PROX1", "CDK6", "GLI3", "SORCS1"],
        "species": ["Homo sapiens"],
        "source": "zuo2024",
    },
    {
        "name": "NRPC",
        "parent": "RPC",
        "markers": ["ATOH7", "PRDM13", "OTX2", "NEUROD1", "HES6", "TFAP2A",
                    "CRX", "NRL", "POU4F2", "ISL1"],
        "species": ["Homo sapiens"],
        "source": "zuo2024",
    },
]

expert_rules = [
    {
        "priority": 10,
        "condition": {
            "markers_present": {
                "VSX2": 1.0,
                "SOX2": 1.0,
                "MKI67": 0.5,
            },
            "markers_absent": ["ATOH7", "OTX2"],
        },
        "action": "PRPC",
        "source": "zuo2024",
        "notes": "VSX2+SOX2+MKI67 without neurogenic TFs specifies proliferating RPCs",
    },
    {
        "priority": 9,
        "condition": {
            "markers_present": {
                "ATOH7": 0.5,
                "POU4F2": 1.0,
                "ISL1": 0.5,
            },
        },
        "action": "RGC",
        "source": "zuo2024",
        "notes": "ATOH7+POU4F2+ISL1: core RGC identity GRN from developing retina",
    },
    {
        "priority": 8,
        "condition": {
            "markers_present": {
                "ONECUT1": 1.0,
                "PROX1": 0.5,
                "CALB1": 0.5,
            },
        },
        "action": "Horizontal_Cell",
        "source": "zuo2024",
        "notes": "ONECUT1+PROX1: experimentally validated HC specification (PMID: 25228773)",
    },
    {
        "priority": 7,
        "condition": {
            "markers_present": {
                "RLBP1": 1.0,
                "GLUL": 1.0,
            },
            "markers_absent": ["MKI67", "ATOH7"],
        },
        "action": "Muller_Glia",
        "source": "zuo2024",
        "notes": "RLBP1+GLUL without progenitor markers distinguishes mature MG from PRPCs",
    },
]

conflicts = []
