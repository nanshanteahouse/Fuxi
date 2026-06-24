"""
menon2019.py — Menon et al. 2019 (Nature Communications, PMID: 31653841)
Single-cell Transcriptomic Atlas of the Human Retina Identifies Cell Types
Associated with Age-Related Macular Degeneration.

Species: Homo sapiens (adult, 3 donors, 42-80 years)
Tissue: retina (macula + peripheral)
Platform: Microfluidics (Fluidigm C1, GSE137537) + Seq-Well (GSE137846)
Cells: ~20,000 (combined across two platforms)
Key findings: First comprehensive single-cell atlas of human retina;
              AMD GWAS enrichment in glia, vascular cells, and cone
              photoreceptors; macroglial heterogeneity identified.
GEO: GSE137847 (SuperSeries: GSE137537 + GSE137846)

Marker source: Supplementary Data 5 (filtered_gene_scores matrix).
    Genes with cell-type-specific score > 2.0 and negative in all other types
    are placed in `confirm`; genes with score > 1.0 are in `add`.
    Negative markers derived from genes with strong negative scores (< -1.0).
"""

source_meta = {
    "id": "menon2019",
    "short_name": "Menon 2019 Nat Commun",
    "pmid": "31653841",
    "journal": "Nature Communications",
    "year": 2019,
    "species": ["Homo sapiens"],
    "tissue": "retina",
    "regions": ["macula", "periphery"],
    "n_cells": 20000,
    "n_subtypes": 12,
    "n_groups": 9,
    "class": "Mammalia",
    "order": "Primates",
}

markers = {
    "Rod_Photoreceptor": {
        "confirm": {
            "PDE6A": ["31653841"],
            "NR2E3": ["31653841"],
            "CNGA1": ["31653841"],
            "PDE6B": ["31653841"],
        },
        "add": {
            "AIPL1": ["31653841"],
            "CNGB1": ["31653841"],
            "UNC119": ["31653841"],
            "MFGE8": ["31653841"],
            "EPB41L2": ["31653841"],
            "SCAPER": ["31653841"],
            "PDC": ["31653841"],
            "SYNE2": ["31653841"],
            "MOK": ["31653841"],
            "OSBP2": ["31653841"],
            "PPEF2": ["31653841"],
            "GNAT1": ["31653841"],
            "SAG": ["31653841"],
            "RCVRN": ["31653841"],
        },
        "refine": {
            "NR2E3": {
                "note": "Rod-fate specification TF; mutation causes Enhanced "
                       "S-Cone Syndrome (ESCS). Strongest rod-specific TF marker.",
                "threshold": "score > 10.0",
                "pmid": "31653841",
            },
        },
    },

    "Cone_Photoreceptor": {
        "confirm": {
            "OPN1LW": ["31653841"],
            "CNGB3": ["31653841"],
            "GRK7": ["31653841"],
            "RAB41": ["31653841"],
        },
        "add": {
            "CRX": ["31653841"],
            "MEGF9": ["31653841"],
            "GNGT2": ["31653841"],
            "ARR3": ["31653841"],
            "GNAT2": ["31653841"],
            "OPN1SW": ["31653841"],
            "PDE6H": ["31653841"],
            "THRB": ["31653841"],
        },
        "refine": {
            "OPN1LW": {
                "note": "L-opsin; strongest cone-specific signal in adult human "
                       "retina. L/M cone subtype; S-cones express OPN1SW instead.",
                "threshold": "score > 8.0",
                "pmid": "31653841",
            },
        },
    },

    "RGC": {
        "confirm": {
            "NEFM": ["31653841"],
            "POU4F2": ["31653841"],
            "RBPMS": ["31653841"],
        },
        "add": {
            "UNC119": ["31653841"],
            "AIPL1": ["31653841"],
            "STMN1": ["31653841"],
            "PCP2": ["31653841"],
            "LAPTM4B": ["31653841"],
            "POU4F1": ["31653841"],
            "NEFL": ["31653841"],
            "SLC17A6": ["31653841"],
            "THY1": ["31653841"],
            "ISL1": ["31653841"],
        },
    },

    "Bipolar_Cell": {
        "confirm": {
            "TRPM1": ["31653841"],
            "GRM6": ["31653841"],
            "VSX1": ["31653841"],
            "TMEM215": ["31653841"],
        },
        "add": {
            "RGS16": ["31653841"],
            "LRTM1": ["31653841"],
            "PRDM8": ["31653841"],
            "DOK6": ["31653841"],
            "EFR3A": ["31653841"],
            "CHN2": ["31653841"],
            "SLC38A1": ["31653841"],
            "OTX2": ["31653841"],
            "VSX2": ["31653841"],
            "ISL1": ["31653841"],
            "PRKCA": ["31653841"],
            "CABP5": ["31653841"],
        },
        "refine": {
            "TRPM1": {
                "note": "ON-bipolar cell marker; signal transduction channel. "
                       "Highest bipolar-specific signal in adult human retina.",
                "threshold": "score > 14.0",
                "pmid": "31653841",
            },
            "GRM6": {
                "note": "Metabotropic glutamate receptor; ON-bipolar synaptic "
                       "marker. Mutations cause congenital stationary night blindness.",
                "threshold": "score > 11.0",
                "pmid": "31653841",
            },
        },
    },

    "Amacrine_Cell": {
        "confirm": {
            "GAD1": ["31653841"],
            "ELAVL3": ["31653841"],
        },
        "add": {
            "C1QL2": ["31653841"],
            "CARTPT": ["31653841"],
            "GABRA2": ["31653841"],
            "MIAT": ["31653841"],
            "PAPPA": ["31653841"],
            "CACNA2D2": ["31653841"],
            "TFAP2A": ["31653841"],
            "GAD2": ["31653841"],
            "CALB2": ["31653841"],
            "SLC32A1": ["31653841"],
        },
    },

    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["31653841"],
            "LHX1": ["31653841"],
        },
        "add": {
            "ONECUT3": ["31653841"],
            "NDFIP1": ["31653841"],
            "FRMPD4": ["31653841"],
            "PROX1": ["31653841"],
            "CALB1": ["31653841"],
            "ONECUT2": ["31653841"],
        },
        "refine": {
            "ONECUT1": {
                "note": "Master TF for horizontal cell specification. "
                       "Strongest HC-specific signal (score > 9.0).",
                "threshold": "score > 9.0",
                "pmid": "31653841",
            },
        },
    },

    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["31653841"],
            "GLUL": ["31653841"],
            "CLU": ["31653841"],
        },
        "add": {
            "TF": ["31653841"],
            "GPX3": ["31653841"],
            "CRABP1": ["31653841"],
            "WIF1": ["31653841"],
            "FRZB": ["31653841"],
            "PTGDS": ["31653841"],
            "DKK3": ["31653841"],
            "RGR": ["31653841"],
            "GPM6B": ["31653841"],
            "SLC1A3": ["31653841"],
            "VIM": ["31653841"],
            "APOE": ["31653841"],
            "AQP4": ["31653841"],
        },
        "refine": {
            "CLU": {
                "note": "Clusterin; highest MG-specific signal (score > 25). "
                       "Complement regulator; AMD-associated.",
                "threshold": "score > 25.0",
                "pmid": "31653841",
            },
        },
    },

    "Microglia": {
        "confirm": {
            "C1QB": ["31653841"],
            "C1QA": ["31653841"],
            "C1QC": ["31653841"],
            "TYROBP": ["31653841"],
        },
        "add": {
            "FCER1G": ["31653841"],
            "LAPTM5": ["31653841"],
            "ITGB2": ["31653841"],
            "AIF1": ["31653841"],
            "P2RY12": ["31653841"],
            "CSF1R": ["31653841"],
            "CX3CR1": ["31653841"],
            "TREM2": ["31653841"],
        },
        "refine": {
            "C1QB": {
                "note": "Complement C1q B chain; highest microglia signal (score > 5.7). "
                       "Complement cascade is key AMD pathway.",
                "threshold": "score > 5.0",
                "pmid": "31653841",
            },
        },
    },

    "Astrocyte": {
        "confirm": {
            "GFAP": ["31653841"],
            "S100B": ["31653841"],
        },
        "add": {
            "ALDH1L1": ["31653841"],
            "AQP4": ["31653841"],
            "GJA1": ["31653841"],
            "SLC1A2": ["31653841"],
        },
    },

    "Vascular_Endothelial": {
        "confirm": {
            "CLDN5": ["31653841"],
            "VWF": ["31653841"],
        },
        "add": {
            "HIGD1B": ["31653841"],
            "GNG11": ["31653841"],
            "PECAM1": ["31653841"],
            "CDH5": ["31653841"],
            "ENG": ["31653841"],
            "FLT1": ["31653841"],
        },
        "refine": {
            "CLDN5": {
                "note": "Claudin-5; tight junction protein. Highest vascular "
                       "endothelial score (> 6.3). Also found in microglia "
                       "(shared expression).",
                "threshold": "score > 2.0",
                "pmid": "31653841",
            },
        },
    },
}

novel_types = []

expert_rules = [
    {
        "priority": 10,
        "condition": {
            "markers_present": {
                "RHO": 1.0,
                "NRL": 1.0,
                "PDE6A": 0.5,
            },
        },
        "action": "Rod_Photoreceptor",
        "source": "menon2019",
        "notes": "PDE6A is the strongest rod-specific marker by gene score (23.0)"
    },
    {
        "priority": 9,
        "condition": {
            "markers_present": {
                "RLBP1": 1.0,
                "GLUL": 1.0,
                "CLU": 1.0,
            },
        },
        "action": "Muller_Glia",
        "source": "menon2019",
        "notes": "CLU is the strongest MG-specific marker in adult human retina (score > 25)"
    },
    {
        "priority": 8,
        "condition": {
            "markers_present": {
                "ONECUT1": 1.0,
                "LHX1": 0.5,
            },
        },
        "action": "Horizontal_Cell",
        "source": "menon2019",
        "notes": "ONECUT1+LHX1 is the canonical HC signature; validated across platforms"
    },
]

conflicts = []
