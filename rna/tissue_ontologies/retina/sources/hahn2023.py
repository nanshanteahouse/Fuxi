"""
hahn2023.py — Hahn et al. 2023 (Nature, PMID: 38092908)
Evolution of neuronal cell classes and types in the vertebrate retina.

Species: 17 vertebrate species across Mammalia, Aves, Reptilia, Teleostei,
         Cyclostomata (see source_meta).
Tissue: retina (multi-species comparison)
Platform: 10X Genomics (various)
Cells: >1M cells across all species
Key findings: All 6 major cell classes are conserved across vertebrates;
              midget RGCs have mammalian orthologs (previously thought
              primate-specific); evolutionary gradient from outer to inner
              retina in cell-type variation.
GEO: GSE237215 (SuperSeries, 13 sub-projects)

NOTE: This source covers MULTIPLE classes — "multi" is used in source_meta.
       Individual cell type entries carry per-type class/order annotations.
"""

source_meta = {
    "id": "hahn2023",
    "short_name": "Hahn 2023 Nature",
    "pmid": "38092908",
    "journal": "Nature",
    "year": 2023,
    "species": [
        "Homo sapiens",
        "Callithrix jacchus",
        "Bos taurus",
        "Ovis aries",
        "Sus scrofa",
        "Mustela putorius furo",
        "Didelphis marsupialis",
        "Tupaia belangeri",
        "Peromyscus maniculatus",
        "Rhabdomys pumilio",
        "Ictidomys tridecemlineatus",
        "Anolis sagrei",
        "Danio rerio",
    ],
    "tissue": "retina",
    "regions": ["whole retina"],
    "n_cells": 1000000,
    "n_subtypes": 60,
    "n_groups": 6,
    "class": "multi",
    "order": "multi",
}

markers = {
    # ── Rod Photoreceptor — conserved across ALL vertebrate classes ─
    "Rod_Photoreceptor": {
        "confirm": {
            "RHO": ["38092908"],
            "GNAT1": ["38092908"],
            "NR2E3": ["38092908"],
        },
        "add": {
            "PDE6A": ["38092908"],
            "PDE6B": ["38092908"],
            "SAG": ["38092908"],
            "RCVRN": ["38092908"],
            "CNGA1": ["38092908"],
            "CNGB1": ["38092908"],
            "RBP3": ["38092908"],
            "NRL": ["38092908"],
            "CRX": ["38092908"],
        },
        "refine": {
            "RHO": {
                "note": "Universally conserved across all vertebrate classes — "
                       "present in Mammalia, Aves, Reptilia, Teleostei, "
                       "Cyclostomata. The most deeply conserved retinal marker.",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
            "NRL": {
                "note": "Rod-fate specifier; conserved in all jawed vertebrates",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
        },
    },

    # ── Cone Photoreceptor — conserved across all vertebrate classes ─
    "Cone_Photoreceptor": {
        "confirm": {
            "ARR3": ["38092908"],
            "GNAT2": ["38092908"],
        },
        "add": {
            "OPN1SW": ["38092908"],
            "OPN1LW": ["38092908"],
            "OPN1MW": ["38092908"],
            "PDE6C": ["38092908"],
            "PDE6H": ["38092908"],
            "CNGB3": ["38092908"],
            "GNB3": ["38092908"],
            "THRB": ["38092908"],
            "RXRB": ["38092908"],
            "GUCY2E": ["38092908"],
            "GUCA1A": ["38092908"],
            "GUCA1B": ["38092908"],
        },
        "refine": {
            "ARR3": {
                "note": "Cone arrestin; conserved across all vertebrate classes",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
        },
    },

    # ── Bipolar Cell — conserved across jawed vertebrates ───────────
    "Bipolar_Cell": {
        "confirm": {
            "VSX2": ["38092908"],
            "GRM6": ["38092908"],
            "OTX2": ["38092908"],
        },
        "add": {
            "TRPM1": ["38092908"],
            "GRIK1": ["38092908"],
            "ISL1": ["38092908"],
            "PRKCA": ["38092908"],
            "CABP5": ["38092908"],
            "NYX": ["38092908"],
            "LRIT3": ["38092908"],
            "GNAO1": ["38092908"],
            "GNG13": ["38092908"],
            "PKC": ["38092908"],
            "SCGN": ["38092908"],
            "SLC4A5": ["38092908"],
        },
        "refine": {
            "VSX2": {
                "note": "Pan-bipolar marker (also expressed in Müller glia); "
                       "experimentally validated via CHX10 immunostaining "
                       "across mouse and human.",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
            "ISL1": {
                "note": "Marks ON-type bipolar cells; conserved across mammals",
                "threshold": "log2FC > 0.5",
                "pmid": "38092908",
            },
            "GRIK1": {
                "note": "Marks OFF-type bipolar cells; conserved across mammals",
                "threshold": "log2FC > 0.5",
                "pmid": "38092908",
            },
        },
    },

    # ── Amacrine Cell — conserved ───────────────────────────────────
    "Amacrine_Cell": {
        "confirm": {
            "TFAP2A": ["38092908"],
            "GAD1": ["38092908"],
        },
        "add": {
            "TFAP2B": ["38092908"],
            "TFAP2C": ["38092908"],
            "GAD2": ["38092908"],
            "SLC32A1": ["38092908"],
            "MEIS2": ["38092908"],
            "TCF4": ["38092908"],
            "GLRA1": ["38092908"],
            "GLRA2": ["38092908"],
            "SLC6A9": ["38092908"],
            "CALB2": ["38092908"],
            "PTH2": ["38092908"],
            "NPY": ["38092908"],
            "VIP": ["38092908"],
        },
        "refine": {
            "TFAP2A": {
                "note": "Pan-amacrine marker; conserved across all vertebrates; "
                       "validated via AP2A IHC in mouse and human.",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
            "GAD1": {
                "note": "GABAergic amacrine subclass marker",
                "threshold": "log2FC > 0.5",
                "pmid": "38092908",
            },
        },
    },

    # ── Horizontal Cell — conserved ─────────────────────────────────
    "Horizontal_Cell": {
        "confirm": {
            "ONECUT1": ["38092908"],
            "PROX1": ["38092908"],
            "CALB1": ["38092908"],
        },
        "add": {
            "ONECUT2": ["38092908"],
            "LHX1": ["38092908"],
            "EBF1": ["38092908"],
            "TSHZ2": ["38092908"],
            "ONECUT3": ["38092908"],
            "PRDM6": ["38092908"],
        },
    },

    # ── Retinal Ganglion Cell — conserved ──────────────────────────
    "RGC": {
        "confirm": {
            "RBPMS": ["38092908"],
            "POU4F2": ["38092908"],
            "NEFM": ["38092908"],
        },
        "add": {
            "POU4F1": ["38092908"],
            "NEFL": ["38092908"],
            "SLC17A6": ["38092908"],
            "THY1": ["38092908"],
            "SNCG": ["38092908"],
            "IRX3": ["38092908"],
            "BNC2": ["38092908"],
            "OPN4": ["38092908"],
            "EOMES": ["38092908"],
            "TBR1": ["38092908"],
            "NEUROD2": ["38092908"],
            "ISL1": ["38092908"],
            "POU4F3": ["38092908"],
            "STMN1": ["38092908"],
        },
        "refine": {
            "RBPMS": {
                "note": "Pan-RGC marker; validated across all vertebrate "
                       "classes via IHC. Most universally conserved RGC-"
                       "specific marker.",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
            "EOMES": {
                "note": "Marks a subset of RGCs including ipRGCs; conserved "
                       "across mammals but NOT found in teleost or lamprey",
                "threshold": "log2FC > 0.5",
                "pmid": "38092908",
            },
        },
    },

    # ── Müller Glia — conserved in jawed vertebrates ────────────────
    "Muller_Glia": {
        "confirm": {
            "RLBP1": ["38092908"],
            "GLUL": ["38092908"],
            "SLC1A3": ["38092908"],
        },
        "add": {
            "VIM": ["38092908"],
            "APOE": ["38092908"],
            "AQP4": ["38092908"],
            "CLU": ["38092908"],
            "DKK3": ["38092908"],
            "SOX2": ["38092908"],
            "NFIA": ["38092908"],
            "CRYM": ["38092908"],
            "GFAP": ["38092908"],
            "CST3": ["38092908"],
        },
        "refine": {
            "RLBP1": {
                "note": "Pan-MG marker conserved across jawed vertebrates; "
                       "not present in lamprey (Cyclostomata)",
                "threshold": "log2FC > 1.0",
                "pmid": "38092908",
            },
        },
    },

    # ── Non-neuronal ────────────────────────────────────────────────
    "Microglia": {
        "confirm": {
            "AIF1": ["38092908"],
            "P2RY12": ["38092908"],
            "CSF1R": ["38092908"],
        },
        "add": {
            "CX3CR1": ["38092908"],
            "ITGAM": ["38092908"],
            "TREM2": ["38092908"],
            "SALL1": ["38092908"],
            "C1QA": ["38092908"],
        },
    },

    "Astrocyte": {
        "confirm": {
            "GFAP": ["38092908"],
            "S100B": ["38092908"],
        },
        "add": {
            "ALDH1L1": ["38092908"],
            "AQP4": ["38092908"],
            "GJA1": ["38092908"],
            "SLC1A2": ["38092908"],
        },
    },

    "Vascular_Endothelial": {
        "confirm": {
            "PECAM1": ["38092908"],
            "CDH5": ["38092908"],
            "CLDN5": ["38092908"],
        },
        "add": {
            "FLT1": ["38092908"],
            "VWF": ["38092908"],
            "ENG": ["38092908"],
            "KDR": ["38092908"],
        },
    },

    # ── Lamprey-specific: Cyclostomata divergent RGC subtypes ──────
    "Lamprey_RGC": {
        "confirm": {
            "POU4F2": ["38092908"],
            "NEFM": ["38092908"],
        },
        "add": {
            "SLC17A6": ["38092908"],
            "THY1": ["38092908"],
            "STMN1": ["38092908"],
        },
        "refine": {
            "POU4F2": {
                "note": "RGC marker even in the most divergent extant vertebrate "
                       "(lamprey, Petromyzon marinus). RBPMS is NOT expressed "
                       "in lamprey RGCs — POU4F2 is the deepest pan-RGC marker.",
                "threshold": "log2FC > 0.5",
                "pmid": "38092908",
            },
        },
    },
}

novel_types = [
    {
        "name": "Lamprey_RGC",
        "parent": "RGC",
        "markers": ["POU4F2", "NEFM", "SLC17A6", "THY1"],
        "species": ["Petromyzon marinus"],
        "source": "hahn2023",
    },
]

expert_rules = [
    {
        "priority": 10,
        "condition": {
            "markers_present": {
                "RBPMS": 1.0,
                "POU4F2": 1.0,
                "SLC17A6": 1.0,
            },
        },
        "action": "RGC",
        "source": "hahn2023",
        "notes": "RBPMS+POU4F2 co-expression universally specifies RGC identity in jawed vertebrates",
    },
    {
        "priority": 9,
        "condition": {
            "markers_present": {
                "VSX2": 1.0,
                "GRM6": 0.5,
            },
            "markers_absent": ["RLBP1"],
        },
        "action": "Bipolar_Cell",
        "source": "hahn2023",
        "notes": "VSX2+GRM6 without RLBP1 distinguishes bipolar cells from Müller glia",
    },
    {
        "priority": 8,
        "condition": {
            "markers_present": {
                "POU4F2": 1.0,
                "NEFM": 1.0,
            },
            "markers_absent": ["RBPMS", "ISL1", "THY1"],
        },
        "action": "Lamprey_RGC",
        "source": "hahn2023",
        "notes": "POU4F2+NEFM without RBPMS, ISL1, or THY1: lamprey/cyclostome RGC"
                " pattern (RBPMS+ISL1+THY1 are gnathostome-specific)",
    },
]

conflicts = [
    {
        "type_a": {
            "cell_type": "RGC",
            "marker": "RBPMS",
        },
        "type_b": {
            "cell_type": "Lamprey_RGC",
            "marker": "RBPMS",
        },
        "notes": "RBPMS is NOT expressed in lamprey RGCs but is universal in "
                "jawed vertebrates. Annotators should use POU4F2+NEFM as the "
                "deepest RGC signature, and treat RBPMS absence as lamprey-"
                "specific rather than non-RGC.",
        "source": {"a": "hahn2023", "b": "hahn2023"},
    },
]
