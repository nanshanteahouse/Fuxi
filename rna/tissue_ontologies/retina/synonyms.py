"""
tissue_ontologies/retina/synonyms.py — Human-readable synonyms for 32 retina KB cell types.

This module provides comprehensive synonym dictionaries covering case variants,
abbreviations, full forms, singular/plural, punctuation variants, alternative
names from literature, and common AI output patterns for each canonical cell
type in the merged retina knowledge base.

Usage::

    from tissue_ontologies.retina.synonyms import RETINA_SYNONYMS

    display = RETINA_SYNONYMS["Muller_Glia"]["display_name"]   # "Muller Glia"
    aliases = RETINA_SYNONYMS["Muller_Glia"]["synonyms"]        # list[str]
"""

RETINA_SYNONYMS: dict[str, dict[str, str | list[str]]] = {
    # ═══════════════════════════════════════════════════════════════════════
    # Major types (no parent) — 5-8+ synonyms each
    # ═══════════════════════════════════════════════════════════════════════

    "ASCL1_Reprogrammed_MG": {
        "display_name": "ASCL1-Reprogrammed Muller Glia",
        "synonyms": [
            "ASCL1-Reprogrammed Muller Glia",
            "ASCL1 reprogrammed Muller glia",
            "ASCL1-reprogrammed MG",
            "ASCL1-reprogrammed Muller glial cell",
            "ASCL1_Reprogrammed_MG",
            "ASCL1-Reprogrammed Müller Glia",
            "ASCL1-reprogrammed Muller cell",
            "Reprogrammed Muller Glia",
            "ASCL1+ Muller Glia",
            "ASCL1-reprogrammed MG cells",
            "ascl1-reprogrammed muller glia",
        ],
    },

    "Amacrine_Cell": {
        "display_name": "Amacrine Cell",
        "synonyms": [
            "Amacrine Cell",
            "amacrine cell",
            "Amacrine",
            "amacrine",
            "amacrine cells",
            "Amacrine cells",
            "AmacrineCells",
            "AC",
            "amacrineCells",
            "Amacrine cell",
            "amacrine Cells",
        ],
    },

    "Amacrine_Precursor": {
        "display_name": "Amacrine Precursor",
        "synonyms": [
            "Amacrine Precursor",
            "Amacrine precursor",
            "amacrine precursor",
            "Amacrine Precursor Cells",
            "amacrine precursor cell",
            "Amacrine progenitor",
            "amacrine progenitor",
            "Amacrine Precursors",
            "Developing Amacrine",
            "developing amacrine",
            "amacrine precursor cells",
        ],
    },

    "Astrocyte": {
        "display_name": "Astrocyte",
        "synonyms": [
            "Astrocyte",
            "astrocyte",
            "Astrocytes",
            "astrocytes",
            "Astrocyte cell",
            "astrocyte cell",
            "Retinal astrocyte",
            "retinal astrocyte",
            "AstrocyteCells",
            "astrocyte cells",
            "astrocytic cell",
        ],
    },

    "Bipolar_Cell": {
        "display_name": "Bipolar Cell",
        "synonyms": [
            "Bipolar Cell",
            "bipolar cell",
            "Bipolar",
            "bipolar",
            "bipolar cells",
            "Bipolar cells",
            "BipolarCells",
            "Retinal bipolar cell",
            "retinal bipolar cell",
            "BC",
            "bipolarCells",
            "bipolar Cells",
            "Retinal Bipolar Cell",
        ],
    },

    "Bipolar_Precursor": {
        "display_name": "Bipolar Precursor",
        "synonyms": [
            "Bipolar Precursor",
            "bipolar precursor",
            "Bipolar precursor",
            "bipolar precursor cell",
            "Bipolar Precursor Cells",
            "Bipolar progenitor",
            "bipolar progenitor",
            "Developing bipolar cell",
            "developing bipolar cell",
            "Bipolar Precursors",
            "bipolar precursor cells",
        ],
    },

    "Cone_Photoreceptor": {
        "display_name": "Cone Photoreceptor",
        "synonyms": [
            "Cone Photoreceptor",
            "cone photoreceptor",
            "Cone photoreceptor",
            "cone",
            "Cone",
            "cone cells",
            "Cone cells",
            "cone photoreceptors",
            "Cone Photoreceptors",
            "CONE",
            "ConePhotoreceptor",
            "Retinal cone",
            "retinal cone",
            "retinal cone cell",
            "cone cell",
            "Cone cell",
        ],
    },

    "Developing_AC_HC_Precursors": {
        "display_name": "Developing Amacrine/Horizontal Cell Precursors",
        "synonyms": [
            "Developing Amacrine/Horizontal Cell Precursors",
            "Developing AC/HC precursors",
            "Developing amacrine horizontal cell precursors",
            "developing AC HC precursors",
            "AC/HC precursors",
            "Amacrine/Horizontal precursors",
            "Developing AC and HC precursors",
            "developing amacrine/horizontal precursors",
            "Developing_AC_HC_Precursors",
            "developing AC/HC precursor cells",
        ],
    },

    "Developing_BC_Photo_Precursors": {
        "display_name": "Developing Bipolar/Photoreceptor Precursors",
        "synonyms": [
            "Developing Bipolar/Photoreceptor Precursors",
            "Developing BC/Photo precursors",
            "Developing bipolar photoreceptor precursors",
            "developing BC photo precursors",
            "BC/Photoreceptor precursors",
            "Bipolar/Photoreceptor precursors",
            "Developing bipolar and photoreceptor precursors",
            "Developing_BC_Photo_Precursors",
            "developing bipolar/photo precursors",
            "BC/Photo precursor cells",
        ],
    },

    "Fetal_RPE": {
        "display_name": "Fetal RPE",
        "synonyms": [
            "Fetal RPE",
            "fetal RPE",
            "Fetal retinal pigment epithelium",
            "fetal retinal pigment epithelium",
            "FetalRPE",
            "Fetal RPE cells",
            "fetal RPE cells",
            "Developing RPE",
            "developing RPE",
            "fetal retinal pigment epithelial cells",
            "Fetal RPE cell",
        ],
    },

    "Fibroblast": {
        "display_name": "Fibroblast",
        "synonyms": [
            "Fibroblast",
            "fibroblast",
            "Fibroblasts",
            "fibroblasts",
            "Fibroblast cell",
            "fibroblast cell",
            "FibroblastCells",
            "Retinal fibroblast",
            "retinal fibroblast",
            "fibroblast cells",
            "Retinal Fibroblast",
        ],
    },

    "Foveal_Cones": {
        "display_name": "Foveal Cones",
        "synonyms": [
            "Foveal Cones",
            "foveal cones",
            "Foveal cone",
            "foveal cone",
            "Foveal Cone",
            "Foveal cone cells",
            "foveal cone cells",
            "FovealCones",
            "Foveal photoreceptors",
            "foveal photoreceptors",
            "Foveal cone photoreceptors",
            "Foveal_Cone",
        ],
    },

    "Foveal_ML_Cones": {
        "display_name": "Foveal M/L Cones",
        "synonyms": [
            "Foveal M/L Cones",
            "foveal M/L cones",
            "Foveal ML cones",
            "foveal ML cones",
            "Foveal medium/long wavelength cones",
            "Foveal M cone",
            "Foveal L cone",
            "Foveal_ML_Cones",
            "Foveal ML Cone",
            "foveal ML cone cells",
            "foveal medium wavelength cones",
            "foveal long wavelength cones",
        ],
    },

    "Foveal_S_Cones": {
        "display_name": "Foveal S Cones",
        "synonyms": [
            "Foveal S Cones",
            "foveal S cones",
            "Foveal short wavelength cones",
            "foveal short wavelength cones",
            "Foveal blue cones",
            "foveal blue cones",
            "Foveal_S_Cones",
            "Foveal S Cone",
            "foveal S cone cells",
            "Foveal short-wave cones",
        ],
    },

    "Horizontal_Cell": {
        "display_name": "Horizontal Cell",
        "synonyms": [
            "Horizontal Cell",
            "horizontal cell",
            "Horizontal",
            "horizontal",
            "horizontal cells",
            "Horizontal cells",
            "HorizontalCells",
            "HC",
            "Retinal horizontal cell",
            "retinal horizontal cell",
            "horizontalCells",
            "Horizontal Cells",
        ],
    },

    "Macula_Specific_Cone": {
        "display_name": "Macula-Specific Cone",
        "synonyms": [
            "Macula-Specific Cone",
            "Macula specific cone",
            "macula-specific cone",
            "Macula-specific cone",
            "Macular cone",
            "macular cone",
            "Macula cone",
            "Macula_Specific_Cone",
            "Macula-specific cone cell",
            "macula specific cone photoreceptor",
            "Macular cone cell",
        ],
    },

    "Marmoset_Foveal_Cones": {
        "display_name": "Marmoset Foveal Cones",
        "synonyms": [
            "Marmoset Foveal Cones",
            "marmoset foveal cones",
            "Marmoset foveal cone",
            "marmoset foveal cone",
            "Marmoset cone",
            "marmoset cone",
            "Marmoset_Foveal_Cones",
            "Marmoset foveal cone cells",
            "marmoset foveal cone photoreceptors",
            "Marmoset Foveal Cone",
        ],
    },

    "Microglia": {
        "display_name": "Microglia",
        "synonyms": [
            "Microglia",
            "microglia",
            "Microglial cell",
            "microglial cell",
            "Microglial cells",
            "microglial cells",
            "Retinal microglia",
            "retinal microglia",
            "MicrogliaCells",
            "microglial",
            "Microglial",
            "retinal microglial cell",
        ],
    },

    "Muller_Glia": {
        "display_name": "Muller Glia",
        "synonyms": [
            "Muller Glia",
            "Muller glia",
            "muller glia",
            "Muller cell",
            "Muller glial cell",
            "Müller glia",
            "Müller cell",
            "Müller glial cell",
            "MG",
            "muller glial cells",
            "Müller glial cells",
            "Muller glial cells",
            "MullerGlia",
            "Muller cells",
            "muller cell",
        ],
    },

    "Neonatal_RPCs": {
        "display_name": "Neonatal RPCs",
        "synonyms": [
            "Neonatal RPCs",
            "neonatal RPCs",
            "Neonatal RPC",
            "Neonatal retinal progenitor cells",
            "neonatal retinal progenitor cells",
            "Neonatal retinal progenitors",
            "Neonatal_RPCs",
            "Neonatal progenitor cells",
            "neonatal retinal progenitor",
            "neonatal RPC",
            "Neonatal Retinal Progenitor Cells",
        ],
    },

    "Oligodendrocyte": {
        "display_name": "Oligodendrocyte",
        "synonyms": [
            "Oligodendrocyte",
            "oligodendrocyte",
            "Oligodendrocytes",
            "oligodendrocytes",
            "Oligodendrocyte cell",
            "oligodendrocyte cell",
            "Oligodendroglia",
            "oligodendroglia",
            "Retinal oligodendrocyte",
            "retinal oligodendrocyte",
            "oligodendrocyte cells",
            "OligodendrocyteCells",
        ],
    },

    "Pericyte": {
        "display_name": "Pericyte",
        "synonyms": [
            "Pericyte",
            "pericyte",
            "Pericytes",
            "pericytes",
            "Pericyte cell",
            "pericyte cell",
            "PericyteCells",
            "Retinal pericyte",
            "retinal pericyte",
            "pericyte cells",
            "Pericyte Cells",
            "retinal pericytes",
        ],
    },

    "Peripheral_Amacrine": {
        "display_name": "Peripheral Amacrine",
        "synonyms": [
            "Peripheral Amacrine",
            "peripheral amacrine",
            "Peripheral amacrine cell",
            "peripheral amacrine cell",
            "Peripheral amacrine cells",
            "peripheral amacrine cells",
            "Peripheral_Amacrine",
            "Peripheral AC",
            "peripheral AC",
            "Peripheral amacrine Cells",
            "peripheral amacrine Cells",
        ],
    },

    "Peripheral_Rods": {
        "display_name": "Peripheral Rods",
        "synonyms": [
            "Peripheral Rods",
            "peripheral rods",
            "Peripheral rod",
            "peripheral rod",
            "Peripheral rod cells",
            "Peripheral rod photoreceptors",
            "peripheral rod photoreceptors",
            "Peripheral_Rods",
            "peripheral rod cell",
            "Peripheral Rod Cells",
            "Peripheral Rod Photoreceptor",
        ],
    },

    "Photoreceptor_Precursor": {
        "display_name": "Photoreceptor Precursor",
        "synonyms": [
            "Photoreceptor Precursor",
            "photoreceptor precursor",
            "Photoreceptor precursor",
            "Photoreceptor precursor cell",
            "photoreceptor precursor cell",
            "Photoreceptor progenitor",
            "photoreceptor progenitor",
            "Photoreceptor Precursors",
            "Photo precursor",
            "photo precursor",
            "Photoreceptor Precursor Cells",
            "photoreceptor precursor cells",
        ],
    },

    "Proliferating_MG": {
        "display_name": "Proliferating Muller Glia",
        "synonyms": [
            "Proliferating Muller Glia",
            "proliferating Muller glia",
            "Proliferating MG",
            "proliferating MG",
            "Proliferating Muller cell",
            "proliferating Muller cell",
            "Proliferating Müller glia",
            "Proliferating_MG",
            "Proliferating Muller glial cells",
            "proliferating Muller glial cell",
            "Proliferating Müller cells",
            "proliferating muller glia",
        ],
    },

    "Proliferating_RPC": {
        "display_name": "Proliferating RPC",
        "synonyms": [
            "Proliferating RPC",
            "proliferating RPC",
            "Proliferating retinal progenitor",
            "proliferating retinal progenitor",
            "Proliferating progenitor cells",
            "Proliferating RPCs",
            "Proliferating_RPC",
            "Proliferating retinal progenitor cells",
            "proliferating retinal progenitor cells",
            "Proliferating RPC cell",
        ],
    },

    "RGC": {
        "display_name": "Retinal Ganglion Cell",
        "synonyms": [
            "RGC",
            "Retinal ganglion cell",
            "retinal ganglion cell",
            "Retinal Ganglion Cell",
            "ganglion cell",
            "RGCs",
            "retinal ganglion cells",
            "Ganglion cells",
            "ganglion cells",
            "RGC cell",
            "RetinalGanglionCell",
            "Retinal Ganglion Cells",
            "rgc",
        ],
    },

    "RPC": {
        "display_name": "Retinal Progenitor Cell",
        "synonyms": [
            "RPC",
            "retinal progenitor cell",
            "Retinal progenitor cell",
            "Retinal Progenitor Cell",
            "retinal progenitor",
            "Retinal progenitor",
            "retinal progenitors",
            "RPCs",
            "RetinalProgenitorCell",
            "Retinal Progenitor Cells",
            "retinal progenitor cells",
            "rpc",
        ],
    },

    "RPE": {
        "display_name": "Retinal Pigment Epithelium",
        "synonyms": [
            "RPE",
            "retinal pigment epithelium",
            "Retinal Pigment Epithelium",
            "RPE cells",
            "retinal pigment epithelial cells",
            "pigment epithelium",
            "RPE cell",
            "Retinal pigment epithelium (RPE)",
            "RPE Cells",
            "retinal pigment epithelial cell",
            "rpe",
            "Pigment Epithelium",
        ],
    },

    "Rod_Photoreceptor": {
        "display_name": "Rod Photoreceptor",
        "synonyms": [
            "Rod Photoreceptor",
            "rod photoreceptor",
            "Rod photoreceptor",
            "rod",
            "Rod",
            "rod cells",
            "Rod cells",
            "rod photoreceptors",
            "Rod Photoreceptors",
            "ROD",
            "RodPhotoreceptor",
            "Retinal rod",
            "retinal rod",
            "rod cell",
            "Rod cell",
            "retinal rod cell",
        ],
    },

    "Vascular_Endothelial": {
        "display_name": "Vascular Endothelial Cell",
        "synonyms": [
            "Vascular Endothelial Cell",
            "vascular endothelial cell",
            "Vascular Endothelial",
            "vascular endothelial",
            "Endothelial",
            "endothelial",
            "Endothelial cell",
            "endothelial cell",
            "VascularEndothelial",
            "Vascular Endothelial Cells",
            "vascular endothelial cells",
            "Endothelial Cells",
        ],
    },

    # ═══════════════════════════════════════════════════════════════════════
    # Developmental / Progenitor types (v3.0.0+, zuo2024 + hahn2023)
    # ═══════════════════════════════════════════════════════════════════════

    "PRPC": {
        "display_name": "Proliferating Retinal Progenitor Cell",
        "synonyms": [
            "Proliferating Retinal Progenitor Cell",
            "PRPC",
            "proliferating RPC",
            "Proliferating RPC",
            "Retinal Progenitor Cell (proliferating)",
            "early RPC",
            "Early RPC",
            "cycling RPC",
            "Cycling RPC",
            "Retinal Progenitor",
            "retinal progenitor cell",
            "proliferating retinal progenitor",
            "RPC-proliferating",
            "PRPCs",
        ],
    },

    "NRPC_RGC_fate": {
        "display_name": "Neurogenic RPC (RGC Fate)",
        "synonyms": [
            "Neurogenic RPC (RGC Fate)",
            "NRPC_RGC_fate",
            "RGC-fated neurogenic RPC",
            "RGC-fated NRPC",
            "RGC progenitor",
            "RGC Precursor",
            "ATOH7+ RPC",
            "RGC-committed RPC",
            "neurogenic retinal progenitor (RGC)",
            "RGC-fated progenitor",
            "developing RGC",
            "RGC neuroblast",
        ],
    },

    "NRPC_AC_HC_fate": {
        "display_name": "Neurogenic RPC (AC/HC Fate)",
        "synonyms": [
            "Neurogenic RPC (AC/HC Fate)",
            "NRPC_AC_HC_fate",
            "AC/HC-fated NRPC",
            "AC/HC progenitor",
            "Amacrine/Horizontal progenitor",
            "amacrine-horizontal precursor",
            "PRDM13+ RPC",
            "AC/HC-committed RPC",
            "neurogenic retinal progenitor (AC/HC)",
            "AC-HC precursor",
        ],
    },

    "NRPC_Cone_BC_fate": {
        "display_name": "Neurogenic RPC (Cone/BC Fate)",
        "synonyms": [
            "Neurogenic RPC (Cone/BC Fate)",
            "NRPC_Cone_BC_fate",
            "Cone/BC-fated NRPC",
            "Cone-Bipolar progenitor",
            "OTX2+ RPC",
            "photoreceptor-bipolar precursor",
            "Cone/BC progenitor",
            "Cone/BC-committed RPC",
            "neurogenic retinal progenitor (Cone/BC)",
            "Cone-BC precursor",
        ],
    },

    "NRPC_Rod_fate": {
        "display_name": "Neurogenic RPC (Rod Fate)",
        "synonyms": [
            "Neurogenic RPC (Rod Fate)",
            "NRPC_Rod_fate",
            "Rod-fated NRPC",
            "Rod progenitor",
            "Rod Precursor",
            "NRL+ RPC",
            "CRX+ NRL+ RPC",
            "Rod-committed RPC",
            "neurogenic retinal progenitor (Rod)",
            "rod photoreceptor precursor",
            "developing rod",
        ],
    },

    # ═══════════════════════════════════════════════════════════════════════
    # Cross-species naming variants (v3.0.0+, hahn2023)
    # ═══════════════════════════════════════════════════════════════════════

    "Lamprey_RGC": {
        "display_name": "Lamprey Retinal Ganglion Cell",
        "synonyms": [
            "Lamprey Retinal Ganglion Cell",
            "Lamprey RGC",
            "lamprey RGC",
            "Cyclostome RGC",
            "cyclostome RGC",
            "Petromyzon RGC",
            "sea lamprey RGC",
            "lamprey retinal ganglion",
            "Agnathan RGC",
            "agnathan RGC",
        ],
    },
}
