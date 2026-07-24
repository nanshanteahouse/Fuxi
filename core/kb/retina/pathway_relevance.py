"""
core/kb/retina/pathway_relevance.py — Curated retina-specific pathway
metadata for tissue-aware enrichment analysis.

This module provides three constant data structures that enable the enrichment
module to (a) highlight retina-relevant pathways even when their p-values are
not the most significant, (b) filter out generic / housekeeping pathways that
carry no tissue-specific information, and (c) map pathway terms to known
retina marker genes for overlap-based relevance scoring.

References
----------
- Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023.
  *Genetics*, 224(1), iyad031. DOI: 10.1093/genetics/iyad031
- Kanehisa, M. et al. (2023). KEGG: Kyoto Encyclopedia of Genes and Genomes.
  *Nucleic Acids Research*, 51(D1), D587–D592. DOI: 10.1093/nar/gkac963
- Swaroop, A., Kim, D., & Forrest, D. (2010). Transcriptional regulation of
  photoreceptor development and homeostasis in the mammalian retina. *Nature
  Reviews Neuroscience*, 11, 563–576. DOI: 10.1038/nrn2880
- Hoon, M., Okawa, H., Della Santina, L., & Wong, R.O.L. (2014). Functional
  architecture of the retina: Development and disease. *Progress in Retinal
  and Eye Research*, 42, 44–84. DOI: 10.1016/j.preteyeres.2014.06.003

Notes
-----
- RETINA_KEY_PATHWAYS:  Enrichment results containing these terms are
  flagged for visual prominence regardless of rank order.
- RETINA_GENERIC_PATHWAYS: Terms on this list are de-emphasised or
  suppressed from top-N output because they reflect ubiquitous cellular
  machinery rather than tissue-specific biology.
- RETINA_KB_PATHWAY_MARKERS: {keyword_lowercase: [gene_symbols,...]}.
  During pathway-relevance scoring, each enriched term is matched
  case-insensitively against the dict keys; if the term *contains* a key,
  the corresponding markers are considered the known-relevant set for
  overlap computation against the enrichment's Overlap column.
"""

RETINA_KEY_PATHWAYS: list[str] = [
    # Phototransduction & visual perception
    "Phototransduction",
    "Phototransduction Cascade",
    "Visual Perception",
    "Sensory Perception of Light Stimulus",
    "Response to Light Stimulus",
    "Cellular Response to Light Stimulus",
    # Retina development & morphogenesis
    "Retina Development in Camera-Type Eye",
    "Eye Development",
    "Visual System Development",
    "Retinal Ganglion Cell Axon Guidance",
    "Axon Guidance",
    # Synaptic structure & function
    "Synaptic Signaling",
    "Synaptic Transmission",
    "Synapse Assembly",
    # Neuron biology
    "Neuron Projection Development",
    # Photoreceptor maintenance & outer segment
    "Photoreceptor Outer Segment",
    "Photoreceptor Cell Maintenance",
    "Retinoid Metabolic Process",
    "Cillium Assembly",
    "Cillium Organization",
]

RETINA_GENERIC_PATHWAYS: list[str] = [
    # Translation & transcription machinery
    "Ribosome",
    "Spliceosome",
    "Basal Transcription Factors",
    "mRNA Surveillance Pathway",
    "RNA Degradation",
    "RNA Transport",
    # Energy metabolism
    "Oxidative Phosphorylation",
    # Protein turnover
    "Proteasome",
    "Protein Export",
    # DNA maintenance
    "Nucleotide Excision Repair",
    "Mismatch Repair",
    "Homologous Recombination",
    "Base Excision Repair",
    "DNA Replication",
    "Fanconi Anemia Pathway",
    "Cell Cycle",
]

RETINA_KB_PATHWAY_MARKERS: dict[str, list[str]] = {
    "phototransduction": [
        "RHO",
        "OPN1SW",
        "OPN1MW",
        "OPN1LW",
        "GNAT1",
        "GNAT2",
        "PDE6A",
        "PDE6B",
        "PDE6G",
        "ARR3",
        "SAG",
        "RCVRN",
        "GRK1",
        "GUCY2D",
        "GUCA1A",
    ],
    "visual_perception": [
        "RHO",
        "OPN1SW",
        "OPN1MW",
        "OPN1LW",
        "CRX",
        "OTX2",
        "NRL",
        "NR2E3",
    ],
    "synapse": [
        "SLC17A6",
        "SLC17A7",
        "GAD1",
        "GAD2",
        "SLC32A1",
        "SYT1",
        "SYT2",
        "DLG4",
        "HOMER1",
    ],
    "retina_development": [
        "VSX2",
        "OTX2",
        "PAX6",
        "RAX",
        "SIX3",
        "SIX6",
        "LHX2",
        "SOX2",
        "VSX1",
        "PRDM1",
    ],
    "axon_guidance": [
        "NETO1",
        "NTN1",
        "DCC",
        "ROBO1",
        "ROBO2",
        "SLIT1",
        "EPHA5",
        "EPHB2",
        "SEMA3A",
    ],
    "photoreceptor_maintenance": [
        "PRPH2",
        "ROM1",
        "ABCA4",
        "CRB1",
        "RP1",
        "RP2",
        "RPGR",
        "USH2A",
    ],
    "cell_cycle": [
        "MKI67",
        "PCNA",
        "TOP2A",
        "CCND1",
        "CDK1",
        "CDK4",
    ],
}
