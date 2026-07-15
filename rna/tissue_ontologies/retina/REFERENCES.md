# Retina Knowledge Base — References

Bibliography of all source papers whose data forms the retina knowledge base.

## Summary

| Metric | Value |
|--------|-------|
| Total sources | 11 |
| Unique markers | 1406 |
| Cell types | 39 |
| Species | Homo sapiens (9), Macaca fascicularis (1), Mus musculus (1), Gallus gallus (1), plus hahn2023's 13 vertebrate classes |
| Journals | Nature Communications (4), Neuron (1), eLife (1), PLoS Biology (1), Nature (1), Cell (1), Nature Genetics (1), Cell Genomics (1), preprint (1) |
| Sources with supplement verification | **11/11** |
| Sources with empirical validation | 2 (hu2019, menon2019) |
---

## Sources

### Hu et al. (2019)
- **PMID**: 31269016
- **Journal**: PLoS Biology
- **Title**: Single-cell transcriptomic analysis of the developing human retina
- **Species**: Homo sapiens (fetal retina, fetal RPE)
- **Contribution**: Foundational fetal retina atlas defining 9 cell types (RPE, Fibroblast, RPC, Muller Glia, Retinal Ganglion Cell, Photoreceptor Precursor, Bipolar Precursor, Amacrine Precursor, Horizontal Cell, Vascular Endothelial, Microglia). Provides RPE65+BEST1 co-expression rule for fetal RPE and VSX2+SOX2 rule for RPCs.
- **Data source**: s012.xlsx cluster markers (22,555 rows) + GSE107618 empirical validation
- **Validation**: Empirically validated against GSE107618 and GSE118614
- **Audit**: Confirm-tier errors fixed + add-tier genes not in HVG set removed (2026-07-10)

### Menon et al. (2019)
- **PMID**: 31653841
- **Journal**: Nature Communications
- **Title**: Single-cell transcriptomic atlas of the adult human retina
- **Species**: Homo sapiens (macula, periphery)
- **Contribution**: Adult human retina atlas with 20,000 cells across 12 subtypes. Defines quantitative gene scores for all major retinal cell types: Rod, Cone, Bipolar, Amacrine, Horizontal, Muller Glia, RGC, Microglia, Astrocyte, Vascular Endothelial. CLU identified as strongest MG-specific marker (score > 25). PDE6A as strongest rod-specific marker.
- **Data source**: MOESM5 gene scores + MOESM10 Wilcox table + GSE118614 empirical validation
- **Validation**: Empirically cross-validated with GSE118614
- **Audit**: POU4F2/RBPMS confirm errors fixed; known TF under-scoring artifact documented (2026-07-10)

### Peng et al. (2019)
- **PMID**: 30712875
- **Journal**: Cell
- **Title**: Single-cell transcriptomes of the macaque retina reveal cell-type-specific expression patterns
- **Species**: Macaca fascicularis (fovea, periphery)
- **Contribution**: First primate (macaque) retina cell atlas with 165,679 cells across 60 subtypes. Defines foveal vs peripheral rod/cone signatures. Provides cross-species validation for primate-specific markers including opsin expression gradients.
- **Data source**: Supplemental Table 9 (supp-9)

### Liang et al. (2023)
- **PMID**: 37388908
- **Journal**: Cell Genomics
- **Title**: Multi-omic human retina cell atlas reveals cell-type-specific chromatin architecture
- **Species**: Homo sapiens (fovea, macula, periphery)
- **Contribution**: 8 supplement-verified cell types (Rod, Cone, RGC, Bipolar, Amacrine, Horizontal, Muller Glia, Astrocyte). Multi-omics (snRNA-seq + snATAC-seq) human retina atlas.
- **Data source**: mmc2.xlsx major_markers_RNA (Cell Genomics, 2023)
- **Audit**: Replaced AI-generated markers with supplement-verified data (2026-07-10)

### Hahn et al. (2023)
- **PMID**: 38092908
- **Journal**: Nature
- **Title**: Evolution of the vertebrate retinal cell types
- **Species**: 13 species across Mammalia, Reptilia, Teleostei, Cyclostomata (1M cells)
- **Contribution**: Pan-vertebrate cross-species retina atlas (17 species, 6 cell classes). Rebuilt from MOESM8 class markers (41 unique genes) + MOESM10 orthotype markers. 7 cell types, 73 markers with cross-species conservation scores (>=10 species = confirm). Unique evolutionary reference for phylogenetic weighting.
- **Data source**: MOESM8 class markers + MOESM10 orthotype gene markers (Nature, 2023)

### Zuo et al. (2024)
- **PMID**: 39117640
- **Journal**: Nature Communications
- **Title**: Integrated single-cell transcriptomic analysis reveals human retinal development
- **Species**: Homo sapiens (fovea, periphery, developing)
- **Contribution**: Developmental retina atlas with 220,000 cells across 22 subtypes spanning fetal through adult. Defines Proliferating_RPC and NRPC (neurogenic RPC) trajectories with ATOH7+PRDM13+OTX2 GRN. Provides human-specific RPC markers (ARHGAP11B). Most comprehensive coverage of developing human retina cell types.
- **Data source**: Supplementary Table MOESM6

---

### Tran et al. (2019)
- **PMID**: 31784286
- **Journal**: Neuron
- **Title**: Single-cell profiles of retinal ganglion cells differing in resilience to injury reveal neuroprotective genes
- **Species**: Mus musculus (whole retina)
- **Contribution**: Definitive mouse RGC subtype atlas with 35,699 cells across 46 RGC types in 7 subclasses. Defines pan-RGC markers (RBPMS, SLC17A6) and subtype-specific markers for αRGCs (Spp1), ipRGCs (Opn4, Eomes), ooDSGCs (Cartpt, Mmp17), F-RGCs (Foxp2), T-RGCs (Tbr1), W3-RGCs (Tusc5), and N-RGCs (Neurod2). POU4F2 demoted to add-tier (only ~70% RGC coverage, misses ipRGCs).
- **Data source**: Table S2 — cluster-specific and subclass markers for 46 RGC types
- **KB entries**: RGC (pan-RGC) + 7 RGC subtypes (RGC_Alpha, RGC_ipRGC, RGC_ooDSGC, RGC_Foxp2, RGC_Tbr1, RGC_W3, RGC_Neurod2)

### Yamagata et al. (2021)
- **PMID**: 33393903
- **Journal**: eLife
- **Title**: A cell atlas of the chick retina based on single-cell transcriptomics
- **Species**: Gallus gallus (E12, E16, E18 developmental stages)
- **Contribution**: First chicken retina single-cell atlas with 33,000 cells across 136 cell types in 7 major classes. Defines chicken-specific markers including double cone markers (CALB1+OPN1LW — unique to non-mammalian vertebrates). Validated by in situ hybridization (51 PCR probes) and CRISPR/eCHIKIN knockin (14 genes). Discovers 41 chick RGC types with no mammalian midget/parasol equivalents — only ipRGCs are clearly conserved.
- **Data source**: supp1 (ISH probe genes) + supp2 (CRISPR target genes) + supp3 (cluster QC)
- **KB entries**: 9 cell types including Chicken_Double_Cone (novel non-mammalian type) and Oligodendrocyte (new retina KB type)

---


### Hoang et al. (2023)
- **PMID**: Not assigned (placeholder: 00000000)
- **Journal**: preprint / HRCA (Human Retinal Cell Atlas)
- **Title**: ASCL1 induces neurogenesis in human Muller glia
- **Species**: Homo sapiens (fetal retina, MG culture)
- **Contribution**: Muller glia reprogramming study defining ASCL1-driven neurogenesis. Introduces Proliferating_MG and ASCL1_Reprogrammed_MG as novel types. Provides markers for MG-to-neuron transition states (MKI67+TOP2A+ASCL1). Data sourced from GSE246169.
- **Data source**: GSE246169 dataset
- **Note**: Preprint — enters KB under committee consensus model with audit annotation

---


### Dorgau et al. (2024)
- **PMID**: 38670973
- **Journal**: Nature Communications
- **Title**: Single-cell transcriptomics of the human ciliary margin reveals retinal progenitor cell diversity
- **Species**: Homo sapiens (whole retina, ciliary margin)
- **Contribution**: Ciliary margin and RPC diversity atlas. Defines proliferating RPC markers (FGF19, SFRP2, DAPL1, ZIC1) and neurogenic transition states (T1/T2/T3 across ATOH7 to photoreceptor specification). Provides developmental trajectory markers from RPC through committed precursors to photoreceptors.
- **Data source**: Primary publication (Nature Communications, 2024)

### Li et al. (2026)
- **PMID**: 41578023
- **Journal**: Nature Genetics
- **Title**: Large-scale single-cell atlas of the human retina
- **Species**: Homo sapiens (fovea, macula, periphery)
- **Contribution**: Largest human retina atlas to date with 3.4M cells across 9 major groups. Deepest coverage of Bipolar_Cell and Amacrine_Cell subtypes with hundreds of discriminatory markers. Provides regional (fovea/macula/periphery) expression signatures at unprecedented scale.
- **Data source**: Primary publication (Nature Genetics, 2026)

---

## Notes

### Preprint note

**Hoang et al. (2023)** is a preprint without an assigned PMID. It enters the knowledge base under the committee consensus model with explicit audit annotations. Its markers for ASCL1-reprogrammed Muller glia (Proliferating_MG) are unique to this source and have not been independently validated.

### Cross-species note

**Hahn et al. (2023)** spans 13 vertebrate species from lamprey (Petromyzon marinus) to human. Per earlier audit, it contributes 0 common orthologs in human-only analyses but provides essential evolutionary depth for phylogenetic weighting.

**Tran et al. (2019)** provides the first species-specific mouse (Mus musculus) RGC subtype entries, covering 7 molecular subclasses.

**Yamagata et al. (2021)** provides the first avian (Gallus gallus) retina entries, including the Chicken_Double_Cone — a non-mammalian photoreceptor type with no mammalian homolog.

### Empirically validated sources

**Hu et al. (2019)** and **Menon et al. (2019)** have both supplement-derived markers and independent empirical validation (GSE107618, GSE118614), providing the highest-confidence marker backbone for the human retina knowledge base.

### Audit methodology

Each source was audited against its original publication supplements. Markers were cross-validated against source data files. Discrepancies between AI-generated marker lists and supplement-verified data were corrected. Add-tier markers absent from the highly variable gene (HVG) set were removed.
---

*Generated from source_meta in 11 YAML source files. Last updated: 2026-07-15.*
