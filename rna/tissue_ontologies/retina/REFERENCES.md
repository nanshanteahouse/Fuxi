# Retina Knowledge Base — References

Bibliography of all source papers whose data forms the retina knowledge base.

## Summary

| Metric | Value |
|--------|-------|
| Total sources | 16 |
| Unique markers | ~2100 |
| Cell types | ~70 |
| Species | Homo sapiens (10), Mus musculus (3), Macaca fascicularis (1), Danio rerio (2), Gallus gallus (1), plus hahn2023's 13 vertebrate classes |
| Journals | Nature Communications (5), Cell (3), Neuron (2), Nature Genetics (1), PLoS Biology (1), PLoS Genetics (1), Nature (1), eLife (1), Cell Genomics (1), Cell Reports (1), preprint (1) |
| Sources with supplement verification | **16/16** |
| Sources with empirical validation | 2 (hu2019, menon2019) |

---

## Sources

### Hu et al. (2019)
- **PMID**: 31269016
- **Journal**: PLoS Biology
- **Title**: Single-cell transcriptomic analysis of the developing human retina
- **Species**: Homo sapiens (fetal retina, fetal RPE)
- **Contribution**: Foundational fetal retina atlas defining 9 cell types. Provides RPE65+BEST1 co-expression rule for fetal RPE and VSX2+SOX2 rule for RPCs.
- **Data source**: s012.xlsx cluster markers (22,555 rows) + GSE107618 empirical validation
- **Validation**: Empirically validated against GSE107618 and GSE118614
- **Audit**: Confirm-tier errors fixed + add-tier genes not in HVG set removed (2026-07-10)

### Menon et al. (2019)
- **PMID**: 31653841
- **Journal**: Nature Communications
- **Title**: Single-cell transcriptomic atlas of the adult human retina
- **Species**: Homo sapiens (macula, periphery)
- **Contribution**: Adult human retina atlas with 20,000 cells across 12 subtypes. CLU identified as strongest MG-specific marker (score > 25).
- **Data source**: MOESM5 gene scores + MOESM10 Wilcox table + GSE118614 empirical validation
- **Validation**: Empirically cross-validated with GSE118614
- **Audit**: POU4F2/RBPMS confirm errors fixed (2026-07-10)

### Peng et al. (2019)
- **PMID**: 30712875
- **Journal**: Cell
- **Title**: Single-cell transcriptomes of the macaque retina reveal cell-type-specific expression patterns
- **Species**: Macaca fascicularis (fovea, periphery)
- **Contribution**: First primate retina cell atlas with 165,679 cells across 60 subtypes.
- **Data source**: Supplemental Table 9 (supp-9)
- **Audit**: 8 null audit dates + 5 empty supplement_verified fields fixed (2026-07-16)

### Shekhar et al. (2016)
- **PMID**: 27565351
- **Journal**: Cell
- **Title**: Comprehensive classification of retinal bipolar neurons by single-cell transcriptomics
- **Species**: Mus musculus (whole retina, P14)
- **Contribution**: First comprehensive mouse bipolar cell atlas. 14 BC subtypes from 26 Louvain-Jaccard clusters. ~228 add-tier markers.
- **Data source**: supplement-10 KnownMarkers + supplement-11 Hi Conf DE genes
- **Audit**: Verified (2026-07-16)

### Macosko et al. (2015) — NEW
- **PMID**: 26000488
- **Journal**: Cell
- **Title**: Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets
- **Species**: Mus musculus (whole retina, P14)
- **Contribution**: Seminal Drop-seq mouse retina atlas with 44,808 cells across 39 clusters. First comprehensive mouse retina single-cell reference covering all major types.
- **Data source**: supp_data_4.xlsx FINAL_MARKERS_FOR_EACH_CLUSTER (4296 markers with AUC/diff/power)
- **KB entries**: 9 cell types, 414 markers

### Chen et al. (2023) / Liang et al. (2023)
- **PMID**: 37388908
- **Journal**: Cell Genomics
- **Title**: Multi-omic human retina cell atlas reveals cell-type-specific chromatin architecture
- **Species**: Homo sapiens (fovea, macula, periphery)
- **Contribution**: 10 cell types. Multi-omics (snRNA-seq + snATAC-seq) human retina atlas.
- **Data source**: mmc2.xlsx major_markers_RNA
- **Audit**: PRKG1/WIF1 intra-source conflicts documented (2026-07-10)

### Hahn et al. (2023)
- **PMID**: 38092908
- **Journal**: Nature
- **Title**: Evolution of the vertebrate retinal cell types
- **Species**: 13 species across Mammalia, Reptilia, Teleostei, Cyclostomata (1M cells)
- **Contribution**: Pan-vertebrate cross-species retina atlas. 7 cell types with cross-species conservation scores.
- **Data source**: MOESM8 class markers + MOESM10 orthotype gene markers

### Zuo et al. (2024)
- **PMID**: 39117640
- **Journal**: Nature Communications
- **Title**: Integrated single-cell transcriptomic analysis reveals human retinal development
- **Species**: Homo sapiens (fovea, periphery, developing)
- **Contribution**: Developmental retina atlas with 220,000 cells across 22 subtypes. Proliferating_RPC and NRPC trajectories.
- **Data source**: Supplementary Table MOESM6
- **Audit**: Microglia/Astrocyte/RPE markers not in MOESM6 documented (2026-07-16)

### Tran et al. (2019)
- **PMID**: 31784286
- **Journal**: Neuron
- **Title**: Single-cell profiles of retinal ganglion cells differing in resilience to injury reveal neuroprotective genes
- **Species**: Mus musculus (whole retina)
- **Contribution**: Mouse RGC subtype atlas with 35,699 cells across 46 RGC types in 7 subclasses.
- **Data source**: Table S2
- **KB entries**: RGC (pan-RGC) + 7 RGC subtypes
- **Audit**: Pan-RGC markers reclassified as domain_knowledge (2026-07-16)

### Yamagata et al. (2021)
- **PMID**: 33393903
- **Journal**: eLife
- **Title**: A cell atlas of the chick retina based on single-cell transcriptomics
- **Species**: Gallus gallus (E12, E16, E18)
- **Contribution**: First chicken retina atlas with 33,000 cells across 136 cell types. Defines Chicken_Double_Cone.
- **Data source**: Paper text and figures (no marker tables in supplements)
- **KB entries**: 9 cell types

### Hoang et al. (2023)
- **PMID**: preprint (HRCA)
- **Journal**: preprint
- **Title**: ASCL1 induces neurogenesis in human Muller glia
- **Species**: Homo sapiens (fetal retina, MG culture)
- **Contribution**: MG reprogramming study. ASCL1 is exogenously overexpressed by retrovirus — NOT endogenous.
- **Data source**: GSE246169
- **Note**: Preprint — enters KB with explicit audit annotations

### Dorgau et al. (2024)
- **PMID**: 38670973
- **Journal**: Nature Communications
- **Title**: Single-cell transcriptomics of the human ciliary margin reveals retinal progenitor cell diversity
- **Species**: Homo sapiens (whole retina, ciliary margin)
- **Contribution**: CMZ and RPC diversity atlas. Neurogenic transition states T1/T2/T3. GN2BL1 and PCPB4 are published symbols (standard HGNC: GNB2L1, PCP4).
- **Data source**: MOESM7 gene signatures + MOESM13 per-type marker lists
- **Audit**: Full rewrite with complete audit blocks (2026-07-16)

### Li et al. (2026)
- **PMID**: 41578023
- **Journal**: Nature Genetics
- **Title**: Large-scale single-cell atlas of the human retina
- **Species**: Homo sapiens (fovea, macula, periphery)
- **Contribution**: Largest human retina atlas to date with 3.4M cells. 6 cell types use canonical markers from add_canonical_markers().
- **Data source**: Supplementary_Tables.xlsx (S8/S9)
- **Audit**: Audit sections added to all 9 cell types (2026-07-16)

### Wang et al. (2026)
- **PMID**: 41528844
- **Journal**: Cell Reports
- **Title**: Single-cell multiome and enhancer connectome of the human retina and choroid
- **Species**: Homo sapiens (RPE/choroid)
- **Contribution**: RPE/choroid atlas. 9 cell types with 95 markers. 100% supplement-verified.
- **Data source**: supplement-6 SupData1 scRNA-seq markers
- **Audit**: Full xlsx cross-verification passed (2026-07-16)

### Kölsch et al. (2021) — NEW
- **PMID**: 33357413
- **Journal**: Neuron
- **Title**: Molecular classification of zebrafish retinal ganglion cells
- **Species**: Danio rerio (adult retina)
- **Contribution**: First zebrafish RGC subtype atlas. 32 RGC subclusters in adult retina with DE gene lists.
- **Data source**: supp-3.xlsx Adult DE tables (9667 rows with avg_logFC)
- **KB entries**: 1 cell type (RGC, 28 markers)

### Liu et al. (2022) — NEW
- **PMID**: 35245286
- **Journal**: PLoS Genetics
- **Title**: Rod genesis driven by mafba in an nrl knockout zebrafish model
- **Species**: Danio rerio (P60 retina)
- **Contribution**: Zebrafish nrl knockout model. Pre-computed marker genes for 25+ retinal clusters covering all major types.
- **Data source**: pgen.1009841.s016.xlsx S4 Dataset (marker genes for each cluster)
- **KB entries**: 10 cell types, 224 markers

---

## Notes

### Cross-species coverage

- **Human/Macaque**: Hu 2019, Menon 2019, Peng 2019, Chen 2023, Zuo 2024, Dorgau 2024, Li 2026, Wang 2026, Hoang 2023
- **Mouse**: Macosko 2015, Shekhar 2016, Tran 2019
- **Zebrafish**: Kölsch 2021, Liu 2022
- **Chicken**: Yamagata 2021
- **Pan-vertebrate**: Hahn 2023 (13 species from lamprey to human)

### Audit methodology

The 2026-07-16 cross-batch audit verified marker data across all 16 sources. Key actions:
- dorgau2024: complete rewrite with audit blocks, gene symbol documentation
- li2026: audit sections added to all 9 cell types
- peng2019/tran2019/yamagata2021/hoang2023: structural fixes (dates, notes, supp verifications)
- wang2026: 100% marker match verified against 11 xlsx supplement files
- merge.py: None-value handling for empty YAML fields
- 3 new sources added: macosko2015, kolsch2021, liu2022

---

*Generated from source_meta in 16 YAML source files. Last updated: 2026-07-16.*
