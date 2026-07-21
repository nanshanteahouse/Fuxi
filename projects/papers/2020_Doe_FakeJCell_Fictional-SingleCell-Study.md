# SYNTHETIC TEST FIXTURE — NOT A REAL PAPER. ALL NAMES FABRICATED.

## Introduction

Single-cell transcriptomics has transformed our understanding of cellular heterogeneity in model organisms. The common laboratory newt (Notophthalmus syntheticus) is a well-established model for studying limb regeneration, but the transcriptional landscape of its regenerative blastema remains incompletely characterized. Understanding the cellular composition and gene expression programs of the regenerating newt limb is essential for elucidating the molecular mechanisms of tissue regeneration.

Figure 1 provides an overview of the experimental design and the major cell types identified in the regenerating newt limb blastema. We performed single-cell RNA sequencing on dissociated blastema tissue at three time points post-amputation to capture the dynamic cellular changes during regeneration.

## Results

Our single-cell RNA sequencing analysis of regenerating newt limb blastema generated a comprehensive atlas comprising over thirty thousand individual cells. After stringent quality control filtering, we retained cells with between five hundred and seven thousand five hundred detected genes and mitochondrial content below fifteen percent.

Figure 2 displays the UMAP embedding of all blastema cells colored by cluster identity. We identified fifteen transcriptionally distinct clusters corresponding to major cell types including mesenchymal progenitor cells, wound epidermis cells, macrophages, endothelial cells, Schwann cells, and muscle satellite cells. Notably, we observed a population of cells co-expressing markers of both mesenchymal and epithelial lineages, suggesting a transitional state during the regeneration process.

Figure 3 presents the expression patterns of key marker genes across all clusters. The mesenchymal progenitor population showed high expression of PRRX1, COL1A1, and FN1, while wound epidermis cells were characterized by KRT5, KRT14, and TP63 expression. Immune cell clusters expressed PTPRC, CD68, and CSF1R at high levels.

## Discussion

This comprehensive single-cell atlas of the regenerating newt limb blastema provides a valuable resource for the regeneration research community. The identification of a transitional mesenchymal-epithelial cell population suggests that lineage plasticity plays an important role in newt limb regeneration. Future studies integrating this atlas with functional perturbation experiments may reveal the gene regulatory networks that govern regenerative capacity in this remarkable model organism.

## Methods

Single-cell RNA sequencing libraries were prepared using the DropLogic Chromium platform according to the manufacturer's protocol. Newt limb blastema tissue was obtained from adult Notophthalmus syntheticus specimens maintained in the laboratory aquatic facility following institutional animal care committee approval. Blastema tissue was dissociated into single-cell suspensions using enzymatic digestion with collagenase and dispase followed by mechanical trituration. Cell viability was assessed by acridine orange and propidium iodide staining, and samples with viability greater than eighty-five percent were processed for library preparation.

Sequencing was performed on an AvantSeq X platform with paired-end reads of one hundred and fifty base pairs. Raw sequencing data were processed using the CellSieve pipeline version three point two with alignment to the custom-annotated newt reference genome Nsyn3.0. Downstream analysis was conducted using the Voxis single-cell analysis pipeline. Quality control filtering removed cells with fewer than five hundred detected genes, more than seven thousand five hundred genes, or mitochondrial gene content exceeding fifteen percent. Doublet detection was performed using DoubletGuard with default parameters.

Data normalization was performed using the pooling-based size factor method implemented in the scnorm package. Highly variable genes were selected using the variance-stabilizing transformation flavor with the top three thousand five hundred genes retained for downstream analysis. Principal component analysis was performed on the scaled data and the top twenty-five principal components were used for uniform manifold approximation and projection dimensionality reduction. Cell clustering was performed using the Louvain algorithm with resolution parameters ranging from zero point three to one point eight, and the optimal resolution was selected based on the Davies-Bouldin index and cluster stability metrics.

Differential expression analysis between clusters was performed using the negative binomial generalized linear model approach with empirical Bayes shrinkage for dispersion estimation. Cell-type annotation was performed using a combination of automated annotation with the Voxis annotation framework and manual curation based on known limb regeneration marker genes compiled from the literature. Marker gene validation was performed against the RepBase regeneration marker database.

Cell-cell communication analysis was performed using the InterCell framework to identify ligand-receptor interactions between blastema cell types. Gene regulatory network analysis was conducted using the Regulonator pipeline to identify transcription factor regulons active in each regenerating cell type. Pseudotime trajectory analysis was performed using the Wanderlust algorithm to reconstruct the differentiation continuum from progenitor cells to mature limb cell types. All statistical analyses were performed in the Hyperion scientific computing environment using the scanpy package version one point ten and associated libraries.
