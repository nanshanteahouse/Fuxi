#!/usr/bin/env python3
"""
ai_prompts.py — AI 注释与解读的提示词模板（RNA + ATAC 统一）
==============================================================

集中存放所有 LLM 提示词模板和构建函数，确保提示词一致、可复用、易维护。
支持 scRNA-seq 和 scATAC-seq 两种组学类型。

设计原则:
  - 提示词与调用逻辑分离（提示词在此模块，调用在 core.ai_caller）
  - RNA 和 ATAC 的 prompt 以前缀区分: ANNOTATION_SYSTEM_PROMPT (RNA), ATAC_ANNOTATION_SYSTEM_PROMPT
  - build_annotation_prompt() 自动运行 marker 基因检测并组装完整提示词

用法 (RNA):
    from core.ai_prompts import ANNOTATION_SYSTEM_PROMPT, build_annotation_prompt
    sys_prompt, user_prompt = build_annotation_prompt(adata, "retina", "human")

用法 (ATAC):
    from core.ai_prompts import ATAC_ANNOTATION_SYSTEM_PROMPT, ATAC_ANNOTATION_USER_PROMPT_TEMPLATE
"""

import json
import os


# ═══════════════════════════════════════════════════════════════════════
#  scRNA-seq — 聚类注释提示词
# ═══════════════════════════════════════════════════════════════════════

ANNOTATION_SYSTEM_PROMPT = """You are an expert single-cell RNA-seq biologist with deep knowledge of cell type identification across tissues and species.

For each cluster ID provided in the user message, return a JSON object mapping cluster IDs to annotations with the following fields:
  - cell_type  : the broad cell type (e.g. "T cell", "Macrophage", "Oligodendrocyte", "Excitatory neuron")
  - state      : activation or functional state (e.g. "resting", "activated", "proliferating", "N/A")
  - subtype    : the most specific subtype (e.g. "CD8+ cytotoxic T cell", "M1 macrophage", "SST+ interneuron", "N/A")
  - confidence : one of "high", "medium", or "low" — based on how specific and well-established the markers are
  - reasoning  : a single sentence citing the key marker genes that support your annotation

Return ONLY a valid JSON object. No explanation, no markdown formatting, no code fences.
Include ALL cluster IDs in the response.

Required format:
{"0":{"cell_type":"T cell","state":"activated","subtype":"CD8+ cytotoxic T cell","confidence":"high","reasoning":"High CD8A, GZMB, PRF1 expression indicates cytotoxic T cells"},"1":{"cell_type":"...","state":"...","subtype":"...","confidence":"...","reasoning":"..."}}

IMPORTANT — Cross-species guidance:
- Gene names have been mapped to human orthologs where possible (prefix "UNMAPPED_" indicates no mapping).
- For retina data, the following non-neuronal cell types may be present and should be considered: Microglia (AIF1, CSF1R, CX3CR1, CD74, P2RY12), Pericytes (PDGFRB, RGS5, CSPG4), Astrocytes (AQP4, SLC1A2, GFAP, ALDH1L1), Retinal Pigment Epithelium (RPE65, RDH5, BEST1), Endothelial cells (PECAM1, CDH5, VWF, CLDN5), Oligodendrocytes (MBP, PLP1, MOG), Vascular Smooth Muscle (ACTA2, MYH11, TAGLN).
- Do NOT assume non-neuronal types are absent — check their markers carefully before classifying a cluster as "Neuron" or "Retinal neuron".
- For clusters with mostly UNMAPPED_ genes, set confidence to "low" and note the limited gene annotation."""


ANNOTATION_USER_PROMPT_TEMPLATE = """Tissue: {tissue}
Species: {species}

Marker genes per cluster (top {n_top} by Wilcoxon score):
{cluster_markers_json}

Return ONLY a valid JSON object mapping each cluster ID to its annotation. Include ALL cluster IDs."""


def build_annotation_prompt(adata, tissue: str, species: str,
                            precomputed_rank: bool = False,
                            extra_context: str = "",
                            compact: bool = False,
                            kb_candidates: list[str] | None = None,
                            unconstrained: bool = False):
    """
    构建 RNA 聚类注释的完整提示词对。

    可自动运行或跳过 Wilcoxon rank-sum 检验。当调用者已经执行过
    rank_genes_groups 时，传入 precomputed_rank=True 避免重复计算，
    直接使用 adata.uns['rank_genes_groups'] 中的已有结果。

    参数:
        adata:  已聚类（leiden 列）的 AnnData 对象
        tissue: 组织名称（如 "hypothalamus", "retina"）
        species: 物种名称（如 "human", "mouse"）
        precomputed_rank: 若为 True，跳过 rank_genes_groups 计算
        extra_context: 额外上下文信息追加到用户提示词尾部
        compact: 若为 True，每聚类仅展示 top 5 而非 top 20 marker 基因
        kb_candidates: 若提供，限制 AI 只能从该列表中选取细胞类型名称
        unconstrained: 若为 True，kb_candidates 作为参考而非约束；
            AI 可以建议列表外的细胞类型，用 ``[NOVEL] `` 前缀标记
            （v3.1.0+ 审计模式 / 新组织类型探测）

    返回:
        (system_prompt, user_prompt) 二元组，可直接传入 ai_query()
    """
    import scanpy as sc

    # ── 计算 marker 基因（如尚未计算）────────────────────────────────
    if not precomputed_rank:
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                            use_raw=True if adata.raw is not None else None)

    # ── 提取每聚类 marker 基因 ───────────────────────────────────────
    n_top = 5 if compact else 20
    clusters = sorted(adata.obs["leiden"].unique(),
                      key=lambda x: int(x))
    cluster_markers: dict = {}
    for cl in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=str(cl))
        top_genes = df.head(n_top)["names"].tolist()
        cluster_markers[cl] = top_genes

    # ── 组装提示词 ────────────────────────────────────────────────────
    user_prompt = ANNOTATION_USER_PROMPT_TEMPLATE.format(
        tissue=tissue,
        species=species,
        n_top=n_top,
        cluster_markers_json=json.dumps(cluster_markers, indent=2),
    )
    if extra_context:
        user_prompt += f"\n\n{extra_context}"

    if kb_candidates:
        candidates_str = "\n".join(f"  - {c}" for c in kb_candidates)
        if unconstrained:
            reference_text = f"""

REFERENCE — Known cell types in the knowledge base:
{candidates_str}

These are known cell types but you are NOT limited to this list.
If a cluster's markers do NOT match any known type, suggest a cell
type outside this list and prefix it with '[NOVEL] ' so it can be
flagged for KB review.  Use this sparingly — only when the markers
clearly indicate a cell type not represented in the reference list.
"""
            user_prompt += reference_text
        else:
            constraint_text = f"""

IMPORTANT — Constrained naming:
You MUST choose cell type names from this list ONLY:
{candidates_str}

Rules:
- Select the type that best matches the cluster's marker genes
- If NO type from the list fits well, output "Uncertain" and explain why
- Do NOT create new type names outside this list
"""
            user_prompt += constraint_text

    return ANNOTATION_SYSTEM_PROMPT, user_prompt


def _filter_candidates_by_category(
    kb_candidates: list[str],
    kb_hierarchy: dict,
    target_category: str | None = None,
) -> list[str]:
    """Filter KB candidates to those belonging to a specific broad category.

    Args:
        kb_candidates: Full list of cell type names from KB
        kb_hierarchy: The _hierarchy dict from KB (with "categories" key)
        target_category: Broad category name (e.g., "Progenitor", "Neuron", ...)
            If None, returns the full (unfiltered) list.

    Returns:
        Filtered list of cell types belonging to target_category.
    """
    if not target_category or not kb_hierarchy:
        return kb_candidates

    categories = kb_hierarchy.get("categories", {})
    if target_category not in categories:
        return kb_candidates  # fallback: return unfiltered

    member_set = set(categories[target_category].get("members", []))
    # Also include the broad type itself
    member_set.add(f"Broad_{target_category}")

    return [c for c in kb_candidates if c in member_set]


def build_hierarchical_annotation_prompt(
    adata,
    tissue: str,
    species: str,
    kb_candidates: list[str] | None = None,
    kb_hierarchy: dict | None = None,
    precomputed_rank: bool = False,
    extra_context: str = "",
    compact: bool = False,
    unconstrained: bool = False,
):
    """Build hierarchical annotation prompt with broad-category constraint.

    Like build_annotation_prompt(), but adds explicit instruction for the
    AI to first determine the broad category (Progenitor / Neuron / Glia /
    Non-neural) and then select a subtype only from within that category.

    Args are identical to build_annotation_prompt() with one addition:
        kb_hierarchy: The _hierarchy section from the KB (contains "categories")

    Returns:
        (system_prompt, user_prompt) tuple
    """
    # Get the base system+user prompts from the existing function
    system_prompt, user_prompt = build_annotation_prompt(
        adata=adata,
        tissue=tissue,
        species=species,
        precomputed_rank=precomputed_rank,
        extra_context=extra_context,
        compact=compact,
        kb_candidates=kb_candidates,
        unconstrained=unconstrained,
    )

    # Add hierarchical instruction to user_prompt
    if kb_hierarchy and kb_candidates:
        cat_names = list(kb_hierarchy.get("categories", {}).keys())
        cat_str = " | ".join(cat_names)
        hierarchy_text = f"""
HIERARCHICAL ANNOTATION — Two-step classification required:
1. FIRST: Determine the broad category for each cluster. Choose from:
   [{cat_str}]

2. SECOND: Within that broad category, select the specific cell type from
   the KB candidate list. A cluster's cell type MUST belong to its broad
   category — cross-category assignments will be flagged as errors.

IMPORTANT: Include "cell_category" AND "cell_type" in your JSON output.
The cell_category field must be one of: {cat_str}.
The cell_type field must be a subtype within that category from the candidate list.
"""
        user_prompt += hierarchy_text

    return system_prompt, user_prompt


# ═══════════════════════════════════════════════════════════════════════
#  scATAC-seq — 染色质状态注释提示词
# ═══════════════════════════════════════════════════════════════════════

ATAC_ANNOTATION_SYSTEM_PROMPT = """You are an expert in single-cell epigenomics and chromatin biology.
Your task is to annotate ATAC-seq clusters based on their marker peaks and genomic features.
Output ONLY valid JSON with cluster IDs as keys and objects containing 'cell_type', 'confidence', and 'reasoning'."""

ATAC_ANNOTATION_USER_PROMPT_TEMPLATE = """Annotate these ATAC-seq clusters from {tissue}:

Cluster summary (JSON):
{cluster_summary}

For each cluster, determine the likely chromatin state or cell type based on:
- Top marker peak regions (chromatin accessibility near specific gene loci)
- Number of cells in the cluster

Return ONLY JSON in format:
{{
    "0": {{"cell_type": "...", "confidence": "high|medium|low", "reasoning": "..."}},
    "1": ...
}}

Possible chromatin states / cell types include (but not limited to):
- Active Progenitors (high accessibility at cell cycle / proliferation genes)
- Primed Neuronal (accessible at neuronal TF loci)
- Photoreceptor lineage (cone/rod)
- Glial lineage (Müller glia, astrocytes)
- Retinal Ganglion Cells
- Excitatory/Inhibitory Neurons
- Interneurons (Amacrine, Horizontal, Bipolar)
- Vascular/Endothelial
- Microglia"""

ATAC_INTERPRETATION_SYSTEM_PROMPT = """You are an expert epigenomics analyst.
Interpret the ATAC-seq analysis results, focusing on biological insights from chromatin accessibility data."""

ATAC_INTERPRETATION_USER_PROMPT_TEMPLATE = """Interpret these ATAC-seq results for {tissue}:

Marker peaks results: {marker_results}
Motif enrichment: {motif_results}
Enrichment analysis: {enrichment_results}

Provide a concise biological interpretation. Focus on:
1. Key cell types/states identified by their chromatin accessibility profiles
2. Important transcription factors (from motif analysis) and their potential roles
3. Notable pathways from enrichment analysis
4. Biological insights from the data
5. Limitations and caveats"""


# ═══════════════════════════════════════════════════════════════════════
#  RNA 后续步骤的提示词桩（TODO: 在对应步骤实现时完善）
# ═══════════════════════════════════════════════════════════════════════

# PARAM_SUGGEST_PROMPT = """..."""
# 用途: 根据数据特征建议 QC 参数阈值

# QC_REVIEW_PROMPT = """..."""
# 用途: 审查 QC 结果并给出质量判断

# DEG_DESIGN_PROMPT = """..."""
# 用途: 建议差异表达分析的对比设计

# INTERPRETATION_PROMPT = """..."""
# 用途: 解读差异表达或富集分析结果


# ═══════════════════════════════════════════════════════════════════════
#  Paper Interpretation Prompts (added 2026-07)
# ═══════════════════════════════════════════════════════════════════════

PAPER_META_SYSTEM_PROMPT = """You are an expert biomedical research analyst who extracts structured metadata from paper abstracts.

Given the abstract of a single-cell genomics paper, extract the following structured information in JSON format:

1. paper_type — the category of the paper: one of "research", "resource", "review", "method"
2. experimental_design — describes the biological and technical setup of the study:
   - species: the model organism (use standard NCBI taxonomy names like homo_sapiens, mus_musculus, macaca_mulatta, danio_rerio, drosophila_melanogaster, etc.)
   - tissue: the tissue/organ studied (e.g. retina, brain, pancreas, liver)
   - tissue_info: a brief description of the tissue context/subregion
   - models: a list of experimental models used, each with a name and brief description
   - conditions: a list of experimental conditions compared, each with name and description
   - modalities: sequencing modalities used (e.g. scRNA-seq, scATAC-seq, multiome, CITE-seq)
   - summary: a 2-3 sentence overview of the experimental design
3. key_findings — an array of 1-6 key biological findings reported in the abstract
4. data_access — accession identifiers for public data repositories:
   - geo_ids: an array of GEO accession IDs (e.g. ["GSE137537"])
   - sra_ids: an array of SRA accession IDs (e.g. ["SRP123456"])
5. data_notes — an array of important data characteristics (e.g. number of cells, sequencing depth, sample origin)

EXAMPLE:

Input: "We performed single-nucleus RNA-seq on postmortem human retinal samples from 6 donors using the 10x Genomics platform. Our analysis identified 58 transcriptionally distinct cell types, including novel subtypes of amacrine cells. The data are available through GEO accession GSE137537."

Expected output:
{
  "paper_type": "research",
  "experimental_design": {
    "species": "homo_sapiens",
    "tissue": "retina",
    "tissue_info": "postmortem human retina from 6 donors",
    "models": [{"name": "postmortem human retina", "description": "Six normal donors"}],
    "conditions": [{"name": "Normal", "description": "Control retinas"}],
    "modalities": ["snRNA-seq"],
    "summary": "Single-nucleus RNA-seq on postmortem human retinal samples from 6 donors, identifying 58 cell types."
  },
  "key_findings": ["Identified 58 transcriptionally distinct cell types", "Novel subtypes of amacrine cells discovered"],
  "data_access": {"geo_ids": ["GSE137537"], "sra_ids": []},
  "data_notes": ["snRNA-seq — use is_nuclei=True", "6 postmortem donors"]
}

The input text may contain minor PDF-conversion artifacts such as word concatenation or watermark text — ignore these and extract only the structured data from the coherent portions.

Return ONLY a valid JSON object. No explanation, no markdown formatting, no code fences.

IMPORTANT: Species must use the homo_sapiens / mus_musculus naming convention.
If the species is unclear, set it to "unknown" and note this in data_notes."""

PAPER_META_USER_TEMPLATE = """Extract experimental design, paper type, key findings, data access, and data notes from this paper abstract:

Abstract:
{abstract_text}

Return ONLY the JSON object as specified."""

PAPER_FIGURE_SYSTEM_PROMPT = """You are an expert in single-cell bioinformatics figure interpretation.

Analyze the given figure legend and extract structured information. Use the following controlled vocabulary for figure types:

  umap, tsne, pca, heatmap, dotplot, violin, barplot, feature_plot, trajectory,
  enrichment, volcano, scatter, lineplot, genome_browser, motif_analysis,
  immunofluorescence, imaging_3d, schematic, electrophysiology, other,
  stacked_violin, ridgeplot, correlation_heatmap, spatial_feature, pseudotime, cellchat

REPRODUCIBILITY: Based on the figure type and description, determine if this figure can be reproduced from raw single-cell data. Provide your reasoning in the reproducibility_reasoning field.

Return a JSON object with the following fields:
  id: the figure identifier (e.g. 'Fig_2b')
  caption: the raw figure legend text verbatim
  type: one of the controlled vocabulary above
  panels: array of panel labels (e.g. ['2B', '2C'])
  parameters.features: array of gene/feature names shown
  parameters.resolution: clustering resolution if mentioned, else null
  parameters.method: computational method if mentioned, else null
  parameters.conditions: array of experimental conditions compared
  parameters.n_value: sample size if mentioned (e.g. "n=6 donors"), else null
  parameters.error_bar_type: type of error bars if mentioned (e.g. "SD", "SEM", "95% CI"), else null
  purpose: one-sentence summary of what this figure shows
  reproducible: boolean indicating whether the figure can be reproduced from raw data
  reproducibility_reasoning: explanation of why the reproducibility decision was made

EXAMPLE:

Input: "Fig. 1: Single-cell transcriptomic analysis of human retina. a, UMAP embedding of 20,091 cells colored by cluster identity (n=6 donors). b, Dotplot showing expression of canonical retinal marker genes across clusters. c, Barplot of cell type proportions per donor."

Expected output:
{
  "id": "Fig_1",
  "caption": "Fig. 1: Single-cell transcriptomic analysis of human retina. a, UMAP embedding of 20,091 cells colored by cluster identity (n=6 donors). b, Dotplot showing expression of canonical retinal marker genes across clusters. c, Barplot of cell type proportions per donor.",
  "type": "umap",
  "panels": ["1a", "1b", "1c"],
  "parameters": {
    "features": null,
    "resolution": null,
    "method": null,
    "conditions": [],
    "n_value": "n=6 donors",
    "error_bar_type": null
  },
  "purpose": "Overview of single-cell transcriptomic analysis of human retina showing clustering, marker gene expression, and cell type proportions.",
  "reproducible": true,
  "reproducibility_reasoning": "UMAP, dotplot, and barplot are computational visualizations generated from scRNA-seq data, reproducible with the same count matrix and parameters."
}

The input text may contain minor PDF-conversion artifacts such as word concatenation or watermark text — ignore these and extract only the structured data from the coherent portions.

Return ONLY a valid JSON object. No explanation, no markdown formatting, no code fences."""

PAPER_FIGURE_USER_TEMPLATE = """Extract structured information from this figure legend:

Figure text:
{figure_text}

Return ONLY the JSON object with figure type, parameters, purpose, caption, and reproducibility status."""

PAPER_METHODS_SYSTEM_PROMPT = """You are an expert in single-cell bioinformatics methods extraction.

Given a paper's Methods section (or relevant portions), extract structured details about the computational and experimental methods used.

Identify:
  1. key_methods — array of specific method/platform names (e.g. '10x Genomics Chromium Single Cell 3\' v3', 'Seurat v4.0', 'Cell Ranger 7.0')
  2. software_versions — object mapping software names to version strings
  3. reference_genome — genome assembly name (e.g. hg38, GRCh38, mm10, GRCm39)
  4. sequencing_platforms — array of sequencing instrument names (e.g. "Illumina NovaSeq 6000", "Illumina NextSeq 2000")
  5. data_notes — array of important notes about data processing (e.g. 'snRNA-seq — use is_nuclei=True', 'data were filtered to remove doublets')

EXAMPLE:

Input: "Single-cell libraries were prepared using the 10x Genomics Chromium Single Cell 3' v3 reagent kit and sequenced on an Illumina NovaSeq 6000 platform. Reads were aligned to the hg38 reference genome using CellRanger v7.0. Downstream analysis was performed with Seurat v4.0."

Expected output:
{
  "key_methods": ["10x Genomics Chromium Single Cell 3' v3", "CellRanger v7.0", "Seurat v4.0"],
  "software_versions": {"CellRanger": "7.0", "Seurat": "4.0"},
  "reference_genome": "hg38",
  "sequencing_platforms": ["Illumina NovaSeq 6000"],
  "data_notes": ["10x Genomics Chromium Single Cell 3' v3 used", "Sequenced on NovaSeq 6000"]
}

The input text may contain minor PDF-conversion artifacts such as word concatenation or watermark text — ignore these and extract only the structured data from the coherent portions.

Return ONLY a valid JSON object with this structure:
{
  "key_methods": ["method 1", "method 2"],
  "software_versions": {"SoftwareName": "version"},
  "reference_genome": "assembly_name",
  "sequencing_platforms": ["platform 1", "platform 2"],
  "data_notes": ["note 1", "note 2"]
}

No explanation, no markdown formatting, no code fences."""

PAPER_METHODS_USER_TEMPLATE = """Extract bioinformatics methods, software versions, reference genome, sequencing platforms, and data processing notes from this text:

Methods text:
{methods_text}

Return ONLY the JSON object."""
