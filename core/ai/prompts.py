#!/usr/bin/env python3
"""
ai_prompts.py — AI 注释与解读的提示词模板（RNA + ATAC 统一）
==============================================================

集中存放所有 LLM 提示词模板和构建函数，确保提示词一致、可复用、易维护。
支持 scRNA-seq 和 scATAC-seq 两种组学类型。

设计原则:
  - 提示词与调用逻辑分离（提示词在 core/prompts/*.yaml，逻辑在此模块）
  - RNA 和 ATAC 的 prompt 以前缀区分: ANNOTATION_SYSTEM_PROMPT (RNA), ATAC_ANNOTATION_SYSTEM_PROMPT
  - build_annotation_prompt() 自动运行 marker 基因检测并组装完整提示词

用法 (RNA):
    from core.ai_prompts import ANNOTATION_SYSTEM_PROMPT, build_annotation_prompt
    sys_prompt, user_prompt = build_annotation_prompt(adata, "retina", "human")

用法 (ATAC):
    from core.ai_prompts import ATAC_ANNOTATION_SYSTEM_PROMPT, ATAC_ANNOTATION_USER_PROMPT_TEMPLATE

提示词源文件:
    core/prompts/annotation.yaml            → scRNA-seq annotation
    core/prompts/atac_annotation.yaml       → scATAC-seq annotation
    core/prompts/atac_interpretation.yaml   → scATAC-seq interpretation
    core/prompts/paper_meta.yaml            → paper metadata extraction
    core/prompts/paper_figure.yaml          → figure legend extraction
    core/prompts/paper_methods.yaml         → methods section extraction
    core/prompts/paper_methodology.yaml     → methodology pattern analysis
"""

import json

from core.ai.templates._loader import load_prompt

# ═══════════════════════════════════════════════════════════════════════
#  scRNA-seq — 聚类注释提示词
# ═══════════════════════════════════════════════════════════════════════

_annotation = load_prompt("annotation")
ANNOTATION_SYSTEM_PROMPT: str = _annotation["system"]
ANNOTATION_USER_PROMPT_TEMPLATE: str = _annotation["user_template"]


def build_annotation_prompt(
    adata,
    tissue: str,
    species: str,
    precomputed_rank: bool = False,
    extra_context: str = "",
    compact: bool = False,
    kb_candidates: list[str] | None = None,
    unconstrained: bool = False,
):
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
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            method="wilcoxon",
            use_raw=True if adata.raw is not None else None,
        )

    # ── 提取每聚类 marker 基因 ───────────────────────────────────────
    n_top = 5 if compact else 20
    clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
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

_atac_ann = load_prompt("atac_annotation")
ATAC_ANNOTATION_SYSTEM_PROMPT: str = _atac_ann["system"]
ATAC_ANNOTATION_USER_PROMPT_TEMPLATE: str = _atac_ann["user_template"]

_atac_interp = load_prompt("atac_interpretation")
ATAC_INTERPRETATION_SYSTEM_PROMPT: str = _atac_interp["system"]
ATAC_INTERPRETATION_USER_PROMPT_TEMPLATE: str = _atac_interp["user_template"]


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
#  Paper Interpretation Prompts (loaded from YAML)
# ═══════════════════════════════════════════════════════════════════════

_paper_meta = load_prompt("paper_meta")
PAPER_META_SYSTEM_PROMPT: str = _paper_meta["system"]
PAPER_META_USER_TEMPLATE: str = _paper_meta["user_template"]

_paper_figure = load_prompt("paper_figure")
PAPER_FIGURE_SYSTEM_PROMPT: str = _paper_figure["system"]
PAPER_FIGURE_USER_TEMPLATE: str = _paper_figure["user_template"]

_paper_methods = load_prompt("paper_methods")
PAPER_METHODS_SYSTEM_PROMPT: str = _paper_methods["system"]
PAPER_METHODS_USER_TEMPLATE: str = _paper_methods["user_template"]

_paper_methodology = load_prompt("paper_methodology")
PAPER_METHODOLOGY_SYSTEM_PROMPT: str = _paper_methodology["system"]
PAPER_METHODOLOGY_USER_TEMPLATE: str = _paper_methodology["user_template"]
