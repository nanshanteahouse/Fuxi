"""core.ai — AI-powered annotation and interpretation.

Provides unified LLM API calling and prompt template management for
cell-type annotation across RNA, ATAC, and Spatial modalities.
"""

from core.ai.caller import ai_query
from core.ai.prompts import (
    ANNOTATION_SYSTEM_PROMPT,
    ATAC_ANNOTATION_SYSTEM_PROMPT,
    ATAC_ANNOTATION_USER_PROMPT_TEMPLATE,
    build_annotation_prompt,
    build_hierarchical_annotation_prompt,
)

__all__ = [
    "ai_query",
    "ANNOTATION_SYSTEM_PROMPT",
    "ATAC_ANNOTATION_SYSTEM_PROMPT",
    "ATAC_ANNOTATION_USER_PROMPT_TEMPLATE",
    "build_annotation_prompt",
    "build_hierarchical_annotation_prompt",
]
