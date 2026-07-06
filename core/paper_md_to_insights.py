#!/usr/bin/env python3
"""
paper_md_to_insights.py — AI-assisted extraction of structured paper insights

Reads a paper's PDF-to-text markdown file and uses LLM (via core.ai_caller)
to generate a structured insights.yaml for reproduction tracking.

Usage:
    python core/paper_md_to_insights.py <paper.md> [--output OUTPUT] [--force]
"""

import re
import json
import os
import sys
import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# Ensure repo root is on sys.path for standalone CLI usage
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from core.ai_caller import ai_query
from core.ai_prompts import (
    PAPER_META_SYSTEM_PROMPT,
    PAPER_META_USER_TEMPLATE,
    PAPER_FIGURE_SYSTEM_PROMPT,
    PAPER_FIGURE_USER_TEMPLATE,
    PAPER_METHODS_SYSTEM_PROMPT,
    PAPER_METHODS_USER_TEMPLATE,
)

logger = logging.getLogger(__name__)

# ── Section / figure regex patterns ────────────────────────────────────────────
_SECTION_RE = re.compile(r'^(?:#+\s*)?(SUMMARY|Abstract|Introduction|Results|Discussion|Methods|Experimental\s*Procedures)\b', re.MULTILINE | re.IGNORECASE)
_FIGURE_RE = re.compile(r'(?:Figure|Fig\.?)\s+\d+[a-z]?', re.IGNORECASE)


@dataclass
class _LLMConfig:
    """Minimal LLM config for standalone CLI use."""
    model: str = "deepseek/deepseek-v4-flash"
    api_base: str = "https://api.deepseek.com/v1"
    api_key: str = ""
    max_tokens: int = 16384
    temperature: float = 0.3
    ai_cache_responses: bool = True
    thinking_enabled: bool = False


def _parse_filename_meta(md_path: str) -> dict:
    """Parse year/author/journal from filename pattern {year}_{author}_{journal}_{title}."""
    stem = Path(md_path).stem
    parts = stem.split('_')
    meta: dict = {
        "filename": Path(md_path).name,
        "stem": stem,
        "year": None,
        "first_author": None,
        "journal": None,
        "title": stem,
    }
    if len(parts) >= 3 and parts[0].isdigit() and len(parts[0]) in (2, 4):
        meta["year"] = parts[0]
        meta["first_author"] = parts[1]
        meta["journal"] = parts[2]
        meta["title"] = '_'.join(parts[3:]) if len(parts) > 3 else parts[2]
    return meta


def _safe_json_parse(raw: str, fallback_label: str = "unknown") -> dict:
    """Parse JSON from LLM response, stripping markdown fences. Returns {} on failure."""
    text = raw.strip()
    if text.startswith('```'):
        idx = text.find('\n')
        if idx != -1:
            text = text[idx + 1:]
        text = text.removesuffix('```').strip()
    if text.startswith('`') and text.endswith('`'):
        text = text[1:-1].strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        logger.warning("Failed to parse JSON for %s: %.200s", fallback_label, text)
        return {}


def _yaml_dump(data: Any, indent: int = 0) -> list[str]:
    """Serialize Python value to YAML lines (no PyYAML dependency)."""
    pfx = "  " * indent
    if data is None:
        return [f"{pfx}null"]
    if isinstance(data, bool):
        return [f"{pfx}{str(data).lower()}"]
    if isinstance(data, (int, float)):
        return [f"{pfx}{data}"]
    if isinstance(data, str):
        if '\n' in data:
            return [f"{pfx}|"] + [f"{pfx}  {l}" for l in data.split('\n')]
        escaped = data.replace('"', '\\"')
        return [f'{pfx}"{escaped}"']
    if isinstance(data, dict):
        if not data:
            return [f"{pfx}{{}}"]
        lines: list[str] = []
        for k, v in data.items():
            ks = str(k)
            sub = _yaml_dump(v, indent + 1)
            if len(sub) == 1 and not sub[0].startswith("  " * (indent + 1) + "-"):
                lines.append(f"{pfx}{ks}: {sub[0].strip()}")
            else:
                lines.append(f"{pfx}{ks}:")
                lines.extend(sub)
        return lines
    if isinstance(data, (list, tuple)):
        if not data:
            return [f"{pfx}[]"]
        lines = []
        for item in data:
            sub = _yaml_dump(item, indent + 1)
            if len(sub) == 1:
                lines.append(f"{pfx}- {sub[0].strip()}")
            elif isinstance(item, dict) and len(sub) > 1:
                # Inline first key after - for readability
                first = sub[0].lstrip()
                remaining = sub[1:] 
                lines.append(f"{pfx}- {first}")
                lines.extend(remaining)
            else:
                lines.append(f"{pfx}-")
                lines.extend(sub)
        return lines
    return [f"{pfx}{data}"]


def _build_cfg_from_env() -> _LLMConfig:
    """Build LLM config from environment variables."""
    return _LLMConfig(
        api_key=os.environ.get("LLM_API_KEY", ""),
        api_base=os.environ.get("LLM_API_BASE", "https://api.deepseek.com/v1"),
        model=os.environ.get("LLM_MODEL", "deepseek/deepseek-v4-flash"),
    )


class PaperMdToInsights:
    """Extract structured insights from paper markdown using AI prompts."""

    @staticmethod
    def split_sections(md_text: str) -> dict[str, str]:
        """Split markdown into sections by header.

        Returns {abstract, introduction, results, discussion, methods} with
        lowercase keys. Missing sections get empty string.
        """
        sections: dict[str, str] = {k: "" for k in ("abstract", "introduction", "results", "discussion", "methods")}
        matches = list(_SECTION_RE.finditer(md_text))
        if not matches:
            sections["results"] = md_text.strip()
            return sections
        _section_key_map = {
            "summary": "abstract",
            "abstract": "abstract",
            "introduction": "introduction",
            "results": "results",
            "discussion": "discussion",
            "methods": "methods",
            "experimental procedures": "methods",
            "experimentalprocedures": "methods",
        }
        for i, m in enumerate(matches):
            raw_key = m.group(1).lower().strip()
            key = _section_key_map.get(raw_key, "unknown")
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(md_text)
            sections[key] = md_text[start:end].strip()
        return sections

    @staticmethod
    def extract_figure_blocks(results_text: str) -> list[str]:
        """Split Results text into figure blocks by Figure/Fig references."""
        if not results_text.strip():
            return []
        blocks = _FIGURE_RE.split(results_text)
        refs = _FIGURE_RE.findall(results_text)
        result: list[str] = []
        for i, block in enumerate(blocks[1:], start=1):
            header = refs[i - 1] if i - 1 < len(refs) else f"Figure_{i}"
            cleaned = block.strip().lstrip(".: ")
            combined = f"{header}. {cleaned}"
            if combined.strip():
                result.append(combined)
        if not result and results_text.strip():
            result.append(results_text.strip())
        return result

    @staticmethod
    def extract_metadata(abstract_text: str, cfg) -> dict:
        """Extract experimental design, key findings, data notes from abstract via LLM.

        Returns parsed dict or {} on failure. Skips gracefully if abstract is empty.
        """
        if not abstract_text.strip():
            logger.info("Abstract empty -- skipping metadata extraction")
            return {}
        raw = ai_query(PAPER_META_SYSTEM_PROMPT, PAPER_META_USER_TEMPLATE.format(abstract_text=abstract_text), cfg, expect_json=True)
        if raw is None:
            return {}
        return _safe_json_parse(raw, "abstract metadata")

    @staticmethod
    def extract_figure(figure_block: str, cfg) -> dict:
        """Extract structured info from figure legend via LLM. Returns {} on failure."""
        if not figure_block.strip():
            return {}
        raw = ai_query(PAPER_FIGURE_SYSTEM_PROMPT, PAPER_FIGURE_USER_TEMPLATE.format(figure_text=figure_block), cfg, expect_json=True)
        if raw is None:
            return {}
        return _safe_json_parse(raw, "figure")

    @staticmethod
    def extract_methods(methods_text: str, cfg) -> dict:
        """Extract methods, software, versions via LLM. Returns {} on failure."""
        if not methods_text.strip():
            logger.info("Methods empty -- skipping methods extraction")
            return {}
        raw = ai_query(PAPER_METHODS_SYSTEM_PROMPT, PAPER_METHODS_USER_TEMPLATE.format(methods_text=methods_text), cfg, expect_json=True)
        if raw is None:
            return {}
        return _safe_json_parse(raw, "methods")

    @staticmethod
    def merge_to_insights(meta: dict, figures: list[dict], methods: dict, md_path: str) -> dict:
        """Combine all extracted data into final insights dict.

        Deduplicates figures by id, sorts them, merges data notes from meta and methods.
        """
        seen: set = set()
        unique_figs: list[dict] = []
        for f in figures:
            fid = f.get("id", "")
            if fid and fid not in seen:
                seen.add(fid)
                unique_figs.append(f)
            elif not fid:
                unique_figs.append(f)
        unique_figs.sort(key=lambda f: f.get("id", ""))

        meta_notes = list(meta.get("data_notes", [])) if meta else []
        methods_notes = list(methods.get("data_notes", [])) if methods else []

        return {
            "paper_meta": _parse_filename_meta(md_path),
            "experimental_design": meta.get("experimental_design", {}) if meta else {},
            "key_findings": list(meta.get("key_findings", [])) if meta else [],
            "data_notes": meta_notes + [n for n in methods_notes if n not in meta_notes],
            "figures": unique_figs,
            "reproduction_status": {
                "pipeline_run": "not_started",
                "overall_match": None,
                "verified_figures": [],
                "notes": "",
            },
        }

    def run(self, md_path: str, cfg, output_path: str | None = None, force: bool = False) -> str:
        """Full pipeline: read, split, extract, merge, write.

        Returns path to written file, or "SKIPPED" if output exists and not force.

        Auto-creates subdirectory for output if needed.
        """
        md_obj = Path(md_path)
        if output_path:
            out_path = Path(output_path)
        else:
            out_path = md_obj.parent / md_obj.stem / "insights.yaml"

        if out_path.exists() and not force:
            logger.info("Output exists at %s, skipping (use --force to overwrite)", out_path)
            return "SKIPPED"

        md_text = md_obj.read_text(encoding="utf-8")
        sections = self.split_sections(md_text)
        logger.info("Found sections: %s", [k for k, v in sections.items() if v])

        meta: dict = {}
        if sections.get("abstract"):
            logger.info("Extracting metadata from abstract...")
            meta = self.extract_metadata(sections["abstract"], cfg)

        figures: list[dict] = []
        if sections.get("results"):
            blocks = self.extract_figure_blocks(sections["results"])
            logger.info("Found %d figure blocks", len(blocks))
            for i, block in enumerate(blocks):
                logger.info("Figure %d/%d...", i + 1, len(blocks))
                fig = self.extract_figure(block, cfg)
                if fig:
                    figures.append(fig)

        methods: dict = {}
        if sections.get("methods"):
            logger.info("Extracting methods...")
            methods = self.extract_methods(sections["methods"], cfg)

        insights = self.merge_to_insights(meta, figures, methods, md_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text('\n'.join(_yaml_dump(insights)) + '\n', encoding="utf-8")
        logger.info("Wrote insights to %s", out_path)
        return str(out_path)


def main() -> None:
    """CLI entry point for paper_md_to_insights."""
    parser = argparse.ArgumentParser(description="Extract structured paper insights using AI.")
    parser.add_argument("md_file", type=str, help="Paper markdown file")
    parser.add_argument("--output", "-o", type=str, default=None, help="Output path (default: <md_dir>/<stem>/insights.yaml)")
    parser.add_argument("--force", "-f", action="store_true", default=False, help="Overwrite existing output")
    parser.add_argument("--verbose", "-v", action="store_true", default=False, help="Enable debug logging")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s [%(name)s] %(message)s", stream=sys.stderr)

    # Use env-based config (LLM_API_KEY, LLM_API_BASE, LLM_MODEL)
    cfg = _build_cfg_from_env()
    logger.info("Using env-based LLM config")

    result = PaperMdToInsights().run(md_path=args.md_file, cfg=cfg, output_path=args.output, force=args.force)
    print("SKIPPED" if result == "SKIPPED" else f"Done: {result}")


if __name__ == "__main__":
    main()
