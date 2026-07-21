"""Prompt YAML loader — loads and caches LLM prompts from YAML files."""

from pathlib import Path

import yaml

_PROMPTS_DIR = Path(__file__).parent
_cache: dict[str, dict] = {}


def load_prompt(name: str) -> dict[str, str]:
    """Load a prompt YAML file. Returns dict with 'system' and 'user_template' keys.

    Caches results in memory — call load_prompt() multiple times with no I/O penalty.
    """
    if name not in _cache:
        path = _PROMPTS_DIR / f"{name}.yaml"
        if not path.is_file():
            raise FileNotFoundError(f"Prompt file not found: {path}")
        with open(path, encoding="utf-8") as f:
            _cache[name] = yaml.safe_load(f)
    return _cache[name]
