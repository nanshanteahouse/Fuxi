"""Extract JSON from LLM responses with possible markdown code fences."""

import re


def extract_json_block(text: str) -> str:
    """Extract a JSON block from text that may contain markdown code fences or prose.

    Strategy:
    1. Try to match a ```json ... ``` code block (multiline, non-greedy)
    2. Fallback: find first '{' to last '}' in the text
    3. Last resort: return raw text (let json.loads raise a readable error)
    """
    # Strategy 1: ```json ... ``` code fence
    pattern = r"```(?:json)?\s*\n?(.*?)```"
    match = re.search(pattern, text, re.DOTALL)
    if match:
        candidate = match.group(1).strip()
        # Verify it looks like JSON (starts with { or [)
        if candidate and (candidate[0] in ("{", "[")):
            return candidate

    # Strategy 2: first { to last }
    start = text.find("{")
    if start != -1:
        end = text.rfind("}")
        if end > start:
            return text[start : end + 1]

    # Strategy 3: try [ for arrays
    start = text.find("[")
    if start != -1:
        end = text.rfind("]")
        if end > start:
            return text[start : end + 1]

    return text
