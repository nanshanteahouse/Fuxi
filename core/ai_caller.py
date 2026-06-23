#!/usr/bin/env python3
"""
ai_caller.py — 统一 LLM 调用模块
==================================

集中管理所有 AI API 调用，提供统一的接口模式。
所有步骤脚本通过本模块调用 LLM，避免重复的客户端创建和错误处理。

设计原则:
  - 使用 OpenAI SDK 作为统一客户端（兼容 OpenAI / Azure / 兼容 API）
  - cfg 为鸭子类型，只需提供 api_base, model, api_key, max_tokens, temperature 属性
  - 不在此模块中导入 scanpy 或任何生物学分析库

用法:
    from core.ai_caller import ai_query
    from core.ai_prompts import ANNOTATION_SYSTEM_PROMPT

    result = ai_query(ANNOTATION_SYSTEM_PROMPT, "User query...", cfg)
"""

import os
import sys
import time
import json
import urllib.request
from functools import lru_cache

@lru_cache(maxsize=1)
def _query_available_models(api_base: str, api_key: str) -> list[str]:
    """Query {api_base}/models and return available model names (cached).

    Uses urllib.request (stdlib) - no new dependencies.
    Returns empty list on any failure (network, auth, non-standard API).
    """
    url = api_base.rstrip('/') + '/models'
    req = urllib.request.Request(url)
    if api_key:
        req.add_header('Authorization', 'Bearer ' + api_key)
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        return [m['id'] for m in data.get('data', []) if isinstance(m, dict) and 'id' in m]
    except Exception:
        return []


def ai_query(system_prompt: str, user_prompt: str, cfg,
             log=None, expect_json: bool = False) -> str:
    """
    统一的 LLM 查询接口。

    使用 OpenAI SDK 构造聊天补全请求。cfg 是任意提供以下属性的对象:
      - api_base (str):      API 端点 URL
      - model (str):         模型名称 (如 'gpt-4o', 'claude-3-opus-20240229')
      - api_key (str):       API 密钥（留空则从 LLM_API_KEY 环境变量读取）
      - max_tokens (int):    最大输出 token 数
      - temperature (float): 采样温度 (0.0 ~ 2.0)

    参数:
        system_prompt: 系统角色提示词
        user_prompt:   用户输入/任务提示词
        cfg:           配置对象（鸭子类型，只需包含上述属性）
        log:           可选的 logger（用于 ATACseq 兼容）
        expect_json:   若为 True，去除 markdown 围栏并验证 JSON

    返回:
        模型生成的文本内容

    示例:
        >>> resp = ai_query("你是一个生物学专家", "注释以下聚类...", llm_cfg)
        >>> print(resp)
        '{"0": {"cell_type": "T cell", ...}}'
    """
    import logging
    if log is None:
        log = logging.getLogger(__name__)

    # ── 磁盘缓存（可配置开关，默认开启） ──────────────────────────────
    if getattr(cfg, 'ai_cache_responses', False):
        import hashlib
        cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', 'cache', 'ai_responses')
        os.makedirs(cache_dir, exist_ok=True)
        cache_key = hashlib.sha256(
            f"{cfg.model}:{system_prompt}:{user_prompt}".encode('utf-8')
        ).hexdigest()[:16]
        cache_path = os.path.join(cache_dir, f"{cache_key}.json")
        if os.path.exists(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    cached = json.load(f)
                return cached['response']
            except (json.JSONDecodeError, KeyError):
                pass  # Corrupted cache — proceed to live call

    from openai import OpenAI
    import openai

    client = OpenAI(
        api_key=cfg.api_key or os.getenv("LLM_API_KEY", "not-needed"),
        base_url=cfg.api_base,
    )

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user",   "content": user_prompt},
    ]

    # Retry loop: some vLLM deployments return content=None transiently
    max_retries = 3
    last_content = None
    for attempt in range(max_retries):
        # Build call kwargs based on thinking mode
        call_kwargs = dict(
            model=cfg.model,
            messages=messages,
            max_tokens=cfg.max_tokens,
            timeout=getattr(cfg, "timeout", None),
        )
        thinking_enabled = getattr(cfg, "thinking_enabled", True)
        if thinking_enabled:
            reasoning_effort = getattr(cfg, "reasoning_effort", "high")
            call_kwargs["extra_body"] = {"thinking": {"type": "enabled"}}
            call_kwargs["reasoning_effort"] = reasoning_effort
            # temperature/top_p/presence_penalty/frequency_penalty ignored in thinking mode
            # Auto-boost: a small max_tokens (e.g. 4096) lets the model burn the
            # entire budget on reasoning and return content=None after the retry
            # loop, failing the annotation silently. Floor at 32768 when thinking
            # is on so the model still has room for the final answer.
            if (call_kwargs["max_tokens"] or 0) < 32768:
                call_kwargs["max_tokens"] = 32768
        else:
            call_kwargs["temperature"] = cfg.temperature

        try:
            resp = client.chat.completions.create(**call_kwargs)
        except openai.APIStatusError as e:
            if e.status_code in (404, 422):
                models = _query_available_models(
                    api_base=cfg.api_base,
                    api_key=cfg.api_key or os.environ.get('LLM_API_KEY', ''),
                )
                if models:
                    msg = '[ai_caller] Model ' + repr(cfg.model) + ' not found. Available models: ' + str(models)
                    print(msg, file=sys.stderr)
                else:
                    msg = '[ai_caller] Model ' + repr(cfg.model) + ' not found (HTTP ' + str(e.status_code) + '). Could not query available models endpoint.'
                    print(msg, file=sys.stderr)
            raise
        content = resp.choices[0].message.content

        # ── JSON extraction (ATACseq compatibility) ──────────────────
        if content is not None and content.strip() and expect_json:
            content = content.strip()
            if content.startswith("```"):
                parts = content.split("```")
                content = parts[1] if len(parts) > 1 else content
                if content.startswith("json"):
                    content = content[4:]
                content = content.strip()
            try:
                json.loads(content)  # Validate
            except json.JSONDecodeError:
                pass  # Keep content as-is even if not valid JSON

        if content is not None and content.strip():
            # ── Write to cache on success ─────────────────────────────────
            if getattr(cfg, 'ai_cache_responses', False):
                try:
                    with open(cache_path, 'w') as f:
                        json.dump({
                            'model': cfg.model,
                            'created': time.time(),
                            'response': content,
                        }, f)
                except Exception:
                    pass  # Cache write failure is non-fatal
            return content

        # Remember reasoning_content as fallback
        rc = getattr(resp.choices[0].message, 'reasoning_content', None)
        if rc and rc.strip():
            last_content = rc.strip()
        if attempt < max_retries - 1:
            wait = 2 ** attempt
            msg = repr(content)[:100] if content is not None else 'None'
            print(f"[ai_caller] Empty response (content={msg}), retrying in {wait}s...", file=sys.stderr)
            time.sleep(wait)

    # Fallback: use reasoning_content if available
    if last_content:
        return last_content
    return None
