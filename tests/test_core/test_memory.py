"""Unit tests for core.utils._memory — unified memory settings + guard rails."""

from __future__ import annotations

import pytest

from core.utils._memory import (
    check_memory_guard,
    estimate_step_peak,
    resolve_memory_budget_bytes,
    resolve_memory_settings,
)


class _Mem:
    def __init__(self, policy="speed", budget="auto", guard="warn"):
        self.policy = policy
        self.budget = budget
        self.guard = guard


class _Exec:
    def __init__(self, mem: _Mem | None = None, policy="speed", limit="auto"):
        self.memory = mem
        # legacy flat fields (survive pre-migration objects)
        self.memory_policy = policy
        self.memory_limit = limit


class _Cfg:
    def __init__(self, exec_cfg):
        self.execution = exec_cfg


def test_budget_parsing() -> None:
    assert resolve_memory_budget_bytes("auto") > 0  # psutil-based
    assert resolve_memory_budget_bytes("64GB") == int(64 * 2**30)
    assert resolve_memory_budget_bytes("128GiB") == int(128 * 2**30)
    assert resolve_memory_budget_bytes("512MB") == int(512 * 2**20)
    assert resolve_memory_budget_bytes("16g") == int(16 * 2**30)
    assert resolve_memory_budget_bytes("garbage") == 0


def test_resolve_memory_settings_nested() -> None:
    cfg = _Cfg(_Exec(_Mem(policy="balanced", budget="64GB", guard="block")))
    policy, budget, guard = resolve_memory_settings(cfg)
    assert policy == "balanced"
    assert budget == int(64 * 2**30)
    assert guard == "block"


def test_resolve_memory_settings_legacy_fallback() -> None:
    cfg = _Cfg(_Exec(mem=None, policy="memory", limit="32GB"))
    policy, budget, guard = resolve_memory_settings(cfg)
    assert policy == "memory"
    assert budget == int(32 * 2**30)
    assert guard == "warn"


def test_estimate_step_peak_monotonic() -> None:
    # Bigger input -> bigger estimate, on every step.
    small = (100_000, 20_000, 20_000_000)
    big = (2_000_000, 35_000, 3_500_000_000)
    for step in (1, 2, 3):
        s = estimate_step_peak(step, *small, policy="balanced", budget_bytes=80 * 2**30)
        b = estimate_step_peak(step, *big, policy="balanced", budget_bytes=80 * 2**30)
        assert b > s, f"step {step}: {s:.1f} -> {b:.1f} not monotonic"


def test_estimate_step03_policy_difference() -> None:
    # speed (dense PCA) must estimate far above balanced on a large matrix.
    n, g, nnz = 2_000_000, 35_000, 3_500_000_000
    speed = estimate_step_peak(3, n, g, nnz, policy="speed")
    bal = estimate_step_peak(3, n, g, nnz, policy="balanced")
    assert speed > bal * 2


def test_check_memory_guard_warn_continues() -> None:
    # warn: over-budget -> logs, returns True (run may proceed).
    assert check_memory_guard({1: 90.0}, 70 * 2**30, "warn") is True


def test_check_memory_guard_block_raises() -> None:
    with pytest.raises(RuntimeError):
        check_memory_guard({1: 90.0}, 70 * 2**30, "block")


def test_check_memory_guard_block_within_budget() -> None:
    assert check_memory_guard({1: 60.0, 2: 40.0}, 70 * 2**30, "block") is True


def test_check_memory_guard_off() -> None:
    assert check_memory_guard({1: 999.0}, 70 * 2**30, "off") is True
