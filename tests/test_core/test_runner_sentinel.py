"""Sentinel completion-tracking tests for ``core/pipeline/runner.py`` (plan
h5ad-incremental-io Item 1.6).

In-place-writeback steps (RNA 06) and steps that may skip their own output
file (RNA 08) cannot be judged by their anchor checkpoint — that file exists
from an earlier step. Completion is tracked by a sentinel file
(``<base_checkpoint>.step{NN}_done``). These tests lock the runner-side
*reading* logic; *writing* the sentinel is the step scripts' job (T5/T6).

Baseline contract: without sentinels passed (or for modalities with no
sentinel steps), ``find_first_incomplete`` behaves exactly as before.
"""

from pathlib import Path

import pytest

from core.pipeline.runner import (
    RNA_CHECKPOINT_FILES,
    RNA_SENTINEL_FILES,
    RNA_STEPS,
    RNA_STEPS_WRITE_CHECKPOINT,
    _remove_anchored_sentinels,
    _sentinel_base,
    find_first_incomplete,
)

# ── Optional dependency guard (mcp) ───────────────────────────────────

_MCP_AVAILABLE: bool = False
try:
    import mcp  # noqa: F401

    _MCP_AVAILABLE = True
except ImportError:
    pass

skipif_no_mcp = pytest.mark.skipif(
    not _MCP_AVAILABLE,
    reason="mcp not installed (pip install 'mcp==2.0.0b2' for MCP server tests)",
)


def _write_nonempty(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("done\n")


def _rna_checkpoints_done(tmp_path: Path, upto: int) -> None:
    """Write non-empty checkpoint files for steps 0..upto-1 (the writing ones)."""
    for i in range(upto):
        _write_nonempty(tmp_path / RNA_CHECKPOINT_FILES[i])


# ── Baseline: no sentinels → pre-change behavior unchanged ──────────────


def test_baseline_no_sentinels_param_unchanged(tmp_path):
    """sentinels 缺省时行为与改造前完全一致：step 06/08 不参与判定。"""
    # 空目录 → 返回 0
    assert (
        find_first_incomplete(
            str(tmp_path), RNA_STEPS, RNA_CHECKPOINT_FILES, RNA_STEPS_WRITE_CHECKPOINT
        )
        == 0
    )
    # 0-5 完成、11 缺失 → 返回 11（跳过 6-10，与现状一致）
    _rna_checkpoints_done(tmp_path, 6)
    assert (
        find_first_incomplete(
            str(tmp_path), RNA_STEPS, RNA_CHECKPOINT_FILES, RNA_STEPS_WRITE_CHECKPOINT
        )
        == 11
    )
    # 0-5 + 11 全部完成 → len(steps)
    _write_nonempty(tmp_path / "11_grn.h5ad")
    assert find_first_incomplete(
        str(tmp_path), RNA_STEPS, RNA_CHECKPOINT_FILES, RNA_STEPS_WRITE_CHECKPOINT
    ) == len(RNA_STEPS)


def test_baseline_empty_sentinels_dict_unchanged(tmp_path):
    """sentinels={} 显式传入（ATAC/SPATIAL/BULK 情形）同样保持原行为。"""
    _rna_checkpoints_done(tmp_path, 6)
    assert (
        find_first_incomplete(
            str(tmp_path),
            RNA_STEPS,
            RNA_CHECKPOINT_FILES,
            RNA_STEPS_WRITE_CHECKPOINT,
            sentinels={},
        )
        == 11
    )


# ── Sentinel steps: presence drives completion ──────────────────────────


def test_missing_sentinel_marks_step_incomplete(tmp_path):
    """05_annotated.h5ad 已存在但无 step06 sentinel → 返回 6（不被误判为已完成）。"""
    _rna_checkpoints_done(tmp_path, 6)
    assert (
        find_first_incomplete(
            str(tmp_path),
            RNA_STEPS,
            RNA_CHECKPOINT_FILES,
            RNA_STEPS_WRITE_CHECKPOINT,
            RNA_SENTINEL_FILES,
        )
        == 6
    )


def test_sentinel_present_marks_step_done(tmp_path):
    """step06 sentinel 存在 → 越过 06；缺 step08 sentinel → 返回 8。"""
    _rna_checkpoints_done(tmp_path, 6)
    _write_nonempty(tmp_path / "05_annotated.h5ad.step06_done")
    assert (
        find_first_incomplete(
            str(tmp_path),
            RNA_STEPS,
            RNA_CHECKPOINT_FILES,
            RNA_STEPS_WRITE_CHECKPOINT,
            RNA_SENTINEL_FILES,
        )
        == 8
    )
    _write_nonempty(tmp_path / "05_final.h5ad.step08_done")
    assert (
        find_first_incomplete(
            str(tmp_path),
            RNA_STEPS,
            RNA_CHECKPOINT_FILES,
            RNA_STEPS_WRITE_CHECKPOINT,
            RNA_SENTINEL_FILES,
        )
        == 11
    )
    # 全链完成（含 11_grn）→ 全部完成
    _write_nonempty(tmp_path / "11_grn.h5ad")
    assert find_first_incomplete(
        str(tmp_path),
        RNA_STEPS,
        RNA_CHECKPOINT_FILES,
        RNA_STEPS_WRITE_CHECKPOINT,
        RNA_SENTINEL_FILES,
    ) == len(RNA_STEPS)


def test_empty_sentinel_counts_incomplete(tmp_path):
    """空 sentinel（0 字节）视为未完成——与 checkpoint 空文件判定一致。"""
    _rna_checkpoints_done(tmp_path, 6)
    (tmp_path / "05_annotated.h5ad.step06_done").touch()
    assert (
        find_first_incomplete(
            str(tmp_path),
            RNA_STEPS,
            RNA_CHECKPOINT_FILES,
            RNA_STEPS_WRITE_CHECKPOINT,
            RNA_SENTINEL_FILES,
        )
        == 6
    )


# ── Naming convention & cleanup ─────────────────────────────────────────


def test_sentinel_naming_convention():
    """``<base>.step{NN}_done`` → base 还原。"""
    assert _sentinel_base("05_annotated.h5ad.step06_done") == "05_annotated.h5ad"
    assert _sentinel_base("05_final.h5ad.step08_done") == "05_final.h5ad"


def test_remove_anchored_sentinels_only_anchor(tmp_path):
    """--cleanup 语义：删除锚定文件上的 sentinel，未锚定的保留。"""
    for f in ("05_annotated.h5ad.step06_done", "05_final.h5ad.step08_done", "05_annotated.h5ad"):
        _write_nonempty(tmp_path / f)
    _remove_anchored_sentinels(str(tmp_path), "05_annotated.h5ad", RNA_SENTINEL_FILES)
    assert not (tmp_path / "05_annotated.h5ad.step06_done").exists()
    assert (tmp_path / "05_final.h5ad.step08_done").exists()  # 未锚定 → 保留
    assert (tmp_path / "05_annotated.h5ad").exists()  # 锚定文件由调用方删除


def test_remove_anchored_sentinels_noop_without_sentinels(tmp_path):
    """无 sentinel 表时是 no-op。"""
    _write_nonempty(tmp_path / "05_annotated.h5ad.step06_done")
    _remove_anchored_sentinels(str(tmp_path), "05_annotated.h5ad", {})
    assert (tmp_path / "05_annotated.h5ad.step06_done").exists()


# ── MCP consistency (core/ai/mcp_tools/pipeline.py) ─────────────────────


@skipif_no_mcp
def test_mcp_sentinel_mirror_in_sync():
    """mcp_tools 的 _SENTINELS 镜像与 runner 的 RNA_SENTINEL_FILES 一致。"""
    from core.ai.mcp_tools.pipeline import _SENTINELS

    assert _SENTINELS["rna"] == RNA_SENTINEL_FILES
    # 其他 modality 无 sentinel → 调用方 fallback 到空表
    assert _SENTINELS.get("atac", {}) == {}


@skipif_no_mcp
def test_mcp_step_completed_uses_sentinel_for_sentinel_steps(tmp_path):
    """MCP 判定：sentinel 步骤只看 sentinel；普通步骤仍按 checkpoint。"""
    from core.ai.mcp_tools.pipeline import _step_completed

    # 锚定文件存在但 sentinel 缺失 → 不算完成（这是误报修复的核心）
    _write_nonempty(tmp_path / "05_annotated.h5ad")
    assert (
        _step_completed(str(tmp_path), "05_annotated.h5ad", "05_annotated.h5ad.step06_done")
        is False
    )
    # sentinel 出现 → 完成
    _write_nonempty(tmp_path / "05_annotated.h5ad.step06_done")
    assert (
        _step_completed(str(tmp_path), "05_annotated.h5ad", "05_annotated.h5ad.step06_done")
        is True
    )
    # 普通步骤：按 checkpoint 文件
    assert _step_completed(str(tmp_path), "00_raw.h5ad") is False
    _write_nonempty(tmp_path / "00_raw.h5ad")
    assert _step_completed(str(tmp_path), "00_raw.h5ad") is True
    # 无 marker → None（no_checkpoint 状态）
    assert _step_completed(str(tmp_path), "") is None
