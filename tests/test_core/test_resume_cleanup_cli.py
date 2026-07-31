"""CLI-level ``--resume`` / ``--cleanup`` behavior tests — plan h5ad-incremental-io
task 11 ②.

Every assertion here is subprocess-level evidence: the pipeline runner is invoked
as ``python -m core.run_pipeline --modality rna ...`` against a tmp project with
hand-staged checkpoints, and the resume/cleanup semantics are judged by the
runner's own stdout.

Scenarios locked:

1. 00-05 checkpoints staged but no step-06 sentinel → ``--resume`` starts at
   step [06] (pre-sentinel behavior would have jumped to [11]). The chain stops
   fast afterwards because the deliberately minimal ``05_annotated.h5ad`` has no
   annotation columns → step 07 exits 1.
2. Plus a non-empty ``05_annotated.h5ad.step06_done`` → ``--resume`` skips 6 and
   starts at step [08]. Step 08 crashes fast on the minimal h5ad (no X_pca).
3. ``--step 11 --cleanup`` with ``grn.run=false``: the runner deletes the
   step-11 dependency ``05_annotated.h5ad`` AND its anchored sentinel
   ``05_annotated.h5ad.step06_done`` together ("同生共死"); an unanchored
   sentinel (``05_final.h5ad.step08_done``) survives.

Sentinel contract used: ``<base_checkpoint>.step{NN}_done``, non-empty content.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Checkpoint files the runner checks for steps 0..5 (RNA).
_CHECKPOINT_00_05 = [
    "00_raw.h5ad",
    "01_doublet.h5ad",
    "02_qc.h5ad",
    "03_integrated.h5ad",
    "04_clustered.h5ad",
    "05_annotated.h5ad",
]


def _stage_checkpoint(h5ad_dir: Path, name: str, *, with_cell_type: bool) -> None:
    """Write a tiny real (gzip) h5ad so the checkpoint passes existence checks."""
    rng = np.random.default_rng(0)
    adata = sc.AnnData(X=rng.random((60, 30)).astype(np.float32))
    if with_cell_type:
        adata.obs["cell_type"] = pd.Categorical([f"T{i % 2}" for i in range(60)])
    adata.write(str(h5ad_dir / name), compression="gzip")


def _build_config(tmp_path: Path) -> dict:
    """Minimal RNA config rooted under tmp_path (dirs created by resolve_config)."""
    return {
        "modality": "rna",
        "tissue": "test",
        "species": "mouse",
        "project_dir": str(tmp_path),
        "h5ad_dir": str(tmp_path / "results" / "h5ad"),
        "figure_dir": str(tmp_path / "results" / "figures"),
        "table_dir": str(tmp_path / "results" / "tables"),
        "log_dir": str(tmp_path / "logs"),
        "execution": {"device": "cpu", "random_seed": 42, "n_jobs": 1},
        "marker": {"subcluster_types": []},  # step 06 skips (exit 2)
        "grn": {"run": False},  # step 11 returns immediately when reached
    }


def _stage_00_05(tmp_path: Path) -> Path:
    """Stage non-empty 00-05 checkpoints; 05_annotated has no annotation columns.

    Returns the h5ad dir. The minimal 05_annotated (X only) makes step 07 fail
    fast (no cell_type/leiden → sys.exit 1) and step 08 fail fast (no X_pca).
    """
    h5ad_dir = tmp_path / "results" / "h5ad"
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    for name in _CHECKPOINT_00_05:
        _stage_checkpoint(h5ad_dir, name, with_cell_type=False)
    return h5ad_dir


def _write_config(tmp_path: Path) -> Path:
    import yaml

    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(_build_config(tmp_path), f)
    return config_path


def _run_cli(argv: list[str], timeout: int = 150):
    """Run the runner as a subprocess; return (returncode_or_None, stdout).

    ``returncode`` is None when the subprocess had to be killed by the timeout
    (stdout is still preserved for assertion diagnostics).
    """
    cmd = [sys.executable, "-m", "core.run_pipeline", *argv]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(_REPO_ROOT),
            timeout=timeout,
        )
        return proc.returncode, proc.stdout
    except subprocess.TimeoutExpired as exc:  # pragma: no cover - safety net
        return None, exc.stdout or ""


def _assert_cleanup_evidence(stdout: str, h5ad_dir: Path) -> None:
    """Both anchored files must be gone; the unanchored sentinel must survive."""
    assert "[run]   Cleaned up: 05_annotated.h5ad" in stdout, stdout[-2000:]
    assert "[run]   Cleaned up: 05_annotated.h5ad.step06_done" in stdout, stdout[-2000:]
    assert not (h5ad_dir / "05_annotated.h5ad").exists(), "anchored file survived cleanup"
    assert not (h5ad_dir / "05_annotated.h5ad.step06_done").exists(), (
        "anchored sentinel survived cleanup"
    )
    assert (h5ad_dir / "05_final.h5ad.step08_done").exists(), (
        "unanchored sentinel must survive cleanup"
    )


def test_resume_without_sentinel_starts_at_step06(tmp_path):
    """00-05 staged, no step-06 sentinel → resume begins at step [06].

    This is the regression the sentinel mechanism fixes: pre-sentinel the runner
    returned 11 here (skipping steps 6-10).
    """
    h5ad_dir = _stage_00_05(tmp_path)
    config_path = _write_config(tmp_path)

    rc, stdout = _run_cli(["--modality", "rna", "--resume", "--config", str(config_path)])

    assert "[run] Resuming from step [06]" in stdout, stdout[-3000:]
    # Step 06 is attempted and skipped (exit 2), then step 07 fails fast on the
    # annotation-less 05_annotated — proving the chain really started at 6.
    assert "[run] Step [06] skipped" in stdout, stdout[-3000:]
    assert "Step [07] failed" in stdout, stdout[-3000:]
    assert rc is not None, "subprocess timed out before the resume decision was provable"
    assert rc != 0, "expected deliberate step-07 failure to stop the run"
    assert not (h5ad_dir / "05_annotated.h5ad.step06_done").exists(), (
        "no sentinel may be written without a writeback"
    )


def test_resume_with_step06_sentinel_skips_to_step08(tmp_path):
    """Non-empty step-06 sentinel → resume skips 6 and starts at step [08]."""
    h5ad_dir = _stage_00_05(tmp_path)
    (h5ad_dir / "05_annotated.h5ad.step06_done").write_text("done\n")
    config_path = _write_config(tmp_path)

    rc, stdout = _run_cli(["--modality", "rna", "--resume", "--config", str(config_path)])

    assert "[run] Resuming from step [08]" in stdout, stdout[-3000:]
    assert "[run] Resuming from step [06]" not in stdout, stdout[-3000:]
    # Step 08 runs and crashes fast (minimal h5ad has no X_pca for neighbors).
    assert "[run] [RNA] Step [08]:" in stdout, stdout[-3000:]
    assert rc is not None, "subprocess timed out before the resume decision was provable"
    assert rc != 0, "expected deliberate step-08 failure to stop the run"


def test_cleanup_removes_sentinel_with_anchored_file(tmp_path):
    """--cleanup: 05_annotated.h5ad.step06_done dies with 05_annotated.h5ad.

    Step 11's dependency is 05_annotated.h5ad (in RNA_STEPS_WRITE_CHECKPOINT), so
    a completing ``--step 11`` run with ``grn.run=false`` triggers the cleanup
    block: the anchor file is removed and ``_remove_anchored_sentinels`` removes
    only the sentinel anchored to it.
    """
    h5ad_dir = tmp_path / "results" / "h5ad"
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    _stage_checkpoint(h5ad_dir, "05_annotated.h5ad", with_cell_type=False)
    (h5ad_dir / "05_annotated.h5ad.step06_done").write_text("done\n")
    (h5ad_dir / "05_final.h5ad.step08_done").write_text("done\n")
    config_path = _write_config(tmp_path)

    rc, stdout = _run_cli(
        ["--modality", "rna", "--step", "11", "--cleanup", "--config", str(config_path)]
    )

    assert rc == 0, f"step 11 (grn.run=false) should complete: {stdout[-3000:]}"
    _assert_cleanup_evidence(stdout, h5ad_dir)
