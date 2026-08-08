"""Unit tests for ``core.tui.backends.pipeline.build_run_command``.

TUI execution is a thin wrapper over ``core/run_pipeline.py``; these tests
verify command construction only — no subprocess is ever spawned.
"""

import os
import sys

from core.tui.backends.pipeline import build_run_command


def test_build_run_command_basic():
    cmd = build_run_command(4, "rna", "config_rna.yaml")

    assert cmd[0] == sys.executable
    assert cmd[1].endswith("core/run_pipeline.py")
    assert cmd[cmd.index("--modality") + 1] == "rna"
    assert cmd[cmd.index("--step") + 1] == "4"
    assert "--config" in cmd
    assert "--cell-type" not in cmd


def test_build_run_command_cell_type_only_rna_step6():
    # cell_type forwarded for rna step 6, value preserved (incl. non-ASCII)
    cmd = build_run_command(6, "rna", "c.yaml", cell_type="Müller Glia")
    assert "--cell-type" in cmd
    assert cmd[cmd.index("--cell-type") + 1] == "Müller Glia"

    # NOT forwarded for rna step 5
    cmd = build_run_command(5, "rna", "c.yaml", cell_type="Müller Glia")
    assert "--cell-type" not in cmd

    # NOT forwarded for atac step 6
    cmd = build_run_command(6, "atac", "c.yaml", cell_type="Müller Glia")
    assert "--cell-type" not in cmd


def test_build_run_command_extra_args_appended():
    cmd = build_run_command(2, "rna", "c.yaml", extra_args=["--plot-only"])
    assert cmd[-1] == "--plot-only"


def test_build_run_command_abs_config():
    cmd = build_run_command(1, "rna", "relative/config.yaml")
    config_value = cmd[cmd.index("--config") + 1]
    assert os.path.isabs(config_value)
    assert config_value.endswith("relative/config.yaml")
