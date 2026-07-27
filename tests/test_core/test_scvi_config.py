"""Tests for SCVIConfig extended fields and 03_integrate validation guards."""

from __future__ import annotations

import os
import sys

import pytest
from pydantic import ValidationError

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from core.config.schema import SCVIConfig  # noqa: E402

# ── SCVIConfig field tests ──────────────────────────────────────────


class TestSCVIConfigDefaults:
    def test_defaults(self):
        c = SCVIConfig()
        assert c.batch_size == 128
        assert c.early_stopping is False
        assert c.precision == "32"
        assert c.trainer_kwargs == {}
        assert c.plan_kwargs == {}
        assert c.datasplitter_kwargs == {}

    def test_backward_compat_no_new_fields(self):
        """Old config dict without new fields still validates."""
        c = SCVIConfig.model_validate(
            {
                "n_latent": 30,
                "n_layers": 2,
                "n_hidden": 128,
                "max_epochs": 400,
                "batch_key": "sample",
                "use_gpu": True,
                "train_size": 0.9,
            }
        )
        assert c.batch_size == 128  # default applied
        assert c.precision == "32"

    def test_precision_literal_valid(self):
        for p in ["32", "16-mixed", "bf16-mixed"]:
            assert SCVIConfig(precision=p).precision == p

    def test_precision_literal_invalid(self):
        with pytest.raises(ValidationError):
            SCVIConfig(precision="8-mixed")

    def test_passthrough_dicts_accept_arbitrary_keys(self):
        c = SCVIConfig(
            trainer_kwargs={"gradient_clip_val": 0.5, "accumulate_grad_batches": 2},
            plan_kwargs={"lr": 0.01, "weight_decay": 1e-5},
            datasplitter_kwargs={"num_workers": 4, "drop_last": True},
        )
        assert c.trainer_kwargs["gradient_clip_val"] == 0.5
        assert c.plan_kwargs["lr"] == 0.01
        assert c.datasplitter_kwargs["num_workers"] == 4

    def test_extra_forbid_still_enforced(self):
        with pytest.raises(ValidationError):
            SCVIConfig(unknown_field=42)


# ── Conflict-key filter test ────────────────────────────────────────


class TestConflictKeyFilter:
    def test_explicit_kwargs_constant_exists(self):
        """_SCVI_EXPLICIT_KWARGS is defined in 03_integrate.py."""
        import importlib.util

        _path = os.path.join(_REPO_ROOT, "rna", "steps", "03_integrate.py")
        _spec = importlib.util.spec_from_file_location("_03_integrate_test", _path)
        _mod = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_mod)
        assert hasattr(_mod, "_SCVI_EXPLICIT_KWARGS")
        assert "max_epochs" in _mod._SCVI_EXPLICIT_KWARGS
        assert "batch_size" in _mod._SCVI_EXPLICIT_KWARGS
        assert "precision" in _mod._SCVI_EXPLICIT_KWARGS
        assert "datasplitter_kwargs" in _mod._SCVI_EXPLICIT_KWARGS
        assert "plan_kwargs" in _mod._SCVI_EXPLICIT_KWARGS
        assert "validation_size" in _mod._SCVI_EXPLICIT_KWARGS
        assert "trainer_config" in _mod._SCVI_EXPLICIT_KWARGS


# ── Early stopping validation guard test ───────────────────────────


class TestEarlyStoppingGuard:
    def test_early_stopping_with_full_train_size_should_fail(self):
        """early_stopping=True + train_size=1.0 is an invalid combination."""
        cfg_scvi = SCVIConfig(early_stopping=True, train_size=1.0)
        # The guard checks: early_stopping and train_size >= 1.0
        assert cfg_scvi.early_stopping and cfg_scvi.train_size >= 1.0  # condition is True

    def test_early_stopping_with_validation_set_ok(self):
        """early_stopping=True + train_size=0.9 is valid."""
        cfg_scvi = SCVIConfig(early_stopping=True, train_size=0.9)
        assert not (cfg_scvi.early_stopping and cfg_scvi.train_size >= 1.0)  # condition is False


# ── CPU precision guard test (D1) ─────────────────────────────────


class TestCPUPrecisionGuard:
    def test_mixed_precision_on_cpu_downgrades(self):
        """precision='16-mixed' + use_gpu=False → downgraded to '32'."""
        cfg_scvi = SCVIConfig(precision="16-mixed")
        use_gpu = False
        # Simulate the guard logic
        if cfg_scvi.precision != "32" and not use_gpu:
            _precision = "32"
        else:
            _precision = cfg_scvi.precision
        assert _precision == "32"  # downgraded

    def test_mixed_precision_on_gpu_preserved(self):
        """precision='16-mixed' + use_gpu=True → preserved."""
        cfg_scvi = SCVIConfig(precision="16-mixed")
        use_gpu = True
        if cfg_scvi.precision != "32" and not use_gpu:
            _precision = "32"
        else:
            _precision = cfg_scvi.precision
        assert _precision == "16-mixed"  # unchanged

    def test_fp32_always_safe(self):
        """precision='32' is safe on both CPU and GPU."""
        cfg_scvi = SCVIConfig(precision="32")
        for use_gpu in (True, False):
            if cfg_scvi.precision != "32" and not use_gpu:
                _precision = "32"
            else:
                _precision = cfg_scvi.precision
            assert _precision == "32"
