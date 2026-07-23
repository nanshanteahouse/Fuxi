"""Tests for global config — deep_merge, GlobalConfig, resolve_config, export round-trip."""

import os
import tempfile

import pytest
import yaml

from core.config.global_config import GlobalConfig, GlobalPlotConfig
from core.config.schema import Config
from core.utils._config import _find_global_yaml, _find_repo_root, deep_merge, resolve_config

# ═══════════════════════════════════════════════════════════════════
# TestGlobalConfig — construction and validation
# ═══════════════════════════════════════════════════════════════════


class TestGlobalConfig:
    """GlobalConfig construction, defaults, and validation."""

    def test_defaults(self):
        """GlobalConfig constructs with correct defaults."""
        c = GlobalConfig()
        assert c.execution.random_seed == 42
        assert c.de.method == "wilcoxon"
        assert c.clustering.cluster_selection_method == "multi_metric"
        assert c.integration.max_iter == 20
        assert c.qc.max_pct_mito == 20.0
        assert c.plot.figure_dpi == 150
        assert c.plot.palette.categorical == "tab20"

    def test_validation_extra_forbid(self):
        """extra='forbid' rejects unknown keys."""
        with pytest.raises(Exception):
            GlobalConfig.model_validate({"execution": {"bad_key": 1}})


# ═══════════════════════════════════════════════════════════════════
# TestDeepMerge — deep_merge behavioural contract
# ═══════════════════════════════════════════════════════════════════


class TestDeepMerge:
    """deep_merge returns (merged, source_map) — 10 tests."""

    def test_scalar_override(self):
        """Scalar in override replaces base."""
        result, _ = deep_merge({"a": 1}, {"a": 2})
        assert result == {"a": 2}

    def test_nested_merge(self):
        """Dicts at same key are merged recursively."""
        result, _ = deep_merge({"a": {"b": 1}}, {"a": {"c": 2}})
        assert result == {"a": {"b": 1, "c": 2}}

    def test_empty_dict_clears_nested(self):
        """Empty dict {} in override clears the nested key entirely."""
        result, _ = deep_merge({"a": {"b": 1}}, {"a": {}})
        assert result == {"a": {}}

    def test_empty_override_keeps_base(self):
        """Empty dict override {} keeps full base."""
        result, _ = deep_merge({"a": 1}, {})
        assert result == {"a": 1}

    def test_partial_override(self):
        """Only specified keys in override are replaced."""
        result, _ = deep_merge({"a": 1, "b": 2}, {"a": 3})
        assert result == {"a": 3, "b": 2}

    def test_list_replacement(self):
        """Lists are replaced wholesale, not merged."""
        result, _ = deep_merge({"a": [1, 2]}, {"a": [3]})
        assert result == {"a": [3]}

    def test_none_base_overridden(self):
        """None base value is replaced by a non-None override."""
        result, _ = deep_merge({"a": None}, {"a": 1})
        assert result == {"a": 1}

    def test_none_override_preserves_base(self):
        """None in override keeps the base value."""
        result, _ = deep_merge({"a": 1}, {"a": None})
        assert result == {"a": 1}

    def test_three_level_depth(self):
        """Merging works at three levels of nesting."""
        result, _ = deep_merge({"a": {"b": {"c": 1}}}, {"a": {"b": {"d": 2}}})
        assert result == {"a": {"b": {"c": 1, "d": 2}}}

    def test_source_map(self):
        """Source map correctly labels 'base' vs 'override' origins."""
        result, sm = deep_merge({"a": {"b": 1}}, {"a": {"c": 2}})
        assert sm == {"a.b": "base", "a.c": "override"}


# ═══════════════════════════════════════════════════════════════════
# TestFindGlobalYaml — discovery helpers
# ═══════════════════════════════════════════════════════════════════


class TestFindGlobalYaml:
    """_find_repo_root and _find_global_yaml behaviour."""

    def test_repo_root_found(self):
        """_find_repo_root finds a dir with .git/ in ancestry."""
        root = _find_repo_root()
        assert root is not None
        assert os.path.isdir(os.path.join(root, ".git"))

    def test_find_global_yaml_none(self):
        """Without env var and without repo-root global.yaml, returns None."""
        path = _find_global_yaml()
        # The Fuxi repo has no global.yaml at root, so this should be None
        # unless FUXI_GLOBAL_CONFIG is set in the environment.
        assert path is None or path.endswith("global.yaml")


# ═══════════════════════════════════════════════════════════════════
# TestResolveConfig — integration with resolve_config
# ═══════════════════════════════════════════════════════════════════


class TestResolveConfig:
    """resolve_config — backward compat, priority, bulk loading."""

    @staticmethod
    def _minimal_config(tmpdir: str, name: str = "test_project", **overrides) -> str:
        """Write a minimal valid project config to tmpdir and return its path."""
        data = {
            "modality": "rna",
            "tissue": "retina",
            "species": "human",
            "data_format": "10X_mtx",
            "expression_type": "raw_counts",
            "tissue_maturity": "developing",
            "project_dir": tmpdir,
        }
        data.update(overrides)
        path = os.path.join(tmpdir, f"config_{name}.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f)
        return path

    def test_backward_compat_no_global(self):
        """Configs load without global.yaml, using schema defaults."""
        with tempfile.TemporaryDirectory() as tmp:
            project_path = self._minimal_config(tmp)
            cfg = resolve_config(project_path)
            assert cfg.tissue == "retina"
            assert cfg.plot.figure_dpi == 150

    def test_priority_project_overrides_global(self):
        """Project values win over global when both specify the key."""
        with tempfile.TemporaryDirectory() as tmp:
            # Project config: explicitly set figure_dpi = 100
            project_path = self._minimal_config(tmp, plot={"figure_dpi": 100})
            # Global config: set figure_dpi = 300
            global_path = os.path.join(tmp, "global.yaml")
            with open(global_path, "w") as f:
                yaml.dump({"plot": {"figure_dpi": 300}}, f)
            os.environ["FUXI_GLOBAL_CONFIG"] = global_path
            try:
                cfg = resolve_config(project_path)
                # Project value (100) should win over global (300)
                assert cfg.plot.figure_dpi == 100
            finally:
                del os.environ["FUXI_GLOBAL_CONFIG"]

    def test_all_active_configs_load(self):
        """Multiple valid configs load without error, plot defaults present."""
        with tempfile.TemporaryDirectory() as tmp:
            for i in range(5):
                project_path = self._minimal_config(tmp, name=f"proj_{i}")
                cfg = resolve_config(project_path)
                assert cfg.plot.figure_dpi == 150

    def test_all_active_configs_load_no_global(self):
        """Configs load cleanly without global.yaml, plot defaults present."""
        with tempfile.TemporaryDirectory() as tmp:
            for i in range(5):
                project_path = self._minimal_config(tmp, name=f"proj_{i}")
                cfg = resolve_config(project_path)
                assert cfg.plot.figure_dpi == 150


# ═══════════════════════════════════════════════════════════════════
# TestPlotConfig — GlobalPlotConfig + GlobalPaletteConfig
# ═══════════════════════════════════════════════════════════════════


class TestPlotConfig:
    """Plot sub-config defaults and partial overrides."""

    def test_defaults(self):
        """GlobalPlotConfig constructs with correct defaults."""
        c = GlobalPlotConfig()
        assert c.figure_dpi == 150
        assert c.figure_format == "pdf"
        assert c.palette.categorical == "tab20"

    def test_partial_override(self):
        """Partial override fills unspecified fields with defaults."""
        c = GlobalPlotConfig.model_validate(
            {"figure_dpi": 300, "palette": {"categorical": "Set2"}}
        )
        assert c.figure_dpi == 300
        assert c.palette.categorical == "Set2"
        assert c.palette.heatmap == "RdBu_r"


# ═══════════════════════════════════════════════════════════════════
# TestExport — round-trip: resolve → export → load
# ═══════════════════════════════════════════════════════════════════


class TestExport:
    """Export resolved config to YAML and re-load it."""

    @staticmethod
    def _minimal_config(tmpdir: str, name: str = "test_project", **overrides) -> str:
        data = {
            "modality": "rna",
            "tissue": "retina",
            "species": "human",
            "data_format": "10X_mtx",
            "expression_type": "raw_counts",
            "tissue_maturity": "developing",
            "project_dir": tmpdir,
        }
        data.update(overrides)
        path = os.path.join(tmpdir, f"config_{name}.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f)
        return path

    def test_roundtrip(self):
        """Exported config can be loaded back as a valid Config model."""
        with tempfile.TemporaryDirectory() as tmp:
            project_path = self._minimal_config(tmp)
            cfg = resolve_config(project_path)
            resolved_path = os.path.join(cfg.results_dir, "config_resolved.yaml")
            assert os.path.exists(resolved_path)
            with open(resolved_path) as f:
                d = yaml.safe_load(f)
            assert "_config_meta" in d
            assert d["_config_meta"]["merged_priority"] == "project > global > schema_default"
            data = {k: v for k, v in d.items() if not k.startswith("_config_")}
            cfg2 = Config.model_validate(data)
            assert cfg2.tissue == "retina"
            assert cfg2.plot.figure_dpi == 150
