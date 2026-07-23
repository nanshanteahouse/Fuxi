"""Test require_*() guards in ``core/utils/_optional.py``.

Uses mocking — no real imports, no network, no I/O.  Verifies that each
``ImportError`` message mentions the correct ``pip install fuxi[...]`` hint
so the hint ↔ extras contract does not drift again.
"""

import types
from unittest.mock import patch

import pytest

# (function_name, expected_pip_hint_substring)
REQUIRE_FUNCTIONS = [
    ("require_scvi", "fuxi[scvi]"),
    ("require_celltypist", "fuxi[celltypist]"),
    ("require_scvelo", "fuxi[scvelo]"),
    ("require_sccoda", "fuxi[sccoda]"),
    ("require_cellbender", "fuxi[cellbender]"),
    ("require_soupx", "fuxi[soupx]"),
]


@pytest.fixture(autouse=True)
def _reset_availability_cache():
    """Reset module-level ``_*_available`` caches before every test.

    Without this, a test that mocks ``find_spec → None`` would permanently
    poison the cache for subsequent parametrised iterations.
    """
    import core.utils._optional as _opt

    for attr in list(vars(_opt)):
        if attr.startswith("_") and attr.endswith("_available"):
            setattr(_opt, attr, None)
    yield


def _load(name: str):
    """Dynamically import a ``require_*`` function by name."""
    import core.utils._optional as _opt

    return getattr(_opt, name)


class TestRequireErrors:
    """ImportError hints must mention the correct pip install command."""

    # ── 1. pip install hint ─────────────────────────────────────────────

    @pytest.mark.parametrize("func_name,pip_hint", REQUIRE_FUNCTIONS)
    def test_hint_mentions_correct_extra(self, func_name: str, pip_hint: str):
        """When the package is missing, the error message must mention ``fuxi[...]``."""
        with patch("importlib.util.find_spec", return_value=None):
            with pytest.raises(ImportError) as exc:
                _load(func_name)()
        assert pip_hint in str(exc.value), (
            f"{func_name}: expected hint to contain '{pip_hint}', got: '{str(exc.value)[:120]}'"
        )

    # ── 2. custom feature string ────────────────────────────────────────

    @pytest.mark.parametrize("func_name,pip_hint", REQUIRE_FUNCTIONS)
    def test_custom_feature_appears_in_message(self, func_name: str, pip_hint: str):
        """A custom ``feature=`` kwarg must show up in the error message."""
        custom = f"_test_feature_{func_name}_"
        with patch("importlib.util.find_spec", return_value=None):
            with pytest.raises(ImportError) as exc:
                _load(func_name)(feature=custom)
        assert custom in str(exc.value), (
            f"{func_name}: custom feature '{custom}' should appear, got: '{str(exc.value)[:120]}'"
        )

    @pytest.mark.parametrize("func_name,pip_hint", REQUIRE_FUNCTIONS)
    def test_passes_when_package_available(self, func_name: str, pip_hint: str):
        """When ``find_spec`` returns a truthy value, no error is raised."""
        with patch("importlib.util.find_spec", return_value=types.ModuleType("fake")):
            # Should not raise
            _load(func_name)()
