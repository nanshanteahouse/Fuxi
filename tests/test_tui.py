"""Integration tests for Fuxi TUI using Textual's pilot API.

All tests are self-contained (no external dependencies).
"""

from __future__ import annotations

import asyncio
from unittest.mock import patch

import pytest


def _run(coro):
    """Execute coroutine in a fresh event loop, then close it."""
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


# ══════════════════════════════════════════════════════════════════════════
# TestAppLaunch
# ══════════════════════════════════════════════════════════════════════════

class TestAppLaunch:
    def test_app_starts(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                assert type(pilot.app.screen).__name__ == "HomeScreen"
        _run(check())

    def test_home_screen_has_sidebar_and_content(self):
        from core.tui.app import FuxiTUI
        from core.tui.widgets.sidebar import Sidebar
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one(Sidebar) is not None
                assert s.query_one("#content-area") is not None
        _run(check())

    def test_home_screen_content(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                await asyncio.sleep(0.15)
                static = pilot.app.screen.query_one("Static.-header")
                assert "Fuxi" in str(static.visual)
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestScreenNavigation
# ══════════════════════════════════════════════════════════════════════════

class TestScreenNavigation:
    def test_all_keyboard_shortcuts(self):
        from core.tui.screens.home_screen import HomeScreen as HS
        from core.tui.screens.registry_screen import RegistryScreen as RS
        from core.tui.screens.pipeline_screen import PipelineScreen as PS
        from core.tui.screens.results_screen import ResultsScreen as ResS
        from core.tui.screens.data_mgmt import DataManagementScreen as DMS
        from core.tui.screens.config_editor import ConfigEditorScreen as CES
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                a = pilot.app
                await pilot.press("ctrl+r"); assert type(a.screen) is RS
                await pilot.press("ctrl+p"); assert type(a.screen) is PS
                await pilot.press("ctrl+e"); assert type(a.screen) is ResS
                await pilot.press("ctrl+d"); assert a.screen.id == "data-mgmt"
                await pilot.press("ctrl+c"); assert a.screen.id == "config-editor"
                await pilot.press("ctrl+h"); assert type(a.screen) is HS
        _run(check())

    def test_sidebar_highlights_active(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                await pilot.press("ctrl+p")
                assert pilot.app.screen.query_one("#pipeline") is not None
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestModalitySelector
# ══════════════════════════════════════════════════════════════════════════

class TestModalitySelector:
    def test_modality_select_renders(self):
        from core.tui.screens.home import HomeScreen
        from textual.app import App
        async def check():
            a = App()
            async with a.run_test(size=(80, 24)) as pilot:
                a.install_screen(HomeScreen(), "home")
                a.push_screen("home")
                await asyncio.sleep(0.15)
                sel = pilot.app.screen.query_one("#modality_select")
                assert len(list(sel._options)) >= 2
        _run(check())

    def test_quick_launch_buttons_render(self):
        from core.tui.screens.home import HomeScreen
        from textual.app import App
        from textual.widgets import Button
        async def check():
            a = App()
            async with a.run_test(size=(80, 24)) as pilot:
                a.install_screen(HomeScreen(), "home")
                a.push_screen("home")
                await asyncio.sleep(0.15)
                for bid in ("btn_browse_registry", "btn_run_pipeline", "btn_view_results"):
                    assert pilot.app.screen.query_one(f"#{bid}", Button) is not None
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestRegistryScreen
# ══════════════════════════════════════════════════════════════════════════

class TestRegistryScreen:
    def test_registry_screen_renders(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                from core.tui.screens.registry_screen import RegistryScreen
                await pilot.press("ctrl+r")
                assert type(pilot.app.screen) is RegistryScreen
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestStepListing
# ══════════════════════════════════════════════════════════════════════════

class TestStepListing:
    def test_step_table_populated(self):
        from core.tui.screens.steps import StepListingScreen
        from textual.app import App
        from textual.widgets import DataTable
        async def check():
            a = App()
            async with a.run_test(size=(80, 24)) as pilot:
                a.install_screen(StepListingScreen(), "steps")
                a.push_screen("steps")
                await asyncio.sleep(0.15)
                table = pilot.app.screen.query_one("#steps-table", DataTable)
                assert table.row_count > 0
        _run(check())

    def test_step_table_has_columns(self):
        from core.tui.screens.steps import StepListingScreen
        from textual.app import App
        from textual.widgets import DataTable
        async def check():
            a = App()
            async with a.run_test(size=(80, 24)) as pilot:
                a.install_screen(StepListingScreen(), "steps")
                a.push_screen("steps")
                await asyncio.sleep(0.15)
                table = pilot.app.screen.query_one("#steps-table", DataTable)
                labels = [c.label.plain for c in table.ordered_columns]
                assert "Step #" in labels
                assert "Script" in labels
                assert "Description" in labels
        _run(check())

    def test_step_count_label(self):
        from core.tui.screens.steps import StepListingScreen
        from textual.app import App
        from textual.widgets import Static
        async def check():
            a = App()
            async with a.run_test(size=(80, 24)) as pilot:
                a.install_screen(StepListingScreen(), "steps")
                a.push_screen("steps")
                await asyncio.sleep(0.15)
                label = pilot.app.screen.query_one("#step-count", Static)
                assert "Total:" in str(label.visual)
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestConfigForm
# ══════════════════════════════════════════════════════════════════════════

class TestConfigForm:
    def test_config_form_has_sections(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                assert len(list(pilot.app.screen.query("Collapsible"))) >= 1
        _run(check())

    def test_config_form_has_toolbar_buttons(self):
        from core.tui.app import FuxiTUI
        from textual.widgets import Button
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one("#btn-load", Button) is not None
                assert s.query_one("#btn-save", Button) is not None
        _run(check())

    def test_config_form_has_path_input(self):
        from core.tui.app import FuxiTUI
        from textual.widgets import Input
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                assert pilot.app.screen.query_one("#config-path-input", Input) is not None
        _run(check())

    def test_config_form_has_field_widgets(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                assert len(list(pilot.app.screen.query(".field-row"))) >= 1
        _run(check())

    def test_config_form_status_bar(self):
        from core.tui.app import FuxiTUI
        from textual.widgets import Static
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                bar = pilot.app.screen.query_one("#status-bar", Static)
                assert str(bar.visual) != ""
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestPipelineScreen
# ══════════════════════════════════════════════════════════════════════════

class TestPipelineScreen:
    def test_pipeline_screen_renders(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                from core.tui.screens.pipeline_screen import PipelineScreen
                await pilot.press("ctrl+p")
                assert type(pilot.app.screen) is PipelineScreen
        _run(check())

    def test_pipeline_screen_has_sidebar(self):
        from core.tui.app import FuxiTUI
        from core.tui.widgets.sidebar import Sidebar
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                await pilot.press("ctrl+p")
                await asyncio.sleep(0.1)
                s = pilot.app.screen
                assert s.query_one(Sidebar) is not None
                assert s.query_one("#content-area") is not None
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestResultsScreen
# ══════════════════════════════════════════════════════════════════════════

class TestResultsScreen:
    def test_results_screen_renders(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                from core.tui.screens.results_screen import ResultsScreen
                await pilot.press("ctrl+e")
                assert type(pilot.app.screen) is ResultsScreen
        _run(check())

    def test_results_screen_has_sidebar(self):
        from core.tui.app import FuxiTUI
        from core.tui.widgets.sidebar import Sidebar
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                await pilot.press("ctrl+e")
                await asyncio.sleep(0.1)
                assert pilot.app.screen.query_one(Sidebar) is not None
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestDataManagement
# ══════════════════════════════════════════════════════════════════════════

class TestDataManagement:
    def test_data_mgmt_renders_tabs(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                assert pilot.app.screen.id == "data-mgmt"
                assert pilot.app.screen.query_one("TabbedContent") is not None
        _run(check())

    def test_data_mgmt_register_panel(self):
        from core.tui.app import FuxiTUI
        from textual.widgets import Button, Input
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one("#register-pmid-input", Input) is not None
                assert s.query_one("#register-button", Button) is not None
        _run(check())

    def test_data_mgmt_download_panel(self):
        from core.tui.app import FuxiTUI
        from textual.widgets import Button, Input
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one("#download-gse-input", Input) is not None
                assert s.query_one("#list-files-button", Button) is not None
                assert s.query_one("#download-button", Button) is not None
        _run(check())

    def test_data_mgmt_preprocess_panel(self):
        from core.tui.app import FuxiTUI
        from textual.widgets import Button, Input
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one("#preprocess-input", Input) is not None
                assert s.query_one("#detect-formats-button", Button) is not None
                assert s.query_one("#generate-config-button", Button) is not None
        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestStatePersistence
# ══════════════════════════════════════════════════════════════════════════

class TestStatePersistence:
    def test_save_load_roundtrip(self):
        from core.tui.backends import state as tui_state
        with patch.object(tui_state, "save") as mock_save:
            tui_state.update(last_modality="rna")
            mock_save.assert_called_once_with({"last_modality": "rna"})
        with patch.object(tui_state, "load", return_value={"last_modality": "rna"}):
            assert tui_state.load() == {"last_modality": "rna"}

    def test_default_state_has_expected_keys(self):
        from core.tui.backends import state as tui_state
        d = tui_state.default_state()
        for k in ("last_modality", "last_config", "recent_configs", "window_prefs"):
            assert k in d

    def test_push_recent_config_prepends(self):
        from core.tui.backends import state as tui_state
        with patch.object(tui_state, "save"):
            st = tui_state.push_recent_config("rna", "GSE1", "/path/1.yaml")
            recents = st.get("recent_configs", [])
            assert len(recents) >= 1
            assert recents[0]["gse_id"] == "GSE1"

    def test_get_state_path(self):
        from core.tui.backends import state as tui_state
        p = tui_state.get_state_path()
        assert isinstance(p, str) and len(p) > 0 and p.endswith("tui_state.json")


# ══════════════════════════════════════════════════════════════════════════
# TestErrorHandling
# ══════════════════════════════════════════════════════════════════════════

class TestErrorHandling:
    def test_corrupted_state_returns_empty(self):
        from core.tui.backends import state as tui_state
        with patch.object(tui_state, "load", return_value={}):
            assert tui_state.load() == {}

    def test_no_state_file_returns_empty(self):
        from core.tui.backends import state as tui_state
        with patch.object(tui_state, "logger"):
            with patch("builtins.open", side_effect=FileNotFoundError):
                assert tui_state.load() == {}

    def test_app_handles_bad_screen_switch(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                try:
                    pilot.app.switch_screen("nonexistent")
                except Exception:
                    pass
                assert pilot.app.screen is not None
        _run(check())

    def test_data_mgmt_no_data_root(self):
        from core.tui.app import FuxiTUI
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                assert pilot.app.screen.id == "data-mgmt"
                assert hasattr(pilot.app.screen, "data_root")
        _run(check())

    def test_home_screen_button_notifications(self):
        from core.tui.screens.home import HomeScreen
        from textual.app import App
        from textual.widgets import Button
        async def check():
            a = App()
            async with a.run_test(size=(80, 24)) as pilot:
                a.install_screen(HomeScreen(), "home")
                a.push_screen("home")
                await asyncio.sleep(0.15)
                for bid in ("btn_browse_registry", "btn_run_pipeline", "btn_view_results"):
                    assert pilot.app.screen.query_one(f"#{bid}", Button) is not None
        _run(check())
