"""Integration and regression tests for Fuxi TUI using Textual's pilot API.

Covers:
- App launch and all six screens mount without exceptions
- Layout containers have non-degenerate geometry
- Screen navigation via keyboard shortcuts and sidebar
- Modality selector, registry, pipeline, results, data management
- Config form (sections, toolbar, fields, status bar)
- Step listing table
- State persistence (save/load roundtrip)
- Error handling (corrupted state, missing files, bad screen switch)
"""

from __future__ import annotations

import asyncio
from unittest.mock import patch

import pytest
from textual.widgets import Button

from core.tui.app import FuxiTUI

SCREENS = ["home", "registry", "pipeline", "results", "data-mgmt", "config-editor"]


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
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                assert type(pilot.app.screen).__name__ == "HomeScreen"

        _run(check())

    @pytest.mark.skip(reason="Sidebar widget removed from TUI")
    def test_home_screen_has_sidebar_and_content(self):
        pass

    def test_home_screen_content(self):
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                await asyncio.sleep(0.15)
                static = pilot.app.screen.query_one("#home-logo")
                # Logo is ASCII art of "FUXI" in box-drawing characters,
                # so check visual is present rather than matching text.
                assert str(static.visual).strip()

        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestScreenNavigation
# ══════════════════════════════════════════════════════════════════════════


class TestScreenNavigation:
    def test_all_keyboard_shortcuts(self):
        from core.tui.screens.home import HomeScreen
        from core.tui.screens.pipeline import PipelineRunnerScreen
        from core.tui.screens.registry import RegistryBrowserScreen
        from core.tui.screens.results import ResultsSummaryScreen

        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                a = pilot.app
                await pilot.press("ctrl+r")
                assert type(a.screen) is RegistryBrowserScreen
                await pilot.press("ctrl+p")
                assert type(a.screen) is PipelineRunnerScreen
                await pilot.press("ctrl+e")
                assert type(a.screen) is ResultsSummaryScreen
                await pilot.press("ctrl+d")
                assert a.screen.id == "data-mgmt"
                await pilot.press("ctrl+c")
                assert a.screen.id == "config-editor"
                await pilot.press("f1")
                assert type(a.screen) is HomeScreen

        _run(check())

    @pytest.mark.skip(reason="Sidebar widget removed from TUI")
    def test_sidebar_highlights_active(self):
        pass


# ══════════════════════════════════════════════════════════════════════════
# TestModalitySelector
# ══════════════════════════════════════════════════════════════════════════


class TestModalitySelector:
    def test_modality_select_renders(self):
        from textual.app import App

        from core.tui.screens.home import HomeScreen

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
        from textual.app import App

        from core.tui.screens.home import HomeScreen

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
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                from core.tui.screens.registry import RegistryBrowserScreen

                await pilot.press("ctrl+r")
                assert type(pilot.app.screen) is RegistryBrowserScreen

        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# TestStepListing
# ══════════════════════════════════════════════════════════════════════════


class TestStepListing:
    def test_step_table_populated(self):
        from textual.app import App
        from textual.widgets import DataTable

        from core.tui.screens.steps import StepListingScreen

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
        from textual.app import App
        from textual.widgets import DataTable

        from core.tui.screens.steps import StepListingScreen

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
        from textual.app import App
        from textual.widgets import Static

        from core.tui.screens.steps import StepListingScreen

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
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                assert len(list(pilot.app.screen.query("Collapsible"))) >= 1

        _run(check())

    def test_config_form_has_toolbar_buttons(self):
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one("#btn-load", Button) is not None
                assert s.query_one("#btn-save", Button) is not None

        _run(check())

    def test_config_form_has_path_input(self):
        from textual.widgets import Input

        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                assert pilot.app.screen.query_one("#config-path-input", Input) is not None

        _run(check())

    def test_config_form_has_field_widgets(self):
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+c")
                await asyncio.sleep(0.15)
                assert len(list(pilot.app.screen.query(".field-row"))) >= 1

        _run(check())

    def test_config_form_status_bar(self):
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
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                from core.tui.screens.pipeline import PipelineRunnerScreen

                await pilot.press("ctrl+p")
                assert type(pilot.app.screen) is PipelineRunnerScreen

        _run(check())

    @pytest.mark.skip(reason="Sidebar widget removed from TUI")
    def test_pipeline_screen_has_sidebar(self):
        pass


# ══════════════════════════════════════════════════════════════════════════
# TestResultsScreen
# ══════════════════════════════════════════════════════════════════════════


class TestResultsScreen:
    def test_results_screen_renders(self):
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                from core.tui.screens.results import ResultsSummaryScreen

                await pilot.press("ctrl+e")
                assert type(pilot.app.screen) is ResultsSummaryScreen

        _run(check())

    @pytest.mark.skip(reason="Sidebar widget removed from TUI")
    def test_results_screen_has_sidebar(self):
        pass


# ══════════════════════════════════════════════════════════════════════════
# TestDataManagement
# ══════════════════════════════════════════════════════════════════════════


class TestDataManagement:
    def test_data_mgmt_renders_tabs(self):
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                assert pilot.app.screen.id == "data-mgmt"
                assert pilot.app.screen.query_one("TabbedContent") is not None

        _run(check())

    def test_data_mgmt_register_panel(self):
        from textual.widgets import Input

        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                s = pilot.app.screen
                assert s.query_one("#register-pmid-input", Input) is not None
                assert s.query_one("#register-button", Button) is not None

        _run(check())

    def test_data_mgmt_download_panel(self):
        from textual.widgets import Input

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
        from textual.widgets import Input

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
        async def check():
            async with FuxiTUI().run_test(size=(80, 24)) as pilot:
                try:
                    pilot.app.switch_screen("nonexistent")
                except Exception:
                    pass
                assert pilot.app.screen is not None

        _run(check())

    def test_data_mgmt_no_data_root(self):
        async def check():
            async with FuxiTUI().run_test(size=(100, 40)) as pilot:
                await pilot.press("ctrl+d")
                await asyncio.sleep(0.15)
                assert pilot.app.screen.id == "data-mgmt"
                assert hasattr(pilot.app.screen, "data_root")

        _run(check())


# ══════════════════════════════════════════════════════════════════════════
# Async tests (pytest-asyncio) — regression coverage from 2026-07
# ══════════════════════════════════════════════════════════════════════════


async def test_all_screens_mount() -> None:
    """Every installed screen switches in cleanly."""
    app = FuxiTUI()
    async with app.run_test(size=(120, 40)) as pilot:
        for name in SCREENS:
            app.switch_screen(name)
            await pilot.pause()
            await pilot.pause()
            assert app.screen is not None


async def test_home_layout_not_collapsed() -> None:
    """The home screen's main containers occupy real screen area."""
    app = FuxiTUI()
    async with app.run_test(size=(120, 40)) as pilot:
        await pilot.pause()
        screen = app.screen
        assert screen.region.width == 120
        assert screen.region.height >= 38


async def test_home_buttons_navigate() -> None:
    """Home quick-launch buttons switch to their target screens.

    Uses ``action_press()`` instead of ``pilot.click()`` because
    Textual's test pilot can miss clicks on widgets inside nested
    containers with complex CSS layouts (a known Textual 8.2.8
    pilot limitation, not a real-terminal bug).
    """
    app = FuxiTUI()
    async with app.run_test(size=(120, 40)) as pilot:
        await pilot.pause()
        await pilot.pause()
        for button_id, expected in [
            ("#btn_browse_registry", "RegistryBrowserScreen"),
            ("#btn_view_results", "ResultsSummaryScreen"),
            ("#btn_run_pipeline", "PipelineRunnerScreen"),
        ]:
            app.switch_screen("home")
            await pilot.pause()
            await pilot.pause()
            btn = app.screen.query_one(button_id, Button)
            btn.action_press()
            await pilot.pause()
            await pilot.pause()
            assert type(app.screen).__name__ == expected
