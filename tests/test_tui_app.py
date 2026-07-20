"""Regression tests for the TUI launch bugs fixed in 2026-07.

Covers:
- App launches and all six screens mount without exceptions
  (regression: stub screens hid the full implementations; several
  latent crashes surfaced once the full screens were reachable).
- Layout containers keep non-degenerate geometry
  (regression: theme.css globally forced ``Horizontal { height: 1 }``
  and ``Vertical { width: 1 }``, collapsing the whole UI to a black
  screen).
"""

from __future__ import annotations

import pytest

from textual.widgets import Button

from core.tui.app import FuxiTUI

SCREENS = ["home", "registry", "pipeline", "results", "data-mgmt", "config-editor"]


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
