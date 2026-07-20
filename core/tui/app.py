"""FuxiTUI — Root Textual application for the Fuxi pipeline."""

from __future__ import annotations

from textual.app import App
from textual.binding import Binding

from core.tui.screens.config_editor import ConfigEditorScreen
from core.tui.screens.data_mgmt import DataManagementScreen
from core.tui.screens.home import HomeScreen
from core.tui.screens.pipeline import PipelineRunnerScreen
from core.tui.screens.registry import RegistryBrowserScreen
from core.tui.screens.results import ResultsSummaryScreen


class FuxiTUI(App):
    """Root TUI application for interacting with the Fuxi pipeline."""

    CSS_PATH = "theme.css"

    def get_css_variables(self) -> dict[str, str]:
        """Extend built-in CSS variables with Fuxi custom palette."""
        variables = super().get_css_variables()
        variables.update({
            "bg-dark": "#1a1a2e",
            "bg-medium": "#16213e",
            "bg-light": "#0f3460",
            "text-primary": "#e0e0e0",
            "text-secondary": "#a0a0a0",
            "text-muted": "#6e7681",
            "accent-hover": "#ff6b6b",
            "border": "#2a2a4a",
            "border-light": "#3a3a5a",
            "highlight": "#2a2a4a",
            "scrollbar-track": "#0f0f1a",
            "scrollbar-thumb": "#4a4a6a",
            "info": "#3498db",
            "warning-bg": "#3a3a2e",
            "warning-fg": "#ffd93d",
        })
        return variables

    BINDINGS = [
        Binding("f1", "switch_screen('home')", "Home", priority=True),
        Binding("ctrl+r", "switch_screen('registry')", "Registry", priority=True),
        Binding("ctrl+p", "switch_screen('pipeline')", "Pipeline", priority=True),
        Binding("ctrl+e", "switch_screen('results')", "Results", priority=True),
        Binding("ctrl+d", "switch_screen('data-mgmt')", "Data Mgmt", priority=True),
        Binding("ctrl+c", "switch_screen('config-editor')", "Config", priority=True),
        Binding("ctrl+q", "quit", "Quit", priority=True),
    ]

    def on_mount(self) -> None:
        """Install all screens, then push the home screen.

        Creates fresh screen instances so each app run is isolated
        (important when running multiple apps in the same process / loop).
        """
        screens = {
            "home": HomeScreen(),
            "registry": RegistryBrowserScreen(),
            "pipeline": PipelineRunnerScreen(),
            "results": ResultsSummaryScreen(),
            "data-mgmt": DataManagementScreen(),
            "config-editor": ConfigEditorScreen(),
        }
        for name, screen in screens.items():
            self.install_screen(screen, name)
        self.push_screen("home")
