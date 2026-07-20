"""TUI screens — full-page views for the Fuxi pipeline."""

from core.tui.screens.config_editor import ConfigEditorScreen
from core.tui.screens.data_mgmt import DataManagementScreen
from core.tui.screens.home import HomeScreen
from core.tui.screens.pipeline import PipelineRunnerScreen
from core.tui.screens.registry import RegistryBrowserScreen
from core.tui.screens.results import ResultsSummaryScreen

__all__ = [
    "ConfigEditorScreen",
    "DataManagementScreen",
    "HomeScreen",
    "PipelineRunnerScreen",
    "RegistryBrowserScreen",
    "ResultsSummaryScreen",
]
