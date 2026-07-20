"""Home screen — modality selector and quick launch buttons."""

from __future__ import annotations

from textual.screen import Screen
from textual.containers import Horizontal
from textual.widgets import Button, Header, Select, Footer
from core.pipeline.runner import MODALITY_MAP


class HomeScreen(Screen):
    """Home screen with modality selector and quick-launch buttons."""

    id = "home"

    def compose(self):
        """Compose the home screen layout."""
        modalities = list(MODALITY_MAP.keys())

        yield Header()

        yield Select(
            [(mod, mod) for mod in modalities],
            value=modalities[0] if modalities else None,
            id="modality_select",
            prompt="Select modality",
        )

        yield Horizontal(
            Button("Browse Registry", id="btn_browse_registry"),
            Button("Run Pipeline", id="btn_run_pipeline"),
            Button("View Results", id="btn_view_results"),
        )

        yield Footer()

    def on_select_changed(self, event: Select.Changed) -> None:
        """Store selected modality in the app."""
        self.app.modality = event.value
        self.app.notify(f"Modality selected: {event.value}")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Navigate to the screen matching the pressed button."""
        button_id = event.button.id

        if button_id == "btn_browse_registry":
            self.app.switch_screen("registry")
        elif button_id == "btn_run_pipeline":
            self.app.switch_screen("pipeline")
        elif button_id == "btn_view_results":
            self.app.switch_screen("results")