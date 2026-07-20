"""Home screen — branded landing page with modality selector and quick-launch buttons."""

from __future__ import annotations

from textual.screen import Screen
from textual.containers import Horizontal, Vertical
from textual.widgets import Button, Header, Select, Footer, Static
from core.pipeline.runner import MODALITY_MAP

LOGO = "\n".join([
    "███████╗██╗   ██╗██╗  ██╗██╗",
    "██╔════╝██║   ██║╚██╗██╔╝██║",
    "█████╗  ██║   ██║ ╚███╔╝ ██║",
    "██╔══╝  ██║   ██║ ██╔██╗ ██║",
    "██║     ╚██████╔╝██╔╝ ██╗██║",
    "╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝",
])

class HomeScreen(Screen):
    """Home screen — branded landing page with modality selector and quick-launch buttons."""

    id = "home"

    DEFAULT_CSS = """
    HomeScreen {
        align: center middle;
    }

    #home-content {
        width: auto;
        max-width: 72;
        min-width: 40;
        height: auto;
        overflow: hidden;
        margin: 1 2;
    }

    #home-logo {
        content-align: center middle;
        color: $accent;
        text-style: bold;
        padding: 1 0;
        margin-bottom: 0;
    }

    #home-subtitle {
        content-align: center middle;
        color: $text-disabled;
        padding: 0;
        margin-bottom: 2;
    }

    #home-divider {
        height: 1;
        background: $accent;
        margin: 0 2 1 2;
    }

    .home-card {
        border: solid $border;
        padding: 1 3;
        margin: 0;
    }

    .home-card > .card-title {
        content-align: center top;
        color: $text-secondary;
        padding: 0;
        margin-bottom: 1;
    }

    .home-card Select {
        margin: 1 0;
        width: 100%;
    }

    .home-card Horizontal {
        width: 100%;
        height: auto;
        align-horizontal: center;
    }

    .home-card Button {
        margin: 0 1;
        min-width: 14;
    }
    """

    def compose(self):
        modalities = list(MODALITY_MAP.keys())

        yield Header()

        with Vertical(id="home-content"):
            yield Static(LOGO, id="home-logo")
            yield Static(
                "scRNA-seq  •  scATAC-seq  •  Spatial Transcriptomics",
                id="home-subtitle",
            )
            yield Static(id="home-divider")

            with Vertical(classes="home-card"):
                yield Static("Select modality to begin", classes="card-title")
                yield Select(
                    [(mod, mod) for mod in modalities],
                    value=modalities[0] if modalities else None,
                    id="modality_select",
                    allow_blank=False,
                )
                with Horizontal():
                    yield Button("Registry", id="btn_browse_registry")
                    yield Button("Run Pipeline", id="btn_run_pipeline")
                    yield Button("Results", id="btn_view_results")

        yield Footer()

    def on_select_changed(self, event: Select.Changed) -> None:
        """Store selected modality in the app."""
        self.app.modality = event.value
        self.app.notify(f"Modality: {event.value}")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Navigate to the screen matching the pressed button."""
        button_id = event.button.id

        if button_id == "btn_browse_registry":
            self.app.switch_screen("registry")
        elif button_id == "btn_run_pipeline":
            self.app.switch_screen("pipeline")
        elif button_id == "btn_view_results":
            self.app.switch_screen("results")
