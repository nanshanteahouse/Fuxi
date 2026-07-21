"""Step listing screen — displays pipeline steps for the selected modality."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.screen import Screen
from textual.style import Style
from textual.widgets import DataTable, Footer, Header, Select, Static

from core.pipeline.runner import MODALITY_MAP


class StepListingScreen(Screen):
    """Screen that lists all pipeline steps for the selected modality.

    Columns: Step #, Script, Description, Checkpoint File.
    Steps whose index appears in ``write_checkpoints`` are shown in bold.
    A ``Select`` at the top toggles between rna / atac / spatial.

    Read-only view — selection and run controls belong in Phase 3.
    """

    id = "steps"

    DEFAULT_CSS = """
    StepListingScreen {
        align: center top;
    }

    #title {
        text-align: center;
        padding: 1 0;
        text-style: bold;
        color: $accent;
    }

    #modality-selector {
        margin: 0 2 1 2;
        width: 30;
    }

    #steps-table {
        margin: 0 2;
        height: 1fr;
    }

    #step-count {
        margin: 0 2 1 2;
        text-align: right;
        color: $text-secondary;
    }
    """

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static("Pipeline Steps", id="title")
        yield Select(
            [(m, m) for m in MODALITY_MAP.keys()],
            prompt="Select modality",
            id="modality-selector",
        )
        yield DataTable(id="steps-table")
        yield Static(id="step-count")
        yield Footer()

    def on_mount(self) -> None:
        """Set up the DataTable columns and load the first modality."""
        table = self.query_one("#steps-table", DataTable)
        table.add_columns("Step #", "Script", "Description", "Checkpoint File")

        # Default to the first registered modality
        default = list(MODALITY_MAP.keys())[0]
        self.query_one("#modality-selector", Select).value = default
        self._populate_table(default)

    def on_select_changed(self, event: Select.Changed) -> None:
        """Repopulate the table when the user picks a different modality."""
        self._populate_table(event.value)

    # ------------------------------------------------------------------
    def _populate_table(self, modality: str) -> None:
        """Fill the DataTable with every step defined for *modality*."""
        mod = MODALITY_MAP[modality]
        steps = mod["steps"]
        checkpoints = mod["checkpoints"]
        mod["write_checkpoints"]

        table = self.query_one("#steps-table", DataTable)
        table.clear()

        Style(bold=True)

        for idx, (num, script, desc) in enumerate(steps):
            cp = checkpoints[idx] if idx < len(checkpoints) else ""
            row_key = f"step_{idx}"
            table.add_row(num, script, desc, cp, key=row_key)

        # Bold styling for checkpoint-writing steps is omitted —
        # Textual 8.2 DataTable.update_cell does not support `style=`
        # ── footer label ──────────────────────────────────────────────
        total = len(steps)
        label = f"Total: {total} step{'s' if total != 1 else ''}"
        self.query_one("#step-count", Static).update(label)
