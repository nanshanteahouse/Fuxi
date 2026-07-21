"""ResultsSummaryScreen — display pipeline output files in tabbed panels.

This screen presents parsed pipeline outputs (QC reports, marker genes,
enrichment results) through tabbed panels with configurable layouts.
"""

from __future__ import annotations

from typing import Any

from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.screen import Screen
from textual.widgets import Button, DataTable, Input, Static, TabbedContent, TabPane

from core.tui.backends.results import (
    list_available_reports,
    parse_enrichment,
    parse_marker_genes,
    parse_qc_report,
)
from core.tui.widgets.config_selector import ConfigSelected, ConfigSelector


class ResultsSummaryScreen(Screen):
    """Screen for browsing parsed pipeline outputs.

    Displays QC reports, marker genes, enrichment results, and a summary
    overview in tabbed panels with a config selector on the left.
    """

    id = "results"

    SCREEN_ID = "results"
    TITLE = "Results Summary"

    DEFAULT_CSS = """
    ResultsSummaryScreen {
        height: 100%;
        layout: horizontal;
    }

    ResultsSummaryScreen > Horizontal {
        height: 100%;
    }

    ResultsSummaryScreen > Horizontal > #config-panel {
        width: 33%;
        height: 100%;
        dock: left;
    }

    ResultsSummaryScreen > Horizontal > #results-panel {
        width: 67%;
        height: 100%;
        layout: vertical;
        padding: 1;
    }

    ResultsSummaryScreen > Horizontal > #results-panel > #toolbar {
        height: 3;
        dock: top;
        padding: 0 1;
        background: $bg-medium;
        border-bottom: solid $border;
        content-align: left middle;
    }

    ResultsSummaryScreen > Horizontal > #results-panel > TabbedContent {
        height: 1fr;
    }

    ResultsSummaryScreen > Horizontal > #results-panel > DataTable {
        height: 1fr;
    }

    #empty-message {
        height: 1fr;
        content-align: center middle;
        color: $text-muted;
        text-style: italic;
    }

    #search-input {
        margin: 0 1;
        width: 1fr;
    }

    #reload-btn {
        margin: 0 1 0 0;
    }
    """

    # ── State ──────────────────────────────────────────────────────────────

    _current_config_path: str | None = None
    _qc_data: dict[str, Any] | None = None
    _marker_data: list[dict[str, Any]] | None = None
    _enrichment_data: list[dict[str, Any]] | None = None
    _available_reports: list[str] | None = None

    # ── Compose ─────────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        """Compose the screen with config selector and tabbed results panels."""
        with Horizontal():
            with Vertical(id="config-panel"):
                yield ConfigSelector()

            with Vertical(id="results-panel"):
                with Horizontal(id="toolbar"):
                    yield Static("Results for: ", id="config-label")
                    yield Button("Reload", id="reload-btn", variant="primary")

                with TabbedContent(initial="qc"):
                    with TabPane("QC Report", id="qc"):
                        yield DataTable(id="qc-table")
                        yield Static("No QC report available", id="qc-empty", classes="hidden")

                    with TabPane("Marker Genes", id="markers"):
                        with Horizontal(id="marker-toolbar"):
                            yield Input(placeholder="Search markers...", id="search-input")
                        yield DataTable(id="marker-table")
                        yield Static(
                            "No marker genes available", id="marker-empty", classes="hidden"
                        )

                    with TabPane("Enrichment", id="enrichment"):
                        yield DataTable(id="enrichment-table")
                        yield Static(
                            "No enrichment results available",
                            id="enrichment-empty",
                            classes="hidden",
                        )

                    with TabPane("Summary", id="summary"):
                        yield Static(id="summary-list")

    # ── Mount and refresh ────────────────────────────────────────────────────

    def on_mount(self) -> None:
        """Initialize with empty state."""
        self._update_config_label()

    # ── Event handlers ───────────────────────────────────────────────────────

    def on_config_selected(self, event: ConfigSelected) -> None:
        """Handle config selection and refresh results."""
        self._current_config_path = event.path
        self._update_config_label()
        self._load_results()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle reload button press."""
        if event.button.id == "reload-btn" and self._current_config_path:
            self._load_results()

    def on_input_changed(self, event: Input.Changed) -> None:
        """Handle search input for marker genes."""
        if event.input.id == "search-input" and self._marker_data:
            self._filter_marker_table(event.value)

    # ── Data loading ────────────────────────────────────────────────────────

    def _load_results(self) -> None:
        """Load all results for the current config."""
        if not self._current_config_path:
            self._clear_all_tables()
            return

        try:
            self._qc_data = parse_qc_report(self._current_config_path)
            self._marker_data = parse_marker_genes(self._current_config_path)
            self._enrichment_data = parse_enrichment(self._current_config_path)
            self._available_reports = list_available_reports(self._current_config_path)

            self._update_qc_table()
            self._update_marker_table()
            self._update_enrichment_table()
            self._update_summary()

        except Exception as e:
            self._show_error(f"Error loading results: {e}")

    def _update_qc_table(self) -> None:
        """Populate the QC report table."""
        table = self.query_one("#qc-table", DataTable)
        empty = self.query_one("#qc-empty", Static)

        if not self._qc_data:
            table.clear()
            table.display = False
            empty.display = True
            empty.update("No QC report available")
            return

        table.display = True
        empty.display = False
        table.clear()

        table.add_column("Metric", width=30)
        table.add_column("Value", width=50)

        if isinstance(self._qc_data, dict):
            # Single-row QC report
            for metric, value in self._qc_data.items():
                table.add_row(str(metric), str(value))
        elif isinstance(self._qc_data, dict) and "summary" in self._qc_data:
            # Multi-row QC report
            for row in self._qc_data["summary"]:
                for metric, value in row.items():
                    table.add_row(str(metric), str(value))

    def _update_marker_table(self) -> None:
        """Populate the marker genes table."""
        table = self.query_one("#marker-table", DataTable)
        empty = self.query_one("#marker-empty", Static)

        if not self._marker_data:
            table.clear()
            table.display = False
            empty.display = True
            empty.update("No marker genes available")
            return

        table.display = True
        empty.display = False
        table.clear()

        # Get columns from first row
        if self._marker_data:
            columns = list(self._marker_data[0].keys())
            for col in columns:
                table.add_column(str(col), width=20 if col != "names" else 40)

            for row in self._marker_data:
                table.add_row(*[str(row.get(col, "")) for col in columns])

    def _filter_marker_table(self, search_text: str) -> None:
        """Filter the marker genes table by search text."""
        table = self.query_one("#marker-table", DataTable)
        if not search_text or not self._marker_data:
            self._update_marker_table()
            return

        table.clear()
        search_lower = search_text.lower()

        if self._marker_data:
            columns = list(self._marker_data[0].keys())
            for col in columns:
                table.add_column(str(col), width=20 if col != "names" else 40)

            for row in self._marker_data:
                row_text = " ".join(str(v).lower() for v in row.values())
                if search_lower in row_text:
                    table.add_row(*[str(row.get(col, "")) for col in columns])

    def _update_enrichment_table(self) -> None:
        """Populate the enrichment results table, sorted by p-value."""
        table = self.query_one("#enrichment-table", DataTable)
        empty = self.query_one("#enrichment-empty", Static)

        if not self._enrichment_data:
            table.clear()
            table.display = False
            empty.display = True
            empty.update("No enrichment results available")
            return

        table.display = True
        empty.display = False
        table.clear()

        # Sort by p-value if available
        sorted_data = sorted(
            self._enrichment_data,
            key=lambda x: float(x.get("pval", float("inf")) or float("inf")),
        )

        # Get columns from first row
        if sorted_data:
            columns = list(sorted_data[0].keys())
            for col in columns:
                table.add_column(str(col), width=15)

            for row in sorted_data:
                table.add_row(*[str(row.get(col, "")) for col in columns])

    def _update_summary(self) -> None:
        """Update the summary panel with available reports."""
        summary_list = self.query_one("#summary-list", Static)

        if not self._available_reports:
            summary_list.update("No reports available for this configuration.")
            return

        lines = ["[bold]Available Reports:[/]\n"]
        for report in self._available_reports:
            lines.append(f"  • {report}")

        summary_list.update("\n".join(lines))

    # ── Helper methods ─────────────────────────────────────────────────────

    def _update_config_label(self) -> None:
        """Update the config label with the selected config basename."""
        label = self.query_one("#config-label", Static)
        if self._current_config_path:
            label.update(f"Results for: [bold]{self._current_config_path.split('/')[-1]}[/]")
        else:
            label.update("Results for: [italic]No config selected[/]")

    def _clear_all_tables(self) -> None:
        """Clear all tables and show empty states."""
        for table_id, empty_id, empty_msg in [
            ("qc-table", "qc-empty", "No QC report available"),
            ("marker-table", "marker-empty", "No marker genes available"),
            ("enrichment-table", "enrichment-empty", "No enrichment results available"),
        ]:
            table = self.query_one(f"#{table_id}", DataTable)
            empty = self.query_one(f"#{empty_id}", Static)
            table.clear()
            table.display = False
            empty.display = True
            empty.update(empty_msg)

        self.query_one("#summary-list", Static).update("No reports available.")

    def _show_error(self, message: str) -> None:
        """Display an error message in the current tab."""
        # Find the active tab and show error there
        try:
            self.query_one(TabbedContent)
            Static(f"[red]{message}[/]")
            # For now, just log the error - could be enhanced to show in each pane
            self.log(message)
        except Exception:
            self.log(message)
