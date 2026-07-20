"""Log panel widget — streaming log display with ANSI color support."""

from __future__ import annotations

from datetime import datetime
from textual.app import ComposeResult
from textual.reactive import reactive
from textual.widgets import RichLog, Static
from textual.widget import Widget


class LogPanel(Widget):
    """Streaming log panel widget for pipeline subprocess output.

    Displays log lines with ANSI color support, auto-scroll, and stream source
    indicators. Each line is prefixed with a timestamp and a stream indicator
    (→ for stdout, ✗ for stderr).
    """

    DEFAULT_CSS = """
    LogPanel {
        height: 100%;
        layout: vertical;
        background: $bg-dark;
    }

    LogPanel > Static#header {
        height: 3;
        content-align: left middle;
        padding: 0 1;
        background: $bg-medium;
        border-bottom: solid $border;
        text-style: bold;
        color: $text;
    }

    LogPanel > RichLog {
        height: 1fr;
        background: $bg-dark;
        border: none;
        scrollbar-background: $bg-dark;
        scrollbar-color: $accent $bg-dark;
    }
    """

    auto_scroll_enabled = reactive(True)

    def compose(self) -> ComposeResult:
        """Compose the widget with header and RichLog."""
        yield Static("Idle", id="header")
        yield RichLog(id="log", wrap=True, markup=True, auto_scroll=True)

    def on_mount(self) -> None:
        """Initialize widget on mount."""
        self._rich_log = self.query_one("#log", RichLog)
        self._header = self.query_one("#header", Static)
        self._rich_log.scroll_end()

    def set_step_info(self, modality: str, step_name: str) -> None:
        """Update the header with current modality and step name.

        Parameters
        ----------
        modality
            The analysis modality (e.g., "rna", "atac", "spatial").
        step_name
            The name of the current step being executed.
        """
        self._header.update(f"Running [{modality}] step {step_name}")

    def write_line(self, line: str, stream: str) -> None:
        """Append a log line with timestamp and stream indicator.

        Parameters
        ----------
        line
            The log line content. May contain ANSI escape sequences.
        stream
            The stream source, either "stdout" or "stderr".
        """
        if not line:
            # Handle empty lines as line breaks
            self._rich_log.write("")
            return

        timestamp = datetime.now().strftime("%H:%M:%S")

        if stream == "stderr":
            # Stderr lines shown in red/dim
            indicator = "✗"
            styled_line = f"[dim red][{timestamp}] {indicator}[/] {line}"
        else:
            # Stdout lines with normal styling
            indicator = "→"
            styled_line = f"[{timestamp}] {indicator} {line}"

        self._rich_log.write(styled_line)

        # Auto-scroll if enabled
        if self.auto_scroll_enabled:
            self._rich_log.scroll_end(animate=False)

    def clear(self) -> None:
        """Reset all log content and header to idle state."""
        self._rich_log.clear()
        self._header.update("Idle")

    def watch_auto_scroll_enabled(self, old_value: bool, new_value: bool) -> None:
        """React to auto-scroll toggle changes.

        Parameters
        ----------
        old_value
            Previous auto-scroll state.
        new_value
            New auto-scroll state.
        """
        if new_value:
            self._rich_log.scroll_end(animate=False)

    def toggle_auto_scroll(self) -> None:
        """Toggle auto-scroll on/off."""
        self.auto_scroll_enabled = not self.auto_scroll_enabled