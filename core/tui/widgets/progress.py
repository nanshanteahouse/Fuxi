"""Progress tracker — compact pipeline progress widget with per-step status badges.

Shows an overall progress percentage bar, an elapsed-time counter,
and a two-column grid of per-step status badges with animated running indicator.

Typical usage::

    tracker = ProgressTracker()
    tracker.set_total_steps(5, ["QC", "Normalize", "HVG", "PCA", "Cluster"])
    tracker.set_step_status(0, "completed")
    tracker.set_step_status(1, "running")
    tracker.set_elapsed(42.0)
    print(tracker.estimated_remaining())   # "01:15"
"""

from __future__ import annotations

from rich.text import Text
from textual.app import ComposeResult
from textual.reactive import reactive
from textual.timer import Timer
from textual.widget import Widget
from textual.widgets import ProgressBar, Static

# ── Status → display mapping ────────────────────────────────────────────
# Hex values match the Fuxi TUI dark-theme variables from theme.css.

_STATUS_COLORS: dict[str, str] = {
    "pending": "#555555",
    "running": "#e94560",  # $accent
    "completed": "#4ecca3",  # $success
    "failed": "#e74c3c",  # $error
    "skipped": "#ffd93d",  # $warning
}

_STATUS_SYMBOLS: dict[str, str] = {
    "pending": "○",
    "running": "●",
    "completed": "✓",
    "failed": "✗",
    "skipped": "△",
}

_ANIMATION_SYMBOLS: tuple[str, ...] = ("●", "◉", "◎", "◌")


class ProgressTracker(Widget):
    """Compact pipeline progress tracker.

    Composed of:
    - A :class:`~textual.widgets.ProgressBar` showing overall completion %
    - A ``Static`` elapsed-time label
    - A ``Static`` grid of per-step status badges (two per row)

    Each status badge displays the step number, a coloured symbol, and the
    step label.  When a step is ``"running"`` the symbol cycles through an
    animation sequence at 0.8 s intervals.
    """

    DEFAULT_CSS = """
    ProgressTracker {
        height: auto;
        padding: 0 1;
        margin: 0 0 1 0;
    }

    ProgressTracker > #progress {
        width: 100%;
        margin-bottom: 0;
    }

    ProgressTracker > #elapsed {
        color: $text-secondary;
        text-style: italic;
        margin-bottom: 1;
        height: 1;
    }

    ProgressTracker > #step-grid {
        height: auto;
        margin-top: 0;
    }
    """

    total_steps: reactive[int] = reactive(0)

    # ── Private state ───────────────────────────────────────────────────

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._labels: list[str] = []
        self._statuses: list[str] = []
        self._elapsed: float = 0.0
        self._pulse: int = 0
        self._anim_timer: Timer | None = None

    # ── Public API ──────────────────────────────────────────────────────

    def set_total_steps(self, total: int, labels: list[str]) -> None:
        """Configure the total number of steps and their display labels.

        Parameters
        ----------
        total
            Number of pipeline steps.
        labels
            Display label for each step (should be at least *total* long).
        """
        self.total_steps = total
        self._labels = list(labels)
        self._statuses = ["pending"] * total
        self._sync_display()

    def set_step_status(self, index: int, status: str) -> None:
        """Update the status badge for the step at *index*.

        Accepted *status* values:
            ``"pending"``, ``"running"``, ``"completed"``,
            ``"failed"``, ``"skipped"``.
        """
        if 0 <= index < len(self._statuses):
            self._statuses[index] = status
            self._sync_display()

    def set_elapsed(self, seconds: float) -> None:
        """Update the elapsed-time display."""
        self._elapsed = seconds
        self._sync_elapsed()

    def estimated_remaining(self) -> str:
        """Return a rough ``MM:SS`` ETA based on elapsed / completed ratio.

        Returns ``"N/A"`` when there are no completed steps or when
        elapsed time is zero or negative.
        """
        if self.total_steps == 0:
            return "N/A"
        completed = sum(1 for s in self._statuses if s == "completed")
        if completed == 0 or self._elapsed <= 0:
            return "N/A"
        rate = self._elapsed / completed
        remaining = rate * (self.total_steps - completed)
        if remaining < 0:
            return "N/A"
        return self._format_time(int(remaining))

    # ── Textual lifecycle ───────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        yield ProgressBar(total=100, show_percentage=True, id="progress")
        yield Static("", id="elapsed")
        yield Static("", id="step-grid")

    def on_mount(self) -> None:
        self._sync_display()

    # ── Internal helpers ────────────────────────────────────────────────

    def _sync_display(self) -> None:
        """Refresh all child widgets."""
        self._sync_progress()
        self._sync_elapsed()
        self._sync_grid()
        self._sync_anim_timer()

    def _sync_progress(self) -> None:
        bar = self.query_one("#progress", ProgressBar)
        if self.total_steps == 0:
            bar.progress = 0.0
            return
        completed = sum(1 for s in self._statuses if s == "completed")
        bar.progress = (completed / self.total_steps) * 100.0

    def _sync_elapsed(self) -> None:
        self.query_one("#elapsed", Static).update(
            f"Elapsed: {self._format_time(int(self._elapsed))}"
        )

    def _sync_grid(self) -> None:
        grid = self.query_one("#step-grid", Static)
        if self.total_steps == 0:
            grid.update("")
            return

        lines: list[Text] = []
        for i in range(0, self.total_steps, 2):
            left = self._badge(i)
            right = self._badge(i + 1) if i + 1 < self.total_steps else None
            line = Text()
            line.append(left)
            if right is not None:
                line.append("  ")
                line.append(right)
            lines.append(line)

        result = Text("\n").join(lines)
        grid.update(result)

    def _badge(self, index: int) -> Text:
        """Build a single step badge as a ``Text`` renderable."""
        status = self._statuses[index]
        label = self._labels[index] if index < len(self._labels) else f"Step {index + 1}"
        color = _STATUS_COLORS.get(status, "#555")

        if status == "running":
            sym = _ANIMATION_SYMBOLS[self._pulse % len(_ANIMATION_SYMBOLS)]
        else:
            sym = _STATUS_SYMBOLS.get(status, "○")

        t = Text()
        t.append(f"{index + 1:02d} ", style="bold")
        t.append(sym, style=color)
        t.append(f" {label}", style=color)
        return t

    def _sync_anim_timer(self) -> None:
        """Start / stop the pulse timer based on whether any step is running."""
        has_running = any(s == "running" for s in self._statuses)
        if has_running and self._anim_timer is None:
            self._anim_timer = self.set_interval(0.8, self._on_anim_tick)
        elif not has_running and self._anim_timer is not None:
            self._anim_timer.stop()
            self._anim_timer = None

    def _on_anim_tick(self) -> None:
        self._pulse = (self._pulse + 1) % 4
        self._sync_grid()

    # ── Formatting ──────────────────────────────────────────────────────

    @staticmethod
    def _format_time(seconds: int) -> str:
        m, s = divmod(max(0, seconds), 60)
        h, m = divmod(m, 60)
        if h > 0:
            return f"{h}:{m:02d}:{s:02d}"
        return f"{m:02d}:{s:02d}"
