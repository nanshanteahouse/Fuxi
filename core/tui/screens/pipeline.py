"""Pipeline runner screen — orchestrates step selection, log display, and progress tracking."""

from __future__ import annotations

import asyncio
import os
from datetime import datetime

from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.reactive import reactive
from textual.screen import Screen
from textual.timer import Timer
from textual.widgets import Button, Footer, Header, Input, Static

from core.pipeline.runner import MODALITY_MAP
from core.tui.backends.pipeline import (
    get_checkpoint_status,
    get_step_dependency,
    run_step,
)
from core.tui.widgets.log_panel import LogPanel
from core.tui.widgets.progress import ProgressTracker
from core.tui.widgets.step_selector import StepSelector


class PipelineRunnerScreen(Screen):
    """Screen for running pipeline steps with step selection, log output, and progress tracking."""

    BINDINGS = [
        ("r", "run_selected_steps", "Run Selected"),
        ("s", "stop_run", "Stop"),
        ("q", "app.pop_screen", "Back"),
        ("ctrl+c", "stop_run", "Stop"),
    ]

    # Reactive state
    modality = reactive("rna", init=False)
    config_path = reactive("", init=False)
    h5ad_dir = reactive("", init=False)
    cell_type = reactive("", init=False)
    annotate_method = reactive("", init=False)
    run_active = reactive(False, init=False)

    def __init__(self, **kwargs) -> None:
        """Initialize the pipeline runner screen."""
        kwargs.setdefault("id", "pipeline")
        super().__init__(**kwargs)
        self._run_task: asyncio.Task[None] | None = None
        self._elapsed_timer: Timer | None = None
        self._start_time: datetime | None = None

    def compose(self) -> ComposeResult:
        """Compose the screen layout."""
        yield Header()
        with Horizontal(id="config-bar"):
            yield Static("Config:", classes="config-label")
            yield Input(
                placeholder="Config file path …",
                id="config_path_input",
            )
            yield Button("Browse", id="btn_browse_config", variant="primary")
        with Horizontal():
            # Left panel: Step selector
            with Vertical(id="selector-panel"):
                yield Static("Step Selection", classes="panel-header")
                yield StepSelector(
                    id="step_selector",
                    modality=self.modality,
                    get_step_dependency=get_step_dependency,
                )
                with Horizontal(id="button-row"):
                    yield Button("Run Selected Steps", id="run_button", variant="primary")
                    yield Button("Stop", id="stop_button", variant="error")

            # Right panel: Log panel and progress tracker
            with Vertical(id="log-panel"):
                yield LogPanel(id="log_panel")
                yield ProgressTracker(id="progress_tracker")
                with Horizontal(id="status-bar"):
                    yield Static("", id="step_count")
                    yield Static("", id="selection_count")

        yield Footer()

    def on_mount(self) -> None:
        """Initialize the screen when mounted."""
        # Sync modality from app (set by HomeScreen)
        if hasattr(self.app, "modality") and self.app.modality:
            self.modality = self.app.modality
        # Auto-populate config path from recent state
        self._auto_populate_config()
        self._load_steps()
        self._update_status_bar()

    def watch_modality(self, old: str, new: str) -> None:
        """Reload steps when modality changes."""
        self._load_steps()

    def watch_config_path(self, old: str, new: str) -> None:
        """Update h5ad_dir and checkpoint status when config path changes."""
        if new and os.path.isfile(new):
            try:
                from core.utils._config import resolve_config

                cfg = resolve_config(new)
                self.h5ad_dir = cfg.h5ad_dir
            except Exception:
                pass
        if new and self.h5ad_dir:
            self._update_checkpoint_status()

    def _auto_populate_config(self) -> None:
        """Auto-populate config_path from recent state or project scanning."""
        import glob

        # Try recent state first
        try:
            from core.tui.backends.state import load

            state = load()
            recent = state.get("recent_configs", [])
            if recent and recent[0].get("modality") == self.modality:
                path = recent[0].get("path", "")
                if path and os.path.isfile(path):
                    self.query_one("#config_path_input", Input).value = path
                    self.config_path = path
                    return
        except Exception:
            pass
        # Scan projects directory
        prefix = f"projects/{self.modality}/"
        if os.path.isdir(prefix):
            configs = sorted(glob.glob(f"{prefix}*/config_*.yaml"))
            if configs:
                # Pick first found
                abs_path = os.path.abspath(configs[0])
                self.query_one("#config_path_input", Input).value = abs_path
                self.config_path = abs_path

    def _browse_config(self) -> None:
        """Handle Browse button — cycle through available configs."""
        import glob

        prefix = f"projects/{self.modality}/"
        configs = sorted(glob.glob(f"{prefix}*/config_*.yaml"))
        if not configs:
            self.app.notify("No configs found", severity="warning")
            return
        current = self.config_path
        try:
            idx = configs.index(os.path.relpath(current)) if current else -1
        except ValueError:
            idx = -1
        next_idx = (idx + 1) % len(configs)
        abs_path = os.path.abspath(configs[next_idx])
        self.query_one("#config_path_input", Input).value = abs_path
        self.config_path = abs_path
        self.app.notify(f"Config: {os.path.basename(abs_path)}")

    def _load_steps(self) -> None:
        """Load step definitions for the current modality."""
        mod = MODALITY_MAP.get(self.modality)
        if not mod:
            return

        selector = self.query_one("#step_selector", StepSelector)
        selector.modality = self.modality
        selector.steps = mod["steps"]
        selector.checkpoints = mod["checkpoints"]

        # Update checkpoint status if we have the data directory
        if self.h5ad_dir:
            self._update_checkpoint_status()

        self._update_status_bar()

    def _update_checkpoint_status(self) -> None:
        """Update the step selector with current checkpoint status."""
        if not self.h5ad_dir:
            return

        try:
            statuses = get_checkpoint_status(self.modality, self.h5ad_dir)
            selector = self.query_one("#step_selector", StepSelector)
            selector.update_checkpoint_status(statuses)
        except Exception as e:
            log_panel = self.query_one("#log_panel", LogPanel)
            log_panel.write_line(f"Warning: Could not load checkpoint status: {e}", "stderr")

    def _update_status_bar(self) -> None:
        """Update the status bar with step and selection counts."""
        selector = self.query_one("#step_selector", StepSelector)
        total_steps = len(selector.steps)
        selected_count = len(selector.selected_indices)

        self.query_one("#step_count", Static).update(f"Total: {total_steps}")
        self.query_one("#selection_count", Static).update(f"Selected: {selected_count}")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        if event.button.id == "run_button":
            self.action_run_selected_steps()
        elif event.button.id == "stop_button":
            self.action_stop_run()
        elif event.button.id == "btn_browse_config":
            self._browse_config()

    def watch_run_active(self, old: bool, new: bool) -> None:
        """Update UI state when running status changes."""
        run_button = self.query_one("#run_button", Button)
        stop_button = self.query_one("#stop_button", Button)

        run_button.disabled = new
        stop_button.disabled = not new

        # Manage elapsed time timer
        if new and self._start_time is None:
            self._start_time = datetime.now()
            self._elapsed_timer = self.set_interval(1.0, self._update_elapsed)
        elif not new and self._elapsed_timer is not None:
            self._elapsed_timer.stop()
            self._elapsed_timer = None

    def _update_elapsed(self) -> None:
        """Update elapsed time display."""
        if self._start_time is None:
            return

        elapsed = (datetime.now() - self._start_time).total_seconds()
        tracker = self.query_one("#progress_tracker", ProgressTracker)
        tracker.set_elapsed(elapsed)

    def action_run_selected_steps(self) -> None:
        """Run all selected steps sequentially."""
        if self.run_active:
            return

        selector = self.query_one("#step_selector", StepSelector)
        selected_indices = selector.selected_indices

        if not selected_indices:
            self.notify("No steps selected. Please select at least one step.", severity="warning")
            return

        log_panel = self.query_one("#log_panel", LogPanel)
        log_panel.clear()
        log_panel.write_line(
            f"Starting pipeline run with {len(selected_indices)} selected steps...", "stdout"
        )

        # Configure progress tracker
        tracker = self.query_one("#progress_tracker", ProgressTracker)
        mod = MODALITY_MAP.get(self.modality)
        if mod:
            labels = [f"{step[0]} {step[1]}" for step in mod["steps"]]
            tracker.set_total_steps(len(mod["steps"]), labels)

        # Reset all step statuses
        if mod:
            for idx in range(len(mod["steps"])):
                tracker.set_step_status(idx, "pending" if idx in selected_indices else "skipped")

        # Start the run task
        self._run_task = asyncio.create_task(self._run_steps(selected_indices))

    async def _run_steps(self, step_indices: list[int]) -> None:
        """Run the specified steps sequentially.

        Parameters
        ----------
        step_indices
            List of step indices to run, in order.
        """
        self.run_active = True
        self._start_time = datetime.now()

        mod = MODALITY_MAP.get(self.modality)
        if not mod:
            self.run_active = False
            return

        log_panel = self.query_one("#log_panel", LogPanel)
        tracker = self.query_one("#progress_tracker", ProgressTracker)

        try:
            for step_idx in step_indices:
                # Check if we were cancelled
                if self._run_task and self._run_task.cancelled():
                    log_panel.write_line("\n\nPipeline run cancelled by user.", "stderr")
                    break

                step_num, script, desc = mod["steps"][step_idx]
                step_name = f"{step_num}_{script}"

                log_panel.set_step_info(self.modality, step_name)
                log_panel.write_line(f"\n{'=' * 60}", "stdout")
                log_panel.write_line(
                    f"Running step {step_idx + 1}/{len(step_indices)}: {step_name}", "stdout"
                )
                log_panel.write_line(f"Description: {desc}", "stdout")
                log_panel.write_line(f"{'=' * 60}\n", "stdout")

                # Update progress tracker
                tracker.set_step_status(step_idx, "running")

                # Run the step
                try:
                    step_succeeded = await self._run_single_step(step_idx)

                    if step_succeeded:
                        tracker.set_step_status(step_idx, "completed")
                        log_panel.write_line(
                            f"\n✓ Step {step_name} completed successfully.", "stdout"
                        )
                    else:
                        tracker.set_step_status(step_idx, "failed")
                        log_panel.write_line(f"\n✗ Step {step_name} failed.", "stderr")

                        # Ask user if they want to continue
                        # For now, we'll continue to the next step
                        log_panel.write_line("Continuing to next step...", "stdout")

                except asyncio.CancelledError:
                    tracker.set_step_status(step_idx, "failed")
                    log_panel.write_line(f"\n✗ Step {step_name} was cancelled.", "stderr")
                    raise
                except Exception as e:
                    tracker.set_step_status(step_idx, "failed")
                    log_panel.write_line(
                        f"\n✗ Step {step_name} encountered an error: {e}", "stderr"
                    )
                    continue

            # All steps done
            log_panel.write_line("\n" + "=" * 60, "stdout")
            log_panel.write_line("Pipeline run completed.", "stdout")
            log_panel.write_line("=" * 60, "stdout")

        except asyncio.CancelledError:
            log_panel.write_line("\n" + "=" * 60, "stderr")
            log_panel.write_line("Pipeline run cancelled.", "stderr")
            log_panel.write_line("=" * 60, "stderr")
        except Exception as e:
            log_panel.write_line(f"\n{'=' * 60}", "stderr")
            log_panel.write_line(f"Pipeline run encountered an unexpected error: {e}", "stderr")
            log_panel.write_line("=" * 60, "stderr")
        finally:
            self.run_active = False
            self._start_time = None
            self._run_task = None

            # Update checkpoint status after run
            self._update_checkpoint_status()

    async def _run_single_step(self, step_idx: int) -> bool:
        """Run a single step and stream its output to the log panel.

        Parameters
        ----------
        step_idx
            Index of the step to run.

        Returns
        -------
        bool
            True if the step succeeded (exit code 0), False otherwise.
        """
        log_panel = self.query_one("#log_panel", LogPanel)

        try:
            async for item in run_step(
                step_idx=step_idx,
                modality=self.modality,
                config_path=self.config_path,
                cell_type=self.cell_type or None,
                annotate_method=self.annotate_method or None,
            ):
                if item["type"] in ("stdout", "stderr"):
                    log_panel.write_line(item["data"], item["type"])
                elif item["type"] == "exit":
                    return item["data"] == 0

        except asyncio.CancelledError:
            raise
        except Exception as e:
            log_panel.write_line(f"Error running step: {e}", "stderr")
            return False

        return False

    def action_stop_run(self) -> None:
        """Stop the currently running pipeline."""
        if not self.run_active or self._run_task is None:
            return

        log_panel = self.query_one("#log_panel", LogPanel)
        log_panel.write_line("\n\nStopping pipeline run...", "stderr")

        self._run_task.cancel()

    # Panel styling
    DEFAULT_CSS = """
    PipelineRunnerScreen {
        layout: vertical;
    }

    #config-bar {
        height: 3;
        padding: 0 1;
        background: $bg-medium;
        border-bottom: solid $border;
    }

    #config-bar > .config-label {
        width: 8;
        content-align: left middle;
        color: $text-secondary;
        text-style: bold;
    }

    #config-bar > #config_path_input {
        width: 1fr;
        margin: 0 1;
    }

    #selector-panel {
        width: 40%;
        height: 1fr;
        overflow: hidden;
        border-right: solid $border;
    }
    #log-panel {
        width: 60%;
        height: 1fr;
        overflow: hidden;
    }
    .panel-header {
        height: 2;
        content-align: left middle;
        padding: 0 1;
        background: $bg-medium;
        border-bottom: solid $border;
        text-style: bold;
        color: $text;
    }

    #button-row {
        height: 3;
        margin-top: 1;
        padding: 0 1;
    }

    #button-row > Button {
        width: 1fr;
        margin-right: 1;
    }

    #button-row > Button:last-child {
        margin-right: 0;
    }

    #status-bar {
        height: 2;
        border-top: solid $border;
        background: $bg-medium;
        padding: 0 1;
    }

    #status-bar > Static {
        content-align: left middle;
        margin-right: 2;
        color: $text-secondary;
        text-style: italic;
    }

    #step-count {
        width: 15;
    }

    #selection-count {
        width: 15;
    }

    StepSelector {
        height: 1fr;
    }

    LogPanel {
        height: 1fr;
    }

    ProgressTracker {
        height: 4;
    }
    """
