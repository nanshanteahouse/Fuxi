"""Step selector — multi-select widget for pipeline steps with dependency visualization.

Usage:
    selector = StepSelector(modality="rna")
    selector.steps = [("00", "00_load.py", "Load raw data → 00_raw.h5ad"), ...]
    selector.checkpoints = ["00_raw.h5ad", "01_doublet.h5ad", ...]
    selector.update_checkpoint_status([True, False, True, ...])
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from textual.reactive import reactive
from textual.widget import Widget
from textual.widgets import Static

if TYPE_CHECKING:
    from typing import Callable


class StepSelector(Widget):
    """Multi-select widget for pipeline steps with dependency visualization.

    Shows all steps for a modality with checkboxes, highlights already-completed
    steps, and indicates required pre-steps for dependency validation.

    Attributes
    ----------
    modality
        The modality identifier (e.g., "rna", "atac", "spatial").
    steps
        List of tuples (num, script, desc) defining each step.
    checkpoints
        List of checkpoint file names for each step.
    completed_steps
        Set of step indices that have been completed.
    selected_steps
        Set of step indices currently selected by the user.
    """

    modality: reactive[str] = reactive("rna", init=False)
    steps: reactive[list[tuple[str, str, str]]] = reactive([], init=False)
    checkpoints: reactive[list[str]] = reactive([], init=False)
    completed_steps: reactive[set[int]] = reactive(set(), init=False)
    selected_steps: reactive[set[int]] = reactive(set(), init=False)

    # Optional external dependency checker (e.g., from pipeline backend)
    _get_step_dependency: Callable[[int, str], str] | None = None

    DEFAULT_CSS = """
    StepSelector {
        height: auto;
        min-height: 10;
        border: solid $border;
        padding: 1;
        background: $bg-dark;
    }

    StepSelector > .step-row {
        height: auto;
        min-height: 2;
        margin: 0 1;
        padding: 0;
    }

    StepSelector > .step-row:hover {
        background: $highlight;
    }

    StepSelector > .step-row.completed {
        color: $success;
    }

    StepSelector > .step-row.selected {
        background: $bg-light;
        border-left: solid $accent;
    }

    StepSelector > .step-row.selected > .checkbox {
        background: $accent;
        color: white;
    }

    StepSelector > .step-number {
        width: 3;
        text-style: bold;
        text-align: right;
        margin-right: 1;
    }

    StepSelector > .checkbox {
        width: 2;
        text-align: center;
        margin-right: 1;
    }

    StepSelector > .script-name {
        width: 20;
        text-style: bold;
        margin-right: 1;
    }

    StepSelector > .description {
        width: 1fr;
        color: $text-secondary;
        text-style: dim;
    }

    StepSelector > .checkpoint-file {
        width: 20;
        color: $text-muted;
        text-style: italic;
    }

    StepSelector > .dependency-warning {
        color: $warning;
        text-style: bold;
        margin-left: 1;
    }

    StepSelector > .completed-indicator {
        color: $success;
        text-style: bold;
    }
    """

    def __init__(
        self,
        modality: str = "rna",
        steps: list[tuple[str, str, str]] | None = None,
        checkpoints: list[str] | None = None,
        get_step_dependency: Callable[[int, str], str] | None = None,
        **kwargs,
    ) -> None:
        """Initialize the step selector.

        Parameters
        ----------
        modality
            The modality identifier (e.g., "rna", "atac", "spatial").
        steps
            List of tuples (num, script, desc) defining each step.
        checkpoints
            List of checkpoint file names for each step.
        get_step_dependency
            Optional function that returns the checkpoint file a step reads from.
            Signature: (step_idx: int, modality: str) -> str
        """
        super().__init__(**kwargs)
        self.modality = modality
        self.steps = steps or []
        self.checkpoints = checkpoints or []
        self._get_step_dependency = get_step_dependency

    def on_mount(self) -> None:
        """Compose the step list when mounted."""
        self._render_steps()

    def watch_modality(self, old: str, new: str) -> None:
        """Re-render when modality changes."""
        self._render_steps()

    def watch_steps(self, old: list, new: list) -> None:
        """Re-render when steps list changes."""
        self._render_steps()

    def watch_checkpoints(self, old: list, new: list) -> None:
        """Re-render when checkpoints list changes."""
        self._render_steps()

    def watch_completed_steps(self, old: set, new: set) -> None:
        """Re-render when completion status changes."""
        self._render_steps()

    def watch_selected_steps(self, old: set, new: set) -> None:
        """Re-render when selection changes."""
        self._render_steps()

    def _render_steps(self) -> None:
        """Render all step rows with current state."""
        # Clear existing children
        self.remove_children()

        for idx, (num, script, desc) in enumerate(self.steps):
            # Determine row state
            is_completed = idx in self.completed_steps
            is_selected = idx in self.selected_steps

            # Check for missing dependencies
            dep_warning = self._check_dependency(idx)

            # Build row classes
            classes = ["step-row"]
            if is_completed:
                classes.append("completed")
            if is_selected:
                classes.append("selected")

            # Get checkpoint file for this step
            ckpt = self.checkpoints[idx] if idx < len(self.checkpoints) else ""

            # Render checkbox or completed indicator
            if is_completed:
                checkbox_content = "✓"
                checkbox_class = "completed-indicator"
            else:
                checkbox_content = "[ ]" if not is_selected else "[X]"
                checkbox_class = "checkbox"

            # Build the row
            row = Static(
                self._render_row_content(
                    num,
                    script,
                    desc,
                    ckpt,
                    checkbox_content,
                    checkbox_class,
                    dep_warning,
                ),
                classes=" ".join(classes),
            )

            # Store row index as data attribute for click handling
            row._step_idx = idx  # type: ignore[attr-defined]

            self.mount(row)

    def _render_row_content(
        self,
        num: str,
        script: str,
        desc: str,
        ckpt: str,
        checkbox_content: str,
        checkbox_class: str,
        dep_warning: str | None,
    ) -> str:
        """Render the content for a single step row.

        Parameters
        ----------
        num
            Step number (e.g., "00", "01").
        script
            Script filename (e.g., "00_load.py").
        desc
            Step description.
        ckpt
            Checkpoint file name.
        checkbox_content
            Checkbox state indicator.
        checkbox_class
            CSS class for the checkbox element.
        dep_warning
            Dependency warning message if applicable.

        Returns
        -------
        str
            Formatted row content.
        """
        completed = checkbox_class == "completed-indicator"
        checkbox = (
            f"[bold green]{checkbox_content}[/]"
            if completed
            else f"[bold]{checkbox_content}[/]"
        )
        parts = [
            checkbox,
            f"[bold cyan]{num}[/]",
            f"[bold]{script}[/]",
            f"[dim]{desc}[/]",
        ]

        if ckpt:
            parts.append(f"[dim italic]({ckpt})[/]")

        if dep_warning:
            parts.append(f"[bold yellow]{dep_warning}[/]")

        return " ".join(parts)

    def _check_dependency(self, step_idx: int) -> str | None:
        """Check if a step has missing dependencies.

        Parameters
        ----------
        step_idx
            Index of the step to check.

        Returns
        -------
        str | None
            Warning message if dependencies are missing, None otherwise.
        """
        if step_idx not in self.selected_steps:
            return None

        # Use external dependency checker if provided
        if self._get_step_dependency is not None:
            required_ckpt = self._get_step_dependency(step_idx, self.modality)
            if required_ckpt:
                # Find which step produces this checkpoint
                for idx, ckpt in enumerate(self.checkpoints):
                    if ckpt == required_ckpt and idx not in self.completed_steps:
                        return "⚠ Dependency missing"
        else:
            # Fallback: check if previous step is completed (simple linear dependency)
            if step_idx > 0 and (step_idx - 1) not in self.completed_steps:
                return "⚠ Dependency missing"

        return None

    def on_click(self, event) -> None:
        """Handle click on step rows to toggle selection."""
        # Check if a step row was clicked
        if hasattr(event.target, "_step_idx"):
            step_idx = event.target._step_idx
            self._toggle_selection(step_idx)
            event.stop()

    def _toggle_selection(self, step_idx: int) -> None:
        """Toggle selection for a step.

        Parameters
        ----------
        step_idx
            Index of the step to toggle.
        """
        if step_idx in self.selected_steps:
            self.selected_steps.discard(step_idx)
        else:
            self.selected_steps.add(step_idx)

        # Trigger re-render to update visual state
        self._render_steps()

    def update_checkpoint_status(self, statuses: list[bool]) -> None:
        """Mark steps as completed based on checkpoint status.

        Parameters
        ----------
        statuses
            List of booleans indicating whether each step's checkpoint exists.
            Index corresponds to step index.
        """
        self.completed_steps = {i for i, s in enumerate(statuses) if s}

    @property
    def selected_indices(self) -> list[int]:
        """Return sorted list of selected step indices.

        Returns
        -------
        list[int]
            Sorted indices of currently selected steps.
        """
        return sorted(self.selected_steps)

    def set_dependency_checker(self, func: Callable[[int, str], str]) -> None:
        """Set an external dependency checking function.

        Parameters
        ----------
        func
            Function that returns the checkpoint file a step reads from.
            Signature: (step_idx: int, modality: str) -> str
        """
        self._get_step_dependency = func
        self._render_steps()

    def select_all(self) -> None:
        """Select all steps."""
        self.selected_steps = set(range(len(self.steps)))
        self._render_steps()

    def deselect_all(self) -> None:
        """Deselect all steps."""
        self.selected_steps = set()
        self._render_steps()

    def select_range(self, start_idx: int, end_idx: int) -> None:
        """Select a contiguous range of steps.

        Parameters
        ----------
        start_idx
            Start of range (inclusive).
        end_idx
            End of range (inclusive).
        """
        for idx in range(start_idx, end_idx + 1):
            self.selected_steps.add(idx)
        self._render_steps()