"""Data management screen — register papers, download GEO data, and generate configs.

Provides three tabbed panels:
- Register: Add papers by PMID to the registry
- Download: List and download GEO supplementary files
- Preprocess: Detect formats and generate pipeline configs

All long-running operations use async backends from core/tui/backends/.
"""

from __future__ import annotations

import asyncio
import logging
import os
from typing import Any

from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.reactive import reactive
from textual.screen import Screen
from textual.widgets import (
    Button,
    DataTable,
    Footer,
    Header,
    Input,
    ProgressBar,
    RichLog,
    Static,
    TabbedContent,
    TabPane,
)

from core.tui.backends.download import download_gse_async, fetch_meta_async, list_suppl_async
from core.tui.backends.preprocess import detect_formats_async, preprocess_async
from core.tui.backends.registry import register_paper_async

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# Screen
# ═══════════════════════════════════════════════════════════════════════════


class DataManagementScreen(Screen):
    """Three-panel screen for paper registration, GEO download, and preprocessing.

    Layout
    ------
    ::

        ┌─ Header ────────────────────────────────────────────┐
        │  Data Management                                     │
        ├─ TabbedContent ─────────────────────────────────────┤
        │  [Register] [Download] [Preprocess]                │
        │  ┌─ Register Panel ─────────────────────────────┐   │
        │  │  PMID Input + Register Button + RichLog      │   │
        │  ├─ Download Panel ──────────────────────────────┤   │
        │  │  GSE Input + List Files + DataTable +        │   │
        │  │  Download Button + ProgressBar + RichLog     │   │
        │  ├─ Preprocess Panel ────────────────────────────┤   │
        │  │  GSE/Dir Input + Detect Formats +            │   │
        │  │  Generate Config Button + RichLog            │   │
        │  └───────────────────────────────────────────────┘   │
        ├─ Footer ────────────────────────────────────────────┤
        └──────────────────────────────────────────────────────┘
    """

    id = "data-mgmt"

    # ── CSS ────────────────────────────────────────────────────────────
    DEFAULT_CSS = """
    DataManagementScreen {
        height: 100%;
    }

    #screen-header {
        padding: 1 2;
        background: $bg-medium;
        border-bottom: solid $accent;
        text-style: bold;
        color: $accent;
        height: 3;
    }

    /* ── TabbedContent ─────────────────────────────────────────────── */
    TabbedContent {
        height: 1fr;
    }

    TabPane {
        padding: 1 2;
        height: 100%;
    }

    /* ── Input rows ─────────────────────────────────────────────────── */
    .input-row {
        dock: top;
        height: auto;
        padding: 0 0 1 0;
    }

    .input-row Input {
        width: 1fr;
    }

    .input-row Button {
        margin-left: 1;
    }

    /* ── Buttons ───────────────────────────────────────────────────── */
    Button {
        width: auto;
    }

    /* ── DataTable (Download panel) ───────────────────────────────── */
    #files-table {
        height: 20;
        margin: 1 0;
    }

    /* ── ProgressBar ───────────────────────────────────────────────── */
    #progress-bar {
        margin: 1 0;
        display: none;
    }

    /* ── RichLog (all panels) ──────────────────────────────────────── */
    RichLog {
        height: 1fr;
        margin-top: 1;
        border: solid $border;
        background: $bg-dark;
    }

    /* ── Format display (Preprocess panel) ─────────────────────────── */
    #format-display {
        height: 15;
        margin: 1 0;
        border: solid $border;
        padding: 1;
        background: $bg-medium;
        overflow-y: auto;
    }

    /* ── SuperSeries warning ───────────────────────────────────────── */
    #superseries-warning {
        margin: 1 0;
        padding: 1;
        background: $warning-bg;
        border: solid $warning;
        color: $warning-fg;
        display: none;
    }
    """

    # ── reactive ──────────────────────────────────────────────────────
    data_root: reactive[str] = reactive("")

    # ── lifecycle -------------------------------------------------------

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._suppl_files: list[dict] = []
        self._selected_files: set[str] = set()
        self._download_gse_id: str = ""
        self._preprocess_gse_id: str = ""
        self._format_result: dict[str, Any] | None = None

    # ── compose ─────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        # Get data root from environment
        self.data_root = os.environ.get("FUXI_DATA_ROOT", "")

        yield Header()
        yield Static(" Data Management", id="screen-header")

        with TabbedContent(id="main-tabs"):
            # ── Register Panel ─────────────────────────────────────────
            with TabPane("Register", id="register-pane"):
                with Horizontal(classes="input-row"):
                    yield Input(
                        placeholder="PMID (e.g., 31493975)",
                        id="register-pmid-input",
                    )
                    yield Button("Register", id="register-button", variant="primary")
                yield RichLog(id="register-log", markup=True)

            # ── Download Panel ────────────────────────────────────────
            with TabPane("Download", id="download-pane"):
                with Horizontal(classes="input-row"):
                    yield Input(
                        placeholder="GSE ID (e.g., GSE123456)",
                        id="download-gse-input",
                    )
                    yield Button("List Files", id="list-files-button", variant="primary")
                yield Static("", id="superseries-warning")
                yield DataTable(id="files-table")
                with Horizontal(classes="input-row"):
                    yield Button(
                        "Download Selected",
                        id="download-button",
                        variant="success",
                        disabled=True,
                    )
                yield ProgressBar(id="progress-bar", show_eta=True)
                yield RichLog(id="download-log", markup=True)

            # ── Preprocess Panel ───────────────────────────────────────
            with TabPane("Preprocess", id="preprocess-pane"):
                with Horizontal(classes="input-row"):
                    yield Input(
                        placeholder="GSE ID or data directory path",
                        id="preprocess-input",
                    )
                    yield Button(
                        "Detect Formats",
                        id="detect-formats-button",
                        variant="primary",
                    )
                yield Static("", id="format-display")
                with Horizontal(classes="input-row"):
                    yield Button(
                        "Generate Config",
                        id="generate-config-button",
                        variant="success",
                        disabled=True,
                    )
                yield RichLog(id="preprocess-log", markup=True)

        yield Footer()

    # ── mount -----------------------------------------------------------

    def on_mount(self) -> None:
        """Initialize the download table when the screen mounts."""
        self._setup_files_table()

    def _setup_files_table(self) -> None:
        """Set up the DataTable columns for the download panel."""
        table = self.query_one("#files-table", DataTable)
        table.add_column("Select", key="select", width=8)
        table.add_column("File Name", key="filename", width=50)
        table.add_column("Size", key="size", width=12)
        table.add_column("Type", key="type", width=15)

    # ── Register Panel ─────────────────────────────────────────────────

    async def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses across all panels."""
        button_id = event.button.id

        if button_id == "register-button":
            await self._handle_register()
        elif button_id == "list-files-button":
            await self._handle_list_files()
        elif button_id == "download-button":
            await self._handle_download()
        elif button_id == "detect-formats-button":
            await self._handle_detect_formats()
        elif button_id == "generate-config-button":
            await self._handle_generate_config()

    async def _handle_register(self) -> None:
        """Register a paper by PMID."""
        pmid_input = self.query_one("#register-pmid-input", Input)
        log = self.query_one("#register-log", RichLog)

        pmid = pmid_input.value.strip()
        if not pmid:
            log.write("[bold red]Error:[/] Please enter a PMID\n")
            return

        log.write(f"[cyan]Registering paper PMID: {pmid}[/]\n")
        log.write("[dim]Fetching metadata from NCBI...[/]\n")

        try:
            result = await register_paper_async(pmid)
            if result:
                log.write(f"[bold green]✓[/] Successfully registered paper {pmid}\n")
            else:
                log.write(f"[bold yellow]⚠[/] Registration returned no result for PMID {pmid}\n")
        except Exception as exc:
            log.write(f"[bold red]Error:[/] {exc}\n")
            logger.exception("Failed to register paper")

    # ── Download Panel ─────────────────────────────────────────────────

    async def _handle_list_files(self) -> None:
        """List supplementary files for a GSE dataset."""
        gse_input = self.query_one("#download-gse-input", Input)
        log = self.query_one("#download-log", RichLog)
        table = self.query_one("#files-table", DataTable)
        warning = self.query_one("#superseries-warning", Static)
        download_button = self.query_one("#download-button", Button)

        gse_id = gse_input.value.strip()
        if not gse_id:
            log.write("[bold red]Error:[/] Please enter a GSE ID\n")
            return

        # Normalize GSE ID
        if not gse_id.upper().startswith("GSE"):
            gse_id = f"GSE{gse_id}"

        self._download_gse_id = gse_id
        log.write(f"[cyan]Fetching metadata for {gse_id}[/]\n")

        # Check for SuperSeries
        try:
            metadata = await fetch_meta_async(gse_id)
            is_superseries = metadata.get("is_superseries", False)

            if is_superseries:
                log.write(f"[bold yellow]⚠ SuperSeries Detected:[/] {gse_id} is a SuperSeries\n")
                subseries = metadata.get("subseries", [])
                if subseries:
                    warning.update(
                        f"[bold]SuperSeries Warning:[/] This is a SuperSeries containing "
                        f"{len(subseries)} sub-series:\n"
                        + "\n".join(f"  • {s}" for s in subseries)
                        + "\n\nConsider downloading individual sub-series instead."
                    )
                    warning.display = True
                else:
                    warning.display = True
                    warning.update(
                        "[bold]SuperSeries Warning:[/] This is a SuperSeries. "
                        "Consider downloading individual sub-series instead."
                    )
            else:
                warning.display = False
        except Exception as exc:
            log.write(f"[bold yellow]Warning:[/] Could not fetch metadata: {exc}\n")
            warning.display = False

        # List supplementary files
        try:
            log.write("[dim]Listing supplementary files...[/]\n")
            files = await list_suppl_async(gse_id)
            self._suppl_files = files

            table.clear()

            if not files:
                log.write(f"[bold yellow]⚠[/] No supplementary files found for {gse_id}\n")
                download_button.disabled = True
                return

            for idx, file_info in enumerate(files):
                filename = file_info.get("name", "")
                size = file_info.get("size_human", "")
                # Infer file type from extension
                ext = os.path.splitext(filename)[1].lower()
                if ext in (".gz",):
                    ext2 = os.path.splitext(os.path.splitext(filename)[0])[1].lower()
                    file_type = ext2 + ext if ext2 else ext
                elif ext in (".csv", ".tsv"):
                    file_type = ext
                elif ext in (".h5", ".h5ad"):
                    file_type = ext
                elif ext in (".mtx",):
                    file_type = ext
                else:
                    file_type = ext.lstrip(".") or "unknown"
                table.add_row(
                    "[ ]",  # Select checkbox
                    filename[:50],  # Truncate long filenames
                    size,
                    file_type,
                    key=f"file_{idx}",
                )

            log.write(f"[bold green]✓[/] Found {len(files)} supplementary file(s)\n")
            download_button.disabled = False

        except Exception as exc:
            log.write(f"[bold red]Error:[/] Failed to list files: {exc}\n")
            logger.exception("Failed to list supplementary files")
            table.clear()
            download_button.disabled = True

    async def _handle_download(self) -> None:
        """Download the selected supplementary files."""
        log = self.query_one("#download-log", RichLog)
        progress = self.query_one("#progress-bar", ProgressBar)

        if not self._download_gse_id:
            log.write("[bold red]Error:[/] No GSE ID specified\n")
            return

        if not self.data_root:
            log.write("[bold red]Error:[/] FUXI_DATA_ROOT environment variable not set\n")
            return

        progress.total = 100
        progress.display = True
        progress.advance(0)

        log.write(f"[cyan]Starting download for {self._download_gse_id}[/]\n")
        log.write(f"[dim]Destination: {self.data_root}[/]\n")

        try:
            progress.advance(10)
            log.write("[dim]Initializing download...[/]\n")

            line_count = 0
            async for line in download_gse_async(self._download_gse_id, self.data_root):
                log.write(f"{line}\n")
                line_count += 1

                # Update progress (simple heuristic)
                if line_count % 10 == 0:
                    progress.advance(1)
                    if progress.progress >= 90:
                        progress.progress = 90

            progress.advance(100 - progress.progress)
            log.write(f"[bold green]✓[/] Download completed for {self._download_gse_id}\n")

        except asyncio.CancelledError:
            log.write("[bold yellow]Download cancelled[/]\n")
            raise
        except Exception as exc:
            log.write(f"[bold red]Error:[/] Download failed: {exc}\n")
            logger.exception("Download failed")
        finally:
            progress.display = False

    # ── Preprocess Panel ────────────────────────────────────────────────

    async def _handle_detect_formats(self) -> None:
        """Detect file formats in the specified directory."""
        input_field = self.query_one("#preprocess-input", Input)
        log = self.query_one("#preprocess-log", RichLog)
        format_display = self.query_one("#format-display", Static)
        generate_button = self.query_one("#generate-config-button", Button)

        input_val = input_field.value.strip()
        if not input_val:
            log.write("[bold red]Error:[/] Please enter a GSE ID or directory path\n")
            return

        # Determine if it's a GSE ID or a path
        if input_val.upper().startswith("GSE"):
            # It's a GSE ID - construct the data path
            gse_id = input_val.upper()
            data_dir = os.path.join(self.data_root, gse_id)
            self._preprocess_gse_id = gse_id
        else:
            # It's a directory path - extract GSE ID from path if possible
            data_dir = input_val
            # Try to extract GSE ID from the directory name
            basename = os.path.basename(os.path.normpath(data_dir))
            if basename.upper().startswith("GSE"):
                self._preprocess_gse_id = basename.upper()
            else:
                self._preprocess_gse_id = basename

        if not os.path.exists(data_dir):
            log.write(f"[bold red]Error:[/] Directory does not exist: {data_dir}\n")
            return

        log.write(f"[cyan]Detecting formats in: {data_dir}[/]\n")
        log.write("[dim]Scanning directory structure...[/]\n")

        try:
            result = await detect_formats_async(data_dir)
            self._format_result = result

            # Format the result for display
            display_lines = []

            modalities = result.get("modalities", {})
            if modalities:
                display_lines.append("[bold]Detected Modalities:[/]")
                for mod, info in modalities.items():
                    display_lines.append(f"  • {mod}: {info}")
                display_lines.append("")

            samples = result.get("samples", {})
            if samples:
                display_lines.append("[bold]Samples:[/]")
                for sample_name, sample_info in samples.items():
                    display_lines.append(f"  • {sample_name}:")
                    if isinstance(sample_info, dict):
                        for key, val in sample_info.items():
                            display_lines.append(f"      {key}: {val}")
                    else:
                        display_lines.append(f"      {sample_info}")
                display_lines.append("")

            unmatched = result.get("unmatched_files", [])
            if unmatched:
                display_lines.append(f"[bold]Unmatched Files:[/] ({len(unmatched)})")
                for f in unmatched[:20]:  # Show first 20
                    display_lines.append(f"  • {f}")
                if len(unmatched) > 20:
                    display_lines.append(f"  … and {len(unmatched) - 20} more")
                display_lines.append("")

            if not display_lines:
                display_lines.append("[dim]No formats detected[/]")

            format_display.update("\n".join(display_lines))
            log.write(
                f"[bold green]✓[/] Format detection complete for {self._preprocess_gse_id}\n"
            )
            generate_button.disabled = False

        except Exception as exc:
            log.write(f"[bold red]Error:[/] Format detection failed: {exc}\n")
            logger.exception("Format detection failed")
            format_display.update("[bold red]Detection failed[/]")
            generate_button.disabled = True

    async def _handle_generate_config(self) -> None:
        """Generate pipeline config for the detected dataset."""
        log = self.query_one("#preprocess-log", RichLog)

        if not self._preprocess_gse_id:
            log.write("[bold red]Error:[/] No GSE ID specified\n")
            return

        if not self.data_root:
            log.write("[bold red]Error:[/] FUXI_DATA_ROOT environment variable not set\n")
            return

        log.write(f"[cyan]Generating config for {self._preprocess_gse_id}[/]\n")
        log.write(f"[dim]Data root: {self.data_root}[/]\n")

        try:
            log.write("[dim]Running preprocessor...[/]\n")

            line_count = 0
            async for line in preprocess_async(
                self._preprocess_gse_id,
                self.data_root,
                modality=None,  # Auto-detect
                skip_extract=False,
                query_ncbi=False,
            ):
                log.write(f"{line}\n")
                line_count += 1

            log.write(
                f"[bold green]✓[/] Config generation complete for {self._preprocess_gse_id}\n"
            )

        except asyncio.CancelledError:
            log.write("[bold yellow]Config generation cancelled[/]\n")
            raise
        except Exception as exc:
            log.write(f"[bold red]Error:[/] Config generation failed: {exc}\n")
            logger.exception("Config generation failed")
