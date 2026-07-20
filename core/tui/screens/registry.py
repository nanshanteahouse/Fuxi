"""Registry browser screen — searchable, sortable DataTable of papers and datasets.

Loads registry data from the real ``MasterRegistry`` API (``load_master_registry``)
via ``asyncio.to_thread`` so the event loop is never blocked.  Read-only browser
for Phase 1 — no edit / register / delete capabilities.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

from rich.text import Text
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.reactive import reactive
from textual.screen import Screen
from textual.widgets import DataTable, Footer, Header, Input, Static

from core.paper.registry import (
    DatasetEntry,
    DatasetStatus,
    InsightStatus,
    PaperEntry,
    load_master_registry,
)

logger = logging.getLogger(__name__)

# ── Status colour mapping ───────────────────────────────────────────────

_STATUS_COLORS: dict[str, str] = {
    # DatasetStatus
    DatasetStatus.DATA_DOWNLOADED: "green",
    DatasetStatus.CONFIG_EXISTS: "blue",
    DatasetStatus.PIPELINE_COMPLETE: "cyan",
    DatasetStatus.DATA_NOT_DOWNLOADED: "grey",
    DatasetStatus.ORPHAN: "yellow",
    DatasetStatus.UNKNOWN: "grey",
    # InsightStatus
    InsightStatus.GENERATED: "green",
    InsightStatus.PENDING: "yellow",
    InsightStatus.FAILED: "red",
    InsightStatus.NO_GEO: "grey",
    InsightStatus.PENDING_REVIEW: "yellow",
    InsightStatus.PREPRINT: "blue",
    InsightStatus.PDF_ONLY: "grey",
}


def _status_badge(status: str) -> Text:
    """Return a coloured Rich ``Text`` badge for the given status string."""
    label = status.lower().replace("_", " ").title()
    colour = _STATUS_COLORS.get(status, "grey")
    return Text(f" {label} ", style=f"bold {colour} on {colour}")


# ═══════════════════════════════════════════════════════════════════════════
# Screen
# ═══════════════════════════════════════════════════════════════════════════


class RegistryBrowserScreen(Screen):
    """Read-only registry browser with search, sort, and detail inspection.

    Layout
    ------
    ::

        ┌─ Header ────────────────────────────────────────────┐
        │  Registry Browser                                     │
        ├─ Search ────────────────────────────────────────────┤
        │  🔍  Search by GSE, PMID, title, species…            │
        ├── DataTable ───────────────┬── Detail Panel ────────┤
        │  Type │ ID    │ Title …   │  Paper / Dataset fields │
        │  ─────┼───────┼───────────┤  …                      │
        │  Paper│ PMID  │ …         │                          │
        │  DS   │ GSE   │ …         │                          │
        ├─ Footer ────────────────────────────────────────────┤
        │  Fuxi …                                              │
        └──────────────────────────────────────────────────────┘
    """
    id = "registry"

    # ── CSS ────────────────────────────────────────────────────────────
    DEFAULT_CSS = """
    RegistryBrowserScreen {
        height: 100%;
    }

    #registry-header {
        padding: 1 2;
        background: $bg-medium;
        border-bottom: solid $accent;
        text-style: bold;
        color: $accent;
        height: 3;
    }

    #search-container {
        padding: 0 2;
        background: $bg-medium;
        border-bottom: solid $border;
        height: auto;
    }

    #search-input {
        width: 100%;
    }

    #registry-body {
        height: 1fr;
    }

    #table-panel {
        height: 100%;
    }

    #registry-table {
        height: 1fr;
        margin: 1 0 0 2;
    }

    #status-bar {
        background: $bg-dark;
        color: $text-muted;
        padding: 0 2;
        height: 1;
        border-top: solid $border;
        text-style: italic;
    }

    #empty-message {
        color: $text-muted;
        text-style: italic;
        padding: 2 4;
        height: 100%;
        content-align: center middle;
    }

    #loading-message {
        color: $text-muted;
        text-style: italic;
        padding: 2 4;
        height: 100%;
        content-align: center middle;
    }

    #error-message {
        color: $error;
        text-style: bold;
        padding: 2 4;
        height: 100%;
        content-align: center middle;
    }

    /* ── detail panel ──────────────────────────────────────────────── */
    #detail-panel {
        width: 42;
        height: 100%;
        background: $bg-medium;
        border-left: solid $border;
        padding: 1 2;
        overflow-y: auto;
    }

    #detail-title {
        text-style: bold;
        color: $accent;
        padding-bottom: 1;
        border-bottom: solid $border-light;
        margin-bottom: 1;
    }

    #detail-content {
        color: $text-primary;
    }

    #detail-content .label {
        text-style: bold;
        color: $text-secondary;
    }

    #detail-content .value {
        color: $text-primary;
    }
    """

    # ── reactive ──────────────────────────────────────────────────────
    filter_text: reactive[str] = reactive("", init=False)

    # ── lifecycle -------------------------------------------------------

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._registry: Any = None
        self._all_rows: list[dict[str, Any]] = []
        self._filtered_rows: list[dict[str, Any]] = []
        self._error: str | None = None

    # ── compose ─────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static(" Registry Browser", id="registry-header")
        with Horizontal(id="search-container"):
            yield Input(
                placeholder="Search by GSE, PMID, title, species…",
                id="search-input",
            )
        with Horizontal(id="registry-body"):
            with Vertical(id="table-panel"):
                yield Static("Loading registry …", id="loading-message")
                yield Static("", id="error-message")
                yield DataTable(id="registry-table")
                yield Static("", id="status-bar")
            with Vertical(id="detail-panel"):
                yield Static("Select a row to inspect", id="detail-title")
                yield Static(
                    "↕  Navigate with arrow keys\n"
                    "🔍  Type to search\n"
                    "⏎  Press Enter to select",
                    id="detail-content",
                )
        yield Footer()

    # ── mount / data loading --------------------------------------------

    async def on_mount(self) -> None:
        """Load registry data when the screen mounts."""
        await self._load_registry()

    async def _load_registry(self) -> None:
        """Load the master registry via ``asyncio.to_thread``."""
        try:
            self._registry = await asyncio.to_thread(load_master_registry)
        except Exception as exc:
            logger.exception("Failed to load registry")
            self._error = f"Failed to load registry: {exc}"
            self._show_error()
        else:
            self._build_rows()

    # ── error display ---------------------------------------------------

    def _show_error(self) -> None:
        """Show the error message and hide other content widgets."""
        loading = self.query_one("#loading-message", Static)
        error = self.query_one("#error-message", Static)
        table = self.query_one("#registry-table", DataTable)
        empty = self.query_one("#status-bar", Static)

        loading.display = False
        error.display = True
        error.update(f"[bold red]Error:[/]\n{self._error}")
        table.display = False
        empty.display = True
        empty.update("⚠ Registry unavailable")

    # ── data building ---------------------------------------------------

    def _build_rows(self) -> None:
        """Build the internal row list from the loaded registry."""
        if self._registry is None:
            return

        self._all_rows = []

        # Papers ---------------------------------------------------------
        for paper in self._registry.papers:
            pmid_display = paper.pmid or paper.paper_id or "?"
            status = (
                paper.insights.status.value
                if paper.insights
                else InsightStatus.PDF_ONLY.value
            )
            self._all_rows.append({
                "type": "Paper",
                "id": pmid_display,
                "title": paper.title or "?",
                "status": status,
                "year": paper.year or "",
                "journal": paper.journal or "",
                "author": paper.first_author or "",
                "doi": paper.doi or "",
                "slug": paper.slug or "",
                "paper_id": paper.paper_id,
                "_entry": paper,
            })

        # Datasets -------------------------------------------------------
        for ds_id, ds in self._registry.datasets.items():
            title = (
                ds.data_root.replace("{FUXI_DATA_ROOT}/", "")
                if ds.data_root
                else ds_id
            )
            species = ds.species or ""
            tissue = ds.tissue or ""
            modality_str = ", ".join(ds.modalities.keys()) if ds.modalities else ""
            links = self._registry.get_paper_links(ds_id)
            linked_pmids = ", ".join(p for p, _ in links) if links else ""
            self._all_rows.append({
                "type": "Dataset",
                "id": ds_id,
                "title": title,
                "status": ds.status or DatasetStatus.UNKNOWN.value,
                "species": species,
                "tissue": tissue,
                "modalities": modality_str,
                "linked_papers": linked_pmids,
                "n_samples": ds.n_samples,
                "n_cells": ds.n_cells,
                "repository": (
                    ds.repository.value if ds.repository else ""
                ),
                "data_root": ds.data_root or "",
                "_entry": ds,
            })

        self._populate_table()

    # ── table population ------------------------------------------------

    def _populate_table(self) -> None:
        """Populate or refresh the DataTable from ``_all_rows`` with the
        current filter applied."""
        table = self.query_one("#registry-table", DataTable)
        loading = self.query_one("#loading-message", Static)
        error = self.query_one("#error-message", Static)
        status_bar = self.query_one("#status-bar", Static)
        detail_title = self.query_one("#detail-title", Static)
        detail_content = self.query_one("#detail-content", Static)

        loading.display = False
        error.display = False

        # Columns — set once
        if not table.columns:
            table.add_column("Type", key="type")
            table.add_column("ID", key="id")
            table.add_column("Title", key="title")
            table.add_column("Status", key="status")
            table.add_column("Year / Species", key="extra")

        table.clear()

        # Filter ---------------------------------------------------------
        q = self.filter_text.strip().lower()
        self._filtered_rows = [
            row
            for row in self._all_rows
            if not q or self._row_matches(row, q)
        ]

        # Empty state ----------------------------------------------------
        if not self._filtered_rows:
            table.display = False
            if q:
                status_bar.update(f"No entries matching “{q}”")
                status_bar.display = True
            else:
                status_bar.display = True
                status_bar.update("No entries in registry")
            detail_title.update("No entries")
            detail_content.update("The registry is empty — no papers or datasets to display.")
            return

        table.display = True

        # Sort: datasets first, then papers; by ID within each group
        self._filtered_rows.sort(
            key=lambda r: (0 if r["type"] == "Dataset" else 1, r["id"]),
        )

        for row in self._filtered_rows:
            badge = _status_badge(row["status"])
            extra = (
                row.get("year")
                or row.get("species")
                or ""
            )
            table.add_row(
                row["type"],
                row["id"],
                (row["title"] or "")[:120],
                badge,
                extra,
            )

        # Status bar
        total = len(self._all_rows)
        shown = len(self._filtered_rows)
        n_papers = sum(1 for r in self._all_rows if r["type"] == "Paper")
        n_datasets = sum(1 for r in self._all_rows if r["type"] == "Dataset")
        if q:
            status_bar.update(
                f"{shown} of {total} entries  ·  {n_papers} papers, {n_datasets} datasets"
                f"  ·  filter: “{q}”"
            )
        else:
            status_bar.update(
                f"{total} entries  ·  {n_papers} papers, {n_datasets} datasets"
            )
        status_bar.display = True

        # Reset detail panel when repopulating
        detail_title.update("Select a row to inspect")
        detail_content.update(
            "↕  Navigate with arrow keys\n🔍  Type to search\n⏎  Press Enter to select"
        )

    @staticmethod
    def _row_matches(row: dict[str, Any], q: str) -> bool:
        """Return ``True`` when *q* appears in any searchable field."""
        fields = [
            row.get("id", ""),
            row.get("title", ""),
            row.get("year", ""),
            row.get("species", ""),
            row.get("journal", ""),
            row.get("author", ""),
            row.get("tissue", ""),
            row.get("modalities", ""),
            row.get("doi", ""),
            row.get("slug", ""),
            row.get("paper_id", ""),
            row.get("linked_papers", ""),
        ]
        return any(q in f.lower() for f in fields if f)

    # ── search ──────────────────────────────────────────────────────────

    def watch_filter_text(self, old: str, new: str) -> None:
        """Re-filter the table whenever the search text changes."""
        if self._all_rows:
            self._populate_table()

    def on_input_changed(self, event: Input.Changed) -> None:
        """Handle live search as the user types."""
        if event.input.id == "search-input":
            self.filter_text = event.value

    # ── detail panel ────────────────────────────────────────────────────

    def on_data_table_row_highlighted(
        self,
        event: DataTable.RowHighlighted,
    ) -> None:
        """Show full entry details when a row is highlighted."""
        self._update_detail(event.cursor_row)

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        """Also update on explicit selection (Enter / click)."""
        self._update_detail(event.cursor_row)

    def _update_detail(self, cursor_row: int) -> None:
        """Refresh the detail panel for the row at *cursor_row*."""
        if not self._filtered_rows:
            return
        if cursor_row < 0 or cursor_row >= len(self._filtered_rows):
            return

        row = self._filtered_rows[cursor_row]
        entry = row.get("_entry")
        if entry is None:
            return

        detail_title = self.query_one("#detail-title", Static)
        detail_content = self.query_one("#detail-content", Static)

        if isinstance(entry, PaperEntry):
            detail_title.update("📄 Paper Details")
            detail_content.update(self._format_paper_detail(entry, row))
        elif isinstance(entry, DatasetEntry):
            detail_title.update("📦 Dataset Details")
            detail_content.update(self._format_dataset_detail(entry, row))

    # ── formatters ──────────────────────────────────────────────────────

    @staticmethod
    def _format_paper_detail(paper: PaperEntry, row: dict[str, Any]) -> str:
        """Build the detail panel content for a paper entry."""
        lines: list[str] = []

        def kv(key: str, val: str, fallback: str = "—") -> None:
            lines.append(f"[bold]{key}[/]")
            lines.append(f"  {val or fallback}")

        kv("PM ID", paper.pmid or row.get("id", "—"))
        if paper.paper_id:
            kv("Paper ID", paper.paper_id)
        if paper.slug:
            kv("Slug", paper.slug)
        lines.append("")
        lines.append(f"[bold]Title[/]")
        lines.append(f"  {paper.title or '—'}")
        lines.append("")
        kv("Journal", paper.journal, fallback="—")
        kv("Year", paper.year, fallback="—")
        if paper.first_author:
            kv("First Author", paper.first_author)
        if paper.doi:
            kv("DOI", paper.doi)

        if paper.insights:
            lines.append("")
            lines.append("[bold]Insights[/]")
            lines.append(
                f"  Status: [bold]{paper.insights.status.value}[/]"
            )
            if paper.insights.insights_path:
                lines.append(f"  Path:   {paper.insights.insights_path}")

        if paper.supplements:
            lines.append("")
            lines.append(f"[bold]Supplements[/]  ({len(paper.supplements)} source(s))")

        if paper.kb_sources:
            lines.append("")
            lines.append(
                f"[bold]KB Sources[/]  ({len(paper.kb_sources)} entry/entries)"
            )

        cross = paper.cross_references or {}
        also_cited = cross.get("also_cited_by", [])
        if also_cited:
            lines.append("")
            lines.append(f"[bold]Cited by[/]  ({len(also_cited)} paper(s))")
            for pmid in also_cited[:8]:
                lines.append(f"  • {pmid}")
            if len(also_cited) > 8:
                lines.append(f"  … and {len(also_cited) - 8} more")

        return "\n".join(lines)

    @staticmethod
    def _format_dataset_detail(ds: DatasetEntry, row: dict[str, Any]) -> str:
        """Build the detail panel content for a dataset entry."""
        lines: list[str] = []

        lines.append(f"[bold]GSE ID[/]")
        lines.append(f"  {row['id']}")
        lines.append(f"[bold]Type[/]")
        lines.append(f"  {ds.type or '—'}")
        lines.append(f"[bold]Status[/]")
        status_val = ds.status or "unknown"
        colour = _STATUS_COLORS.get(status_val, "grey")
        lines.append(f"  [{colour}]{status_val}[/]")
        if ds.data_root:
            lines.append(f"[bold]Data Root[/]")
            lines.append(f"  {ds.data_root}")

        lines.append("")
        lines.append("[bold]Biology[/]")
        lines.append(f"  Species: {ds.species or '—'}")
        lines.append(f"  Tissue:  {ds.tissue or '—'}")
        if ds.n_samples:
            lines.append(f"  Samples: {ds.n_samples}")
        if ds.n_cells is not None:
            lines.append(f"  Cells:   {ds.n_cells:,}" if ds.n_cells > 0 else "  Cells:   0")
        if ds.data_format:
            lines.append(f"  Format:  {ds.data_format}")

        if ds.modalities:
            lines.append("")
            lines.append("[bold]Modalities[/]")
            for mod_key, mod_info in ds.modalities.items():
                mod_colour = _STATUS_COLORS.get(mod_info.status, "grey")
                lines.append(f"  {mod_key}: [{mod_colour}]{mod_info.status}[/]")

        if ds.subseries:
            lines.append("")
            lines.append("[bold]Sub-series[/]")
            for sub in ds.subseries:
                sid = sub.get("id", "?")
                note = sub.get("note", "")
                line = f"  {sid}"
                if note:
                    line += f"  —  {note}"
                lines.append(line)

        if ds.relationships:
            lines.append("")
            lines.append("[bold]Relationships[/]")
            for rel in ds.relationships:
                lines.append(f"  {rel.type.value}  →  {rel.dataset_id}")

        linked = row.get("linked_papers", "")
        if linked:
            lines.append("")
            lines.append(f"[bold]Linked Papers[/]  ({linked})")

        # Associations from dataset entry
        if hasattr(ds, "paper_pmids") and ds.paper_pmids:
            lines.append("")
            lines.append(f"[bold]paper_pmids[/]  ({', '.join(ds.paper_pmids[:5])})")
            if len(ds.paper_pmids) > 5:
                lines.append(f"  … and {len(ds.paper_pmids) - 5} more")

        if ds.notes:
            lines.append("")
            lines.append(f"[bold]Notes[/]")
            lines.append(f"  {ds.notes}")

        return "\n".join(lines)
