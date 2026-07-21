"""Config selector — config-file browser widget for the Fuxi pipeline TUI.

Lists available YAML config files found in the ``projects/`` directory
tree, grouped by modality with collapsible sections, and reports checkpoint
pipeline status visually.

Usage
-----
    selector = ConfigSelector()
    selector.reload_configs()                       # scan all modalities
    selector.select("rna", "GSE123456")      # programmatic selection
    path = selector.selected_path            # → str | None

    # Or filter to a single modality
    selector2 = ConfigSelector(modality="atac")
    selector2.reload_configs()

Events
------
    :class:`ConfigSelected`
        Posted when the user clicks a config leaf node.
    :class:`ConfigConfirmed`
        Posted when the user double-clicks a config leaf node.
"""

from __future__ import annotations

import glob as _glob
import os as _os
from typing import Any

from textual.app import ComposeResult
from textual.message import Message
from textual.reactive import reactive
from textual.widget import Widget
from textual.widgets import Static, Tree
from textual.widgets._tree import TreeNode

from core.utils._path import repo_root

# ═══════════════════════════════════════════════════════════════════════
#  Custom messages
# ═══════════════════════════════════════════════════════════════════════


class ConfigSelected(Message):
    """Posted when a config leaf is clicked (single-click).

    Attributes
    ----------
    path : str
        Absolute path to the selected YAML config file.
    """

    def __init__(self, path: str) -> None:
        self.path = path
        super().__init__()


class ConfigConfirmed(Message):
    """Posted when a config leaf is double-clicked.

    Attributes
    ----------
    path : str
        Absolute path to the confirmed YAML config file.
    """

    def __init__(self, path: str) -> None:
        self.path = path
        super().__init__()


# ═══════════════════════════════════════════════════════════════════════
#  Widget
# ═══════════════════════════════════════════════════════════════════════

_CHECKPOINT_MARKER = "●"
_NO_CHECKPOINT_MARKER = "○"


class ConfigSelector(Widget):
    """Config-file browser widget for selecting pipeline configs.

    Scans the ``projects/`` directory tree for ``config_*.yaml`` files
    (expected layout: ``projects/{modality}/{GSE_ID}/config_{GSE_ID}.yaml``),
    groups results by modality, and presents them as a collapsible tree.

    *Single-click* a leaf to select it (updates :attr:`selected_path` and
    posts :class:`ConfigSelected`).  *Double-click* a leaf to confirm,
    which additionally posts :class:`ConfigConfirmed`.

    Attributes
    ----------
    modality : str | None
        When set, only this modality subdirectory is scanned.
    selected_path : str | None
        Absolute path of the currently selected config file (read-only).
    available_configs : dict[str, dict[str, str]]
        The scanned results: ``{modality: {gse_id: config_path}}``.
    """

    DEFAULT_CSS = """
    ConfigSelector {
        height: 100%;
        layout: vertical;
    }

    ConfigSelector > #cs-header {
        height: 3;
        content-align: left middle;
        padding: 0 1;
        background: $bg-medium;
        border-bottom: solid $border;
        text-style: bold;
        color: $text;
    }

    ConfigSelector > Tree {
        height: 1fr;
        background: $bg-dark;
        border: none;
        padding: 0 1;
        margin: 0;
    }

    ConfigSelector > Tree:focus {
        border: none;
    }

    ConfigSelector > #cs-empty {
        height: 1fr;
        content-align: center middle;
        color: $text-muted;
        text-style: italic;
    }

    ConfigSelector > #cs-status {
        height: 1;
        padding: 0 1;
        color: $text-secondary;
        background: $bg-medium;
        border-top: solid $border;
    }
    """

    # ── Reactive ──────────────────────────────────────────────────────

    modality: reactive[str | None] = reactive(None, init=False)

    # ── Private state ──────────────────────────────────────────────────

    _projects_dir: str = ""
    _configs: dict[str, dict[str, str]] = {}
    _leaf_nodes: list[TreeNode[dict[str, Any]]] = []
    _selected_path: str | None = None

    # ── Init ───────────────────────────────────────────────────────────

    def __init__(self, modality: str | None = None, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.modality = modality

    def on_mount(self) -> None:
        """Resolve ``projects/`` and perform the initial scan."""
        self._projects_dir = _os.path.join(repo_root(), "projects")
        self._scan_and_render()

    def watch_modality(self, old: str | None, new: str | None) -> None:
        """Re-scan when the modality filter changes."""
        self._scan_and_render()

    # ── Public API ─────────────────────────────────────────────────────

    def reload_configs(self, modality: str | None = None) -> None:
        """Re-scan the ``projects/`` directory for config files.

        Parameters
        ----------
        modality
            If given, restrict the scan to this modality's subdirectory
            (e.g. ``"rna"``).  When *None*, uses the current filter.
        """
        if modality is not None:
            self.modality = modality
        else:
            self._scan_and_render()

    def select(self, modality: str, gse_id: str) -> bool:
        """Programmatically select a config by modality and GSE ID.

        Parameters
        ----------
        modality
            Modality directory name (e.g. ``"rna"``, ``"atac"``).
        gse_id
            GSE accession ID (e.g. ``"GSE123456"``).

        Returns
        -------
        bool
            ``True`` if the config was found and selected, ``False`` if
            no matching leaf exists.
        """
        if modality not in self._configs or gse_id not in self._configs[modality]:
            return False

        self._selected_path = self._configs[modality][gse_id]
        self._highlight_selected()
        self._update_status_bar()
        return True

    @property
    def selected_path(self) -> str | None:
        """Return the absolute path of the selected config file, or None."""
        return self._selected_path

    @property
    def available_configs(self) -> dict[str, dict[str, str]]:
        """Return the scanned configs: ``{modality: {gse_id: config_path}}``."""
        return {mod: dict(cfgs) for mod, cfgs in self._configs.items()}

    # ── Scan helpers ───────────────────────────────────────────────────

    def _scan_and_render(self) -> None:
        """Scan ``projects/`` and rebuild the tree display."""
        self._configs.clear()
        self._leaf_nodes.clear()
        self._selected_path = None
        self._scan_directory()
        self._rebuild_tree()

    def _scan_directory(self) -> None:
        """Walk ``projects/`` and collect every matching ``config_*.yaml``.

        Expected layout::

            projects/{modality}/{GSE_ID}/config_{GSE_ID}.yaml

        The glob pattern ``<modality>/*/config_*.yaml`` naturally skips
        directories without config files.
        """
        if not _os.path.isdir(self._projects_dir):
            return

        # Determine which modality sub-directories to walk.
        candidates: list[str] = []
        if self.modality:
            candidates = [self.modality]
        else:
            try:
                candidates = sorted(
                    e.name
                    for e in _os.scandir(self._projects_dir)
                    if e.is_dir() and not e.name.startswith(".")
                )
            except OSError:
                return

        for mod_name in candidates:
            mod_path = _os.path.join(self._projects_dir, mod_name)
            if not _os.path.isdir(mod_path):
                continue

            for config_path in sorted(_glob.glob(_os.path.join(mod_path, "*", "config_*.yaml"))):
                gse_id = self._gse_id_from_path(config_path)
                if gse_id:
                    self._configs.setdefault(mod_name, {})[gse_id] = config_path

    @staticmethod
    def _gse_id_from_path(config_path: str) -> str | None:
        """Derive the GSE ID from the parent directory name.

        The expected path is ``.../{modality}/{GSE_ID}/config_{GSE_ID}.yaml``,
        so the parent directory of the config file *is* the GSE ID.
        """
        parent = _os.path.basename(_os.path.dirname(config_path))
        if parent.startswith("."):
            return None
        return parent

    @staticmethod
    def _has_checkpoints(config_path: str) -> bool:
        """Return ``True`` if the pipeline generated any h5ad checkpoint.

        Checks for any ``.h5ad`` files inside
        ``<project-dir>/results/h5ad/``.
        """
        h5ad_dir = _os.path.join(_os.path.dirname(config_path), "results", "h5ad")
        if not _os.path.isdir(h5ad_dir):
            return False
        try:
            return any(f.endswith(".h5ad") for f in _os.listdir(h5ad_dir))
        except OSError:
            return False

    # ── Tree rendering ─────────────────────────────────────────────────

    def _rebuild_tree(self) -> None:
        """Build the :class:`~textual.widgets.Tree` from scanned configs."""
        tree = self.query_one(Tree)
        tree.clear()
        is_empty = not self._configs
        self.query_one("#cs-empty", Static).display = is_empty

        if is_empty:
            tree.display = False
            self._update_status_bar()
            return
        tree.display = True

        for mod_name in sorted(self._configs):
            datasets = self._configs[mod_name]
            mod_label = f"[bold]{mod_name.upper()}[/bold]  ({len(datasets)})"
            mod_node: TreeNode[dict[str, Any]] = tree.root.add(mod_label, expand=True)

            for gse_id in sorted(datasets):
                config_path = datasets[gse_id]
                has_ckpt = self._has_checkpoints(config_path)

                # Choose colour and marker based on checkpoint status.
                if has_ckpt:
                    marker = f"[green]{_CHECKPOINT_MARKER}[/green]"
                else:
                    marker = _NO_CHECKPOINT_MARKER

                label = f"{marker} {gse_id}   [dim]{config_path}[/dim]"
                leaf = mod_node.add_leaf(label)
                leaf.data = {
                    "modality": mod_name,
                    "gse_id": gse_id,
                    "config_path": config_path,
                    "has_checkpoints": has_ckpt,
                }
                self._leaf_nodes.append(leaf)

        tree.root.expand()

    def _highlight_selected(self) -> None:
        """Highlight the tree node matching ``_selected_path``."""
        if not self._selected_path:
            return
        tree = self.query_one(Tree)
        for node in self._leaf_nodes:
            if node.data and node.data.get("config_path") == self._selected_path:
                tree.select_node(node)
                # Scroll to the selected node so it's visible.
                tree.scroll_to_node(node, animate=False)
                return

    def _update_status_bar(self) -> None:
        """Update the status bar with current selection or summary info."""
        sb = self.query_one("#cs-status", Static)
        if self._selected_path:
            sb.update(f"Selected: [bold]{_os.path.basename(self._selected_path)}[/bold]")
        elif self._configs:
            total = sum(len(ds) for ds in self._configs.values())
            mods = ", ".join(sorted(self._configs.keys()))
            sb.update(f"{total} configs across [italic]{mods}[/italic]")
        else:
            sb.update("No config files found")

    # ── Event handlers ─────────────────────────────────────────────────

    def on_tree_node_selected(self, event: Tree.NodeSelected[dict[str, Any]]) -> None:
        """Handle single-click node selection.

        For leaf nodes (config files) the path is stored and a
        :class:`ConfigSelected` message is posted.  Branch (modality)
        nodes are ignored.
        """
        if event.node.data is None:
            return

        self._selected_path = event.node.data.get("config_path", "")
        self._update_status_bar()
        self.post_message(ConfigSelected(self._selected_path or ""))

    def on_click(self, event) -> None:
        """Detect double-click on leaf nodes; post :class:`ConfigConfirmed`."""

        # Let single-click handling complete first (selected_path is set by
        # on_tree_node_selected).  Then check for double-click.
        def _post_if_double() -> None:
            num_clicks = getattr(event, "num_clicks", 0)
            if num_clicks >= 2 and self._selected_path:
                self.post_message(ConfigConfirmed(self._selected_path))

        self.call_later(_post_if_double)

    # ── Compose ────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        yield Static("Config Selector", id="cs-header")
        yield Tree[dict[str, Any]]("Projects", id="cs-tree")
        yield Static("No config files found", id="cs-empty")
        yield Static("", id="cs-status")
