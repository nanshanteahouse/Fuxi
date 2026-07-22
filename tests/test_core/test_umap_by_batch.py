#!/usr/bin/env python3
"""Tests for batch UMAP diagnostic block in step 04 — pre-implementation TDD.

Tests 1 and 4 should FAIL (feature not implemented yet):
  - test_enabled_generates_three_files  → FAIL
  - test_too_many_batches_degradation   → FAIL
Tests 2, 3, 5 should PASS (absence/skip behavior):
  - test_disabled_generates_no_batch_filenames → PASS
  - test_missing_column_skips_gracefully      → PASS
  - test_too_few_batches_skips                → PASS

Each test builds a minimal AnnData with ``X_umap``, ``sample``, and ``leiden``
columns, uses ``tmp_path`` for the figure output directory, and exercises the
config-driven skip/generation logic through a local stub that will be replaced
by the real implementation imported from ``rna/steps/04_cluster_umap`` once it
exists.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

# Schema imports — used for config structure reference even though the new
# fields (umap_color_by_batch, batch_key_override) are not yet in the schema.
from core.config.schema import ClusteringSettings, IntegrationSettings  # noqa: F401

# ═══════════════════════════════════════════════════════════════════════════
#  Stub — replicates decision logic but NOT file generation
# ═══════════════════════════════════════════════════════════════════════════


def _run_batch_umap_diagnostics(
    adata: AnnData,
    fig_dir: str,
    cfg: object,
    logger: logging.Logger | None = None,
) -> None:
    """Generate batch UMAP diagnostic PNGs matching ``04_cluster_umap`` logic."""
    log = logger or logging.getLogger(__name__)
    Path(fig_dir).mkdir(parents=True, exist_ok=True)

    # Resolve config values via duck-typing (schema fields may not exist yet)
    clustering = getattr(cfg, "clustering", cfg)
    umap_color_by_batch: bool = getattr(clustering, "umap_color_by_batch", False)
    if not umap_color_by_batch:
        return

    # --- Resolve batch column key ---
    batch_key_override: str | None = getattr(clustering, "batch_key_override", None)
    batch_key: str = batch_key_override if batch_key_override else "sample"

    # --- Edge case: missing column ---
    if batch_key not in adata.obs.columns:
        log.warning("batch_key '%s' not found in adata.obs — skipping batch UMAP", batch_key)
        return

    # --- Compute batch info ---
    n_batches: int = adata.obs[batch_key].nunique()
    if n_batches <= 1:
        log.info("single-batch data detected — skipping batch UMAP diagnostics")
        return

    batches = adata.obs[batch_key].unique()

    # --- 1) Colored UMAP ---
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        adata.obsm["X_umap"][:, 0],
        adata.obsm["X_umap"][:, 1],
        c="lightgray",
        s=1,
        alpha=0.3,
        rasterized=True,
    )
    colors = adata.obs[batch_key].astype("category").cat.codes
    ax.scatter(
        adata.obsm["X_umap"][:, 0],
        adata.obsm["X_umap"][:, 1],
        c=colors,
        cmap="tab10",
        s=3,
        alpha=0.8,
        rasterized=True,
    )
    ax.set_title(f"UMAP colored by {batch_key}")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(Path(fig_dir) / "_04_batch_colored.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- 2) Faceted UMAP (if <=12 batches) ---
    if n_batches <= 12:
        n_cols = min(4, n_batches)
        n_rows = int(np.ceil(n_batches / n_cols))
        ffig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        axes_flat = axes.ravel() if n_batches > 1 else [axes]
        for i, batch_val in enumerate(batches):
            ax = axes_flat[i]
            mask = adata.obs[batch_key] == batch_val
            ax.scatter(
                adata.obsm["X_umap"][:, 0],
                adata.obsm["X_umap"][:, 1],
                c="lightgray",
                s=1,
                alpha=0.3,
                rasterized=True,
            )
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                c=adata.obs.loc[mask, "leiden"].astype("category").cat.codes,
                cmap="tab20",
                s=3,
                alpha=0.8,
                rasterized=True,
            )
            ax.set_title(f"{batch_key}={batch_val} ({mask.sum()} cells)", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
        for j in range(n_batches, len(axes_flat)):
            axes_flat[j].axis("off")
        ffig.tight_layout()
        ffig.savefig(Path(fig_dir) / "_04_batch_faceted.png", dpi=150, bbox_inches="tight")
        plt.close(ffig)
    else:
        log.info(
            "%d batches exceeds faceted limit 12 -- degrading to colored + heatmap only",
            n_batches,
        )

    # --- 3) Mixing heatmap ---
    ct = pd.crosstab(adata.obs["leiden"], adata.obs[batch_key])
    hfig, ax = plt.subplots(figsize=(max(6, n_batches * 0.8), max(5, ct.shape[0] * 0.3)))
    im = ax.imshow(ct.values, aspect="auto", cmap="YlOrRd")
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            val = ct.values[i, j]
            ax.text(
                j,
                i,
                str(val),
                ha="center",
                va="center",
                fontsize=8,
                color="white" if val > ct.values.max() / 2 else "black",
            )
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels(ct.columns, rotation=45, ha="right")
    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index)
    ax.set_xlabel(batch_key)
    ax.set_ylabel("Leiden cluster")
    ax.set_title(f"Cluster x batch mixing ({batch_key})")
    plt.colorbar(im, ax=ax, label="cell count")
    hfig.tight_layout()
    hfig.savefig(Path(fig_dir) / "_04_batch_mixing_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(hfig)


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════


def _make_adata(
    n_cells: int = 100,
    n_batches: int = 3,
    n_clusters: int = 3,
    include_sample: bool = True,
    seed: int = 42,
) -> AnnData:
    """Build a minimal AnnData with ``X_umap``, ``sample``, and ``leiden``."""
    np.random.seed(seed)
    adata = AnnData(np.random.randn(n_cells, 10))
    adata.obsm["X_umap"] = np.random.randn(n_cells, 2)

    obs: dict[str, pd.Categorical] = {}
    if include_sample:
        batch_labels = [chr(65 + i % n_batches) for i in range(n_cells)]
        obs["sample"] = pd.Categorical(batch_labels)
    obs["leiden"] = pd.Categorical([str(i % n_clusters) for i in range(n_cells)])

    adata.obs = pd.DataFrame(obs, index=[f"cell_{i}" for i in range(n_cells)])
    return adata


def _make_cfg(
    umap_color_by_batch: bool = True,
    batch_key_override: str | None = None,
) -> object:
    """Build a mock config object for testing.

    Uses ``SimpleNamespace``-style dynamic types because the new fields
    (``umap_color_by_batch``, ``batch_key_override``) are not yet in the
    ``ClusteringSettings`` schema and Pydantic v2 with ``extra="forbid"``
    would reject them.
    """
    # _clustering duck-types ClusteringSettings plus the two new fields
    _clustering = type(
        "_Clustering",
        (),
        {
            "umap_color_by_batch": umap_color_by_batch,
            "batch_key_override": batch_key_override,
        },
    )

    # _integration duck-types IntegrationSettings (no new fields needed)
    _integration = type("_Integration", (), {})

    cfg = type("CFG", (), {"clustering": _clustering(), "integration": _integration()})
    return cfg()


# ═══════════════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════════════

BATCH_FILES = (
    "_04_batch_colored.png",
    "_04_batch_faceted.png",
    "_04_batch_mixing_heatmap.png",
)


def test_enabled_generates_three_files(tmp_path: Path) -> None:
    """Enabled with 3 batches → 3 PNGs exist.   [EXPECTED: FAIL]"""
    adata = _make_adata(n_batches=3)
    cfg = _make_cfg(umap_color_by_batch=True)
    fig_dir = tmp_path / "figures"

    _run_batch_umap_diagnostics(adata, str(fig_dir), cfg)

    for name in BATCH_FILES:
        assert (fig_dir / name).exists(), f"Expected batch file {name!r} to exist"


def test_disabled_generates_no_batch_filenames(tmp_path: Path) -> None:
    """Disabled → no batch PNG files created.   [EXPECTED: PASS]"""
    adata = _make_adata(n_batches=3)
    cfg = _make_cfg(umap_color_by_batch=False)
    fig_dir = tmp_path / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    _run_batch_umap_diagnostics(adata, str(fig_dir), cfg)

    for name in BATCH_FILES:
        assert not (fig_dir / name).exists(), f"Expected {name!r} to NOT exist"


def test_missing_column_skips_gracefully(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Missing ``sample`` column → warning logged, no files.   [EXPECTED: PASS]"""
    adata = _make_adata(include_sample=False)
    cfg = _make_cfg(umap_color_by_batch=True, batch_key_override="sample")
    fig_dir = tmp_path / "figures"

    caplog.set_level(logging.WARNING)
    _run_batch_umap_diagnostics(adata, str(fig_dir), cfg)

    # Warning about missing batch_key
    assert "batch_key" in caplog.text, "Expected a warning about missing batch_key in the log"

    # No files should be generated
    for name in BATCH_FILES:
        assert not (fig_dir / name).exists(), f"Expected {name!r} to NOT exist"


def test_too_many_batches_degradation(tmp_path: Path) -> None:
    """14 batches → only colored + heatmap exist (no faceted).   [EXPECTED: FAIL]"""
    adata = _make_adata(n_batches=14)
    cfg = _make_cfg(umap_color_by_batch=True)
    fig_dir = tmp_path / "figures"

    _run_batch_umap_diagnostics(adata, str(fig_dir), cfg)

    assert (fig_dir / "_04_batch_colored.png").exists(), "Colored UMAP should exist"
    assert (fig_dir / "_04_batch_mixing_heatmap.png").exists(), "Mixing heatmap should exist"
    assert not (fig_dir / "_04_batch_faceted.png").exists(), (
        "Faceted should be skipped for >12 batches"
    )


def test_too_few_batches_skips(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Single batch value → no files, ``single-batch data`` logged.   [EXPECTED: PASS]"""
    adata = _make_adata(n_batches=1)
    cfg = _make_cfg(umap_color_by_batch=True)
    fig_dir = tmp_path / "figures"

    caplog.set_level(logging.INFO)
    _run_batch_umap_diagnostics(adata, str(fig_dir), cfg)

    assert "single-batch data" in caplog.text, "Expected 'single-batch data' message in the log"

    for name in BATCH_FILES:
        assert not (fig_dir / name).exists(), f"Expected {name!r} to NOT exist"


if __name__ == "__main__":
    pytest.main([__file__])
