#!/usr/bin/env python3
"""Project-wide global configuration defaults — shared across all modalities.

Contains GlobalPaletteConfig, GlobalPlotConfig, Execution/DE/Clustering/
Integration/QC sub-configs, and a top-level GlobalConfig that aggregates them
all.  Every field has a sensible default so no project is required to specify
them; projects override only what differs via global.yaml.
"""

from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class GlobalPaletteConfig(BaseModel):
    """Matplotlib/seaborn color palette defaults for every plot type."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    categorical: str = "tab20"
    heatmap: str = "RdBu_r"
    pseudotime: str = "plasma"
    dotplot_fill: str = "YlOrRd"
    gsea_dotplot: str = "RdBu_r"
    batch_heatmap: str = "YlOrRd"
    composition: str = "tab20"
    interaction_heatmap: str = "RdYlBu_r"
    qc_hist: str = "steelblue"
    qc_second: str = "indianred"
    qc_third: str = "darkorange"
    qc_threshold: str = "red"
    batch_gini_good: str = "#2ecc71"
    batch_gini_bad: str = "#e74c3c"
    batch_gini_ambiguous: str = "#f39c12"
    significance_edge: str = "grey"
    grn_facecolor: str = "white"


class GlobalPlotConfig(BaseModel):
    """Global visualization defaults controlling figure output.

    Every field has a sensible default so projects only override what they
    need in their global.yaml.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    figure_dpi: int = 150
    figure_format: str = "pdf"
    figure_transparent: bool = True
    scatter_dot_size: int = 4
    legend_fontsize: int = 8
    legend_fontsize_dense: int = 5
    title_fontsize: int = 12
    axis_label_fontsize: int = 10
    qc_figure_size: list[int] = [8, 5]
    umap_panel_size: list[int] = [6, 5]
    palette: GlobalPaletteConfig = Field(default_factory=GlobalPaletteConfig)


class GlobalExecutionConfig(BaseModel):
    """Execution environment defaults — parallelization, memory, float precision.

    These values control runtime behaviour of every pipeline step across
    all modalities. Projects override only what differs from the common
    baseline in their global.yaml.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    random_seed: int = 42
    memory: dict = {"policy": "speed", "budget": "auto", "guard": "warn"}
    n_jobs: int = 0  # 0 = auto-detect
    limit_blas_threads: bool = True
    force_csr: bool = True
    use_float32: bool = True


class GlobalDEConfig(BaseModel):
    """Differential expression analysis defaults — method, gene count, cutoffs."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    method: str = "wilcoxon"
    n_genes: int = 50
    pval_cutoff: float = 0.05
    logfc_cutoff: float = 0.25


class GlobalClusteringConfig(BaseModel):
    """Clustering and UMAP defaults — selection method, parameter grids.

    Controls how the pipeline picks cluster resolution (multi-metric) and the
    UMAP layout (convex-hull) without requiring per-project config.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    cluster_selection_method: Optional[str] = "multi_metric"
    umap_selection_method: str = "convex_hull"
    umap_selection_metric: str = "trustworthiness"
    param_grid_spread: list[float] = [1.0]
    umap_min_dist: float = 0.3
    umap_spread: float = 1.0


class GlobalIntegrationConfig(BaseModel):
    """Integration / batch-correction defaults — Harmony, diagnosis thresholds."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    max_iter: int = 20
    diagnose: bool = True
    diagnose_report: bool = True
    gini_batch_threshold: float = 0.3
    gini_biology_threshold: float = 0.6


class GlobalQCConfig(BaseModel):
    """Quality-control thresholds — mito, gene detection, expression pattern."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    max_pct_mito: float = 20.0
    min_cells_per_gene: int = 3
    mt_gene_pattern: str = "MT-"


class GlobalConfig(BaseModel):
    """Top-level global config aggregating every sub-config namespace.

    This is the single import consumers use to access all project-wide
    defaults.  Any sub-config can be overridden via global.yaml without
    touching the remainder.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    execution: GlobalExecutionConfig = Field(default_factory=GlobalExecutionConfig)
    de: GlobalDEConfig = Field(default_factory=GlobalDEConfig)
    clustering: GlobalClusteringConfig = Field(default_factory=GlobalClusteringConfig)
    integration: GlobalIntegrationConfig = Field(default_factory=GlobalIntegrationConfig)
    qc: GlobalQCConfig = Field(default_factory=GlobalQCConfig)
    plot: GlobalPlotConfig = Field(default_factory=GlobalPlotConfig)
