#!/usr/bin/env python3
"""Global visualization configuration defaults — shared across all modalities.

Contains GlobalPaletteConfig (color palettes) and GlobalPlotConfig (DPI, size,
format) that serve as the single source of truth for all plot rendering in
the Fuxi pipeline.

These models are designed to be overridable per project via global.yaml;
every field has a sensible default so no project is required to specify them.
"""

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
