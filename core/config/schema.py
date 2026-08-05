#!/usr/bin/env python3
"""
config.py — Fuxi (伏羲) 统一配置 (Pydantic v2)

25+ Pydantic BaseModel classes:
  23 topic sub-models + several nested + 1 top-level Config

设计原则:
  - 所有参数集中在一个 Config(BaseModel) 中
  - modality 字段区分组学类型: 'rna' | 'atac' | 'spatial'
  - 向后兼容现有项目配置文件
  - 路径自动解析在 Todo 1.2 (model_post_init) 中

使用方法:
    from core.config import Config
    CFG = Config(project_dir=...)
"""

import os
from typing import Any, Dict, List, Literal, Optional

# ── Auto-load .env from repo root ────────────────────────────────────
# This runs before any data_root() call, so FUXI_DATA_ROOT in .env
# is available to the pipeline without manual sourcing.
from dotenv import load_dotenv
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from core.config.global_config import GlobalPlotConfig

_env_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), ".env"
)

if os.path.isfile(_env_path):
    load_dotenv(_env_path)


# ── Named constants ───────────────────────────────────────────────────
SILHOUETTE_SAMPLE_THRESHOLD: int = 10000


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 1 — DataInputConfig
# ═══════════════════════════════════════════════════════════════════════
class DataInputConfig(BaseModel):
    """RNA / ATAC data input paths and format settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    mtx_prefix: str = ""
    mtx_dir: str = ""
    # Multi-sample 10X MTX: glob under mtx_dir matching sample subdirs.
    # Empty = legacy single-directory behavior (mtx_prefix applies directly).
    mtx_dir_pattern: str = ""
    # Optional regex to extract sample name from matched subdir basename.
    # Empty = use basename as-is.
    mtx_sample_regex: str = ""
    # >0: batch (tree) concat for the slow path; 0 = one-shot concat.
    mtx_concat_batch: int = 0
    matrix_file: str = ""
    metadata_file: str = ""
    barcodes_file: str = ""
    features_file: str = ""
    csv_sep: Optional[str] = None
    csv_decimal: str = "."
    gene_symbol_column: str = ""
    input_h5ad: str = ""
    backed: str = ""
    h5_file_pattern: str = "*filtered_feature_bc_matrix.h5"
    h5_dir: str = ""
    fragment_file: str = ""
    # ── Preprocessed (embedded metadata columns) format ──
    file_pattern: str = "*.tsv.gz"
    separator: str = ""  # empty = auto-detect (tab vs comma)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 2 — SampleMetaConfig
# ═══════════════════════════════════════════════════════════════════════
class SampleMetaConfig(BaseModel):
    """Sample / stage metadata mapping."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    sample_map: Dict[int | str, str] = Field(default_factory=dict)
    stage_map: Dict[int | str, str] = Field(default_factory=dict)
    stage_order: List[str] = Field(default_factory=list)
    meta_columns: Dict[str, str] = Field(default_factory=dict)
    barcode_parse_regex: str = ""
    barcode_parse_groups: Dict[str, str] = Field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 3 — QCSettings
# ═══════════════════════════════════════════════════════════════════════
class QCSettings(BaseModel):
    """Quality control thresholds."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    min_genes: int = 500
    max_genes: int = 7500
    max_pct_mito: float = 20.0
    mt_gene_pattern: str = "MT-"
    mt_gene_list: List[str] = Field(default_factory=list)
    min_genes_per_umi: float = 0.7
    min_cells_per_gene: int = 3
    use_adaptive_thresholds: bool = False
    mad_n_mads: float = 3.0
    ncount_max_mad: float = 5.0
    min_mad_upper_genes: int = 4000
    min_mad_upper_genes_nuclei: int = 3000
    is_nuclei: bool = False
    max_pct_mito_nuclei: float = 5.0
    # 流式块大小：int 固定值（如 200000），或 "auto" 按可用内存 × 0.4 / (每行nnz×12B×(prefetch+1))
    # 反推并 clamp [50k, 500k]。时间对块大小不敏感（平台期），故内存约束优先。
    block_size: int | str = "auto"

    @field_validator("block_size")
    @classmethod
    def _validate_block_size(cls, v):
        if isinstance(v, int):
            if v <= 0:
                raise ValueError("block_size must be a positive int")
            return v
        if v != "auto":
            raise ValueError('block_size must be a positive int or "auto"')
        return v

    max_pct_mito_nuclei: float = 5.0


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 4 — ScrubletSettings
# ═══════════════════════════════════════════════════════════════════════
class ScrubletSettings(BaseModel):
    """Doublet detection (Scrublet) settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    run: bool = True
    expected_doublet_rate: Optional[float] = None
    batch_key: str = "sample"
    min_counts: int = 2
    min_cells: int = 3
    min_gene_var_pctl: int = 85
    n_prin_comps: int = 30
    serial_threshold: int = 15000  # cells: >threshold → serial, ≤threshold → parallel
    # PCA/SVD backend: arpack (exact, default) vs randomized (~20-50% faster,
    # doublet labels ~97-99% consistent with arpack)
    svd_solver: Literal["arpack", "randomized"] = "arpack"
    # zscore output dtype: float32 halves zscore peak memory (~20 GiB per
    # 157k-cell group) but drifts doublet labels ~+9%; default float64 keeps
    # bit-exact parity with scrublet's original implementation
    zscore_float32: bool = False

    on_non_counts: Literal["skip_warn", "skip_silent", "abort"] = (
        "skip_warn"  # policy when expression_type != raw_counts
    )


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 4b — PearsonConfig
# ═══════════════════════════════════════════════════════════════════════
class PearsonConfig(BaseModel):
    """Pearson residuals normalization configuration."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_top_genes: int = 4000
    clip: Optional[float] = None


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 5 — NormalizationSettings
# ═══════════════════════════════════════════════════════════════════════
class NormalizationSettings(BaseModel):
    """Normalization, cell-cycle regression, sex detection."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    normalize_target_sum: float = 1e4
    use_regress_out: bool = True
    score_cell_cycle: bool = False
    regress_out_genes: List[str] = Field(default_factory=list)
    detect_sex: bool = True

    # ── Normalization method ──
    method: Literal["log_cpm", "pearson_residuals"] = "log_cpm"
    pearson_residuals: PearsonConfig = Field(default_factory=PearsonConfig)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 6 — HVGSettings
# ═══════════════════════════════════════════════════════════════════════
class HVGSettings(BaseModel):
    """Highly variable gene selection."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_top_genes: int = 4000
    flavor: str = "seurat_v3"
    batch_key: str = "sample"
    forced_genes: list[str] = Field(default_factory=list)
    auto_forced_genes: bool = Field(
        default=False,
        description="When True and forced_genes is empty, auto-fill from the tissue KB "
        "(build_forced_genes) at the 'high' consensus threshold.  Default False "
        "keeps legacy behaviour.",
    )


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 7 — PCASettings
# ═══════════════════════════════════════════════════════════════════════
class PCASettings(BaseModel):
    """Principal component analysis settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_pcs_full: int = 100
    n_pcs_use: int = 50


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 8 — SCVIConfig
# ═══════════════════════════════════════════════════════════════════════
class SCVIConfig(BaseModel):
    """scVI integration model configuration."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_latent: int = 30
    n_layers: int = 2
    n_hidden: int = 128
    max_epochs: int = 400
    batch_key: str = "sample"
    use_gpu: bool = True
    train_size: float = 0.9

    # ── Training control (high-frequency overrides) ──
    batch_size: int = 128
    early_stopping: bool = False
    precision: Literal["32", "16-mixed", "bf16-mixed"] = "32"

    # ── Passthrough dicts (advanced params via global.yaml) ──
    trainer_kwargs: dict[str, Any] = {}
    plan_kwargs: dict[str, Any] = {}
    datasplitter_kwargs: dict[str, Any] = {}


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 9 — IntegrationSettings
# ═══════════════════════════════════════════════════════════════════════
class IntegrationSettings(BaseModel):
    """Integration / batch correction settings (Harmony, scVI, Combat)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    method: Literal["harmony", "combat", "scvi"] = "harmony"
    use_harmony: bool = True
    batch_key: str = "sample"
    max_iter: int = 20

    diagnose: bool = True
    diagnose_report: bool = True
    diagnose_max_cells: int = 50_000  # batch diagnosis 子采样上限（秒级 vs 全量分钟级）
    diagnose_exclude_patterns: list[str] = Field(
        default_factory=lambda: ["*leiden*", "*cell_type*", "*annotation*", "*annotated*"]
    )
    gini_batch_threshold: float = 0.3
    gini_biology_threshold: float = 0.6
    collinearity_guard: bool = True
    # 流式写 .raw：03 写盘时从 02_qc.h5ad 分块读取 counts → normalize+log1p
    # → 直写输出文件 raw 组（省内存峰值，避免全基因矩阵常驻）。
    # 默认 False（保持原行为）；大数据集（>100 万细胞）建议开启。
    stream_raw: bool = False
    scvi: SCVIConfig = Field(default_factory=SCVIConfig)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model — AmbientSettings
# ═══════════════════════════════════════════════════════════════════════
class AmbientSettings(BaseModel):
    """Ambient RNA correction / removal settings (CellBender, SoupX)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    run: bool = False
    method: Literal["cellbender", "soupx", "none"] = "none"
    raw_matrix_path: str = ""
    expected_cells: int = 0
    fallback_fraction: float = 0.1


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 10 — ClusteringSettings
# ═══════════════════════════════════════════════════════════════════════
class ClusteringSettings(BaseModel):
    """Clustering, UMAP, and parameter grid search settings."""

    model_config = ConfigDict(extra="ignore", validate_assignment=True)

    n_neighbors: int = 30  # single-value KNN for non-clustering uses (trajectory, ATAC)
    use_rep: str = Field(
        default="X_pca",
        description="Representation in adata.obsm used for nearest-neighbor graph and clustering.",
    )
    param_grid_n_neighbors: list = Field(default_factory=lambda: [15, 20, 30])
    param_grid_n_neighbors_adaptive: bool = Field(
        default=True,
        description="When True, auto-compute param_grid_n_neighbors based on adata.n_obs (adaptive scaling). When False, use the explicit param_grid_n_neighbors value.",
    )
    param_grid_resolutions: list = Field(default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    leiden_flavor: str = "igraph"
    leiden_n_iterations: int = Field(
        default=2,
        description="Leiden n_iterations for main clustering paths (grid search, annotate, subcluster, batch diagnostics). Use -1 for full convergence.",
    )
    best_resolution: float = 1.0
    best_n_neighbors: int = 0
    cluster_selection_method: Optional[str] = "multi_metric"
    multi_metric_weights: dict = Field(
        default_factory=lambda: {
            "silhouette": 0.2,
            "stability": 0.2,
            "cluster_coherence": 0.3,
            "splitting_gain": 0.2,
            "kb_annotatable_rate": 0.1,
        }
    )
    stability_n_seeds: int = Field(
        default=12,
        description="Number of random seeds for stability evaluation in multi-metric clustering.",
    )
    multi_metric_stability_top_k: int = Field(
        default=4,
        description="Only the top-K entries (by silhouette) per n_neighbors group get the expensive 5-seed stability scan; the rest get stability_score=None (treated as 0 in multi-metric). Measured: stability is 50-77% of Step 04 wall time; n_neighbors is nearly flat vs resolution.",
    )
    stability_leiden_n_iterations: int = Field(
        default=-1,
        description="Leiden n_iterations for the stability sub-path in multi-metric selection. -1 = full convergence for accurate ARI.",
    )
    leiden_gpu_min_cells: int = Field(
        default=20_000,
        description="Below this n_obs, Leiden runs on CPU (igraph) even when a GPU is present. cuGraph's graph-build + anndata_to_GPU copy overhead dominates small datasets (measured: 5-20k cells GPU is 1.5-3x slower than CPU).",
    )
    multi_metric_adaptive_resolution: bool = True
    multi_metric_coverage_ratio_threshold: float = 2.5
    multi_metric_granularity_cv_threshold: float = 0.05
    multi_metric_granularity_min_clusters: int = 10
    multi_metric_de_gate_threshold: int = 25
    umap_selection_method: Optional[str] = "convex_hull"
    umap_selection_metric: Literal["trustworthiness", "convex_hull"] = "trustworthiness"
    umap_paga_init: bool = Field(
        default=False,
        description="Use PAGA-initialized UMAP positions during parameter selection.",
    )
    # min_dist sweep is an O(n^2) trustworthiness cost and makes silhouette
    # non-comparable across grid combos (different embedding geometry).
    # Default: single value (production, reproducible). Multi-value grids are
    # a research-mode escape hatch — step 04 degrades to the first value with a
    # warning when n_obs > 30k (76k cells ≈ 92GB pairwise matrix).
    param_grid_min_dist: Optional[list] = Field(default_factory=lambda: [0.3])
    param_grid_spread: Optional[list] = Field(default_factory=lambda: [1.0])
    umap_min_dist: float = 0.5
    de_pairwise_max_clusters: int = Field(
        default=30,
        description="Max clusters for pairwise DE in _select_de_gated. Above this, falls back to one-vs-rest for performance. 0 = always pairwise (no cap).",
    )
    umap_spread: float = 1.0
    umap_maxiter: Optional[int] = Field(
        default=None,
        description="Max UMAP iterations (None = auto/unlimited). Set to 200-500 on large datasets to bound convergence time.",
    )
    umap_n_epochs: Optional[int] = Field(
        default=None,
        description="Number of UMAP training epochs (None = auto). Lower values (e.g., 500) speed up at slight quality cost.",
    )
    umap_color_by_batch: bool = False
    batch_key_override: Optional[str] = None
    # Performance: matplotlib scatter on >100k cells is the dominant cost in
    # step 04 (Li2026 1M cells spent ~4h just drawing UMAP scatter). Default
    # 'auto' uses full rendering under the threshold and subsamples above it,
    # preserving backward compatibility for small datasets while making large
    # datasets actually finish. 'skip' disables plotting entirely.
    umap_plot_mode: Literal["auto", "full", "subsample", "skip"] = "auto"
    umap_plot_max_cells: int = 50000
    # Stratified scatter subsample caps (step 04 figures): label sizes <=
    # cap_small are kept in full (rare populations stay visible); 500-50k
    # labels capped at cap_medium; > umap_plot_stratum_large at cap_large.
    # The hard ceiling remains umap_plot_max_cells.
    umap_plot_cap_small: int = 500
    umap_plot_cap_medium: int = 500
    umap_plot_cap_large: int = 1000
    umap_plot_stratum_large: int = 50000
    # plot_per_combo: skip per-(n_neighbors, resolution) UMAP subplot loop.
    # When False, skips the per-combo individual UMAP scatter plots (the dominant
    # matplotlib cost at scale — Li2026 1M cells: ~2h for this loop alone).
    # Summary grid figure, best-params UMAP, batch-colored UMAP, and other
    # aggregate plots still generate regardless of this flag.
    plot_per_combo: bool = Field(
        default=False,
        description="Plot per-(n_neighbors, resolution) UMAP subplots. Disable to save ~80% plot wall time at 1M scale.",
    )

    # ── Funnel mode ──
    funnel_enabled: bool = Field(
        default=False,
        description="Progressive funnel grid search (subsample grid + full re-validation) for datasets > funnel_threshold. 2026-08-05 A/B: full grid now faster (134k 182s vs 6613s; 1.05M 907s vs 6731s) and unbiased — funnel subsampling changed best param on both test sets. Keep off by default; escape hatch only for >2M-cell datasets.",
    )
    funnel_threshold: int = Field(
        default=100_000, description="Minimum n_obs to trigger funnel mode."
    )
    funnel_subsample_size: int = Field(
        default=50_000, description="Subsample size for funnel grid search."
    )
    funnel_top_k: int = Field(
        default=3, description="Number of top candidates to re-validate on full data."
    )
    funnel_kmeans_dim: int | None = Field(
        default=None,
        description="KMeans dims for stratified subsampling (None=full embedding; dim reduction rejected after measurement: 100->50 changes 70% of picks for only 2x).",
    )
    funnel_kmeans_n_init: int = Field(
        default=1,
        description="KMeans init runs for stratified subsampling (single init suffices for coverage).",
    )

    # ── Target mode ──
    target_n_clusters: int | None = Field(
        default=None,
        description="Target cluster count. When set, binary-search resolution instead of grid.",
    )
    target_search_max_iters: int = Field(
        default=10, description="Max iterations for binary search in target mode."
    )


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 11 — MarkerSettings
# ═══════════════════════════════════════════════════════════════════════
class MarkerSettings(BaseModel):
    """Cell-type marker / annotation settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    marker_dict: Dict[str, List[str]] = Field(default_factory=dict)
    subcluster_types: List[str] = Field(default_factory=list)
    subcluster_resolution: float = 0.4
    min_cells_subcluster: int = 50
    expert_rule_strictness: str = "default"
    expert_rule_top_n: int = 0
    candidate_pool_expand_steps: List[int] = Field(default_factory=lambda: [50, 100, 200])
    validation_n_top_genes: int = 15
    validation_min_overlap: float = 0.5
    validation_marginal_threshold: float = 0.25
    validation_pass_rate_min: float = 0.1
    quality_gate_min_pass_rate: float = 0.10
    developmental_mode: bool = False
    step10_groupby: List[str] = Field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 12 — DESettings
# ═══════════════════════════════════════════════════════════════════════
class DESettings(BaseModel):
    """Differential expression analysis settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    method: Literal["wilcoxon", "pseudobulk"] = "pseudobulk"
    n_genes: int = 50
    pval_cutoff: float = 0.05
    logfc_cutoff: float = 0.25
    pairwise_method: Literal["wilcoxon", "t-test"] = "wilcoxon"
    branch_method: Literal["wilcoxon", "t-test"] = "t-test"
    stage_pairwise: bool = True
    use_raw: bool = True
    auto_switch_on_low_quality: bool = False
    pseudobulk: "PseudobulkDESettings" = Field(default_factory=lambda: PseudobulkDESettings())


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 12b — PseudobulkDESettings (nested under de)
# ═══════════════════════════════════════════════════════════════════════
class PseudobulkDESettings(BaseModel):
    """Pseudobulk differential expression via PyDESeq2."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    celltype_col: str = "cell_type"
    sample_col: str = "sample"
    design: str = "~batch + condition"
    contrast_column: str = "condition"
    contrast_treatment: str = ""
    contrast_baseline: str = ""
    alpha: float = 0.05
    min_cells_per_sample: int = 10
    min_cells_per_group: int = 3
    lfc_shrink: bool = True
    n_jobs: int = 0
    output_dir: str = ""


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 12c — ScVeloConfig
# ═══════════════════════════════════════════════════════════════════════
class ScVeloConfig(BaseModel):
    """scVelo RNA velocity configuration."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    loom_path: str = ""
    mode: Literal["stochastic", "dynamical"] = "stochastic"


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 13 — TrajectorySettings
# ═══════════════════════════════════════════════════════════════════════
class TrajectorySettings(BaseModel):
    """Pseudotime / trajectory analysis settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    method: Literal["paga_dpt", "scvelo_cellrank"] = "paga_dpt"
    scvelo: ScVeloConfig = Field(default_factory=ScVeloConfig)

    root_cell_types: List[str] = Field(default_factory=list)
    root_markers: List[str] = Field(default_factory=list)
    n_diffmap_comps: int = 15
    n_branchings: int = 2
    pseudotime_genes: List[str] = Field(default_factory=list)
    pseudotime_n_branch_de: int = 10
    pseudotime_n_correlated: int = 10
    pseudotime_cor_pval: float = 0.05
    save_final_h5ad: bool = Field(
        default=True,
        description="是否输出 05_final.h5ad。该文件全仓库无下游消费者，允许关闭以减少中间产物占用。",
    )


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 14 — EnrichmentSettings
# ═══════════════════════════════════════════════════════════════════════
class EnrichmentSettings(BaseModel):
    """Gene-set enrichment analysis settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    run: bool = True
    method: str = "both"
    gene_sets: list = Field(
        default_factory=lambda: [
            "GO_Biological_Process_2023",
            "KEGG_2021_Human",
        ]
    )
    organism: str = "human"
    n_top_genes: int = 200
    pval_cutoff: float = 0.05
    min_size: int = 10
    max_size: int = 500
    permutations: int = 1000
    tissue_mode: str = "off"
    tissue_pathways_whitelist: list = Field(default_factory=list)
    tissue_pathways_blacklist: list = Field(default_factory=list)
    redundancy_cluster: bool = False
    redundancy_threshold: float = 0.6
    use_kb_relevance: bool = False
    gene_sets_tissue: list = Field(default_factory=list)
    background_restrict: bool = False
    peak_gene_distance: int = 100000
    gene_annotation_bed: str = ""


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 15 — GRNSettings
# ═══════════════════════════════════════════════════════════════════════
class GRNSettings(BaseModel):
    """Gene regulatory network analysis settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    run: bool = True
    method: str = "decoupler"
    species: str = "human"
    n_top_regulons: int = 50
    min_regulon_size: int = 5
    confidence_levels: list = Field(default_factory=lambda: ["A", "B", "C"])
    tissue_mode: str = "off"
    use_kb_relevance: bool = False
    export_filtered: bool = False


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 16 — CCISettings
# ═══════════════════════════════════════════════════════════════════════
class CCISettings(BaseModel):
    """Cell-cell interaction analysis settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    run: bool = True
    method: str = "liana"
    lr_database: str = "consensus"
    permutations: int = 100
    n_top_interactions: int = 50
    spatial_method: str = "liana_spatial"
    spatial_distance: float = 0.0
    lr_cache_dir: str = ""
    adjacency: str = "off"
    tissue: str = ""
    adjacency_file: str = ""
    adjacency_types: list = Field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 17 — DownsampleSettings
# ═══════════════════════════════════════════════════════════════════════
class DownsampleSettings(BaseModel):
    """Downsampling and subset filtering settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    target: Optional[int] = None
    strategy: str = "stratified"
    max_per_sample: Optional[int] = None
    random_seed: int = 42
    sample_keep: List[str] = Field(default_factory=list)
    obs_filter: str = ""
    subset_suffix: str = ""


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 18 — SpatialConfig
# ═══════════════════════════════════════════════════════════════════════
class SpatialConfig(BaseModel):
    """Spatial transcriptomics platform and processing settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    platform: str = "visium"
    library_id: str = ""
    img_path: str = ""
    spot_diameter: float = 0.0
    crop_image: bool = True
    img_rescale: float = 1.0
    neighbors_n: int = 6
    neighbors_radius: float = 0.0
    run_autocorr: bool = True
    moran_percentile: int = 90
    svg_n_top: int = 2000
    run_segmentation: bool = False


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 19 — ATACConfig
# ═══════════════════════════════════════════════════════════════════════
class ATACConfig(BaseModel):
    """ATAC-specific configuration fields."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    genome: str = "hg38"
    chrom_sizes: str = ""
    blacklist_bed: str = ""
    tss_bed: str = ""
    min_fragments: int = 1000
    max_fragments: int = 50000
    min_tsse: float = 7.0
    max_blacklist_ratio: float = 0.05
    min_peak_region_fragments: int = 300
    peak_qval: float = 0.05
    peak_width: int = 500
    use_macs3: bool = True
    n_features: int = 50000
    n_spectral: int = 30
    marker_peaks_log2fc: float = 0.5
    marker_peaks_fdr: float = 0.05
    motif_db: str = "JASPAR2024"
    terminal_cell_types: List[str] = Field(default_factory=list)
    max_cells: Optional[int] = None
    harmony_use_harmony: bool = False
    harmony_batch_key: str = "sample"
    multi_metric_enabled: bool = False


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 20 — ExecutionConfig
# ═══════════════════════════════════════════════════════════════════════


class MemoryConfig(BaseModel):
    """Memory budget / policy / guard rail for steps 01-03."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    # speed   = no tradeoffs (dense PCA, full regress_out, fastest)
    # balanced = avoid dense copies (arpack PCA, skip regress_out)
    # memory  = maximum savings (arpack + subsample UMAP train)
    policy: Literal["speed", "balanced", "memory"] = "speed"
    budget: str = "auto"  # auto = psutil detect (80% phys RAM) | e.g. 64GB | 128GiB | 32000MB
    guard: Literal["warn", "block", "off"] = "warn"  # pre-run peak estimate vs budget


class ExecutionConfig(BaseModel):
    """Execution environment settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_jobs: int = 0
    limit_blas_threads: bool = True
    random_seed: int = 42
    scanpy_verbosity: int = 2
    force_csr: bool = True
    use_float32: bool = True
    # ── Device selection (GPU acceleration via rapids-singlecell) ──
    # auto = detect at runtime, fall back to CPU if RAPIDS unavailable
    # gpu  = force GPU, raise on missing RAPIDS
    # cpu  = force CPU (skip detection entirely)
    device: Literal["auto", "cpu", "gpu"] = "auto"
    # ── Memory budget / policy / guard (steps 01-03) ──
    memory: MemoryConfig = Field(default_factory=MemoryConfig)
    # ── Statistical approximation (step 05 fast mode) ──
    # orthogonal to use_float32: that flag trades numeric dtype (ulp-level
    # representation), while approximation trades statistical precision
    # (downsampled cluster statistics are estimates of the full-data ones).
    approximation: Literal["exact", "fast"] = "exact"
    fast_sampling: int = Field(default=5000, ge=100)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 21 — AIConfig
# ═══════════════════════════════════════════════════════════════════════
class AIConfig(BaseModel):
    """AI / LLM configuration — all AI features controlled here."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    enabled: bool = False
    api_base: str = ""
    model: str = "deepseek-v4-flash"
    api_key: str = ""
    max_tokens: int = 4096
    temperature: float = 0.1
    thinking_enabled: bool = True
    reasoning_effort: str = "high"
    timeout: Optional[int] = None
    qc_review: bool = False
    param_suggest: bool = False
    annotation: bool = True
    ai_annotation: bool = False  # AI fallback for low-confidence clusters (Unified KB mode only)
    subcluster: bool = True
    deg_design: bool = False
    interpretation: bool = True
    cache_responses: bool = True
    unconstrained_annotation: bool = False
    subcluster_kb_constrained: bool = True
    """Constrain Step 06 AI subcluster naming to KB subtype space when a hierarchy exists."""


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 22 — BulkConfig
# ═══════════════════════════════════════════════════════════════════════
class BulkConfig(BaseModel):
    """Bulk RNA-seq specific configuration fields."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    design: str = "~condition"
    contrast_column: str = "condition"
    contrast_treatment: str = ""
    contrast_baseline: str = ""
    alpha: float = 0.05
    lfc_shrink: bool = True
    normalization_method: str = "deseq2_median_ratios"
    min_counts_per_gene: int = 10
    min_samples_per_group: int = 2
    n_jobs: int = 0
    output_dir: str = ""
    batch_correct: bool = False
    batch_column: str = "batch"


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 23 — ScCODAConfig (nested under exploratory)
# ═══════════════════════════════════════════════════════════════════════
class ScCODAConfig(BaseModel):
    """scCODA compositional analysis configuration."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    sample_col: str = ""
    condition_col: str = ""
    reference_cell_type: str = ""


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 24 — ExploratorySettings
# ═══════════════════════════════════════════════════════════════════════
class ExploratorySettings(BaseModel):
    """Exploratory analysis settings (composition test, etc.)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    composition_test: Literal["none", "sccoda"] = "none"
    sccoda: ScCODAConfig = Field(default_factory=ScCODAConfig)


# ═══════════════════════════════════════════════════════════════════════
# CellTypistConfig — nested sub-model for AnnotationSettings
# ═══════════════════════════════════════════════════════════════════════
class CellTypistConfig(BaseModel):
    """CellTypist annotation model configuration."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    model: str = ""
    """Model name (e.g. \"Immune_All_Low.pkl\"). Empty = auto-select disabled."""
    majority_voting: bool = True
    enabled: bool = False


# ═══════════════════════════════════════════════════════════════════════
# AnnotationSettings — annotation method selector + supplementary evidence
# ═══════════════════════════════════════════════════════════════════════
class AnnotationSettings(BaseModel):
    """Annotation method and supplementary evidence settings."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    method: Literal["kb_unified", "celltypist", "ai", "score_genes"] = "kb_unified"
    celltypist: CellTypistConfig = Field(default_factory=CellTypistConfig)
    multi_peak_min_types: int = 3
    """多峰歧义降级阈值:过渡路径上 ≥N 个类型 score≥floor 即降级 ambiguous。
    (D1) min_types=3 而非 4:3 型并列与 4 型并列同属 Top-2 任意截断失败模式,
    第 3 名会被 _is_transition_state 静默丢弃。"""
    multi_peak_score_floor: float = 0.9
    """多峰歧义降级分数下限:score≥floor 的类型计入并列集合。"""

    canonical_pct_floor: float = 0.05
    """canonical 表达兜底的 pct 表达下限 (D3):融合后, confidence ∈ (high, medium)
    的 marker_scoring 决策, winning 类型的 top-consensus 标记 (consensus≥2, 取前 3)
    在簇中 pct 表达全部低于此值 → 降级 confidence=low + 计入 review_queue
    (reason=no_canonical_expression)。"""

    # ── Layer-3 KADP developmental potency (plan annotation-kadp-metc todo 5) ──
    # Thresholds default off; the exact values are locked by calibration
    # (todo 6/7) on GSE246169 fetal.
    kadp_enabled: bool = False
    """Enable the KADP potency axis: ambiguous multi-peak clusters are named as
    differentiating precursors when the Progenitor pole dominates."""
    kadp_ratio_threshold: float = 2.0
    """KADP ratio variant pass threshold: max_prog / max(max_term, epsilon)."""
    kadp_abs_threshold: float = 0.6
    """KADP abs variant pass threshold: max_prog (with max_prog > max_term guard)."""
    kadp_gap_threshold: float = 0.1
    """KADP gap variant pass threshold: max_prog - max_term."""
    use_gap_criterion: bool = False
    """When True, the gap variant joins the ratio/abs pass combination."""

    # ── Layer-4 METC multi-source transition voting (plan annotation-kadp-metc todo 10) ──
    # Default off; thresholds mirror METCConfig in rna/utils/evidence_fusion.py.
    metc_enabled: bool = False
    """Enable METC: multi-source voting (marker / expert / AI / CellTypist)
    arbitrates ambiguous/transition_state candidates that KADP did not name.
    Requires ``allows_transitions`` context (developing tissue or _dev_mode)."""
    metc_min_sources: int = 3
    """Minimum number of speaking evidence sources for METC arbitration."""
    metc_min_distinct_transition: int = 3
    """Minimum number of distinct votes to emit a transitional (vs a 2-way
    ambiguous split) decision."""


# ═══════════════════════════════════════════════════════════════════════
# Top-level Config
# ═══════════════════════════════════════════════════════════════════════
class Config(BaseModel):
    """Fuxi unified config — contains all RNA + ATAC + Spatial fields."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    # ═══════════════════════════════════════════════════════════════════
    # 组学类型 & 元信息
    # ═══════════════════════════════════════════════════════════════════
    modality: str = "rna"
    tissue: str = "unknown"
    species: str = "human"
    tissue_maturity: str = "unknown"
    expression_type: str = "raw_counts"
    data_format: str = "10X_mtx"

    # ═══════════════════════════════════════════════════════════════════
    # 路径设置
    # ═══════════════════════════════════════════════════════════════════
    data_dir: str = ""
    results_dir: str = "results"
    h5ad_dir: str = "results/h5ad"
    figure_dir: str = "results/figures"
    table_dir: str = "results/tables"
    log_dir: str = "logs"
    project_dir: str = ""

    # ═══════════════════════════════════════════════════════════════════
    # 执行环境相关（留在根级，与 ExecutionConfig 配合使用）
    # ═══════════════════════════════════════════════════════════════════
    h5ad_compression: str = "gzip"
    per_step_h5ad_compression: dict[str, str] = Field(
        default_factory=lambda: {"integrated": "gzip"}
    )
    h5ad_compression_opts: Optional[int] = Field(
        default=None,
        description="h5py 压缩级别（如 gzip level 1=快速、文件略大；level 4=默认平衡；level 9=最小文件、最慢）。None 时用 anndata/h5py 默认。仅对 gzip 压缩生效。大数据集（>100 万细胞）建议 level 1 提速写盘。",
    )
    h5ad_tempdir: str = Field(
        default="/tmp/Fuxi",
        description="Temporary directory for h5ad writes. Default /tmp/Fuxi uses tmpfs (RAM-backed, fast but volatile). On WSL2 with /mnt/ data, cross-filesystem rename costs ~1.5-2x. For large datasets (>4 GB), set to the same filesystem as the output directory to avoid cross-device copy overhead.",
    )
    cleanup_intermediates: bool = False
    perf_monitoring: bool = True
    verify_write_integrity: bool = True
    incremental_io: bool = Field(
        default=True,
        description="增量写入开关（默认开）。True → 步骤对 obs/obsm/obsp/uns 走 in-place 追加写；False → 全部回退全量 safe_write。默认开的理由：copy+append 崩溃即删源文件兜底、in-place 写回有 .bak 兜底；逃生口：WSL /mnt 文件系统不稳时置 False 回退全量写入。",
    )

    # ═══════════════════════════════════════════════════════════════════
    # ATAC → RNA 整合 / Spatial RNA 参考
    # ═══════════════════════════════════════════════════════════════════
    rna_h5ad: str = ""
    rna_ref: str = ""
    rna_marker_top_n: int = 10
    rna_marker_pval_threshold: float = 0.05
    rna_marker_logfc_min: float = 0.0

    # ═══════════════════════════════════════════════════════════════════
    # 组织知识库 (tissue_kb / tissue_ontology 保留在根级，由下游模块使用)
    # ═══════════════════════════════════════════════════════════════════
    tissue_kb: str = ""
    tissue_ontology: str = ""

    # ═══════════════════════════════════════════════════════════════════
    # 系统发育过滤 (phylogenetic filtering for KB)
    # ═══════════════════════════════════════════════════════════════════
    @model_validator(mode="after")
    def _normalize_species_validator(self):
        """Defence-in-depth: normalise species to canonical pipeline key.

        The canonical normalisation happens in ``resolve_config`` (single
        source of truth), but this validator catches any path that
        instantiates Config directly (tests, scripts, ad-hoc usage).
        """
        from core.preprocess.format_detector import _SPECIES_NORMALISE

        raw = self.species
        norm = _SPECIES_NORMALISE.get(raw)
        if norm is None:
            norm = _SPECIES_NORMALISE.get(raw.lower(), raw)
        if norm != raw:
            object.__setattr__(self, "species", norm)
        return self

    target_class: str = ""
    target_order: str = ""

    # ═══════════════════════════════════════════════════════════════════
    # 23 个主题子模型
    # ═══════════════════════════════════════════════════════════════════
    data_input: DataInputConfig = Field(default_factory=DataInputConfig)
    sample_meta: SampleMetaConfig = Field(default_factory=SampleMetaConfig)
    qc: QCSettings = Field(default_factory=QCSettings)
    scrublet: ScrubletSettings = Field(default_factory=ScrubletSettings)
    normalization: NormalizationSettings = Field(default_factory=NormalizationSettings)
    hvg: HVGSettings = Field(default_factory=HVGSettings)
    pca: PCASettings = Field(default_factory=PCASettings)
    integration: IntegrationSettings = Field(default_factory=IntegrationSettings)
    clustering: ClusteringSettings = Field(default_factory=ClusteringSettings)
    marker: MarkerSettings = Field(default_factory=MarkerSettings)
    de: DESettings = Field(default_factory=DESettings)
    trajectory: TrajectorySettings = Field(default_factory=TrajectorySettings)
    enrichment: EnrichmentSettings = Field(default_factory=EnrichmentSettings)
    exploratory: ExploratorySettings = Field(default_factory=ExploratorySettings)
    grn: GRNSettings = Field(default_factory=GRNSettings)
    cci: CCISettings = Field(default_factory=CCISettings)
    downsample: DownsampleSettings = Field(default_factory=DownsampleSettings)
    spatial: SpatialConfig = Field(default_factory=SpatialConfig)
    atac: ATACConfig = Field(default_factory=ATACConfig)
    execution: ExecutionConfig = Field(default_factory=ExecutionConfig)
    ai: AIConfig = Field(default_factory=AIConfig)
    bulk: BulkConfig = Field(default_factory=BulkConfig)
    plot: GlobalPlotConfig = Field(default_factory=GlobalPlotConfig)
    ambient: AmbientSettings = Field(default_factory=AmbientSettings)
    annotation: AnnotationSettings = Field(default_factory=AnnotationSettings)

    def model_post_init(self, __context):
        """Resolve relative paths after construction.

        Replaces the old resolve_paths() with Pydantic-native
        post-init hook.
        """
        base = self.project_dir if self.project_dir else os.path.dirname(os.path.abspath(__file__))

        # Treat '.' as "not set" for mtx_dir
        if self.data_input.mtx_dir == ".":
            self.data_input.mtx_dir = ""

        # Resolve top-level relative paths to absolute
        for attr in (
            "data_dir",
            "results_dir",
            "h5ad_dir",
            "figure_dir",
            "table_dir",
            "log_dir",
        ):
            val = getattr(self, attr)
            if val and not os.path.isabs(val):
                setattr(self, attr, os.path.join(base, val))

        # Resolve sub-model relative paths to absolute
        sub = self.data_input
        if sub.mtx_dir and not os.path.isabs(sub.mtx_dir):
            sub.mtx_dir = os.path.join(base, sub.mtx_dir)
        if sub.h5_dir and not os.path.isabs(sub.h5_dir):
            sub.h5_dir = os.path.join(base, sub.h5_dir)

        # Auto-resolve data_dir from FUXI_DATA_ROOT when empty
        if not self.data_dir:
            _data_root = os.environ.get("FUXI_DATA_ROOT") or os.environ.get("SCRNA_DATA_ROOT")
            if _data_root:
                dataset_id = os.path.basename(self.project_dir or base)
                self.data_dir = os.path.join(_data_root, dataset_id)
            else:
                self.data_dir = base

        # Resolve individual data files relative to data_dir
        if sub.matrix_file and not os.path.isabs(sub.matrix_file):
            sub.matrix_file = os.path.join(self.data_dir, sub.matrix_file)
        if sub.barcodes_file and not os.path.isabs(sub.barcodes_file):
            sub.barcodes_file = os.path.join(self.data_dir, sub.barcodes_file)
        if sub.features_file and not os.path.isabs(sub.features_file):
            sub.features_file = os.path.join(self.data_dir, sub.features_file)
        if sub.fragment_file and not os.path.isabs(sub.fragment_file):
            sub.fragment_file = os.path.join(self.data_dir, sub.fragment_file)
        if sub.input_h5ad and not os.path.isabs(sub.input_h5ad):
            sub.input_h5ad = os.path.join(self.data_dir, sub.input_h5ad)
        if sub.metadata_file and not os.path.isabs(sub.metadata_file):
            sub.metadata_file = os.path.join(self.data_dir, sub.metadata_file)

        # Auto-fill mtx_dir and h5_dir from data_dir
        if not self.data_input.mtx_dir:
            self.data_input.mtx_dir = self.data_dir
        if not self.data_input.h5_dir:
            self.data_input.h5_dir = self.data_dir

        # Subset filter: auto-append suffix to output dirs
        ds = self.downsample
        if ds.sample_keep or (ds.obs_filter and ds.obs_filter.strip()):
            suffix = ds.subset_suffix if ds.subset_suffix else "_subset"
            self.h5ad_dir = self.h5ad_dir.rstrip("/\\") + suffix
            self.figure_dir = self.figure_dir.rstrip("/\\") + suffix
            self.table_dir = self.table_dir.rstrip("/\\") + suffix
            self.log_dir = self.log_dir.rstrip("/\\") + suffix
            print(f"[Config] Subset active → output dir suffix: '{suffix}'")

        # tissue_kb auto-inference from tissue
        if not self.tissue_kb and self.tissue not in ("unknown", ""):
            _repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            _kb_dir = os.path.join(_repo_root, "kb", self.tissue)
            if os.path.isdir(_kb_dir):
                self.tissue_kb = self.tissue
                self.tissue_ontology = self.tissue_ontology or self.tissue

    # ═══════════════════════════════════════════════════════════════════
    # RNA checkpoint 路径（属性）
    # ═══════════════════════════════════════════════════════════════════
    @property
    def raw_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "00_raw.h5ad")

    @property
    def qc_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "02_qc.h5ad")

    @property
    def doublet_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "01_doublet.h5ad")

    @property
    def norm_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_normalized.h5ad")

    @property
    def integrated_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_integrated.h5ad")

    @property
    def cluster_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "04_clustered.h5ad")

    @property
    def annotated_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "05_annotated.h5ad")

    @property
    def final_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "05_final.h5ad")

    @property
    def grn_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "11_grn.h5ad")

    # ── ATAC: checkpoint 路径 ──
    @property
    def filtered_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "02_filtered.h5ad")

    @property
    def processed_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_processed.h5ad")

    @property
    def clustered_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "04_clustered.h5ad")

    @property
    def trajectory_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "10_trajectory.h5ad")

    @property
    def peak_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "05_peaks.h5ad")

    # ── Spatial: checkpoint paths ──
    @property
    def sq_image_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "02_image.h5ad")

    @property
    def sq_processed_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_processed.h5ad")

    # ═══════════════════════════════════════════════════════════════════
    #  方法
    # ═══════════════════════════════════════════════════════════════════
    def has_sample_mapping(self) -> bool:
        return len(self.sample_meta.sample_map) > 0

    def has_stage_mapping(self) -> bool:
        return len(self.sample_meta.stage_map) > 0

    def has_markers(self) -> bool:
        return len(self.marker.marker_dict) > 0

    def has_rna_data(self) -> bool:
        """ATAC: check if RNA data is available for integration"""
        return bool(self.rna_h5ad) and os.path.exists(self.rna_h5ad)
