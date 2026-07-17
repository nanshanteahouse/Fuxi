#!/usr/bin/env python3
"""
config.py — Fuxi (伏羲) 统一配置 (Pydantic v2)

21 Pydantic BaseModel classes:
  20 topic sub-models + 1 top-level Config

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
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, ConfigDict

# ── Auto-load .env from repo root ────────────────────────────────────
# This runs before any data_root() call, so FUXI_DATA_ROOT in .env
# is available to the pipeline without manual sourcing.
from dotenv import load_dotenv

_env_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '.env')
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
    matrix_file: str = ""
    barcodes_file: str = ""
    features_file: str = ""
    csv_sep: Optional[str] = None
    csv_decimal: str = '.'
    gene_symbol_column: str = ''
    input_h5ad: str = ""
    backed: str = ""
    h5_file_pattern: str = "*filtered_feature_bc_matrix.h5"
    h5_dir: str = ""
    fragment_file: str = ""
    # ── Preprocessed (embedded metadata columns) format ──
    file_pattern: str = "*.tsv.gz"
    separator: str = ""          # empty = auto-detect (tab vs comma)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 2 — SampleMetaConfig
# ═══════════════════════════════════════════════════════════════════════
class SampleMetaConfig(BaseModel):
    """Sample / stage metadata mapping."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    sample_map: Dict[int, str] = Field(default_factory=dict)
    stage_map: Dict[int, str] = Field(default_factory=dict)
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


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 7 — PCASettings
# ═══════════════════════════════════════════════════════════════════════
class PCASettings(BaseModel):
    """Principal component analysis settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_pcs_full: int = 100
    n_pcs_use: int = 50


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 8 — HarmonySettings
# ═══════════════════════════════════════════════════════════════════════
class HarmonySettings(BaseModel):
    """Harmony batch correction settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    use_harmony: bool = True
    batch_key: str = "sample"
    max_iter: int = 20

    diagnose: bool = True
    diagnose_report: bool = True
    diagnose_exclude_patterns: list[str] = Field(
        default_factory=lambda: ['*leiden*', '*cell_type*', '*annotation*', '*annotated*']
    )
    gini_batch_threshold: float = 0.3
    gini_biology_threshold: float = 0.6
    collinearity_guard: bool = True


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 9 — ClusteringSettings
# ═══════════════════════════════════════════════════════════════════════
class ClusteringSettings(BaseModel):
    """Clustering, UMAP, and parameter grid search settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_neighbors: int = 30
    leiden_resolutions: List[float] = Field(
        default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    )
    param_grid_n_neighbors: list = Field(default_factory=lambda: [15, 20, 30])
    param_grid_resolutions: list = Field(default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    leiden_flavor: str = "igraph"
    best_resolution: float = 1.0
    best_n_neighbors: int = 0
    cluster_selection_method: Optional[str] = "multi_metric"
    multi_metric_weights: dict = Field(default_factory=lambda: {
        "silhouette": 0.2,
        "stability": 0.2,
        "cluster_coherence": 0.3,
        "splitting_gain": 0.2,
        "kb_annotatable_rate": 0.1
    })
    multi_metric_n_stability_seeds: int = 5
    multi_metric_adaptive_resolution: bool = True
    multi_metric_coverage_ratio_threshold: float = 1.5
    multi_metric_coherence_dominance: float = 1.5
    multi_metric_granularity_cv_threshold: float = 0.05
    multi_metric_granularity_min_clusters: int = 10
    multi_metric_de_gate_threshold: int = 25
    umap_selection_method: Optional[str] = "convex_hull"
    param_grid_min_dist: Optional[list] = Field(default_factory=lambda: [0.1, 0.3, 0.5])
    param_grid_spread: Optional[list] = Field(default_factory=lambda: [1.0])
    umap_min_dist: float = 0.3
    umap_spread: float = 1.0


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 10 — MarkerSettings
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
    expert_rule_pval_cutoff: float = 0.0
    validation_n_top_genes: int = 15
    validation_min_overlap: float = 0.5
    validation_marginal_threshold: float = 0.25
    validation_pass_rate_min: float = 0.1
    quality_gate_min_pass_rate: float = 0.10
    developmental_mode: bool = False
    step10_groupby: List[str] = Field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 11 — DESettings
# ═══════════════════════════════════════════════════════════════════════
class DESettings(BaseModel):
    """Differential expression analysis settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    method: str = "wilcoxon"
    n_genes: int = 50
    pval_cutoff: float = 0.05
    logfc_cutoff: float = 0.25
    stage_pairwise: bool = True
    auto_switch_on_low_quality: bool = False


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 12 — TrajectorySettings
# ═══════════════════════════════════════════════════════════════════════
class TrajectorySettings(BaseModel):
    """Pseudotime / trajectory analysis settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    root_cell_types: List[str] = Field(default_factory=list)
    root_markers: List[str] = Field(default_factory=list)
    n_diffmap_comps: int = 15
    n_branchings: int = 2
    pseudotime_genes: List[str] = Field(default_factory=list)
    pseudotime_n_branch_de: int = 10
    pseudotime_n_correlated: int = 10
    pseudotime_cor_pval: float = 0.05


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 13 — EnrichmentSettings
# ═══════════════════════════════════════════════════════════════════════
class EnrichmentSettings(BaseModel):
    """Gene-set enrichment analysis settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    run: bool = True
    method: str = "both"
    gene_sets: list = Field(default_factory=lambda: [
        'GO_Biological_Process_2023',
        'KEGG_2021_Human',
    ])
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
# Sub-model 14 — GRNSettings
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
# Sub-model 15 — CCISettings
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
# Sub-model 16 — DownsampleSettings
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
# Sub-model 17 — SpatialConfig
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
# Sub-model 18 — ATACConfig
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


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 19 — ExecutionConfig
# ═══════════════════════════════════════════════════════════════════════
class ExecutionConfig(BaseModel):
    """Execution environment settings."""
    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    n_jobs: int = 0
    limit_blas_threads: bool = True
    random_seed: int = 42
    scanpy_verbosity: int = 2
    force_csr: bool = True
    use_float32: bool = True


# ═══════════════════════════════════════════════════════════════════════
# Sub-model 20 — AIConfig
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
    subcluster: bool = True
    deg_design: bool = False
    interpretation: bool = True
    cache_responses: bool = True
    unconstrained_annotation: bool = False


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
    h5ad_tempdir: str = "/tmp/Fuxi"
    cleanup_intermediates: bool = False
    perf_monitoring: bool = True

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
    target_class: str = ""
    target_order: str = ""

    # ═══════════════════════════════════════════════════════════════════
    # 20 个主题子模型
    # ═══════════════════════════════════════════════════════════════════
    data_input: DataInputConfig = Field(default_factory=DataInputConfig)
    sample_meta: SampleMetaConfig = Field(default_factory=SampleMetaConfig)
    qc: QCSettings = Field(default_factory=QCSettings)
    scrublet: ScrubletSettings = Field(default_factory=ScrubletSettings)
    normalization: NormalizationSettings = Field(default_factory=NormalizationSettings)
    hvg: HVGSettings = Field(default_factory=HVGSettings)
    pca: PCASettings = Field(default_factory=PCASettings)
    harmony: HarmonySettings = Field(default_factory=HarmonySettings)
    clustering: ClusteringSettings = Field(default_factory=ClusteringSettings)
    marker: MarkerSettings = Field(default_factory=MarkerSettings)
    de: DESettings = Field(default_factory=DESettings)
    trajectory: TrajectorySettings = Field(default_factory=TrajectorySettings)
    enrichment: EnrichmentSettings = Field(default_factory=EnrichmentSettings)
    grn: GRNSettings = Field(default_factory=GRNSettings)
    cci: CCISettings = Field(default_factory=CCISettings)
    downsample: DownsampleSettings = Field(default_factory=DownsampleSettings)
    spatial: SpatialConfig = Field(default_factory=SpatialConfig)
    atac: ATACConfig = Field(default_factory=ATACConfig)
    execution: ExecutionConfig = Field(default_factory=ExecutionConfig)
    ai: AIConfig = Field(default_factory=AIConfig)


    def model_post_init(self, __context):
        """Resolve relative paths after construction.

        Replaces the old resolve_paths() with Pydantic-native
        post-init hook.
        """
        base = self.project_dir if self.project_dir else os.path.dirname(os.path.abspath(__file__))

        # Treat '.' as "not set" for mtx_dir
        if self.data_input.mtx_dir == '.':
            self.data_input.mtx_dir = ""

        # Resolve top-level relative paths to absolute
        for attr in (
            "data_dir", "results_dir", "h5ad_dir",
            "figure_dir", "table_dir", "log_dir",
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
            _data_root = os.environ.get('FUXI_DATA_ROOT') or os.environ.get('SCRNA_DATA_ROOT')
            if _data_root:
                dataset_id = os.path.basename(self.project_dir or base)
                self.data_dir = os.path.join(_data_root, dataset_id)
            else:
                self.data_dir = base

        # Auto-fill mtx_dir and h5_dir from data_dir
        if not self.data_input.mtx_dir:
            self.data_input.mtx_dir = self.data_dir
        if not self.data_input.h5_dir:
            self.data_input.h5_dir = self.data_dir

        # Subset filter: auto-append suffix to output dirs
        ds = self.downsample
        if ds.sample_keep or (ds.obs_filter and ds.obs_filter.strip()):
            suffix = ds.subset_suffix if ds.subset_suffix else "_subset"
            self.h5ad_dir = self.h5ad_dir.rstrip('/\\') + suffix
            self.figure_dir = self.figure_dir.rstrip('/\\') + suffix
            self.table_dir = self.table_dir.rstrip('/\\') + suffix
            self.log_dir = self.log_dir.rstrip('/\\') + suffix
            print(f"[Config] Subset active → output dir suffix: '{suffix}'")

        # tissue_kb auto-inference from tissue
        if not self.tissue_kb and self.tissue not in ("unknown", ""):
            _repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            _kb_dir = os.path.join(_repo_root, "rna", "tissue_ontologies", self.tissue)
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
    def harmony_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_harmony.h5ad")

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
        return os.path.join(self.h5ad_dir, "01_filtered.h5ad")

    @property
    def processed_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "02_processed.h5ad")

    @property
    def clustered_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "03_clustered.h5ad")

    @property
    def trajectory_h5ad(self) -> str:
        return os.path.join(self.h5ad_dir, "07_trajectory.h5ad")

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
