#!/usr/bin/env python3
"""
config.py — Fuxi (伏羲) 统一配置
==================================

合并 scRNAseq_pipeline 和 ATACseq_pipeline 的配置系统。
一个 Config dataclass 包含所有组学类型的字段。

设计原则:
  - 所有参数集中在一个 dataclass 中
  - modality 字段区分组学类型: 'rna' | 'atac' | 'spatial'
  - 向后兼容现有项目配置文件
  - 路径自动解析: 默认所有路径相对于 config.py 所在目录

使用方法:
    from core.config import Config, AIConfig, CFG
    CFG.resolve_paths()
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class AIConfig:
    """AI configuration — all AI features controlled here"""
    # Global switch
    enabled: bool = False

    # Inference endpoint
    api_base: str = ""
    model: str = "deepseek-v4-flash"
    api_key: str = ""              # 留空则从环境变量 LLM_API_KEY 读取
    max_tokens: int = 4096
    temperature: float = 0.1
    thinking_enabled: bool = True
    reasoning_effort: str = "high"
    timeout: Optional[int] = None

    # RNA task-level switches
    ai_qc_review: bool = False
    ai_param_suggest: bool = False
    ai_annotation: bool = True
    ai_subcluster: bool = True
    ai_deg_design: bool = False
    ai_interpretation: bool = True
    ai_cache_responses: bool = True

    # Unconstrained annotation mode (v3.1.0+).
    # When True, AI is NOT constrained to KB candidates and can freely
    # suggest cell types.  Useful for audit mode (KB blind-spot detection)
    # or novel tissue types.  Default False (backward compatible).
    unconstrained_annotation: bool = False


@dataclass
class Config:
    """Fuxi unified config — 包含 RNA + ATAC 所有字段"""

    # ═══════════════════════════════════════════════════════════════════
    #  组学类型
    # ═══════════════════════════════════════════════════════════════════
    modality: str = "rna"                # 'rna' | 'atac' | 'spatial'

    # ═══════════════════════════════════════════════════════════════════
    #  路径设置（通用）
    # ═══════════════════════════════════════════════════════════════════
    data_dir: str = ""
    results_dir: str = "results"
    h5ad_dir: str = "results/h5ad"
    figure_dir: str = "results/figures"
    table_dir: str = "results/tables"
    log_dir: str = "logs"
    project_dir: str = ""

    # ═══════════════════════════════════════════════════════════════════
    #  数据输入格式（通用）
    # ═══════════════════════════════════════════════════════════════════
    # RNA: '10X_mtx' | 'csv_matrix' | 'h5ad' | '10X_h5'
    # ATAC: '10x_fragments' | 'h5ad'
    data_format: str = "10X_mtx"

    # ── RNA: 10X MTX 格式 ──
    mtx_prefix: str = ""
    mtx_dir: str = ""

    # ── RNA: CSV 矩阵格式 ──
    matrix_file: str = ""
    barcodes_file: str = ""
    features_file: str = ""
    csv_sep: Optional[str] = None
    csv_decimal: str = '.'
    gene_symbol_column: str = ''

    # ── 通用: h5ad 格式 ──
    input_h5ad: str = ""
    backed: str = ""                     # ''=全量, 'r'=只读backed

    # ── RNA: 10X HDF5 格式 ──
    h5_file_pattern: str = "*filtered_feature_bc_matrix.h5"
    h5_dir: str = ""

    # ── ATAC: 10X fragment 模式 ──
    fragment_file: str = ""
    barcodes_file: str = ""

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 样本元数据映射
    # ═══════════════════════════════════════════════════════════════════
    sample_map: Dict[int, str] = field(default_factory=dict)
    stage_map: Dict[int, str] = field(default_factory=dict)
    stage_order: List[str] = field(default_factory=list)
    meta_columns: Dict[str, str] = field(default_factory=dict)
    barcode_parse_regex: str = ""
    barcode_parse_groups: Dict[str, str] = field(default_factory=dict)

    # ═══════════════════════════════════════════════════════════════════
    #  数据集元信息（通用）
    # ═══════════════════════════════════════════════════════════════════
    tissue: str = "unknown"
    species: str = "human"

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: QC 阈值
    # ═══════════════════════════════════════════════════════════════════
    min_genes: int = 500
    max_genes: int = 7500
    max_pct_mito: float = 20.0
    mt_gene_pattern: str = "MT-"
    mt_gene_list: List[str] = field(default_factory=list)
    min_genes_per_umi: float = 0.7
    min_cells_per_gene: int = 3
    use_adaptive_thresholds: bool = False
    mad_n_mads: float = 3.0

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: 参考基因组
    # ═══════════════════════════════════════════════════════════════════
    genome: str = "hg38"
    chrom_sizes: str = ""
    blacklist_bed: str = ""
    tss_bed: str = ""

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: QC 阈值
    # ═══════════════════════════════════════════════════════════════════
    min_fragments: int = 1000
    max_fragments: int = 50000
    min_tsse: float = 7.0
    max_blacklist_ratio: float = 0.05
    min_peak_region_fragments: int = 300

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: Peak calling
    # ═══════════════════════════════════════════════════════════════════
    peak_qval: float = 0.05
    peak_width: int = 500
    use_macs3: bool = True

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: Scrublet 双细胞检测
    # ═══════════════════════════════════════════════════════════════════
    run_scrublet: bool = True
    scrublet_expected_doublet_rate: float = 0.06
    scrublet_batch_key: str = "sample"
    scrublet_min_counts: int = 2
    scrublet_min_cells: int = 3
    scrublet_min_gene_var_pctl: int = 85
    scrublet_n_prin_comps: int = 30

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 归一化与 HVG
    # ═══════════════════════════════════════════════════════════════════
    normalize_target_sum: float = 1e4
    n_top_genes: int = 4000
    hvg_flavor: str = "seurat_v3"
    hvg_batch_key: str = "sample"
    use_regress_out: bool = True

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: PCA
    # ═══════════════════════════════════════════════════════════════════
    n_pcs_full: int = 100
    n_pcs_use: int = 50

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: Harmony 批次校正
    # ═══════════════════════════════════════════════════════════════════
    use_harmony: bool = True
    harmony_batch_key: str = "sample"
    harmony_max_iter: int = 20

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 聚类与 UMAP
    # ═══════════════════════════════════════════════════════════════════
    n_neighbors: int = 30
    leiden_resolutions: List[float] = field(
        default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    )
    param_grid_n_neighbors: list = field(default_factory=lambda: [15, 20, 30])
    param_grid_resolutions: list = field(default_factory=lambda: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
    leiden_flavor: str = "igraph"
    best_resolution: float = 1.0

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: 降维
    # ═══════════════════════════════════════════════════════════════════
    n_features: int = 50000
    n_spectral: int = 30

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 细胞类型注释
    # ═══════════════════════════════════════════════════════════════════
    marker_dict: Dict[str, List[str]] = field(default_factory=dict)
    subcluster_types: List[str] = field(default_factory=list)
    subcluster_resolution: float = 0.4
    min_cells_subcluster: int = 50
    tissue_kb: str = ""
    tissue_ontology: str = ""

    # Phylogenetic filtering (v3.0.0+ KB feature).
    # When non-empty, marker_scoring applies taxonomic weighting.
    #   target_class  (str): Desired class (纲), e.g. "Mammalia".
    #   target_order  (str): Desired order (目), e.g. "Primates".
    # Both default to "" (no filtering — all KB sources used at full weight).
    target_class: str = ""
    target_order: str = ""

    # Expert-rule constraint parameters (v3.0.0+ self-audit).
    #
    # Layer 1 — Precise overrides (0 / 0.0 = use template default):
    #   expert_rule_top_n (int):
    #       Only examine the top-N DE genes for rule-matching.
    #       0 = use the value from the strictness template.
    #   expert_rule_pval_cutoff (float):
    #       Only consider genes with pvals_adj < this value.
    #       0.0 = use the value from the strictness template.
    #
    # Layer 2 — Convenience template:
    #   expert_rule_strictness (str):
    #       "strict"   → top_n=50,   pval=0.01  (mature tissue, high-confidence)
    #       "default"  → top_n=50,   pval=0.05  (general purpose, built-in default)
    #       "deep"     → top_n=200,  pval=0.05  (KB markers rank deep but significant)
    #       "wide"     → top_n=1000, pval=0.05  (developmental data, first triggers appear)
    #       "relaxed"  → top_n=5000, pval=0.05  (developmental/organoid sweet spot)
    #       "manual"   → requires both top_n (>0) AND pval_cutoff (>0.0)
    #
    #   pval_cutoff stays ≤ 0.05 in ALL presets — this is the last line of
    #   defence against noise-triggered rules.  To relax pval beyond 0.05,
    #   you MUST use "manual" + explicit expert_rule_pval_cutoff.
    expert_rule_strictness: str = "default"
    expert_rule_top_n: int = 0
    expert_rule_pval_cutoff: float = 0.0

    # Marker validation thresholds (v3.1.0+).
    # Controls how StandardOntology.validate() cross-checks assigned cell types
    # against KB markers using top DE genes.
    #   marker_validation_n_top_genes (int):
    #       Number of top DE genes per cluster to compare.  Default 15.
    #   marker_validation_min_overlap (float):
    #       Minimum overlap ratio (found/KB_total) for PASS status.
    #       Default 0.5 (backward compatible).
    #   marker_validation_marginal_threshold (float):
    #       Threshold for MARGINAL tier (PASS > MARGINAL > LOW > FAIL).
    #       Default 0.25.  Set to 0 to disable MARGINAL tier.
    #   marker_validation_pass_rate_min (float):
    #       Minimum PASS cell-rate for trajectory quality gate (Step 08).
    #       Also used by Steps 07/09 for quality warnings.  Default 0.1.
    marker_validation_n_top_genes: int = 15
    marker_validation_min_overlap: float = 0.5
    marker_validation_marginal_threshold: float = 0.25
    marker_validation_pass_rate_min: float = 0.1

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 差异表达分析
    # ═══════════════════════════════════════════════════════════════════
    de_method: str = "wilcoxon"
    de_n_genes: int = 50
    de_pval_cutoff: float = 0.05
    de_logfc_cutoff: float = 0.25
    de_stage_pairwise: bool = True
    # When True and marker_validation pass_rate < threshold, Step 07
    # automatically uses leiden-based grouping instead of cell_type labels.
    # Default False: only warns, user decides.  (v3.1.0+)
    de_auto_switch_on_low_quality: bool = False

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: 差异分析
    # ═══════════════════════════════════════════════════════════════════
    marker_peaks_log2fc: float = 0.5
    marker_peaks_fdr: float = 0.05

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: Motif
    # ═══════════════════════════════════════════════════════════════════
    motif_db: str = "JASPAR2024"

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 轨迹分析
    # ═══════════════════════════════════════════════════════════════════
    root_cell_types: List[str] = field(default_factory=list)
    root_markers: List[str] = field(default_factory=list)
    n_diffmap_comps: int = 15
    n_branchings: int = 2

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: 轨迹
    # ═══════════════════════════════════════════════════════════════════
    terminal_cell_types: List[str] = field(default_factory=list)

    # ═══════════════════════════════════════════════════════════════════
    #  富集分析（通用）
    # ═══════════════════════════════════════════════════════════════════
    run_enrichment: bool = True
    enrichment_method: str = "both"      # 'ora' | 'prerank' | 'both'
    enrichment_gene_sets: list = field(
        default_factory=lambda: [
            'GO_Biological_Process_2023',
            'KEGG_2021_Human',
        ]
    )
    enrichment_organism: str = "human"
    enrichment_n_top_genes: int = 200
    enrichment_pval_cutoff: float = 0.05
    enrichment_min_size: int = 10
    enrichment_max_size: int = 500
    enrichment_permutations: int = 1000
    peak_gene_distance: int = 100000    # ATAC: peak-to-gene 映射距离

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: 降采样 (Step 01)
    # ═══════════════════════════════════════════════════════════════════
    downsample_target: Optional[int] = None
    downsample_strategy: str = "stratified"
    downsample_max_per_sample: Optional[int] = None
    downsample_random_seed: int = 42

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: 降采样
    # ═══════════════════════════════════════════════════════════════════
    max_cells: Optional[int] = None

    # ═══════════════════════════════════════════════════════════════════
    #  执行环境（通用）
    # ═══════════════════════════════════════════════════════════════════
    n_jobs: int = 0                     # 0 = os.cpu_count()；项目 config 可覆写
    limit_blas_threads: bool = True     # 控制 OMP/MKL 线程数，锁定为 n_jobs 值
    random_seed: int = 42
    scanpy_verbosity: int = 2
    force_csr: bool = True
    use_float32: bool = False
    h5ad_compression: str = "gzip"       # 'gzip' | 'lzf' | 'zstd'
    h5ad_tempdir: str = "/tmp/Fuxi"      # safe_write 临时目录
    cleanup_intermediates: bool = False  # ATAC: 自动删除上游中间 checkpoint

    # ═══════════════════════════════════════════════════════════════════
    #  AI 配置（通用）
    # ═══════════════════════════════════════════════════════════════════
    ai: AIConfig = field(default_factory=AIConfig)

    # ═══════════════════════════════════════════════════════════════════
    #  ATAC: RNA 整合 (Step 09)
    # ═══════════════════════════════════════════════════════════════════
    rna_h5ad: str = ""

    # ═══════════════════════════════════════════════════════════════════
    #  RNA: checkpoint 路径（属性）
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
        return os.path.join(self.h5ad_dir, "02_normalized.h5ad")

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

    # ──────────────────────────────────────────────────────────────────
    #  方法
    # ──────────────────────────────────────────────────────────────────
    def resolve_paths(self):
        """解析所有路径。非绝对路径视为相对于 project_dir 或 config.py 所在目录。"""
        base = self.project_dir if self.project_dir else os.path.dirname(os.path.abspath(__file__))
        for attr in [
            "data_dir", "results_dir", "h5ad_dir",
            "figure_dir", "table_dir", "log_dir",
            "mtx_dir", "h5_dir",
        ]:
            val = getattr(self, attr)
            if val and not os.path.isabs(val):
                setattr(self, attr, os.path.join(base, val))

        if not self.data_dir:
            self.data_dir = base
        if not self.mtx_dir:
            self.mtx_dir = self.data_dir
        if not self.h5_dir:
            self.h5_dir = self.data_dir

        for d in [self.results_dir, self.h5ad_dir,
                  self.figure_dir, self.table_dir, self.log_dir]:
            os.makedirs(d, exist_ok=True)

    def has_sample_mapping(self) -> bool:
        return len(self.sample_map) > 0

    def has_stage_mapping(self) -> bool:
        return len(self.stage_map) > 0

    def has_markers(self) -> bool:
        return len(self.marker_dict) > 0

    def has_rna_data(self) -> bool:
        """ATAC: check if RNA data is available for integration"""
        return bool(self.rna_h5ad) and os.path.exists(self.rna_h5ad)


# ═══════════════════════════════════════════════════════════════════════
#  全局实例
# ═══════════════════════════════════════════════════════════════════════
CFG = Config()
