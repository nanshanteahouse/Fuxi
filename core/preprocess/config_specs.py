"""Per-format config specifications — the single source for starter configs.

Every format's starter config (fields, defaults, placeholders, comments)
is declared here as data.  The scaffold CLI (``core/config/scaffold.py``)
renders these specs on demand (``--format``) as human starting points;
``generate_config()`` assembles project configs directly from the specs
at runtime — no template files exist on disk.

Field names/defaults come from :mod:`core.config.schema` (single source of
truth) — ``validate_specs()`` fails on any dotted path that is not a real
schema field, so schema renames surface immediately instead of silently
dropping template entries.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, Union

from pydantic import BaseModel

from core.config.introspect import _resolve_base_model
from core.config.schema import Config

# ═══════════════════════════════════════════════════════════════════
# Data model
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class SpecComment:
    """Raw comment lines (and blank lines) rendered verbatim.

    Each line is rendered with ``"# "`` prefix unless it already starts
    with ``"#"``; ``""`` renders as a blank line.  ``indent`` is the YAML
    depth (2 spaces per level).
    """

    lines: Tuple[str, ...]
    indent: int = 0


@dataclass(frozen=True)
class SpecField:
    """One config field with a dotted schema path.

    ``placeholder`` (e.g. ``"TISSUE"``) renders as ``{{TISSUE}}`` and is
    filled from detected values at assembly time; ``value`` is the typed
    literal used when no placeholder is set.
    """

    path: str
    value: Any = None
    placeholder: str = ""
    comment: str = ""


@dataclass(frozen=True)
class FormatSpec:
    """Starter-config definition for one data format."""

    key: str
    modality: str
    data_format: str
    items: Tuple[Union[SpecComment, SpecField], ...]


# ═══════════════════════════════════════════════════════════════════
# Placeholders (filled by generate_config from detected values)
# ═══════════════════════════════════════════════════════════════════

SUPPORTED_PLACEHOLDERS = frozenset(
    {
        "MTX_PREFIX",
        "MTX_DIR",
        "MTX_DIR_PATTERN",
        "MATRIX_FILE",
        "BARCODES_FILE",
        "FEATURES_FILE",
        "FRAGMENT_FILE",
        "LIBRARY_ID",
        "TISSUE",
        "SPECIES",
        "GENOME",
        "EXPRESSION_TYPE",
    }
)


# ═══════════════════════════════════════════════════════════════════
# Shared comment blocks (RNA templates)
# ═══════════════════════════════════════════════════════════════════

_GLOBAL_NOTE = SpecComment(
    (
        "公共默认参数已移至 global.yaml (参考 global.example.yaml)",
        "本模板仅包含数据集特异性和生物学相关参数",
    )
)


def _header(data_format: str) -> Tuple[SpecComment, SpecField, SpecComment]:
    return (
        SpecComment(("── Data format ──",)),
        SpecField("data_format", value=data_format),
        _GLOBAL_NOTE,
    )


_SCVI_BLOCK = (
    SpecField("integration.scvi.n_latent", value=30),
    SpecField("integration.scvi.n_layers", value=2),
    SpecField("integration.scvi.n_hidden", value=128),
    SpecField("integration.scvi.max_epochs", value=400),
    SpecField("integration.scvi.batch_key", value="sample"),
    SpecField("integration.scvi.use_gpu", value=True),
    SpecField("integration.scvi.train_size", value=0.9),
    SpecField(
        "integration.scvi.batch_size",
        value=1024,
        comment="128→1024: scverse 官方线性加速建议（census 35M 用 1024），实测 0.51→0.32 s/epoch",
    ),
    SpecField(
        "integration.scvi.early_stopping",
        value=True,
        comment="ELBO 早停（Lightning patience=3 默认），避免傻跑满 max_epochs",
    ),
    SpecField(
        "integration.scvi.precision",
        value="16-mixed",
        comment='3090 AMP 加速（实测 ~1.5-2x）；CPU 自动降级 "32"',
    ),
    SpecComment(("── WSL2 关键（实测 GPU 饿死根因修复，勿改）──",), indent=4),
    SpecField(
        "integration.scvi.datasplitter_kwargs.num_workers",
        value=0,
        comment="WSL2 IPC 阻塞 → GPU 0% 利用率；0 = 主进程加载",
    ),
    SpecField(
        "integration.scvi.trainer_kwargs.enable_progress_bar",
        value=False,
        comment="去 tqdm 每 batch 管道写（加剧 pipe_write 阻塞）",
    ),
    SpecComment(("Advanced params: see global.example.yaml → integration.scvi",), indent=4),
)

_DIAGNOSE_BLOCK = (
    SpecComment(
        ("-- v4.x+: batch diagnosis (auto-classify batch vs biology columns) --",), indent=2
    ),
    SpecComment(
        ("diagnose: true                        # 全量诊断（大细胞数深拷贝，百万级建议 false）",),
        indent=2,
    ),
    SpecComment(
        ("diagnose_max_cells: 50000            # 诊断子采样上限（默认 5 万，秒级）",), indent=2
    ),
    SpecComment(
        (
            "stream_raw: true                     # 流式写 .raw（从 02_qc 分块变换直写；>100 万细胞建议 true，峰值内存 -24%）",
        ),
        indent=2,
    ),
)

_QC_COMMON = (
    SpecField("qc.min_genes", value=200),
    SpecField("qc.max_genes", value=7500),
    SpecComment(
        ("max_pct_mito: 20.0  # global default; snRNA-seq→5-10%, heart/liver→15%",), indent=1
    ),
    SpecField("qc.min_genes_per_umi", value=0.7),
    SpecComment(("── Adaptive thresholds (managed by global.yaml) ──",), indent=1),
    SpecComment(("use_adaptive_thresholds: true   # ← global.yaml default",), indent=1),
    SpecComment(
        ("mad_n_mads: 3.0                  # MAD multiplier (developing tissue auto-uses 5.0)",),
        indent=1,
    ),
    SpecComment(
        ("ncount_max_mad: 5.0              # nCount upper bound MAD multiplier",), indent=1
    ),
    SpecComment(("── Single-nucleus (snRNA-seq) ──",), indent=1),
    SpecComment(("is_nuclei: true",), indent=1),
    SpecComment(("max_pct_mito_nuclei: 5.0",), indent=1),
)

_SCRUBLET = (
    SpecComment(
        (
            "serial_threshold: 15000             # cells: >threshold → serial, ≤threshold → parallel",
        ),
        indent=1,
    ),
    SpecField(
        "scrublet.run", value=True, comment="auto-disabled when expression_type != raw_counts"
    ),
)

_HVG = (
    SpecField("hvg.n_top_genes", value=4000),
    SpecField("hvg.batch_key", value="sample"),
    SpecField("hvg.flavor", value="seurat_v3"),
    SpecComment(("use_regress_out: false",), indent=1),
)

_CLUSTERING_COMMON = (
    SpecField("clustering.param_grid_resolutions", value=[0.3, 0.5, 0.8, 1.0, 1.5, 2.0]),
    SpecField("clustering.best_resolution", value=1.0),
    SpecField("clustering.best_n_neighbors", value=0),
    SpecField("clustering.umap_color_by_batch", value=False),
    SpecComment(
        ("batch_key_override: null  # optionally override integration.batch_key",), indent=1
    ),
    SpecComment(("── Multi-metric evaluation ──",), indent=1),
    SpecComment(
        ("stability_n_seeds: 12                  # random seeds for stability scoring",), indent=1
    ),
    SpecComment(
        ("stability_leiden_n_iterations: 1000    # Leiden iterations for stability",), indent=1
    ),
    SpecComment(
        ("umap_selection_metric: trustworthiness # UMAP quality metric (trustworthiness /",),
        indent=1,
    ),
    SpecComment(
        (
            "umap_maxiter: null                    # max UMAP iterations (null=auto). 200-500 for 1M+ cells.",
        ),
        indent=1,
    ),
    SpecComment(
        ("umap_n_epochs: null                   # UMAP training epochs (null=auto)",), indent=1
    ),
    SpecComment(("  continuity / distance_mse / sil_score)",), indent=1),
    SpecComment(
        ("param_grid_n_neighbors_adaptive: true  # auto-scale neighbors from n_obs",), indent=1
    ),
    SpecComment(
        ("plot_per_combo: false                 # per-(n_neighbors, resolution) UMAP plots",),
        indent=1,
    ),
    SpecComment(("── Funnel mode (large datasets) ──",), indent=1),
    SpecComment(
        ("funnel_enabled: true                   # progressive grid search for >100k cells",),
        indent=1,
    ),
    SpecComment(
        ("funnel_threshold: 100000               # min n_obs to trigger funnel",), indent=1
    ),
    SpecComment(("funnel_subsample_size: 50000           # subsample for funnel grid",), indent=1),
    SpecComment(
        ("funnel_top_k: 3                        # top candidates re-validated on full data",),
        indent=1,
    ),
    SpecComment(("── Target cluster mode ──",), indent=1),
    SpecComment(
        ("target_n_clusters: null                # target cluster count (binary-search",), indent=1
    ),
    SpecComment(("  resolution instead of grid)",), indent=1),
    SpecComment(
        ("target_search_max_iters: 10            # max binary-search iterations",), indent=1
    ),
)

_MARKER = (
    SpecField("marker.marker_dict", value={}),
    SpecField("marker.subcluster_types", value=[]),
    SpecComment(
        ("candidate_pool_expand_steps: [50, 100, 200]   # adaptive scoring tiers → n_genes",),
        indent=1,
    ),
    SpecField("marker.subcluster_resolution", value=0.4),
    SpecField("marker.min_cells_subcluster", value=50),
)

_TRAJECTORY = (
    SpecField("trajectory.root_cell_types", value=[]),
    SpecComment(("root_markers: [SOX2, PAX6, NES]",), indent=1),
    SpecComment(("n_diffmap_comps: 15",), indent=1),
    SpecComment(("n_branchings: 2",), indent=1),
)

_DE = (
    SpecField("de.stage_pairwise", value=True),
    SpecComment(("auto_switch_on_low_quality: true",), indent=1),
)

_PSEUDOBULK_BLOCK = SpecComment(
    (
        "── Pseudobulk DE (PyDESeq2) ──",
        "de.pseudobulk:",
        "  sample_col: sample",
        '  design: "~condition"',
        "  contrast_column: condition",
        "  contrast_treatment: treated",
        "  contrast_baseline: control",
        "  alpha: 0.05",
    )
)

_ENRICHMENT_BLOCK = SpecComment(
    (
        "── Enrichment (uncomment to enable) ──",
        "enrichment:",
        "  run: true",
        "  method: both          # ora | prerank | both",
        "  gene_sets: [GO_Biological_Process_2023, KEGG_2021_Human]",
        "  organism: human",
    )
)

_GRN = (
    SpecField("grn.run", value=True),
    SpecField("grn.method", value="decoupler"),
    SpecField("grn.species", value="human"),
    SpecField("grn.n_top_regulons", value=50),
    SpecField("grn.min_regulon_size", value=5),
)

_CCI = (
    SpecField("cci.run", value=True),
    SpecComment(("method: liana",), indent=1),
    SpecComment(("lr_database: consensus",), indent=1),
    SpecComment(("permutations: 1000",), indent=1),
    SpecComment(("n_top_interactions: 50",), indent=1),
    SpecComment(("spatial_method: liana_spatial",), indent=1),
    SpecComment(("spatial_distance: 0.0",), indent=1),
    SpecComment(('lr_cache_dir: ""',), indent=1),
)

_AI_BLOCK = SpecComment(
    (
        "── AI settings (uncomment to enable) ──",
        "ai:",
        "  enabled: true",
        "  api_base: https://api.deepseek.com/v1",
        "  model: deepseek-v4-pro",
        "  max_tokens: 32768",
        "  temperature: 0.1",
        "  thinking_enabled: true",
        "  reasoning_effort: high",
        "  annotation: true",
        "  subcluster: true",
        "  interpretation: true",
        "  cache_responses: true",
    )
)

_DOWNSAMPLE_BLOCK = SpecComment(
    (
        "── Downsampling (optional) ──",
        "downsample:",
        "  target: 5000",
        "  strategy: random",
    )
)

_TRAILER = (
    SpecComment(("",), indent=0),
    SpecComment(
        (
            "h5ad_compression: zstd        # gzip | zstd | lzf。03 浮点高熵数据：gzip 文件小、zstd 写盘快（-32% 写时，文件 +40%）",
        ),
        indent=0,
    ),
    SpecComment(
        ("h5ad_compression_opts: 1      # gzip 压缩级别（默认 4→1：写盘 -44%，文件略大）",),
        indent=0,
    ),
)


# ═══════════════════════════════════════════════════════════════════
# RNA specs
# ═══════════════════════════════════════════════════════════════════

_10X_H5 = FormatSpec(
    key="10X_h5",
    modality="rna",
    data_format="10X_h5",
    items=(
        *_header("10X_h5"),
        SpecComment(("── Data input ──",)),
        SpecComment((".h5 file matching pattern (glob syntax)",), indent=1),
        SpecComment(
            ('Default matches Cell Ranger filtered output: "*filtered_feature_bc_matrix.h5"',),
            indent=1,
        ),
        SpecComment(('For raw matrix: "*raw_feature_bc_matrix.h5"',), indent=1),
        SpecComment(('Custom: "*_counts.h5" etc.',), indent=1),
        SpecField("data_input.h5_file_pattern", value="*filtered_feature_bc_matrix.h5"),
        SpecComment(('h5_dir: "data"',), indent=1),
        SpecComment(("",)),
        SpecComment(("── Sample metadata ──",)),
        SpecComment(
            ("When multiple .h5 files exist, sample names are auto-extracted from filenames.",),
            indent=1,
        ),
        SpecComment(
            ("For 10X multi-channel aggregation, barcodes use -1, -2, ... suffixes.",), indent=1
        ),
        SpecComment(("sample_map:",), indent=1),
        SpecComment(("  1: sample1",), indent=1),
        SpecComment(("  2: sample2",), indent=1),
        SpecField("sample_meta.sample_map", value={}),
        SpecComment(("",)),
        SpecComment(("barcode suffix -> developmental stage / experimental group",), indent=1),
        SpecComment(("stage_map:",), indent=1),
        SpecComment(("  1: Day3",), indent=1),
        SpecComment(("  2: Day7",), indent=1),
        SpecField("sample_meta.stage_map", value={}),
        SpecField("sample_meta.stage_order", value=[]),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", value="unknown"),
        SpecComment(
            ("tissue_kb: \"\"          # auto-inferred from 'tissue' when matching KB exists",),
            indent=0,
        ),
        SpecComment(
            ('tissue_maturity: "unknown"   # "developing" | "adult" | "unknown"',), indent=0
        ),
        SpecField("species", value="human"),
        SpecField("expression_type", value="raw_counts"),
        SpecComment(("raw_counts | log1p_counts | TPM | CPM | FPKM",), indent=0),
        SpecComment(
            ("TPM/FPKM/CPM: total_counts & complexity filters are auto-skipped.",), indent=0
        ),
        SpecComment(("Scrublet is auto-disabled for non-raw_counts data.",), indent=0),
        SpecComment(("",)),
        SpecComment(("── QC ──",)),
        *_QC_COMMON,
        SpecComment(("",)),
        SpecComment(("── Scrublet (doublet detection) ──",)),
        *_SCRUBLET,
        SpecComment(("",)),
        SpecComment(("── HVG ──",)),
        *_HVG,
        SpecComment(("",)),
        SpecComment(("── Batch correction / Integration ──",)),
        SpecField("integration.method", value="harmony"),
        SpecField("integration.batch_key", value="sample"),
        SpecComment(("",)),
        SpecComment(("── Neighbor persistence (skip re-computation) ──",), indent=1),
        SpecComment(
            ("persist_neighbors: false               # reuse neighbors from clustering step",),
            indent=1,
        ),
        SpecComment(
            ("persisted_neighbors_k: 50              # k for persisted neighbor graph",), indent=1
        ),
        *_DIAGNOSE_BLOCK,
        SpecComment(("-- scVI (deep generative model, method: scvi) --",), indent=1),
        *_SCVI_BLOCK,
        SpecComment(("",)),
        SpecComment(("── PCA ──",)),
        SpecComment(
            (
                "n_pcs_full = PCA 计算量（elbow 图横轴）；n_pcs_use = 下游消费维数（04+ 用 min() 兑底）",
            ),
            indent=0,
        ),
        SpecComment(
            ("method: scvi 时 n_pcs_use 建议 = 30（scVI latent 维，与 scvi.n_latent 一致）",),
            indent=0,
        ),
        SpecField("pca.n_pcs_full", value=100),
        SpecField("pca.n_pcs_use", value=50),
        SpecComment(("",)),
        SpecComment(("── Clustering ──",)),
        *_CLUSTERING_COMMON,
        SpecComment(("",)),
        SpecComment(("── Cell type markers ──",)),
        SpecComment(("Format: CellType: [marker1, marker2, ...]",), indent=0),
        SpecComment(("Replace with actual markers for your tissue type.",), indent=0),
        *_MARKER,
        SpecComment(("",)),
        SpecComment(("── Knowledge base ──",)),
        SpecComment(('tissue_kb: ""   # Tissue KB name from core/kb/, e.g. "retina"',), indent=0),
        SpecComment(("",)),
        SpecComment(("── Trajectory ──",)),
        *_TRAJECTORY,
        SpecComment(
            (
                "save_final_h5ad: true   # 输出 05_final.h5ad（默认开；该文件无下游消费者，可关以省空间）",
            ),
            indent=1,
        ),
        SpecComment(("",)),
        SpecComment(("── Differential expression ──",)),
        *_DE,
        _PSEUDOBULK_BLOCK,
        _ENRICHMENT_BLOCK,
        SpecComment(
            (
                "",
                "Tissue-specific gene set libraries (v4.0+):",
                "Available: CellMarker_Augmented_2021, PanglaoDB_Augmented_2021,",
                "           Tabula_Sapiens, Azimuth_Cell_Types, Allen_Brain_Atlas",
                "  gene_sets_tissue:",
                "    - CellMarker_Augmented_2021",
                "    - PanglaoDB_Augmented_2021",
                "",
                "Tissue-aware ranking/filtering (v4.0+):",
                "  'off'  = pure statistical ranking (default, current behavior)",
                "  'soft' = annotate tissue_relevant column + output _tissue_relevant.csv",
                "  'hard' = keep only tissue-relevant pathways",
                "  tissue_mode: soft",
                "  use_kb_relevance: true",
                "  redundancy_cluster: true",
            )
        ),
        SpecComment(("",)),
        SpecComment(("── GRN regulatory network analysis ──",)),
        *_GRN,
        SpecComment(("confidence_levels: [A, B, C]",), indent=1),
        SpecComment(("",)),
        SpecComment(("── CCI cell-cell interaction analysis ──",)),
        *_CCI,
        _AI_BLOCK,
        *_TRAILER,
        SpecComment(("",)),
        SpecComment(("── h5ad 增量写入 ──",)),
        SpecComment(
            (
                "incremental_io: true   # 增量写开关（默认开）：obs/obsm/obsp/uns 走 in-place 追加写",
            ),
            indent=0,
        ),
        SpecComment(
            (
                "                       # 关闭(false) → 全部回退全量 safe_write（WSL /mnt 不稳时的逃生口）",
            ),
            indent=0,
        ),
        SpecComment(("",)),
        _DOWNSAMPLE_BLOCK,
        SpecComment(("",)),
        SpecComment(
            ("── Subset filtering (run on subset of samples/cells, suffix auto-appended) ──",)
        ),
        SpecComment(("sample_keep: [GSM9292434_SCR205, GSM9292436_SCR206]",), indent=0),
        SpecComment(("obs_filter: \"stage == 'PCW8'\"",), indent=0),
        SpecComment(('subset_suffix: "_pcw8"',), indent=0),
    ),
)

_10X_MTX = FormatSpec(
    key="10X_mtx",
    modality="rna",
    data_format="10X_mtx",
    items=(
        *_header("10X_mtx"),
        SpecComment(("── Data input ──",)),
        SpecField("data_input.mtx_prefix", placeholder="MTX_PREFIX"),
        SpecField("data_input.mtx_dir", placeholder="MTX_DIR"),
        SpecComment(
            ('多样本: glob 匹配 mtx_dir 下的样本子目录 (如 "*/" 或 "GSM*/")。',), indent=1
        ),
        SpecComment(("Step 00 自动合并匹配到的目录; 空 = 传统单目录加载。",), indent=1),
        SpecField("data_input.mtx_dir_pattern", placeholder="MTX_DIR_PATTERN"),
        SpecComment(("可选: 正则从子目录名提取样本名 (空 = 用目录 basename)",), indent=1),
        SpecComment(('mtx_sample_regex: ""',), indent=1),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", placeholder="TISSUE"),
        SpecComment(
            ('tissue_kb: ""          # auto-inferred from "tissue" when matching KB exists',),
            indent=0,
        ),
        SpecComment(
            ('tissue_maturity: "unknown"   # "developing" | "adult" | "unknown"',), indent=0
        ),
        SpecField("species", placeholder="SPECIES"),
        SpecField("expression_type", placeholder="EXPRESSION_TYPE"),
        SpecComment(("raw_counts | log1p_counts | TPM | CPM | FPKM",), indent=0),
        SpecComment(
            ("TPM/FPKM/CPM: total_counts & complexity filters are auto-skipped.",), indent=0
        ),
        SpecComment(("Scrublet is auto-disabled for non-raw_counts data.",), indent=0),
        SpecComment(("Adjust min_genes / max_pct_mito accordingly.",), indent=0),
        SpecComment(("",)),
        SpecComment(("── Sample metadata ──",)),
        SpecComment(("Map 10X barcode suffixes (-1, -2, ...) to sample names.",), indent=1),
        SpecComment(("sample_map:",), indent=1),
        SpecComment(("  1: sample1",), indent=1),
        SpecComment(("  2: sample2",), indent=1),
        SpecField("sample_meta.sample_map", value={}),
        SpecComment(("",)),
        SpecComment(("Stage mapping (if developmental data)",), indent=1),
        SpecComment(("stage_map:",), indent=1),
        SpecComment(("  1: Day3",), indent=1),
        SpecComment(("  2: Day7",), indent=1),
        SpecField("sample_meta.stage_map", value={}),
        SpecField("sample_meta.stage_order", value=[]),
        SpecComment(("",)),
        SpecComment(("── QC ──",)),
        *_QC_COMMON,
        SpecComment(("",)),
        SpecComment(("── Scrublet (doublet detection) ──",)),
        *_SCRUBLET,
        SpecComment(("",)),
        SpecComment(("── HVG ──",)),
        *_HVG,
        SpecComment(("",)),
        SpecComment(("── Batch correction / Integration ──",)),
        SpecField("integration.method", value="harmony"),
        SpecField("integration.batch_key", value="sample"),
        SpecComment(("",)),
        *_DIAGNOSE_BLOCK,
        SpecComment(("-- scVI (deep generative model, method: scvi) --",), indent=1),
        *_SCVI_BLOCK,
        SpecComment(("",)),
        SpecComment(("── PCA ──",)),
        SpecComment(
            ("n_pcs_full = PCA 计算量；n_pcs_use = 下游消费维数（04+ 用 min() 兑底）",), indent=0
        ),
        SpecComment(
            ("method: scvi 时 n_pcs_use 建议 = 30（scVI latent 维，与 scvi.n_latent 一致）",),
            indent=0,
        ),
        SpecField("pca.n_pcs_full", value=100),
        SpecField("pca.n_pcs_use", value=50),
        SpecComment(("",)),
        SpecComment(("── Clustering ──",)),
        *_CLUSTERING_COMMON,
        SpecComment(("",)),
        SpecComment(("── Cell type markers ──",)),
        SpecComment(("TODO: Add known marker genes for {{TISSUE}} tissue.",), indent=0),
        SpecComment(
            (
                "umap_maxiter: null                    # max UMAP iterations (null=auto). 200-500 for 1M+ cells.",
            ),
            indent=1,
        ),
        SpecComment(
            ("umap_n_epochs: null                   # UMAP training epochs (null=auto)",), indent=1
        ),
        *_MARKER,
        SpecComment(("",)),
        SpecComment(("── Knowledge base ──",)),
        SpecComment(("TODO: Set to a tissue KB name if one exists in core/kb/.",), indent=0),
        SpecComment(('tissue_kb: ""   # e.g. "retina", "hypothalamus"',), indent=0),
        SpecComment(("",)),
        SpecComment(("── Trajectory ──",)),
        *_TRAJECTORY,
        SpecComment(("",)),
        SpecComment(("── Differential expression ──",)),
        *_DE,
        _PSEUDOBULK_BLOCK,
        _ENRICHMENT_BLOCK,
        SpecComment(("",)),
        SpecComment(("── GRN regulatory network analysis ──",)),
        *_GRN,
        SpecComment(("confidence_levels: [A, B, C]",), indent=1),
        SpecComment(("tissue_mode: off",), indent=1),
        SpecComment(("use_kb_relevance: false",), indent=1),
        SpecComment(("export_filtered: false",), indent=1),
        SpecComment(("",)),
        SpecComment(("── CCI cell-cell interaction analysis ──",)),
        *_CCI,
        _AI_BLOCK,
        *_TRAILER,
        SpecComment(("",)),
        _DOWNSAMPLE_BLOCK,
    ),
)

_CSV = FormatSpec(
    key="csv_matrix",
    modality="rna",
    data_format="csv_matrix",
    items=(
        *_header("csv_matrix"),
        SpecComment(("── Data input ──",)),
        SpecField("data_input.matrix_file", placeholder="MATRIX_FILE"),
        SpecField("data_input.barcodes_file", placeholder="BARCODES_FILE"),
        SpecField("data_input.features_file", placeholder="FEATURES_FILE"),
        SpecComment(('csv_sep: "\\t"',), indent=1),
        SpecComment(('csv_decimal: "."',), indent=1),
        SpecComment(("",)),
        SpecComment(("Metadata column mapping",), indent=1),
        SpecComment(("If barcodes/file CSV contains sample/stage/tissue info,",), indent=1),
        SpecComment(("map those columns to pipeline obs fields.",), indent=1),
        SpecComment(("meta_columns:",), indent=1),
        SpecComment(("  sample: sample",), indent=1),
        SpecComment(("  stage: age",), indent=1),
        SpecComment(("  tissue: tissue",), indent=1),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", placeholder="TISSUE"),
        SpecComment(
            ('tissue_maturity: "unknown"   # "developing" | "adult" | "unknown"',), indent=0
        ),
        SpecField("species", placeholder="SPECIES"),
        SpecField("expression_type", placeholder="EXPRESSION_TYPE"),
        SpecComment(("raw_counts | log1p_counts | TPM | CPM | FPKM",), indent=0),
        SpecComment(
            ("TPM/FPKM/CPM: total_counts & complexity filters are auto-skipped.",), indent=0
        ),
        SpecComment(("Scrublet is auto-disabled for non-raw_counts data.",), indent=0),
        SpecComment(("Adjust min_genes / max_pct_mito accordingly.",), indent=0),
        SpecComment(("",)),
        SpecComment(("── Sample metadata ──",)),
        SpecField("sample_meta.sample_map", value={}),
        SpecComment(("",)),
        SpecComment(("Stage mapping (if developmental data)",), indent=1),
        SpecComment(("stage_map:",), indent=1),
        SpecComment(("  Week5: EarlyFetal",), indent=1),
        SpecComment(("  Week12: MidFetal",), indent=1),
        SpecComment(("  Week24: LateFetal",), indent=1),
        SpecField("sample_meta.stage_map", value={}),
        SpecField("sample_meta.stage_order", value=[]),
        SpecComment(("",)),
        SpecComment(("── QC ──",)),
        *_QC_COMMON,
        SpecComment(("",)),
        SpecComment(("── Scrublet (doublet detection) ──",)),
        *_SCRUBLET,
        SpecComment(("",)),
        SpecComment(("── HVG ──",)),
        *_HVG,
        SpecComment(("",)),
        SpecComment(("── Batch correction / Integration ──",)),
        SpecField("integration.method", value="harmony"),
        SpecField("integration.batch_key", value="sample"),
        SpecComment(("",)),
        *_DIAGNOSE_BLOCK,
        SpecComment(("-- scVI (deep generative model, method: scvi) --",), indent=1),
        *_SCVI_BLOCK,
        SpecComment(("",)),
        SpecComment(("── PCA ──",)),
        SpecComment(
            ("n_pcs_full = PCA 计算量；n_pcs_use = 下游消费维数（04+ 用 min() 兑底）",), indent=0
        ),
        SpecComment(
            ("method: scvi 时 n_pcs_use 建议 = 30（scVI latent 维，与 scvi.n_latent 一致）",),
            indent=0,
        ),
        SpecField("pca.n_pcs_full", value=100),
        SpecField("pca.n_pcs_use", value=50),
        SpecComment(("",)),
        SpecComment(("── Clustering ──",)),
        *_CLUSTERING_COMMON,
        SpecComment(("",)),
        SpecComment(("── Cell type markers ──",)),
        SpecComment(("TODO: Add known marker genes for {{TISSUE}} tissue.",), indent=0),
        SpecComment(
            (
                "umap_maxiter: null                    # max UMAP iterations (null=auto). 200-500 for 1M+ cells.",
            ),
            indent=1,
        ),
        SpecComment(
            ("umap_n_epochs: null                   # UMAP training epochs (null=auto)",), indent=1
        ),
        *_MARKER,
        SpecComment(("",)),
        SpecComment(("── Knowledge base ──",)),
        SpecComment(("TODO: Set to a tissue KB name if one exists in core/kb/.",), indent=0),
        SpecComment(('tissue_kb: ""   # e.g. "retina", "hypothalamus"',), indent=0),
        SpecComment(("",)),
        SpecComment(("── Trajectory ──",)),
        *_TRAJECTORY,
        SpecComment(("",)),
        SpecComment(("── Differential expression ──",)),
        *_DE,
        _PSEUDOBULK_BLOCK,
        _ENRICHMENT_BLOCK,
        SpecComment(("",)),
        SpecComment(("── GRN regulatory network analysis ──",)),
        *_GRN,
        SpecComment(("confidence_levels: [A, B, C]",), indent=1),
        SpecComment(("",)),
        SpecComment(("── CCI cell-cell interaction analysis ──",)),
        *_CCI,
        _AI_BLOCK,
        *_TRAILER,
        SpecComment(("",)),
        _DOWNSAMPLE_BLOCK,
    ),
)

_PREPROCESSED = FormatSpec(
    key="preprocessed",
    modality="rna",
    data_format="preprocessed",
    items=(
        *_header("preprocessed"),
        SpecField(
            "expression_type",
            placeholder="EXPRESSION_TYPE",
            comment="log1p_counts / raw_counts / TPM / CPM / FPKM",
        ),
        SpecComment(("",)),
        SpecComment(("── Data input (optional overrides) ──",)),
        SpecComment(("All auto-detectable by default; only override if needed.",), indent=0),
        SpecComment(("data_input:",), indent=0),
        SpecComment(
            ('  file_pattern: "*.tsv.gz"   # glob pattern (default: *.tsv.gz)',), indent=0
        ),
        SpecComment(
            ('  separator: "\\t"            # delimiter (default: auto-detect tab vs comma)',),
            indent=0,
        ),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", placeholder="TISSUE"),
        SpecField("species", placeholder="SPECIES"),
        SpecComment(("",)),
        SpecComment(("── QC (relaxed — data already QC'd by original authors) ──",)),
        SpecComment(
            ("Complexity filter (min_genes_per_umi) and total_counts upper bound are",), indent=0
        ),
        SpecComment(
            ("auto-skipped for non-raw_counts data. MAD adaptive thresholds are NOT",), indent=0
        ),
        SpecComment(
            ("recommended here — the distribution of log-normalized values differs",), indent=0
        ),
        SpecComment(("from raw counts and would produce misleading MAD bounds.",), indent=0),
        SpecField("qc.min_genes", value=200),
        SpecField("qc.max_genes", value=7500),
        SpecComment(
            ("max_pct_mito: 20.0  # global default; snRNA-seq→5-10%, heart/liver→15%",), indent=1
        ),
        SpecComment(("",)),
        SpecComment(("── Scrublet (doublet detection) ──",)),
        *_SCRUBLET,
        SpecComment(("",)),
        SpecComment(("── HVG ──",)),
        *_HVG,
        SpecComment(("",)),
        SpecComment(("── Batch correction / Integration ──",)),
        SpecField("integration.method", value="harmony"),
        SpecField("integration.batch_key", value="sample"),
        SpecComment(("",)),
        *_DIAGNOSE_BLOCK,
        SpecComment(("-- scVI (deep generative model, method: scvi) --",), indent=1),
        *_SCVI_BLOCK,
        SpecComment(("",)),
        SpecComment(("── PCA ──",)),
        SpecComment(
            ("n_pcs_full = PCA 计算量；n_pcs_use = 下游消费维数（04+ 用 min() 兑底）",), indent=0
        ),
        SpecComment(
            ("method: scvi 时 n_pcs_use 建议 = 30（scVI latent 维，与 scvi.n_latent 一致）",),
            indent=0,
        ),
        SpecField("pca.n_pcs_full", value=100),
        SpecField("pca.n_pcs_use", value=50),
        SpecComment(("",)),
        SpecComment(("── Clustering ──",)),
        *_CLUSTERING_COMMON,
        SpecComment(("",)),
        SpecComment(("── Cell type markers ──",)),
        SpecComment(
            (
                "umap_maxiter: null                    # max UMAP iterations (null=auto). 200-500 for 1M+ cells.",
            ),
            indent=1,
        ),
        SpecComment(
            ("umap_n_epochs: null                   # UMAP training epochs (null=auto)",), indent=1
        ),
        *_MARKER,
        SpecComment(("",)),
        SpecComment(("── Trajectory ──",)),
        *_TRAJECTORY,
        SpecComment(("",)),
        SpecComment(("── Differential expression ──",)),
        *_DE,
        _PSEUDOBULK_BLOCK,
        SpecComment(("",)),
        SpecComment(("── GRN regulatory network ──",)),
        *_GRN,
        SpecComment(("",)),
        SpecComment(("── CCI cell-cell interaction ──",)),
        SpecField("cci.run", value=True),
        SpecComment(("",)),
        _DOWNSAMPLE_BLOCK,
    ),
)


# ═══════════════════════════════════════════════════════════════════
# ATAC spec
# ═══════════════════════════════════════════════════════════════════

_FRAGMENTS = FormatSpec(
    key="10x_fragments",
    modality="atac",
    data_format="10x_fragments",
    items=(
        SpecComment(("── Modality & data format ──",)),
        SpecField("modality", value="atac"),
        SpecField("data_format", value="10x_fragments"),
        _GLOBAL_NOTE,
        SpecComment(("",)),
        SpecComment(("── Data input ──",)),
        SpecField("data_input.fragment_file", placeholder="FRAGMENT_FILE"),
        SpecField("data_input.barcodes_file", placeholder="BARCODES_FILE"),
        SpecComment(
            (
                "sorted_by_barcode: true  # fragments 是否按 barcode 排序；10x 原始文件按位置排序需设 false",
            ),
            indent=1,
        ),
        SpecComment(("",)),
        SpecComment(("── ATAC-specific config ──",)),
        SpecComment(("Reference genome",), indent=1),
        SpecField("atac.genome", placeholder="GENOME"),
        SpecComment(('chrom_sizes: ""',), indent=1),
        SpecComment(('blacklist_bed: ""',), indent=1),
        SpecComment(("",)),
        SpecComment(("QC",), indent=1),
        SpecField("atac.min_fragments", value=1000),
        SpecField("atac.max_fragments", value=50000),
        SpecField("atac.min_tsse", value=7.0),
        SpecComment(("",)),
        SpecComment(("Harmony batch correction",), indent=1),
        SpecField("atac.harmony_use_harmony", value=False),
        SpecField("atac.harmony_batch_key", value="sample"),
        SpecComment(("",)),
        SpecComment(("Multi-metric clustering (stability + gain)",), indent=1),
        SpecField("atac.multi_metric_enabled", value=False),
        SpecComment(
            ("max_blacklist_ratio: 0.05   # 细胞级黑名单比例，当前未实现（TODO）",), indent=1
        ),
        SpecComment(("",)),
        SpecComment(("Peak calling",), indent=1),
        SpecField("atac.peak_qval", value=0.05),
        SpecComment(
            (
                "use_pseudo_replicates: true   # 伪重复峰验证（replicate + replicate_qvalue），false 回退",
            ),
            indent=1,
        ),
        SpecComment(("",)),
        SpecComment(("Feature selection / dimensionality reduction",), indent=1),
        SpecField("atac.n_features", value=50000),
        SpecField("atac.n_spectral", value=30),
        SpecComment(
            ("spectral_sample_size: null   # Nyström 近似阈值；null=禁用（默认），设整数值启用",),
            indent=1,
        ),
        SpecComment(("",)),
        SpecComment(("Differential accessibility",), indent=1),
        SpecField("atac.marker_peaks_log2fc", value=0.5),
        SpecField("atac.marker_peaks_fdr", value=0.05),
        SpecComment(
            (
                'marker_peaks_method: "quick"   # "quick"=snap.marker_regions | "bpc"=pseudobulk+背景匹配Wilcoxon',
            ),
            indent=1,
        ),
        SpecComment(("",)),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", placeholder="TISSUE"),
        SpecComment(
            ('tissue_maturity: "unknown"   # "developing" | "adult" | "unknown"',), indent=0
        ),
        SpecField("species", placeholder="SPECIES"),
        SpecComment(("",)),
        SpecComment(("── Clustering ──",)),
        SpecField("clustering.n_neighbors", value=15),
        SpecField("clustering.param_grid_n_neighbors", value=[15, 20, 30]),
        SpecField("clustering.param_grid_resolutions", value=[0.3, 0.5, 0.8, 1.0, 1.5, 2.0]),
        SpecField("clustering.umap_color_by_batch", value=False),
        SpecComment(
            ("batch_key_override: null  # optionally override integration.batch_key",), indent=1
        ),
        SpecComment(("",)),
        SpecComment(("── Trajectory ──",)),
        SpecComment(
            (
                "umap_maxiter: null                    # max UMAP iterations (null=auto). 200-500 for 1M+ cells.",
            ),
            indent=1,
        ),
        SpecComment(
            ("umap_n_epochs: null                   # UMAP training epochs (null=auto)",), indent=1
        ),
        SpecField("trajectory.root_cell_types", value=[]),
        SpecComment(("terminal_cell_types: []",), indent=1),
        SpecComment(("",)),
        SpecComment(("── Peak-to-gene enrichment ──",)),
        SpecField("enrichment.peak_gene_distance", value=100000),
        SpecComment(('gene_annotation_bed: ""',), indent=1),
        SpecField("enrichment.gene_sets", value=["GO_Biological_Process_2023", "KEGG_2021_Human"]),
        SpecField("enrichment.organism", value="human"),
        SpecComment(("method: both",), indent=1),
        SpecComment(("n_top_genes: 2000",), indent=1),
        SpecComment(("pval_cutoff: 0.05",), indent=1),
        SpecComment(("min_size: 15",), indent=1),
        SpecComment(("max_size: 500",), indent=1),
        SpecComment(("permutations: 1000",), indent=1),
        SpecComment(("",)),
        SpecComment(("── RNA integration (multiome / paired RNA+ATAC) ──",)),
        SpecComment(('rna_h5ad: ""',), indent=0),
        SpecComment(("rna_marker_top_n: 200",), indent=0),
        SpecComment(("rna_marker_pval_threshold: 0.05",), indent=0),
        SpecComment(("rna_marker_logfc_min: 0.25",), indent=0),
        SpecComment(("",)),
        _AI_BLOCK,
        SpecComment(("",)),
        SpecField("h5ad_tempdir", value="/tmp/Fuxi"),
        SpecComment(("h5ad_compression: gzip",), indent=0),
        SpecComment(("cleanup_intermediates: true",), indent=0),
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Spatial spec
# ═══════════════════════════════════════════════════════════════════

_VISIUM = FormatSpec(
    key="visium",
    modality="spatial",
    data_format="visium",
    items=(
        SpecComment(("── Modality & data format ──",)),
        SpecField("modality", value="spatial"),
        SpecField("data_format", value="visium"),
        _GLOBAL_NOTE,
        SpecComment(("",)),
        SpecComment(("── Spatial platform ──",)),
        SpecField("spatial.platform", value="visium"),
        SpecField("spatial.library_id", placeholder="LIBRARY_ID"),
        SpecField("spatial.img_path", value=""),
        SpecField("spatial.crop_image", value=True),
        SpecField("spatial.img_rescale", value=1.0),
        SpecComment(("",)),
        SpecComment(("Spatial graph",), indent=1),
        SpecField("spatial.neighbors_n", value=6),
        SpecField("spatial.neighbors_radius", value=0.0),
        SpecComment(("",)),
        SpecComment(("Spatially variable genes",), indent=1),
        SpecField("spatial.run_autocorr", value=True),
        SpecField("spatial.svg_n_top", value=2000),
        SpecField("spatial.moran_percentile", value=90),
        SpecComment(("run_segmentation: false",), indent=1),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", placeholder="TISSUE"),
        SpecComment(
            ("tissue_kb: \"\"  # auto-inferred from 'tissue' when matching KB exists",), indent=0
        ),
        SpecComment(
            ('tissue_maturity: "unknown"   # "developing" | "adult" | "unknown"',), indent=0
        ),
        SpecField("species", placeholder="SPECIES"),
        SpecComment(("",)),
        SpecComment(("── QC ──",)),
        SpecField("qc.min_genes", value=500),
        SpecField("qc.max_genes", value=7500),
        SpecComment(
            ("max_pct_mito: 20.0  # global default; snRNA-seq→5-10%, heart/liver→15%",), indent=1
        ),
        SpecField("qc.min_genes_per_umi", value=0.7),
        SpecComment(("── Adaptive thresholds (managed by global.yaml) ──",), indent=1),
        SpecComment(("use_adaptive_thresholds: true   # ← global.yaml default",), indent=1),
        SpecComment(("mad_n_mads: 3.0                  # MAD multiplier",), indent=1),
        SpecComment(
            ("ncount_max_mad: 5.0              # nCount upper bound MAD multiplier",), indent=1
        ),
        SpecComment(("",)),
        SpecComment(("── HVG ──",)),
        SpecField("hvg.n_top_genes", value=4000),
        SpecField("hvg.flavor", value="seurat_v3"),
        SpecComment(("use_regress_out: false",), indent=1),
        SpecComment(("forced_genes: []",), indent=1),
        SpecComment(
            ("auto_forced_genes: false        # fill forced_genes from tissue KB when empty",),
            indent=1,
        ),
        SpecComment(("",)),
        SpecComment(("── PCA ──",)),
        SpecField("pca.n_pcs_full", value=100),
        SpecField("pca.n_pcs_use", value=50),
        SpecComment(("",)),
        SpecComment(("── Clustering ──",)),
        SpecField("clustering.n_neighbors", value=30),
        SpecField("clustering.param_grid_resolutions", value=[0.3, 0.5, 0.8, 1.0, 1.5, 2.0]),
        SpecField("clustering.best_resolution", value=1.0),
        SpecField("clustering.best_n_neighbors", value=0),
        SpecField("clustering.leiden_flavor", value="igraph"),
        SpecField("clustering.umap_color_by_batch", value=False),
        SpecComment(
            ("batch_key_override: null  # optionally override integration.batch_key",), indent=1
        ),
        SpecComment(("umap_paga_init: false",), indent=1),
        SpecComment(("",)),
        SpecComment(("── Cell type markers ──",)),
        SpecComment(("TODO: Add known marker genes for {{TISSUE}} tissue.",), indent=0),
        SpecComment(
            (
                "umap_maxiter: null                    # max UMAP iterations (null=auto). 200-500 for 1M+ cells.",
            ),
            indent=1,
        ),
        SpecComment(
            ("umap_n_epochs: null                   # UMAP training epochs (null=auto)",), indent=1
        ),
        SpecField("marker.marker_dict", value={}),
        SpecComment(("",)),
        SpecComment(("── Subcluster (fine-resolution annotation) ──",)),
        SpecComment(("subcluster:",), indent=0),
        SpecComment(("  types: []",), indent=0),
        SpecComment(("  resolution: 0.4",), indent=0),
        SpecComment(("  min_cells: 50",), indent=0),
        SpecComment(("── Differential expression ──",)),
        SpecComment(
            (
                "candidate_pool_expand_steps: [50, 100, 200]   # adaptive scoring tiers → n_genes for DE",
            ),
            indent=1,
        ),
        SpecField("de.stage_pairwise", value=True),
        SpecComment(("auto_switch_on_low_quality: true",), indent=1),
        SpecComment(("",)),
        _PSEUDOBULK_BLOCK,
        SpecComment(("── Trajectory ──",)),
        SpecField("trajectory.root_cell_types", value=[]),
        SpecComment(("root_markers: [SOX2, PAX6, NES]",), indent=1),
        SpecField("trajectory.n_diffmap_comps", value=15),
        SpecField("trajectory.n_branchings", value=2),
        SpecComment(("",)),
        SpecComment(("── Enrichment ──",)),
        SpecField("enrichment.run", value=True),
        SpecField("enrichment.method", value="both"),
        SpecField("enrichment.gene_sets", value=["GO_Biological_Process_2023", "KEGG_2021_Human"]),
        SpecField("enrichment.organism", value="human"),
        SpecComment(("n_top_genes: 2000",), indent=1),
        SpecComment(("pval_cutoff: 0.05",), indent=1),
        SpecComment(("min_size: 15",), indent=1),
        SpecComment(("max_size: 500",), indent=1),
        SpecComment(("permutations: 1000",), indent=1),
        SpecComment(("",)),
        SpecComment(("── GRN regulatory network analysis ──",)),
        SpecComment(("grn:",), indent=0),
        SpecComment(("  run: true",), indent=0),
        SpecComment(("  species: human",), indent=0),
        SpecComment(("  min_regulon_size: 5",), indent=0),
        SpecComment(("  n_top_regulons: 50",), indent=0),
        SpecComment(("",)),
        SpecComment(("── CCI spatial cell-cell interaction analysis ──",)),
        *_CCI,
        SpecComment(("",)),
        _AI_BLOCK,
        SpecComment(("",)),
        SpecComment(("h5ad_compression: gzip",), indent=0),
        SpecComment(("h5ad_tempdir: /tmp/Fuxi",), indent=0),
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Bulk spec
# ═══════════════════════════════════════════════════════════════════

_BULK = FormatSpec(
    key="bulk",
    modality="bulk",
    data_format="count_matrix",
    items=(
        SpecComment(("── Data format ──",)),
        SpecField("data_format", value="count_matrix"),
        SpecComment(("",)),
        SpecComment(("── Modality ──",)),
        SpecField("modality", value="bulk"),
        _GLOBAL_NOTE,
        SpecComment(("",)),
        SpecComment(("── Data input ──",)),
        SpecField("data_input.matrix_file", value=""),
        SpecComment(('metadata_file: "{{METADATA_FILE}}"',), indent=1),
        SpecComment(('csv_sep: ","',), indent=1),
        SpecComment(('gene_symbol_column: ""',), indent=1),
        SpecComment(("",)),
        SpecComment(("── Dataset metadata ──",)),
        SpecField("tissue", placeholder="TISSUE"),
        SpecField("species", placeholder="SPECIES"),
        SpecField("expression_type", value="raw_counts"),
        SpecComment(("raw_counts | tpm | fpkm | cpm",), indent=0),
        SpecComment(("",)),
        SpecComment(("── Bulk-specific settings ──",)),
        SpecField("bulk.design", value="~condition"),
        SpecComment(('contrast_column: "condition"',), indent=1),
        SpecComment(('contrast_treatment: "treated"',), indent=1),
        SpecComment(('contrast_baseline: "control"',), indent=1),
        SpecField("bulk.alpha", value=0.05),
        SpecField("bulk.lfc_shrink", value=True),
        SpecField("bulk.normalization_method", value="deseq2_median_ratios"),
        SpecComment(("min_counts_per_gene: 10",), indent=1),
        SpecComment(("min_samples_per_group: 2",), indent=1),
        SpecComment(("n_jobs: 0",), indent=1),
        SpecComment(('output_dir: ""',), indent=1),
        SpecComment(("batch_correct: false",), indent=1),
        SpecComment(('batch_column: "batch"',), indent=1),
        SpecComment(("",)),
        SpecComment(("── Enrichment ──",)),
        SpecField("enrichment.run", value=True),
        SpecField("enrichment.method", value="both"),
        SpecField("enrichment.gene_sets", value=["GO_Biological_Process_2023", "KEGG_2021_Human"]),
        SpecField("enrichment.organism", value="human"),
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Modality → format → spec lookup (two-level)
# ═══════════════════════════════════════════════════════════════════

MODALITY_SPECS: Dict[str, Dict[str, Union[FormatSpec, str]]] = {
    "rna": {
        "default": "10X_h5",
        "10X_h5": _10X_H5,
        "10X_mtx": _10X_MTX,
        "csv_matrix": _CSV,
        "preprocessed": _PREPROCESSED,
        "h5ad": "10X_h5",  # reuse
    },
    "atac": {
        "default": "10x_fragments",
        "10x_fragments": _FRAGMENTS,
        "10x_peak_h5": "10x_fragments",  # reuse
    },
    "spatial": {
        "default": "visium",
        "visium": _VISIUM,
    },
    "bulk": {
        "default": "count_matrix",
        "count_matrix": _BULK,
        "tpm_matrix": _BULK,  # reuse
        "bulk_h5ad": _BULK,  # reuse
    },
    # multiome inherits the rna mapping (RNA is the primary modality)
}


def lookup_spec(modality: str, data_format: str) -> Optional[FormatSpec]:
    """Resolve the spec for *(modality, data_format)*.

    - ``data_format == "unknown"`` → the modality's ``default`` spec.
    - Unknown modality (e.g. ``multiome``) → inherits the rna mapping.
    - Format not present in the modality map → falls back to any other
      modality that defines it (preserves old cross-format behaviour,
      e.g. fragments-only multiome).
    - Returns None when nothing matches.
    """
    mod = MODALITY_SPECS.get(modality) or MODALITY_SPECS["rna"]
    if data_format == "unknown":
        data_format = str(mod["default"])
    target: Optional[Union[FormatSpec, str]] = mod.get(data_format)
    if target is None:
        for other in MODALITY_SPECS.values():
            cand = other.get(data_format)
            if cand is None:
                continue
            while isinstance(cand, str):
                cand = other.get(cand)
            if isinstance(cand, FormatSpec):
                return cand
        return None
    while isinstance(target, str):
        target = mod.get(target)
    return target


def materialized_specs() -> list[FormatSpec]:
    """Distinct specs (one per committed template file), in file order."""
    seen: set[str] = set()
    out: list[FormatSpec] = []
    for mod in MODALITY_SPECS.values():
        for target in mod.values():
            if isinstance(target, FormatSpec) and target.key not in seen:
                seen.add(target.key)
                out.append(target)
    return out


# ═══════════════════════════════════════════════════════════════════
# Validation (schema drift guard)
# ═══════════════════════════════════════════════════════════════════


def _path_is_valid(path: str) -> bool:
    """True if *path* resolves against the Config model tree.

    Traversal stops at the first non-model (scalar / dict / list) field —
    any remaining components are free-form keys inside that field (e.g.
    ``integration.scvi.datasplitter_kwargs.num_workers``).
    """
    parts = path.split(".")
    model: type[BaseModel] = Config
    for i, part in enumerate(parts):
        if part not in model.model_fields:
            return False
        field_info = model.model_fields[part]
        inner = _resolve_base_model(field_info.annotation)
        if inner is None:
            return True  # free-form sub-keys allowed below a non-model field
        model = inner
    return True


def validate_specs() -> list[str]:
    """Return a list of spec/schema inconsistencies (empty = healthy).

    Checks: every ``SpecField.path`` exists in the schema; every
    placeholder is supported; every reuse reference resolves; every
    ``FormatSpec.modality`` matches its owning modality map.
    """
    errors: list[str] = []
    for mod_key, mod in MODALITY_SPECS.items():
        if "default" not in mod:
            errors.append(f"[{mod_key}] missing 'default' format")
            continue
        if mod["default"] not in mod:
            errors.append(f"[{mod_key}] default '{mod['default']}' not in map")
        for fmt, target in mod.items():
            if fmt == "default":
                continue
            if isinstance(target, str):
                if target not in mod:
                    errors.append(f"[{mod_key}] reuse '{fmt}' → missing '{target}'")
                continue
            if not isinstance(target, FormatSpec):
                errors.append(f"[{mod_key}] entry '{fmt}' is not a FormatSpec")
                continue
            if target.modality != mod_key:
                errors.append(
                    f"[{mod_key}] spec '{target.key}' modality='{target.modality}' mismatch"
                )
            for item in target.items:
                if not isinstance(item, SpecField):
                    continue
                if not _path_is_valid(item.path):
                    errors.append(f"[{target.key}] unknown schema path: {item.path}")
                if item.placeholder and item.placeholder not in SUPPORTED_PLACEHOLDERS:
                    errors.append(f"[{target.key}] unsupported placeholder: {item.placeholder}")
    return errors
