#!/usr/bin/env python3
"""
10X HDF5 (.h5) 配置模板 — 通用 scRNA-seq
==========================================
这是针对 10X Genomics HDF5 格式（.h5 文件）的通用配置模板。
适用于 Cell Ranger 输出的 filtered_feature_bc_matrix.h5 或 raw_feature_bc_matrix.h5。

使用方法:
    cp config_10X_h5.py ../../config_myproject.py
    编辑 config_myproject.py，填写样本映射、阶段映射和标记基因
    python ../../run_pipeline.py --config config_myproject.py

数据准备:
    将所有 .h5 文件放在同一目录下（默认与 config 文件同级目录）。
    文件名示例: sample1_filtered_feature_bc_matrix.h5, sample2_filtered_feature_bc_matrix.h5
"""
import os
from core.config import CFG

# ── 数据格式 ──
CFG.data_format = '10X_h5'

# .h5 文件匹配模式（glob 语法）
# 默认匹配 Cell Ranger 过滤后的输出: "*filtered_feature_bc_matrix.h5"
# 如需使用原始矩阵，改为: "*raw_feature_bc_matrix.h5"
# 如需自定义，改为: "*_counts.h5" 等
CFG.h5_file_pattern = "*filtered_feature_bc_matrix.h5"

# .h5 文件所在目录（留空则默认使用 data_dir，即 config.py 所在目录）
# CFG.h5_dir = "data"

# ── 样本映射 ──
# 当存在多个 .h5 文件时，管线会自动从文件名提取样本名。
# 如果文件名不符合 pattern，或需要自定义样本名，可在此处设置 barcode 后缀映射。
# 10X 多通道聚合时，barcode 以 -1, -2, ... 后缀区分样本。
# 示例: {1: 'sample1', 2: 'sample2'}
CFG.sample_map = {
    # 1: 'sample1',
    # 2: 'sample2',
}

# ── 阶段映射 ──
# barcode 后缀 → 发育阶段 / 实验分组
# 示例: {1: 'Day3', 2: 'Day7'}
STAGE_MAP = {
    # 1: 'Day3',
    # 2: 'Day7',
}
CFG.stage_map = STAGE_MAP
CFG.stage_order = []  # 阶段顺序，用于图例排序

# ── 数据集元信息 ──
CFG.tissue = "unknown"   # 组织类型: brain, retina, liver, etc.
CFG.species = "human"    # 物种: human, mouse, rat
CFG.expression_type = "raw_counts"   # 10X H5 默认是 raw UMI counts
# raw_counts | log1p_counts | TPM | CPM | FPKM
# TPM/FPKM/CPM: total_counts & complexity filters are auto-skipped.
# Scrublet is auto-disabled for non-raw_counts data.

# ── QC ──
CFG.min_genes = 500
CFG.max_genes = 7500
CFG.max_pct_mito = 20.0
CFG.min_cells_per_gene = 3
CFG.run_scrublet = True           # auto-disabled when expression_type != "raw_counts"
CFG.min_genes_per_umi = 0.7       # complexity filter — only applied when expression_type="raw_counts"
# CFG.use_adaptive_thresholds = True   # 替代固定阈值，基于 MAD
# CFG.mad_n_mads = 3.0
# CFG.qc_ncount_max_mad = 5.0

# ── HVG ──
CFG.n_top_genes = 4000
CFG.hvg_batch_key = 'sample'
CFG.hvg_flavor = 'seurat_v3'
# CFG.use_regress_out = False     # 回归细胞周期/n_counts 等协变量

# ── 批次校正 / Harmony ──
CFG.harmony_batch_key = 'sample'
CFG.use_harmony = True
CFG.harmony_max_iter = 20

# ── PCA ──
CFG.n_pcs_full = 100
CFG.n_pcs_use = 50

# ── 聚类 ──
CFG.leiden_resolutions = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
CFG.best_resolution = 1.0
CFG.umap_min_dist = 0.3    # 增大（如0.5）可增强 UMAP 散点延展性
CFG.umap_spread = 1.0

# ── 已知细胞类型标记 ──
# 格式: {'CellType': ['marker1', 'marker2', ...]}
# 请根据组织类型替换为实际标记基因
CFG.marker_dict = {
    # 'Progenitor': ['SOX2', 'PAX6', 'NES', 'VIM'],
    # 'Neuron': ['RBFOX3', 'MAP2', 'DCX', 'TUBB3'],
    # 'Astrocyte': ['GFAP', 'AQP4', 'S100B'],
    # 'Microglia': ['PTPRC', 'CX3CR1', 'CSF1R'],
    # 'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
}

# ── 知识库 ──
# CFG.tissue_kb = ''   # rna/tissue_ontologies/ 中的组织 KB 名称，如 'retina'

# ── 子聚类 ──
# 对哪些细胞类型进行进一步子聚类
CFG.subcluster_types = []
CFG.subcluster_resolution = 0.4
CFG.min_cells_subcluster = 50

# ── 轨迹 ──
# 根细胞类型（发育起点）
CFG.root_cell_types = []
# CFG.root_markers = ['SOX2', 'PAX6', 'NES']   # 通过标记基因自动检测根
# CFG.n_diffmap_comps = 15
# CFG.n_branchings = 2

# ── 差异表达 ──
CFG.de_method = 'wilcoxon'
CFG.de_n_genes = 50
CFG.de_pval_cutoff = 0.05
CFG.de_logfc_cutoff = 0.25
CFG.de_stage_pairwise = True
# CFG.de_auto_switch_on_low_quality = True   # 低质量数据自动切换方法

# ── 富集分析 (取消注释以启用) ──
# CFG.run_enrichment = True
# CFG.enrichment_method = 'both'       # 'ora' | 'prerank' | 'both'
# CFG.enrichment_gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
# CFG.enrichment_organism = 'human'

# ── GRN 调控网络分析 (取消注释以启用) ──
CFG.run_grn = True
CFG.grn_method = "decoupler"           # 'decoupler' only for now (pySCENIC TBD)
CFG.grn_species = "human"              # 'human' | 'mouse'
CFG.grn_n_top_regulons = 50            # top N variable TFs for heatmap
CFG.grn_min_regulon_size = 5           # minimum target genes per regulon
# CFG.grn_confidence_levels = ["A", "B", "C"]  # DoRothEA confidence levels

# ── AI 设置 (取消注释以启用) ──
# CFG.ai.enabled = True
# CFG.ai.api_base = 'https://api.deepseek.com/v1'
# CFG.ai.model = 'deepseek-v4-pro'
# CFG.ai.api_key = os.environ.get('LLM_API_KEY', '')
# CFG.ai.max_tokens = 32768
# CFG.ai.temperature = 0.1
# CFG.ai.thinking_enabled = True
# CFG.ai.reasoning_effort = 'high'
# CFG.ai.ai_annotation = True
# CFG.ai.ai_subcluster = True
# CFG.ai.ai_interpretation = True
# CFG.ai.ai_cache_responses = True

# ── 执行 ──
CFG.n_jobs = 0  # auto-detect (override in project config if needed)
CFG.random_seed = 42
# CFG.force_csr = True           # 强制 CSR 稀疏矩阵格式
# CFG.use_float32 = False        # 使用 float32 减少内存
# CFG.limit_blas_threads = True  # 限制 BLAS 线程数
# CFG.scanpy_verbosity = 2       # Scanpy 日志级别 (0=error, 1=warning, 2=info)
# CFG.h5ad_compression = 'gzip'  # h5ad 写入压缩

# ── 降采样 (可选) ──
# CFG.downsample_target = 5000         # 每样本最多保留 N 细胞
# CFG.downsample_strategy = 'sample'   # 'sample' | 'random'
