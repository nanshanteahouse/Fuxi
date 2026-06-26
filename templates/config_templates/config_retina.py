#!/usr/bin/env python3
"""
配置模板 — 人视网膜发育 scRNA-seq (CSV 矩阵)
================================================
可复制到项目目录并修改各参数适配具体数据集。

使用方法:
    cp config_retina.py ../../config_myproject.py
    python ../../run_pipeline.py --config config_myproject.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.config import CFG

# ── 数据格式 (CSV 矩阵) ──
CFG.data_format = 'csv_matrix'
CFG.matrix_file = '<GSE_ID>_matrix.mtx.gz'
CFG.barcodes_file = '<GSE_ID>_barcodes.csv.gz'
CFG.features_file = '<GSE_ID>_genes.csv.gz'

# ── CSV 中已有的元数据列映射 ──
CFG.meta_columns = {
    'sample': 'sample',
    'stage': 'age',
    'tissue': 'sample_type',
}

# ── 阶段分组映射 ──
CFG.stage_map = {
    # 注意: 此处 stage 来自 meta_columns 中 'stage' 列的原始值，
    # 实际映射在步骤 04 中使用 pandas map 手动完成。
    # 参见 rna/steps/04_cluster_umap.py 中的 _map_stages() 函数。
}
CFG.stage_order = ['Organoid', 'EarlyFetal', 'MidFetal', 'LateFetal', 'Postnatal', 'Adult']

# ── QC ──
CFG.min_genes = 200
CFG.max_genes = 7500
CFG.max_pct_mito = 20.0
CFG.min_cells_per_gene = 3
CFG.run_scrublet = True
# CFG.use_adaptive_thresholds = True   # 替代固定阈值，基于 MAD
# CFG.mad_n_mads = 3.0

# ── HVG ──
CFG.n_top_genes = 5000
CFG.hvg_batch_key = 'sample'
CFG.hvg_flavor = 'seurat_v3'
# CFG.use_regress_out = False

# ── 批次校正 / Harmony ──
CFG.harmony_batch_key = 'age'
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

# ── 已知视网膜细胞类型标记 ──
CFG.marker_dict = {
    'RPCs':              ['VSX2', 'PAX6', 'SOX2', 'HES1', 'NOTCH1'],
    'Neurogenic Cells':  ['ASCL1', 'NEUROG2', 'TUBB3', 'DCX', 'STMN2'],
    'RGCs':              ['POU4F2', 'POU4F1', 'NEFM', 'NEFL', 'ELAVL4'],
    'Cones':             ['ARR3', 'OPN1SW', 'GNAT2', 'PDE6C', 'RCVRN'],
    'Rods':              ['NRL', 'RHO', 'GNAT1', 'PDE6B', 'SAG'],
    'Horizontal':        ['ONECUT1', 'ONECUT2', 'PROX1', 'CALB1', 'LHX1'],
    'Bipolar':           ['VSX1', 'PRKCA', 'TRPM1', 'GRM6', 'CABP5'],
    'Amacrine':          ['GAD1', 'GAD2', 'SLC6A9', 'TFAP2A', 'TFAP2B'],
    'Muller Glia':       ['RLBP1', 'GFAP', 'VIM', 'SLC1A3', 'GLUL'],
    'BC/Photo_Precurs':  ['CRX', 'OTX2', 'VSX1'],
    'AC/HC_Precurs':     ['PTF1A', 'TFAP2A', 'PROX1', 'ONECUT1'],
    'Microglia':         ['AIF1', 'P2RY12', 'CSF1R', 'CX3CR1', 'CD74'],
    'Astrocytes':        ['GFAP', 'AQP4', 'ALDH1L1', 'S100B', 'SOX9'],
    'Pericytes':         ['PDGFRB', 'CSPG4', 'RGS5', 'ANPEP', 'ACTA2'],
    'RPE':               ['RPE65', 'BEST1', 'LRAT', 'TIMP3', 'RDH10'],
    'Endothelial':       ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'EGFL7'],
    'Oligodendrocytes':  ['MBP', 'PLP1', 'MOG', 'OLIG2', 'SOX10'],
    'Fibroblast':        ['COL1A1', 'DCN', 'LUM', 'COL3A1', 'COL1A2'],
}

# ── 子聚类 ──
# CFG.subcluster_types = ['RPCs', 'Muller Glia']
# CFG.subcluster_resolution = 0.4
# CFG.min_cells_subcluster = 50

# ── 知识库 ──
# CFG.tissue_kb = ''   # rna/tissue_ontologies/ 中的组织 KB 名称

# ── 差异表达 ──
CFG.de_method = 'wilcoxon'
CFG.de_n_genes = 50
CFG.de_pval_cutoff = 0.05
CFG.de_logfc_cutoff = 0.25
CFG.de_stage_pairwise = True
# CFG.de_auto_switch_on_low_quality = True

# ── 轨迹 ──
CFG.root_markers = ['VSX2', 'PAX6', 'SOX2', 'HES1', 'NOTCH1']
# CFG.root_cell_types = []
# CFG.n_diffmap_comps = 15
# CFG.n_branchings = 2

# ── 富集分析 (取消注释以启用) ──
# CFG.run_enrichment = True
# CFG.enrichment_method = 'both'
# CFG.enrichment_gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
# CFG.enrichment_organism = 'human'

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
# CFG.force_csr = True
# CFG.use_float32 = False
# CFG.limit_blas_threads = True
# CFG.scanpy_verbosity = 2
# CFG.h5ad_compression = 'gzip'

# ── 降采样 (可选) ──
# CFG.downsample_target = 5000
# CFG.downsample_strategy = 'sample'
