#!/usr/bin/env python3
"""
配置模板 — 人胎下丘脑发育 scRNA-seq (10X MTX)
================================================
可复制到项目目录并修改各参数适配具体数据集。

使用方法:
    cp config_hypothalamus.py ../../config_myproject.py
    python ../../run_pipeline.py --config config_myproject.py
"""
import os
from core.config import CFG

# ── 数据格式 ──
CFG.data_format = '10X_mtx'
CFG.mtx_prefix = '<GSE_ID>_'

# ── 样本映射 (示例: 19 个样本, GW7~GW20) ──
CFG.sample_map = {
    1: 'GW7-lane1',   2: 'GW7-lane2',
    3: 'GW8-1',       4: 'GW8-2',
    5: 'GW10',
    6: 'GW12_01',     7: 'GW12_02',
    8: 'GW15-A',      9: 'GW15-M',     10: 'GW15-P',
    11: 'GW18-01-A',  12: 'GW18-01-M', 13: 'GW18-01-P',
    14: 'GW18-02-lane1', 15: 'GW18-02-lane2', 16: 'GW18-02-lane3',
    17: 'GW20-A',     18: 'GW20-M',     19: 'GW20-P',
}

# ── 阶段映射 ──
STAGE_MAP = {
    1: 'GW7',  2: 'GW7',
    3: 'GW8',  4: 'GW8',
    5: 'GW10',
    6: 'GW12', 7: 'GW12',
    8: 'GW15', 9: 'GW15', 10: 'GW15',
    11: 'GW18', 12: 'GW18', 13: 'GW18',
    14: 'GW18', 15: 'GW18', 16: 'GW18',
    17: 'GW20', 18: 'GW20', 19: 'GW20',
}
CFG.stage_map = STAGE_MAP
CFG.stage_order = ['GW7', 'GW8', 'GW10', 'GW12', 'GW15', 'GW18', 'GW20']

# ── QC ──
CFG.min_genes = 500
CFG.max_genes = 7000
CFG.max_pct_mito = 20.0
CFG.min_genes_per_umi = 0.7
CFG.run_scrublet = True

# ── HVG ──
CFG.n_top_genes = 4000
CFG.hvg_batch_key = 'sample'
CFG.hvg_flavor = 'seurat_v3'
# CFG.use_regress_out = False

# ── 批次校正 / Harmony ──
CFG.harmony_batch_key = 'sample'
CFG.use_harmony = True
CFG.harmony_max_iter = 20

# ── PCA ──
CFG.n_pcs_full = 100
CFG.n_pcs_use = 50

# ── 聚类 ──
CFG.leiden_resolutions = [0.3, 0.5, 0.8]
CFG.best_resolution = 0.8   # Only used when cluster_selection_method is None
CFG.best_n_neighbors = 0     # Only used when cluster_selection_method is None; 0 = auto-pick
CFG.cluster_selection_method = "pareto_elbow"  # "pareto_elbow" | "silhouette" | None
# UMAP visualization parameter sweep — tried AFTER best cluster params selected.
# Selection method:  "convex_hull" (auto, default) | None (manual: use umap_min_dist/spread)
CFG.umap_selection_method = "convex_hull"
CFG.param_grid_min_dist = [0.1, 0.3, 0.5]  # values to sweep
CFG.param_grid_spread = [1.0]
CFG.umap_min_dist = 0.3    # 增大（如0.5）可增强 UMAP 散点延展性
CFG.umap_spread = 1.0

# ── 已知下丘脑细胞类型标记 ──
CFG.marker_dict = {
    'NE':           ['HES1', 'HES5', 'SOX2', 'NES', 'NOTCH1'],
    'NP':           ['SOX2', 'PAX6', 'ASCL1', 'MKI67', 'TOP2A'],
    'Neuron':       ['RBFOX3', 'MAP2', 'SYN1', 'DCX', 'STMN2', 'TUBB3'],
    'Astrocyte':    ['GFAP', 'AQP4', 'S100B', 'ALDH1L1', 'GJA1'],
    'OPC':          ['PDGFRA', 'CSPG4', 'SOX10', 'OLIG1', 'OLIG2'],
    'OL':           ['MBP', 'MOG', 'PLP1', 'MAG', 'CLDN11'],
    'Ependymocyte': ['FOXJ1', 'RSPH1', 'DNAH9', 'CFAP53', 'TUBB4B'],
    'Microglia':    ['PTPRC', 'CSF1R', 'CX3CR1', 'TREM2', 'ITGAM'],
    'Endothelial':  ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'FLT1'],
    'Mural':        ['RGS5', 'PDGFRB', 'ACTA2', 'MYH11', 'CSPG4'],
    'VLMC':         ['LUM', 'DCN', 'COL1A1', 'COL3A1', 'FN1'],
}

# ── 子聚类 ──
CFG.subcluster_types = ['NE', 'NP', 'Neuron']
CFG.subcluster_resolution = 0.4
CFG.min_cells_subcluster = 50

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
CFG.root_cell_types = ['NE', 'NP']
# CFG.root_markers = ['SOX2', 'PAX6', 'NES']
# CFG.n_diffmap_comps = 15
# CFG.n_branchings = 2

# ── 富集分析 (取消注释以启用) ──
# CFG.run_enrichment = True
# CFG.enrichment_method = 'both'
# CFG.enrichment_gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
# CFG.enrichment_organism = 'human'

# ── GRN 调控网络分析 (取消注释以启用) ──
CFG.run_grn = True
CFG.grn_method = "decoupler"           # 目前仅 'decoupler' (pySCENIC 待定)
CFG.grn_species = "human"              # 'human' | 'mouse'
CFG.grn_n_top_regulons = 50            # 热图显示的方差最高 TF 数量
CFG.grn_min_regulon_size = 5           # 每个 regulon 最少靶基因数
# CFG.grn_confidence_levels = ["A", "B", "C"]  # DoRothEA 置信度等级

# ── CCI 细胞通讯分析 (取消注释以启用) ──
CFG.run_cci = True
# CFG.cci_method = "liana"                  # 'liana' (LIANA+ rank_aggregate)
# CFG.cci_lr_database = "consensus"         # 'consensus' | 'cellphonedb' | 'cellchat' | 'celltalkdb' | 'ramilowski'
# CFG.cci_permutations = 1000               # permutation 迭代次数
# CFG.cci_n_top_interactions = 50           # 热图展示的 top 互作对数量
# CFG.cci_spatial_method = "liana_spatial"  # spatial 方法 (reserves 'commot')
# CFG.cci_spatial_distance = 0.0            # 0 = use existing spatial_connectivities
# CFG.cci_lr_cache_dir = ""                 # LIANA cache 目录; 空 = auto (~/.cache/liana)
# # CFG.cci_multi_condition = False         # future: multi-condition CCI

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
