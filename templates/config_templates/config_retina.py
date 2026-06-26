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
STAGE_MAP = {
    'Hgw9':  'EarlyFetal',  'Hgw11': 'EarlyFetal',  'Hgw12': 'EarlyFetal',
    'Hgw13': 'EarlyFetal',
    'Hgw14': 'MidFetal',    'Hgw15': 'MidFetal',    'Hgw16': 'MidFetal',
    'Hgw17': 'MidFetal',    'Hgw18': 'MidFetal',    'Hgw19': 'MidFetal',
    'Hgw20': 'MidFetal',
    'Hgw22': 'LateFetal',   'Hgw24': 'LateFetal',   'Hgw27': 'LateFetal',
    'Hpnd8': 'Postnatal',
    '24_Day': 'Organoid',   '30_Day': 'Organoid',   '42_Day': 'Organoid',
    '59_Day': 'Organoid',
    'Adult': 'Adult',
}
# 注意: 此处 stage_map 为 dict[int, str] 格式，而实际使用时 age 值
# 可能是字符串，因此使用 meta_columns + 步骤 04 中的手动映射。
CFG.stage_map = {}
CFG.stage_order = ['Organoid', 'EarlyFetal', 'MidFetal', 'LateFetal', 'Postnatal', 'Adult']

# ── QC ──
CFG.min_genes = 200
CFG.max_genes = 7500
CFG.max_pct_mito = 20.0
CFG.min_cells_per_gene = 3
CFG.run_scrublet = True

# ── HVG ──
CFG.n_top_genes = 5000
CFG.hvg_batch_key = 'sample'

# ── 批次校正 ──
CFG.harmony_batch_key = 'age'

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

# ── 轨迹根细胞使用标记基因自动检测 ──
CFG.root_markers = ['VSX2', 'PAX6', 'SOX2', 'HES1', 'NOTCH1']

# ── 执行 ──
CFG.n_jobs = 0  # auto-detect (override in project config if needed)
CFG.random_seed = 42
