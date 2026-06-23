#!/usr/bin/env python3
"""
GSE169109 配置模板 — 人胎下丘脑发育 scRNA-seq
================================================
这是通用管线的具体配置示例，对应 GSE169109 项目。

使用方法:
    cp config_hypothalamus.py ../../config_myproject.py
    python ../../run_pipeline.py --config config_myproject.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.config import CFG

# ── 数据格式 ──
CFG.data_format = '10X_mtx'
CFG.mtx_prefix = 'GSE169109_'

# ── 样本映射 (19 个样本, GW7~GW20) ──
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

# ── 批次校正 ──
CFG.harmony_batch_key = 'sample'

# ── 聚类 ──
CFG.leiden_resolutions = [0.3, 0.5, 0.8]
CFG.best_resolution = 0.8

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

# ── 轨迹 ──
CFG.root_cell_types = ['NE', 'NP']

# ── 执行 ──
CFG.n_jobs = 0  # auto-detect (override in project config if needed)
CFG.random_seed = 42
