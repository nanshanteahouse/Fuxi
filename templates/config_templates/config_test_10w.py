#!/usr/bin/env python3
"""config_test_10w.py — 测试配置模板 (ATAC)"""

from core.config import Config, CFG
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import data_root

# ── 数据路径 ──
CFG.data_dir = os.path.join(data_root(), "<GSE_ID>")
CFG.data_format = "10x_fragments"
CFG.fragment_file = os.path.join(CFG.data_dir, "<fragments.tsv.gz>")
CFG.barcodes_file = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "<barcodes.txt>",
)

# ── 参考基因组 ──
CFG.genome = "hg38"

# ── 降采样（测试用，只取前500细胞） ──
CFG.max_cells = 500

# ── QC ──
CFG.min_fragments = 1000
CFG.max_fragments = 50000
CFG.min_tsse = 7.0

# ── Peak calling ──
CFG.peak_qval = 0.05
CFG.peak_width = 500

# ── 降维/聚类 ──
CFG.n_features = 10000
CFG.n_spectral = 30
CFG.n_neighbors = 15
CFG.param_grid_n_neighbors = [15, 20]
CFG.param_grid_resolutions = [0.3, 0.8, 1.5]
CFG.umap_min_dist = 0.3    # ATAC 也使用 UMAP
CFG.umap_spread = 1.0

# ── 差异分析 ──
CFG.marker_peaks_log2fc = 0.25
CFG.marker_peaks_fdr = 0.1

# ── 轨迹 ──
CFG.root_cell_types = []

# ── 富集 ──
CFG.peak_gene_distance = 50000

# ── AI ──
CFG.ai.enabled = False

# ── 执行 ──
CFG.n_jobs = 4
CFG.random_seed = 42
