#!/usr/bin/env python3
"""config_test_10w.py — 测试配置模板 (ATAC)"""

from core.config import Config, CFG
import os
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
# CFG.chrom_sizes = ''
# CFG.blacklist_bed = ''
# CFG.tss_bed = ''

# ── 降采样（测试用，只取前500细胞） ──
CFG.max_cells = 500

# ── QC ──
CFG.min_fragments = 1000
CFG.max_fragments = 50000
CFG.min_tsse = 7.0
# CFG.max_blacklist_ratio = 0.05
# CFG.min_peak_region_fragments = 300

# ── Peak calling ──
CFG.peak_qval = 0.05
CFG.peak_width = 500
# CFG.use_macs3 = True

# ── 降维/聚类 ──
CFG.n_features = 10000
CFG.n_spectral = 30
CFG.n_neighbors = 15
CFG.param_grid_n_neighbors = [15, 20]
CFG.param_grid_resolutions = [0.3, 0.8, 1.5]
# UMAP sweep (only meaningful for RNA/Spatial; ATAC umap is unexposed)
CFG.umap_selection_method = "convex_hull"
CFG.param_grid_min_dist = [0.1, 0.3, 0.5]
CFG.param_grid_spread = [1.0]
CFG.umap_min_dist = 0.3    # ATAC 也使用 UMAP
CFG.umap_spread = 1.0

# ── 差异分析 ──
CFG.marker_peaks_log2fc = 0.25
CFG.marker_peaks_fdr = 0.1

# ── 轨迹 ──
CFG.root_cell_types = []
# CFG.terminal_cell_types = []

# ── 富集 ──
CFG.peak_gene_distance = 50000
# CFG.gene_annotation_bed = ''
# CFG.enrichment_method = 'both'
# CFG.enrichment_gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
# CFG.enrichment_organism = 'human'
# CFG.enrichment_n_top_genes = 2000
# CFG.enrichment_pval_cutoff = 0.05
# CFG.enrichment_min_size = 15
# CFG.enrichment_max_size = 500
# CFG.enrichment_permutations = 1000

# ── CCI 细胞通讯分析 (ATAC 暂无对应步骤，占位) ──
# CFG.run_cci = True
# CFG.cci_method = "liana"                  # 'liana' (LIANA+ rank_aggregate)
# CFG.cci_lr_database = "consensus"         # 'consensus' | 'cellphonedb' | 'cellchat' | 'celltalkdb' | 'ramilowski'
# CFG.cci_permutations = 1000               # permutation 迭代次数
# CFG.cci_n_top_interactions = 50           # 热图展示的 top 互作对数量

# ── RNA 整合 ──
# CFG.rna_h5ad = ''
# CFG.rna_marker_top_n = 200
# CFG.rna_marker_pval_threshold = 0.05
# CFG.rna_marker_logfc_min = 0.25

# ── AI ──
CFG.ai.enabled = False

# ── 执行 ──
CFG.n_jobs = 4
CFG.random_seed = 42
# CFG.force_csr = True
# CFG.use_float32 = False
# CFG.limit_blas_threads = True
# CFG.scanpy_verbosity = 2
# CFG.h5ad_compression = 'gzip'
# CFG.h5ad_tempdir = '/tmp/Fuxi'
# CFG.cleanup_intermediates = True
