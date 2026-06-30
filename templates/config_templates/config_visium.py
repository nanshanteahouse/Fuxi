#!/usr/bin/env python3
"""
10X Visium (SpaceRanger output) config template — spatial transcriptomics
==========================================================================
This template is for SpaceRanger output directories containing:
  - filtered_feature_bc_matrix.h5  (gene expression)
  - spatial/tissue_positions_list.csv  (spot coordinates)
  - spatial/tissue_hires_image.png     (H&E image)
  - spatial/scalefactors_json.json     (scaling factors)

Usage:
    python core/run_pipeline.py --modality spatial --config this_config.py
"""
import os
from core.config import CFG

# ── Modality ──
CFG.modality = 'spatial'
CFG.data_format = 'visium'

# ── Spatial platform ──
CFG.spatial_platform = 'visium'
CFG.library_id = '{{LIBRARY_ID}}'   # e.g. 'V1_Adult_Mouse_Brain'
CFG.img_path = ''                    # Auto-detected from Visium directory

# ── Dataset metadata ──
CFG.tissue = '{{TISSUE}}'   # TODO: verify tissue type
CFG.species = '{{SPECIES}}'

# ── QC ──
CFG.min_genes = 500
CFG.max_genes = 7500
CFG.max_pct_mito = 20.0
CFG.min_genes_per_umi = 0.7
CFG.min_cells_per_gene = 3
# CFG.use_adaptive_thresholds = True
# CFG.mad_n_mads = 3.0

# ── HVG ──
CFG.n_top_genes = 4000
CFG.hvg_flavor = 'seurat_v3'
# CFG.use_regress_out = False

# ── Spatial graph ──
CFG.spatial_neighbors_n = 6
CFG.spatial_neighbors_radius = 0.0   # Use n_neighbors mode

# ── PCA ──
CFG.n_pcs_full = 100
CFG.n_pcs_use = 50

# ── Clustering ──
CFG.n_neighbors = 30
CFG.leiden_resolutions = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
CFG.best_resolution = 1.0   # Only used when cluster_selection_method is None
CFG.best_n_neighbors = 0     # Only used when cluster_selection_method is None; 0 = auto-pick
CFG.cluster_selection_method = "pareto_elbow"  # "pareto_elbow" | "silhouette" | None
CFG.leiden_flavor = 'igraph'
# UMAP visualization parameter sweep — tried AFTER best cluster params selected.
# Selection method:  "convex_hull" (auto, default) | None (manual: use umap_min_dist/spread)
CFG.umap_selection_method = "convex_hull"
CFG.param_grid_min_dist = [0.1, 0.3, 0.5]  # values to sweep
CFG.param_grid_spread = [1.0]
CFG.umap_min_dist = 0.3    # increase (e.g. 0.5) for more UMAP spread
CFG.umap_spread = 1.0

# ── Cell type markers ──
# TODO: Add known marker genes for {{TISSUE}} tissue.
CFG.marker_dict = {
    # 'CellTypeA': ['GENE1', 'GENE2', 'GENE3'],
    # 'CellTypeB': ['GENE4', 'GENE5', 'GENE6'],
}

# ── Knowledge Base ──
CFG.tissue_kb = ''   # e.g. 'retina', 'hypothalamus'

# ── Spatially variable genes ──
CFG.run_spatial_autocorr = True
CFG.svg_n_top = 2000
CFG.moran_percentile = 90

# ── Differential expression ──
CFG.de_method = 'wilcoxon'
CFG.de_n_genes = 50
CFG.de_pval_cutoff = 0.05
CFG.de_logfc_cutoff = 0.25

# ── Trajectory ──
CFG.root_cell_types = []   # TODO: developmental root cell types
# CFG.root_markers = ['SOX2', 'PAX6', 'NES']
CFG.n_diffmap_comps = 15
CFG.n_branchings = 2

# ── Differential expression ──
CFG.de_method = 'wilcoxon'
CFG.de_n_genes = 50
CFG.de_pval_cutoff = 0.05
CFG.de_logfc_cutoff = 0.25
CFG.de_stage_pairwise = True
# CFG.de_auto_switch_on_low_quality = True

# ── Enrichment ──
CFG.run_enrichment = True
CFG.enrichment_method = 'both'   # 'ora' | 'prerank' | 'both'
CFG.enrichment_gene_sets = [
    'GO_Biological_Process_2023',
    'KEGG_2021_Human',
]
CFG.enrichment_organism = 'human'
# CFG.enrichment_n_top_genes = 2000
# CFG.enrichment_pval_cutoff = 0.05
# CFG.enrichment_min_size = 15
# CFG.enrichment_max_size = 500
# CFG.enrichment_permutations = 1000

# ── CCI spatial cell-cell interaction analysis (uncomment to enable) ──
CFG.run_cci = True
# CFG.cci_method = "liana"                  # 'liana' (LIANA+ rank_aggregate)
# CFG.cci_lr_database = "consensus"         # 'consensus' | 'cellphonedb' | 'cellchat' | 'celltalkdb' | 'ramilowski'
# CFG.cci_permutations = 1000               # permutation test iterations
# CFG.cci_n_top_interactions = 50           # top N interactions for heatmap
# CFG.cci_spatial_method = "liana_spatial"  # spatial method (reserves 'commot')
# CFG.cci_spatial_distance = 0.0            # 0 = use existing spatial_connectivities
# CFG.cci_lr_cache_dir = ""                 # LIANA cache dir; empty = auto (~/.cache/liana)
# # CFG.cci_multi_condition = False         # future: multi-condition CCI

# ── Image processing ──
CFG.crop_image = True
CFG.img_rescale = 1.0

# ── AI settings (uncomment to enable) ──
# CFG.ai.enabled = True
# CFG.ai.api_base = 'https://api.deepseek.com/v1'
# CFG.ai.model = 'deepseek-v4-pro'
# CFG.ai.api_key = os.environ.get('LLM_API_KEY', '')
# CFG.ai.max_tokens = 32768
# CFG.ai.temperature = 0.1
# CFG.ai.thinking_enabled = True
# CFG.ai.reasoning_effort = 'high'
# CFG.ai.ai_annotation = True  # Enable AI cell type annotation
# CFG.ai.ai_interpretation = True
# CFG.ai.ai_cache_responses = True

# ── Execution ──
CFG.n_jobs = 0   # auto-detect
CFG.random_seed = 42
# CFG.force_csr = True
# CFG.use_float32 = False
# CFG.limit_blas_threads = True
# CFG.scanpy_verbosity = 2
# CFG.h5ad_compression = 'gzip'
# CFG.h5ad_tempdir = '/tmp/Fuxi'
