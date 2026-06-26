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
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
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

# ── HVG ──
CFG.n_top_genes = 4000
CFG.hvg_flavor = 'seurat_v3'

# ── Spatial graph ──
CFG.spatial_neighbors_n = 6
CFG.spatial_neighbors_radius = 0.0   # Use n_neighbors mode

# ── PCA ──
CFG.n_pcs_full = 100
CFG.n_pcs_use = 50

# ── Clustering ──
CFG.n_neighbors = 30
CFG.leiden_resolutions = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
CFG.best_resolution = 1.0
CFG.leiden_flavor = 'igraph'
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
# Alternative: auto-detect root via marker genes
# CFG.root_markers = ['SOX2', 'PAX6', 'NES']
CFG.n_diffmap_comps = 15
CFG.n_branchings = 2

# ── Enrichment ──
CFG.run_enrichment = True
CFG.enrichment_method = 'both'   # 'ora' | 'prerank' | 'both'
CFG.enrichment_gene_sets = [
    'GO_Biological_Process_2023',
    'KEGG_2021_Human',
]
CFG.enrichment_organism = 'human'

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
