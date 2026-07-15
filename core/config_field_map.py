#!/usr/bin/env python3
"""
config_field_map.py — Field-to-topic mapping for Config Pydantic migration

Defines three exports consumed by:
  - Todo 1.1: Pydantic sub-model definitions (FIELD_MAP → topic boundaries)
  - Todo 3.1: .py → .yaml config migration script (FIELD_MAP → nested paths)
  - Todo 4.2: Consumer rename CFG.xxx → CFG.topic.xxx (FIELD_MAP → target paths)

Usage:
    from core.config_field_map import FIELD_MAP, DUPLICATE_FIELDS, TOPIC_NAMES

    new_path = FIELD_MAP["n_top_genes"]        # → "hvg.n_top_genes"
    dup_info = DUPLICATE_FIELDS["csv_sep"]     # → {"sources": {...}, "winner": "Config", ...}
    model_cls_name = TOPIC_NAMES["hvg"]        # → "HVGSettings"
"""

from typing import Dict, Any


# ═══════════════════════════════════════════════════════════════════════
# TOPIC_NAMES — Topic key → Pydantic sub-model class name
# ═══════════════════════════════════════════════════════════════════════
# Each topic key is the dotted prefix used in FIELD_MAP values.
# Each value is the Pydantic BaseModel class name to be defined in
# core/config.py Todo 1.1.
TOPIC_NAMES: Dict[str, str] = {
    # ── RNA data input (formats, file paths) ──
    "data_input": "DataInputConfig",
    # ── Sample / stage metadata mapping ──
    "sample_meta": "SampleMetaConfig",
    # ── Quality control thresholds ──
    "qc": "QCSettings",
    # ── Doublet detection (Scrublet) ──
    "scrublet": "ScrubletSettings",
    # ── Normalization, cell-cycle, sex detection ──
    "normalization": "NormalizationSettings",
    # ── Highly variable gene selection ──
    "hvg": "HVGSettings",
    # ── Principal component analysis ──
    "pca": "PCASettings",
    # ── Harmony batch correction ──
    "harmony": "HarmonySettings",
    # ── Clustering, UMAP, parameter sweep ──
    "clustering": "ClusteringSettings",
    # ── Cell-type marker / annotation ──
    "marker": "MarkerSettings",
    # ── Differential expression ──
    "de": "DESettings",
    # ── Pseudotime / trajectory ──
    "trajectory": "TrajectorySettings",
    # ── Gene-set enrichment ──
    "enrichment": "EnrichmentSettings",
    # ── Gene regulatory network ──
    "grn": "GRNSettings",
    # ── Cell-cell interaction ──
    "cci": "CCISettings",
    # ── Downsampling / subset filtering ──
    "downsample": "DownsampleSettings",
    # ── Spatial transcriptomics ──
    "spatial": "SpatialConfig",
    # ── ATAC-specific ──
    "atac": "ATACConfig",
    # ── Execution environment ──
    "execution": "ExecutionConfig",
    # ── AI / LLM configuration ──
    "ai": "AIConfig",
}


# ═══════════════════════════════════════════════════════════════════════
# FIELD_MAP — Every current flat field name → new dotted topic path
# ═══════════════════════════════════════════════════════════════════════
#
# Naming convention:
#   If a field name has a `<topic>_` prefix that matches its target topic,
#   the prefix is dropped in the new path (e.g. `de_method` → `de.method`).
#   Otherwise the field name is preserved as-is.
#
# Top-level fields (no topic prefix) map to themselves — they stay on the
# Config root in the new model.
FIELD_MAP: Dict[str, str] = {
    # ═══════════════════════════════════════════════════════════════════
    # Root-level fields (stay on Config root)
    # ═══════════════════════════════════════════════════════════════════
    "modality": "modality",
    "tissue": "tissue",
    "species": "species",
    "tissue_maturity": "tissue_maturity",
    "expression_type": "expression_type",
    "data_format": "data_format",
    "data_dir": "data_dir",
    "results_dir": "results_dir",
    "h5ad_dir": "h5ad_dir",
    "figure_dir": "figure_dir",
    "table_dir": "table_dir",
    "log_dir": "log_dir",
    "project_dir": "project_dir",
    "h5ad_compression": "h5ad_compression",
    "h5ad_tempdir": "h5ad_tempdir",
    "cleanup_intermediates": "cleanup_intermediates",
    "perf_monitoring": "perf_monitoring",
    "rna_h5ad": "rna_h5ad",
    "rna_ref": "rna_ref",
    "rna_marker_top_n": "rna_marker_top_n",
    "rna_marker_pval_threshold": "rna_marker_pval_threshold",
    "rna_marker_logfc_min": "rna_marker_logfc_min",
    "tissue_kb": "tissue_kb",
    "tissue_ontology": "tissue_ontology",
    "target_class": "target_class",
    "target_order": "target_order",

    # ═══════════════════════════════════════════════════════════════════
    # data_input — DataInputConfig
    # ═══════════════════════════════════════════════════════════════════
    "mtx_prefix": "data_input.mtx_prefix",
    "mtx_dir": "data_input.mtx_dir",
    "matrix_file": "data_input.matrix_file",
    "barcodes_file": "data_input.barcodes_file",
    "features_file": "data_input.features_file",
    "csv_sep": "data_input.csv_sep",
    "csv_decimal": "data_input.csv_decimal",
    "gene_symbol_column": "data_input.gene_symbol_column",
    "input_h5ad": "data_input.input_h5ad",
    "backed": "data_input.backed",
    "h5_file_pattern": "data_input.h5_file_pattern",
    "h5_dir": "data_input.h5_dir",
    "fragment_file": "data_input.fragment_file",
    "file_pattern": "data_input.file_pattern",
    "separator": "data_input.separator",

    # ═══════════════════════════════════════════════════════════════════
    # sample_meta — SampleMetaConfig
    # ═══════════════════════════════════════════════════════════════════
    "sample_map": "sample_meta.sample_map",
    "stage_map": "sample_meta.stage_map",
    "stage_order": "sample_meta.stage_order",
    "meta_columns": "sample_meta.meta_columns",
    "barcode_parse_regex": "sample_meta.barcode_parse_regex",
    "barcode_parse_groups": "sample_meta.barcode_parse_groups",

    # ═══════════════════════════════════════════════════════════════════
    # qc — QCSettings
    # ═══════════════════════════════════════════════════════════════════
    "min_genes": "qc.min_genes",
    "max_genes": "qc.max_genes",
    "max_pct_mito": "qc.max_pct_mito",
    "mt_gene_pattern": "qc.mt_gene_pattern",
    "mt_gene_list": "qc.mt_gene_list",
    "min_genes_per_umi": "qc.min_genes_per_umi",
    "min_cells_per_gene": "qc.min_cells_per_gene",
    "use_adaptive_thresholds": "qc.use_adaptive_thresholds",
    "mad_n_mads": "qc.mad_n_mads",
    "qc_ncount_max_mad": "qc.ncount_max_mad",
    "min_mad_upper_genes": "qc.min_mad_upper_genes",
    "min_mad_upper_genes_nuclei": "qc.min_mad_upper_genes_nuclei",
    "is_nuclei": "qc.is_nuclei",
    "max_pct_mito_nuclei": "qc.max_pct_mito_nuclei",

    # ═══════════════════════════════════════════════════════════════════
    # scrublet — ScrubletSettings
    # ═══════════════════════════════════════════════════════════════════
    "run_scrublet": "scrublet.run",
    "scrublet_expected_doublet_rate": "scrublet.expected_doublet_rate",
    "scrublet_batch_key": "scrublet.batch_key",
    "scrublet_min_counts": "scrublet.min_counts",
    "scrublet_min_cells": "scrublet.min_cells",
    "scrublet_min_gene_var_pctl": "scrublet.min_gene_var_pctl",
    "scrublet_n_prin_comps": "scrublet.n_prin_comps",

    # ═══════════════════════════════════════════════════════════════════
    # normalization — NormalizationSettings
    # ═══════════════════════════════════════════════════════════════════
    "normalize_target_sum": "normalization.normalize_target_sum",
    "use_regress_out": "normalization.use_regress_out",
    "score_cell_cycle": "normalization.score_cell_cycle",
    "regress_out_genes": "normalization.regress_out_genes",
    "detect_sex": "normalization.detect_sex",

    # ═══════════════════════════════════════════════════════════════════
    # hvg — HVGSettings
    # ═══════════════════════════════════════════════════════════════════
    "n_top_genes": "hvg.n_top_genes",
    "hvg_flavor": "hvg.flavor",
    "hvg_batch_key": "hvg.batch_key",
    "hvg_forced_genes": "hvg.forced_genes",

    # ═══════════════════════════════════════════════════════════════════
    # pca — PCASettings
    # ═══════════════════════════════════════════════════════════════════
    "n_pcs_full": "pca.n_pcs_full",
    "n_pcs_use": "pca.n_pcs_use",

    # ═══════════════════════════════════════════════════════════════════
    # harmony — HarmonySettings
    # ═══════════════════════════════════════════════════════════════════
    "use_harmony": "harmony.use_harmony",
    "harmony_batch_key": "harmony.batch_key",
    "harmony_max_iter": "harmony.max_iter",
    "harmony_collinearity_guard": "harmony.collinearity_guard",

    # ═══════════════════════════════════════════════════════════════════
    # clustering — ClusteringSettings
    # ═══════════════════════════════════════════════════════════════════
    "n_neighbors": "clustering.n_neighbors",
    "leiden_resolutions": "clustering.leiden_resolutions",
    "param_grid_n_neighbors": "clustering.param_grid_n_neighbors",
    "param_grid_resolutions": "clustering.param_grid_resolutions",
    "leiden_flavor": "clustering.leiden_flavor",
    "best_resolution": "clustering.best_resolution",
    "best_n_neighbors": "clustering.best_n_neighbors",
    "cluster_selection_method": "clustering.cluster_selection_method",
    "multi_metric_weights": "clustering.multi_metric_weights",
    "multi_metric_n_stability_seeds": "clustering.multi_metric_n_stability_seeds",
    "multi_metric_adaptive_resolution": "clustering.multi_metric_adaptive_resolution",
    "multi_metric_coverage_ratio_threshold": "clustering.multi_metric_coverage_ratio_threshold",
    "multi_metric_coherence_dominance": "clustering.multi_metric_coherence_dominance",
    "multi_metric_granularity_cv_threshold": "clustering.multi_metric_granularity_cv_threshold",
    "multi_metric_granularity_min_clusters": "clustering.multi_metric_granularity_min_clusters",
    "multi_metric_de_gate_threshold": "clustering.multi_metric_de_gate_threshold",
    "umap_selection_method": "clustering.umap_selection_method",
    "param_grid_min_dist": "clustering.param_grid_min_dist",
    "param_grid_spread": "clustering.param_grid_spread",
    "umap_min_dist": "clustering.umap_min_dist",
    "umap_spread": "clustering.umap_spread",

    # ═══════════════════════════════════════════════════════════════════
    # marker — MarkerSettings
    # ═══════════════════════════════════════════════════════════════════
    "marker_dict": "marker.marker_dict",
    "subcluster_types": "marker.subcluster_types",
    "subcluster_resolution": "marker.subcluster_resolution",
    "min_cells_subcluster": "marker.min_cells_subcluster",
    "expert_rule_strictness": "marker.expert_rule_strictness",
    "expert_rule_top_n": "marker.expert_rule_top_n",
    "expert_rule_pval_cutoff": "marker.expert_rule_pval_cutoff",
    "marker_validation_n_top_genes": "marker.validation_n_top_genes",
    "marker_validation_min_overlap": "marker.validation_min_overlap",
    "marker_validation_marginal_threshold": "marker.validation_marginal_threshold",
    "marker_validation_pass_rate_min": "marker.validation_pass_rate_min",
    "step10_groupby": "marker.step10_groupby",
    "quality_gate_min_pass_rate": "marker.quality_gate_min_pass_rate",

    # ═══════════════════════════════════════════════════════════════════
    # de — DESettings
    # ═══════════════════════════════════════════════════════════════════
    "de_method": "de.method",
    "de_n_genes": "de.n_genes",
    "de_pval_cutoff": "de.pval_cutoff",
    "de_logfc_cutoff": "de.logfc_cutoff",
    "de_stage_pairwise": "de.stage_pairwise",
    "de_auto_switch_on_low_quality": "de.auto_switch_on_low_quality",

    # ═══════════════════════════════════════════════════════════════════
    # trajectory — TrajectorySettings
    # ═══════════════════════════════════════════════════════════════════
    "root_cell_types": "trajectory.root_cell_types",
    "root_markers": "trajectory.root_markers",
    "n_diffmap_comps": "trajectory.n_diffmap_comps",
    "n_branchings": "trajectory.n_branchings",
    "pseudotime_genes": "trajectory.pseudotime_genes",
    "pseudotime_n_branch_de": "trajectory.pseudotime_n_branch_de",
    "pseudotime_n_correlated": "trajectory.pseudotime_n_correlated",
    "pseudotime_cor_pval": "trajectory.pseudotime_cor_pval",

    # ═══════════════════════════════════════════════════════════════════
    # enrichment — EnrichmentSettings
    # ═══════════════════════════════════════════════════════════════════
    "run_enrichment": "enrichment.run",
    "enrichment_method": "enrichment.method",
    "enrichment_gene_sets": "enrichment.gene_sets",
    "enrichment_organism": "enrichment.organism",
    "enrichment_n_top_genes": "enrichment.n_top_genes",
    "enrichment_pval_cutoff": "enrichment.pval_cutoff",
    "enrichment_min_size": "enrichment.min_size",
    "enrichment_max_size": "enrichment.max_size",
    "enrichment_permutations": "enrichment.permutations",
    "enrichment_tissue_mode": "enrichment.tissue_mode",
    "enrichment_tissue_pathways_whitelist": "enrichment.tissue_pathways_whitelist",
    "enrichment_tissue_pathways_blacklist": "enrichment.tissue_pathways_blacklist",
    "enrichment_redundancy_cluster": "enrichment.redundancy_cluster",
    "enrichment_redundancy_threshold": "enrichment.redundancy_threshold",
    "enrichment_use_kb_relevance": "enrichment.use_kb_relevance",
    "enrichment_gene_sets_tissue": "enrichment.gene_sets_tissue",
    "enrichment_background_restrict": "enrichment.background_restrict",
    "peak_gene_distance": "enrichment.peak_gene_distance",
    "gene_annotation_bed": "enrichment.gene_annotation_bed",

    # ═══════════════════════════════════════════════════════════════════
    # grn — GRNSettings
    # ═══════════════════════════════════════════════════════════════════
    "run_grn": "grn.run",
    "grn_method": "grn.method",
    "grn_species": "grn.species",
    "grn_n_top_regulons": "grn.n_top_regulons",
    "grn_min_regulon_size": "grn.min_regulon_size",
    "grn_confidence_levels": "grn.confidence_levels",
    "grn_tissue_mode": "grn.tissue_mode",
    "grn_use_kb_relevance": "grn.use_kb_relevance",
    "grn_export_filtered": "grn.export_filtered",

    # ═══════════════════════════════════════════════════════════════════
    # cci — CCISettings
    # ═══════════════════════════════════════════════════════════════════
    "run_cci": "cci.run",
    "cci_method": "cci.method",
    "cci_lr_database": "cci.lr_database",
    "cci_permutations": "cci.permutations",
    "cci_n_top_interactions": "cci.n_top_interactions",
    "cci_spatial_method": "cci.spatial_method",
    "cci_spatial_distance": "cci.spatial_distance",
    "cci_lr_cache_dir": "cci.lr_cache_dir",
    "cci_adjacency": "cci.adjacency",
    "cci_tissue": "cci.tissue",
    "cci_adjacency_file": "cci.adjacency_file",
    "cci_adjacency_types": "cci.adjacency_types",

    # ═══════════════════════════════════════════════════════════════════
    # downsample — DownsampleSettings
    # ═══════════════════════════════════════════════════════════════════
    "downsample_target": "downsample.target",
    "downsample_strategy": "downsample.strategy",
    "downsample_max_per_sample": "downsample.max_per_sample",
    "downsample_random_seed": "downsample.random_seed",
    "sample_keep": "downsample.sample_keep",
    "obs_filter": "downsample.obs_filter",
    "subset_suffix": "downsample.subset_suffix",

    # ═══════════════════════════════════════════════════════════════════
    # spatial — SpatialConfig
    # ═══════════════════════════════════════════════════════════════════
    "spatial_platform": "spatial.platform",
    "library_id": "spatial.library_id",
    "img_path": "spatial.img_path",
    "spot_diameter": "spatial.spot_diameter",
    "crop_image": "spatial.crop_image",
    "img_rescale": "spatial.img_rescale",
    "spatial_neighbors_n": "spatial.neighbors_n",
    "spatial_neighbors_radius": "spatial.neighbors_radius",
    "run_spatial_autocorr": "spatial.run_autocorr",
    "moran_percentile": "spatial.moran_percentile",
    "svg_n_top": "spatial.svg_n_top",

    # ═══════════════════════════════════════════════════════════════════
    # atac — ATACConfig
    # ═══════════════════════════════════════════════════════════════════
    "genome": "atac.genome",
    "chrom_sizes": "atac.chrom_sizes",
    "blacklist_bed": "atac.blacklist_bed",
    "tss_bed": "atac.tss_bed",
    "min_fragments": "atac.min_fragments",
    "max_fragments": "atac.max_fragments",
    "min_tsse": "atac.min_tsse",
    "max_blacklist_ratio": "atac.max_blacklist_ratio",
    "min_peak_region_fragments": "atac.min_peak_region_fragments",
    "peak_qval": "atac.peak_qval",
    "peak_width": "atac.peak_width",
    "use_macs3": "atac.use_macs3",
    "n_features": "atac.n_features",
    "n_spectral": "atac.n_spectral",
    "marker_peaks_log2fc": "atac.marker_peaks_log2fc",
    "marker_peaks_fdr": "atac.marker_peaks_fdr",
    "motif_db": "atac.motif_db",
    "terminal_cell_types": "atac.terminal_cell_types",
    "max_cells": "atac.max_cells",

    # ═══════════════════════════════════════════════════════════════════
    # execution — ExecutionConfig
    # ═══════════════════════════════════════════════════════════════════
    "n_jobs": "execution.n_jobs",
    "limit_blas_threads": "execution.limit_blas_threads",
    "random_seed": "execution.random_seed",
    "scanpy_verbosity": "execution.scanpy_verbosity",
    "force_csr": "execution.force_csr",
    "use_float32": "execution.use_float32",

    # ═══════════════════════════════════════════════════════════════════
    # ai — AIConfig
    # ═══════════════════════════════════════════════════════════════════
    # These are already nested under Config.ai in the current dataclass;
    # they remain under Config.ai in the Pydantic model.
    "enabled": "ai.enabled",
    "api_base": "ai.api_base",
    "model": "ai.model",
    "api_key": "ai.api_key",
    "max_tokens": "ai.max_tokens",
    "temperature": "ai.temperature",
    "thinking_enabled": "ai.thinking_enabled",
    "reasoning_effort": "ai.reasoning_effort",
    "timeout": "ai.timeout",
    "ai_qc_review": "ai.qc_review",
    "ai_param_suggest": "ai.param_suggest",
    "ai_annotation": "ai.annotation",
    "ai_subcluster": "ai.subcluster",
    "ai_deg_design": "ai.deg_design",
    "ai_interpretation": "ai.interpretation",
    "ai_cache_responses": "ai.cache_responses",
    "unconstrained_annotation": "ai.unconstrained_annotation",
}


# ═══════════════════════════════════════════════════════════════════════
# DUPLICATE_FIELDS — Fields that exist in BOTH Config and a nested config
# ═══════════════════════════════════════════════════════════════════════
#
# These are fields defined on the top-level Config dataclass AND also on
# RNAConfig, ATACConfig, or SpatialConfig.  In the current system Config's
# value shadows the nested config's value (since Python MRO checks the
# class's own __dict__ first).  In the new Pydantic model each duplicated
# field lives in exactly ONE topic sub-model.
#
# Fields with identical defaults on both sides are marked accordingly.
# For fields with DIFFERENT defaults, the rationale for the winner is
# documented.
DUPLICATE_FIELDS: Dict[str, Dict[str, Any]] = {
    # ── DataInputConfig duplicates (Config ↔ RNAConfig) ──
    "mtx_prefix": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.mtx_prefix",
        "rationale": "Same default in both sources; Config-level takes precedence as user-visible value.",
    },
    "mtx_dir": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.mtx_dir",
        "rationale": "Same default in both sources.",
    },
    "matrix_file": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.matrix_file",
        "rationale": "Same default in both sources.",
    },
    "barcodes_file": {
        "sources": {"Config": "", "RNAConfig": "", "ATACConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.barcodes_file",
        "rationale": "Same default in all three sources; Config is the user-facing layer.",
    },
    "features_file": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.features_file",
        "rationale": "Same default in both sources.",
    },
    "csv_sep": {
        "sources": {"Config": None, "RNAConfig": ","},
        "winner": "Config",
        "winner_default": None,
        "target_field": "data_input.csv_sep",
        "rationale": (
            "Config default None makes csv_sep optional — not all projects use CSV input. "
            "RNAConfig default ',' would force comma-delimited assumption. "
            "None allows DataInputConfig to distinguish 'not set' from 'explicit comma'."
        ),
    },
    "csv_decimal": {
        "sources": {"Config": ".", "RNAConfig": "."},
        "winner": "Config",
        "winner_default": ".",
        "target_field": "data_input.csv_decimal",
        "rationale": "Same default in both sources.",
    },
    "gene_symbol_column": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.gene_symbol_column",
        "rationale": "Same default in both sources.",
    },
    "input_h5ad": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.input_h5ad",
        "rationale": "Same default in both sources.",
    },
    "backed": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.backed",
        "rationale": "Same default in both sources.",
    },
    "h5_file_pattern": {
        "sources": {"Config": "*filtered_feature_bc_matrix.h5", "RNAConfig": "*filtered_feature_bc_matrix.h5"},
        "winner": "Config",
        "winner_default": "*filtered_feature_bc_matrix.h5",
        "target_field": "data_input.h5_file_pattern",
        "rationale": "Same default in both sources.",
    },
    "h5_dir": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.h5_dir",
        "rationale": "Same default in both sources.",
    },

    # ── SampleMetaConfig duplicates (Config ↔ RNAConfig) ──
    "sample_map": {
        "sources": {"Config": {}, "RNAConfig": {}},
        "winner": "Config",
        "winner_default": {},
        "target_field": "sample_meta.sample_map",
        "rationale": "Same default in both sources.",
    },
    "stage_map": {
        "sources": {"Config": {}, "RNAConfig": {}},
        "winner": "Config",
        "winner_default": {},
        "target_field": "sample_meta.stage_map",
        "rationale": "Same default in both sources.",
    },
    "stage_order": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "sample_meta.stage_order",
        "rationale": "Same default in both sources.",
    },
    "meta_columns": {
        "sources": {"Config": {}, "RNAConfig": {}},
        "winner": "Config",
        "winner_default": {},
        "target_field": "sample_meta.meta_columns",
        "rationale": "Same default in both sources.",
    },
    "barcode_parse_regex": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "sample_meta.barcode_parse_regex",
        "rationale": "Same default in both sources.",
    },
    "barcode_parse_groups": {
        "sources": {"Config": {}, "RNAConfig": {}},
        "winner": "Config",
        "winner_default": {},
        "target_field": "sample_meta.barcode_parse_groups",
        "rationale": "Same default in both sources.",
    },

    # ── QCSettings duplicates (Config ↔ RNAConfig) ──
    "min_genes": {
        "sources": {"Config": 500, "RNAConfig": 500},
        "winner": "Config",
        "winner_default": 500,
        "target_field": "qc.min_genes",
        "rationale": "Same default in both sources.",
    },
    "max_genes": {
        "sources": {"Config": 7500, "RNAConfig": 7500},
        "winner": "Config",
        "winner_default": 7500,
        "target_field": "qc.max_genes",
        "rationale": "Same default in both sources.",
    },
    "max_pct_mito": {
        "sources": {"Config": 20.0, "RNAConfig": 20.0},
        "winner": "Config",
        "winner_default": 20.0,
        "target_field": "qc.max_pct_mito",
        "rationale": "Same default in both sources.",
    },
    "mt_gene_pattern": {
        "sources": {"Config": "MT-", "RNAConfig": "MT-"},
        "winner": "Config",
        "winner_default": "MT-",
        "target_field": "qc.mt_gene_pattern",
        "rationale": "Same default in both sources.",
    },
    "mt_gene_list": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "qc.mt_gene_list",
        "rationale": "Same default in both sources.",
    },
    "min_genes_per_umi": {
        "sources": {"Config": 0.7, "RNAConfig": 0.7},
        "winner": "Config",
        "winner_default": 0.7,
        "target_field": "qc.min_genes_per_umi",
        "rationale": "Same default in both sources.",
    },
    "min_cells_per_gene": {
        "sources": {"Config": 3, "RNAConfig": 3},
        "winner": "Config",
        "winner_default": 3,
        "target_field": "qc.min_cells_per_gene",
        "rationale": "Same default in both sources.",
    },
    "use_adaptive_thresholds": {
        "sources": {"Config": False, "RNAConfig": False},
        "winner": "Config",
        "winner_default": False,
        "target_field": "qc.use_adaptive_thresholds",
        "rationale": "Same default in both sources.",
    },
    "mad_n_mads": {
        "sources": {"Config": 3.0, "RNAConfig": 3.0},
        "winner": "Config",
        "winner_default": 3.0,
        "target_field": "qc.mad_n_mads",
        "rationale": "Same default in both sources.",
    },
    "qc_ncount_max_mad": {
        "sources": {"Config": 5.0, "RNAConfig": 5.0},
        "winner": "Config",
        "winner_default": 5.0,
        "target_field": "qc.ncount_max_mad",
        "rationale": "Same default in both sources.",
    },

    # ── ScrubletSettings duplicates (Config ↔ RNAConfig) ──
    "run_scrublet": {
        "sources": {"Config": True, "RNAConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "scrublet.run",
        "rationale": "Same default in both sources.",
    },
    "scrublet_expected_doublet_rate": {
        "sources": {"Config": None, "RNAConfig": None},
        "winner": "Config",
        "winner_default": None,
        "target_field": "scrublet.expected_doublet_rate",
        "rationale": "Same default in both sources.",
    },
    "scrublet_batch_key": {
        "sources": {"Config": "sample", "RNAConfig": "sample"},
        "winner": "Config",
        "winner_default": "sample",
        "target_field": "scrublet.batch_key",
        "rationale": "Same default in both sources.",
    },
    "scrublet_min_counts": {
        "sources": {"Config": 2, "RNAConfig": 2},
        "winner": "Config",
        "winner_default": 2,
        "target_field": "scrublet.min_counts",
        "rationale": "Same default in both sources.",
    },
    "scrublet_min_cells": {
        "sources": {"Config": 3, "RNAConfig": 3},
        "winner": "Config",
        "winner_default": 3,
        "target_field": "scrublet.min_cells",
        "rationale": "Same default in both sources.",
    },
    "scrublet_min_gene_var_pctl": {
        "sources": {"Config": 85, "RNAConfig": 85},
        "winner": "Config",
        "winner_default": 85,
        "target_field": "scrublet.min_gene_var_pctl",
        "rationale": "Same default in both sources.",
    },
    "scrublet_n_prin_comps": {
        "sources": {"Config": 30, "RNAConfig": 30},
        "winner": "Config",
        "winner_default": 30,
        "target_field": "scrublet.n_prin_comps",
        "rationale": "Same default in both sources.",
    },

    # ── NormalizationSettings duplicates (Config ↔ RNAConfig) ──
    "normalize_target_sum": {
        "sources": {"Config": 10000.0, "RNAConfig": 10000.0},
        "winner": "Config",
        "winner_default": 10000.0,
        "target_field": "normalization.normalize_target_sum",
        "rationale": "Same default in both sources.",
    },
    "use_regress_out": {
        "sources": {"Config": True, "RNAConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "normalization.use_regress_out",
        "rationale": "Same default in both sources.",
    },
    "score_cell_cycle": {
        "sources": {"Config": False, "RNAConfig": False},
        "winner": "Config",
        "winner_default": False,
        "target_field": "normalization.score_cell_cycle",
        "rationale": "Same default in both sources.",
    },

    # ── HVGSettings duplicates (Config ↔ RNAConfig) ──
    "n_top_genes": {
        "sources": {"Config": 4000, "RNAConfig": 4000},
        "winner": "Config",
        "winner_default": 4000,
        "target_field": "hvg.n_top_genes",
        "rationale": "Same default in both sources.",
    },
    "hvg_flavor": {
        "sources": {"Config": "seurat_v3", "RNAConfig": "seurat_v3"},
        "winner": "Config",
        "winner_default": "seurat_v3",
        "target_field": "hvg.flavor",
        "rationale": "Same default in both sources.",
    },
    "hvg_batch_key": {
        "sources": {"Config": "sample", "RNAConfig": "sample"},
        "winner": "Config",
        "winner_default": "sample",
        "target_field": "hvg.batch_key",
        "rationale": "Same default in both sources.",
    },

    # ── PCASettings duplicates (Config ↔ RNAConfig) ──
    "n_pcs_full": {
        "sources": {"Config": 100, "RNAConfig": 100},
        "winner": "Config",
        "winner_default": 100,
        "target_field": "pca.n_pcs_full",
        "rationale": "Same default in both sources.",
    },
    "n_pcs_use": {
        "sources": {"Config": 50, "RNAConfig": 50},
        "winner": "Config",
        "winner_default": 50,
        "target_field": "pca.n_pcs_use",
        "rationale": "Same default in both sources.",
    },

    # ── HarmonySettings duplicates (Config ↔ RNAConfig) ──
    "use_harmony": {
        "sources": {"Config": True, "RNAConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "harmony.use_harmony",
        "rationale": "Same default in both sources.",
    },
    "harmony_batch_key": {
        "sources": {"Config": "sample", "RNAConfig": "sample"},
        "winner": "Config",
        "winner_default": "sample",
        "target_field": "harmony.batch_key",
        "rationale": "Same default in both sources.",
    },
    "harmony_max_iter": {
        "sources": {"Config": 20, "RNAConfig": 20},
        "winner": "Config",
        "winner_default": 20,
        "target_field": "harmony.max_iter",
        "rationale": "Same default in both sources.",
    },

    # ── ClusteringSettings duplicates (Config ↔ RNAConfig) ──
    "n_neighbors": {
        "sources": {"Config": 30, "RNAConfig": 30},
        "winner": "Config",
        "winner_default": 30,
        "target_field": "clustering.n_neighbors",
        "rationale": "Same default in both sources.",
    },
    "leiden_resolutions": {
        "sources": {"Config": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0], "RNAConfig": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]},
        "winner": "Config",
        "winner_default": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0],
        "target_field": "clustering.leiden_resolutions",
        "rationale": "Same default in both sources.",
    },
    "param_grid_n_neighbors": {
        "sources": {"Config": [15, 20, 30], "RNAConfig": [15, 20, 30]},
        "winner": "Config",
        "winner_default": [15, 20, 30],
        "target_field": "clustering.param_grid_n_neighbors",
        "rationale": "Same default in both sources.",
    },
    "param_grid_resolutions": {
        "sources": {"Config": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0], "RNAConfig": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]},
        "winner": "Config",
        "winner_default": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0],
        "target_field": "clustering.param_grid_resolutions",
        "rationale": "Same default in both sources.",
    },
    "leiden_flavor": {
        "sources": {"Config": "igraph", "RNAConfig": "igraph"},
        "winner": "Config",
        "winner_default": "igraph",
        "target_field": "clustering.leiden_flavor",
        "rationale": "Same default in both sources.",
    },
    "best_resolution": {
        "sources": {"Config": 1.0, "RNAConfig": 1.0},
        "winner": "Config",
        "winner_default": 1.0,
        "target_field": "clustering.best_resolution",
        "rationale": "Same default in both sources.",
    },
    "best_n_neighbors": {
        "sources": {"Config": 0, "RNAConfig": 0},
        "winner": "Config",
        "winner_default": 0,
        "target_field": "clustering.best_n_neighbors",
        "rationale": "Same default in both sources.",
    },
    "cluster_selection_method": {
        "sources": {"Config": "multi_metric", "RNAConfig": "pareto_elbow"},
        "winner": "Config",
        "winner_default": "multi_metric",
        "target_field": "clustering.cluster_selection_method",
        "rationale": (
            "Config default 'multi_metric' is the newer, more comprehensive method. "
            "RNAConfig default 'pareto_elbow' is the older RNA-specific default. "
            "Config wins because Config is the user-facing layer and multi_metric is "
            "the current recommended approach."
        ),
    },
    "umap_selection_method": {
        "sources": {"Config": "convex_hull", "RNAConfig": "convex_hull"},
        "winner": "Config",
        "winner_default": "convex_hull",
        "target_field": "clustering.umap_selection_method",
        "rationale": "Same default in both sources.",
    },
    "param_grid_min_dist": {
        "sources": {"Config": [0.1, 0.3, 0.5], "RNAConfig": [0.1, 0.3, 0.5]},
        "winner": "Config",
        "winner_default": [0.1, 0.3, 0.5],
        "target_field": "clustering.param_grid_min_dist",
        "rationale": "Same default in both sources.",
    },
    "param_grid_spread": {
        "sources": {"Config": [1.0], "RNAConfig": [1.0]},
        "winner": "Config",
        "winner_default": [1.0],
        "target_field": "clustering.param_grid_spread",
        "rationale": "Same default in both sources.",
    },
    "umap_min_dist": {
        "sources": {"Config": 0.3, "RNAConfig": 0.3},
        "winner": "Config",
        "winner_default": 0.3,
        "target_field": "clustering.umap_min_dist",
        "rationale": "Same default in both sources.",
    },
    "umap_spread": {
        "sources": {"Config": 1.0, "RNAConfig": 1.0},
        "winner": "Config",
        "winner_default": 1.0,
        "target_field": "clustering.umap_spread",
        "rationale": "Same default in both sources.",
    },

    # ── MarkerSettings duplicates (Config ↔ RNAConfig) ──
    "marker_dict": {
        "sources": {"Config": {}, "RNAConfig": {}},
        "winner": "Config",
        "winner_default": {},
        "target_field": "marker.marker_dict",
        "rationale": "Same default in both sources.",
    },
    "subcluster_types": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "marker.subcluster_types",
        "rationale": "Same default in both sources.",
    },
    "subcluster_resolution": {
        "sources": {"Config": 0.4, "RNAConfig": 0.4},
        "winner": "Config",
        "winner_default": 0.4,
        "target_field": "marker.subcluster_resolution",
        "rationale": "Same default in both sources.",
    },
    "min_cells_subcluster": {
        "sources": {"Config": 50, "RNAConfig": 50},
        "winner": "Config",
        "winner_default": 50,
        "target_field": "marker.min_cells_subcluster",
        "rationale": "Same default in both sources.",
    },
    "expert_rule_strictness": {
        "sources": {"Config": "default", "RNAConfig": "default"},
        "winner": "Config",
        "winner_default": "default",
        "target_field": "marker.expert_rule_strictness",
        "rationale": "Same default in both sources.",
    },
    "expert_rule_top_n": {
        "sources": {"Config": 0, "RNAConfig": 0},
        "winner": "Config",
        "winner_default": 0,
        "target_field": "marker.expert_rule_top_n",
        "rationale": "Same default in both sources.",
    },
    "expert_rule_pval_cutoff": {
        "sources": {"Config": 0.0, "RNAConfig": 0.0},
        "winner": "Config",
        "winner_default": 0.0,
        "target_field": "marker.expert_rule_pval_cutoff",
        "rationale": "Same default in both sources.",
    },
    "marker_validation_n_top_genes": {
        "sources": {"Config": 15, "RNAConfig": 15},
        "winner": "Config",
        "winner_default": 15,
        "target_field": "marker.validation_n_top_genes",
        "rationale": "Same default in both sources.",
    },
    "marker_validation_min_overlap": {
        "sources": {"Config": 0.5, "RNAConfig": 0.5},
        "winner": "Config",
        "winner_default": 0.5,
        "target_field": "marker.validation_min_overlap",
        "rationale": "Same default in both sources.",
    },
    "marker_validation_marginal_threshold": {
        "sources": {"Config": 0.25, "RNAConfig": 0.25},
        "winner": "Config",
        "winner_default": 0.25,
        "target_field": "marker.validation_marginal_threshold",
        "rationale": "Same default in both sources.",
    },
    "marker_validation_pass_rate_min": {
        "sources": {"Config": 0.1, "RNAConfig": 0.1},
        "winner": "Config",
        "winner_default": 0.1,
        "target_field": "marker.validation_pass_rate_min",
        "rationale": "Same default in both sources.",
    },
    "step10_groupby": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "marker.step10_groupby",
        "rationale": "Same default in both sources.",
    },

    # ── DESettings duplicates (Config ↔ RNAConfig) ──
    "de_method": {
        "sources": {"Config": "wilcoxon", "RNAConfig": "wilcoxon"},
        "winner": "Config",
        "winner_default": "wilcoxon",
        "target_field": "de.method",
        "rationale": "Same default in both sources.",
    },
    "de_n_genes": {
        "sources": {"Config": 50, "RNAConfig": 50},
        "winner": "Config",
        "winner_default": 50,
        "target_field": "de.n_genes",
        "rationale": "Same default in both sources.",
    },
    "de_pval_cutoff": {
        "sources": {"Config": 0.05, "RNAConfig": 0.05},
        "winner": "Config",
        "winner_default": 0.05,
        "target_field": "de.pval_cutoff",
        "rationale": "Same default in both sources.",
    },
    "de_logfc_cutoff": {
        "sources": {"Config": 0.25, "RNAConfig": 0.25},
        "winner": "Config",
        "winner_default": 0.25,
        "target_field": "de.logfc_cutoff",
        "rationale": "Same default in both sources.",
    },
    "de_stage_pairwise": {
        "sources": {"Config": True, "RNAConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "de.stage_pairwise",
        "rationale": "Same default in both sources.",
    },
    "de_auto_switch_on_low_quality": {
        "sources": {"Config": False, "RNAConfig": False},
        "winner": "Config",
        "winner_default": False,
        "target_field": "de.auto_switch_on_low_quality",
        "rationale": "Same default in both sources.",
    },

    # ── TrajectorySettings duplicates (Config ↔ RNAConfig) ──
    "root_cell_types": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "trajectory.root_cell_types",
        "rationale": "Same default in both sources.",
    },
    "root_markers": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "trajectory.root_markers",
        "rationale": "Same default in both sources.",
    },
    "n_diffmap_comps": {
        "sources": {"Config": 15, "RNAConfig": 15},
        "winner": "Config",
        "winner_default": 15,
        "target_field": "trajectory.n_diffmap_comps",
        "rationale": "Same default in both sources.",
    },
    "n_branchings": {
        "sources": {"Config": 2, "RNAConfig": 2},
        "winner": "Config",
        "winner_default": 2,
        "target_field": "trajectory.n_branchings",
        "rationale": "Same default in both sources.",
    },

    # ── GRNSettings duplicates (Config ↔ RNAConfig) ──
    "run_grn": {
        "sources": {"Config": True, "RNAConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "grn.run",
        "rationale": "Same default in both sources.",
    },
    "grn_method": {
        "sources": {"Config": "decoupler", "RNAConfig": "decoupler"},
        "winner": "Config",
        "winner_default": "decoupler",
        "target_field": "grn.method",
        "rationale": "Same default in both sources.",
    },
    "grn_species": {
        "sources": {"Config": "human", "RNAConfig": "human"},
        "winner": "Config",
        "winner_default": "human",
        "target_field": "grn.species",
        "rationale": "Same default in both sources.",
    },
    "grn_n_top_regulons": {
        "sources": {"Config": 50, "RNAConfig": 50},
        "winner": "Config",
        "winner_default": 50,
        "target_field": "grn.n_top_regulons",
        "rationale": "Same default in both sources.",
    },
    "grn_min_regulon_size": {
        "sources": {"Config": 5, "RNAConfig": 5},
        "winner": "Config",
        "winner_default": 5,
        "target_field": "grn.min_regulon_size",
        "rationale": "Same default in both sources.",
    },
    "grn_confidence_levels": {
        "sources": {"Config": ["A", "B", "C"], "RNAConfig": ["A", "B", "C"]},
        "winner": "Config",
        "winner_default": ["A", "B", "C"],
        "target_field": "grn.confidence_levels",
        "rationale": "Same default in both sources.",
    },
    "grn_tissue_mode": {
        "sources": {"Config": "off", "RNAConfig": "off"},
        "winner": "Config",
        "winner_default": "off",
        "target_field": "grn.tissue_mode",
        "rationale": "Same default in both sources.",
    },
    "grn_use_kb_relevance": {
        "sources": {"Config": False, "RNAConfig": False},
        "winner": "Config",
        "winner_default": False,
        "target_field": "grn.use_kb_relevance",
        "rationale": "Same default in both sources.",
    },
    "grn_export_filtered": {
        "sources": {"Config": False, "RNAConfig": False},
        "winner": "Config",
        "winner_default": False,
        "target_field": "grn.export_filtered",
        "rationale": "Same default in both sources.",
    },

    # ── CCISettings duplicates (Config ↔ RNAConfig) ──
    "run_cci": {
        "sources": {"Config": True, "RNAConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "cci.run",
        "rationale": "Same default in both sources.",
    },
    "cci_method": {
        "sources": {"Config": "liana", "RNAConfig": "liana"},
        "winner": "Config",
        "winner_default": "liana",
        "target_field": "cci.method",
        "rationale": "Same default in both sources.",
    },
    "cci_lr_database": {
        "sources": {"Config": "consensus", "RNAConfig": "consensus"},
        "winner": "Config",
        "winner_default": "consensus",
        "target_field": "cci.lr_database",
        "rationale": "Same default in both sources.",
    },
    "cci_permutations": {
        "sources": {"Config": 100, "RNAConfig": 1000},
        "winner": "Config",
        "winner_default": 100,
        "target_field": "cci.permutations",
        "rationale": (
            "Config default 100 is more conservative (faster) for broad CCI screening. "
            "RNAConfig default 1000 provides higher statistical power at 10x compute cost. "
            "Config wins because 100 is the safer default that users can increase if needed; "
            "most CCI analyses in this pipeline use 100 permutations."
        ),
    },
    "cci_n_top_interactions": {
        "sources": {"Config": 50, "RNAConfig": 50},
        "winner": "Config",
        "winner_default": 50,
        "target_field": "cci.n_top_interactions",
        "rationale": "Same default in both sources.",
    },
    "cci_spatial_method": {
        "sources": {"Config": "liana_spatial", "RNAConfig": "liana_spatial"},
        "winner": "Config",
        "winner_default": "liana_spatial",
        "target_field": "cci.spatial_method",
        "rationale": "Same default in both sources.",
    },
    "cci_spatial_distance": {
        "sources": {"Config": 0.0, "RNAConfig": 0.0},
        "winner": "Config",
        "winner_default": 0.0,
        "target_field": "cci.spatial_distance",
        "rationale": "Same default in both sources.",
    },
    "cci_lr_cache_dir": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "cci.lr_cache_dir",
        "rationale": "Same default in both sources.",
    },

    # ── DownsampleSettings duplicates (Config ↔ RNAConfig) ──
    "downsample_target": {
        "sources": {"Config": None, "RNAConfig": None},
        "winner": "Config",
        "winner_default": None,
        "target_field": "downsample.target",
        "rationale": "Same default in both sources.",
    },
    "downsample_strategy": {
        "sources": {"Config": "stratified", "RNAConfig": "stratified"},
        "winner": "Config",
        "winner_default": "stratified",
        "target_field": "downsample.strategy",
        "rationale": "Same default in both sources.",
    },
    "downsample_max_per_sample": {
        "sources": {"Config": None, "RNAConfig": None},
        "winner": "Config",
        "winner_default": None,
        "target_field": "downsample.max_per_sample",
        "rationale": "Same default in both sources.",
    },
    "downsample_random_seed": {
        "sources": {"Config": 42, "RNAConfig": 42},
        "winner": "Config",
        "winner_default": 42,
        "target_field": "downsample.random_seed",
        "rationale": "Same default in both sources.",
    },
    "sample_keep": {
        "sources": {"Config": [], "RNAConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "downsample.sample_keep",
        "rationale": "Same default in both sources.",
    },
    "obs_filter": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "downsample.obs_filter",
        "rationale": "Same default in both sources.",
    },
    "subset_suffix": {
        "sources": {"Config": "", "RNAConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "downsample.subset_suffix",
        "rationale": "Same default in both sources.",
    },

    # ── ATAC duplicates (Config ↔ ATACConfig) ──
    "fragment_file": {
        "sources": {"Config": "", "ATACConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.fragment_file",
        "rationale": "Same default in both sources. Maps to data_input because fragment_file is an input format path.",
    },
    "barcodes_file": {
        "sources": {"Config": "", "ATACConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "data_input.barcodes_file",
        "rationale": "Same default in both sources. Already covered above (also duplicated with RNAConfig).",
    },
    "genome": {
        "sources": {"Config": "hg38", "ATACConfig": "hg38"},
        "winner": "Config",
        "winner_default": "hg38",
        "target_field": "atac.genome",
        "rationale": "Same default in both sources.",
    },
    "chrom_sizes": {
        "sources": {"Config": "", "ATACConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "atac.chrom_sizes",
        "rationale": "Same default in both sources.",
    },
    "blacklist_bed": {
        "sources": {"Config": "", "ATACConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "atac.blacklist_bed",
        "rationale": "Same default in both sources.",
    },
    "tss_bed": {
        "sources": {"Config": "", "ATACConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "atac.tss_bed",
        "rationale": "Same default in both sources.",
    },
    "min_fragments": {
        "sources": {"Config": 1000, "ATACConfig": 1000},
        "winner": "Config",
        "winner_default": 1000,
        "target_field": "atac.min_fragments",
        "rationale": "Same default in both sources.",
    },
    "max_fragments": {
        "sources": {"Config": 50000, "ATACConfig": 50000},
        "winner": "Config",
        "winner_default": 50000,
        "target_field": "atac.max_fragments",
        "rationale": "Same default in both sources.",
    },
    "min_tsse": {
        "sources": {"Config": 7.0, "ATACConfig": 7.0},
        "winner": "Config",
        "winner_default": 7.0,
        "target_field": "atac.min_tsse",
        "rationale": "Same default in both sources.",
    },
    "max_blacklist_ratio": {
        "sources": {"Config": 0.05, "ATACConfig": 0.05},
        "winner": "Config",
        "winner_default": 0.05,
        "target_field": "atac.max_blacklist_ratio",
        "rationale": "Same default in both sources.",
    },
    "min_peak_region_fragments": {
        "sources": {"Config": 300, "ATACConfig": 300},
        "winner": "Config",
        "winner_default": 300,
        "target_field": "atac.min_peak_region_fragments",
        "rationale": "Same default in both sources.",
    },
    "peak_qval": {
        "sources": {"Config": 0.05, "ATACConfig": 0.05},
        "winner": "Config",
        "winner_default": 0.05,
        "target_field": "atac.peak_qval",
        "rationale": "Same default in both sources.",
    },
    "peak_width": {
        "sources": {"Config": 500, "ATACConfig": 500},
        "winner": "Config",
        "winner_default": 500,
        "target_field": "atac.peak_width",
        "rationale": "Same default in both sources.",
    },
    "use_macs3": {
        "sources": {"Config": True, "ATACConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "atac.use_macs3",
        "rationale": "Same default in both sources.",
    },
    "n_features": {
        "sources": {"Config": 50000, "ATACConfig": 50000},
        "winner": "Config",
        "winner_default": 50000,
        "target_field": "atac.n_features",
        "rationale": "Same default in both sources.",
    },
    "n_spectral": {
        "sources": {"Config": 30, "ATACConfig": 30},
        "winner": "Config",
        "winner_default": 30,
        "target_field": "atac.n_spectral",
        "rationale": "Same default in both sources.",
    },
    "marker_peaks_log2fc": {
        "sources": {"Config": 0.5, "ATACConfig": 0.5},
        "winner": "Config",
        "winner_default": 0.5,
        "target_field": "atac.marker_peaks_log2fc",
        "rationale": "Same default in both sources.",
    },
    "marker_peaks_fdr": {
        "sources": {"Config": 0.05, "ATACConfig": 0.05},
        "winner": "Config",
        "winner_default": 0.05,
        "target_field": "atac.marker_peaks_fdr",
        "rationale": "Same default in both sources.",
    },
    "motif_db": {
        "sources": {"Config": "JASPAR2024", "ATACConfig": "JASPAR2024"},
        "winner": "Config",
        "winner_default": "JASPAR2024",
        "target_field": "atac.motif_db",
        "rationale": "Same default in both sources.",
    },
    "terminal_cell_types": {
        "sources": {"Config": [], "ATACConfig": []},
        "winner": "Config",
        "winner_default": [],
        "target_field": "atac.terminal_cell_types",
        "rationale": "Same default in both sources.",
    },
    "max_cells": {
        "sources": {"Config": None, "ATACConfig": None},
        "winner": "Config",
        "winner_default": None,
        "target_field": "atac.max_cells",
        "rationale": "Same default in both sources.",
    },

    # ── Spatial duplicates (Config ↔ SpatialConfig) ──
    "spatial_platform": {
        "sources": {"Config": "visium", "SpatialConfig": "visium"},
        "winner": "Config",
        "winner_default": "visium",
        "target_field": "spatial.platform",
        "rationale": "Same default in both sources.",
    },
    "library_id": {
        "sources": {"Config": "", "SpatialConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "spatial.library_id",
        "rationale": "Same default in both sources.",
    },
    "img_path": {
        "sources": {"Config": "", "SpatialConfig": ""},
        "winner": "Config",
        "winner_default": "",
        "target_field": "spatial.img_path",
        "rationale": "Same default in both sources.",
    },
    "spot_diameter": {
        "sources": {"Config": 0.0, "SpatialConfig": 0.0},
        "winner": "Config",
        "winner_default": 0.0,
        "target_field": "spatial.spot_diameter",
        "rationale": "Same default in both sources.",
    },
    "crop_image": {
        "sources": {"Config": True, "SpatialConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "spatial.crop_image",
        "rationale": "Same default in both sources.",
    },
    "img_rescale": {
        "sources": {"Config": 1.0, "SpatialConfig": 1.0},
        "winner": "Config",
        "winner_default": 1.0,
        "target_field": "spatial.img_rescale",
        "rationale": "Same default in both sources.",
    },
    "spatial_neighbors_n": {
        "sources": {"Config": 6, "SpatialConfig": 6},
        "winner": "Config",
        "winner_default": 6,
        "target_field": "spatial.neighbors_n",
        "rationale": "Same default in both sources.",
    },
    "spatial_neighbors_radius": {
        "sources": {"Config": 0.0, "SpatialConfig": 0.0},
        "winner": "Config",
        "winner_default": 0.0,
        "target_field": "spatial.neighbors_radius",
        "rationale": "Same default in both sources.",
    },
    "run_spatial_autocorr": {
        "sources": {"Config": True, "SpatialConfig": True},
        "winner": "Config",
        "winner_default": True,
        "target_field": "spatial.run_autocorr",
        "rationale": "Same default in both sources.",
    },
    "moran_percentile": {
        "sources": {"Config": 90, "SpatialConfig": 90},
        "winner": "Config",
        "winner_default": 90,
        "target_field": "spatial.moran_percentile",
        "rationale": "Same default in both sources.",
    },
    "svg_n_top": {
        "sources": {"Config": 2000, "SpatialConfig": 2000},
        "winner": "Config",
        "winner_default": 2000,
        "target_field": "spatial.svg_n_top",
        "rationale": "Same default in both sources.",
    },

    # ── Execution duplicates (Config ↔ RNAConfig) ──
    # n_jobs exists on Config (line 727) and RNAConfig has NO n_jobs.
    # However n_jobs is listed in the plan's discrepancy list — checking Config vs ATACConfig?
    # Actually, n_jobs is Config-only. But the plan mentions it as "documented discrepancy"
    # because RNAConfig/ATACConfig might reference it during __getattr__. Not an actual duplicate.
    "n_jobs": {
        "sources": {"Config": 0},
        "winner": "Config",
        "winner_default": 0,
        "target_field": "execution.n_jobs",
        "rationale": "Config-only field. 0 = auto-detect os.cpu_count(). Not duplicated in any nested config.",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# RNAConfig-only fields (not duplicated on Config — accessed via __getattr__)
# ═══════════════════════════════════════════════════════════════════════
# These fields exist ONLY on RNAConfig, not on Config. In the current
# system they are accessed through Config.__getattr__ delegation.
# In the new Pydantic model they live in the appropriate topic sub-model.
RNAONLY_FIELDS: Dict[str, str] = {
    "regress_out_genes": "normalization.regress_out_genes",
    "detect_sex": "normalization.detect_sex",
    "pseudotime_genes": "trajectory.pseudotime_genes",
    "pseudotime_n_branch_de": "trajectory.pseudotime_n_branch_de",
    "pseudotime_n_correlated": "trajectory.pseudotime_n_correlated",
    "pseudotime_cor_pval": "trajectory.pseudotime_cor_pval",
}

# ═══════════════════════════════════════════════════════════════════════
# Completeness check helpers
# ═══════════════════════════════════════════════════════════════════════

# All field names known to exist in the current dataclass system.
# Used by the schema-parity test (Todo 4.5.2) to verify completeness.
ALL_KNOWN_FIELDS: list[str] = sorted(FIELD_MAP.keys())

# Topic membership: maps each field name → topic key
FIELD_TOPIC: Dict[str, str] = {}
for _name, _path in FIELD_MAP.items():
    _parts = _path.split(".")
    FIELD_TOPIC[_name] = _parts[0] if len(_parts) > 1 else "__root__"
