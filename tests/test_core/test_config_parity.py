"""Verify new Pydantic Config model defaults and YAML round-trip.

This test suite ensures the Pydantic v2 Config model produces identical
defaults to the old dataclass and that YAML serialization/deserialization
works correctly.
"""

import pathlib

import pytest
import yaml

from core.config.schema import Config
from core.utils._config import resolve_config

# ═══════════════════════════════════════════════════════════════════════
#  Defaults parity
# ═══════════════════════════════════════════════════════════════════════


class TestConfigDefaults:
    """Verify every sub-model has the correct default values."""

    def test_top_level_defaults(self) -> None:
        """Root-level string/primitive defaults match the dataclass era."""
        cfg = Config()
        assert cfg.modality == "rna"
        assert cfg.tissue == "unknown"
        assert cfg.species == "human"
        assert cfg.tissue_maturity == "unknown"
        assert cfg.expression_type == "raw_counts"
        assert cfg.data_format == "10X_mtx"
        assert cfg.h5ad_compression == "gzip"
        assert cfg.h5ad_tempdir == "/tmp/Fuxi"
        assert cfg.cleanup_intermediates is False
        assert cfg.perf_monitoring is True

    def test_data_input_defaults(self) -> None:
        """DataInputConfig defaults (skip path fields resolved by model_post_init)."""
        cfg = Config()
        assert cfg.data_input.mtx_prefix == ""
        assert cfg.data_input.matrix_file == ""
        assert cfg.data_input.barcodes_file == ""
        assert cfg.data_input.features_file == ""
        assert cfg.data_input.csv_sep is None
        assert cfg.data_input.csv_decimal == "."
        assert cfg.data_input.gene_symbol_column == ""
        assert cfg.data_input.input_h5ad == ""
        assert cfg.data_input.backed == ""
        assert cfg.data_input.h5_file_pattern == "*filtered_feature_bc_matrix.h5"
        assert cfg.data_input.fragment_file == ""
        # mtx_dir, h5_dir, data_dir are auto-resolved by model_post_init
        # from FUXI_DATA_ROOT env var — skip in pure-defaults test

    def test_qc_defaults(self) -> None:
        """QC settings match the dataclass era."""
        cfg = Config()
        assert cfg.qc.min_genes == 500
        assert cfg.qc.max_genes == 7500
        assert cfg.qc.max_pct_mito == 20.0
        assert cfg.qc.mt_gene_pattern == "MT-"
        assert cfg.qc.min_genes_per_umi == 0.7
        assert cfg.qc.min_cells_per_gene == 3
        assert cfg.qc.use_adaptive_thresholds is False
        assert cfg.qc.mad_n_mads == 3.0
        assert cfg.qc.ncount_max_mad == 5.0
        assert cfg.qc.min_mad_upper_genes == 4000
        assert cfg.qc.min_mad_upper_genes_nuclei == 3000
        assert cfg.qc.is_nuclei is False
        assert cfg.qc.max_pct_mito_nuclei == 5.0

    def test_scrublet_defaults(self) -> None:
        """Scrublet defaults."""
        cfg = Config()
        assert cfg.scrublet.run is True
        assert cfg.scrublet.expected_doublet_rate is None
        assert cfg.scrublet.batch_key == "sample"
        assert cfg.scrublet.min_counts == 2
        assert cfg.scrublet.min_cells == 3
        assert cfg.scrublet.min_gene_var_pctl == 85
        assert cfg.scrublet.n_prin_comps == 30

    def test_normalization_defaults(self) -> None:
        """Normalization defaults."""
        cfg = Config()
        assert cfg.normalization.normalize_target_sum == 1e4
        assert cfg.normalization.use_regress_out is True
        assert cfg.normalization.score_cell_cycle is False
        assert cfg.normalization.detect_sex is True

    def test_hvg_defaults(self) -> None:
        """HVG defaults."""
        cfg = Config()
        assert cfg.hvg.n_top_genes == 4000
        assert cfg.hvg.flavor == "seurat_v3"
        assert cfg.hvg.batch_key == "sample"

    def test_pca_defaults(self) -> None:
        """PCA defaults."""
        cfg = Config()
        assert cfg.pca.n_pcs_full == 100
        assert cfg.pca.n_pcs_use == 50

    def test_integration_defaults(self) -> None:
        """Integration defaults."""
        cfg = Config()
        assert cfg.integration.method == "harmony"
        assert cfg.integration.batch_key == "sample"
        assert cfg.integration.max_iter == 20

    def test_clustering_defaults(self) -> None:
        """Clustering defaults."""
        cfg = Config()
        assert cfg.clustering.n_neighbors == 30
        assert cfg.clustering.param_grid_resolutions == [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
        assert cfg.clustering.param_grid_n_neighbors == [15, 20, 30]
        assert cfg.clustering.param_grid_resolutions == [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
        assert cfg.clustering.leiden_flavor == "igraph"
        assert cfg.clustering.best_resolution == 1.0
        assert cfg.clustering.best_n_neighbors == 0
        assert cfg.clustering.cluster_selection_method == "multi_metric"
        assert cfg.clustering.umap_selection_method == "convex_hull"
        assert cfg.clustering.umap_min_dist == 0.5
        assert cfg.clustering.umap_spread == 1.0

    def test_marker_defaults(self) -> None:
        """Marker settings defaults."""
        cfg = Config()
        assert cfg.marker.marker_dict == {}
        assert cfg.marker.subcluster_types == []
        assert cfg.marker.subcluster_resolution == 0.4
        assert cfg.marker.min_cells_subcluster == 50
        assert cfg.marker.validation_n_top_genes == 15
        assert cfg.marker.validation_min_overlap == 0.5
        assert cfg.marker.validation_pass_rate_min == 0.1
        assert cfg.marker.quality_gate_min_pass_rate == 0.10

    def test_de_defaults(self) -> None:
        """DE settings defaults."""
        cfg = Config()
        assert cfg.de.method == "pseudobulk"
        assert cfg.de.n_genes == 50
        assert cfg.de.pval_cutoff == 0.05
        assert cfg.de.logfc_cutoff == 0.25
        assert cfg.de.stage_pairwise is True

    def test_trajectory_defaults(self) -> None:
        """Trajectory defaults."""
        cfg = Config()
        assert cfg.trajectory.root_cell_types == []
        assert cfg.trajectory.root_markers == []
        assert cfg.trajectory.n_diffmap_comps == 15
        assert cfg.trajectory.n_branchings == 2
        assert cfg.trajectory.pseudotime_genes == []
        assert cfg.trajectory.pseudotime_n_branch_de == 10
        assert cfg.trajectory.pseudotime_n_correlated == 10
        assert cfg.trajectory.pseudotime_cor_pval == 0.05

    def test_enrichment_defaults(self) -> None:
        """Enrichment defaults."""
        cfg = Config()
        assert cfg.enrichment.run is True
        assert cfg.enrichment.method == "both"
        assert cfg.enrichment.organism == "human"
        assert cfg.enrichment.n_top_genes == 200
        assert cfg.enrichment.pval_cutoff == 0.05
        assert cfg.enrichment.min_size == 10
        assert cfg.enrichment.max_size == 500
        assert cfg.enrichment.permutations == 1000

    def test_grn_defaults(self) -> None:
        """GRN defaults."""
        cfg = Config()
        assert cfg.grn.run is True
        assert cfg.grn.method == "decoupler"
        assert cfg.grn.species == "human"
        assert cfg.grn.n_top_regulons == 50
        assert cfg.grn.min_regulon_size == 5

    def test_cci_defaults(self) -> None:
        """CCI defaults."""
        cfg = Config()
        assert cfg.cci.run is True
        assert cfg.cci.method == "liana"
        assert cfg.cci.permutations == 100
        assert cfg.cci.n_top_interactions == 50

    def test_execution_defaults(self) -> None:
        """Execution defaults."""
        cfg = Config()
        assert cfg.execution.n_jobs == 0
        assert cfg.execution.limit_blas_threads is True
        assert cfg.execution.random_seed == 42
        assert cfg.execution.scanpy_verbosity == 2
        assert cfg.execution.force_csr is True
        assert cfg.execution.use_float32 is True

    def test_atac_defaults(self) -> None:
        """ATAC-specific defaults."""
        cfg = Config()
        assert cfg.atac.genome == "hg38"
        assert cfg.atac.min_fragments == 1000
        assert cfg.atac.max_fragments == 50000
        assert cfg.atac.min_tsse == 7.0
        assert cfg.atac.n_features == 50000
        assert cfg.atac.use_macs3 is True

    def test_spatial_defaults(self) -> None:
        """Spatial defaults."""
        cfg = Config()
        assert cfg.spatial.platform == "visium"
        assert cfg.spatial.neighbors_n == 6
        assert cfg.spatial.run_autocorr is True
        assert cfg.spatial.svg_n_top == 2000

    def test_ai_defaults(self) -> None:
        """AI config defaults."""
        cfg = Config()
        assert cfg.ai.enabled is False
        assert cfg.ai.model == "deepseek-v4-flash"
        assert cfg.ai.max_tokens == 4096
        assert cfg.ai.temperature == 0.1
        assert cfg.ai.thinking_enabled is True
        assert cfg.ai.annotation is True
        assert cfg.ai.subcluster is True
        assert cfg.ai.interpretation is True
        assert cfg.ai.cache_responses is True
        assert cfg.ai.unconstrained_annotation is False

    def test_annotation_multi_peak_defaults(self) -> None:
        """AnnotationSettings 多峰歧义降级阈值默认值 (D1)."""
        cfg = Config()
        assert cfg.annotation.multi_peak_min_types == 3
        assert cfg.annotation.multi_peak_score_floor == 0.9

    def test_annotation_canonical_pct_floor_default(self) -> None:
        """AnnotationSettings canonical 表达兜底 pct 下限默认值 (D3)."""
        cfg = Config()
        assert cfg.annotation.canonical_pct_floor == 0.05

    def test_annotation_kadp_defaults(self) -> None:
        """AnnotationSettings KADP developmental-potency defaults (plan todo 5)."""
        cfg = Config()
        assert cfg.annotation.kadp_enabled is False
        assert cfg.annotation.kadp_ratio_threshold == 2.0
        assert cfg.annotation.kadp_abs_threshold == 0.6
        assert cfg.annotation.kadp_gap_threshold == 0.1
        assert cfg.annotation.use_gap_criterion is False

    def test_annotation_metc_defaults(self) -> None:
        """AnnotationSettings METC multi-source voting defaults (plan todo 10)."""
        cfg = Config()
        assert cfg.annotation.metc_enabled is False
        assert cfg.annotation.metc_min_sources == 3
        assert cfg.annotation.metc_min_distinct_transition == 3

    def test_downsample_defaults(self) -> None:
        """Downsample defaults."""
        cfg = Config()
        assert cfg.downsample.target is None
        assert cfg.downsample.strategy == "stratified"
        assert cfg.downsample.max_per_sample is None
        assert cfg.downsample.random_seed == 42

    def test_all_sub_models_exist(self) -> None:
        """All 20 expected sub-models are present on Config()."""
        cfg = Config()
        expected = [
            "data_input",
            "sample_meta",
            "qc",
            "scrublet",
            "normalization",
            "hvg",
            "pca",
            "integration",
            "clustering",
            "marker",
            "de",
            "trajectory",
            "enrichment",
            "grn",
            "cci",
            "downsample",
            "spatial",
            "atac",
            "execution",
            "ai",
        ]
        for name in expected:
            assert hasattr(cfg, name), f"Missing sub-model: {name}"
        assert len(expected) == 20, "Expected exactly 20 sub-models"


# ═══════════════════════════════════════════════════════════════════════
#  YAML round-trip
# ═══════════════════════════════════════════════════════════════════════


class TestConfigYamlRoundTrip:
    """Config.model_validate() with a YAML dict produces correct values."""

    SAMPLE_YAML = """
modality: rna
data_format: 10X_h5
tissue: retina
species: human

data_input:
  h5_file_pattern: "*filtered_feature_bc_matrix.h5"

qc:
  min_genes: 200
  max_genes: 6000
  max_pct_mito: 25.0

scrublet:
  run: true

hvg:
  n_top_genes: 3000

pca:
  n_pcs_use: 40

integration:
  method: harmony

clustering:
  param_grid_resolutions: [0.5, 1.0, 2.0]
  cluster_selection_method: "multi_metric"

de:
  method: "wilcoxon"
  n_genes: 100

trajectory:
  root_cell_types: ["RGC"]
  n_diffmap_comps: 10

enrichment:
  run: false

grn:
  run: false

cci:
  run: false

execution:
  n_jobs: 8
  random_seed: 123
"""

    def test_model_validate_from_yaml_dict(self) -> None:
        """Config.model_validate() works with a parsed YAML dict."""
        data = yaml.safe_load(self.SAMPLE_YAML)
        cfg = Config.model_validate(data)

        # Top-level
        assert cfg.modality == "rna"
        assert cfg.tissue == "retina"
        assert cfg.species == "human"
        assert cfg.data_format == "10X_h5"

        # Sub-model fields
        assert cfg.qc.min_genes == 200
        assert cfg.qc.max_genes == 6000
        assert cfg.qc.max_pct_mito == 25.0
        assert cfg.hvg.n_top_genes == 3000
        assert cfg.pca.n_pcs_use == 40
        assert cfg.clustering.param_grid_resolutions == [0.5, 1.0, 2.0]
        assert cfg.de.n_genes == 100
        assert cfg.trajectory.root_cell_types == ["RGC"]
        assert cfg.trajectory.n_diffmap_comps == 10
        assert cfg.execution.n_jobs == 8
        assert cfg.execution.random_seed == 123

        # Overridden non-default values
        assert cfg.enrichment.run is False
        assert cfg.grn.run is False
        assert cfg.cci.run is False

        # Unset values still have defaults
        assert cfg.qc.mt_gene_pattern == "MT-"
        assert cfg.de.method == "wilcoxon"
        assert cfg.pca.n_pcs_full == 100

    def test_model_validate_with_defaults(self) -> None:
        """Parsing empty YAML yields the same defaults as Config()."""
        data = {}
        cfg = Config.model_validate(data)

        # Spot-check a few defaults
        assert cfg.modality == "rna"
        assert cfg.hvg.n_top_genes == 4000
        assert cfg.qc.min_genes == 500
        assert cfg.de.method == "pseudobulk"
        assert cfg.clustering.cluster_selection_method == "multi_metric"
        assert cfg.execution.n_jobs == 0

    def test_model_validate_rejects_extra_keys(self) -> None:
        """Extra fields in YAML raise ValidationError (extra='forbid')."""
        data = {"nonexistent_field": 123}
        with pytest.raises(Exception, match="extra"):
            Config.model_validate(data)

    def test_model_validate_rejects_extra_subkeys(self) -> None:
        """Extra sub-keys in a sub-model raise ValidationError."""
        data = {"qc": {"nonexistent_sub_field": 42}}
        with pytest.raises(Exception, match="extra"):
            Config.model_validate(data)

    def test_yaml_template_round_trip(self, tmp_path: pathlib.Path) -> None:
        """Loading the 10X H5 starter config YAML works."""
        from core.config.scaffold import render_template_text
        from core.preprocess.config_specs import materialized_specs

        spec = next(s for s in materialized_specs() if s.key == "10X_h5")
        cfg_path = tmp_path / "config_10X_h5.yaml"
        cfg_path.write_text(render_template_text(spec), encoding="utf-8")

        with open(cfg_path) as f:
            data = yaml.safe_load(f)

        cfg = Config.model_validate(data)

        # Spot-check the template values match
        assert cfg.data_format == "10X_h5"
        assert cfg.data_input.h5_file_pattern == "*filtered_feature_bc_matrix.h5"
        assert cfg.qc.min_genes == 200
        assert cfg.hvg.n_top_genes == 4000
        assert cfg.integration.method == "harmony"
        assert cfg.clustering.cluster_selection_method == "multi_metric"
        assert cfg.de.method == "pseudobulk"
        assert cfg.execution.n_jobs == 0


# ═══════════════════════════════════════════════════════════════════════
#  Model serialization
# ═══════════════════════════════════════════════════════════════════════


class TestConfigSerialization:
    """Config can be serialized to dict/JSON and re-loaded."""

    def test_model_dump_round_trip(self) -> None:
        """Config.model_dump() → model_validate recovers the same values."""
        cfg = Config.model_validate(
            {
                "modality": "atac",
                "tissue": "brain",
                "hvg": {"n_top_genes": 2000},
                "qc": {"min_genes": 300},
            }
        )
        dumped = cfg.model_dump()
        restored = Config.model_validate(dumped)

        assert restored.modality == "atac"
        assert restored.tissue == "brain"
        assert restored.hvg.n_top_genes == 2000
        assert restored.qc.min_genes == 300

        # Unset defaults survive the round-trip
        assert restored.de.method == "pseudobulk"
        assert restored.execution.random_seed == 42

    def test_model_dump_json_round_trip(self) -> None:
        """JSON serialization preserves nested defaults."""
        cfg = Config.model_validate({"clustering": {"umap_min_dist": 0.5, "umap_spread": 1.5}})
        json_str = cfg.model_dump_json()
        restored = Config.model_validate_json(json_str)

        assert restored.clustering.umap_min_dist == 0.5
        assert restored.clustering.umap_spread == 1.5


# ═══════════════════════════════════════════════════════════════════════
#  Edge cases
# ═══════════════════════════════════════════════════════════════════════


class TestConfigEdgeCases:
    """Edge-case behavior for the Config model."""

    def test_nested_dict_override(self) -> None:
        """Nested dict overrides only specified sub-fields, keeping defaults."""
        cfg = Config.model_validate({"hvg": {"n_top_genes": 2000}})
        assert cfg.hvg.n_top_genes == 2000
        assert cfg.hvg.flavor == "seurat_v3"  # default preserved
        assert cfg.hvg.batch_key == "sample"  # default preserved

    def test_partial_sub_model(self) -> None:
        """Partial sub-model dicts fill unspecified fields with defaults."""
        cfg = Config.model_validate({"qc": {"min_genes": 200}})
        assert cfg.qc.min_genes == 200
        assert cfg.qc.max_genes == 7500  # default preserved
        assert cfg.qc.max_pct_mito == 20.0  # default preserved

    def test_empty_string_tissue_kb(self) -> None:
        """Empty tissue_kb is accepted."""
        cfg = Config(tissue_kb="")
        assert cfg.tissue_kb == ""


class TestAIConfigSubclusterKbConstrained:
    """ai.subcluster_kb_constrained — Step 06 KB-constrained subtype naming flag."""

    def test_absent_defaults_to_true(self) -> None:
        """Key absent from YAML → default True."""
        cfg = Config()
        assert cfg.ai.subcluster_kb_constrained is True

    def test_yaml_false_round_trips(self, tmp_path) -> None:
        """YAML ai.subcluster_kb_constrained: false parses and round-trips via resolve_config."""
        cfg_path = tmp_path / "config_t4.yaml"
        cfg_path.write_text("modality: rna\nai:\n  subcluster_kb_constrained: false\n")
        cfg = resolve_config(str(cfg_path))
        assert cfg.ai.subcluster_kb_constrained is False

    def test_yaml_true_round_trips(self, tmp_path) -> None:
        """Explicit YAML true value survives round-trip."""
        cfg_path = tmp_path / "config_t4.yaml"
        cfg_path.write_text("modality: rna\nai:\n  subcluster_kb_constrained: true\n")
        cfg = resolve_config(str(cfg_path))
        assert cfg.ai.subcluster_kb_constrained is True
