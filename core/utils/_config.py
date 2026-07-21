"""Config loading — resolve_config, species validation, dataset.yaml helpers."""

import logging
import os
from typing import Optional

import yaml

from core.config.dataset import load_dataset
from core.config.schema import Config

# ── Dataset.yaml helpers ───────────────────────────────────────────


def _find_dataset_yaml(cfg) -> Optional[str]:
    """Search project_dir → data_dir for dataset.yaml."""
    for base in [cfg.project_dir, cfg.data_dir]:
        if base:
            path = os.path.join(base, "dataset.yaml")
            if os.path.exists(path):
                return path
    return None


def resolve_config(config_path: Optional[str] = None) -> Config:
    """
    Load a YAML config file and return a resolved Config instance.

    Replaces the old importlib-based .py loading with Pydantic v2 YAML loading.
    """
    if config_path is None:
        # Default: look for config.yaml in parent directory
        config_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            "config.yaml",
        )

    config_path = os.path.abspath(config_path)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    if config_path.endswith(".py"):
        raise ValueError(
            f"Python configs are no longer supported: {config_path}\n"
            f"Migrate to YAML format using: python core/migrate_configs.py --gse <GSE_ID>"
        )
    if not config_path.endswith(".yaml"):
        raise ValueError(f"Config file must be .yaml format: {config_path}")

    # ── Load YAML and validate ──
    with open(config_path) as f:
        data = yaml.safe_load(f)

    # Auto-detect project_dir from config file location BEFORE validation
    # (must precede model_validate so model_post_init sees correct project_dir)
    if not data.get("project_dir"):
        data["project_dir"] = os.path.dirname(config_path)

    cfg = Config.model_validate(data)

    # ── Resolve n_jobs ──
    if cfg.execution.n_jobs == 0:
        cfg.execution.n_jobs = os.cpu_count() or 1

    # ── Auto-fill is_nuclei from dataset.yaml ──
    # Check if is_nuclei is explicitly set in YAML config first
    if not cfg.qc.is_nuclei:
        _yaml = _find_dataset_yaml(cfg)
        if _yaml:
            try:
                _ds = load_dataset(_yaml)
                if getattr(_ds, "assay_type", None) == "snRNAseq":
                    cfg.qc.is_nuclei = True
                    print(f"[Config] Auto-set is_nuclei=True from {_yaml}")
            except Exception:
                pass  # Graceful: no crash on invalid yaml

    # ── Species sanity check ──
    _validate_species(cfg)

    # ── Create output directories ──
    for d in [cfg.results_dir, cfg.h5ad_dir, cfg.figure_dir, cfg.table_dir, cfg.log_dir]:
        os.makedirs(d, exist_ok=True)

    return cfg


# Species pipeline-keys known to rna/ortholog.py and rna/utils/marker_scoring.py.
# Keep in sync when adding new species.
_KNOWN_SPECIES_KEYS = frozenset(
    {
        "human",
        "mouse",
        "macaque",
        "cynomolgus",
        "marmoset",
        "zebrafish",
        "chicken",
        "lamprey",
        "frog",
        "pig",
        "cow",
        "sheep",
        "ferret",
        "squirrel",
        "opossum",
        "tree_shrew",
        "treeshrew",
        "lizard",
        "anolis",
        "peromyscus",
        "deer_mouse",
        "rhabdomys",
        "striped_mouse",
        # mus_musculus is used by some legacy mouse configs
        "mus_musculus",
    }
)


def _validate_species(cfg) -> None:
    """Log a warning if *cfg.species* is not a recognised pipeline key."""
    species = getattr(cfg, "species", "")
    if not species:
        return  # empty species = user hasn't set it yet; not necessarily wrong
    if species not in _KNOWN_SPECIES_KEYS:
        logging.getLogger("core").warning(
            "CFG.species=%r is not a recognised species key. "
            "KB marker scoring and phylogenetic weighting may silently fail. "
            "Add this species to _SPECIES_SYNONYMS in rna/utils/marker_scoring.py "
            "if it is a valid species.",
            species,
        )
