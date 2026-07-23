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


def _find_global_yaml() -> str | None:
    """Discover the global.yaml file path.

    Priority:
    1. FUXI_GLOBAL_CONFIG env var (explicit override, full path)
    2. Walk up from CWD looking for repo root (first dir containing .git/)
       and check for <repo_root>/global.yaml
    3. Fall back to None

    Never raises — returns None on any error.
    """
    # 1. Explicit env var override
    env_path = os.environ.get("FUXI_GLOBAL_CONFIG")
    if env_path and os.path.isfile(env_path):
        return env_path

    # 2. Walk up from CWD looking for repo root
    try:
        current = os.path.abspath(os.getcwd())
        for _ in range(20):  # safety limit
            if os.path.isdir(os.path.join(current, ".git")):
                global_yaml = os.path.join(current, "global.yaml")
                if os.path.isfile(global_yaml):
                    return global_yaml
                break  # found repo root, but no global.yaml
            parent = os.path.dirname(current)
            if parent == current:
                break
            current = parent
    except Exception:
        pass

    return None


# ── Deep merge ──────────────────────────────────────────────────────


def deep_merge(base: dict, override: dict, _prefix: str = "") -> tuple[dict, dict]:
    """Recursively merge override into base, returning (merged, source_map).

    source_map is a dict of dot-separated paths -> "base" | "override":
    - "override" if the value came from override dict
    - "base" if retained from base dict

    Rules:
    - Scalar override: override[key] replaces base[key]
    - Dict recursion: if both values are dicts, recurse into them
    - Key retention: keys in base but not in override are kept as-is
    - List replacement: override[key] replaces base[key] (no merging)
    - None handling: None in override does NOT override non-None base values
    - Empty dict handling: {} in override DOES clear nested (logs WARNING)
    """
    merged: dict = {}
    source_map: dict = {}

    all_keys = set(base.keys()) | set(override.keys())

    for key in all_keys:
        full_key = f"{_prefix}.{key}" if _prefix else key
        in_base = key in base
        in_override = key in override

        if in_base and in_override:
            base_val = base[key]
            override_val = override[key]

            if override_val is None:
                # None in override keeps base value
                merged[key] = base_val
                source_map[full_key] = "base"
            elif isinstance(base_val, dict) and isinstance(override_val, dict):
                if not override_val:
                    # Empty dict in override clears nested base values
                    logging.getLogger("core").warning(
                        "deep_merge: override key '%s' is empty dict, clearing nested base values",
                        full_key,
                    )
                    merged[key] = {}
                    source_map[full_key] = "override"
                else:
                    sub_merged, sub_sources = deep_merge(base_val, override_val, full_key)
                    merged[key] = sub_merged
                    source_map.update(sub_sources)
            else:
                # Scalar or list: override wins
                merged[key] = override_val
                source_map[full_key] = "override"

        elif in_base:
            merged[key] = base[key]
            source_map[full_key] = "base"

        else:  # in_override only
            override_val = override[key]
            if override_val is not None:
                merged[key] = override_val
                source_map[full_key] = "override"
            # None-only keys in override are silently ignored

    return merged, source_map


# ── Resolve config ──────────────────────────────────────────────────


def resolve_config(config_path: Optional[str] = None) -> Config:
    """
    Load a YAML config file and return a resolved Config instance.

    Replaces the old importlib-based .py loading with Pydantic v2 YAML loading.
    Supports global.yaml merging: project config values always win over global.
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

    # ── Load project config ──
    with open(config_path) as f:
        project_data = yaml.safe_load(f)

    # Auto-detect project_dir from config file location BEFORE validation
    # (must precede model_validate so model_post_init sees correct project_dir)
    if not project_data.get("project_dir"):
        project_data["project_dir"] = os.path.dirname(config_path)

    # ── Load global.yaml (if available) and merge ──
    global_path = _find_global_yaml()
    global_data: dict = {}
    if global_path:
        with open(global_path) as f:
            global_data = yaml.safe_load(f) or {}

    # Merge: global as base, project as override (project always wins)
    merged, _source_map = deep_merge(global_data, project_data)

    # ── Validate merged config ──
    cfg = Config.model_validate(merged)

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

    # ── Auto-export resolved config (Todo 5) ──
    # _export_resolved(cfg, merged, _source_map, global_path, config_path)

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
