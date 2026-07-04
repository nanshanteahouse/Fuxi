"""Config loading — resolve_config, species validation, dataset.yaml helpers."""

import logging
import os
import sys
from typing import Optional

from core.dataset_schema import load_dataset


# ── Dataset.yaml helpers ───────────────────────────────────────────


def _find_dataset_yaml(cfg) -> Optional[str]:
    """Search project_dir → data_dir for dataset.yaml."""
    for base in [cfg.project_dir, cfg.data_dir]:
        if base:
            path = os.path.join(base, 'dataset.yaml')
            if os.path.exists(path):
                return path
    return None


def _has_explicit_is_nuclei(config_path: str) -> bool:
    """Check if config.py explicitly references 'is_nuclei'."""
    try:
        with open(config_path) as f:
            return 'is_nuclei' in f.read()
    except OSError:
        return False


def resolve_config(config_path: Optional[str] = None):
    """
    解析 --config 参数，返回配置模块的 CFG 对象。

    所有步骤脚本统一使用本函数加载配置。
    """
    if config_path is None:
        # 默认寻找父目录的 config.py
        config_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            "config.py",
        )

    config_path = os.path.abspath(config_path)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"配置文件不存在: {config_path}")

    import importlib.util
    spec = importlib.util.spec_from_file_location("pipeline_config", config_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_config"] = mod
    spec.loader.exec_module(mod)

    # Auto-detect project_dir from config file location if not explicitly set.
    if not mod.CFG.project_dir:
        mod.CFG.project_dir = os.path.dirname(config_path)

    # ── Resolve n_jobs ───────────────────────────────────────────────
    # 0 means "auto-detect" but joblib.Parallel rejects 0 outright.
    # Resolve here so both standalone step scripts and subprocess steps
    # get a usable value (run_pipeline.py also does this, but steps
    # shouldn't depend on the launcher).
    _nc = getattr(mod.CFG, 'n_jobs', 0)
    if _nc == 0:
        mod.CFG.n_jobs = os.cpu_count() or 1

    mod.CFG.resolve_paths()

    # ── Auto-fill is_nuclei from dataset.yaml ──
    if not _has_explicit_is_nuclei(config_path):
        _yaml = _find_dataset_yaml(mod.CFG)
        if _yaml:
            try:
                _ds = load_dataset(_yaml)
                if getattr(_ds, 'assay_type', None) == 'snRNAseq' and not mod.CFG.is_nuclei:
                    mod.CFG.is_nuclei = True
                    print(f"[Config] Auto-set is_nuclei=True from {_yaml}")
            except Exception:
                pass  # Graceful: no crash on invalid yaml

    # ── Species sanity check ───────────────────────────────────────────
    _validate_species(mod.CFG)

    return mod.CFG


# Species pipeline-keys known to rna/ortholog.py and rna/utils/marker_scoring.py.
# Keep in sync when adding new species.
_KNOWN_SPECIES_KEYS = frozenset({
    "human", "mouse", "macaque", "cynomolgus", "marmoset",
    "zebrafish", "chicken", "lamprey", "frog",
    "pig", "cow", "sheep", "ferret",
    "squirrel", "opossum", "tree_shrew", "treeshrew",
    "lizard", "anolis", "peromyscus", "deer_mouse",
    "rhabdomys", "striped_mouse",
    # mus_musculus is used by some legacy mouse configs
    "mus_musculus",
})


def _validate_species(cfg) -> None:
    """Log a warning if *cfg.species* is not a recognised pipeline key."""
    species = getattr(cfg, 'species', '')
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
