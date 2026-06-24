#!/usr/bin/env python3
"""
Step 04: AI-assisted chromatin state annotation
=================================================
  - Reads clustered AnnData
  - Computes marker regions per cluster
  - AI annotation with disk caching to avoid redundant LLM calls

Input:  03_clustered.h5ad
Output: 04_annotated.h5ad
"""

import sys, os, time, argparse, json, hashlib
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from core.utils import setup_logger, resolve_config, safe_write
import numpy as np
import pandas as pd
import snapatac2 as snap


def _cache_key(cfg, cluster_summary):
    """Generate a deterministic cache key from cluster summary content."""
    payload = json.dumps(cluster_summary, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _load_cache(cache_dir, cache_key):
    """Return cached annotation dict or None."""
    if not cache_dir:
        return None
    cache_path = os.path.join(cache_dir, f"ai_annot_{cache_key}.json")
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            return json.load(f)
    return None


def _save_cache(cache_dir, cache_key, annotations):
    """Persist annotation result to disk cache."""
    if not cache_dir:
        return
    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, f"ai_annot_{cache_key}.json")
    with open(cache_path, 'w') as f:
        json.dump(annotations, f, indent=2)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    CFG = resolve_config(args.config)
    log = setup_logger("04_annotate", os.path.join(CFG.log_dir, "04_annotate.log"))
    log.info("Step 04: Chromatin state annotation")

    if os.path.exists(CFG.annotated_h5ad):
        log.info("Skip: %s exists.", CFG.annotated_h5ad)
        return

    # In-memory mode for pandas compatibility
    data = snap.read(CFG.clustered_h5ad)
    if data.isbacked:
        data = data.to_memory()
    log.info("Loaded: %d cells, %d clusters", data.n_obs, data.obs.get('leiden').nunique() if 'leiden' in data.obs else 0)

    cluster_col = None
    for c in data.obs.columns:
        if c.startswith('leiden'):
            cluster_col = c
            break
    if cluster_col is None:
        log.error("No leiden cluster found.")
        sys.exit(1)
    log.info("Using: %s (%d clusters)", cluster_col, data.obs[cluster_col].nunique())

    # ── Marker regions per cluster (SnapATAC2 2.9 returns dict, not stored in uns) ──
    markers = snap.tl.marker_regions(data, groupby=cluster_col, pvalue=CFG.marker_peaks_fdr)
    # markers is Dict[str, pd.Index] — cluster -> peak names

    cluster_summary = []
    for c in sorted(data.obs[cluster_col].unique()):
        mask = data.obs[cluster_col] == c
        top_peaks = []
        if str(c) in markers:
            top_peaks = list(markers[str(c)][:10])
        cluster_summary.append({
            'cluster': int(c),
            'n_cells': int(mask.sum()),
            'top_peaks': top_peaks,
        })
    log.info("Marker regions found for %d clusters", len(cluster_summary))

    # ── AI annotation (with disk cache) ──
    annotations = {}
    if CFG.ai.enabled and CFG.ai.ai_annotation:
        use_cache = getattr(CFG.ai, 'ai_cache_responses', True)
        cache_dir = None
        ck = None
        if use_cache:
            cache_dir = os.path.join(CFG.table_dir if CFG.table_dir else CFG.h5ad_dir,
                                     ".ai_cache")
            ck = _cache_key(CFG, cluster_summary)
            cached = _load_cache(cache_dir, ck)
            if cached:
                log.info("AI annotation loaded from cache (%s)", ck)
                annotations = cached

        if not annotations:
            try:
                from core.ai_caller import ai_query
                from core.ai_prompts import ATAC_ANNOTATION_SYSTEM_PROMPT, ATAC_ANNOTATION_USER_PROMPT_TEMPLATE
                log.info("AI annotation...")
                user_prompt = ATAC_ANNOTATION_USER_PROMPT_TEMPLATE.format(
                    tissue=CFG.tissue,
                    cluster_summary=json.dumps(cluster_summary, indent=2),
                )
                response = ai_query(ATAC_ANNOTATION_SYSTEM_PROMPT, user_prompt, cfg=CFG.ai, log=log)
                annotations = json.loads(response)
                if use_cache and ck:
                    _save_cache(cache_dir, ck, annotations)
            except Exception as e:
                log.warning("AI annotation failed: %s", e)

    # ── Fallback ──
    if not annotations:
        for c in sorted(data.obs[cluster_col].unique()):
            annotations[str(c)] = {'cell_type': f'Cluster_{c}', 'confidence': 'medium'}

    data.obs['cell_type'] = data.obs[cluster_col].astype(str).map(
        lambda x: annotations.get(x, {}).get('cell_type', f'Cluster_{x}'))
    data.obs['annot_confidence'] = data.obs[cluster_col].astype(str).map(
        lambda x: annotations.get(x, {}).get('confidence', 'medium'))

    safe_write(data, CFG.annotated_h5ad, cfg=CFG, compression_override=None)
    log.info("Step 04 complete, took %.1fs", time.time() - t0)


if __name__ == '__main__':
    main()
