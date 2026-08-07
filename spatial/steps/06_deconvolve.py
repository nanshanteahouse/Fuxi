#!/usr/bin/env python3
"""
Step 06: cell2location spatial deconvolution
=============================================
  Infer per-spot cell-type proportions from a scRNA-seq reference.

  The spot data is paired with a scRNA reference (``cfg.rna_ref`` or
  auto-discovered via ``find_rna_h5ad``).  The reference's per-cell-type
  expression signatures are estimated with ``RegressionModel``, then
  ``Cell2location`` deconvolves each spot into cell-type abundances which
  are normalised into a spot × cell-type proportion matrix.

  Key conventions (verified against scvi-tools 1.5 — see
  ``.omo/evidence/phase2-env-compat.md``):
    * ``train()``/``export_posterior()`` use ``accelerator=`` and NEVER
      ``use_gpu=`` / ``train_aggressive`` (removed in scvi-tools 1.x).
    * ``cell_state_df`` must be the DataFrame sliced from
      ``ref.varm["means_per_cluster_mu_fg"]`` — never a model object.

  Gene space: spatial var_names may be ENSEMBL ids (ENSG) while the
  reference uses symbols — a mygene lookup maps ENSG → symbol before
  alignment.  If the shared-gene overlap drops below 50% the step records a
  failed status (state preserved) and exits 1 (RNA quality-gate convention).

  Optional NMF tissue zones: when ``cfg.spatial.deconv_n_factors > 0`` a
  ``CoLocatedGroupsSklearnNMF`` over the abundance matrix assigns each spot
  a ``tissue_zone`` and exports the per-zone cell-type composition.

Input:  05_annotated.h5ad  (+ optional scRNA reference h5ad)
Output: 06_deconvolved.h5ad (obsm['deconv_proportions'], obs['tissue_zone'],
        uns['deconvolution']), proportions.csv, tissue_zone_composition.csv
"""

import argparse
import os
import re
import sys
import time

# Add repo root so `from core.*` resolves correctly
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
import pandas as pd
import scanpy as sc

from core.utils import find_rna_h5ad, resolve_config, safe_write, setup_logger, timed_substep

# Minimum fraction of spot genes that must be shared with the reference.
# Below this the step treats the pair as incompatible (RNA quality gate).
MIN_GENE_OVERLAP = 0.5

# Gene set retained when the reference has an HVG column (training speed).
MIN_HVG_FOR_CAP = 1000

_ENSEMBL_RE = re.compile(r"^ENSG\d+(\.\d+)?$", re.IGNORECASE)


def _cfg_val(cfg, section: str, name: str, default):
    """Backward-compatible config read (schema fields land in Wave 3).

    Pydantic configs raise AttributeError for missing fields, so a plain
    ``getattr(..., default)`` returns the default.  None/"" collapse to the
    default as well.
    """
    section_obj = getattr(cfg, section, None) if section else cfg
    if section_obj is None:
        return default
    val = getattr(section_obj, name, default)
    return default if val is None or val == "" else val


def _get_models(log=None):
    """Lazily import cell2location.models; ``None`` when unavailable.

    Imported lazily (not at module top) so the step can run and report a
    graceful skip on machines without cell2location.
    """
    try:
        from cell2location import models

        return models
    except ImportError:
        if log is not None:
            log.warning("cell2location not importable — spatial deconvolution will be skipped")
        return None


def _resolve_accelerator(cfg) -> str:
    """Map cfg.execution.device ('gpu'|'cpu'|'auto') to a lightning accelerator."""
    device = _cfg_val(cfg, "execution", "device", "auto")
    if device == "gpu":
        return "gpu"
    if device == "cpu":
        return "cpu"
    return "auto"


def _resolve_label_col(ref_adata, log) -> str | None:
    """Pick the reference's cell-type label column, tolerating aliases.

    Seurat-derived references commonly store the label under ``CellType``
    (or ``celltype``/``cell_types``) rather than scanpy's ``cell_type``.
    Returns ``None`` when no usable label column exists."""
    for candidate in ("cell_type", "CellType", "celltype", "cell_types", "leiden"):
        if candidate in ref_adata.obs:
            log.info("Reference cell-type label column: '%s'", candidate)
            return candidate
    return None


def _resolve_rna_ref(cfg, log):
    """Resolve the scRNA reference h5ad path.

    Search order:
      1. ``cfg.rna_ref`` — existing file path, or a bare dataset ID /
         directory resolved via ``find_rna_h5ad``.
      2. ``find_rna_h5ad`` auto-discovery (projects/rna/{dataset_id}).
    Returns ``None`` when no reference can be found.
    """
    ref = _cfg_val(cfg, None, "rna_ref", "") or ""
    if ref:
        if os.path.isfile(ref):
            log.info("RNA reference (explicit path): %s", ref)
            return ref
        dataset_id = os.path.basename(os.path.normpath(ref)) if os.path.isdir(ref) else ref
        found = find_rna_h5ad(cfg=cfg, dataset_id=dataset_id, log=log)
        if found:
            log.info("RNA reference resolved from rna_ref='%s': %s", ref, found)
            return found
        log.warning("rna_ref='%s' did not resolve to an existing RNA h5ad", ref)
        return None

    found = find_rna_h5ad(cfg=cfg, log=log)
    if found:
        log.info("RNA reference auto-discovered: %s", found)
        return found
    log.warning("No scRNA reference found (cfg.rna_ref empty and auto-discovery found nothing)")
    return None


def _looks_like_ensembl(var_names) -> bool:
    """True when >90% of var_names look like ENSEMBL gene ids (ENSG...)."""
    if len(var_names) == 0:
        return False
    hits = sum(1 for v in var_names if _ENSEMBL_RE.match(str(v)))
    return hits / len(var_names) > 0.9


def _mygene_query(ensg_ids, log, chunk_size: int = 500):
    """Batch ENSG → symbol lookup via mygene (Ensembl gene scope).

    Returns ``{ensg: symbol}`` (first symbol per id, no duplicates) or
    ``None`` when mygene is unavailable / returns nothing.  Query failures
    on a chunk are logged and skipped (genes in that chunk are dropped).
    """
    try:
        import mygene
    except ImportError:
        log.warning(
            "mygene not installed — cannot map ENSG→symbol "
            "(install with: pip install 'mygene>=3.2')"
        )
        return None

    mg = mygene.MyGeneInfo()
    mapping: dict = {}
    for start in range(0, len(ensg_ids), chunk_size):
        chunk = list(ensg_ids[start : start + chunk_size])
        try:
            results = mg.querymany(
                chunk,
                scopes="ensembl.gene",
                fields="symbol",
                species="human",
                verbose=False,
            )
        except Exception as e:  # noqa: BLE001 - any query failure drops this chunk
            log.warning("mygene querymany failed for %d genes: %s", len(chunk), e)
            continue
        for r in results:
            if not isinstance(r, dict):
                continue
            q = r.get("query")
            if q and not r.get("notfound", False):
                sym = r.get("symbol")
                if isinstance(sym, str) and sym:
                    mapping.setdefault(q, sym)
    if not mapping:
        log.warning("mygene returned no ENSG→symbol mappings")
        return None
    return mapping


def _align_genes(spot_adata, ref_adata, log):
    """Align spot & reference to a common symbol gene set.

    Returns ``(spot_counts, ref_counts, n_genes_mapped)`` or ``None`` when
    the overlap gate (<50%) fails.  ``n_genes_mapped`` counts spot genes
    that survive alignment (unmapped / absent-from-reference genes dropped).
    """
    spot_var = spot_adata.var_names
    ref_symbols = set(ref_adata.var_names)

    if _looks_like_ensembl(spot_var):
        log.info("Spatial var_names look like ENSEMBL ids — mapping ENSG→symbol via mygene")
        mapping = _mygene_query(list(spot_var), log)
        if not mapping:
            log.error("ENSG→symbol mapping failed (mygene returned nothing)")
            return None
        # First symbol wins per gene; drop unmapped; dedupe symbols
        seen = set()
        symbol_map = {}
        for g in spot_var:
            s = mapping.get(g)
            if s and s not in seen:
                seen.add(s)
                symbol_map[g] = s
        mapped_df = pd.DataFrame(
            [(g, s) for g, s in symbol_map.items()],
            columns=["ensg", "symbol"],
        )
        present = mapped_df[mapped_df["symbol"].isin(ref_symbols)]
        overlap = len(present) / len(spot_var)
        log.info(
            "  ENSG→symbol: %d/%d mapped, %d present in reference (%.1f%% overlap)",
            len(mapped_df),
            len(spot_var),
            len(present),
            overlap * 100,
        )
        if overlap < MIN_GENE_OVERLAP:
            log.error(
                "Gene overlap %.1f%% < %.0f%% — cannot run deconvolution",
                overlap * 100,
                MIN_GENE_OVERLAP * 100,
            )
            return None
        spot_sub = spot_adata[:, present["ensg"].tolist()].copy()
        spot_sub.var_names = present["symbol"].tolist()
        n_mapped = len(spot_sub.var_names)
    else:
        common = spot_var.intersection(ref_symbols)
        overlap = len(common) / len(spot_var)
        log.info(
            "  Symbol alignment: %d/%d spot genes shared with reference (%.1f%% overlap)",
            len(common),
            len(spot_var),
            overlap * 100,
        )
        if overlap < MIN_GENE_OVERLAP:
            log.error(
                "Gene overlap %.1f%% < %.0f%% — cannot run deconvolution",
                overlap * 100,
                MIN_GENE_OVERLAP * 100,
            )
            return None
        spot_sub = spot_adata[:, common.tolist()].copy()
        n_mapped = len(common)

    ref_sub = ref_adata[:, spot_sub.var_names].copy()
    ref_sub, spot_sub = _maybe_hvg_subset(ref_sub, spot_sub, log)

    spot_counts = _extract_counts(spot_sub, log)
    ref_counts = _extract_counts(ref_sub, log)
    log.info(
        "  Aligned gene set: %d genes (%d spot genes mapped)",
        spot_counts.n_vars,
        n_mapped,
    )
    return spot_counts, ref_counts, n_mapped


def _maybe_hvg_subset(ref_sub, spot_sub, log):
    """Cap the shared gene set to reference HVGs when available (speed)."""
    if "highly_variable" not in ref_sub.var:
        return ref_sub, spot_sub
    hvg = ref_sub.var_names[ref_sub.var["highly_variable"].values.astype(bool)]
    if len(hvg) >= MIN_HVG_FOR_CAP and len(hvg) < len(ref_sub.var_names):
        log.info(
            "  HVG subset: using %d highly-variable genes (of %d shared)",
            len(hvg),
            len(ref_sub.var_names),
        )
        return ref_sub[:, hvg], spot_sub[:, hvg]
    return ref_sub, spot_sub


def _extract_counts(adata, log):
    """Return an AnnData carrying count-like data for cell2location.

    Priority: ``adata.raw`` (full-gene counts; spatial keeps raw counts) >
    ``layers['counts']`` > ``adata.X``.  The result is subset to the current
    var_names and keeps obs (batch/label columns) intact.
    """
    if adata.raw is not None:
        try:
            out = adata.raw.to_adata()
            if out.n_vars != adata.n_vars:
                out = out[:, adata.var_names]
            log.info("  Counts source: adata.raw")
            return out
        except Exception as e:  # noqa: BLE001 - fall back to the next source
            log.warning("  adata.raw extraction failed (%s) — falling back", e)
    if "counts" in adata.layers:
        out = adata.copy()
        out.X = out.layers["counts"]
        log.info("  Counts source: adata.layers['counts']")
        return out
    log.info("  Counts source: adata.X (no raw counts layer found)")
    return adata


def _run_reference_estimation(models, ref_counts, accelerator, cfg, log):
    """RegressionModel reference estimation → inf_aver (gene × cell-type)."""
    labels_key = _resolve_label_col(ref_counts, log) or "leiden"
    n_before = ref_counts.n_obs
    labeled = ref_counts.obs[labels_key].notna()
    if not labeled.all():
        ref_counts = ref_counts[labeled].copy()
        log.warning(
            "Dropped %d reference cells with missing '%s' label (%d remain)",
            n_before - ref_counts.n_obs,
            labels_key,
            ref_counts.n_obs,
        )
    batch_key = "sample" if "sample" in ref_counts.obs else None
    max_epochs = int(_cfg_val(cfg, "spatial", "deconv_ref_max_epochs", 500))
    n_samples = int(_cfg_val(cfg, "spatial", "deconv_ref_n_samples", 1000))
    batch_size = int(_cfg_val(cfg, "spatial", "deconv_ref_batch_size", 2500))
    log.info(
        "Reference estimation: %d cells × %d genes, %d cell types (labels_key='%s')",
        ref_counts.n_obs,
        ref_counts.n_vars,
        ref_counts.obs[labels_key].nunique(),
        labels_key,
    )
    models.RegressionModel.setup_anndata(ref_counts, labels_key=labels_key, batch_key=batch_key)
    sc_model = models.RegressionModel(ref_counts)
    with timed_substep("RegressionModel reference training", log=log):
        sc_model.train(max_epochs=max_epochs, accelerator=accelerator)
    with timed_substep("Reference posterior export", log=log):
        ref_out = sc_model.export_posterior(
            ref_counts,
            sample_kwargs={"num_samples": n_samples, "batch_size": batch_size},
        )
    return _extract_inf_aver(ref_out, log)


def _extract_inf_aver(ref_out, log):
    """Slice per-cell-type expression signatures from the reference posterior.

    Standard cell2location shape: genes × cell types
    (``ref_out.varm['means_per_cluster_mu_fg']`` sliced by the factor_names
    recorded in ``ref_out.uns['mod']``).  Returns ``None`` when the expected
    slots are missing.
    """
    try:
        varm = ref_out.varm["means_per_cluster_mu_fg"]
        factors = list(ref_out.uns["mod"]["factor_names"])
    except (KeyError, AttributeError, TypeError) as e:
        log.error(
            "Reference posterior missing expected slots "
            "(varm['means_per_cluster_mu_fg'] / uns['mod']['factor_names']): %s",
            e,
        )
        return None
    cols = [f"means_per_cluster_mu_fg_{f}" for f in factors]
    cols = [c for c in cols if c in varm.columns]
    if not cols:
        log.error(
            "No means_per_cluster_mu_fg columns matched factor_names — "
            "reference estimation produced no signatures"
        )
        return None
    inf_aver = varm[cols].copy()
    inf_aver.columns = [c[len("means_per_cluster_mu_fg_") :] for c in cols]
    log.info("  Reference signatures: %d genes × %d cell types", *inf_aver.shape)
    return inf_aver


def _run_deconvolution(models, spot_counts, inf_aver, accelerator, cfg, log):
    """Cell2location spatial deconvolution → absolute cell-type abundance.

    Returns a spot × cell-type abundance DataFrame (posterior 5% quantile,
    ``obsm['q05_cell_abundance_w_sf']``) or ``None`` on failure.
    """
    max_epochs = int(_cfg_val(cfg, "spatial", "deconv_max_epochs", 30000))
    batch_size = int(_cfg_val(cfg, "spatial", "deconv_batch_size", 256))
    n_cells = int(_cfg_val(cfg, "spatial", "deconv_n_cells_per_location", 8))
    alpha = float(_cfg_val(cfg, "spatial", "deconv_detection_alpha", 20.0))
    batch_key = "sample" if "sample" in spot_counts.obs else None
    log.info(
        "Spatial deconvolution: %d spots × %d genes (%d cell types), "
        "N_cells_per_location=%d, detection_alpha=%.0f",
        spot_counts.n_obs,
        spot_counts.n_vars,
        inf_aver.shape[1],
        n_cells,
        alpha,
    )
    models.Cell2location.setup_anndata(spot_counts, batch_key=batch_key)
    mod = models.Cell2location(
        spot_counts,
        cell_state_df=inf_aver,
        N_cells_per_location=n_cells,
        detection_alpha=alpha,
    )
    with timed_substep("Cell2location training", log=log):
        mod.train(max_epochs=max_epochs, accelerator=accelerator, batch_size=batch_size)
    with timed_substep("Spatial posterior export", log=log):
        spot_out = mod.export_posterior(
            spot_counts,
            sample_kwargs={"num_samples": 1000, "batch_size": batch_size},
        )
    try:
        abundance = spot_out.obsm["q05_cell_abundance_w_sf"]
    except (KeyError, AttributeError) as e:
        log.error("Deconvolution posterior missing obsm['q05_cell_abundance_w_sf']: %s", e)
        return None
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(
            abundance,
            index=spot_counts.obs_names,
            columns=inf_aver.columns,
        )
    log.info("  Estimated abundance: %d spots × %d cell types", *abundance.shape)
    return abundance


def _proportions_from_abundance(abundance) -> pd.DataFrame:
    """Normalise abundance rows to proportions (each row sums to 1).

    All-zero rows are left as 0 (never fabricated).
    """
    row_sum = abundance.sum(axis=1).replace(0, np.nan)
    return abundance.div(row_sum, axis=0).fillna(0.0)


def _run_nmf_zones(models, adata, abundance, cfg, log):
    """CoLocatedGroupsSklearnNMF over abundances → obs['tissue_zone'].

    Only runs when ``cfg.spatial.deconv_n_factors > 0``.  Exports the
    per-zone cell-type composition CSV (cell types × factors).
    """
    n_factors = int(_cfg_val(cfg, "spatial", "deconv_n_factors", 0))
    if n_factors <= 0:
        return None
    if abundance.shape[0] < n_factors or abundance.shape[1] < 2:
        log.warning(
            "NMF tissue zones skipped: need ≥%d spots and ≥2 cell types (have %d×%d)",
            n_factors,
            abundance.shape[0],
            abundance.shape[1],
        )
        return None
    random_state = int(_cfg_val(cfg, "execution", "random_seed", 42))
    log.info(
        "NMF tissue zones: %d factors on %d spots × %d cell types",
        n_factors,
        abundance.shape[0],
        abundance.shape[1],
    )
    with timed_substep("CoLocatedGroupsSklearnNMF", log=log):
        nmf_model = models.CoLocatedGroupsSklearnNMF(
            n_fact=n_factors,
            X_data=abundance.values,
            n_iter=10000,
            random_state=random_state,
            var_names=abundance.columns.tolist(),
            obs_names=abundance.index.tolist(),
            fact_names=[f"factor{i}" for i in range(n_factors)],
        )
        nmf_model.fit(n=3, n_type="restart")
        nmf_model.sample2df(node_name="location_factors", ct_node_name="cell_type_factors")

    zone = nmf_model.location_factors_df.idxmax(axis=1)
    zone.name = "tissue_zone"
    adata.obs["tissue_zone"] = zone.astype("category")

    composition = nmf_model.cell_type_fractions  # cell types × factors
    comp_csv = os.path.join(cfg.table_dir, "tissue_zone_composition.csv")
    composition.to_csv(comp_csv)
    log.info("  tissue_zone assigned to %d spots; composition → %s", len(zone), comp_csv)
    return zone


def _record_deconvolution(adata, *, status: str, method: str, **extra) -> None:
    rec = {
        "method": method,
        "status": status,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    rec.update(extra)
    adata.uns["deconvolution"] = rec


def _finish_skipped(adata, cfg, log, output_path, method, reason) -> None:
    log.warning("Deconvolution skipped: %s", reason)
    _record_deconvolution(adata, status="skipped", method=method, reason=reason)
    safe_write(adata, output_path, cfg=cfg)
    log.info("State recorded (status='skipped') → %s", output_path)


def _finish_failed(adata, cfg, log, output_path, method, reason) -> None:
    log.error("Deconvolution failed: %s", reason)
    _record_deconvolution(adata, status="failed", method=method, reason=reason)
    safe_write(adata, output_path, cfg=cfg)
    log.info("State preserved (status='failed') → %s", output_path)
    sys.exit(1)


def main():
    t0 = time.time()
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--config", default="../config.py")
    args = args_parser.parse_args()

    cfg = resolve_config(args.config)
    log = setup_logger("06_deconvolve", os.path.join(cfg.log_dir, "06_deconvolve.log"))
    log.info("Step 06: cell2location spatial deconvolution")

    input_path = os.path.join(cfg.h5ad_dir, "05_annotated.h5ad")
    output_path = os.path.join(cfg.h5ad_dir, "06_deconvolved.h5ad")
    if os.path.exists(output_path):
        log.info("Skip: %s already exists.", output_path)
        return
    if not os.path.exists(input_path):
        log.error("Input not found: %s. Run Step 05 first.", input_path)
        sys.exit(1)

    adata = sc.read(input_path)
    log.info("Loaded: %d spots × %d genes", adata.n_obs, adata.n_vars)

    method = _cfg_val(cfg, "spatial", "deconv_method", "cell2location") or "none"
    if method in ("none", ""):
        _finish_skipped(adata, cfg, log, output_path, method, reason=f"deconv_method='{method}'")
        return

    models = _get_models(log=log)
    if models is None:
        _finish_skipped(
            adata,
            cfg,
            log,
            output_path,
            method,
            reason="cell2location not installed (pip install fuxi[scvi])",
        )
        return

    ref_path = _resolve_rna_ref(cfg, log)
    if not ref_path:
        _finish_skipped(
            adata,
            cfg,
            log,
            output_path,
            method,
            reason="no scRNA reference found (cfg.rna_ref / auto-discovery)",
        )
        return

    ref_adata = sc.read(ref_path)
    log.info(
        "Reference loaded: %d cells × %d genes from %s",
        ref_adata.n_obs,
        ref_adata.n_vars,
        ref_path,
    )
    label_col = _resolve_label_col(ref_adata, log)
    if label_col is None:
        _finish_failed(
            adata,
            cfg,
            log,
            output_path,
            method,
            reason="reference lacks a cell-type label column (cell_type/CellType/leiden)",
        )

    accelerator = _resolve_accelerator(cfg)
    log.info("accelerator='%s' (from cfg.execution.device)", accelerator)

    aligned = _align_genes(adata, ref_adata, log)
    if aligned is None:
        _finish_failed(
            adata,
            cfg,
            log,
            output_path,
            method,
            reason="gene overlap < 50% (see log)",
        )
    spot_counts, ref_counts, n_genes_mapped = aligned

    with timed_substep("Reference estimation (RegressionModel)", log=log):
        inf_aver = _run_reference_estimation(models, ref_counts, accelerator, cfg, log)
    if inf_aver is None:
        _finish_failed(
            adata,
            cfg,
            log,
            output_path,
            method,
            reason="reference estimation produced no cell-state signatures",
        )

    with timed_substep("Spatial deconvolution (Cell2location)", log=log):
        abundance = _run_deconvolution(models, spot_counts, inf_aver, accelerator, cfg, log)
    if abundance is None:
        _finish_failed(
            adata,
            cfg,
            log,
            output_path,
            method,
            reason="spatial deconvolution produced no abundance estimates",
        )

    proportions = _proportions_from_abundance(abundance)
    adata.obsm["deconv_proportions"] = proportions.values
    adata.obsm["q05_cell_abundance"] = abundance.values

    os.makedirs(cfg.table_dir, exist_ok=True)
    proportions.to_csv(os.path.join(cfg.table_dir, "proportions.csv"))
    log.info(
        "Proportion matrix: %d spots × %d cell types → proportions.csv",
        *proportions.shape,
    )

    n_factors = int(_cfg_val(cfg, "spatial", "deconv_n_factors", 0))
    if n_factors > 0:
        _run_nmf_zones(models, adata, abundance, cfg, log)

    _record_deconvolution(
        adata,
        status="completed",
        method=method,
        rna_ref=ref_path,
        accelerator=accelerator,
        n_ref_cells=int(ref_adata.n_obs),
        n_cell_types=int(abundance.shape[1]),
        n_spots=int(adata.n_obs),
        n_genes_mapped=n_genes_mapped,
        n_genes_used=int(spot_counts.n_vars),
        max_epochs=int(_cfg_val(cfg, "spatial", "deconv_max_epochs", 30000)),
        n_factors=n_factors,
        wall_sec=round(time.time() - t0, 1),
    )
    safe_write(adata, output_path, cfg=cfg)
    log.info("Deconvolution results → %s", output_path)
    log.info("Step 06 complete, took %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
