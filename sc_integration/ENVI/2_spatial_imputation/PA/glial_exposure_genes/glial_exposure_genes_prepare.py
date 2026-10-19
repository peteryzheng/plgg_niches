#!/usr/bin/env python3
# Gene-level myeloid glial-adjacency contrast (Python streaming stage) -- pilocytic astrocytoma.
#
# Purpose: the companion notebook PA_myeloid_mapk_neighborhood.qmd asks whether a
# myeloid cell's programs shift when it sits next to PA-niche (MAPK-high) glia, but
# it only ever tests 56 *signatures*. A collaborator wants the underlying GENES, so
# this stage runs the identical per-sample contrast one gene at a time over the
# ENVI-imputed transcriptome and emits per-sample standardized effect sizes.
#
# WHY CLOSED-FORM (no lm per gene): Model 1 is `y ~ exposure_group`, a single
# categorical predictor, so the pa_niche_glia coefficient IS mean(pa) - mean(other),
# and its SE is sqrt(s^2 * (1/n_pa + 1/n_other)) with s^2 the residual MS pooled over
# all three groups (df = n_myeloid - 3). Both follow from per-group sums and
# sums-of-squares, so one streaming pass reproduces lm exactly -- no 18,077 model fits.
# Model 2 (`+ mean_glial_mapk`) has a design matrix fixed per sample, so (X'X)^-1 is
# precomputed once and each gene needs only three dot products (Sy, S(g*y), S(mapk*y)).
#
# MEMORY NOTE: obsm['imputation'] is stored by anndata as a DataFrame group -- one
# HDF5 dataset per gene -- so we stream it gene-by-gene with h5py (idiom lifted from
# CL/lipid_state_gsea/lipid_state_prepare.py). Peak memory is one gene column plus the
# dense 266-gene panel, so this runs in a few hundred MB and does NOT need a big-RAM node.
#
# Inputs
#   - PA/myeloid_glial_exposure.tsv -- the FINISHED per-cell design from the notebook
#     (cell_id, sample_id, n_glial, n_pa_glial, mean_glial_mapk, exposure_group).
#     The kNN graph, the programmatic PA-glial-niche derivation and the group
#     assignment are already done there; this stage re-derives none of it.
#   - Phase-0 cell metadata TSV -- cell_type, used only to locate glial cells for the
#     spillover control.
#   - PA ENVI spatial h5ads: obsm['imputation'] (imputed transcriptome, per-gene
#     datasets) and X (the 266-gene MEASURED panel, already depth-normalized).
#
# Outputs (written next to this script, consumed by glial_exposure_genes.R)
#   - *_per_sample.tsv         sample x gene x source long table of effect sizes
#   - *_sample_diagnostics.tsv per-sample cell counts, join rate, gene universe size
#
# Key assumptions / decisions
#   - sd_M (the standardization denominator) is computed over ALL myeloid cells in the
#     sample INCLUDING no_glial, because the notebook z-scores per patient before
#     Model 2 drops those cells. Both models divide by that same sd_M, so the
#     Model-1 -> Model-2 attenuation stays comparable to the notebook's.
#   - Spillover is reported as glial_share = mean_glial / (mean_glial + mean_myeloid),
#     bounded in [0,1], rather than a log ratio: imputed values are tiny (many gene
#     means ~1e-4) so a hardcoded epsilon in a log ratio would silently reorder most of
#     the 18k genes depending on its magnitude.
#   - mean_glial_mapk (the Model-2 covariate) was computed on the OLD imputation. It is
#     reused as-is for a sensitivity arm; the provenance mismatch is documented.

import argparse
import glob
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# Glial lineages, finalized lowercase labels -- must match the notebook's glial_types.
GLIAL_TYPES = ("astro_glial", "opc", "radial_glia")

# Reference level is other_glia, so the reported beta is the PA-specific effect
# (near PA-niche glia vs near OTHER glia), exactly as in the notebook's Model 1.
GROUP_PA = "pa_niche_glia"
GROUP_OTHER = "other_glia"
GROUP_NONE = "no_glial"
GROUPS = (GROUP_PA, GROUP_OTHER, GROUP_NONE)

MIN_N_GRP = 30  # min cells per contrasted group; same gate as the notebook


def detect_workdir() -> str:
    home = os.path.expanduser("~")
    if home in ["/Users/youyun", "/Users/youyunzheng"]:
        return home + "/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
    if home == "/home/yz762":
        return "/mnt/storage/dept/medonc/beroukhim/"
    if home == "/PHShome/yz762":
        return "/data/beroukhim1/"
    return "/xchip/beroukhimlab/"


def _decode(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def read_obs_ids_and_genes(h5path: str) -> tuple:
    # CHEAP pass: obs ids, imputed gene names (+ their HDF5 keys) and panel var names.
    # Touches no matrix data.
    with h5py.File(h5path, "r") as h:
        obs = h["obs"]
        iname = _decode(obs.attrs.get("_index", "_index"))
        obs_ids = np.array([_decode(x) for x in obs[iname][:]])

        imp = h["obsm"]["imputation"]
        if isinstance(imp, h5py.Dataset):
            raise TypeError(
                f"obsm['imputation'] is a plain matrix (no gene names) in {h5path}; "
                "expected an anndata DataFrame group with per-gene datasets."
            )
        if "column-order" in imp.attrs:
            raw_keys = [_decode(g) for g in imp.attrs["column-order"]]
        else:
            raw_keys = [k for k in imp.keys() if k != "_index"]

        vname = _decode(h["var"].attrs.get("_index", "_index"))
        panel_genes = np.array([_decode(x).upper() for x in h["var"][vname][:]])

    # Uppercased gene -> original HDF5 dataset key, keeping the first on any
    # case-duplicate (mirrors the dedup in the sibling pipelines).
    gene_key_map = {}
    for k in raw_keys:
        gene_key_map.setdefault(k.upper(), k)
    return obs_ids, gene_key_map, panel_genes


def build_sample_context(h5path: str, exposure: pd.DataFrame, cell_meta: pd.DataFrame) -> dict:
    """Resolve every row index / design quantity this sample needs, before touching expression."""
    sample_id = Path(h5path).name.split("_")[0]
    obs_ids, gene_key_map, panel_genes = read_obs_ids_and_genes(h5path)

    exp_s = exposure[exposure["sample_id"] == sample_id]
    if exp_s.empty:
        raise ValueError(f"No exposure-table rows for sample {sample_id}.")

    # Gather myeloid rows in h5ad obs order. This join is the single most dangerous
    # step in the pipeline (a scrambled gather still yields plausible betas), so the
    # join rate is asserted at 100% and reported in the diagnostics.
    obs_pos = pd.Index(obs_ids)
    row_of_cell = obs_pos.get_indexer(exp_s["cell_id"].to_numpy())
    n_missing = int((row_of_cell < 0).sum())
    if n_missing:
        raise ValueError(
            f"{sample_id}: {n_missing}/{len(row_of_cell)} exposure cells absent from h5ad obs; "
            "refusing to proceed on a partial join."
        )

    groups = exp_s["exposure_group"].to_numpy()
    mapk = exp_s["mean_glial_mapk"].to_numpy(dtype=float)

    # Group ids over the myeloid rows: 0 = pa_niche_glia, 1 = other_glia, 2 = no_glial.
    gid = np.full(groups.shape[0], -1, dtype=np.int64)
    for j, g in enumerate(GROUPS):
        gid[groups == g] = j
    if (gid < 0).any():
        bad = sorted(set(groups[gid < 0]))
        raise ValueError(f"{sample_id}: unexpected exposure_group values {bad}.")
    n_by_group = np.bincount(gid, minlength=len(GROUPS))

    # Model-2 subset: glia-adjacent cells only (mean_glial_mapk is undefined for
    # no_glial). Guard against any NaN covariate rather than letting it poison a dot product.
    adj = (gid <= 1) & np.isfinite(mapk)
    n_dropped_mapk = int(((gid <= 1) & ~np.isfinite(mapk)).sum())
    # 1 for pa_niche_glia, 0 for other_glia -> the coefficient reads as pa-vs-other.
    g_adj = (gid[adj] == 0).astype(float)
    m_adj = mapk[adj]

    # (X'X)^-1 for X = [1, group, mapk] over the adjacent subset. Fixed per sample, so
    # each gene only needs X'y. g is 0/1, hence sum(g^2) == sum(g).
    n_adj = g_adj.shape[0]
    xtx = np.array([
        [n_adj,        g_adj.sum(),        m_adj.sum()],
        [g_adj.sum(),  g_adj.sum(),        float(g_adj @ m_adj)],
        [m_adj.sum(),  float(g_adj @ m_adj), float(m_adj @ m_adj)],
    ])
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = None  # singular (e.g. constant covariate) -> Model 2 unavailable

    # Glial rows for the spillover control, from the finalized cell-type labels.
    ct_vec = pd.Index(obs_ids).map(cell_meta["cell_type"]).to_numpy()
    rows_glial = np.flatnonzero(np.isin(ct_vec, GLIAL_TYPES))

    fittable = (n_by_group[0] >= MIN_N_GRP) and (n_by_group[1] >= MIN_N_GRP)

    return {
        "sample_id": sample_id,
        "h5path": h5path,
        "gene_key_map": gene_key_map,
        "panel_genes": panel_genes,
        "rows_myeloid": row_of_cell,
        "gid": gid,
        "adj": adj,
        "g_adj": g_adj,
        "m_adj": m_adj,
        "xtx_inv": xtx_inv,
        "rows_glial": rows_glial,
        "n_by_group": n_by_group,
        "n_adj": n_adj,
        "n_dropped_mapk": n_dropped_mapk,
        "fittable": bool(fittable),
    }


def gene_stats(y_myeloid: np.ndarray, y_glial_mean: float, ctx: dict) -> dict:
    """Model 1 + Model 2 effect sizes for one gene, from per-group sums only.

    Reproduces lm() exactly for Model 1: with a single 3-level factor the
    pa_niche_glia coefficient is mean(pa) - mean(other) and its SE uses the residual
    MS pooled over all three groups (df = n_myeloid - 3).
    """
    gid = ctx["gid"]
    n_grp = ctx["n_by_group"]
    n_m = y_myeloid.shape[0]

    sums = np.bincount(gid, weights=y_myeloid, minlength=len(GROUPS))
    sumsq = np.bincount(gid, weights=y_myeloid * y_myeloid, minlength=len(GROUPS))
    means = np.divide(sums, n_grp, out=np.full(len(GROUPS), np.nan), where=n_grp > 0)

    mean_m = y_myeloid.sum() / n_m
    var_m = (float((y_myeloid * y_myeloid).sum()) - n_m * mean_m ** 2) / (n_m - 1)
    sd_m = float(np.sqrt(var_m)) if var_m > 0 else 0.0

    # Residual MS pooled over the occupied groups == lm's sigma^2.
    occupied = n_grp > 0
    rss = float((sumsq[occupied] - n_grp[occupied] * means[occupied] ** 2).sum())
    df_resid = n_m - int(occupied.sum())
    s2 = rss / df_resid if df_resid > 0 else np.nan

    n_pa, n_other = int(n_grp[0]), int(n_grp[1])
    if n_pa > 0 and n_other > 0:
        diff = float(means[0] - means[1])
        se_raw = float(np.sqrt(s2 * (1.0 / n_pa + 1.0 / n_other))) if np.isfinite(s2) else np.nan
    else:
        diff, se_raw = np.nan, np.nan

    # Standardize by sd_M: dividing y by a constant scales beta and SE identically, so
    # t is unchanged -- this is exactly the notebook's per-patient z-scoring.
    if sd_m > 0:
        beta_m1 = diff / sd_m
        se_m1 = se_raw / sd_m
    else:
        beta_m1 = se_m1 = np.nan
    t_m1 = diff / se_raw if (se_raw and np.isfinite(se_raw) and se_raw > 0) else np.nan

    # Model 2: MAPK-dose-adjusted, glia-adjacent cells only.
    beta_m2 = np.nan
    if ctx["xtx_inv"] is not None and sd_m > 0:
        y_adj = y_myeloid[ctx["adj"]]
        xty = np.array([y_adj.sum(), float(ctx["g_adj"] @ y_adj), float(ctx["m_adj"] @ y_adj)])
        beta_m2 = float((ctx["xtx_inv"] @ xty)[1]) / sd_m

    return {
        "n_pa": n_pa,
        "n_other": n_other,
        "n_no_glial": int(n_grp[2]),
        "n_adj": ctx["n_adj"],
        "mean_pa": float(means[0]),
        "mean_other": float(means[1]),
        "mean_expr_myeloid": float(mean_m),
        "sd_myeloid": sd_m,
        "mean_expr_glial": float(y_glial_mean),
        # bounded spillover index; 0.5 = equally expressed in glia and myeloid
        "glial_share": float(y_glial_mean / (y_glial_mean + mean_m))
        if (y_glial_mean + mean_m) > 0 else np.nan,
        "beta_m1": beta_m1,
        "se_m1": se_m1,
        "t_m1": t_m1,
        "df_m1": df_resid,
        "beta_m2": beta_m2,
    }


def run_sample(ctx: dict, genes: list, do_panel: bool) -> list:
    """Stream the imputed genes (and optionally the measured panel) for one sample."""
    rows_m = ctx["rows_myeloid"]
    rows_g = ctx["rows_glial"]
    out = []

    with h5py.File(ctx["h5path"], "r") as h:
        imp = h["obsm"]["imputation"]
        for i, gene in enumerate(genes):
            col = imp[ctx["gene_key_map"][gene]][:]
            y_m = col[rows_m].astype(np.float64, copy=False)
            g_mean = float(col[rows_g].mean()) if rows_g.size else np.nan
            rec = gene_stats(y_m, g_mean, ctx)
            rec.update(sample_id=ctx["sample_id"], gene=gene, source="imputed")
            out.append(rec)
            if (i + 1) % 2500 == 0:
                print(f"  {ctx['sample_id']}: {i + 1}/{len(genes)} imputed genes", flush=True)

        # On-panel validation arm: the 266 MEASURED genes, immune to imputation
        # artifacts. X is dense and already depth-normalized (constant row sums).
        if do_panel:
            X = h["X"]
            panel = np.asarray(X[:, :]) if isinstance(X, h5py.Dataset) else None
            if panel is None:
                print(f"  {ctx['sample_id']}: X is not dense; skipping panel arm.", flush=True)
            else:
                for j, gene in enumerate(ctx["panel_genes"]):
                    col = panel[:, j]
                    y_m = col[rows_m].astype(np.float64, copy=False)
                    g_mean = float(col[rows_g].mean()) if rows_g.size else np.nan
                    rec = gene_stats(y_m, g_mean, ctx)
                    rec.update(sample_id=ctx["sample_id"], gene=gene, source="panel")
                    out.append(rec)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--envi-dir", default="ENVI_results_091126",
                    help="ENVI result dir under data/sc_integration (the 2026-09-11 full-transcriptome export).")
    ap.add_argument("--out-tag", default="glial_exposure_genes",
                    help="Output filename prefix.")
    ap.add_argument("--genes", default=None,
                    help="Optional comma-separated gene subset (verification runs).")
    ap.add_argument("--shared-only", action="store_true",
                    help="Restrict to genes imputed in ALL samples (a no-op on the current "
                         "export, whose per-sample universes are identical).")
    ap.add_argument("--no-panel", action="store_true", help="Skip the measured-panel arm.")
    args = ap.parse_args()

    workdir = detect_workdir()
    out_dir = Path(__file__).resolve().parent
    pa_dir = out_dir.parent

    data_root = os.path.join(workdir, "youyun/plgg/data/sc_integration", args.envi_dir)
    if not os.path.isdir(data_root):
        raise FileNotFoundError(f"ENVI directory not found: {data_root}")
    spatial_files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_spatial_envi.h5ad")))
    if not spatial_files:
        raise FileNotFoundError(f"No PA spatial h5ads in {data_root}")

    # The finished per-cell design from PA_myeloid_mapk_neighborhood.qmd.
    exposure = pd.read_csv(pa_dir / "myeloid_glial_exposure.tsv", sep="\t",
                           dtype={"sample_id": str, "cell_id": str})
    cell_meta = pd.read_csv(
        os.path.join(workdir, "youyun/plgg/data/xenium_cell_types/"
                              "cell_metadata_resct0.2_resni0.5_20260826.tsv"),
        sep="\t", usecols=["cell_id", "cell_type"]).set_index("cell_id")

    print(f"ENVI dir: {data_root}")
    print(f"Exposure design: {len(exposure)} myeloid cells, "
          f"{exposure['sample_id'].nunique()} samples", flush=True)

    contexts = [build_sample_context(f, exposure, cell_meta) for f in spatial_files]

    diag = pd.DataFrame([{
        "sample_id": c["sample_id"],
        "n_myeloid": int(c["n_by_group"].sum()),
        "n_pa_niche_glia": int(c["n_by_group"][0]),
        "n_other_glia": int(c["n_by_group"][1]),
        "n_no_glial": int(c["n_by_group"][2]),
        "n_glia_adjacent_model2": c["n_adj"],
        "n_dropped_nonfinite_mapk": c["n_dropped_mapk"],
        "n_glial_cells": int(c["rows_glial"].size),
        "n_imputed_genes": len(c["gene_key_map"]),
        "n_panel_genes": int(c["panel_genes"].size),
        "join_rate_pct": 100.0,  # asserted in build_sample_context
        "fittable": c["fittable"],
    } for c in contexts])
    print("\nPer-sample diagnostics:")
    print(diag.to_string(index=False), flush=True)

    contexts = [c for c in contexts if c["fittable"]]
    if not contexts:
        raise ValueError(f"No sample has >= {MIN_N_GRP} cells in both contrasted groups.")

    # Gene universe. The old ENVI run imputed a different ~2.2k HVG set per sample, so
    # a cohort-level cross-sample comparison there requires the intersection.
    gene_sets = [set(c["gene_key_map"].keys()) for c in contexts]
    if args.genes:
        requested = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
        universe = [g for g in requested if all(g in s for s in gene_sets)]
        print(f"\nRequested {len(requested)} genes; {len(universe)} present in all samples.")
    elif args.shared_only:
        universe = sorted(set.intersection(*gene_sets))
        print(f"\nShared imputed universe: {len(universe)} genes")
    else:
        sizes = {len(s) for s in gene_sets}
        if len(sizes) > 1:
            raise ValueError(
                f"Per-sample imputed universes differ ({sorted(sizes)}); pass --shared-only "
                "to restrict to the intersection.")
        universe = sorted(set.intersection(*gene_sets))
        print(f"\nImputed universe (identical across samples): {len(universe)} genes")
    if not universe:
        raise ValueError("Empty gene universe.")

    records = []
    for c in contexts:
        print(f"\nStreaming {c['sample_id']} ({len(universe)} genes)...", flush=True)
        records.extend(run_sample(c, universe, do_panel=not args.no_panel))

    res = pd.DataFrame.from_records(records)
    col_order = ["sample_id", "gene", "source", "n_pa", "n_other", "n_no_glial", "n_adj",
                 "mean_pa", "mean_other", "mean_expr_myeloid", "sd_myeloid",
                 "mean_expr_glial", "glial_share", "beta_m1", "se_m1", "t_m1", "df_m1",
                 "beta_m2"]
    res = res[col_order]

    res_path = out_dir / f"{args.out_tag}_per_sample.tsv"
    diag_path = out_dir / f"{args.out_tag}_sample_diagnostics.tsv"
    res.to_csv(res_path, sep="\t", index=False)
    diag.to_csv(diag_path, sep="\t", index=False)

    print(f"\nWrote {len(res)} rows -> {res_path}")
    print(f"Wrote diagnostics -> {diag_path}")
    print("Rows by source:", res["source"].value_counts().to_dict())


if __name__ == "__main__":
    sys.exit(main())
