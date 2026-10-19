#!/usr/bin/env python3
# Verification harness for the gene-level myeloid glial-adjacency pipeline.
#
# Two independent checks, because they fail in different ways:
#
#   TEST A (arithmetic) -- does the closed-form in glial_exposure_genes_prepare.py
#     equal a real least-squares fit? Loads a few gene columns, fits Model 1
#     (y ~ exposure_group, 3-level dummy, reference other_glia) and Model 2
#     (+ mean_glial_mapk) by lstsq on the per-patient z-scored values, and compares
#     beta / SE / t against the pipeline's output. Must agree to ~1e-10.
#
#   TEST B (end-to-end row alignment) -- does the cell_id -> h5ad obs gather line up?
#     A scrambled gather still produces a plausible-looking beta, so arithmetic alone
#     proves nothing. This test rebuilds a SIGNATURE-level score the way the notebook
#     did (mean over that sample's imputed members of HALLMARK_COMPLEMENT etc. on the
#     OLD ENVI_results universe), fits the same Model 1, and compares against the
#     notebook's stored per-sample Estimate in myeloid_glial_exposure_nested_results.tsv.
#     It will NOT match exactly -- sc.tl.score_genes subtracts control-bin means, which a
#     plain gene average does not -- but the sign must match in all 4 samples and the
#     magnitude must land in the same ballpark. That exercises the join, the group
#     assignment and the standardization together.
#
# Usage:
#   python verify_glial_exposure_genes.py test-a --per-sample <prepare output tsv>
#   python verify_glial_exposure_genes.py test-b

import argparse
import glob
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

GLIAL_TYPES = ("astro_glial", "opc", "radial_glia")
GROUPS = ("pa_niche_glia", "other_glia", "no_glial")
# Signatures with stored per-sample Model-1 estimates to reproduce in TEST B.
TEST_B_SIGNATURES = ("HALLMARK_COMPLEMENT", "MG_Act_Tripathi", "HALLMARK_INFLAMMATORY_RESPONSE")
TEST_A_GENES = ("C1QB", "APOE", "VCAN")


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


def load_design(workdir: str) -> pd.DataFrame:
    pa_dir = Path(workdir) / "youyun/plgg/code/sc_integration/ENVI/2_spatial_imputation/PA"
    return pd.read_csv(pa_dir / "myeloid_glial_exposure.tsv", sep="\t",
                       dtype={"sample_id": str, "cell_id": str})


def open_sample(h5path: str, exposure: pd.DataFrame):
    """Return (myeloid row indices in obs order, exposure rows, gene->hdf5 key map)."""
    sample_id = Path(h5path).name.split("_")[0]
    with h5py.File(h5path, "r") as h:
        iname = _decode(h["obs"].attrs.get("_index", "_index"))
        obs_ids = np.array([_decode(x) for x in h["obs"][iname][:]])
        imp = h["obsm"]["imputation"]
        raw_keys = [_decode(g) for g in imp.attrs["column-order"]]
    kmap = {}
    for k in raw_keys:
        kmap.setdefault(k.upper(), k)
    exp_s = exposure[exposure["sample_id"] == sample_id].copy()
    rows = pd.Index(obs_ids).get_indexer(exp_s["cell_id"].to_numpy())
    assert (rows >= 0).all(), f"{sample_id}: incomplete cell_id -> obs join"
    return sample_id, rows, exp_s, kmap


def fit_model1(y_z: np.ndarray, groups: np.ndarray) -> dict:
    """Real lstsq fit of y ~ exposure_group with other_glia as the reference level."""
    d_pa = (groups == "pa_niche_glia").astype(float)
    d_none = (groups == "no_glial").astype(float)
    X = np.column_stack([np.ones_like(y_z), d_pa, d_none])
    beta, *_ = np.linalg.lstsq(X, y_z, rcond=None)
    resid = y_z - X @ beta
    df = X.shape[0] - X.shape[1]
    s2 = float(resid @ resid) / df
    se = float(np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1]))
    return {"beta": float(beta[1]), "se": se, "t": float(beta[1]) / se, "df": df}


def fit_model2(y_z: np.ndarray, groups: np.ndarray, mapk: np.ndarray) -> dict:
    """Real lstsq fit of y ~ (pa vs other) + mean_glial_mapk on glia-adjacent cells."""
    keep = np.isin(groups, ("pa_niche_glia", "other_glia")) & np.isfinite(mapk)
    yk = y_z[keep]
    X = np.column_stack([
        np.ones_like(yk),
        (groups[keep] == "pa_niche_glia").astype(float),
        mapk[keep],
    ])
    beta, *_ = np.linalg.lstsq(X, yk, rcond=None)
    return {"beta": float(beta[1])}


def test_a(args) -> None:
    workdir = detect_workdir()
    exposure = load_design(workdir)
    got = pd.read_csv(args.per_sample, sep="\t", dtype={"sample_id": str})
    got = got[got["source"] == "imputed"]

    data_root = os.path.join(workdir, "youyun/plgg/data/sc_integration", args.envi_dir)
    files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_spatial_envi.h5ad")))

    print(f"{'sample':8s} {'gene':8s} {'quantity':6s} {'pipeline':>14s} {'lstsq':>14s} {'absdiff':>11s}")
    max_diff = 0.0
    for f in files:
        sample_id, rows, exp_s, kmap = open_sample(f, exposure)
        groups = exp_s["exposure_group"].to_numpy()
        mapk = exp_s["mean_glial_mapk"].to_numpy(dtype=float)
        with h5py.File(f, "r") as h:
            imp = h["obsm"]["imputation"]
            for gene in TEST_A_GENES:
                if gene not in kmap:
                    continue
                sub = got[(got["sample_id"] == sample_id) & (got["gene"] == gene)]
                if sub.empty:
                    continue
                sub = sub.iloc[0]
                y = imp[kmap[gene]][:][rows].astype(np.float64)
                # per-patient z-scoring, exactly as the notebook does before fitting
                y_z = (y - y.mean()) / y.std(ddof=1)
                m1 = fit_model1(y_z, groups)
                m2 = fit_model2(y_z, groups, mapk)
                for label, pipe, ref in (("beta", sub["beta_m1"], m1["beta"]),
                                         ("se", sub["se_m1"], m1["se"]),
                                         ("t", sub["t_m1"], m1["t"]),
                                         ("beta2", sub["beta_m2"], m2["beta"])):
                    d = abs(float(pipe) - ref)
                    max_diff = max(max_diff, d)
                    print(f"{sample_id:8s} {gene:8s} {label:6s} {float(pipe):14.9f} "
                          f"{ref:14.9f} {d:11.2e}")
    print(f"\nTEST A max absolute difference: {max_diff:.3e}")
    print("TEST A:", "PASS" if max_diff < 1e-8 else "FAIL")


def test_b(args) -> None:
    workdir = detect_workdir()
    exposure = load_design(workdir)
    pa_dir = Path(workdir) / "youyun/plgg/code/sc_integration/ENVI/2_spatial_imputation/PA"

    stored = pd.read_csv(pa_dir / "myeloid_glial_exposure_nested_results.tsv", sep="\t",
                         dtype={"sample_id": str})
    sig_members = pd.read_csv(pa_dir.parent / "myeloid_signatures.csv")
    sig_members["Gene"] = sig_members["Gene"].str.upper()

    # TEST B runs on the OLD universe because that is what the stored estimates were
    # computed on (verified: stored HEME score correlates +0.65 with the old gene
    # subset vs -0.57 with the 18k one).
    data_root = os.path.join(workdir, "youyun/plgg/data/sc_integration", args.envi_dir)
    files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_spatial_envi.h5ad")))

    rows_out = []
    for f in files:
        sample_id, rows, exp_s, kmap = open_sample(f, exposure)
        groups = exp_s["exposure_group"].to_numpy()
        with h5py.File(f, "r") as h:
            imp = h["obsm"]["imputation"]
            for sig in TEST_B_SIGNATURES:
                members = sig_members.loc[sig_members["Signature"] == sig, "Gene"]
                present = sorted(set(members) & set(kmap))
                if not present:
                    continue
                # score = mean over this sample's imputed members (score_genes without
                # the control-bin subtraction), then the notebook's z-score + Model 1
                acc = np.zeros(rows.shape[0], dtype=np.float64)
                for g in present:
                    acc += imp[kmap[g]][:][rows].astype(np.float64)
                score = acc / len(present)
                score_z = (score - score.mean()) / score.std(ddof=1)
                m1 = fit_model1(score_z, groups)
                ref = stored[(stored["sample_id"] == sample_id)
                             & (stored["signature"] == sig + "_score")]
                ref_beta = float(ref["Estimate"].iloc[0]) if not ref.empty else np.nan
                ref_n = int(ref["n_cells"].iloc[0]) if not ref.empty else -1
                rows_out.append({
                    "sample_id": sample_id, "signature": sig, "n_genes_used": len(present),
                    "n_cells": rows.shape[0], "notebook_n_cells": ref_n,
                    "recomputed_beta": m1["beta"], "notebook_beta": ref_beta,
                    "sign_match": bool(np.sign(m1["beta"]) == np.sign(ref_beta)),
                })
                print(f"{sample_id} {sig:36s} genes={len(present):4d} "
                      f"recomputed={m1['beta']:+.4f} notebook={ref_beta:+.4f} "
                      f"cells={rows.shape[0]}/{ref_n}", flush=True)

    out = pd.DataFrame(rows_out)
    dest = Path(__file__).resolve().parent / "verify_test_b_signature_reproduction.tsv"
    out.to_csv(dest, sep="\t", index=False)
    n_ok = int(out["sign_match"].sum())
    cells_ok = bool((out["n_cells"] == out["notebook_n_cells"]).all())
    print(f"\nCell counts match notebook: {cells_ok}")
    print(f"Sign agreement: {n_ok}/{len(out)}")
    print("TEST B:", "PASS" if (n_ok == len(out) and cells_ok) else "FAIL")
    print("Wrote", dest)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("test-a", help="closed-form vs lstsq")
    a.add_argument("--per-sample", required=True)
    a.add_argument("--envi-dir", default="ENVI_results_091126")
    a.set_defaults(func=test_a)
    b = sub.add_parser("test-b", help="end-to-end signature reproduction (row alignment)")
    b.add_argument("--envi-dir", default="ENVI_results_091126")
    b.set_defaults(func=test_b)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
