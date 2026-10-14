#!/usr/bin/env python3
# Replay the ORIGINAL STAlign notebook pipeline (annotations/STalign/pathology_<sample>_proseg.qmd) with
# different landmarks, to test what limited the old 030226 map.
#
# The pipeline is copied from the notebook cells: whole-slide JPG (normalised RGB, 1.76 um/px) as source;
# a 10 um cell-density raster of the ProSeg centroids (row = -x_centroid, col = -y_centroid) as target;
# L, T from the landmarks; LDDMM with niter 2000, sigmaM 0.15, sigmaB 0.10, sigmaA 0.11, epV 10,
# muB black, muA white, everything else STalign's default (a = 500, sigmaP = 20, float64), landmarks
# also passed to LDDMM. Only the landmarks change:
#   original  : the notebook's hand-placed pointsI / pointsJ (should reproduce 030226)
#   true      : same Xenium-side points (pointsJ); pointsI moved to where the nucleus-verified new map
#               puts them on the H&E (exact correspondences)
#   true_grid : the `true` pointsI rounded to the 50 JPG px grid the notebooks used (tests the rounding)
# Output (NEW_DIR/qc/replay_old/): replay_<sample>_<landmarks>.npz with every cell's H&E (row, col) in JPG
# px (cells.csv.gz order, the 030226 row order), the landmarks used, and the energy history;
# a one-line summary is appended to replay_summary.tsv: median / p90 distance (um) of each cell from
# the new map (= error) and from the 030226 map (replication check for `original`).
# Usage: python replay_old_stalign.py <sample> {original,true,true_grid} [--niter 2000] [--threads 2] [--tag x]
#        (--tag only renames the output, e.g. original_1thread for a thread-count check)
#        python replay_old_stalign.py --figure   -> replay_old/fig_replay_old.png + replay_start_end.tsv: for every
#        finished run, the error of the starting affine (L, T from its landmarks) and of the final map
import argparse, ast, re, sys, time
from pathlib import Path
import numpy as np, pandas as pd, torch
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stalign_nucleus import NEW_DIR, NOTEBOOKS, no_lddmm_figures, _Sink
from nucleus_windows import GROUP, DATA, meta_path, folder, MPP0, JPG_TO_L0
from STalign import STalign

Image.MAX_IMAGE_PIXELS = None
JPG_UM = JPG_TO_L0 * MPP0
OUT = NEW_DIR / "qc/replay_old"


def notebook_landmarks(s):
    """pointsI (H&E JPG px, row/col) and pointsJ (um, row = -x_centroid, col = -y_centroid) from the qmd."""
    txt = (NOTEBOOKS / f"pathology_{s}_proseg.qmd").read_text()
    get = lambda k: np.array(ast.literal_eval(re.search(k + r" = np\.array\((\[.*?\])\)", txt, re.S).group(1)), float)
    return get("pointsI"), get("pointsJ")


def new_map_jpg(s, xy):
    """Xenium (x, y) um -> H&E (row, col) JPG px through the new meta: piecewise-linear over the cells, with a
    global least-squares affine for points outside the cells' hull (landmarks can sit on the tissue edge)."""
    from scipy.interpolate import LinearNDInterpolator
    m = pd.read_csv(meta_path(s, NEW_DIR), usecols=["x_centroid", "y_centroid", "STAlign_x", "STAlign_y"]).dropna()
    src, dst = m[["x_centroid", "y_centroid"]].values, m[["STAlign_x", "STAlign_y"]].values
    out = LinearNDInterpolator(src, dst)(xy)
    M, *_ = np.linalg.lstsq(np.c_[src, np.ones(len(src))], dst, rcond=None)
    bad = np.isnan(out).any(1)
    out[bad] = np.c_[xy[bad], np.ones(bad.sum())] @ M
    return out


def figure():
    """Start (landmark affine) -> end (after LDDMM) median cell error for every finished replay, vs 030226."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for f in sorted(OUT.glob("replay_*.npz")):
        s, name = f.stem.split("_", 2)[1:]
        z = np.load(f)
        c = folder(s)
        cells = pd.read_csv(DATA / f"segmentation/proseg_run_121024/{c}/{c}_proseg/outs/cells.csv.gz",
                            usecols=["cell_id", "x_centroid", "y_centroid"])
        new = pd.read_csv(meta_path(s, NEW_DIR), usecols=["cell_id", "STAlign_x", "STAlign_y"]).set_index("cell_id").loc[cells.cell_id].values
        old = pd.read_csv(meta_path(s), usecols=["cell_id", "STAlign_x", "STAlign_y"]).set_index("cell_id").loc[cells.cell_id].values
        L, T = STalign.L_T_from_points(z["pointsI"], z["pointsJ"])
        # the notebook's affine-only placement: inv(L) @ ([row, col] - T), with row = -x, col = -y
        start = np.linalg.solve(L, np.vstack([-cells.x_centroid.values - T[0], -cells.y_centroid.values - T[1]])).T
        med = lambda a: np.nanmedian(np.hypot(*(a - new).T)) * JPG_UM
        rows.append(dict(sample=s, landmarks=name, start_um=med(start), end_um=med(z["rowcol_jpg"]), map_030226_um=med(old)))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "replay_start_end.tsv", sep="\t", index=False)
    lab = {"original": "hand landmarks", "original_1thread": "hand landmarks, 1 thread", "true": "exact landmarks",
           "true_grid": "exact, rounded to 50 grid"}
    col = {"original": "#d1242f", "original_1thread": "#d1242f", "true": "#2da44e", "true_grid": "#e5a50a"}
    order = ["true", "true_grid", "original", "original_1thread"]
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(d) + 1.8))
    y, ticks = 0, []
    for s in d["sample"].unique():
        x = d[d["sample"] == s].set_index("landmarks")
        for k in [o for o in order if o in x.index]:
            r = x.loc[k]
            ax.annotate("", xy=(r.end_um, y), xytext=(r.start_um, y), arrowprops=dict(arrowstyle="->", color=col[k], lw=1.6))
            ax.scatter(r.start_um, y, s=22, color=col[k], zorder=3)
            ax.text(max(r.start_um, r.end_um) * 1.12, y, f"{r.start_um:.1f} -> {r.end_um:.1f} um", va="center", fontsize=7.5, color="#52514e")
            ticks.append((y, f"{s}: {lab[k]}")); y += 1
        ax.scatter(x.map_030226_um.iloc[0], y, marker="D", s=26, color="#1f5fae", zorder=3)
        ax.text(x.map_030226_um.iloc[0] * 1.12, y, f"{x.map_030226_um.iloc[0]:.1f} um", va="center", fontsize=7.5, color="#52514e")
        ticks.append((y, f"{s}: stored 030226 map")); y += 1.6
    ax.set_yticks([p for p, _ in ticks], [q for _, q in ticks], fontsize=8); ax.invert_yaxis()
    ax.set_xscale("log"); ax.set_xlim(0.1, 2000)
    ax.set_xlabel("median cell error vs the nucleus-verified map (um, log)")
    ax.set_title("Replaying the original STAlign notebooks with different landmarks\n"
                 "dot = starting affine from the landmarks; arrow head = after the old LDDMM", loc="left", fontsize=10)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "fig_replay_old.png", dpi=150); plt.close(fig)
    print(d.round(1).to_string(index=False))


def main():
    if sys.argv[1:] == ["--figure"]:
        return figure()
    ap = argparse.ArgumentParser()
    ap.add_argument("sample"); ap.add_argument("landmarks", choices=["original", "true", "true_grid"])
    ap.add_argument("--niter", type=int, default=2000); ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--tag", default="", help="suffix for the output name only")
    args = ap.parse_args(); s = args.sample
    name = args.landmarks + (f"_{args.tag}" if args.tag else "")
    torch.set_num_threads(args.threads); torch.set_default_device("cpu")
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- notebook cells, verbatim in substance
    V = np.asarray(Image.open(DATA / f"Xenium_annotations/images/230918_Xenium_CytAssist_{GROUP[s]}.jpg").convert("RGB"))
    I = STalign.normalize(V).transpose(2, 0, 1)
    YI, XI = np.arange(I.shape[1]) * 1., np.arange(I.shape[2]) * 1.
    f = folder(s)
    cells = pd.read_csv(DATA / f"segmentation/proseg_run_121024/{f}/{f}_proseg/outs/cells.csv.gz")
    xM, yM = -cells.y_centroid.values, -cells.x_centroid.values
    with no_lddmm_figures():
        XJ, YJ, M, _ = STalign.rasterize(xM, yM, dx=10)
    J = STalign.normalize(np.vstack((M, M, M)))

    # ---- landmarks for this condition
    pI, pJ = notebook_landmarks(s)
    if args.landmarks != "original":
        pI = new_map_jpg(s, np.c_[-pJ[:, 0], -pJ[:, 1]])            # true H&E position of each pointsJ
        if args.landmarks == "true_grid":
            pI = np.round(pI / 50) * 50
    L, T = STalign.L_T_from_points(pI, pJ)
    params = dict(L=L, T=T, niter=args.niter, pointsI=pI, pointsJ=pJ, device="cpu", sigmaM=0.15, sigmaB=0.10,
                  sigmaA=0.11, epV=10, muB=torch.tensor([0, 0, 0]), muA=torch.tensor([1, 1, 1]))
    t0 = time.time()
    with no_lddmm_figures():
        out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **params)
    secs = time.time() - t0
    tp = STalign.transform_points_target_to_source(out["xv"], out["v"], out["A"], np.stack([yM, xM], -1))
    tp = tp.detach().cpu().numpy()                                    # (row, col) JPG px, cells.csv order

    # ---- score: distance from the new map (error) and from 030226 (replication of `original`)
    new = pd.read_csv(meta_path(s, NEW_DIR), usecols=["cell_id", "STAlign_x", "STAlign_y"]).set_index("cell_id")
    old = pd.read_csv(meta_path(s), usecols=["cell_id", "STAlign_x", "STAlign_y"]).set_index("cell_id")
    ids = cells.cell_id.values
    err = np.hypot(*(tp - new.loc[ids].values).T) * JPG_UM
    rep = np.hypot(*(tp - old.loc[ids].values).T) * JPG_UM
    lm_err = np.hypot(*(pI - new_map_jpg(s, np.c_[-pJ[:, 0], -pJ[:, 1]])).T) * JPG_UM
    np.savez_compressed(OUT / f"replay_{s}_{name}.npz", rowcol_jpg=tp.astype(np.float32), pointsI=pI,
                        pointsJ=pJ, energies=np.array(_Sink.energies), A=out["A"].cpu().numpy())
    row = dict(sample=s, landmarks=name, niter=args.niter, seconds=round(secs),
               landmark_err_median_um=np.median(lm_err), err_vs_new_median_um=np.nanmedian(err),
               err_vs_new_p90_um=np.nanquantile(err, 0.9), dist_to_030226_median_um=np.nanmedian(rep),
               dist_to_030226_p90_um=np.nanquantile(rep, 0.9))
    summ = OUT / "replay_summary.tsv"
    pd.DataFrame([row]).to_csv(summ, sep="\t", index=False, mode="a", header=not summ.exists())
    print(pd.DataFrame([row]).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
