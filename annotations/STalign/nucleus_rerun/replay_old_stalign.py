#!/usr/bin/env python3
# Replay the ORIGINAL STAlign notebook pipeline (annotations/STalign/pathology_<sample>_proseg.qmd) with
# different landmarks, to test what limited the old 030226 map.
#
# The pipeline is copied from the notebook cells: whole-slide JPG (normalised RGB, 1.76 um/px) as source;
# a 10 um cell-density raster of the ProSeg centroids (row = -x_centroid, col = -y_centroid) as target;
# L, T from the landmarks; LDDMM with niter 2000, sigmaM 0.15, sigmaB 0.10, sigmaA 0.11, epV 10,
# muB black, muA white, everything else STalign's default (a = 500 JPG px, sigmaP = 20 um, float64), landmarks
# also passed to LDDMM. Only the landmarks change:
#   original  : the notebook's hand-placed pointsI / pointsJ (should reproduce 030226)
#   true      : same Xenium-side points (pointsJ); pointsI moved to where the nucleus-verified new map
#               puts them on the H&E (exact correspondences)
#   true_grid : the `true` pointsI rounded to the 50 JPG px grid the notebooks used (tests the rounding)
# Step sizes: the notebooks set epV = 10 and leave epL / epT at STalign's defaults (2e-8 / 2e-1; LDDMM divides both
# by 10 from iteration diffeo_start = 0 on). --epL / --epT / --epV override them (e.g. --epL 0 --epT 0 = affine held
# at the landmark fit; 10x smaller steps = --epL 2e-9 --epT 2e-2 --epV 1) to test what moves the map. The affine
# (L, T) of every iteration is recorded by wrapping STalign.to_A, which LDDMM calls once per iteration.
# Pixel re-weighting: LDDMM re-estimates its per-pixel weights (matched / artifact / background, WM / WA / WB) every 5
# iterations from a hard-coded iteration 50. --reweight-start n recompiles LDDMM with n in its place; n > niter switches
# the re-weighting off (every pixel keeps weight 0.5), to test whether the drift that starts at iteration 50 comes from it.
# Output (NEW_DIR/qc/replay_old/): replay_<sample>_<landmarks>.npz with every cell's H&E (row, col) in JPG
# px (cells.csv.gz order, the 030226 row order), the landmarks used, the energy history and the affine per iteration;
# a one-line summary is appended to replay_summary.tsv: median / p90 distance (um) of each cell from
# the new map (= error) and from the 030226 map (replication check for `original`).
# Usage: python replay_old_stalign.py <sample> {original,true,true_grid} [--niter 2000] [--threads 2] [--tag x]
#        [--epL e --epT e --epV e] [--reweight-start n]
#        (--tag only renames the output, e.g. original_1thread for a thread-count check)
#        python replay_old_stalign.py --figure   -> replay_old/fig_replay_old.png + replay_start_end.tsv: for every
#        finished run, the error of the starting affine (L, T from its landmarks) and of the final map
#        python replay_old_stalign.py --trajectory -> replay_old/fig_replay_trajectory.png + replay_trajectory.tsv: for
#        runs from exact landmarks that recorded their affine, per iteration: how far the landmarks have moved (from
#        the landmark energy EP, sigmaP = 20 um in the Xenium frame), how far the affine alone has moved the section centre, and EM
#        python replay_old_stalign.py --conditions -> replay_old/fig_replay_conditions.png + replay_conditions.tsv: for
#        CL42112 and PA320784, the final median cell error from each start (hand / exact landmarks) under notebook
#        steps, 10x smaller steps, and 10x smaller steps with the pixel re-weighting off
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


def set_reweight_start(start):
    """Recompile STalign.LDDMM in its own module namespace with the weight re-estimation starting at `start` instead of
    the hard-coded 50 (the only change to its source)."""
    import inspect
    src = inspect.getsource(STalign.LDDMM)
    assert src.count("if it >= 50:") == 1, "STalign.LDDMM changed; re-check the patch"
    exec(compile(src.replace("if it >= 50:", f"if it >= {start}:"), STalign.__file__, "exec"), vars(STalign))


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


def trajectory():
    """Runs from exact landmarks that recorded their affine, per iteration:
    - the RMS landmark displacement of the full map, from EP = sum |d|^2 / (2 sigmaP^2). LDDMM compares A(pointsI) with
      pointsJ in the target frame, so d and sigmaP = 20 are Xenium um (checked: at iteration 0 this equals the RMS
      residual of the landmark affine). The landmarks start exact, so this is the distance moved from the truth there;
    - the section centre's displacement under the affine alone (A^-1 applied to the cells' centre in J);
    - the image term EM relative to iteration 0 (weights fixed until iteration 50)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": "#52514e", "xtick.color": "#52514e",
                         "ytick.color": "#52514e", "axes.spines.top": False, "axes.spines.right": False})
    rows, runs = [], {}
    for f in sorted(OUT.glob("replay_*_true*.npz")):
        s, name = f.stem.split("_", 2)[1:]
        if name == "true_grid":
            continue
        z = np.load(f)
        E = np.array(z["energies"])
        lm = np.sqrt(2 * 20.0 ** 2 * E[:, 3] / len(z["pointsI"]))       # Xenium um (no JPG px conversion)
        aff = np.full(len(lm), np.nan)                                 # runs before the affine was recorded: EP only
        if "affine_traj" in z:
            c = folder(s)
            cells = pd.read_csv(DATA / f"segmentation/proseg_run_121024/{c}/{c}_proseg/outs/cells.csv.gz", usecols=["x_centroid", "y_centroid"])
            ctr = np.r_[-cells.x_centroid.mean(), -cells.y_centroid.mean()]      # J frame: row = -x, col = -y
            src = np.array([np.linalg.solve(a[:4].reshape(2, 2), ctr - a[4:]) for a in z["affine_traj"]])
            aff = np.hypot(*(src - src[0]).T) * JPG_UM
        runs[(s, name)] = (lm, aff)
        rows += [dict(sample=s, run=name, iteration=i, landmark_rms_um=lm[i], affine_centre_shift_um=aff[i], EM_rel=E[i, 1] / E[0, 1])
                 for i in range(len(lm))]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "replay_trajectory.tsv", sep="\t", index=False)
    # per run: how much the landmarks jump per iteration, and where they are over the last 1000 iterations
    t["jump_um"] = t.groupby(["sample", "run"]).landmark_rms_um.diff().abs()
    late = t[t.iteration >= 1000]
    print(pd.concat([t.groupby(["sample", "run"]).jump_um.median().rename("median_jump_per_iteration_um"),
                     t[t.iteration <= 10].groupby(["sample", "run"]).EM_rel.max().rename("EM_rel_max_first_10_iterations"),
                     late.groupby(["sample", "run"]).landmark_rms_um.median().rename("landmark_rms_median_last_1000_um"),
                     late.groupby(["sample", "run"]).landmark_rms_um.max().rename("landmark_rms_max_last_1000_um")], axis=1).round(2).to_string())
    samples = sorted({s for s, _ in runs})
    # notebook settings: the 1-thread rerun where it exists, else the original 2-thread replay ("true")
    has_1thread = {s for s, n in runs if n == "true_traj_1thread"}
    style = {"true_traj_1thread": ("#eb6834", "notebook step sizes"), "true": ("#eb6834", "notebook step sizes"),
             "true_frozen_affine_1thread": ("#e87ba4", "affine frozen, notebook epV"),
             "true_smallstep_1thread": ("#2a78d6", "all steps 10x smaller"),
             "true_smallstep_noreweight_1thread": ("#1baf7a", "10x smaller, re-weighting off")}
    import matplotlib.ticker
    fig, ax = plt.subplots(1, len(samples), figsize=(4.4 * len(samples), 3.9), sharey=True, squeeze=False)
    for a, s in zip(ax[0], samples):
        for name, (col, lab) in style.items():
            if (s, name) in runs and not (name == "true" and s in has_1thread):
                a.plot(runs[(s, name)][0], color=col, lw=0.9, label=lab)
        a.axvline(50, color="#8a8984", lw=0.8, ls=(0, (3, 3)))
        a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        a.set_yscale("log"); a.set_ylim(0.1, 1000); a.set_xlim(0, 2000)
        a.set_title(f"{s}", loc="left", fontsize=9.5); a.set_xlabel("LDDMM iteration")
    ax[0, 0].set_ylabel("landmarks' distance from their true position\n(RMS of 4-5 landmarks, µm, log)")
    ax[0, 0].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.suptitle("Old STAlign pipeline (JPG vs 10 µm cell-count map), started from exact landmarks\n"
                 "dashed line: STalign starts re-estimating its pixel weights", x=0.01, ha="left", fontsize=9.5)
    fig.tight_layout(); fig.savefig(OUT / "fig_replay_trajectory.png", dpi=150); plt.close(fig)


def conditions():
    """Final median cell error (vs the new map) per start and optimizer setting, from replay_summary.tsv; the start is
    the landmark affine's error (as in figure()). Notebook steps were run on 1 and 2 threads where both exist."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = pd.read_csv(OUT / "replay_summary.tsv", sep="\t")
    setting = {"": "notebook steps", "_1thread": "notebook steps", "_traj_1thread": "notebook steps",
               "_smallstep_1thread": "10x smaller steps", "_smallstep_noreweight_1thread": "10x smaller steps, re-weighting off"}
    rows = []
    for s in ("CL42112", "PA320784"):
        c = folder(s)
        cells = pd.read_csv(DATA / f"segmentation/proseg_run_121024/{c}/{c}_proseg/outs/cells.csv.gz", usecols=["cell_id", "x_centroid", "y_centroid"])
        new = pd.read_csv(meta_path(s, NEW_DIR), usecols=["cell_id", "STAlign_x", "STAlign_y"]).set_index("cell_id").loc[cells.cell_id].values
        for start in ("original", "true"):
            z = np.load(OUT / f"replay_{s}_{start}.npz")
            L, T = STalign.L_T_from_points(z["pointsI"], z["pointsJ"])
            st = np.linalg.solve(L, np.vstack([-cells.x_centroid.values - T[0], -cells.y_centroid.values - T[1]])).T
            for suffix, lab in setting.items():
                r = t[(t["sample"] == s) & (t.landmarks == start + suffix)]
                if len(r):
                    rows.append(dict(sample=s, start=start, start_um=np.nanmedian(np.hypot(*(st - new).T)) * JPG_UM,
                                     setting=lab, run=start + suffix, end_um=r.err_vs_new_median_um.iloc[-1]))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "replay_conditions.tsv", sep="\t", index=False)
    print(d.round(1).to_string(index=False))
    col = {"notebook steps": "#eb6834", "10x smaller steps": "#2a78d6", "10x smaller steps, re-weighting off": "#1baf7a"}
    lab = {"original": "hand landmarks", "true": "exact landmarks"}
    import matplotlib.ticker
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    ticks = []
    for i, ((s, start), g) in enumerate(d.groupby(["sample", "start"], sort=False)):
        y = i * 1.25                                                   # sub-rows: start, then one per setting
        ax.scatter(g.start_um.iloc[0], y - 0.36, s=40, facecolor="white", edgecolor="#8a8984", lw=1.2, zorder=3)
        ax.text(g.start_um.iloc[0] * 1.15, y - 0.36, f"{g.start_um.iloc[0]:.1f}", va="center", fontsize=7.5, color="#52514e")
        for k, (setname, h) in enumerate(g.groupby("setting", sort=False)):
            yy = y + (k - 0.5) * 0.24
            ax.plot([h.end_um.min(), h.end_um.max()], [yy, yy], color=col[setname], lw=1.2)
            ax.scatter(h.end_um, [yy] * len(h), s=34, color=col[setname], edgecolor="white", lw=1, zorder=4,
                       label=setname if i == 0 else None)
            ax.text(h.end_um.max() * 1.15, yy, " / ".join(f"{v:.1f}" for v in h.end_um), va="center", fontsize=7.5, color="#52514e")
        ticks.append((y, f"{s}\n{lab[start]}"))
    ax.scatter([], [], s=40, facecolor="white", edgecolor="#8a8984", lw=1.2, label="start (landmark affine)")
    ax.set_yticks([p for p, _ in ticks], [q for _, q in ticks], fontsize=8); ax.invert_yaxis()
    ax.set_xscale("log"); ax.set_xlim(0.1, 600)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("median cell error vs the nucleus-verified map (µm, log)")
    ax.grid(axis="x", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Old STAlign pipeline: where 2000 iterations end, by optimizer setting\n"
                 "same old images and landmarks; notebook steps on 1 and 2 CPU threads", loc="left", fontsize=9.5)
    fig.tight_layout(); fig.savefig(OUT / "fig_replay_conditions.png", dpi=150); plt.close(fig)


def main():
    if sys.argv[1:] == ["--conditions"]:
        return conditions()
    if sys.argv[1:] == ["--figure"]:
        return figure()
    if sys.argv[1:] == ["--trajectory"]:
        return trajectory()
    ap = argparse.ArgumentParser()
    ap.add_argument("sample"); ap.add_argument("landmarks", choices=["original", "true", "true_grid"])
    ap.add_argument("--niter", type=int, default=2000); ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--tag", default="", help="suffix for the output name only")
    ap.add_argument("--epL", type=float, default=None); ap.add_argument("--epT", type=float, default=None)
    ap.add_argument("--epV", type=float, default=10.0, help="the notebooks' 10 unless overridden")
    ap.add_argument("--reweight-start", type=int, default=50, help="STalign's 50 unless overridden; > niter = off")
    args = ap.parse_args(); s = args.sample
    name = args.landmarks + (f"_{args.tag}" if args.tag else "")
    torch.set_num_threads(args.threads); torch.set_default_device("cpu")
    if args.reweight_start != 50:
        set_reweight_start(args.reweight_start)
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
                  sigmaA=0.11, epV=args.epV, muB=torch.tensor([0, 0, 0]), muA=torch.tensor([1, 1, 1]))
    params.update({k: getattr(args, k) for k in ("epL", "epT") if getattr(args, k) is not None})
    traj, to_A = [], STalign.to_A
    def to_A_logged(L_, T_):                                          # record the affine of each iteration
        traj.append(np.r_[L_.detach().numpy().ravel(), T_.detach().numpy()])
        return to_A(L_, T_)
    STalign.to_A = to_A_logged
    t0 = time.time()
    with no_lddmm_figures():
        out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **params)
    STalign.to_A = to_A
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
                        pointsJ=pJ, energies=np.array(_Sink.energies), A=out["A"].cpu().numpy(), affine_traj=np.array(traj))
    row = dict(sample=s, landmarks=name, niter=args.niter, epL=params.get("epL", "default"), epT=params.get("epT", "default"), epV=args.epV,
               reweight_start=args.reweight_start, seconds=round(secs),
               landmark_err_median_um=np.median(lm_err), err_vs_new_median_um=np.nanmedian(err),
               err_vs_new_p90_um=np.nanquantile(err, 0.9), dist_to_030226_median_um=np.nanmedian(rep),
               dist_to_030226_p90_um=np.nanquantile(rep, 0.9))
    summ = OUT / "replay_summary.tsv"
    t = pd.DataFrame([row])
    if summ.exists():                                                 # rows written before these were recorded: the defaults
        t = pd.concat([pd.read_csv(summ, sep="\t"), t]).fillna({"epL": "default", "epT": "default", "epV": 10.0, "reweight_start": 50})
    t.to_csv(summ, sep="\t", index=False)
    print(pd.DataFrame([row]).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
