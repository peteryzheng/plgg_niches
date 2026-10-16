#!/usr/bin/env python3
# Nucleus-guided STAlign rerun for one sample (H&E <-> Xenium), replacing the 030226 registration.
#
# Why    : the 030226 map (whole-slide JPG vs a 10 um cell-density raster, a = 500, 4-6 hand landmarks)
#          is off by 16-75 um in 10/11 samples (README.md). Here STAlign itself is rerun with images
#          that show nuclei on both sides and landmarks measured from nucleus patterns.
# Steps  : 1. windows : 150 um windows across the section; the nucleus check on the OLD map gives each
#                       a shift (nucleus_windows.py). Confident windows (z >= Z_CONF) are split in half
#                       at random: one half become landmarks, the other half is held out for scoring.
#          2. images  : I = H&E hematoxylin (SVS, not the JPG), J = Xenium DAPI morphology_focus, both in
#                       um so that STAlign's smoothness `a` is in um. J is put in H&E orientation
#                       (J row = Xenium x, J col = -Xenium y; checked per sample from the old map).
#          3. affine  : STalign.L_T_from_points on the landmark half (iterative outlier rejection).
#          4. LDDMM   : coarse (~4 um/px, a = 200 um) then fine (~2 um/px, a = 150 um), landmark term kept
#                       tight (sigmaP 0.5 um) so the weak image term cannot pull the map off the nuclei.
#                       Kept for consistency with the original STAlign method (user decision, 2026-10-14).
#                       Pilot note: on CL42112 the landmark affine alone scored 0.71 um on held-out
#                       windows vs 1.06 um with LDDMM -- same slide, so the true map is near-affine;
#                       --niter1 0 --niter2 0 gives the affine-only map (archived in qc/trial_affine_all/).
#          5. output  : every cell of the old meta mapped into the H&E (STAlign_x = row, STAlign_y = col,
#                       JPG px = level-0 / 7) and relabelled with the sample notebook's own polygon /
#                       recode logic; written to NEW_DIR. The 030226 folder is only read.
#          6. validate: the nucleus check on the NEW map at the same window centres (before vs after,
#                       landmark vs held-out windows reported separately).
# Usage  : python stalign_nucleus.py CL42112 [--threads 4]                    (affine + LDDMM; conda env spatial)
#          python stalign_nucleus.py CL42112 --niter1 0 --niter2 0 --tag affine  (affine-only variant)
#          python stalign_nucleus.py CL42112 --out-dir <empty dir>  (from-scratch rerun; nothing cached is reused)
# Outputs: NEW_DIR/proseg_pathology_<s>_proseg_meta.csv; NEW_DIR/qc/{windows_before,windows_after}_<s>.tsv,
#          transform_<s>.npz, summary_<s>.tsv, run_<s>.log
import argparse, ast, json, re, sys, time
from contextlib import contextmanager
from pathlib import Path
import numpy as np, pandas as pd, tifffile, torch
import matplotlib
matplotlib.use("Agg")
import geopandas as gpd
from skimage import color
from STalign import STalign

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from nucleus_windows import (Mapping, draw_windows, measure, meta_path, svs_for, dapi_path, DAPI, SVS, MPP0,
                             JPG_TO_L0, DATA, GROUP, DAPI_UM, Z_CONF)

NOTEBOOKS = HERE.parent                      # annotations/STalign/pathology_<s>_proseg.qmd
NEW_DIR = DATA / "Xenium_annotations/pathology_annotations_101326_nucleus"
SVS_PAGE, DAPI_LEVEL = 2, 2                  # SVS page 2 = level 0 / 4 (~1.01 um/px); DAPI level 2 = 0.85 um/px
N_WIN, WIN_SPACING_UM = 100, 250.0           # window budget per sample and minimum centre spacing


# ---------------------------------------------------------------- images
def block_mean(a, b):
    """Integer binning (mean) of a 2D image; trailing rows/cols that do not fill a block are dropped."""
    h, w = a.shape[0] // b * b, a.shape[1] // b * b
    return a[:h, :w].reshape(h // b, b, w // b, b).mean((1, 3))


def scale01(a, lo_q=0.5, hi_q=99.8):
    """Robust 0-1 scaling (percentile clip) so the two modalities have comparable contrast ranges;
    STAlign additionally fits a linear contrast map from I to J."""
    lo, hi = np.percentile(a, [lo_q, hi_q])
    return np.clip((a - lo) / (hi - lo), 0, 1)


def he_image(s, mp, margin_um):
    """Hematoxylin channel (nuclei bright) of this sample's section at ~1 um/px. The crop is the old
    map's extent of the sample's cells plus a margin larger than the old error (< 150 um), so the
    whole Xenium field lands inside I. Returns (image, row centres um, col centres um) where um =
    level-0 px * MPP0."""
    svs = SVS(DATA / f"Xenium_annotations/images/230918_Xenium_CytAssist_{GROUP[s]}.svs", page=SVS_PAGE)
    lo = np.floor((mp.l0_rc.min(0) - margin_um / MPP0) / svs.down).astype(int)
    hi = np.ceil((mp.l0_rc.max(0) + margin_um / MPP0) / svs.down).astype(int)
    rgb = svs.read(lo[0], lo[1], hi[0] - lo[0], hi[1] - lo[1])
    hem = color.rgb2hed(rgb)[..., 0].astype(np.float32)
    return hem, lo * svs.down * MPP0, svs.down * MPP0              # origin (edge, um) and px size (um)


def dapi_image(s, xy_um, margin_um):
    """DAPI morphology_focus at 0.85 um/px cropped to the cells' extent plus a margin.
    Returns (image[row = y, col = x], origin (y, x) um, px size um)."""
    lev = tifffile.TiffFile(dapi_path(s)).series[0].levels[DAPI_LEVEL]
    px = DAPI_UM * 2 ** DAPI_LEVEL
    a = lev.asarray(maxworkers=4)
    lo = np.floor((xy_um.min(0)[::-1] - margin_um) / px).clip(0).astype(int)        # (row, col)
    hi = np.ceil((xy_um.max(0)[::-1] + margin_um) / px).astype(int)
    return a[lo[0]:hi[0], lo[1]:hi[1]].astype(np.float32), lo * px, px


def binned(img, origin, px, b):
    """Bin by b and return the image with pixel-centre coordinates (um) along rows and cols."""
    im = block_mean(img, b)
    return im, [origin[i] + (np.arange(im.shape[i]) + 0.5) * px * b for i in range(2)]


def to_j_frame(img, coords):
    """DAPI (row = y, col = x) -> J frame (row = x, col = -y): H&E rows follow Xenium x and H&E cols
    run against Xenium y, so this makes the initial affine close to a pure translation."""
    y, x = coords
    return img.T[:, ::-1].copy(), [x, -y[::-1]]


# ---------------------------------------------------------------- landmarks
def landmark_pairs(t, mp):
    """Each confident window gives a pair: H&E window centre (um) <-> the Xenium point (J frame) that
    truly sits under it = old-map preimage of (centre + measured shift)."""
    c = t[["row_l0", "col_l0"]].values.astype(float)
    q = mp.to_um(c + t[["shift_row_um", "shift_col_um"]].values / MPP0)
    return c * MPP0, np.c_[q[:, 0], -q[:, 1]]


def robust_affine(pI, pJ, min_tol_um=15.0):
    """L_T_from_points with iterative rejection of pairs whose affine residual exceeds
    max(min_tol_um, median + 3 robust SD): a wrong FFT peak can pass the z threshold, and a bad
    landmark would pull the landmark term of LDDMM. Non-rigid tissue error up to ~min_tol_um is
    expected and kept."""
    keep = np.ones(len(pI), bool)
    for _ in range(10):
        L, T = STalign.L_T_from_points(pI[keep], pJ[keep])
        r = np.linalg.norm(pI @ L.T + T - pJ, axis=1)
        med = np.median(r[keep]); mad = 1.4826 * np.median(np.abs(r[keep] - med))
        new = r <= max(min_tol_um, med + 3 * mad)
        if (new == keep).all():
            break
        keep = new
    return L, T, keep, r


# ---------------------------------------------------------------- LDDMM
class _Sink:
    """Stand-in for the matplotlib figures STalign.LDDMM draws every 10 iterations (slow on full-size
    images, useless headless). It also keeps a reference to the energy history LDDMM passes to plot(),
    which LDDMM does not return."""
    energies = None
    def __getattr__(self, k): return self
    def __getitem__(self, k): return self
    def __call__(self, *a, **k): return self
    def plot(self, *a, **k):
        if a and isinstance(a[0], list) and a[0] and isinstance(a[0][0], list) and len(a[0][0]) >= 3:
            _Sink.energies = a[0]
        return self


@contextmanager
def no_lddmm_figures():
    orig = STalign.plt.subplots
    STalign.plt.subplots = lambda *a, **k: (_Sink(), _Sink())
    try:
        yield
    finally:
        STalign.plt.subplots = orig


def lddmm(xI, I, xJ, J, pI, pJ, init, log, **kw):
    """One STalign.LDDMM stage on single-channel images; `init` holds either L, T or A, v, xv."""
    t0 = time.time()
    with no_lddmm_figures():
        out = STalign.LDDMM(xI, I[None], xJ, J[None], pointsI=pI, pointsJ=pJ, device="cpu",
                            dtype=torch.float32, **init, **kw)
    E = np.array(_Sink.energies)
    log(f"    {kw.get('niter')} iterations in {time.time() - t0:.0f} s; energy E, EM, ER, EP: "
        f"start {np.round(E[0], 1).tolist()} -> end {np.round(E[-1], 1).tolist()}")
    return {k: (v.double() if torch.is_tensor(v) else [x.double() for x in v]) for k, v in out.items() if k in ("A", "v", "xv")}, E


def to_he(tr, pJ):
    """J-frame points (row = x, col = -y, um) -> H&E (row, col) um through the fitted map."""
    p = STalign.transform_points_target_to_source(tr["xv"], tr["v"], tr["A"], torch.as_tensor(pJ, dtype=torch.float64))
    return p.numpy()


# ---------------------------------------------------------------- labels
def notebook_spec(s):
    """GeoJSON path and colour->label recode map exactly as in this sample's STAlign notebook."""
    txt = (NOTEBOOKS / f"pathology_{s}_proseg.qmd").read_text()
    gj = re.search(r'GEOJSON_PATH = XENIUM_ANNOT_ROOT / "geojsons" / f?"([^"]+)"', txt).group(1)
    rm = ast.literal_eval(re.search(r"recode_map = (\{.*?\})", txt, re.S).group(1))
    return DATA / "Xenium_annotations/geojsons" / gj.replace("{IMAGE_GROUP}", GROUP[s]), rm


def label_cells(s, row_jpg, col_jpg):
    """Pathology label per cell, reproducing the notebook's geojson_contains_xy: first polygon (in
    file order) containing the point (col, row) in JPG px; its hex colour recoded by recode_map;
    region_id = the feature's `id`, else its row index; 'unannotated' / NA when no polygon hits."""
    gj, recode = notebook_spec(s)
    poly = gpd.read_file(gj)
    rgb = poly["classification"].apply(lambda t: json.loads(t).get("color"))
    poly["color"] = rgb.apply(lambda c: "#{:02x}{:02x}{:02x}".format(int(c[0]), int(c[1]), int(c[2])))
    poly.geometry = poly.geometry.scale(xfact=1 / 7, yfact=1 / 7, origin=(0, 0))
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(col_jpg, row_jpg), crs=poly.crs)   # same (pixel) frame
    hit = gpd.sjoin(pts, poly[["geometry"]], predicate="within", how="left")["index_right"]
    first = hit.groupby(level=0).min()                          # lowest feature index = iloc[0] in the notebook
    ann = pd.Series("unannotated", index=pts.index, dtype=object)
    rid = pd.Series(pd.NA, index=pts.index, dtype="string")
    has = first.notna()
    j = first[has].astype(int).values
    ann[has] = poly["color"].values[j]
    ids = poly["id"].values[j] if "id" in poly else np.full(len(j), None)
    rid[has] = [str(i) if i is not None and not pd.isna(i) else str(k) for i, k in zip(ids, j)]
    return ann.replace(recode), rid


def write_meta(s, old, row_jpg, col_jpg, path):
    """Old meta with the registration-dependent columns replaced (same columns, same order)."""
    new = old.copy()
    new["STAlign_x"], new["STAlign_y"] = row_jpg, col_jpg
    new["adjust_xy"] = [f"({c!r}, {r!r})" for c, r in zip(col_jpg.tolist(), row_jpg.tolist())]
    new["annotation"], new["region_id"] = label_cells(s, row_jpg, col_jpg)
    new.to_csv(path, index=False)
    return new


# ---------------------------------------------------------------- scoring
def summarise(t, label):
    """Shift distribution of confident windows (z >= Z_CONF) plus the share of H&E nuclei with a DAPI
    partner within MATCH_UM at zero shift (direct, threshold-free view of nucleus-level agreement)."""
    c = t[t.peak_z >= Z_CONF]
    return dict(set=label, n_windows=len(t), n_confident=len(c), median_um=c.shift_um.median(),
                p90_um=c.shift_um.quantile(0.9), match_at_zero=t.match_at_zero.median(),
                match_at_best=t.match_at_best.median())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample")
    ap.add_argument("--b1", type=int, default=4, help="H&E binning for stage 1 (x ~1.01 um)")
    ap.add_argument("--b2", type=int, default=2, help="H&E binning for stage 2")
    # LDDMM defaults from the CL42112 step-size scan: epV >= 10 drifts tens of um within 100 iterations,
    # and with sigmaM 0.2 / sigmaP 5 the (weak, ~4 um) image term pulls landmarks off; epL = epT = 0 keep
    # the landmark affine fixed so only the velocity field adapts. niter 0 = skip that stage.
    ap.add_argument("--niter1", type=int, default=150); ap.add_argument("--niter2", type=int, default=100)
    ap.add_argument("--a1", type=float, default=200.0); ap.add_argument("--a2", type=float, default=150.0)
    ap.add_argument("--epV", type=float, default=0.3); ap.add_argument("--epL", type=float, default=0.0)
    ap.add_argument("--epT", type=float, default=0.0)
    ap.add_argument("--sigmaM", type=float, default=0.5); ap.add_argument("--sigmaP", type=float, default=0.5)
    ap.add_argument("--sigmaR", type=float, default=5e5)
    ap.add_argument("--tag", default="", help="suffix for a trial run (outputs go to qc/trial_<tag>/)")
    # output root (default: the production folder). The old-map window table is cached in <out-dir>/qc, so an
    # empty out-dir reruns everything from the raw inputs without touching the production files.
    ap.add_argument("--out-dir", default=str(NEW_DIR))
    ap.add_argument("--threads", type=int, default=4, help="torch CPU threads (match the SGE slots per job)")
    args = ap.parse_args()
    s = args.sample
    torch.set_num_threads(args.threads)
    root = Path(args.out_dir)
    qc = root / "qc" / (f"trial_{args.tag}" if args.tag else ""); qc.mkdir(parents=True, exist_ok=True)
    out_meta = (qc if args.tag else root) / f"proseg_pathology_{s}_proseg_meta.csv"
    logf = open(qc / f"run_{s}.log", "a")
    def log(m):
        print(m, flush=True); logf.write(m + "\n"); logf.flush()
    log(f"== {s} {time.strftime('%F %T')} {vars(args)}")

    # 1. windows on the old map (cached: the same centres are reused for the after-measurement)
    mp_old = Mapping(meta_path(s))
    old = pd.read_csv(meta_path(s))
    wb = root / "qc" / f"windows_before_{s}.tsv"
    if not wb.exists():
        c = draw_windows(mp_old, N_WIN, WIN_SPACING_UM)
        measure(s, mp_old, c, svs_for(s), DAPI(dapi_path(s))).to_csv(wb, sep="\t", index=False)
    t = pd.read_csv(wb, sep="\t")
    # sanity: the frame convention (H&E row ~ +x, col ~ -y) must hold for this sample
    cr = np.corrcoef(np.c_[mp_old.um_xy, mp_old.l0_rc].T)
    log(f"  frame check: corr(x, row) {cr[0, 2]:+.3f}, corr(y, col) {cr[1, 3]:+.3f}")
    assert cr[0, 2] > 0.99 and cr[1, 3] < -0.99, "unexpected H&E/Xenium orientation for this sample"

    # split confident windows into landmark / held-out halves (fixed seed, recorded in the table)
    conf = np.flatnonzero((t.peak_z >= Z_CONF).values)
    rng = np.random.default_rng(1)
    fit_idx = np.sort(rng.permutation(conf)[: len(conf) // 2])
    t["role"] = np.where(t.index.isin(fit_idx), "landmark", np.where(t.peak_z >= Z_CONF, "heldout", "unconfident"))
    pI, pJ = landmark_pairs(t.loc[fit_idx], mp_old)
    L, T, keep, r = robust_affine(pI, pJ)
    t.loc[fit_idx[~keep], "role"] = "landmark_rejected"
    pI, pJ = pI[keep], pJ[keep]
    log(f"  windows {len(t)}, confident {len(conf)}: {keep.sum()} landmarks (+{(~keep).sum()} rejected), "
        f"{(t.role == 'heldout').sum()} held out; affine residual median {np.median(r[keep]):.1f} um, "
        f"max {r[keep].max():.1f} um\n  L = {np.round(L, 4).tolist()}, T = {np.round(T, 1).tolist()}")
    # held-out check of the affine alone (no LDDMM yet): how far do held-out pairs miss?
    ho = t[t.role == "heldout"]
    hI, hJ = landmark_pairs(ho, mp_old)
    log(f"  affine-only error on held-out pairs: median {np.median(np.linalg.norm(hI @ L.T + T - hJ, axis=1)):.1f} um")

    # affine-only transform (zero velocity field) = the starting point, and the result if LDDMM is skipped
    A0 = np.eye(3); A0[:2, :2], A0[:2, 2] = L, T
    tr = dict(A=torch.tensor(A0), v=torch.zeros(1, 2, 2, 2, dtype=torch.float64),
              xv=[torch.tensor([-1e6, 1e6], dtype=torch.float64)] * 2)
    init = dict(L=L, T=T)

    # 2. images (um frames), only needed for LDDMM
    if args.niter1 or args.niter2:
        hem, oI, pxI = he_image(s, mp_old, margin_um=400)
        dap, oJ, pxJ = dapi_image(s, mp_old.um_xy, margin_um=100)
        hem, dap = scale01(hem), scale01(np.log1p(dap))           # log: DAPI is heavy-tailed
        log(f"  images: H&E {hem.shape} at {pxI:.3f} um, DAPI {dap.shape} at {pxJ:.3f} um")

    # 4. LDDMM coarse -> fine; J binning chosen to match the H&E pixel size
    for stage, (b, niter, a) in enumerate([(args.b1, args.niter1, args.a1), (args.b2, args.niter2, args.a2)], 1):
        if niter == 0:
            continue
        I, xI = binned(hem, oI, pxI, b)
        J, xJ = to_j_frame(*binned(dap, oJ, pxJ, max(1, round(b * pxI / pxJ))))
        log(f"  stage {stage}: I {I.shape} at {b * pxI:.2f} um, J {J.shape} at {xJ[0][1] - xJ[0][0]:.2f} um, a {a}")
        tr, _ = lddmm(xI, I, xJ, J, pI, pJ, init, log, a=a, niter=niter, epV=args.epV, epL=args.epL,
                      epT=args.epT, sigmaM=args.sigmaM, sigmaB=args.sigmaM, sigmaA=2 * args.sigmaM,
                      sigmaP=args.sigmaP, sigmaR=args.sigmaR, nt=3, diffeo_start=0)
        lm = np.linalg.norm(to_he(tr, pJ) - pI, axis=1)
        hm = np.linalg.norm(to_he(tr, hJ) - hI, axis=1)
        log(f"    landmark error median {np.median(lm):.1f} um; held-out pair error median {np.median(hm):.1f} um, "
            f"p90 {np.quantile(hm, .9):.1f} um")
        init = dict(A=tr["A"].numpy(), v=tr["v"].numpy(), xv=[x.numpy() for x in tr["xv"]])
    np.savez(qc / f"transform_{s}.npz", A=tr["A"].numpy(), v=tr["v"].numpy(),
             xv0=tr["xv"][0].numpy(), xv1=tr["xv"][1].numpy(), args=json.dumps(vars(args)))

    # 5. map every cell and relabel
    he_um = to_he(tr, np.c_[old.x_centroid.values, -old.y_centroid.values])
    row_jpg, col_jpg = he_um[:, 0] / MPP0 / JPG_TO_L0, he_um[:, 1] / MPP0 / JPG_TO_L0
    new = write_meta(s, old, row_jpg, col_jpg, out_meta)
    moved = np.hypot(row_jpg - old.STAlign_x, col_jpg - old.STAlign_y) * JPG_TO_L0 * MPP0
    changed = (new.annotation != old.annotation)[old.annotation != "unannotated"].mean()
    log(f"  cells moved median {np.median(moved):.1f} um; annotated cells whose label changed {changed:.1%}")

    # 6. validation at the same window centres on the new map
    mp_new = Mapping(out_meta)
    ta = measure(s, mp_new, t[["row_l0", "col_l0"]].values, svs_for(s), DAPI(dapi_path(s)))
    ta["role"] = t.role.values
    ta.to_csv(qc / f"windows_after_{s}.tsv", sep="\t", index=False)
    t.to_csv(qc / f"windows_before_roles_{s}.tsv", sep="\t", index=False)
    rows = []
    for role in ["heldout", "landmark", "unconfident"]:
        for when, tt in (("before", t), ("after", ta)):
            rows.append(dict(sample=s, role=role, when=when, **summarise(tt[tt.role == role], role)))
    sm = pd.DataFrame(rows).drop(columns="set")
    sm.to_csv(qc / f"summary_{s}.tsv", sep="\t", index=False)
    log(sm.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
