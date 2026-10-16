#!/usr/bin/env python3
# What does each image pair's matching term reward? Similarity of the H&E and Xenium images as a function of
# a rigid offset from the nucleus-verified (new) map.
#
# Why    : the old STAlign (whole-slide JPG vs a cell-density raster) ends 3.5-79 um from the truth even when
#          started from exact landmarks, and where it ends changes with the start and the thread count
#          (replay_old_stalign.py). Its landmark term is at most 1.3% of the energy, so the image term and the
#          optimizer decide. This scan asks what the image term itself rewards, without the optimizer.
#          STAlign's image term is the squared error of a linear contrast fit, Xenium image ~ b0 + b . H&E
#          channels (EM; mixture weights aside), i.e. it rewards a high R^2 of that fit. This scan places every
#          Xenium pixel on the H&E through the new map plus an offset d (um, H&E row/col) and computes that R^2.
#          A pair whose R^2 peaks sharply at d = 0 pins the map there; a broad or off-centre peak lets LDDMM
#          settle anywhere inside it.
# Pairs  : (H&E side vs Xenium side; Xenium pixels inside the notebooks' raster area = the cells' bounding box
#          x 1.1, as STalign.rasterize draws it)
#   old    : JPG RGB (1.76 um/px)      vs cell-density raster: STalign.rasterize(dx = 10) of the ProSeg centroids,
#                                        i.e. the notebooks' J (Gaussian kernel sd = 2 dx = 20 um)
#   dapi10s: JPG RGB                   vs DAPI (log, robust 0-1, as the rerun prepares it), block mean to 10.2 um,
#                                        Gaussian sd 20 um like the raster -> swaps only the Xenium image's content
#   dapi10 : JPG RGB                   vs the same DAPI without the 20 um smoothing
#   dapi2  : JPG RGB                   vs DAPI at 1.7 um            -> also the pixel size
#   new    : SVS hematoxylin (2.02 um) vs DAPI at 1.7 um            -> the rerun's stage-2 images
#          Each Xenium image uses a fixed random subset of 150,000 of its pixels (seed 0). Each R^2 is computed
#          over all of them ("all") and over tissue pixels only ("tissue": a ProSeg centroid within 20 um), the
#          latter because STAlign's mixture weights can discount background.
# Offsets: +-150 um in 10 um steps, plus +-6 um in 1 um steps around 0.
# Output : NEW_DIR/qc/image_term/scan_<sample>.npz (R^2 per offset, pair and area),
#          --figure: scan_summary.tsv (one row per sample x pair x area: R^2 at the truth, the best offset and the
#          share of R^2 lost 2 / 5 / 10 / 30 / 100 um from the truth, plus the old map's median offset),
#          fig_image_term.png (profiles + best offsets) and fig_image_term_maps.png (R^2 maps for one sample).
# Usage  : python diag_image_term.py <sample> [...]      (conda env spatial)
#          python diag_image_term.py --figure
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, torch
from PIL import Image
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stalign_nucleus import NEW_DIR, he_image, dapi_image, binned, scale01, to_he
from nucleus_windows import Mapping, meta_path, folder, GROUP, DATA, MPP0, JPG_TO_L0, Z_CONF
from STalign import STalign

Image.MAX_IMAGE_PIXELS = None
OUT = NEW_DIR / "qc/image_term"
STEP = np.arange(-150, 150.1, 10.0)
FINE = np.arange(-6, 6.1, 1.0)
OFFSETS = np.unique(np.vstack([np.array(np.meshgrid(STEP, STEP, indexing="ij")).reshape(2, -1).T,
                               np.array(np.meshgrid(FINE, FINE, indexing="ij")).reshape(2, -1).T]), axis=0)
N_SUB, TISSUE_UM, RASTER_SD_UM = 150_000, 20.0, 20.0
PAIRS = {"old": ("jpg", "dens10"), "dapi10s": ("jpg", "dapi10s"), "dapi10": ("jpg", "dapi10"), "dapi2": ("jpg", "dapi2"),
         "new": ("hem", "dapi2")}


def cells_of(s):
    f = folder(s)
    return pd.read_csv(DATA / f"segmentation/proseg_run_121024/{f}/{f}_proseg/outs/cells.csv.gz",
                       usecols=["cell_id", "x_centroid", "y_centroid"])


def truth_map(s):
    """Xenium (x, y) um -> H&E (row, col) level-0 um through the new map's stored transform."""
    z = np.load(NEW_DIR / f"qc/transform_{s}.npz")
    tr = dict(A=torch.tensor(z["A"]), v=torch.tensor(z["v"]), xv=[torch.tensor(z["xv0"]), torch.tensor(z["xv1"])])
    return lambda xy: to_he(tr, np.c_[xy[:, 0], -xy[:, 1]])


def xenium_images(s, cells):
    """Xenium-side images inside the raster area: key -> (pixel values, pixel centres (x, y) um, tissue flag)."""
    x, y = cells.x_centroid.values, cells.y_centroid.values
    # the notebooks' J: rasterize(xM = -y, yM = -x, dx = 10); rows run along Y_ (= -x), cols along X_ (= -y)
    X_, Y_, M = STalign.rasterize(-y, -x, dx=10, draw=0)
    gy, gx = np.meshgrid(Y_, X_, indexing="ij")
    out = {"dens10": (M[0].ravel(), np.c_[-gy.ravel(), -gx.ravel()])}
    xlo, xhi, ylo, yhi = -Y_.max(), -Y_.min(), -X_.max(), -X_.min()
    # DAPI as the rerun prepares it (stalign_nucleus.main): log, robust 0-1 over the crop, then block means
    dap, oJ, pxJ = dapi_image(s, np.c_[x, y], margin_um=100)
    dap = scale01(np.log1p(dap))
    for key, b, sd_um in (("dapi10s", 12, RASTER_SD_UM), ("dapi10", 12, 0), ("dapi2", 2, 0)):
        im, (cy, cx) = binned(dap, oJ, pxJ, b)
        if sd_um:
            im = ndi.gaussian_filter(im, sd_um / (cy[1] - cy[0]))
        Y, X = np.meshgrid(cy, cx, indexing="ij")
        keep = ((X >= xlo) & (X <= xhi) & (Y >= ylo) & (Y <= yhi)).ravel()
        out[key] = (im.ravel()[keep], np.c_[X.ravel()[keep], Y.ravel()[keep]])
    for key, (v, xy) in out.items():
        k = np.random.default_rng(0).choice(len(v), min(N_SUB, len(v)), replace=False)
        k = k[np.lexsort((xy[k, 1], np.floor(xy[k, 0] / 8)))]       # 8 um bands of x = H&E rows: cache-friendly sampling
        out[key] = (v[k], xy[k])
    tree = cKDTree(np.c_[x, y])
    return {key: (v, xy, tree.query(xy, distance_upper_bound=TISSUE_UM)[0] < np.inf) for key, (v, xy) in out.items()}


def he_images(s, mp_old):
    """H&E-side images: key -> (channels x H x W float32, function H&E um (row, col) -> pixel index)."""
    # slide JPG (what the notebooks registered), cropped to the old map's extent of the cells + 400 um
    lo = np.floor((mp_old.l0_rc.min(0) - 400 / MPP0) / JPG_TO_L0).astype(int).clip(0)
    hi = np.ceil((mp_old.l0_rc.max(0) + 400 / MPP0) / JPG_TO_L0).astype(int)
    jpg = Image.open(DATA / f"Xenium_annotations/images/230918_Xenium_CytAssist_{GROUP[s]}.jpg")
    rgb = np.asarray(jpg.crop((lo[1], lo[0], hi[1], hi[0])), np.float32).transpose(2, 0, 1).copy()
    # JPG px i covers level-0 px [7i, 7i + 7), so its centre sits at (7i + 3.5) * MPP0 um
    jpg_idx = lambda um: um / MPP0 / JPG_TO_L0 - 0.5 - lo
    # SVS hematoxylin as the rerun's stage 2 sees it (page 2, binned 2)
    hem, oI, pxI = he_image(s, mp_old, margin_um=400)
    hem, (cr, cc) = binned(scale01(hem), oI, pxI, 2)
    hem_idx = lambda um: (um - [cr[0], cc[0]]) / (cr[1] - cr[0])
    return {"jpg": (rgb, jpg_idx), "hem": (hem[None].astype(np.float32), hem_idx)}


def r2_scan(I, idx, vals, he_um, masks):
    """R^2 of vals ~ b0 + b . I(he_um + d) for every offset d in OFFSETS, one column per mask. I is sampled
    bilinearly with border padding, as STalign.interp does (torch grid_sample, align_corners=True)."""
    W = np.stack(masks, 1).astype(np.float64)                         # one 0/1 weight column per mask
    n, Sy = W.sum(0), W.T @ vals
    Syy = W.T @ vals ** 2 - Sy ** 2 / n                                # total sum of squares per mask
    It = torch.from_numpy(I)[None]
    size = np.array(I.shape[1:]) - 1.0
    out = np.empty((len(OFFSETS), len(masks)))
    for k, d in enumerate(OFFSETS):
        g = idx(he_um + d) / size * 2 - 1                                   # (row, col) -> [-1, 1]
        g = torch.from_numpy(g[:, ::-1].astype(np.float32).copy())[None, None]      # grid_sample wants (x = col, y = row)
        X = torch.nn.functional.grid_sample(It, g, align_corners=True, padding_mode="border")[0, :, 0].T.double().numpy()
        Sx, Sxy = W.T @ X, (W * vals[:, None]).T @ X                    # per mask: sums of X and of X * y
        for j in range(W.shape[1]):
            XtX = (X * W[:, j:j + 1]).T @ X - np.outer(Sx[j], Sx[j]) / n[j]      # centred cross-products
            Xty = Sxy[j] - Sx[j] * Sy[j] / n[j]
            beta = np.linalg.solve(XtX + 1e-9 * np.eye(len(Xty)), Xty)
            out[k, j] = beta @ Xty / Syy[j]                             # explained / total sum of squares
    return out


RADII = np.r_[1:7, 10:151:10]                       # profile radii: full rings of the 1 um grid (<= 6 um), then the 10 um grid


def ring_mean(v, dist, r):
    """Mean R^2 over the offsets at distance ~r from the truth (+-0.5 um on the 1 um grid, +-3.6 um on the 10 um grid)."""
    return v[np.abs(dist - r) <= (0.5 if r < 10 else 3.6)].mean()


def summarise(s, z):
    """One row per pair x area: R^2 at the truth, the best offset, and the share of R^2 lost 2 / 5 / 10 / 30 / 100 um away."""
    rows, off = [], z["offsets"]
    dist = np.hypot(*off.T)
    # the stored map's typical error: median offset length of its confident windows (old_map = median vector,
    # which understates it where the error is non-rigid and changes direction, e.g. CL87352)
    w = pd.read_csv(NEW_DIR / f"qc/windows_before_{s}.tsv", sep="\t").query("peak_z >= @Z_CONF")
    for p in PAIRS:
        for j, area in enumerate(("all", "tissue")):
            v = z[p][:, j]; v0 = v[dist == 0][0]; b = v.argmax()
            rows.append(dict(sample=s, pair=p, area=area, r2_truth=v0, best_row_um=off[b, 0], best_col_um=off[b, 1],
                             best_um=dist[b], r2_best=v[b],
                             **{f"loss_{r}um": (v0 - ring_mean(v, dist, r)) / v0 for r in (2, 5, 10, 30, 100)},
                             old_map_row_um=z["old_map"][0], old_map_col_um=z["old_map"][1], old_map_median_um=w.shift_um.median()))
    return rows


def run(s):
    t0 = time.time()
    cells = cells_of(s)
    mp_old = Mapping(meta_path(s))
    xen = xenium_images(s, cells)
    hes = he_images(s, mp_old)
    to_um = truth_map(s)
    r2 = {}
    for p, (hk, xk) in PAIRS.items():
        vals, xy, tissue = xen[xk]
        I, idx = hes[hk]
        r2[p] = r2_scan(I, idx, vals.astype(np.float64), to_um(xy), [np.ones(len(vals), bool), tissue])
        print(f"  {s} {p}: {len(vals)} px, R^2 at truth {r2[p][OFFSETS.tolist().index([0.0, 0.0])].round(3)}, "
              f"{time.time() - t0:.0f} s", flush=True)
    # the old 030226 map's error: median offset of its confident windows (same sign convention as d)
    w = pd.read_csv(NEW_DIR / f"qc/windows_before_{s}.tsv", sep="\t").query("peak_z >= @Z_CONF")
    old_vec = w[["shift_row_um", "shift_col_um"]].median().values
    np.savez_compressed(OUT / f"scan_{s}.npz", offsets=OFFSETS, old_map=old_vec, **{p: a for p, a in r2.items()})
    print(pd.DataFrame(summarise(s, np.load(OUT / f"scan_{s}.npz"))).round(3).to_string(index=False))


# ---------------------------------------------------------------- figures
NAMES = {"old": "H&E colour vs cell-count map\n(the old pipeline)", "dapi10s": "H&E colour vs DAPI,\nblurred like the map",
         "dapi10": "H&E colour vs DAPI,\n10 µm pixels, unblurred", "dapi2": "H&E colour vs DAPI,\n1.7 µm pixels",
         "new": "hematoxylin vs DAPI, 1.7 µm\n(the new pipeline)"}
# dataviz reference palette, slots 1-5 (validated: adjacent CVD dE >= 9.1; the low-contrast slots get direct labels)
COLOR = {"new": "#2a78d6", "old": "#eb6834", "dapi2": "#1baf7a", "dapi10s": "#eda100", "dapi10": "#e87ba4"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": MUTED, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    scans = {f.stem[5:]: np.load(f) for f in sorted(OUT.glob("scan_*.npz"))}
    summ = pd.concat([pd.DataFrame(summarise(s, z)) for s, z in scans.items()])
    summ.to_csv(OUT / "scan_summary.tsv", sep="\t", index=False)

    # 1. tissue-area profile (R^2 left at distance r, relative to the truth) + where each sample's match peaks
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    for p in PAIRS:
        prof = []
        for z in scans.values():
            v = z[p][:, 1]; dist = np.hypot(*z["offsets"].T)
            prof.append([ring_mean(v, dist, r) / v[dist == 0][0] for r in RADII])
        prof = np.array(prof)
        a.plot(RADII, prof.T, color=COLOR[p], lw=0.6, alpha=0.25)
        med = np.median(prof, 0)
        a.plot(RADII, med, color=COLOR[p], lw=2, label=NAMES[p].replace("\n", " "))
    a.set_xscale("log"); a.set_xlim(1, 150); a.set_ylim(0, 1.03)
    a.set_xticks([1, 2, 5, 10, 20, 50, 100, 150], ["1", "2", "5", "10", "20", "50", "100", "150"])
    a.grid(color=GRID, lw=0.6); a.set_axisbelow(True)
    a.legend(frameon=False, fontsize=7.5, loc="lower left")
    a.set_xlabel("map moved away from the true position (µm, log)")
    a.set_ylabel("image match left (R², tissue pixels; 1 = at the truth)")
    a.set_title("How sharply each image pair pins the map\n(thick = median of samples, thin = each sample)", loc="left", fontsize=9.5)
    b = ax[1]
    t = summ[(summ.area == "tissue") & summ.pair.isin(["old", "new"])].pivot(index="sample", columns="pair", values="best_um")
    om = summ.drop_duplicates("sample").set_index("sample").old_map_median_um.loc[t.index]
    y = np.arange(len(t))[::-1]
    b.hlines(y, t.old, om, color=GRID, lw=2, zorder=1)
    b.scatter(om, y, marker="D", s=30, color=INK, zorder=3, label="stored 030226 map (median window offset)")
    b.scatter(t.old, y, s=40, color=COLOR["old"], edgecolor="white", lw=1.2, zorder=4, label="best match, old image pair")
    b.scatter(t.new, y, s=22, facecolor="none", edgecolor=COLOR["new"], lw=1.2, zorder=5, label="best match, new image pair")
    b.set_yticks(y, t.index); b.set_xlim(-2, 100); b.set_ylim(-0.7, len(t) + 1.3)
    b.set_xlabel("distance from the true position (µm)")
    b.set_title("Best image match vs the stored map", loc="left", fontsize=9.5)
    b.legend(frameon=False, fontsize=7.5, loc="upper right", ncol=1)
    fig.tight_layout(); fig.savefig(OUT / "fig_image_term.png", dpi=150); plt.close(fig)

    # 2. R^2 maps over the 10 um offset grid for the worst old-map sample (tissue pixels, relative to the truth)
    s = "CL42112" if "CL42112" in scans else next(iter(scans))
    z = scans[s]; off = z["offsets"]; coarse = (np.abs(off) % 10 == 0).all(1)
    g = np.unique(off[coarse, 0]); n = len(g)
    fig, ax = plt.subplots(1, len(PAIRS), figsize=(3.0 * len(PAIRS), 3.5))
    for a, p in zip(ax, PAIRS):
        v = z[p][:, 1]; v0 = v[(off == 0).all(1)][0]
        m = np.full((n, n), np.nan)
        m[np.searchsorted(g, off[coarse, 0]), np.searchsorted(g, off[coarse, 1])] = v[coarse] / v0
        im = a.imshow(m, extent=(g[0] - 5, g[-1] + 5, g[-1] + 5, g[0] - 5), cmap="Blues", vmin=0, vmax=1)
        a.scatter(0, 0, marker="+", s=60, color=COLOR["old"], lw=1.5)
        a.scatter(z["old_map"][1], z["old_map"][0], marker="D", s=26, facecolor="none", edgecolor=INK, lw=1.2)
        a.set_title(NAMES[p], fontsize=8.5, loc="left")
        a.set_xlabel("offset, H&E col (µm)")
    ax[0].set_ylabel("offset, H&E row (µm)")
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.01); cb.set_label("image match left (1 = at the truth)", color=MUTED)
    fig.suptitle(f"{s}: image match vs map offset, tissue pixels (+ = true position, diamond = stored 030226 map)",
                 x=0.01, ha="left", fontsize=9.5)
    fig.savefig(OUT / "fig_image_term_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    print(summ[summ.area == "tissue"].round(3).to_string(index=False))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1:] == ["--figure"]:
        return figure()
    torch.set_num_threads(1)
    for s in sys.argv[1:]:
        run(s)


if __name__ == "__main__":
    main()
