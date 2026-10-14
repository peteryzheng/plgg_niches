#!/usr/bin/env python3
# Why H&E finds more "nuclei" than DAPI: a one-to-one accounting of the peaks the window check calls.
#
# For every held-out window of every sample, under the NEW map at zero shift, peaks are called exactly as in
# nucleus_windows.measure() (same reads, padding, smoothing, Otsu threshold and spacing; checked against
# peaks() on every window). Each peak inside the window (8 um border excluded, so partners just outside
# still count) is put in exactly one class:
#   paired      : H&E and DAPI peaks that are each other's nearest neighbour, <= 2.5 um apart (one-to-one,
#                 so #paired H&E == #paired DAPI)
#   duplicate   : the nearest other-modality peak is <= 2.5 um away but is paired with another peak
#                 (two peaks of one modality on one nucleus of the other)
#   near_miss   : nearest other-modality peak 2.5-6 um away
#   no_partner  : nothing from the other modality within 6 um
# Per peak it also records strength = smoothed intensity at the peak / the image's Otsu threshold (> 1 by
# construction), and the OTHER modality's signal at the same spot: other_signal = its smoothed intensity / its
# own Otsu threshold, and other_rel = other_signal / the median other_signal at paired peaks of the same window
# (1 = as strong as under a typical matched nucleus, 0 = nothing). other_rel is the readable one: in tissue the
# Otsu threshold mostly separates tissue from glass, not nucleus from background.
# --figure-only redraws the figure and summary from detection_peaks.tsv.gz without re-reading images.
# Outputs (NEW_DIR/qc/detection/): detection_peaks.tsv.gz (one row per peak), detection_summary.tsv,
# fig_detection.png, detect_example_<sample>.jpg (patient H&E: SendUserFile / Notion only).
# Usage: python diag_detection.py [sample ...]     (default: all 11)
#        python diag_detection.py --figure-only
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import color, feature, filters
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nucleus_core import peaks
from nucleus_windows import (GROUP, DATA, Mapping, meta_path, svs_for, dapi_path, DAPI, MPP0, DAPI_UM, WIN_UM,
                             SHIFT_UM, NUC_MIN_DIST_UM, F, UM)

NEW_DIR = DATA / "Xenium_annotations/pathology_annotations_101326_nucleus"
OUT = NEW_DIR / "qc/detection"
MATCH, NEAR, BORDER = 2.5, 6.0, 8.0          # um
CLASSES = ["paired", "duplicate", "near_miss", "no_partner"]


def called(img, sigma_px, min_dist_px):
    """peaks() plus what it computes internally: smoothed image and Otsu threshold."""
    sm = ndi.gaussian_filter(img, sigma_px)
    thr = filters.threshold_otsu(sm)
    pk = feature.peak_local_max(sm, min_distance=int(min_dist_px), threshold_abs=thr, exclude_border=False)
    return pk, sm, thr


def sample_at(sm, thr, rc):
    """Smoothed intensity / threshold at fractional pixel positions (NaN outside the image)."""
    v = ndi.map_coordinates(sm, rc.T, order=1, cval=np.nan)
    return v / thr


def classify(a, b):
    """Class of each point of a (vs b) and of b (vs a); one-to-one pairs = mutual nearest within MATCH."""
    da, ia = cKDTree(b).query(a)
    db, ib = cKDTree(a).query(b)
    pa = (da <= MATCH) & (ib[ia] == np.arange(len(a)))
    pb = (db <= MATCH) & (ia[ib] == np.arange(len(b)))
    lab = lambda d, p: np.where(p, "paired", np.where(d <= MATCH, "duplicate", np.where(d <= NEAR, "near_miss", "no_partner")))
    return lab(da, pa), da, lab(db, pb), db, ib


def run_sample(s):
    mp = Mapping(meta_path(s, NEW_DIR))
    w = pd.read_csv(NEW_DIR / f"qc/windows_after_{s}.tsv", sep="\t")
    w = w[(w.role == "heldout") & w.n_he.ge(20)]
    svs, dapi = svs_for(s), DAPI(dapi_path(s))
    half0 = int(WIN_UM / 2 / MPP0)
    rows = []
    for _, r in w.iterrows():
        cy, cx = int(r.row_l0), int(r.col_l0)
        # H&E: exactly measure(): 2x2-averaged RGB (0.504 um/px) -> rgb2hed hematoxylin -> peaks
        rgb = svs.read_scaled(cy - half0, cx - half0, 2 * half0, 2 * half0, F).astype(np.uint8)
        hem = color.rgb2hed(rgb)[..., 0]
        hp, hsm, hthr = called(hem, 1.0 / UM, NUC_MIN_DIST_UM / UM)
        assert np.array_equal(hp, peaks(hem, 1.0 / UM, NUC_MIN_DIST_UM / UM))
        he = hp * UM
        # DAPI: exactly measure(): the mapped window padded by SHIFT_UM, raw morphology_focus at 0.2125 um/px
        corners = np.array([[cy - half0, cx - half0], [cy + half0, cx + half0], [cy - half0, cx + half0], [cy + half0, cx - half0]])
        cu = mp.to_um(corners)
        if np.isnan(cu).any():
            continue
        x0, x1 = cu[:, 0].min() - SHIFT_UM, cu[:, 0].max() + SHIFT_UM
        y0, y1 = cu[:, 1].min() - SHIFT_UM, cu[:, 1].max() + SHIFT_UM
        r0, c0 = int(y0 / DAPI_UM), int(x0 / DAPI_UM)
        img = dapi.read(r0, c0, int((y1 - y0) / DAPI_UM), int((x1 - x0) / DAPI_UM))
        dp, dsm, dthr = called(img, 1.0 / DAPI_UM, NUC_MIN_DIST_UM / DAPI_UM)
        d_l0 = mp.to_l0(np.c_[(dp[:, 1] + c0) * DAPI_UM, (dp[:, 0] + r0) * DAPI_UM])
        ok = ~np.isnan(d_l0).any(1)
        dw, dstr = (d_l0[ok] - [cy - half0, cx - half0]) * MPP0, dsm[dp[ok, 0], dp[ok, 1]] / dthr
        near = ((dw >= -NEAR - 1) & (dw <= WIN_UM + NEAR + 1)).all(1)      # only DAPI that can partner window peaks
        dw, dstr = dw[near], dstr[near]
        hlab, hdist, dlab, ddist, partner = classify(he, dw)
        inner = lambda p: ((p >= BORDER) & (p <= WIN_UM - BORDER)).all(1)
        # the other modality's signal at each peak: H&E peak -> Xenium um -> DAPI px; DAPI peak -> H&E px
        he_um = mp.to_um(he / MPP0 + [cy - half0, cx - half0])
        dapi_at_he = sample_at(dsm, dthr, np.c_[he_um[:, 1] / DAPI_UM - r0, he_um[:, 0] / DAPI_UM - c0])
        hem_at_dapi = sample_at(hsm, hthr, dw / UM)
        nn = lambda p: cKDTree(p).query(p, k=2)[0][:, 1]
        # a pair is counted where its H&E member sits, so both modalities count the same pairs
        m_he = inner(he)
        m_dapi = np.where(dlab == "paired", m_he[partner], inner(dw))
        for mod, pts, m, lab, dist, st, other, nns in (
                ("H&E", he, m_he, hlab, hdist, hsm[hp[:, 0], hp[:, 1]] / hthr, dapi_at_he, nn(he)),
                ("DAPI", dw, m_dapi, dlab, ddist, dstr, hem_at_dapi, nn(dw))):
            rows.append(pd.DataFrame(dict(sample=s, window=int(r.window), modality=mod, cls=lab[m],
                                          nearest_other_um=dist[m], strength=st[m], other_signal=other[m],
                                          nearest_same_um=nns[m], row_um=pts[m, 0], col_um=pts[m, 1])))
    t = pd.concat(rows, ignore_index=True)
    print(f"{s}: {t.window.nunique()} windows, {len(t)} peaks", flush=True)
    return t


def example_window(s):
    """The held-out window whose new-map offset is closest to the sample's median (a typical window)."""
    w = pd.read_csv(NEW_DIR / f"qc/windows_after_{s}.tsv", sep="\t")
    w = w[(w.role == "heldout") & w.n_he.ge(20)].reset_index(drop=True)
    return int(w.window[(w.shift_um - w.shift_um.median()).abs().argmin()])


def draw_example(s, k, t, out):
    """One window, two panels with the same peak classes: the H&E, and the Xenium DAPI resampled into the H&E
    window through the new map (so both show the same tissue, pixel for pixel)."""
    mp = Mapping(meta_path(s, NEW_DIR))
    r = pd.read_csv(NEW_DIR / f"qc/windows_after_{s}.tsv", sep="\t").set_index("window").loc[k]
    cy, cx, half0 = int(r.row_l0), int(r.col_l0), int(WIN_UM / 2 / MPP0)
    rgb = svs_for(s).read_scaled(cy - half0, cx - half0, 2 * half0, 2 * half0, F).astype(np.uint8)
    # every H&E pixel centre -> Xenium um -> DAPI px; read that DAPI area once and sample it
    rr, cc = np.mgrid[0:rgb.shape[0], 0:rgb.shape[1]]
    um = mp.to_um(np.c_[rr.ravel(), cc.ravel()] * F + 0.5 * F + [cy - half0, cx - half0])
    y0, x0 = np.nanmin(um[:, 1]) / DAPI_UM - 5, np.nanmin(um[:, 0]) / DAPI_UM - 5
    y1, x1 = np.nanmax(um[:, 1]) / DAPI_UM + 5, np.nanmax(um[:, 0]) / DAPI_UM + 5
    img = DAPI(dapi_path(s)).read(int(y0), int(x0), int(y1 - y0), int(x1 - x0))
    warped = ndi.map_coordinates(img, [um[:, 1] / DAPI_UM - int(y0), um[:, 0] / DAPI_UM - int(x0)], order=1,
                                 cval=np.nan).reshape(rgb.shape[:2])
    lo, hi = np.nanpercentile(warped, [1, 99.5])
    x = t[(t["sample"] == s) & (t.window == k)]
    sty = {"paired": ("#1f5fae", "paired"), "duplicate": ("#7fb2e5", "duplicate"),
           "near_miss": ("#e5a50a", "near miss 2.5-6 um"), "no_partner": ("#d1242f", "no partner within 6 um")}
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2))
    for ax, im, ttl in ((axes[0], rgb, "H&E"), (axes[1], np.clip((warped - lo) / (hi - lo), 0, 1), "Xenium DAPI (morphology_focus), same window")):
        ax.imshow(im, cmap="gray" if im.ndim == 2 else None)
        for mod, mk in (("H&E", "o"), ("DAPI", "+")):
            for c, (col, name) in sty.items():
                p = x[(x.modality == mod) & (x.cls == c)]
                ax.scatter(p.col_um / UM, p.row_um / UM, s=34, marker=mk, lw=1.3, label=f"{mod}: {name} ({len(p)})",
                           **({"facecolors": "none", "edgecolors": col} if mk == "o" else {"color": col}))
        ax.set_xlim(0, rgb.shape[1]); ax.set_ylim(rgb.shape[0], 0); ax.axis("off"); ax.set_title(ttl, fontsize=10)
        ax.plot([8, 8 + 25 / UM], [rgb.shape[0] - 8] * 2, color="k" if im.ndim == 3 else "w", lw=3)
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)
    fig.suptitle(f"{s}, held-out window {k} (150 um, new map): circles = H&E peaks, + = DAPI peaks, colour = class; "
                 "8 um border not scored; bar 25 um", fontsize=10)
    fig.savefig(out, dpi=100, bbox_inches="tight", pil_kwargs={"quality": 85}); plt.close(fig)


def figure(t, out):
    col = {"paired": "#1f5fae", "duplicate": "#7fb2e5", "near_miss": "#e5a50a", "no_partner": "#d1242f"}
    samples = [s for s in GROUP if s in set(t["sample"])]
    frac = t.groupby(["sample", "modality"]).cls.value_counts(normalize=True).unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0)
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1], hspace=0.45, wspace=0.35)
    # A: every peak in exactly one class -> each bar sums to 100%
    a = fig.add_subplot(gs[0, :])
    ypos = []
    for i, s in enumerate(samples):
        for j, mod in enumerate(("H&E", "DAPI")):
            y = i * 2.6 + j; ypos.append((y, f"{s}  {mod}")); left = 0
            n = int((t["sample"].eq(s) & t.modality.eq(mod)).sum())
            for c in CLASSES:
                v = frac.loc[(s, mod), c]
                a.barh(y, v, left=left, color=col[c], height=0.85, edgecolor="white", lw=1)
                if v > 0.06:
                    a.text(left + v / 2, y, f"{v:.0%}", ha="center", va="center", fontsize=7, color="white")
                left += v
            a.text(1.01, y, f"n = {n:,}", va="center", fontsize=7, color="#52514e")
    a.set_yticks([p[0] for p in ypos], [p[1] for p in ypos], fontsize=7); a.invert_yaxis()
    a.set_xlim(0, 1); a.set_xlabel("share of peaks (each bar = 100% of that modality's peaks)")
    a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col[c]) for c in CLASSES],
             labels=["paired one-to-one (<= 2.5 um, mutual nearest)", "duplicate (partner taken by another peak)",
                     "near miss (nearest other peak 2.5-6 um)", "no partner within 6 um"],
             ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=8)
    a.set_title("A. Every called peak, by what lies under it in the other modality (held-out windows, new map)", loc="left",
                fontsize=10, pad=26)
    # B / C: per class, pooled distribution (violin) + one dot per sample (median)
    for k, (mod, var, ttl, ylab) in enumerate([
            ("H&E", "strength", "B. H&E peak strength", "hematoxylin at peak / window Otsu threshold"),
            ("DAPI", "strength", "B. DAPI peak strength", "DAPI at peak / Otsu threshold"),
            ("H&E", "other_rel", "C. DAPI under each H&E peak", "DAPI there / DAPI at matched nuclei"),
            ("DAPI", "other_rel", "C. Hematoxylin under each DAPI peak", "hematoxylin there / at matched nuclei")]):
        ax = fig.add_subplot(gs[1, k]); x = t[t.modality == mod]
        data = [x.loc[x.cls == c, var].dropna().values for c in CLASSES]
        keep = [i for i, d in enumerate(data) if len(d) > 5]
        vp = ax.violinplot([data[i] for i in keep], positions=keep, showextrema=False, widths=0.8)
        for b, i in zip(vp["bodies"], keep):
            b.set_facecolor(col[CLASSES[i]]); b.set_alpha(0.35)
        med = x.groupby(["sample", "cls"])[var].median().unstack().reindex(columns=CLASSES)
        for i in keep:
            ax.scatter(np.full(len(med), i) + np.linspace(-0.18, 0.18, len(med)), med[CLASSES[i]], s=12,
                       color=col[CLASSES[i]], edgecolor="white", lw=0.5, zorder=3)
        if var == "other_rel":
            ax.axhline(1, color="#52514e", lw=1, ls=(0, (3, 3)))
            ax.set_ylim(-0.1, 2.5)
        ax.set_xticks(range(4), ["paired", "dup.", "near\nmiss", "no\npartner"], fontsize=8)
        ax.set_title(ttl, loc="left", fontsize=9); ax.set_ylabel(ylab, fontsize=8)
    fig.text(0.5, 0.005, "B, C: violins = all peaks pooled; dots = one sample's median. B: 1 = the detection threshold. "
             "C: 1 (dashed) = the other stain's level under a typical matched nucleus in the same window; 0 = none.",
             ha="center", fontsize=8, color="#52514e")
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1:] == ["--figure-only"]:
        t = pd.read_csv(OUT / "detection_peaks.tsv.gz", sep="\t").drop(columns="other_rel", errors="ignore")
    else:
        parts = []
        for s in sys.argv[1:] or list(GROUP):
            parts.append(run_sample(s))
        t = pd.concat(parts, ignore_index=True)
    # other-modality signal relative to its level under the paired (matched) peaks of the same window
    ref = t[t.cls == "paired"].groupby(["sample", "window", "modality"]).other_signal.median().rename("ref")
    t = t.join(ref, on=["sample", "window", "modality"])
    t["other_rel"] = t.other_signal / t.pop("ref")
    t.to_csv(OUT / "detection_peaks.tsv.gz", sep="\t", index=False)
    for s in ("CL42112", "PA188570"):                    # a dense CL and a looser PA example
        if s in set(t["sample"]):
            draw_example(s, example_window(s), t, OUT / f"detect_example_{s}.jpg")
    # per sample x modality: counts per class (they add up), and medians per class
    g = t.groupby(["sample", "modality"])
    summ = g.cls.value_counts().unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0).add_prefix("n_")
    summ["n_total"] = summ.sum(1)
    summ["validation_match_share"] = (summ.n_paired + summ.n_duplicate) / summ.n_total   # the 'matched' metric
    for v in ("strength", "other_signal", "other_rel", "nearest_same_um"):
        summ = summ.join(t.groupby(["sample", "modality", "cls"])[v].median().unstack().reindex(columns=CLASSES).add_prefix(f"{v}_med_"))
    summ.to_csv(OUT / "detection_summary.tsv", sep="\t")
    figure(t, OUT / "fig_detection.png")
    print(summ[[c for c in summ.columns if c.startswith("n_")]].to_string())


if __name__ == "__main__":
    main()
