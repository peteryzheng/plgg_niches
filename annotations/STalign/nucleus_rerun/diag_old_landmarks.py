#!/usr/bin/env python3
# How good were the original notebooks' hand-placed STAlign landmarks, and is the 50-unit grid they sit on
# the problem?
#
# Each notebook (annotations/STalign/pathology_<sample>_proseg.qmd) pairs pointsI (H&E JPG px, row/col)
# with pointsJ (Xenium um, row = -x_centroid, col = -y_centroid). For every landmark:
#   clicked  : pointsI, where the landmark was put on the H&E
#   truth    : where the new, nucleus-verified map puts pointsJ on the H&E
#   old_map  : where the final 030226 map puts pointsJ on the H&E
# A landmark's error is |clicked - truth|; it mixes misplacement on either image (the two can't be told
# apart). Rounding bound: if both points had been placed perfectly and then rounded to the 50 grid, the pair
# would be off by at most 25*sqrt(2) JPG px (62 um) + 25*sqrt(2) um (35 um) = 97 um; the typical value comes
# from simulating uniform rounding on both sides.
# Outputs (NEW_DIR/qc/landmarks/): old_landmarks.tsv, fig_old_landmarks_errors.png (no image data),
# fig_old_landmarks_map.jpg (H&E thumbnails: SendUserFile / Notion only).
# Usage: python diag_old_landmarks.py
import sys
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image
from scipy.interpolate import LinearNDInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from replay_old_stalign import notebook_landmarks
from nucleus_windows import GROUP, DATA, meta_path, MPP0, JPG_TO_L0

Image.MAX_IMAGE_PIXELS = None
NEW_DIR = DATA / "Xenium_annotations/pathology_annotations_101326_nucleus"
OUT = NEW_DIR / "qc/landmarks"
JPG_UM = JPG_TO_L0 * MPP0                  # um per JPG px (1.765)


def map_jpg(meta_csv, xy):
    """Xenium (x, y) um -> H&E (row, col) JPG px through a meta file: piecewise-linear over its cells, global
    least-squares affine outside their hull (landmarks often sit on the tissue edge)."""
    m = pd.read_csv(meta_csv, usecols=["x_centroid", "y_centroid", "STAlign_x", "STAlign_y"]).dropna()
    src, dst = m[["x_centroid", "y_centroid"]].values, m[["STAlign_x", "STAlign_y"]].values
    out = LinearNDInterpolator(src, dst)(xy)
    M, *_ = np.linalg.lstsq(np.c_[src, np.ones(len(src))], dst, rcond=None)
    bad = np.isnan(out).any(1)
    out[bad] = np.c_[xy[bad], np.ones(bad.sum())] @ M
    return out, m


def rounding_only(n=200_000, seed=0):
    """Pair error (um) from rounding a perfect H&E point to 50 JPG px and a perfect Xenium point to 50 um."""
    r = np.random.default_rng(seed)
    e = r.uniform(-25, 25, (n, 2)) * JPG_UM + r.uniform(-25, 25, (n, 2))
    return np.hypot(*e.T)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, cells = [], {}
    for s in GROUP:
        pI, pJ = notebook_landmarks(s)
        xy = np.c_[-pJ[:, 0], -pJ[:, 1]]                                   # back to (x_centroid, y_centroid)
        truth, m = map_jpg(meta_path(s, NEW_DIR), xy)
        oldm, _ = map_jpg(meta_path(s), xy)
        cells[s] = m[["STAlign_x", "STAlign_y"]].values
        on_grid = bool(((pI % 50) == 0).all() and ((pJ % 50) == 0).all())
        for i in range(len(pI)):
            rows.append(dict(sample=s, landmark=i, clicked_row=pI[i, 0], clicked_col=pI[i, 1],
                             pointsJ_row=pJ[i, 0], pointsJ_col=pJ[i, 1], notebook_on_50_grid=on_grid,
                             truth_row=truth[i, 0], truth_col=truth[i, 1], old_map_row=oldm[i, 0], old_map_col=oldm[i, 1],
                             landmark_err_um=np.hypot(*(pI[i] - truth[i])) * JPG_UM,
                             old_map_err_um=np.hypot(*(oldm[i] - truth[i])) * JPG_UM,
                             old_map_to_click_um=np.hypot(*(oldm[i] - pI[i])) * JPG_UM,
                             truth_to_nearest_grid_um=np.hypot(*(np.round(truth[i] / 50) * 50 - truth[i])) * JPG_UM))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "old_landmarks.tsv", sep="\t", index=False)
    ro = rounding_only()
    print(f"rounding-only pair error: median {np.median(ro):.0f} um, 95th pct {np.quantile(ro, .95):.0f}, max 97 um")
    print(t.groupby("sample")[["landmark_err_um", "old_map_err_um", "old_map_to_click_um"]].median().round(0).join(
        t.groupby("sample").notebook_on_50_grid.first()).to_string())

    # ---- errors (no image data)
    samples = list(GROUP)[::-1]
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.axvspan(0, 97, color="#e9e8e3", zorder=0)
    ax.axvline(np.median(ro), color="#8a8984", lw=1, ls=(0, (3, 3)), zorder=1)
    ax.text(48, len(samples) - 0.35, "range 50-grid rounding\nalone can produce", ha="center", va="bottom", fontsize=7.5, color="#52514e")
    for y, s in enumerate(samples):
        x = t[t["sample"] == s]
        ax.scatter(x.landmark_err_um, np.full(len(x), y + 0.12), marker="x", s=36, color="#d1242f", lw=1.4, zorder=3,
                   label="hand landmark vs truth" if y == 0 else None)
        ax.scatter(x.old_map_err_um, np.full(len(x), y - 0.12), s=30, facecolor="none", edgecolor="#1f5fae", lw=1.3, zorder=3,
                   label="old final map vs truth, same spot" if y == 0 else None)
        grid = "" if x.notebook_on_50_grid.iloc[0] else "   (not on grid)"
        ax.text(705, y, f"{x.landmark_err_um.median():.0f} -> {x.old_map_err_um.median():.0f}{grid}", va="center", fontsize=7.5, color="#52514e")
    ax.set_yticks(range(len(samples)), samples); ax.set_xlim(0, 820); ax.set_ylim(-0.6, len(samples) + 0.4)
    ax.set_xlabel("distance on the H&E from the true position (um)")
    ax.text(705, len(samples) - 0.35, "median um:\nlandmark -> old map", fontsize=7.5, color="#52514e", va="bottom")
    ax.legend(loc="lower right", bbox_to_anchor=(0.86, 0.0), frameon=False, fontsize=8)
    ax.set_title("Original STAlign landmarks: each x is one landmark; dashed = typical rounding error "
                 f"({np.median(ro):.0f} um)", loc="left", fontsize=9)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "fig_old_landmarks_errors.png", dpi=150); plt.close(fig)

    # ---- where they were: per sample, a zoom on its median-error landmark (H&E JPG, the frame clicked in)
    imgs = {g: np.asarray(Image.open(DATA / f"Xenium_annotations/images/230918_Xenium_CytAssist_{g}.jpg").convert("RGB"))
            for g in sorted(set(GROUP.values()))}
    fig, axes = plt.subplots(3, 4, figsize=(15, 11.5))
    for a, s in zip(axes.ravel(), GROUP):
        x = t[t["sample"] == s]
        l = x.iloc[(x.landmark_err_um - x.landmark_err_um.median()).abs().argmin()]
        pts = np.array([[l.clicked_row, l.clicked_col], [l.truth_row, l.truth_col], [l.old_map_row, l.old_map_col]])
        ctr = pts.mean(0); half = max(170, 0.5 * np.ptp(pts, 0).max() + 90)       # JPG px
        r0, c0 = np.round(ctr - half).astype(int); r1, c1 = np.round(ctr + half).astype(int)
        a.imshow(imgs[GROUP[s]][max(r0, 0):r1, max(c0, 0):c1], extent=(max(c0, 0), c1, r1, max(r0, 0)))
        for g in range(int(np.ceil(c0 / 50)) * 50, c1, 50):                      # the notebooks' 50 px grid
            a.axvline(g, color="white", lw=0.4, alpha=0.6)
        for g in range(int(np.ceil(r0 / 50)) * 50, r1, 50):
            a.axhline(g, color="white", lw=0.4, alpha=0.6)
        a.annotate("", xy=(l.truth_col, l.truth_row), xytext=(l.clicked_col, l.clicked_row),
                   arrowprops=dict(arrowstyle="->", color="#d1242f", lw=1.6))
        a.scatter(l.clicked_col, l.clicked_row, marker="x", s=70, color="#d1242f", lw=2, label="clicked (pointsI)", zorder=4)
        a.scatter(l.truth_col, l.truth_row, s=60, color="#2da44e", edgecolor="white", lw=1, label="true position", zorder=4)
        a.scatter(l.old_map_col, l.old_map_row, s=70, facecolor="none", edgecolor="#1f5fae", lw=2, label="old final map", zorder=4)
        a.plot([c0 + 15, c0 + 15 + 100 / JPG_UM], [r1 - 15] * 2, color="k", lw=3)
        a.set_xlim(c0, c1); a.set_ylim(r1, r0); a.set_xticks([]); a.set_yticks([])
        a.set_title(f"{s} landmark {int(l.landmark)}: clicked {l.landmark_err_um:.0f} um off;\n"
                    f"old map {l.old_map_err_um:.0f} um off; nearest grid point to truth {l.truth_to_nearest_grid_um:.0f} um", fontsize=8.5)
    axes.ravel()[-1].axis("off")
    h, lab = axes.ravel()[0].get_legend_handles_labels()
    axes.ravel()[-1].legend(h, lab, loc="center", frameon=False, fontsize=10, title="white lines = 50 JPG px grid (88 um)")
    fig.suptitle("Each sample's median hand landmark on the H&E JPG it was clicked in. Arrow: clicked -> true position "
                 "of the paired Xenium point. Bar = 100 um.", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(OUT / "fig_old_landmarks_map.jpg", dpi=90, pil_kwargs={"quality": 85})
    plt.close(fig)


if __name__ == "__main__":
    main()
