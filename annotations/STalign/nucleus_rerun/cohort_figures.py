#!/usr/bin/env python3
# Cohort figures for the nucleus-guided STAlign rerun (documentation / Notion).
#
#   fig_offsets.png      : per sample, held-out windows only -- median nucleus offset (log scale) and the
#                          share of H&E nuclei with a DAPI partner within 2.5 um, old map -> new map;
#                          the affine-only variant is drawn as a reference marker, the 5 um target and
#                          the per-sample matching ceiling (best possible shift) as ticks.
#   fig_error_field.png  : per sample, how far each cell's H&E position moves old map -> new map (= the old
#                          map's error, since the new map is verified to ~1 um), median per 200 um square on
#                          a shared log colour scale, with direction arrows (x3) on an 800 um grid. Uses every
#                          cell, so it also shows regions the window check cannot measure on the old map
#                          (errors beyond its +-150 um search come out as unconfident windows).
#   fig_overlay_<slide>.jpg : one held-out window per sample (the median old offset), H&E with Xenium
#                          DAPI nuclei mapped through the old map (left) and the new map (right).
# Usage : python cohort_figures.py <new meta dir> <qc dir> <out dir> [affine-only qc dir]
# Output: PNG/JPG figures; the overlays contain patient H&E -> share only via SendUserFile / Notion.
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nucleus_windows import (DATA, Mapping, meta_path, svs_for, dapi_path, DAPI, peaks, GROUP, Z_CONF, MPP0, DAPI_UM,
                             WIN_UM, NUC_MIN_DIST_UM, F, UM, JPG_TO_L0)

# old-map window measurements (shared by every variant; written once by stalign_nucleus.py)
BEFORE = DATA / "Xenium_annotations/pathology_annotations_101326_nucleus/qc"

meta_dir, qc_dir, out_dir = (Path(a) for a in sys.argv[1:4])
aff_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else None
out_dir.mkdir(parents=True, exist_ok=True)
SAMPLES = [s for s in GROUP if (qc_dir / f"windows_after_{s}.tsv").exists()]
# one hue, two shades for before -> after (dataviz reference palette, slot 1 blue); text in neutral ink
OLD, NEW, REF, INK, MUTED = "#9cc3ef", "#1f5fae", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": MUTED,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})


def held(qc, s):
    """Held-out windows: old-map and new-map rows at the same centres."""
    b = pd.read_csv(BEFORE / f"windows_before_{s}.tsv", sep="\t")
    a = pd.read_csv(qc / f"windows_after_{s}.tsv", sep="\t")
    m = (a.role == "heldout").values
    return b[m].reset_index(drop=True), a[m].reset_index(drop=True)


def stats(t):
    c = t[t.peak_z >= Z_CONF]
    return c.shift_um.median(), t.match_at_zero.median(), t.match_at_best.median()


# ---------------------------------------------------------------- fig_offsets
rows = []
for s in SAMPLES:
    b, a = held(qc_dir, s)
    r = dict(sample=s, old=stats(b)[0], new=stats(a)[0], m_old=stats(b)[1], m_new=stats(a)[1], ceiling=stats(a)[2])
    # windows unmeasurable on the old map (z < Z_CONF; never landmarks): how many, and the new map's offset there
    u = pd.read_csv(qc_dir / f"windows_after_{s}.tsv", sep="\t").query("role == 'unconfident'")
    r.update(n_unmeasurable_old=len(u), n_confident_new=int((u.peak_z >= Z_CONF).sum()),
             unmeasurable_new_median=u[u.peak_z >= Z_CONF].shift_um.median())
    if aff_dir is not None and (aff_dir / f"windows_after_{s}.tsv").exists():
        _, aa = held(aff_dir, s); r.update(aff=stats(aa)[0], m_aff=stats(aa)[1])
    rows.append(r)
d = pd.DataFrame(rows).iloc[::-1].reset_index(drop=True)          # first sample at the top
d.to_csv(out_dir / "fig_offsets_data.tsv", sep="\t", index=False)
fig, ax = plt.subplots(1, 2, figsize=(10, 0.42 * len(d) + 1.6), sharey=True)
y = np.arange(len(d))
for k, (o, n, f, lab, lim) in enumerate([("old", "new", "aff", "median nucleus offset, held-out windows (um, log)", (0.3, 200)),
                                         ("m_old", "m_new", "m_aff", "H&E nuclei with a DAPI nucleus within 2.5 um", (0, 1))]):
    a_ = ax[k]
    a_.hlines(y, d[o], d[n], color="#d9d8d2", lw=2, zorder=1)
    a_.scatter(d[o], y, s=60, color=OLD, edgecolor="white", lw=2, zorder=3, label="old map (030226)")
    a_.scatter(d[n], y, s=60, color=NEW, edgecolor="white", lw=2, zorder=4, label="new map (landmark affine + LDDMM)")
    if f in d:
        a_.scatter(d[f], y, s=36, facecolor="none", edgecolor=REF, lw=1.2, zorder=5, label="affine only (reference)")
    if k == 0:
        a_.set_xscale("log"); a_.axvline(5, color=REF, lw=1, ls=(0, (3, 3)), zorder=0)
        a_.text(5, len(d) - 0.4, " 5 um target", color=MUTED, fontsize=8, va="bottom")
        for yi, v in zip(y, d[n]):                                     # above the dot: the affine marker sits beside it
            a_.text(v, yi + 0.2, f"{v:.1f}", color=INK, fontsize=8, ha="center", va="bottom")
        for yi, v in zip(y, d[o]):
            a_.text(v * 1.25, yi, f"{v:.0f}", color=MUTED, fontsize=8, ha="left", va="center")
    else:
        a_.scatter(d["ceiling"], y, marker="|", s=120, color=INK, lw=1.2, zorder=2, label="ceiling (best shift)")
    a_.set_xlim(*lim); a_.set_xlabel(lab); a_.grid(axis="x", color="#ecebe6", lw=0.8); a_.set_axisbelow(True)
ax[0].set_yticks(y, d["sample"]); ax[0].tick_params(axis="y", length=0, labelcolor=INK)
h, l = ax[1].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=len(l), frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(out_dir / "fig_offsets.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- fig_error_field
from matplotlib.colors import LogNorm
L0_MM = JPG_TO_L0 * MPP0 / 1000                                            # mm per JPG px
FINE, COARSE, EXAG = 0.2, 0.8, 3                                           # mm, mm, arrow exaggeration
n = len(SAMPLES); nc = 4; nr = -(-n // nc)
fig, axes = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.3 * nr))
norm = LogNorm(1, 600)
for a_, s in zip(axes.ravel(), SAMPLES):
    cols = ["cell_id", "STAlign_x", "STAlign_y"]
    o = pd.read_csv(meta_path(s), usecols=cols).set_index("cell_id")
    m = pd.read_csv(meta_dir / f"proseg_pathology_{s}_proseg_meta.csv", usecols=cols).set_index("cell_id").loc[o.index]
    ok = o.STAlign_x.notna() & m.STAlign_x.notna()
    # position = new-map H&E (mm); displacement = old minus new (mm), i.e. where the old map put the cell
    c = pd.DataFrame(dict(r=m.STAlign_x[ok] * L0_MM, c=m.STAlign_y[ok] * L0_MM,
                          dr=(o.STAlign_x - m.STAlign_x)[ok] * L0_MM, dc=(o.STAlign_y - m.STAlign_y)[ok] * L0_MM))
    c["d_um"] = np.hypot(c.dr, c.dc) * 1000
    g = c.groupby([(c.r // FINE).astype(int), (c.c // FINE).astype(int)]).d_um.median()
    ri, ci = g.index.get_level_values(0), g.index.get_level_values(1)
    img = np.full((ri.max() - ri.min() + 1, ci.max() - ci.min() + 1), np.nan)
    img[ri - ri.min(), ci - ci.min()] = g.values
    im = a_.imshow(img, cmap="Blues", norm=norm, interpolation="nearest",
                   extent=(ci.min() * FINE, (ci.max() + 1) * FINE, (ri.max() + 1) * FINE, ri.min() * FINE))
    q = c.groupby([(c.r // COARSE).astype(int), (c.c // COARSE).astype(int)]).filter(lambda x: len(x) >= 200) \
         .groupby([(c.r // COARSE).astype(int), (c.c // COARSE).astype(int)])[["r", "c", "dr", "dc"]].median()
    a_.quiver(q.c, q.r, q.dc * EXAG, q.dr * EXAG, angles="xy", scale_units="xy", scale=1, color=INK, width=0.006,
              headwidth=3.5, headlength=4)
    a_.set_aspect("equal"); a_.set_title(f"{s}  (median {c.d_um.median():.0f} um)", fontsize=9, color=INK)
    a_.tick_params(labelsize=7)
for a_ in axes.ravel()[n:]:
    a_.axis("off")
fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5, label="old-map error (um, median per 200 um square)")
fig.suptitle(f"Where the old map (030226) put each cell vs the new map; arrows x{EXAG} point to the old position; "
             "axes = H&E position (mm)", fontsize=10, color=INK)
fig.savefig(out_dir / "fig_error_field.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- fig_overlay_<slide>
half0 = int(WIN_UM / 2 / MPP0)
for slide in sorted(set(GROUP[s] for s in SAMPLES)):
    ss = [s for s in SAMPLES if GROUP[s] == slide]
    fig, axes = plt.subplots(len(ss), 2, figsize=(7, 3.6 * len(ss)), squeeze=False)
    for i, s in enumerate(ss):
        b, a = held(qc_dir, s)
        j = int(np.argsort(b.shift_um.values)[len(b) // 2])            # median old offset among held-out windows
        cy, cx = int(a.row_l0[j]), int(a.col_l0[j])
        svs, dapi = svs_for(s), DAPI(dapi_path(s))
        rgb = svs.read_scaled(cy - half0, cx - half0, 2 * half0, 2 * half0, F).astype(np.uint8)
        for k, (name, mp, sh) in enumerate([("old map", Mapping(meta_path(s)), b.shift_um[j]),
                                            ("new map", Mapping(meta_dir / f"proseg_pathology_{s}_proseg_meta.csv"), a.shift_um[j])]):
            cu = mp.to_um(np.array([[cy - half0, cx - half0], [cy + half0, cx + half0], [cy - half0, cx + half0], [cy + half0, cx - half0]]))
            x0, x1, y0, y1 = cu[:, 0].min() - 10, cu[:, 0].max() + 10, cu[:, 1].min() - 10, cu[:, 1].max() + 10
            r0, c0 = int(y0 / DAPI_UM), int(x0 / DAPI_UM)
            dp = peaks(dapi.read(r0, c0, int((y1 - y0) / DAPI_UM), int((x1 - x0) / DAPI_UM)), 1.0 / DAPI_UM, NUC_MIN_DIST_UM / DAPI_UM)
            p = (mp.to_l0(np.c_[(dp[:, 1] + c0) * DAPI_UM, (dp[:, 0] + r0) * DAPI_UM]) - [cy - half0, cx - half0]) * MPP0 / UM
            ax_ = axes[i, k]
            ax_.imshow(rgb, extent=(0, rgb.shape[1], rgb.shape[0], 0))
            ax_.scatter(p[:, 1], p[:, 0], s=12, facecolors="none", edgecolors="lime", linewidths=0.8)
            ax_.set_xlim(0, rgb.shape[1]); ax_.set_ylim(rgb.shape[0], 0); ax_.set_xticks([]); ax_.set_yticks([])
            ax_.set_title(f"{s}, {name}: offset {sh:.1f} um", fontsize=9)
        axes[i, 0].plot([8, 8 + 50 / UM], [rgb.shape[0] - 8] * 2, color="k", lw=3)   # 50 um bar
    fig.suptitle(f"{slide}: Xenium DAPI nuclei (green) on H&E, one held-out 150 um window per sample", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig(out_dir / f"fig_overlay_{slide}.jpg", dpi=100, pil_kwargs={"quality": 85})
    plt.close(fig)
print("wrote figures to", out_dir)
