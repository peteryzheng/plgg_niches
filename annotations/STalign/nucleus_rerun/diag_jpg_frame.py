#!/usr/bin/env python3
# Is the whole-slide JPG exactly the SVS at level 0 / 7, and do the pathology polygons sit in that frame?
#
# Every frame in this pipeline assumes it: the GeoJSON polygons (SVS level-0 px) are scaled by 1/7 to label
# cells, STAlign_x / _y are "JPG px", and the original notebooks registered the JPG image itself. If the JPG
# were shifted or rescaled against the SVS, the old map's error and the replays (replay_old_stalign.py) would
# partly reflect that mismatch rather than STAlign.
# Test: at 6 random cells of each sample (030226 positions, i.e. inside each section), a 192 JPG px patch
# (~340 um) of the JPG is compared with the same area of the SVS (page 2 = level 0 / 4, block-averaged and
# resampled onto the JPG grid with JPG px (i, j) <-> level-0 (7i + 3, 7j + 3), the centre of its 7x7 block).
# The offset comes from phase correlation at 1/20 px; corr = pixel correlation of the two patches.
# Polygons: per sample, the outlines of the notebook's GeoJSON polygons (scaled 1/7, as for labelling) near the
# section are compared with the tissue edge of the JPG (grey < 215 after a 5 px median; opened twice). Score =
# mean over outline points of min(distance to the nearest tissue edge, 10 px), for polygon shifts of +-40 JPG px;
# a frame offset would put the minimum away from zero, and the same offset would recur in every sample of a
# GeoJSON file. Hand-drawn outlines follow the edge only in places (share_near_edge = outline points within
# 15 px of an edge), so the test is informative only where that share is high.
# Output: NEW_DIR/qc/frame/jpg_vs_svs.tsv (one row per patch; offsets in JPG px and um),
#         NEW_DIR/qc/frame/polygon_edges.tsv (one row per sample).
# Usage: python diag_jpg_frame.py
import sys
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image
import geopandas as gpd
from scipy import ndimage as ndi
from shapely.geometry import Polygon, MultiPolygon, box
from skimage.registration import phase_cross_correlation

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nucleus_core import SVS, MPP0, JPG_TO_L0
from nucleus_windows import GROUP, DATA, meta_path
from stalign_nucleus import notebook_spec

Image.MAX_IMAGE_PIXELS = None
NEW_DIR = DATA / "Xenium_annotations/pathology_annotations_101326_nucleus"
OUT = NEW_DIR / "qc/frame"
IMG = DATA / "Xenium_annotations/images/230918_Xenium_CytAssist_"
P = 192                                     # patch side, JPG px
CAP = 10                                    # chamfer distance cap, JPG px (an outline far from any edge adds a constant)


def outline_points(geoms):
    """Points every ~1 JPG px along all rings (exteriors and holes) of the polygons, as (x = col, y = row)."""
    for g in geoms:
        for p in [g] if isinstance(g, Polygon) else list(g.geoms) if isinstance(g, MultiPolygon) else []:
            for r in [p.exterior, *p.interiors]:
                yield np.array([r.interpolate(t).coords[0] for t in np.linspace(0, r.length, max(int(r.length), 2), endpoint=False)])


def polygon_edges():
    rows = []
    for g in ("LGG1", "LGG2"):
        im = np.asarray(Image.open(f"{IMG}{g}.jpg").convert("L"))
        tissue = ndi.binary_opening(ndi.median_filter(im, 5) < 215, iterations=2)
        dist = ndi.distance_transform_edt(~(tissue ^ ndi.binary_erosion(tissue)))     # px to the nearest tissue edge
        for s in [k for k, v in GROUP.items() if v == g]:
            gj, _ = notebook_spec(s)
            p = gpd.read_file(gj); p = p[p.geometry.notna()]
            p.geometry = p.geometry.scale(1 / JPG_TO_L0, 1 / JPG_TO_L0, origin=(0, 0))
            # this section's area: its cells' (new map) bounding box + 100 JPG px; outline points inside it only
            m = pd.read_csv(meta_path(s, NEW_DIR), usecols=["STAlign_x", "STAlign_y"]).dropna().values
            (r0, c0), (r1, c1) = m.min(0) - 100, m.max(0) + 100
            sel = p[p.intersects(box(c0, r0, c1, r1))]
            pts = np.vstack(list(outline_points(sel.geometry)))
            pts = pts[(pts[:, 0] > c0) & (pts[:, 0] < c1) & (pts[:, 1] > r0) & (pts[:, 1] < r1)]
            score = lambda sh: np.minimum(ndi.map_coordinates(dist, [pts[:, 1] + sh[0], pts[:, 0] + sh[1]], order=1), CAP).mean()
            coarse = min(((dy, dx) for dy in range(-40, 41, 4) for dx in range(-40, 41, 4)), key=score)
            best = min([(coarse[0] + dy, coarse[1] + dx) for dy in range(-4, 5) for dx in range(-4, 5)] + [(0, 0)], key=score)
            rows.append(dict(sample=s, geojson=gj.name, n_polygons=len(sel), n_outline_pts=len(pts),
                             share_near_edge=(ndi.map_coordinates(dist, [pts[:, 1], pts[:, 0]], order=1) <= 15).mean(),
                             chamfer_at_zero=score((0, 0)), best_row_px=best[0], best_col_px=best[1], chamfer_best=score(best),
                             best_um=np.hypot(*best) * JPG_TO_L0 * MPP0))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "polygon_edges.tsv", sep="\t", index=False)
    print(t.round(3).to_string(index=False))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for g in ("LGG1", "LGG2"):
        jpg = np.asarray(Image.open(f"{IMG}{g}.jpg").convert("L"), np.float32)
        svs = SVS(f"{IMG}{g}.svs", page=2)
        down = svs.down                                        # level-0 px per page-2 px (4)
        ii, jj = np.mgrid[0:P, 0:P]
        for s in [k for k, v in GROUP.items() if v == g]:
            m = pd.read_csv(meta_path(s), usecols=["STAlign_x", "STAlign_y"]).dropna().values
            for r, c in m[np.random.default_rng(0).choice(len(m), 6, replace=False)].astype(int):
                r0, c0 = r - P // 2, c - P // 2
                a = jpg[r0:r0 + P, c0:c0 + P]
                # page-2 px coordinate of each JPG px centre; page-2 px k covers level 0 [4k, 4k + 4), centre 4k + 1.5
                pr = (JPG_TO_L0 * (r0 + ii) + 3 - (down - 1) / 2) / down
                pc = (JPG_TO_L0 * (c0 + jj) + 3 - (down - 1) / 2) / down
                lo = np.floor([pr.min(), pc.min()]).astype(int) - 4
                n = int(np.ceil(P * JPG_TO_L0 / down)) + 16             # page-2 px covering the patch + margin
                blk = svs.read(lo[0], lo[1], n, n).astype(np.float32).mean(2)
                blk = ndi.uniform_filter(blk, JPG_TO_L0 / down)  # ~ the JPG's 7x7 level-0 block average
                b = ndi.map_coordinates(blk, [pr - lo[0], pc - lo[1]], order=1)
                sh, _, _ = phase_cross_correlation(b, a, upsample_factor=20)
                rows.append(dict(slide=g, sample=s, jpg_row=r, jpg_col=c, d_row_px=sh[0], d_col_px=sh[1],
                                 d_um=np.hypot(*sh) * JPG_TO_L0 * MPP0, corr=np.corrcoef(a.ravel(), b.ravel())[0, 1]))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "jpg_vs_svs.tsv", sep="\t", index=False)
    print(t.groupby("sample")[["d_um", "corr"]].agg(["median", "max", "min"]).round(3).to_string())
    print(f"all patches: max offset {t.d_um.max():.2f} um, min corr {t['corr'].min():.3f}")
    polygon_edges()


if __name__ == "__main__":
    main()
