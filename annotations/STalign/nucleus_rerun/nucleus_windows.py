#!/usr/bin/env python3
# Nucleus-scale H&E <-> Xenium registration check, parameterised by the registration (meta) file.
#
# Purpose : measure, in 150 um windows, how far Xenium DAPI nuclei land from H&E nuclei when mapped
#           through a given per-cell registration (STAlign_x = H&E row, STAlign_y = H&E col, JPG px).
#           Used twice by stalign_nucleus.py: on the old 030226 map (to build landmarks + a baseline)
#           and on the new map (held-out validation), always at the SAME window centres.
# Inputs  : a proseg_pathology_<sample>_proseg_meta.csv (any folder), the SVS H&E, the Xenium
#           morphology_focus DAPI image.
# Output  : one row per window (see measure()); shift d is such that DAPI_mapped ~= H&E + d, in um,
#           as (row, col) of the H&E frame -- the convention of registration_nuclei_check.best_shift.
# Assumes : the H&E is a post-Xenium stain of the same section, so nuclei correspond one-to-one.
#
# The measurement primitives (nucleus peaks, FFT shift search, image readers) come from the CL lipid-map
# session's registration_nuclei_check.py / he_image.py, vendored verbatim into nucleus_core.py.
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
from skimage import color

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nucleus_core import (DAPI, peaks, best_shift, folder, GROUP, DATA, DAPI_UM, WIN_UM, SHIFT_UM, NUC_MIN_DIST_UM,
                          SVS, MPP0, JPG_TO_L0, F, UM)   # F = 2 level-0 px per H&E analysis px -> UM = 0.504 um/px

OLD_DIR = DATA / "Xenium_annotations/pathology_annotations_030226"
MIN_CELLS = 100      # a window must hold >= this many mapped cells (skips glass / sparse edges)
Z_CONF = 5.0         # FFT peak z-score at/above which a window's shift is trusted


def meta_path(s, root=OLD_DIR):
    return Path(root) / f"proseg_pathology_{s}_proseg_meta.csv"


def svs_for(s):
    return SVS(DATA / f"Xenium_annotations/images/230918_Xenium_CytAssist_{GROUP[s]}.svs")


def dapi_path(s):
    f = folder(s)
    return DATA / f"segmentation/proseg_run_121024/{f}/{f}_proseg/outs/morphology_focus/morphology_focus_0000.ome.tif"


class Mapping:
    """Xenium um (x, y) <-> H&E level-0 px (row, col), piecewise-linear over a meta file's cell pairs
    (cells are ~10 um apart, far finer than the errors being measured)."""
    def __init__(self, path):
        m = pd.read_csv(path, usecols=["x_centroid", "y_centroid", "STAlign_x", "STAlign_y"]).dropna()
        self.um_xy = m[["x_centroid", "y_centroid"]].values
        self.l0_rc = m[["STAlign_x", "STAlign_y"]].values * JPG_TO_L0
        self.to_l0 = LinearNDInterpolator(self.um_xy, self.l0_rc)
        self.to_um = LinearNDInterpolator(self.l0_rc, self.um_xy)


def draw_windows(mp, n, spacing_um, seed=0):
    """Up to n window centres (level-0 px) at random mapped cells whose window holds >= MIN_CELLS
    cells, kept >= spacing_um apart (Chebyshev) so windows are distinct tissue."""
    tree = cKDTree(mp.l0_rc)
    half0, sp0 = WIN_UM / 2 / MPP0, spacing_um / MPP0
    out = []
    for j in np.random.default_rng(seed).permutation(len(mp.l0_rc)):
        if len(out) >= n:
            break
        c = mp.l0_rc[j]
        if out and np.abs(np.asarray(out) - c).max(1).min() < sp0:
            continue
        if len(tree.query_ball_point(c, half0, p=np.inf)) < MIN_CELLS:
            continue
        out.append(c)
    return np.round(np.asarray(out)).astype(int)


def measure(s, mp, centres, svs, dapi):
    """Nucleus offset in each window centred at `centres` (level-0 row, col); the body of
    registration_nuclei_check.run_cells with the mapping passed in. Windows whose corners fall
    outside the mapped area, or with < 20 nuclei on either side, are reported with NaN shift."""
    half0 = int(WIN_UM / 2 / MPP0)
    rows = []
    for k, (cy, cx) in enumerate(centres):
        row = dict(sample=s, window=k, row_l0=cy, col_l0=cx)
        rgb = svs.read_scaled(cy - half0, cx - half0, 2 * half0, 2 * half0, F).astype(np.uint8)
        he = peaks(color.rgb2hed(rgb)[..., 0], 1.0 / UM, NUC_MIN_DIST_UM / UM) * UM   # um, window frame
        # DAPI read region = the window's corners mapped into Xenium um, padded by the search range
        corners = np.array([[cy - half0, cx - half0], [cy + half0, cx + half0], [cy - half0, cx + half0], [cy + half0, cx - half0]])
        cu = mp.to_um(corners)
        if np.isnan(cu).any():
            rows.append(row); continue
        x0, x1 = cu[:, 0].min() - SHIFT_UM, cu[:, 0].max() + SHIFT_UM
        y0, y1 = cu[:, 1].min() - SHIFT_UM, cu[:, 1].max() + SHIFT_UM
        r0, c0 = int(y0 / DAPI_UM), int(x0 / DAPI_UM)
        img = dapi.read(r0, c0, int((y1 - y0) / DAPI_UM), int((x1 - x0) / DAPI_UM))
        dp = peaks(img, 1.0 / DAPI_UM, NUC_MIN_DIST_UM / DAPI_UM)
        # DAPI nuclei -> H&E level 0 through the mapping under test -> um in the window frame
        d_l0 = mp.to_l0(np.c_[(dp[:, 1] + c0) * DAPI_UM, (dp[:, 0] + r0) * DAPI_UM])
        dw = (d_l0[~np.isnan(d_l0).any(1)] - [cy - half0, cx - half0]) * MPP0
        dw = dw[(dw >= -SHIFT_UM).all(1) & (dw <= WIN_UM + SHIFT_UM).all(1)]
        row.update(n_he=len(he), n_dapi=len(dw))
        if len(he) >= 20 and len(dw) >= 20:
            dy, dx, z, m0, mb, res = best_shift(he, dw)
            row.update(shift_row_um=dy, shift_col_um=dx, shift_um=np.hypot(dy, dx), peak_z=z,
                       match_at_zero=m0, match_at_best=mb, residual_um=res)
        rows.append(row)
        print(f"  {s} window {k}: " + (f"shift {row['shift_um']:.1f} um, z {row['peak_z']:.1f}"
              if "shift_um" in row else "skipped"), flush=True)
    t = pd.DataFrame(rows)
    for c in ["shift_row_um", "shift_col_um", "shift_um", "peak_z", "match_at_zero", "match_at_best", "residual_um"]:
        if c not in t:
            t[c] = np.nan
    return t
