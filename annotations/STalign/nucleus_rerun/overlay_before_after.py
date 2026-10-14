#!/usr/bin/env python3
# Before/after picture of the H&E <-> Xenium registration at nucleus scale.
#
# For three HELD-OUT windows (never used as landmarks; picked at the 25th/50th/75th percentile of the
# old offset), shows the H&E crop with Xenium DAPI nuclei mapped on top through the old 030226 map
# (left) and the new map (right). A correct map puts a marker on every hematoxylin nucleus.
# Usage : python overlay_before_after.py <sample> <new meta csv> <windows_after tsv> <out.jpg>
# Output: a JPG of patient H&E -- share only via SendUserFile / Notion (README standing rules).
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nucleus_windows import DATA, Mapping, meta_path, svs_for, dapi_path, DAPI, peaks, MPP0, DAPI_UM, WIN_UM, NUC_MIN_DIST_UM, F, UM

s, new_meta, win_tsv, out = sys.argv[1:5]
w = pd.read_csv(win_tsv, sep="\t")
# old-map measurement at the same centres (rows align with windows_after)
before = pd.read_csv(DATA / "Xenium_annotations/pathology_annotations_101326_nucleus/qc" / f"windows_before_{s}.tsv", sep="\t")
# held-out windows ordered by their old offset; take the quartile windows so the panel is representative
ho = w[w.role == "heldout"].assign(old_shift=before.loc[w.role == "heldout", "shift_um"].values).sort_values("old_shift")
pick = ho.iloc[[len(ho) // 4, len(ho) // 2, 3 * len(ho) // 4]]

maps = {"old map (030226)": Mapping(meta_path(s)), "new map (nucleus landmarks + LDDMM)": Mapping(new_meta)}
svs, dapi = svs_for(s), DAPI(dapi_path(s))
half0 = int(WIN_UM / 2 / MPP0)
fig, axes = plt.subplots(len(pick), 2, figsize=(10, 5 * len(pick)))
for i, (_, r) in enumerate(pick.iterrows()):
    cy, cx = int(r.row_l0), int(r.col_l0)
    rgb = svs.read_scaled(cy - half0, cx - half0, 2 * half0, 2 * half0, F).astype(np.uint8)
    for j, (name, mp) in enumerate(maps.items()):
        # DAPI nuclei from the Xenium area this map says lies under the window (+10 um border)
        cu = mp.to_um(np.array([[cy - half0, cx - half0], [cy + half0, cx + half0], [cy - half0, cx + half0], [cy + half0, cx - half0]]))
        x0, x1, y0, y1 = cu[:, 0].min() - 10, cu[:, 0].max() + 10, cu[:, 1].min() - 10, cu[:, 1].max() + 10
        r0, c0 = int(y0 / DAPI_UM), int(x0 / DAPI_UM)
        dp = peaks(dapi.read(r0, c0, int((y1 - y0) / DAPI_UM), int((x1 - x0) / DAPI_UM)), 1.0 / DAPI_UM, NUC_MIN_DIST_UM / DAPI_UM)
        d = (mp.to_l0(np.c_[(dp[:, 1] + c0) * DAPI_UM, (dp[:, 0] + r0) * DAPI_UM]) - [cy - half0, cx - half0]) * MPP0 / UM   # window px
        ax = axes[i, j]
        ax.imshow(rgb, extent=(0, rgb.shape[1], rgb.shape[0], 0))
        ax.scatter(d[:, 1], d[:, 0], s=14, facecolors="none", edgecolors="lime", linewidths=0.8)
        ax.set_xlim(0, rgb.shape[1]); ax.set_ylim(rgb.shape[0], 0); ax.set_xticks([]); ax.set_yticks([])
        sh = r.old_shift if j == 0 else r.shift_um
        ax.set_title(f"{name}\nwindow {int(r.window)}: measured offset {sh:.1f} um", fontsize=9)
    axes[i, 0].plot([10, 10 + 50 / UM], [rgb.shape[0] - 10] * 2, color="k", lw=3)   # 50 um scale bar
fig.suptitle(f"{s}: Xenium DAPI nuclei (green circles) on H&E, held-out windows (150 um)", fontsize=11)
fig.tight_layout(); fig.savefig(out, dpi=110, pil_kwargs={"quality": 85})
print("wrote", out)
