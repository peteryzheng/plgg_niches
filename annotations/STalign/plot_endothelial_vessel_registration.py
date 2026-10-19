#!/usr/bin/env python
"""
STalign registration figure (deck slide 14 inset): the same endothelial cells placed on the PA258482 H&E by
the original map and by the nucleus-guided rerun, with an arrow from each cell's old position to its new one.

Endothelial cells (cell type) are drawn on a full-resolution SVS crop of a vessel-dense ROI. Both maps put them
on the vessel, but the nucleus map moves every cell by ~20 µm in this ROI. A vessel overlay alone therefore
cannot tell a 1 µm map from a 20 µm one; the accuracy evidence is the held-out nucleus check in
annotations/STalign/nucleus_rerun/ (~1 µm for the new map).

Frames: the meta files' `STAlign_x`/`STAlign_y` are cell centroids in H&E JPG pixels (STAlign_x = row,
STAlign_y = column). The JPG is the SVS level 0 downsampled exactly 7x (checked by nucleus_rerun/diag_jpg_frame.py),
so level-0 px = JPG px * 7, at MPP0 = 0.2522 µm per level-0 px.

Inputs (under WORKDIR = detect_workdir()/youyun/plgg/data):
  - Xenium_annotations/images/230918_Xenium_CytAssist_LGG1.svs              (H&E, read at level 0)
  - Xenium_annotations/pathology_annotations_030226/proseg_pathology_PA258482_proseg_meta.csv        (old map)
  - Xenium_annotations/pathology_annotations_101326_nucleus/proseg_pathology_PA258482_proseg_meta.csv (nucleus map)
  - xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826_pathology_nucleus.tsv  (per-cell cell_type;
        cell_id is `<base>__<sample_id>`, sample_id numeric)

Output:
  - <repo>/presentation/assets/fig_endothelial_vessel.png

Assumptions: PA258482 (numeric id 258482) is on the LGG1 slide; the ROI center below is "endo hotspot 3", the
clearest isolated vessel found in a density sweep (chosen on the JPG, so given in JPG px). Patient H&E: share the
output only through the deck, SendUserFile or Notion.
Run with the spatial env python:
  python annotations/STalign/plot_endothelial_vessel_registration.py
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REPO = Path(__file__).resolve().parents[2]                 # annotations/STalign/ -> repo root
sys.dont_write_bytecode = True                              # keep __pycache__ out of the nucleus_rerun folder
sys.path.insert(0, str(REPO / "annotations/STalign/nucleus_rerun"))
from nucleus_core import SVS, MPP0, JPG_TO_L0              # tiled level-0 SVS reader and the frame constants


# path resolution: derive the shared-data root from HOME (AGENTS.md convention)
def detect_workdir() -> Path:
    home = os.path.expanduser("~")
    if home in ("/Users/youyun", "/Users/youyunzheng"):
        return Path(home) / "Documents/HMS/PhD/beroukhimlab/dfci_mount"
    if home == "/home/yz762":
        return Path("/mnt/storage/dept/medonc/beroukhim")
    if home == "/PHShome/yz762":
        return Path("/data/beroukhim1")
    return Path("/xchip/beroukhimlab")


DATA = detect_workdir() / "youyun/plgg/data"
ASSETS = REPO / "presentation" / "assets"

# --- sample / ROI configuration ---
SAMPLE_NUM = "258482"           # PA258482, numeric id used in the metadata TSV
IMAGE_GROUP = "LGG1"            # PA258482 sits on the LGG1 CytAssist slide
CR, CC, HALF = 4992, 7405, 110  # ROI center and half-width in JPG px (row=STAlign_x, col=STAlign_y); "endo hotspot 3"
OLD_COLOR, NEW_COLOR = "#f28e2b", "#1f6feb"  # hollow rings so the H&E shows through each marker

SVS_PATH = DATA / "Xenium_annotations/images" / f"230918_Xenium_CytAssist_{IMAGE_GROUP}.svs"
META = "Xenium_annotations/{}/proseg_pathology_PA" + SAMPLE_NUM + "_proseg_meta.csv"
CELL_TSV = DATA / "xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826_pathology_nucleus.tsv"

# --- per-cell positions on both maps, in level-0 px ---
cols = ["cell_id", "STAlign_x", "STAlign_y"]
old = pd.read_csv(DATA / META.format("pathology_annotations_030226"), usecols=cols)
new = pd.read_csv(DATA / META.format("pathology_annotations_101326_nucleus"), usecols=cols)
df = old.merge(new, on="cell_id", suffixes=("_old", "_new"))
for m in ("old", "new"):
    df[f"r_{m}"] = df[f"STAlign_x_{m}"] * JPG_TO_L0
    df[f"c_{m}"] = df[f"STAlign_y_{m}"] * JPG_TO_L0
df["shift_um"] = np.hypot(df.r_new - df.r_old, df.c_new - df.c_old) * MPP0

# --- cell type from the metadata TSV (chunked + filtered: the file is ~200 MB) ---
# TSV cell_id is `<base>__<sample_id>`; strip the suffix to match the meta files' `<base>`.
meta = pd.concat(
    chunk[chunk["sample_id"].astype(str) == SAMPLE_NUM]
    for chunk in pd.read_csv(CELL_TSV, sep="\t", usecols=["cell_id", "sample_id", "cell_type"],
                             dtype=str, chunksize=500_000)
)
meta["cell_id"] = meta["cell_id"].str.split("__").str[0]
df = df.merge(meta[["cell_id", "cell_type"]], on="cell_id", how="left")

# --- ROI crop at level 0 + the endothelial cells with either position inside it ---
r0, c0, side = (CR - HALF) * JPG_TO_L0, (CC - HALF) * JPG_TO_L0, 2 * HALF * JPG_TO_L0
crop = SVS(SVS_PATH).read(r0, c0, side, side)
inside = lambda r, c: r.between(r0, r0 + side) & c.between(c0, c0 + side)
endo = df[(df.cell_type == "endothelial") & (inside(df.r_old, df.c_old) | inside(df.r_new, df.c_new))]
print(f"ROI {side} level-0 px (~{side * MPP0:.0f} µm) | endothelial cells: {len(endo)} | "
      f"median shift {endo.shift_um.median():.1f} µm (whole sample {df.shift_um.median():.1f} µm)")

# --- figure: arrows old -> new, hollow rings on both ends, 50 µm scale bar ---
fig, ax = plt.subplots(figsize=(7.5, 7.5))
ax.imshow(crop)
ax.quiver(endo.c_old - c0, endo.r_old - r0, endo.c_new - endo.c_old, endo.r_new - endo.r_old,
          angles="xy", scale_units="xy", scale=1, color="white", width=0.0025, headwidth=4, alpha=0.9)
ax.scatter(endo.c_old - c0, endo.r_old - r0, s=60, facecolors="none", edgecolors=OLD_COLOR, linewidths=1.3,
           label="Original STalign map")
ax.scatter(endo.c_new - c0, endo.r_new - r0, s=60, facecolors="none", edgecolors=NEW_COLOR, linewidths=1.3,
           label="Nucleus-guided map")
bar_px = 50 / MPP0
ax.add_patch(Rectangle((0.05 * side, 0.94 * side), bar_px, 0.008 * side, color="black"))
ax.text(0.05 * side + bar_px / 2, 0.935 * side, "50 µm", ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.set_xlim(0, side); ax.set_ylim(side, 0)                 # cells near the edge may point outside the crop
ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
ax.set_title(f"Endothelial cells on the PA258482 H&E, original vs nucleus-guided map\n"
             f"arrows old → new; median shift {endo.shift_um.median():.0f} µm here, "
             f"{df.shift_um.median():.0f} µm across the sample", fontsize=11)
ax.axis("off")
fig.tight_layout()
out = ASSETS / "fig_endothelial_vessel.png"
fig.savefig(out, dpi=220, bbox_inches="tight")
print(f"wrote {out}")
