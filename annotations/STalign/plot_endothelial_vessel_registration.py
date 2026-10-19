#!/usr/bin/env python
"""
STalign registration positive control (deck Edge 2): show that reverse-mapped
Xenium cells land on the right H&E structure.

Endothelial cells (cell type) for PA258482 are overlaid on the co-registered H&E
jpg at a vessel-dense ROI; they trace the wall of a visible blood vessel (lumen +
pooled RBCs). This is the qualitative control that the compact-vs-myxoid niche
null (fig_niche_pa_morphology) rests on: the alignment is accurate, so the null is
not a registration artifact.

Key idea: STalign was already applied per sample -- the `STAlign_x`/`STAlign_y`
columns in the pathology CSV ARE the cell centroids in H&E jpg PIXEL space
(STAlign_x = row, STAlign_y = column). No LDDMM re-run is needed to overlay.

Inputs (canonical, under WORKDIR = detect_workdir()/youyun/plgg/data):
  - Xenium_annotations/images/230918_Xenium_CytAssist_LGG1.jpg      (co-registered H&E)
  - Xenium_annotations/pathology_annotations_030226/proseg_pathology_PA258482_proseg_meta.csv
        (per-cell STAlign_x/STAlign_y + centroids, from the STalign pathology qmd)
  - xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv  (per-cell cell_type;
        cell_id is suffixed `<base>__<sample_id>`, sample_id is numeric)

Output:
  - <repo>/presentation/assets/fig_endothelial_vessel.png           (deck slide-14 inset)

Assumptions: PA258482 (numeric id 258482) is on the LGG1 image; the ROI center
below is "endo hotspot 3", the clearest isolated vessel found in a density sweep.
Run with the spatial env python, e.g.:
  /Users/youyun/miniforge3/envs/spatial/bin/python \
      annotations/STalign/plot_endothelial_vessel_registration.py
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

Image.MAX_IMAGE_PIXELS = None  # the CytAssist jpg is large; disable the decompression guard


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


WORKDIR = detect_workdir()
DATA = WORKDIR / "youyun/plgg/data"
REPO = Path(__file__).resolve().parents[2]                 # annotations/STalign/ -> repo root
ASSETS = REPO / "presentation" / "assets"

# --- sample / ROI configuration ---
SAMPLE_NUM = "258482"          # PA258482, numeric id used in the metadata TSV
IMAGE_GROUP = "LGG1"           # PA258482 sits on the LGG1 CytAssist image
CR, CC, HALF = 4992, 7405, 110  # ROI center in jpg px (row=STAlign_x, col=STAlign_y); "endo hotspot 3"
CELL_COLOR = "#1f6feb"         # blue hollow rings so the H&E shows through each marker

JPG = DATA / "Xenium_annotations/images" / f"230918_Xenium_CytAssist_{IMAGE_GROUP}.jpg"
STALIGN_CSV = (DATA / "Xenium_annotations/pathology_annotations_030226"
               / f"proseg_pathology_PA{SAMPLE_NUM}_proseg_meta.csv")
META_TSV = DATA / "xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv"

# --- load H&E + per-cell STalign coords ---
V = np.asarray(Image.open(JPG).convert("RGB"))
H, W = V.shape[:2]
st = pd.read_csv(STALIGN_CSV)  # has cell_id, x_centroid/y_centroid (µm), STAlign_x/y (jpg px)

# --- per-cell cell type from the metadata TSV (chunked + filtered: the file is ~200MB) ---
# TSV cell_id is `<base>__<sample_id>`; strip the suffix to match the STalign CSV's `<base>`.
meta = pd.concat(
    chunk[chunk["sample_id"].astype(str) == SAMPLE_NUM]
    for chunk in pd.read_csv(META_TSV, sep="\t", usecols=["cell_id", "sample_id", "cell_type"],
                             dtype=str, chunksize=500_000)
)
meta["cell_id"] = meta["cell_id"].str.split("__").str[0]

df = st.merge(meta[["cell_id", "cell_type"]], on="cell_id", how="left")
df = df[df.STAlign_x.between(0, H - 1) & df.STAlign_y.between(0, W - 1)]
endo = df[df.cell_type == "endothelial"]              # the actual vessel-lining cells
print(f"coverage: {df.cell_type.notna().mean()*100:.1f}% cells typed | endothelial: {len(endo)}")

# --- empirical px/µm scale for the scale bar (global affine fit; LDDMM is near-rigid) ---
A = np.c_[df.x_centroid, df.y_centroid, np.ones(len(df))]
bx, *_ = np.linalg.lstsq(A, df.STAlign_x, rcond=None)
by, *_ = np.linalg.lstsq(A, df.STAlign_y, rcond=None)
px_per_um = float(np.sqrt(abs(np.linalg.det(np.array([[bx[0], bx[1]], [by[0], by[1]]])))))

# --- crop the ROI + the endothelial cells inside it ---
rr0, rr1, cc0, cc1 = CR - HALF, CR + HALF, CC - HALF, CC + HALF
crop = V[rr0:rr1, cc0:cc1]
s = endo[endo.STAlign_x.between(rr0, rr1) & endo.STAlign_y.between(cc0, cc1)]
print(f"ROI {2*HALF}px (~{2*HALF/px_per_um:.0f} µm) | endothelial cells in ROI: {len(s)}")

# --- figure: hollow blue rings on the H&E crop + a 100 µm scale bar ---
bar_um, side = 100, 2 * HALF
fig, ax = plt.subplots(figsize=(5.2, 5.2))
ax.imshow(crop)
ax.scatter(s.STAlign_y - cc0, s.STAlign_x - rr0, s=32, facecolors="none",
           edgecolors=CELL_COLOR, linewidths=1.1)                      # cell overlay
ax.add_patch(Rectangle((0.06 * side, 0.93 * side), bar_um * px_per_um, 4, color="black"))
ax.text(0.06 * side + bar_um * px_per_um / 2, 0.93 * side - 6, f"{bar_um} µm",
        ha="center", va="bottom", fontsize=9, fontweight="bold")       # scale bar
ax.scatter([], [], s=44, facecolors="none", edgecolors=CELL_COLOR, linewidths=1.1,
           label="Endothelial cell")
ax.legend(loc="upper right", fontsize=8, framealpha=0.9, markerscale=1.1)
ax.set_title("Endothelial cells map onto a blood vessel\nPA258482 H&E (STalign registration)",
             fontsize=11)
ax.axis("off")
fig.tight_layout()
out = ASSETS / "fig_endothelial_vessel.png"
fig.savefig(out, dpi=220, bbox_inches="tight")
print(f"wrote {out}")
