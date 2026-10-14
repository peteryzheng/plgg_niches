#!/usr/bin/env python3
# How much the nucleus-guided registration changes each cell's H&E position and pathology label.
#
# Compares, cell by cell (same cell_id), the old 030226 meta with the new meta (and optionally an
# alternative variant such as the affine-only trial):
#   label_impact.tsv    : per sample -- cells, H&E displacement old -> new (median / p90, um), share of
#                         cells whose annotation changed old -> new, over all cells and over the cells the old
#                         map placed inside an annotated polygon (and new vs the alternative)
#   label_retention.tsv : per sample x old annotation -- cells with that label under the old map, share
#                         that keep it under the new map, cells with it under the new map, and the share
#                         of those that are new entrants (had another label before)
# Usage : python label_impact.py <new meta dir> <out dir> [alternative meta dir]
# Assumes: metas share cell_id order-free; STAlign_x/_y are JPG px (7 level-0 px, MPP0 um each).
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nucleus_windows import GROUP, meta_path, MPP0, JPG_TO_L0

JPG_UM = JPG_TO_L0 * MPP0          # um per JPG px (~1.77)
COLS = ["cell_id", "STAlign_x", "STAlign_y", "annotation"]

new_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
alt_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else None
out_dir.mkdir(parents=True, exist_ok=True)

impact, keep = [], []
for s in GROUP:
    newp = meta_path(s, new_dir)
    if not newp.exists():
        continue
    o = pd.read_csv(meta_path(s), usecols=COLS).set_index("cell_id")
    n = pd.read_csv(newp, usecols=COLS).set_index("cell_id").loc[o.index]
    # cells the old map never placed (NaN coordinates) carry no label to compare
    ok = o.STAlign_x.notna() & n.STAlign_x.notna()
    o, n = o[ok], n[ok]
    d = np.hypot(n.STAlign_x - o.STAlign_x, n.STAlign_y - o.STAlign_y) * JPG_UM
    r = dict(sample=s, n_cells=len(o), disp_median_um=d.median(), disp_p90_um=d.quantile(0.9),
             label_changed_vs_old=(n.annotation != o.annotation).mean())
    # annotated-only rate (the figure stalign_nucleus.py logs): small denominators in sparsely annotated sections
    ann = o.annotation != "unannotated"
    r.update(n_annotated_old=int(ann.sum()), label_changed_annotated=(n.annotation != o.annotation)[ann].mean())
    if alt_dir is not None and meta_path(s, alt_dir).exists():
        a = pd.read_csv(meta_path(s, alt_dir), usecols=["cell_id", "annotation"]).set_index("cell_id").loc[o.index]
        r["label_changed_vs_alt"] = (n.annotation != a.annotation).mean()
    impact.append(r)
    # per old label: retention, and how much of the new label population is newly recruited
    for lab, g in o.groupby("annotation"):
        new_lab = n.annotation == lab
        keep.append(dict(sample=s, annotation=lab, n_old=len(g), kept=(n.annotation[g.index] == lab).mean(),
                         n_new=int(new_lab.sum()), new_entrants=(o.annotation[new_lab] != lab).mean() if new_lab.any() else np.nan,
                         disp_median_um=d[g.index].median()))

pd.DataFrame(impact).to_csv(out_dir / "label_impact.tsv", sep="\t", index=False)
pd.DataFrame(keep).to_csv(out_dir / "label_retention.tsv", sep="\t", index=False)
print(pd.DataFrame(impact).round(3).to_string(index=False))
