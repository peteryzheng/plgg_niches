#!/usr/bin/env python3
# Measured-panel pseudobulk builder -- the ground-truth arm of the imputation QC.
#
# Purpose: ENVI imputation attenuates within-sample regional differences. To quantify by
# how much, we recompute the SAME region contrast on the 266 genes that are actually
# MEASURED (adata.X), over identical cells and identical sample x cell_type x region
# groups, and compare against the imputed result gene-for-gene.
#
# SCALE NOTE (why this is not run through voom like sample_comparisons/dlgnt12_*):
#   adata.X is NOT raw counts. Per-cell totals are exactly constant (65.4 in every cell of
#   every sample), i.e. the matrix is library-size normalized linear expression. There are
#   no `layers` and no `raw` group, so counts are not recoverable from these objects.
#   voom/TMM require genuine library sizes, so they do not apply. Instead the measured arm
#   gets the IDENTICAL treatment as the imputed arm (expression filter -> log2 -> limma
#   trend), which is what makes the two directly comparable: any residual difference is
#   attributable to imputation, not to the statistical method.
#
# Outputs (per cohort, consumed by measured_panel_fit.R):
#   measured_<COH>_pseudobulk_matrix.csv   gene x (sample__celltype__region), mean of X
#   measured_<COH>_sample_design.csv       column_id, sample, cell_type, region
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd, h5py

ROOT = Path("/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data")
ENVI = ROOT / "sc_integration/ENVI_results_091126"
OUT = Path(__file__).resolve().parent
MIN_CELLS = 20   # must match morphology_state_prepare.py

COHORTS = {
    "CL": dict(glob="*liponeurocytoma_spatial_envi.h5ad",
               regmap={"unannotated": "nonlipid", "Lipidized cells": "lipid"},
               cts=("neuronal", "astro_glial", "opc", "radial_glia", "oligo_mature", "myeloid")),
    "PA": dict(glob="*astrocytoma_spatial_envi.h5ad",
               regmap={"Compact. fibrillary component": "compact", "Loose, myxoid component": "loose"},
               cts=("opc", "astro_glial", "myeloid")),
}

def _dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

def main():
    meta = pd.read_csv(
        # pathology labels from the nucleus STAlign rerun (pathology_annotations_101326_nucleus)
        ROOT / "xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826_pathology_nucleus.tsv",
        sep="\t", usecols=["cell_id", "cell_type", "pathology_annotation"]).set_index("cell_id")

    # optional cohort names on the command line (e.g. `PA`) rebuild only those cohorts; the
    # other cohort's files stay as they are, so its measured/imputed pairing is unchanged.
    cohorts = {k: COHORTS[k] for k in sys.argv[1:]} if len(sys.argv) > 1 else COHORTS
    for coh, cfg in cohorts.items():
        cols, design = {}, []
        for f in sorted(glob.glob(str(ENVI / cfg["glob"]))):
            sid = Path(f).name.replace("_spatial_envi.h5ad", "")
            sample = sid.replace("_", " ")            # matches the imputed pipeline's column ids
            with h5py.File(f, "r") as h:
                ix = h["obs"].attrs.get("_index", "_index")
                ix = ix.decode() if isinstance(ix, bytes) else ix
                obs = np.array(_dec(h["obs"][ix][:]))
                vx = h["var"].attrs.get("_index", "_index")
                vx = vx.decode() if isinstance(vx, bytes) else vx
                genes = np.array([g.upper() for g in _dec(h["var"][vx][:])])
                ct = pd.Index(obs).map(meta["cell_type"]).to_numpy()
                reg = pd.Index(obs).map(meta["pathology_annotation"]).map(cfg["regmap"]).to_numpy()
                X = h["X"][:]                          # dense, 266 genes, normalized linear
            regions = tuple(dict.fromkeys(cfg["regmap"].values()))
            for c in cfg["cts"]:
                # a (sample, cell_type) pair is only usable when BOTH regions clear the gate,
                # so the paired ~ sample + region design stays estimable -- same rule as the
                # imputed pipeline, so the two arms are built on exactly the same groups.
                counts = {r: int(((ct == c) & (reg == r)).sum()) for r in regions}
                if any(counts[r] < MIN_CELLS for r in regions):
                    continue
                for r in regions:
                    m = (ct == c) & (reg == r)
                    cid = f"{sample}__{c}__{r}"
                    cols[cid] = pd.Series(X[m].mean(0), index=genes)
                    design.append({"column_id": cid, "sample": sample, "cell_type": c, "region": r})
            print(f"  {coh} {sid}: done", flush=True)

        pb = pd.DataFrame(cols)
        pb.index.name = "gene"
        pb = pb[~pb.index.duplicated(keep="first")]
        des = pd.DataFrame(design).sort_values(["cell_type", "sample", "region"]).reset_index(drop=True)
        pb.to_csv(OUT / f"measured_{coh}_pseudobulk_matrix.csv")
        des.to_csv(OUT / f"measured_{coh}_sample_design.csv", index=False)
        print(f"{coh}: measured pseudobulk {pb.shape[0]} genes x {pb.shape[1]} columns")

if __name__ == "__main__":
    main()
