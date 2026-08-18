#!/usr/bin/env python3

import os
import glob
from pathlib import Path
import itertools

import anndata as ad
import numpy as np
import pandas as pd


TARGET_GROUPS = ("Myeloid 1", "Myeloid 2")
MIN_CELLS_PER_GROUP = 20


def detect_workdir() -> str:
    home = os.path.expanduser("~")
    if home in ["/Users/youyun", "/Users/youyunzheng"]:
        return home + "/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
    if home == "/home/yz762":
        return "/mnt/storage/dept/medonc/beroukhim/"
    if home == "/PHShome/yz762":
        return "/data/beroukhim1/"
    return "/data/beroukhim1/"


def get_data_root(workdir: str) -> str:
    candidates = [
        os.path.join(workdir, "youyun/plgg/data/sc_integration/ENVI_results_022626"),
        os.path.join(workdir, "youyun/plgg/data/sc_integration/ENVI_results"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"No ENVI result directory found among: {candidates}")


def load_sample(spatial_h5ad: str) -> dict:
    adata = ad.read_h5ad(spatial_h5ad)
    sample_id = Path(spatial_h5ad).name.replace("_spatial_envi.h5ad", "")
    sample_name = sample_id.replace("_", " ")

    if "Fine.Cell.Type.UMAP" not in adata.obs.columns:
        raise KeyError(f"Missing Fine.Cell.Type.UMAP in obs: {sample_name}")

    if "imputation" not in adata.obsm_keys():
        raise KeyError(f"Missing obsm['imputation']: {sample_name}")

    imp = adata.obsm["imputation"]
    if isinstance(imp, pd.DataFrame):
        imputation_mat = imp.to_numpy()
        genes = imp.columns.astype(str).str.upper().tolist()
    else:
        imputation_mat = np.asarray(imp)
        if adata.n_vars == imputation_mat.shape[1]:
            genes = pd.Index(adata.var_names).astype(str).str.upper().tolist()
        else:
            genes = [f"GENE_{i+1}" for i in range(imputation_mat.shape[1])]

    genes = pd.Index(genes)
    keep_cols = ~genes.duplicated(keep="first")
    genes = genes[keep_cols]
    imputation_mat = imputation_mat[:, keep_cols]

    group_vec = adata.obs["Fine.Cell.Type.UMAP"].astype("string")
    myeloid_mask = group_vec.isin(TARGET_GROUPS).fillna(False).to_numpy()
    myeloid_groups = group_vec[myeloid_mask].astype(str).to_numpy()
    myeloid_mat = imputation_mat[myeloid_mask, :]

    n_m1 = int(np.sum(myeloid_groups == "Myeloid 1"))
    n_m2 = int(np.sum(myeloid_groups == "Myeloid 2"))

    return {
        "sample_id": sample_id,
        "sample": sample_name,
        "genes": genes.to_numpy(),
        "myeloid_groups": myeloid_groups,
        "myeloid_mat": myeloid_mat,
        "n_cells_myeloid1": n_m1,
        "n_cells_myeloid2": n_m2,
        "used_in_spatial_imputation": n_m1 >= MIN_CELLS_PER_GROUP and n_m2 >= MIN_CELLS_PER_GROUP,
    }


def main() -> None:
    workdir = detect_workdir()
    # Keep handoff files local to the spatial-imputation pipeline directory.
    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)

    data_root = get_data_root(workdir)
    spatial_files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_spatial_envi.h5ad")))
    if not spatial_files:
        raise FileNotFoundError(f"No PA spatial files found in {data_root}")

    sample_objs = [load_sample(f) for f in spatial_files]
    diag_df = pd.DataFrame({
        "sample": [x["sample"] for x in sample_objs],
        "n_tested_genes": [len(x["genes"]) for x in sample_objs],
        "n_cells_myeloid1": [x["n_cells_myeloid1"] for x in sample_objs],
        "n_cells_myeloid2": [x["n_cells_myeloid2"] for x in sample_objs],
        "used_in_spatial_imputation": [x["used_in_spatial_imputation"] for x in sample_objs],
    })

    # Require minimum cells in both groups so downstream limma design is stable.
    used_objs = [x for x in sample_objs if x["used_in_spatial_imputation"]]
    if len(used_objs) < 2:
        raise ValueError("Need at least two samples with >=20 cells in each myeloid group.")

    gene_sets = {x["sample"]: set(x["genes"].tolist()) for x in used_objs}
    shared_universe = sorted(set.intersection(*gene_sets.values()))
    if not shared_universe:
        raise ValueError("Shared gene universe is empty.")

    overlap_rows = []
    for a, b in itertools.combinations(sorted(gene_sets.keys()), 2):
        a_set = gene_sets[a]
        b_set = gene_sets[b]
        inter = len(a_set.intersection(b_set))
        union = len(a_set.union(b_set))
        overlap_rows.append({
            "sample_a": a,
            "sample_b": b,
            "n_intersection": inter,
            "n_union": union,
            "jaccard": inter / union if union > 0 else np.nan,
        })
    overlap_df = pd.DataFrame(overlap_rows)

    # Build sample-level pseudobulk means on a shared tested-gene universe.
    pb_columns = {}
    design_rows = []
    for obj in used_objs:
        idx = pd.Index(obj["genes"]).get_indexer(shared_universe)
        sample_mat = obj["myeloid_mat"][:, idx]
        sample_grp = obj["myeloid_groups"]

        for group in TARGET_GROUPS:
            col_id = f"{obj['sample']}__{group}"
            if np.sum(sample_grp == group) == 0:
                continue
            pb_columns[col_id] = np.nanmean(sample_mat[sample_grp == group, :], axis=0)
            design_rows.append({
                "column_id": col_id,
                "sample": obj["sample"],
                "group": group,
            })

    pb_df = pd.DataFrame(pb_columns, index=shared_universe)
    pb_df.index.name = "gene"
    design_df = pd.DataFrame(design_rows)
    design_df = design_df[design_df["column_id"].isin(pb_df.columns)].copy()
    design_df = design_df.sort_values(["sample", "group"]).reset_index(drop=True)

    coverage_df = diag_df[["sample", "n_tested_genes"]].copy()
    shared_df = pd.DataFrame({"gene": shared_universe})

    coverage_df.to_csv(output_dir / "myeloid_spatial_imputation_py_coverage_summary.csv", index=False)
    overlap_df.to_csv(output_dir / "myeloid_spatial_imputation_py_pairwise_gene_overlap.csv", index=False)
    shared_df.to_csv(output_dir / "myeloid_spatial_imputation_py_shared_universe.csv", index=False)
    pb_df.to_csv(output_dir / "myeloid_spatial_imputation_py_pseudobulk_matrix.csv")
    design_df.to_csv(output_dir / "myeloid_spatial_imputation_py_sample_design.csv", index=False)
    diag_df.to_csv(output_dir / "myeloid_spatial_imputation_py_sample_diagnostics.csv", index=False)

    print("Saved Python handoff files to:", output_dir)
    print("n_samples_total:", diag_df.shape[0])
    print("n_samples_used:", int(diag_df["used_in_spatial_imputation"].sum()))
    print("shared_universe_n:", len(shared_universe))
    print("pseudobulk_shape:", pb_df.shape)


if __name__ == "__main__":
    main()
