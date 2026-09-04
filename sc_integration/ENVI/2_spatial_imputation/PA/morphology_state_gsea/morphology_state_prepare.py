#!/usr/bin/env python3
# Morphology-state pseudobulk builder (Python prepare stage).
#
# Purpose: for the PA cohort, test whether ENVI-imputed cell state differs
# between COMPACT (fibrillary) and LOOSE (myxoid) pathology regions, holding
# cell type constant. This stage builds sample x region pseudobulk means for
# each target cell type on a shared imputed-gene universe; the R stage then
# runs limma (~ sample + region) + fgsea. Pseudobulk-per-sample keeps the
# inferential unit at the sample level (n = usable PA samples), not per-cell.
#
# Inputs
#   - PA ENVI spatial h5ads: */ENVI_results_022626/*astrocytoma_spatial_envi.h5ad
#     imputed full transcriptome lives in obsm['imputation'], genes = var_names.
#   - Phase-0 cell metadata TSV (single source of truth for finalized labels):
#     data/xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv
#     joined on cell_id (== h5ad obs_names == "barcode-1__sampleID").
#
# Outputs (handoff CSVs written next to this script, consumed by morphology_state.R)
#   - *_py_pseudobulk_matrix.csv   gene x (sample__celltype__region) mean imputed expr
#   - *_py_sample_design.csv       column_id, sample, cell_type, region
#   - *_py_shared_universe.csv     shared imputed-gene universe
#   - *_py_coverage_summary.csv    per-sample tested-gene count
#   - *_py_pairwise_gene_overlap.csv  jaccard of imputed genes between samples
#   - *_py_cell_counts.csv         per sample x cell_type x region cell counts + usable flag
#
# Key assumptions
#   - Region is read from pathology_annotation; only the two morphology regions
#     of interest are kept (compact = reference, loose = myxoid).
#   - A sample is "usable" for a cell type only if it has >= MIN_CELLS_PER_GROUP
#     cells in BOTH regions, so the paired limma design stays estimable. Sample
#     258482 carries no compact annotation and is expected to drop everywhere.

import os
import glob
from pathlib import Path
import itertools

import anndata as ad
import numpy as np
import pandas as pd


# Target lineages (raw snake_case keys from the phase-0 object; display names are
# applied only at plot time in the notebook via helpers/cell_type_display.R).
CELL_TYPES = ("opc", "astro_glial", "myeloid")

# Morphology regions; COMPACT is the reference level so the loose coefficient in
# the R stage reads directly as loose-vs-compact.
REGION_MAP = {
    "Compact. fibrillary component": "compact",
    "Loose, myxoid component": "loose",
}
REGIONS = ("compact", "loose")

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


def load_cell_meta(workdir: str) -> pd.DataFrame:
    # One row per cell in the finalized BANKSY object; index by cell_id so we can
    # map h5ad obs_names -> cell_type + region. pathology_annotation -> region.
    meta_path = os.path.join(
        workdir, "youyun/plgg/data/xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv"
    )
    meta = pd.read_csv(
        meta_path, sep="\t", usecols=["cell_id", "cell_type", "pathology_annotation"]
    ).set_index("cell_id")
    meta["region"] = meta["pathology_annotation"].map(REGION_MAP)
    return meta


def load_sample(spatial_h5ad: str, cell_meta: pd.DataFrame) -> dict:
    adata = ad.read_h5ad(spatial_h5ad)
    sample_id = Path(spatial_h5ad).name.replace("_spatial_envi.h5ad", "")
    sample_name = sample_id.replace("_", " ")

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

    # Drop duplicated gene columns (keep first), same as the myeloid prepare stage.
    genes = pd.Index(genes)
    keep_cols = ~genes.duplicated(keep="first")
    genes = genes[keep_cols]
    imputation_mat = imputation_mat[:, keep_cols]

    # Map obs_names (== cell_id) to finalized cell_type and morphology region.
    obs_names = pd.Index(adata.obs_names)
    cell_type_vec = obs_names.map(cell_meta["cell_type"]).to_numpy()
    region_vec = obs_names.map(cell_meta["region"]).to_numpy()

    # Per (cell_type, region) usable-cell counts drive the gating below.
    counts = {}
    for ct in CELL_TYPES:
        for reg in REGIONS:
            counts[(ct, reg)] = int(np.sum((cell_type_vec == ct) & (region_vec == reg)))

    return {
        "sample_id": sample_id,
        "sample": sample_name,
        "genes": genes.to_numpy(),
        "cell_type_vec": cell_type_vec,
        "region_vec": region_vec,
        "imputation_mat": imputation_mat,
        "counts": counts,
    }


def main() -> None:
    workdir = detect_workdir()
    # Keep handoff files local to this pipeline directory.
    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)

    data_root = get_data_root(workdir)
    spatial_files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_spatial_envi.h5ad")))
    if not spatial_files:
        raise FileNotFoundError(f"No PA spatial files found in {data_root}")

    cell_meta = load_cell_meta(workdir)
    sample_objs = [load_sample(f, cell_meta) for f in spatial_files]

    # Long-form per sample x cell_type x region cell counts + usable flag (both
    # regions >= MIN_CELLS_PER_GROUP). A (sample, cell_type) is usable only if
    # both regions clear the threshold, so the within-sample contrast exists.
    count_rows = []
    usable_pairs = set()  # (sample_name, cell_type)
    for obj in sample_objs:
        for ct in CELL_TYPES:
            n_comp = obj["counts"][(ct, "compact")]
            n_loose = obj["counts"][(ct, "loose")]
            usable = (n_comp >= MIN_CELLS_PER_GROUP) and (n_loose >= MIN_CELLS_PER_GROUP)
            if usable:
                usable_pairs.add((obj["sample"], ct))
            count_rows.append({
                "sample": obj["sample"],
                "cell_type": ct,
                "n_compact": n_comp,
                "n_loose": n_loose,
                "usable": usable,
            })
    count_df = pd.DataFrame(count_rows)

    # Samples contributing to at least one cell type define the shared gene
    # universe (imputed genes are per-sample, so intersect across used samples).
    used_sample_names = sorted({s for (s, _ct) in usable_pairs})
    if len(used_sample_names) < 2:
        raise ValueError("Need at least two samples usable in at least one cell type.")
    used_objs = [o for o in sample_objs if o["sample"] in used_sample_names]

    gene_sets = {o["sample"]: set(o["genes"].tolist()) for o in used_objs}
    shared_universe = sorted(set.intersection(*gene_sets.values()))
    if not shared_universe:
        raise ValueError("Shared gene universe is empty.")

    # Pairwise imputed-gene overlap (diagnostic, mirrors the myeloid stage).
    overlap_rows = []
    for a, b in itertools.combinations(sorted(gene_sets.keys()), 2):
        inter = len(gene_sets[a].intersection(gene_sets[b]))
        union = len(gene_sets[a].union(gene_sets[b]))
        overlap_rows.append({
            "sample_a": a,
            "sample_b": b,
            "n_intersection": inter,
            "n_union": union,
            "jaccard": inter / union if union > 0 else np.nan,
        })
    overlap_df = pd.DataFrame(overlap_rows)

    # Build sample x cell_type x region pseudobulk means on the shared universe.
    pb_columns = {}
    design_rows = []
    for obj in used_objs:
        # Map this sample's genes onto the shared universe once.
        idx = pd.Index(obj["genes"]).get_indexer(shared_universe)
        mat = obj["imputation_mat"][:, idx]
        ct_vec = obj["cell_type_vec"]
        reg_vec = obj["region_vec"]
        for ct in CELL_TYPES:
            if (obj["sample"], ct) not in usable_pairs:
                continue  # skip cell types this sample cannot pair across regions
            for reg in REGIONS:
                mask = (ct_vec == ct) & (reg_vec == reg)
                if int(np.sum(mask)) < MIN_CELLS_PER_GROUP:
                    continue
                col_id = f"{obj['sample']}__{ct}__{reg}"
                pb_columns[col_id] = np.nanmean(mat[mask, :], axis=0)
                design_rows.append({
                    "column_id": col_id,
                    "sample": obj["sample"],
                    "cell_type": ct,
                    "region": reg,
                })

    pb_df = pd.DataFrame(pb_columns, index=shared_universe)
    pb_df.index.name = "gene"
    design_df = pd.DataFrame(design_rows)
    design_df = design_df[design_df["column_id"].isin(pb_df.columns)].copy()
    design_df = design_df.sort_values(["cell_type", "sample", "region"]).reset_index(drop=True)

    coverage_df = pd.DataFrame({
        "sample": [o["sample"] for o in sample_objs],
        "n_tested_genes": [len(o["genes"]) for o in sample_objs],
    })
    shared_df = pd.DataFrame({"gene": shared_universe})

    coverage_df.to_csv(output_dir / "morphology_state_py_coverage_summary.csv", index=False)
    overlap_df.to_csv(output_dir / "morphology_state_py_pairwise_gene_overlap.csv", index=False)
    shared_df.to_csv(output_dir / "morphology_state_py_shared_universe.csv", index=False)
    pb_df.to_csv(output_dir / "morphology_state_py_pseudobulk_matrix.csv")
    design_df.to_csv(output_dir / "morphology_state_py_sample_design.csv", index=False)
    count_df.to_csv(output_dir / "morphology_state_py_cell_counts.csv", index=False)

    print("Saved Python handoff files to:", output_dir)
    print("n_samples_total:", len(sample_objs))
    print("n_samples_used:", len(used_sample_names))
    print("shared_universe_n:", len(shared_universe))
    print("pseudobulk_shape:", pb_df.shape)
    for ct in CELL_TYPES:
        n_used = design_df.loc[design_df["cell_type"] == ct, "sample"].nunique()
        print(f"  cell_type={ct}: usable samples (both regions)={n_used}")


if __name__ == "__main__":
    main()
