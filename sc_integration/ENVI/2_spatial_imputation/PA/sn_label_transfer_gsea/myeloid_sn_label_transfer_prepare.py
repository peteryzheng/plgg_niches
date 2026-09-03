#!/usr/bin/env python3

import glob
import itertools
import os
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

TARGET_GROUPS = ("Myeloid 1", "Myeloid 2")
SN_MYELOID_SOURCE_LABELS = ("Monocytes", "Microglia", "Macrophage")
SN_MYELOID_SOURCE_LABELS_NORM = {x.lower() for x in SN_MYELOID_SOURCE_LABELS}
MIN_CELLS_PER_GROUP = 20
FALLBACK_MIN_CELLS_PER_GROUP = 5
K_NEIGHBORS = 25
CHUNK_SIZE = 256
CONFIDENCE_LOW_THRESHOLD = 0.6


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


def to_dense(arr):
    if sparse.issparse(arr):
        return arr.toarray()
    return np.asarray(arr)


def normalize_log1p_rows(x: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    row_sums = x.sum(axis=1)
    scale = np.zeros_like(row_sums, dtype=float)
    nonzero = row_sums > 0
    scale[nonzero] = target_sum / row_sums[nonzero]
    x = x * scale[:, None]
    return np.log1p(x)


def make_coembedding_plot(
    sample_name: str,
    sn_latent_all: np.ndarray,
    sn_original_labels: np.ndarray,
    sn_candidate_mask: np.ndarray,
    sn_transferred_labels: np.ndarray,
    sp_latent_all: np.ndarray,
    sp_original_labels: np.ndarray,
    out_path: Path,
) -> None:
    myeloid_colors = {
        "Myeloid 1": "#1f77b4",
        "Myeloid 2": "#d62728",
    }
    sn_source_colors = {
        "Monocytes": "#2ca02c",
        "Microglia": "#ff7f0e",
        "Macrophage": "#9467bd",
    }
    other_color = "#d0d0d0"

    lim_arr = np.concatenate([sn_latent_all, sp_latent_all], axis=0)
    delta = 1.0
    pre = 0.1
    xmin = np.percentile(lim_arr[:, 0], pre) - delta
    xmax = np.percentile(lim_arr[:, 0], 100 - pre) + delta
    ymin = np.percentile(lim_arr[:, 1], pre) - delta
    ymax = np.percentile(lim_arr[:, 1], 100 - pre) + delta

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    sn_original_series = pd.Series(sn_original_labels, dtype="string")
    sn_original_norm = sn_original_series.str.strip().str.lower().fillna("unlabeled")
    sn_transfer_full = np.full(sn_latent_all.shape[0], "Other", dtype=object)
    sn_transfer_full[sn_candidate_mask] = sn_transferred_labels
    sp_series = pd.Series(sp_original_labels, dtype="string")
    sp_norm = sp_series.str.strip()

    # Panel 1: SN latent with original source labels.
    sn_other_mask = ~sn_original_norm.isin(SN_MYELOID_SOURCE_LABELS_NORM).fillna(False).to_numpy(dtype=bool)
    if np.any(sn_other_mask):
        axes[0].scatter(
            sn_latent_all[sn_other_mask, 0],
            sn_latent_all[sn_other_mask, 1],
            s=5,
            alpha=0.35,
            c=other_color,
            label="Other",
            linewidths=0,
        )
    for label in SN_MYELOID_SOURCE_LABELS:
        mask = (sn_original_norm == label.lower()).fillna(False).to_numpy(dtype=bool)
        if np.any(mask):
            axes[0].scatter(
                sn_latent_all[mask, 0],
                sn_latent_all[mask, 1],
                s=6,
                alpha=0.8,
                c=sn_source_colors[label],
                label=label,
                linewidths=0,
            )
    axes[0].set_title("SN Latent (Original labels)")
    axes[0].set_xlim([xmin, xmax])
    axes[0].set_ylim([ymin, ymax])
    axes[0].axis("off")
    axes[0].legend(loc="upper right", fontsize=8)

    # Panel 2: SN latent with transferred Myeloid 1/2 labels.
    transfer_other_mask = sn_transfer_full == "Other"
    if np.any(transfer_other_mask):
        axes[1].scatter(
            sn_latent_all[transfer_other_mask, 0],
            sn_latent_all[transfer_other_mask, 1],
            s=5,
            alpha=0.35,
            c=other_color,
            label="Other",
            linewidths=0,
        )
    for group in TARGET_GROUPS:
        mask = sn_transfer_full == group
        if np.any(mask):
            axes[1].scatter(
                sn_latent_all[mask, 0],
                sn_latent_all[mask, 1],
                s=6,
                alpha=0.8,
                c=myeloid_colors[group],
                label=group,
                linewidths=0,
            )
    axes[1].set_title("SN Latent (Transferred labels)")
    axes[1].set_xlim([xmin, xmax])
    axes[1].set_ylim([ymin, ymax])
    axes[1].axis("off")
    axes[1].legend(loc="upper right", fontsize=8)

    # Panel 3: Spatial latent with original Myeloid 1/2 labels.
    sp_target_mask = sp_norm.isin(TARGET_GROUPS).fillna(False).to_numpy(dtype=bool)
    sp_other_mask = ~sp_target_mask
    if np.any(sp_other_mask):
        axes[2].scatter(
            sp_latent_all[sp_other_mask, 0],
            sp_latent_all[sp_other_mask, 1],
            s=5,
            alpha=0.35,
            c=other_color,
            label="Other",
            linewidths=0,
        )
    for group in TARGET_GROUPS:
        mask = (sp_norm == group).fillna(False).to_numpy(dtype=bool)
        if np.any(mask):
            axes[2].scatter(
                sp_latent_all[mask, 0],
                sp_latent_all[mask, 1],
                s=6,
                alpha=0.8,
                c=myeloid_colors[group],
                label=group,
                linewidths=0,
            )
    axes[2].set_title("Spatial Latent (Original labels)")
    axes[2].set_xlim([xmin, xmax])
    axes[2].set_ylim([ymin, ymax])
    axes[2].axis("off")
    axes[2].legend(loc="upper right", fontsize=8)

    n_candidates = int(sn_candidate_mask.sum())
    n_m1 = int(np.sum(sn_transferred_labels == "Myeloid 1"))
    n_m2 = int(np.sum(sn_transferred_labels == "Myeloid 2"))
    fig.suptitle(
        f"{sample_name}\nSN candidates={n_candidates}, transferred Myeloid 1={n_m1}, Myeloid 2={n_m2}"
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def transfer_labels(sn_latent: np.ndarray, sp_latent: np.ndarray, sp_labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_train = sp_latent.shape[0]
    if n_train == 0:
        raise ValueError("No spatial training cells available for label transfer.")

    k_eff = min(K_NEIGHBORS, n_train)
    tree = cKDTree(sp_latent)
    dists, neigh_idx = tree.query(sn_latent, k=k_eff)

    if k_eff == 1:
        dists = dists[:, None]
        neigh_idx = neigh_idx[:, None]

    neigh_labels = sp_labels[neigh_idx]
    weights = 1.0 / (dists + 1e-6)
    weight_sum = weights.sum(axis=1)
    weight_sum[weight_sum == 0] = 1.0

    w_m1 = (weights * (neigh_labels == "Myeloid 1")).sum(axis=1)
    score_m1 = w_m1 / weight_sum

    pred = np.where(score_m1 >= 0.5, "Myeloid 1", "Myeloid 2")
    conf = np.maximum(score_m1, 1.0 - score_m1)
    return pred.astype("U16"), conf.astype(float)


def summarize_sample(sample_id: str, sn_path: str, spatial_path: str, output_dir: Path, coembed_dir: Path) -> dict:
    sample_name = sample_id.replace("_", " ")

    # Use backed mode to avoid loading full h5ads into memory.
    sn_adata = ad.read_h5ad(sn_path, backed="r")
    sp_adata = ad.read_h5ad(spatial_path, backed="r")

    try:
        if "Fine.Cell.Type.UMAP" not in sp_adata.obs.columns:
            raise KeyError("missing spatial obs['Fine.Cell.Type.UMAP']")
        if "latent_umap" not in sp_adata.obsm.keys():
            raise KeyError("missing spatial obsm['latent_umap']")
        if "annot_v1" not in sn_adata.obs.columns:
            raise KeyError("missing SN obs['annot_v1']")
        if "latent_umap" not in sn_adata.obsm.keys():
            raise KeyError("missing SN obsm['latent_umap']")

        sp_group = sp_adata.obs["Fine.Cell.Type.UMAP"].astype("string")
        sp_latent_all = np.asarray(sp_adata.obsm["latent_umap"])[:, :2]
        sp_mask = sp_group.isin(TARGET_GROUPS).fillna(False).to_numpy(dtype=bool)
        sp_labels = sp_group[sp_mask].astype(str).to_numpy()
        sp_latent = sp_latent_all[sp_mask, :]

        if sp_latent.shape[0] == 0:
            raise ValueError("no spatial Myeloid 1/2 cells")

        sn_annot = sn_adata.obs["annot_v1"].astype("string")
        sn_annot_norm = sn_annot.str.strip().str.lower()
        sn_mask = sn_annot_norm.isin(SN_MYELOID_SOURCE_LABELS_NORM).fillna(False).to_numpy(dtype=bool)
        sn_idx = np.where(sn_mask)[0]
        if sn_idx.size == 0:
            raise ValueError(
                "no SN candidates from annot_v1 labels "
                f"{sorted(SN_MYELOID_SOURCE_LABELS)}"
            )

        sn_latent_all = np.asarray(sn_adata.obsm["latent_umap"])[:, :2]
        sn_latent = sn_latent_all[sn_idx, :]

        # Transfer Myeloid 1/2 labels in latent space via weighted kNN.
        pred_labels, conf = transfer_labels(sn_latent=sn_latent, sp_latent=sp_latent, sp_labels=sp_labels)

        plot_path = coembed_dir / f"{sample_id}_coembedding_post_transfer.png"
        make_coembedding_plot(
            sample_name=sample_name,
            sn_latent_all=sn_latent_all,
            sn_original_labels=sn_annot.to_numpy(),
            sn_candidate_mask=sn_mask,
            sn_transferred_labels=pred_labels,
            sp_latent_all=sp_latent_all,
            sp_original_labels=sp_group.to_numpy(),
            out_path=plot_path,
        )

        genes = pd.Index(sn_adata.var_names).astype(str).str.upper()
        keep_cols = ~genes.duplicated(keep="first")
        keep_col_idx = np.where(keep_cols)[0]
        unique_genes = genes[keep_cols].to_numpy()

        sums = {
            "Myeloid 1": np.zeros(len(unique_genes), dtype=float),
            "Myeloid 2": np.zeros(len(unique_genes), dtype=float),
        }
        counts = {
            "Myeloid 1": int(np.sum(pred_labels == "Myeloid 1")),
            "Myeloid 2": int(np.sum(pred_labels == "Myeloid 2")),
        }

        # Aggregate normalized observed SN expression in chunks for memory safety.
        for start in range(0, len(sn_idx), CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, len(sn_idx))
            idx_chunk = sn_idx[start:end]
            label_chunk = pred_labels[start:end]

            x_chunk = to_dense(sn_adata.X[idx_chunk, :])
            x_chunk = np.asarray(x_chunk, dtype=float)
            x_chunk = x_chunk[:, keep_col_idx]
            x_chunk = normalize_log1p_rows(x_chunk)

            for group in TARGET_GROUPS:
                gmask = label_chunk == group
                if np.any(gmask):
                    sums[group] += x_chunk[gmask, :].sum(axis=0)

        used_in_analysis = counts["Myeloid 1"] >= MIN_CELLS_PER_GROUP and counts["Myeloid 2"] >= MIN_CELLS_PER_GROUP

        mean_expr = {}
        for group in TARGET_GROUPS:
            if counts[group] > 0:
                mean_expr[group] = sums[group] / counts[group]
            else:
                mean_expr[group] = np.zeros(len(unique_genes), dtype=float)

        pb_cols = {
            f"{sample_name}__Myeloid 1": mean_expr["Myeloid 1"],
            f"{sample_name}__Myeloid 2": mean_expr["Myeloid 2"],
        }

        design_rows = [
            {"column_id": f"{sample_name}__Myeloid 1", "sample": sample_name, "group": "Myeloid 1"},
            {"column_id": f"{sample_name}__Myeloid 2", "sample": sample_name, "group": "Myeloid 2"},
        ]

        result = {
            "sample_id": sample_id,
            "sample": sample_name,
            "genes": unique_genes,
            "pb_cols": pb_cols,
            "design_rows": design_rows,
            "n_tested_genes": len(unique_genes),
            "n_spatial_train_m1": int(np.sum(sp_labels == "Myeloid 1")),
            "n_spatial_train_m2": int(np.sum(sp_labels == "Myeloid 2")),
            "n_sn_myeloid_candidates": int(len(sn_idx)),
            "n_pred_m1": counts["Myeloid 1"],
            "n_pred_m2": counts["Myeloid 2"],
            "mean_confidence": float(np.mean(conf)),
            "median_confidence": float(np.median(conf)),
            "frac_conf_lt_0_6": float(np.mean(conf < CONFIDENCE_LOW_THRESHOLD)),
            "sn_subset_rule": "annot_v1 exact in {Monocytes, Microglia, Macrophage}",
            "used_in_sn_label_transfer": used_in_analysis,
            "coembedding_plot_path": str(plot_path),
            "error": "",
        }
        return result
    finally:
        if sn_adata.isbacked:
            sn_adata.file.close()
        if sp_adata.isbacked:
            sp_adata.file.close()


def main() -> None:
    workdir = detect_workdir()
    # Write handoff files alongside this script to keep pipeline-local I/O.
    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)
    coembed_dir = output_dir / "coembedding"
    coembed_dir.mkdir(parents=True, exist_ok=True)

    data_root = get_data_root(workdir)
    sn_files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_sn_envi.h5ad")))
    sp_files = sorted(glob.glob(os.path.join(data_root, "*astrocytoma_spatial_envi.h5ad")))

    sn_by_sample = {
        re.sub("_sn_envi.h5ad$", "", Path(x).name): x
        for x in sn_files
    }
    sp_by_sample = {
        re.sub("_spatial_envi.h5ad$", "", Path(x).name): x
        for x in sp_files
    }

    sample_ids = sorted(set(sn_by_sample.keys()) & set(sp_by_sample.keys()))
    if not sample_ids:
        raise FileNotFoundError("No matched SN/spatial ENVI samples found.")

    sample_results = []
    for sample_id in sample_ids:
        try:
            sample_results.append(
                summarize_sample(
                    sample_id=sample_id,
                    sn_path=sn_by_sample[sample_id],
                    spatial_path=sp_by_sample[sample_id],
                    output_dir=output_dir,
                    coembed_dir=coembed_dir,
                )
            )
        except Exception as err:
            sample_name = sample_id.replace("_", " ")
            sample_results.append({
                "sample_id": sample_id,
                "sample": sample_name,
                "genes": np.array([], dtype=str),
                "pb_cols": {},
                "design_rows": [],
                "n_tested_genes": 0,
                "n_spatial_train_m1": 0,
                "n_spatial_train_m2": 0,
                "n_sn_myeloid_candidates": 0,
                "n_pred_m1": 0,
                "n_pred_m2": 0,
                "mean_confidence": np.nan,
                "median_confidence": np.nan,
                "frac_conf_lt_0_6": np.nan,
                "sn_subset_rule": "annot_v1 exact in {Monocytes, Microglia, Macrophage}",
                "used_in_sn_label_transfer": False,
                "coembedding_plot_path": "",
                "error": str(err),
            })

    diag_df = pd.DataFrame([
        {
            "sample": x["sample"],
            "n_tested_genes": x["n_tested_genes"],
            "n_spatial_train_m1": x["n_spatial_train_m1"],
            "n_spatial_train_m2": x["n_spatial_train_m2"],
            "n_sn_myeloid_candidates": x["n_sn_myeloid_candidates"],
            "n_pred_m1": x["n_pred_m1"],
            "n_pred_m2": x["n_pred_m2"],
            "used_in_sn_label_transfer": x["used_in_sn_label_transfer"],
            "error": x["error"],
        }
        for x in sample_results
    ])

    transfer_diag_df = pd.DataFrame([
        {
            "sample": x["sample"],
            "n_spatial_train_m1": x["n_spatial_train_m1"],
            "n_spatial_train_m2": x["n_spatial_train_m2"],
            "n_sn_myeloid_candidates": x["n_sn_myeloid_candidates"],
            "n_pred_m1": x["n_pred_m1"],
            "n_pred_m2": x["n_pred_m2"],
            "mean_confidence": x["mean_confidence"],
            "median_confidence": x["median_confidence"],
            "frac_conf_lt_0_6": x["frac_conf_lt_0_6"],
            "sn_subset_rule": x["sn_subset_rule"],
            "coembedding_plot_path": x["coembedding_plot_path"],
            "error": x["error"],
        }
        for x in sample_results
    ])
    print(transfer_diag_df[["sample", "n_spatial_train_m1", "n_spatial_train_m2", "n_sn_myeloid_candidates", "n_pred_m1", "n_pred_m2", "mean_confidence", "median_confidence", "frac_conf_lt_0_6", "error"]])
    print(diag_df[["sample", "n_tested_genes", "n_spatial_train_m1", "n_spatial_train_m2", "n_sn_myeloid_candidates", "n_pred_m1", "n_pred_m2", "used_in_sn_label_transfer", "error"]])

    # Keep only samples with enough transferred cells in both groups for limma design.
    primary_used = [x for x in sample_results if x["used_in_sn_label_transfer"]]
    fallback_used = [
        x for x in sample_results
        if x["n_pred_m1"] >= FALLBACK_MIN_CELLS_PER_GROUP and x["n_pred_m2"] >= FALLBACK_MIN_CELLS_PER_GROUP
    ]
    minimal_used = [
        x for x in sample_results
        if x["n_pred_m1"] >= 1 and x["n_pred_m2"] >= 1
    ]

    if len(primary_used) >= 2:
        used_objs = primary_used
        selection_rule = f">={MIN_CELLS_PER_GROUP} cells per transferred group"
        analysis_mode = "cohort"
    elif len(primary_used) == 1:
        used_objs = primary_used
        selection_rule = f">={MIN_CELLS_PER_GROUP} cells per transferred group (single-sample exploratory)"
        analysis_mode = "single_sample_exploratory"
    elif len(fallback_used) >= 2:
        used_objs = fallback_used
        selection_rule = f">={FALLBACK_MIN_CELLS_PER_GROUP} cells per transferred group (fallback)"
        analysis_mode = "cohort"
    elif len(fallback_used) == 1:
        used_objs = fallback_used
        selection_rule = f">={FALLBACK_MIN_CELLS_PER_GROUP} cells per transferred group (single-sample exploratory fallback)"
        analysis_mode = "single_sample_exploratory"
    elif len(minimal_used) >= 2:
        used_objs = minimal_used
        selection_rule = ">=1 cell per transferred group (minimal fallback)"
        analysis_mode = "cohort"
    elif len(minimal_used) == 1:
        used_objs = minimal_used
        selection_rule = ">=1 cell per transferred group (single-sample exploratory minimal fallback)"
        analysis_mode = "single_sample_exploratory"
    else:
        raise ValueError(
            "Need at least one sample with both transferred groups present. "
            "No threshold produced an eligible sample."
        )

    selected_samples = {x["sample"] for x in used_objs}
    diag_df["selected_for_analysis"] = diag_df["sample"].isin(selected_samples)
    diag_df["analysis_mode"] = analysis_mode
    transfer_diag_df["selected_for_analysis"] = transfer_diag_df["sample"].isin(selected_samples)
    transfer_diag_df["selection_rule"] = selection_rule
    transfer_diag_df["analysis_mode"] = analysis_mode

    gene_sets = {x["sample"]: set(x["genes"].tolist()) for x in used_objs}
    shared_universe = sorted(set.intersection(*gene_sets.values()))
    if not shared_universe:
        raise ValueError("Shared gene universe is empty for SN label-transfer analysis.")

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
    overlap_df = pd.DataFrame(
        overlap_rows,
        columns=["sample_a", "sample_b", "n_intersection", "n_union", "jaccard"]
    )

    pb_columns = {}
    design_rows = []
    for obj in used_objs:
        idx = pd.Index(obj["genes"]).get_indexer(shared_universe)
        for row in obj["design_rows"]:
            col = row["column_id"]
            group = row["group"]
            pb_columns[col] = obj["pb_cols"][col][idx]
            design_rows.append(row)

    pb_df = pd.DataFrame(pb_columns, index=shared_universe)
    pb_df.index.name = "gene"
    design_df = pd.DataFrame(design_rows)
    design_df = design_df[design_df["column_id"].isin(pb_df.columns)].copy()
    design_df = design_df.sort_values(["sample", "group"]).reset_index(drop=True)

    coverage_df = diag_df[["sample", "n_tested_genes"]].copy()
    shared_df = pd.DataFrame({"gene": shared_universe})

    coverage_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_coverage_summary.csv", index=False)
    overlap_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_pairwise_gene_overlap.csv", index=False)
    shared_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_shared_universe.csv", index=False)
    pb_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_pseudobulk_matrix.csv")
    design_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_sample_design.csv", index=False)
    diag_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_sample_diagnostics.csv", index=False)
    transfer_diag_df.to_csv(output_dir / "myeloid_sn_label_transfer_py_transfer_diagnostics.csv", index=False)

    print("Saved SN label-transfer handoff files to:", output_dir)
    print("n_samples_total:", diag_df.shape[0])
    print("n_samples_used:", len(used_objs))
    print("analysis_mode:", analysis_mode)
    print("selection_rule:", selection_rule)
    print("shared_universe_n:", len(shared_universe))
    print("pseudobulk_shape:", pb_df.shape)


if __name__ == "__main__":
    main()
