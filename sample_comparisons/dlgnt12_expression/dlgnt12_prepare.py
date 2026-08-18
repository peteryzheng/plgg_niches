#!/usr/bin/env python3
"""Prepare DLGNT_1 vs DLGNT_2 pseudobulk handoffs from a backed H5AD.

This script reads single-cell metadata and raw counts from an AnnData HDF5 file,
summarizes cell-type support across annot_v0/annot_v1/annot_v2, and writes one
raw-count pseudobulk handoff per retained cell type.

Inputs
------
- H5AD at the repo-standard workdir-relative single-cell path
- cell-type column: annot_v0, annot_v1, or annot_v2

Outputs
-------
- support summaries at all annotation levels
- one raw-count pseudobulk matrix per retained cell type
- aligned design and sample diagnostics tables

Assumptions
-----------
- The analysis is exploratory because DLGNT subtype and location are confounded.
- The input H5AD is always stored relative to the standard runtime-dependent workdir.
- layers["counts"] is stored as a CSC sparse matrix with cells as rows.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

VALID_CELL_TYPE_COLUMNS = ("annot_v0", "annot_v1", "annot_v2")
DLGNT_LEVELS = ("DLGNT_1", "DLGNT_2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare DLGNT_1 vs DLGNT_2 raw pseudobulk handoffs."
    )
    parser.add_argument(
        "--cell-type-column",
        default="annot_v1",
        choices=VALID_CELL_TYPE_COLUMNS,
        help="Cell-type definition used for the per-cell-type pseudobulk analysis.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Analysis directory. Defaults to the directory containing this script.",
    )
    parser.add_argument(
        "--min-samples-per-subtype",
        type=int,
        default=3,
        help="Minimum retained samples in each subtype for a cell type to be analyzed.",
    )
    parser.add_argument(
        "--min-cells-per-sample",
        type=int,
        default=20,
        help="Minimum cells in a sample x cell_type pseudobulk column.",
    )
    return parser.parse_args()


def detect_workdir() -> str:
    home = os.path.expanduser("~")
    # Match the repo-wide mounted-data convention so the same script works on
    # local machines and cluster environments without path edits.
    if home in ["/Users/youyun", "/Users/youyunzheng"]:
        return home + "/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
    if home == "/PHShome/yz762":
        return "/data/beroukhim1/"
    if home == "/home/yz762":
        return "/mnt/storage/dept/medonc/beroukhim/"
    return "/xchip/beroukhimlab/"


def default_h5ad_path(workdir: str) -> Path:
    return Path(workdir) / "youyun/plgg/data/single_cell/Extended_sndata_filtered.csv_subset.h5ad"


def decode_scalar(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_obs_column(obs_group: h5py.Group, key: str) -> np.ndarray:
    obj = obs_group[key]
    if isinstance(obj, h5py.Group):
        # AnnData categoricals are stored as integer codes plus category labels.
        # Decode them here so downstream filtering works on plain strings.
        categories = np.array([decode_scalar(x) for x in obj["categories"][:]], dtype=object)
        codes = obj["codes"][:]
        out = np.full(codes.shape[0], "", dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    return np.array([decode_scalar(x) for x in obj[:]], dtype=object)


def make_safe_names(values: list[str]) -> dict[str, str]:
    used: set[str] = set()
    out: dict[str, str] = {}
    for value in sorted(values):
        # File stems need to stay stable across reruns even when cell-type labels
        # contain spaces or punctuation.
        base = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
        if not base:
            base = "cell_type"
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate)
        out[value] = candidate
    return out


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def build_support_tables(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    detail_frames: list[pd.DataFrame] = []

    for level in VALID_CELL_TYPE_COLUMNS:
        # Count distinct samples per subtype rather than cells because support is
        # defined at the replicate level for the downstream pseudobulk analysis.
        level_df = metadata.loc[metadata[level] != "", ["orig.ident", "IntDx", level]].drop_duplicates()
        counts = (
            level_df.groupby([level, "IntDx"])["orig.ident"]
            .nunique()
            .unstack(fill_value=0)
            .reindex(columns=list(DLGNT_LEVELS), fill_value=0)
            .reset_index()
            .rename(
                columns={
                    level: "cell_type",
                    "DLGNT_1": "n_samples_dlgnt1",
                    "DLGNT_2": "n_samples_dlgnt2",
                }
            )
            .sort_values("cell_type")
            .reset_index(drop=True)
        )
        counts["annotation_level"] = level
        detail_frames.append(counts)

        for threshold in range(1, 5):
            # These threshold summaries mirror the planning table used to decide
            # which annotation level has enough subtype support.
            summary_rows.append(
                {
                    "annotation_level": level,
                    "sample_threshold_per_subtype": threshold,
                    "n_cell_types": int(
                        (
                            (counts["n_samples_dlgnt1"] >= threshold)
                            & (counts["n_samples_dlgnt2"] >= threshold)
                        ).sum()
                    ),
                }
            )

    detail_df = pd.concat(detail_frames, ignore_index=True)
    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["annotation_level", "sample_threshold_per_subtype"]
    )
    return summary_df, detail_df


def aggregate_pseudobulk(
    counts_group: h5py.Group,
    row_to_group: np.ndarray,
    n_genes: int,
    n_groups: int,
) -> np.ndarray:
    indptr = counts_group["indptr"][:]
    indices_ds = counts_group["indices"]
    data_ds = counts_group["data"]
    agg = np.zeros((n_genes, n_groups), dtype=np.int64)

    # Stream one CSC column (gene) at a time so the full sparse matrix never
    # needs to be materialized in memory.
    for gene_idx in range(n_genes):
        if gene_idx % 1000 == 0:
            print(f"Aggregating gene {gene_idx:,} / {n_genes:,}")
        start = int(indptr[gene_idx])
        end = int(indptr[gene_idx + 1])
        if start == end:
            continue

        rows = indices_ds[start:end]
        groups = row_to_group[rows]
        keep = groups >= 0
        if not np.any(keep):
            continue

        values = np.rint(np.asarray(data_ds[start:end])[keep]).astype(np.int64, copy=False)
        agg[gene_idx, :] = np.bincount(groups[keep], weights=values, minlength=n_groups).astype(
            np.int64
        )
    return agg


def main() -> None:
    args = parse_args()

    script_dir = Path(__file__).resolve().parent
    analysis_dir = Path(args.outdir).resolve() if args.outdir else script_dir
    handoff_dir = analysis_dir / "handoff"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    handoff_dir.mkdir(parents=True, exist_ok=True)
    workdir = detect_workdir()
    h5ad_path = default_h5ad_path(workdir)

    if not h5ad_path.exists():
        raise FileNotFoundError(
            f"Expected H5AD at {h5ad_path} derived from workdir {workdir}, but it was not found."
        )

    with h5py.File(h5ad_path, "r") as handle:
        obs_group = handle["obs"]
        var_names = np.array([decode_scalar(x) for x in handle["var"]["_index"][:]], dtype=object)

        orig_ident = read_obs_column(obs_group, "orig.ident")
        intdx = read_obs_column(obs_group, "IntDx")
        sex = read_obs_column(obs_group, "Sex")
        location_key = "Location_standard" if "Location_standard" in obs_group else "Location"
        location = read_obs_column(obs_group, location_key)
        annot_cols = {col: read_obs_column(obs_group, col) for col in VALID_CELL_TYPE_COLUMNS}

        # Restrict immediately to the two DLGNT subtypes of interest so all later
        # summaries and pseudobulks operate on the same cohort.
        dlgnt_mask = np.isin(intdx, DLGNT_LEVELS)
        cell_index = np.flatnonzero(dlgnt_mask)
        if cell_index.size == 0:
            raise ValueError("No DLGNT_1 / DLGNT_2 cells were found in the input H5AD.")

        metadata = pd.DataFrame(
            {
                "cell_index": cell_index,
                "orig.ident": orig_ident[dlgnt_mask],
                "IntDx": intdx[dlgnt_mask],
                "Sex": sex[dlgnt_mask],
                "Location": location[dlgnt_mask],
                "annot_v0": annot_cols["annot_v0"][dlgnt_mask],
                "annot_v1": annot_cols["annot_v1"][dlgnt_mask],
                "annot_v2": annot_cols["annot_v2"][dlgnt_mask],
            }
        )

        # Normalize obvious missing-value encodings once so support counts and
        # output tables do not split on whitespace or literal "NA" strings.
        for col in ["orig.ident", "IntDx", "Sex", "Location", *VALID_CELL_TYPE_COLUMNS]:
            metadata[col] = metadata[col].astype(str).str.strip()
        metadata.loc[metadata["Location"] == "NA", "Location"] = ""
        metadata = metadata.loc[metadata["IntDx"].isin(DLGNT_LEVELS)].reset_index(drop=True)

        sample_metadata = (
            metadata.loc[:, ["orig.ident", "IntDx", "Sex", "Location"]]
            .drop_duplicates()
            .sort_values(["IntDx", "orig.ident"])
            .reset_index(drop=True)
        )
        if sample_metadata["orig.ident"].duplicated().any():
            raise ValueError("Each sample must map to one IntDx / Sex / Location combination.")

        # Save the all-level support tables before choosing one annotation level
        # for pseudobulk so the notebook can report the tradeoff transparently.
        support_summary, support_detail = build_support_tables(metadata)
        write_tsv(sample_metadata, analysis_dir / "dlgnt_sample_metadata.tsv")
        write_tsv(support_summary, analysis_dir / "annotation_level_support_summary.tsv")
        write_tsv(support_detail, analysis_dir / "annotation_level_support_detail.tsv")

        analysis_meta = metadata.loc[:, ["cell_index", "orig.ident", "IntDx", "Sex", "Location"]].copy()
        analysis_meta["cell_type"] = metadata[args.cell_type_column].astype(str).str.strip()
        analysis_meta = analysis_meta.loc[analysis_meta["cell_type"] != ""].reset_index(drop=True)

        # First require enough cells to form a reasonable pseudobulk column for a
        # given sample x cell_type combination.
        sample_cell_counts = (
            analysis_meta.groupby(["orig.ident", "IntDx", "Sex", "Location", "cell_type"], as_index=False)
            .size()
            .rename(columns={"size": "n_cells"})
        )
        retained_sample_cell = sample_cell_counts.loc[
            sample_cell_counts["n_cells"] >= args.min_cells_per_sample
        ].copy()

        # Then require enough retained samples in each subtype so the DE model has
        # actual sample-level replication within that cell type.
        eligibility = (
            retained_sample_cell.groupby(["cell_type", "IntDx"])["orig.ident"]
            .nunique()
            .unstack(fill_value=0)
            .reindex(columns=list(DLGNT_LEVELS), fill_value=0)
        )
        totals = (
            retained_sample_cell.groupby(["cell_type", "IntDx"])["n_cells"]
            .sum()
            .unstack(fill_value=0)
            .reindex(columns=list(DLGNT_LEVELS), fill_value=0)
        )
        eligible_cell_types = (
            eligibility.rename(
                columns={"DLGNT_1": "n_samples_dlgnt1", "DLGNT_2": "n_samples_dlgnt2"}
            )
            .join(
                totals.rename(
                    columns={"DLGNT_1": "total_cells_dlgnt1", "DLGNT_2": "total_cells_dlgnt2"}
                )
            )
            .reset_index()
        )
        eligible_cell_types = eligible_cell_types.loc[
            (eligible_cell_types["n_samples_dlgnt1"] >= args.min_samples_per_subtype)
            & (eligible_cell_types["n_samples_dlgnt2"] >= args.min_samples_per_subtype)
        ].copy()
        eligible_cell_types["n_samples_total"] = (
            eligible_cell_types["n_samples_dlgnt1"] + eligible_cell_types["n_samples_dlgnt2"]
        )

        if eligible_cell_types.empty:
            raise ValueError(
                "No cell types met the requested sample support threshold after filtering."
            )

        safe_name_map = make_safe_names(eligible_cell_types["cell_type"].tolist())
        eligible_cell_types["cell_type_safe"] = eligible_cell_types["cell_type"].map(safe_name_map)
        eligible_cell_types = eligible_cell_types.sort_values(
            ["cell_type", "n_samples_total"], ascending=[True, False]
        ).reset_index(drop=True)
        write_tsv(eligible_cell_types, analysis_dir / "eligible_cell_types.tsv")

        # Keep only cells that belong to sample x cell_type combinations that
        # survive both filters. Those are the cells that contribute to pseudobulk.
        retained_keys = retained_sample_cell.merge(
            eligible_cell_types.loc[:, ["cell_type"]],
            on="cell_type",
            how="inner",
        ).loc[:, ["orig.ident", "cell_type"]].drop_duplicates()
        analysis_cells = analysis_meta.merge(retained_keys, on=["orig.ident", "cell_type"], how="inner")

        # This design table becomes the source of truth for pseudobulk column
        # ordering and for the aligned metadata consumed by the R stage.
        group_design = (
            retained_sample_cell.merge(
                eligible_cell_types.loc[:, ["cell_type", "cell_type_safe"]],
                on="cell_type",
                how="inner",
            )
            .sort_values(["cell_type", "IntDx", "orig.ident"])
            .reset_index(drop=True)
        )
        group_design["column_id"] = (
            group_design["orig.ident"].astype(str) + "__" + group_design["cell_type_safe"].astype(str)
        )
        if group_design["column_id"].duplicated().any():
            raise ValueError("column_id values must be unique across retained pseudobulk columns.")
        group_design["group_index"] = np.arange(group_design.shape[0], dtype=np.int32)

        group_lookup = group_design.loc[:, ["orig.ident", "cell_type", "column_id", "group_index"]]
        analysis_cells = analysis_cells.merge(group_lookup, on=["orig.ident", "cell_type"], how="left")
        if analysis_cells["group_index"].isna().any():
            raise ValueError("Failed to map selected cells to retained pseudobulk groups.")

        # row_to_group lets the sparse-count streamer jump directly from each cell
        # index to the output pseudobulk column it should be added into.
        row_to_group = np.full(int(obs_group["_index"].shape[0]), -1, dtype=np.int32)
        row_to_group[analysis_cells["cell_index"].to_numpy(dtype=np.int64)] = analysis_cells[
            "group_index"
        ].to_numpy(dtype=np.int32)

        agg = aggregate_pseudobulk(
            counts_group=handle["layers"]["counts"],
            row_to_group=row_to_group,
            n_genes=var_names.shape[0],
            n_groups=group_design.shape[0],
        )

    # Record the resolved runtime paths and cohort sizes so downstream summaries
    # do not need to infer how this run was configured.
    prepare_metadata = pd.DataFrame(
        {
            "key": [
                "workdir",
                "h5ad_path",
                "cell_type_column",
                "location_column",
                "min_samples_per_subtype",
                "min_cells_per_sample",
                "n_dlgnt_cells",
                "n_dlgnt_samples",
                "n_genes",
                "n_eligible_cell_types",
                "analysis_note",
            ],
            "value": [
                workdir,
                str(h5ad_path.resolve()),
                args.cell_type_column,
                location_key,
                str(args.min_samples_per_subtype),
                str(args.min_cells_per_sample),
                str(metadata.shape[0]),
                str(sample_metadata["orig.ident"].nunique()),
                str(var_names.shape[0]),
                str(eligible_cell_types.shape[0]),
                "Exploratory subtype-only analysis; subtype and location are confounded.",
            ],
        }
    )
    write_tsv(prepare_metadata, analysis_dir / "prepare_metadata.tsv")

    for row in eligible_cell_types.itertuples(index=False):
        cell_type = row.cell_type
        cell_type_safe = row.cell_type_safe
        group_subset = group_design.loc[group_design["cell_type"] == cell_type].copy()
        group_indices = group_subset["group_index"].to_numpy(dtype=int)
        # These diagnostics make it easy to spot pseudobulks dominated by very
        # low counts even when they passed the cell-count filter.
        group_subset["n_counts_total"] = agg[:, group_indices].sum(axis=0).astype(np.int64)
        group_subset["n_genes_nonzero"] = (agg[:, group_indices] > 0).sum(axis=0).astype(np.int64)

        counts_df = pd.DataFrame(agg[:, group_indices], columns=group_subset["column_id"].tolist())
        counts_df.insert(0, "gene", var_names)

        sample_diag = group_subset.loc[
            :,
            [
                "column_id",
                "orig.ident",
                "IntDx",
                "Sex",
                "Location",
                "cell_type",
                "cell_type_safe",
                "n_cells",
                "n_counts_total",
                "n_genes_nonzero",
            ],
        ].rename(columns={"orig.ident": "sample"})

        write_tsv(counts_df, handoff_dir / f"{cell_type_safe}_counts.tsv")
        write_tsv(
            sample_diag.loc[
                :,
                [
                    "column_id",
                    "sample",
                    "IntDx",
                    "Sex",
                    "Location",
                    "cell_type",
                    "cell_type_safe",
                    "n_cells",
                ],
            ],
            handoff_dir / f"{cell_type_safe}_design.tsv",
        )
        write_tsv(sample_diag, handoff_dir / f"{cell_type_safe}_sample_diagnostics.tsv")

    print(f"Wrote prepare outputs to: {analysis_dir}")
    print(f"DLGNT cells: {metadata.shape[0]:,}")
    print(f"DLGNT samples: {sample_metadata['orig.ident'].nunique():,}")
    print(f"Genes: {var_names.shape[0]:,}")
    print(f"Eligible cell types ({args.cell_type_column}): {eligible_cell_types.shape[0]:,}")


if __name__ == "__main__":
    main()
