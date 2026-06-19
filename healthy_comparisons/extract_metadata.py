#!/usr/bin/env python3
"""
Extract .obs metadata from the integrated healthy+tumor h5ad object into a
parquet file for downstream cell-type composition modeling.

Input:  combined_healthyTumor_scvi.h5ad (opened in backed mode; X not loaded)
Output: combined_healthyTumor_obs_metadata.parquet (all obs columns + cell_id)
"""

import os
import anndata

# Resolve data root from HOME — matches repo path-resolution convention in AGENTS.md.
home = os.path.expanduser("~")
if home.startswith(("/Users/youyun", "/Users/youyunzheng")):
    data_root = os.path.join(home, "Documents/HMS/PhD/beroukhimlab/dfci_mount")
elif home == "/PHShome/yz762":
    data_root = "/data/beroukhim1"
elif home == "/home/yz762":
    data_root = "/mnt/storage/dept/medonc/beroukhim"
else:
    data_root = "/xchip/beroukhimlab"

sc_dir = os.path.join(data_root, "youyun/plgg/data/single_cell")
h5ad_path = os.path.join(sc_dir, "combined_healthyTumor_scvi.h5ad")
out_path = os.path.join(sc_dir, "combined_healthyTumor_obs_metadata.parquet")

print(f"Reading:  {h5ad_path}")
print(f"Output:   {out_path}")

# backed='r' memory-maps X so only .obs is fully loaded into RAM.
adata = anndata.read_h5ad(h5ad_path, backed="r")

print(f"\nShape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
print(f"\nobs columns ({len(adata.obs.columns)}):")
print(adata.obs.dtypes.to_string())

# Promote cell barcode index to a plain column before saving.
obs = adata.obs.copy()
obs.index.name = "cell_id"
obs = obs.reset_index()

obs.to_parquet(out_path, index=False)
print(f"\nSaved {len(obs):,} rows → {out_path}")
