# Leiden resolution/k sweep on the saved BANKSY subset object.
#
# Purpose: pick data-backed defaults for cell-type (lam0.2) and niche (lam0.8)
# Leiden clustering by running a (k_leiden x resolution) grid on the
# Harmony-corrected BANKSY PCs that are already stored in the subset object.
# We do not rerun BANKSY/Harmony here -- only Leiden -- so this is cheap.
#
# Inputs : the 90k-cell subset RDS from the prior banksy_cohort run.
# Outputs: a wide table of n_clusters per (lambda, k_leiden, resolution)
#          printed to stdout and written as TSV next to this script.

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(SpatialExperiment)
    library(data.table)
    library(Banksy)
})

# Path resolution mirrors the AGENTS.md convention so the script works on
# both the local mac and the cluster.
home <- Sys.getenv("HOME")
if (home %in% c("/Users/youyun", "/Users/youyunzheng")) {
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
    # The subset object lives outside dfci_mount on local; fall back to the
    # explicit local-only path used in banksy_clusters_proseg.qmd.
    subset_path <- "/Users/youyun/Documents/HMS/PhD/beroukhimlab/plgg/data/banksy_clusters_connected_subset_k_geom_15_30_pc_20_lam_0.2_0.8_k_leiden_30_50_res_0.75_1_20241230_232215.rds"
} else if (home == "/PHShome/yz762") {
    workdir <- "/data/beroukhim1/"
    subset_path <- paste0(workdir, "youyun/plgg/data/banksy_param_search/k1_15_k2_30_lambda1_0.2_lambda2_0.8_npcs_20_kc1_30_kc2_50_res1_0.75_res2_1/banksy_clusters_connected_subset_k_geom_15_30_pc_20_lam_0.2_0.8_k_leiden_30_50_res_0.75_1_20241230_232215.rds")
} else if (home == "/home/yz762") {
    workdir <- "/mnt/storage/dept/medonc/beroukhim/"
    subset_path <- paste0(workdir, "youyun/plgg/data/banksy_param_search/k1_15_k2_30_lambda1_0.2_lambda2_0.8_npcs_20_kc1_30_kc2_50_res1_0.75_res2_1/banksy_clusters_connected_subset_k_geom_15_30_pc_20_lam_0.2_0.8_k_leiden_30_50_res_0.75_1_20241230_232215.rds")
} else {
    workdir <- "/xchip/beroukhimlab/"
    subset_path <- paste0(workdir, "youyun/plgg/data/banksy_param_search/k1_15_k2_30_lambda1_0.2_lambda2_0.8_npcs_20_kc1_30_kc2_50_res1_0.75_res2_1/banksy_clusters_connected_subset_k_geom_15_30_pc_20_lam_0.2_0.8_k_leiden_30_50_res_0.75_1_20241230_232215.rds")
}

cat(sprintf("[%s] Loading: %s\n", format(Sys.time()), subset_path))
o <- readRDS(subset_path)
cat(sprintf("Cells: %d ; reducedDims: %s\n",
            ncol(o), paste(reducedDimNames(o), collapse = ", ")))

# Drop existing clust_* columns so the new clustering grid lands in clean
# colData slots and parsing column names is unambiguous.
existing_clust <- grep("^clust_", colnames(colData(o)), value = TRUE)
if (length(existing_clust) > 0) colData(o)[, existing_clust] <- NULL

# Sweep grid. Same grid applied to both lambdas; analysis below splits them.
k_vec   <- c(50, 80, 120)
res_vec <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
seed_val <- 55555

set.seed(seed_val)

# Run Leiden once per lambda; clusterBanksy expands the k x res grid internally.
for (lam in c(0.2, 0.8)) {
    cat(sprintf("[%s] Clustering lam%s (k x res grid = %d combos) ...\n",
                format(Sys.time()), lam, length(k_vec) * length(res_vec)))
    o <- Banksy::clusterBanksy(
        o, dimred = paste0("Harmony_BANKSY_lam", lam),
        k_neighbors = k_vec,
        resolution  = res_vec,
        algo = "leiden",
        seed = seed_val
    )
}

# Parse k and res from column names like 'clust_Harmony_BANKSY_lam0.2_k50_res0.4'.
clust_cols <- grep("^clust_Harmony_BANKSY_lam(0\\.2|0\\.8)_k\\d+_res",
                   colnames(colData(o)), value = TRUE)
results <- rbindlist(lapply(clust_cols, function(cc) {
    data.table(
        cluster_col = cc,
        lambda  = as.numeric(sub(".*_lam([0-9.]+)_k.*",  "\\1", cc)),
        k       = as.integer(sub(".*_k([0-9]+)_res.*",    "\\1", cc)),
        res     = as.numeric(sub(".*_res([0-9.]+)$",      "\\1", cc)),
        n_clusters = length(unique(colData(o)[[cc]]))
    )
}))

cat("\n=== Cluster counts ===\n")
# Wide table per lambda for human reading.
wide_ct <- dcast(results[lambda == 0.2], k ~ res, value.var = "n_clusters")
wide_ni <- dcast(results[lambda == 0.8], k ~ res, value.var = "n_clusters")
cat("\nlam0.2 (cell type):\n"); print(wide_ct)
cat("\nlam0.8 (niche):\n");      print(wide_ni)

# Resolve the directory of this script so the output lands next to it,
# regardless of where Rscript was invoked from.
script_dir <- tryCatch(
    dirname(normalizePath(sys.frames()[[1]]$ofile)),
    error = function(e) getwd()
)
out_tsv <- file.path(script_dir, "leiden_resolution_sweep_counts.tsv")
fwrite(results, out_tsv, sep = "\t")
cat(sprintf("\n[%s] Wrote %s\n", format(Sys.time()), out_tsv))
