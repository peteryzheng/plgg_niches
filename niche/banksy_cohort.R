library(data.table)
library(ggplot2)

library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
library(harmony)
library(scran)

# local vs UGER
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}
source(paste0(workdir, "youyun/plgg/code/niche/spatial_helper.R"))

SEED = 55555

# LOAD DATA ===========================================================
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Loading data...'))
xenium_dirs = list.files(paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Analyzer/'
), full.names = TRUE)

metadata = data.table(readxl::read_excel(paste0(
    workdir,'coja/Spatial_PLGG/data/metadata/Xenium_PS.xlsx'
)))

total_se = do.call(
    'cbind',
    lapply(xenium_dirs, function(x){
        tmp_se = xenium2SPE(x)
        file_name = basename(x)
        sample_id = metadata[file == file_name]$id
        sample_idat = metadata[file == file_name]$idat
        histology = metadata[file == file_name]$mc
        alteration = metadata[file == file_name]$alt
        tmp_se$sample_id = sample_id
        tmp_se$sample_idat = sample_idat
        tmp_se$histology = histology
        tmp_se$alteration = alteration
        colnames(tmp_se) = paste0(colnames(tmp_se), '__', sample_id)
        return(tmp_se)
    }
))


# STAGGER SPATIAL COORDINATES ==========================================
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Staggering spatial coordinates...'))
total_se_staggered = stagger_spatial_coords(total_se)

# SUBSET to only gene expression and not control probes ===============
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Subsetting to gene expression...'))
total_se_staggered = total_se_staggered[rowData(total_se_staggered)$feature_type == 'Gene Expression',]

# QC and NORMALIZATION =================================================
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','QC and normalization...'))
qcstats <- perCellQCMetrics(total_se_staggered)
# QC based on total counts
# inspired by the MAD idea from single cell best practice 
# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-cells
is_outlier = function(metric, nmad = 3){
    metric = as.numeric(metric)
    mad = median(abs(metric - median(metric)))
    abs(metric - median(metric)) > nmad * mad
}
outliers = is_outlier(log10(qcstats$total + 1)) | is_outlier(log10(qcstats$detected + 1))
table(outliers)
keep <- !outliers
total_se_staggered <- total_se_staggered[, keep]

# Normalization to mean library size
total_se_staggered <- computeLibraryFactors(total_se_staggered)
aname <- "normcounts"
assay(total_se_staggered, aname) <- normalizeCounts(total_se_staggered, log = FALSE)

#  BANKSY ==============================================================
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Running BANKSY...'))
# default k_geom
k_geom <- c(15, 30)
# calculate mean neighborhood and AGF matricies
total_se_staggered <- Banksy::computeBanksy(total_se_staggered, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

# weighted embedding
# 0 for non spatial clustering, 0.2 for cell typing, and 0.8 for domain segmentation
lambda = c(0, 0.2, 0.8)
total_se_staggered <- Banksy::runBanksyPCA(total_se_staggered, use_agf = TRUE, lambda = lambda, seed = SEED)

# Batch correction =====================================================
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Batch correction...'))
set.seed(SEED)
harmony_embedding_PCA_M1_lam0 <- HarmonyMatrix(
    data_mat = reducedDim(total_se_staggered, "PCA_M1_lam0"),
    meta_data = colData(total_se_staggered),
    vars_use = c("sample_id"),
    do_pca = FALSE,
    max_iter = 50,
    verbose = TRUE
)
reducedDim(total_se_staggered, "Harmony_BANKSY_lam0") <- harmony_embedding_PCA_M1_lam0

harmony_embedding_PCA_M1_lam0.2 <- HarmonyMatrix(
    data_mat = reducedDim(total_se_staggered, "PCA_M1_lam0.2"),
    meta_data = colData(total_se_staggered),
    vars_use = c("sample_id"),
    do_pca = FALSE,
    max_iter = 50,
    verbose = TRUE
)
reducedDim(total_se_staggered, "Harmony_BANKSY_lam0.2") <- harmony_embedding_PCA_M1_lam0.2

harmony_embedding_PCA_M1_lam0.8 <- HarmonyMatrix(
    data_mat = reducedDim(total_se_staggered, "PCA_M1_lam0.8"),
    meta_data = colData(total_se_staggered),
    vars_use = c("sample_id"),
    do_pca = FALSE,
    max_iter = 50,
    verbose = TRUE
)
reducedDim(total_se_staggered, "Harmony_BANKSY_lam0.8") <- harmony_embedding_PCA_M1_lam0.8

# UMAP ================================================================
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Running UMAP...'))
# uncorrected
total_se_staggered <- Banksy::runBanksyUMAP(total_se_staggered, use_agf = TRUE, lambda = lambda)
# corrected
total_se_staggered <- runBanksyUMAP(total_se_staggered, dimred = "Harmony_BANKSY_lam0")
total_se_staggered <- runBanksyUMAP(total_se_staggered, dimred = "Harmony_BANKSY_lam0.2")
total_se_staggered <- runBanksyUMAP(total_se_staggered, dimred = "Harmony_BANKSY_lam0.8")

total_se_staggered_subset = total_se_staggered[, sample(
    c(1:ncol(total_se_staggered)),ncol(total_se_staggered) * 0.05
)]

current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

saveRDS(total_se_staggered, paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_embeddings_',
    current_time,'.rds'
))
saveRDS(total_se_staggered_subset, paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_embeddings_subset_',
    current_time,'.rds'
))

# Leiden and connect clusters =====================================
# total_se_staggered <- Banksy::clusterBanksy(total_se_staggered, use_agf = TRUE, lambda = lambda, resolution = 1.2, seed = SEED)
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Clustering...'))
total_se_staggered <- Banksy::clusterBanksy(
    total_se_staggered, dimred = "Harmony_BANKSY_lam0", 
    resolution = c(1, 0.5), seed = SEED
)
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finished clustering lam0 ...'))
total_se_staggered <- Banksy::clusterBanksy(
    total_se_staggered, dimred = "Harmony_BANKSY_lam0.2", 
    resolution = c(1, 0.5), seed = SEED
)
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finished clustering lam0.2 ...'))
total_se_staggered <- Banksy::clusterBanksy(
    total_se_staggered, dimred = "Harmony_BANKSY_lam0.8", 
    resolution = c(1, 0.5), seed = SEED
)
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finished clustering lam0.8 ...'))


saveRDS(total_se_staggered, paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_clusters_',
    current_time,'.rds'
))

# Different clustering runs can be relabeled to minimise their differences with connectClusters:
total_se_staggered <- Banksy::connectClusters(total_se_staggered, map_to = clusterNames(total_se_staggered)[1])
print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finished connecting clusters'))

saveRDS(total_se_staggered, paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_clusters_',
    current_time,'.rds'
))

# total_se_staggered = readRDS('/xchip/beroukhimlab/coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_clusters_20240711_135041.rds')
total_se_staggered_subset = total_se_staggered[, sample(
    c(1:ncol(total_se_staggered)),ncol(total_se_staggered) * 0.05
)]
# saveRDS(
#     total_se_staggered_subset,
#     '/xchip/beroukhimlab/coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_clusters_subset_20240711_135041.rds'
# )
saveRDS(total_se_staggered_subset, paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_clusters_subset_',
    current_time,'.rds'
))

cell_type_markers = findMarkers(
    assay(total_se_staggered, "counts"),
    groups = total_se_staggered$clust_Harmony_BANKSY_lam0.2_k50_res1.2,
    test.type="wilcox"
)

saveRDS(cell_type_markers, paste0(
    workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Objects/total_spatial_banksy_clusters_markers_',
    current_time,'.rds'
))