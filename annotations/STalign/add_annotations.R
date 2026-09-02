library(data.table)
library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)


# local vs remote
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
} else if (Sys.getenv("HOME") == "/PHShome/yz762"){
    # in eristwo, the home directory is at /PHShome/yz762
    workdir <- "/data/beroukhim1/"
} else if (Sys.getenv("HOME") == "/home/yz762"){
    # on the new server, the home directory is at /home/yz762
    workdir <- "/mnt/storage/dept/medonc/beroukhim/"
} else {
    # in dipg or uger, the home directory is '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}

banksy_rds_path <- file.path(
    workdir,
    "youyun/plgg/data/banksy_param_search/k1_15_k2_30_lambda1_0.2_lambda2_0.8_npcs_20_kct_50_resct_0.5,1_kni_50_resni_0.5,1",
    "banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134.rds"
)
annotation_input_dir <- file.path(
    workdir,
    "youyun/plgg/data/Xenium_annotations/pathology_annotations_030226"
)
annotated_rds_path <- file.path(
    workdir,
    "youyun/plgg/data/Xenium_annotations",
    "banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134_pathology.rds"
)
annotated_subset_rds_path <- file.path(
    workdir,
    "youyun/plgg/data/Xenium_annotations",
    "banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134_pathology_subset.rds"
)

parse_sample_id <- function(path) {
    sample_id <- sub(
        "^proseg_pathology_[A-Z]{2}([0-9]+)_proseg_meta\\.csv$",
        "\\1",
        basename(path)
    )
    if (identical(sample_id, basename(path))) {
        stop("Unable to parse sample_id from file: ", basename(path))
    }
    sample_id
}

normalize_pathology_annotation <- function(annotation_values) {
    fcase(
        # some sample specific encoding
        annotation_values == "#000000", "Blood",
        annotation_values == "#800000", "Compact. fibrillary component",
        annotation_values == "#ff0000", "Microvascular proliferation",
        annotation_values == "Compact, eosinophilic component", "Compact. fibrillary component",
        default = annotation_values
    )
}

banksy_embeddings <- readRDS(banksy_rds_path)
original_metadata <- as.data.frame(colData(banksy_embeddings))
original_order <- rownames(original_metadata)
original_metadata$sample_id_chr <- as.character(original_metadata$sample_id)
original_metadata$merge_id <- paste0(
    original_metadata$cell_id, "__", original_metadata$sample_id_chr
)

annotation_cols_to_replace <- intersect(
    c("pathology_annotation", "region_id"),
    colnames(original_metadata)
)
if (length(annotation_cols_to_replace) > 0) {
    original_metadata[, annotation_cols_to_replace] <- NULL
}

proseg_files <- list.files(
    annotation_input_dir,
    full.names = TRUE,
    pattern = "^proseg_pathology_.*_proseg_meta\\.csv$"
)
if (length(proseg_files) == 0) {
    stop("No annotation CSVs found in: ", annotation_input_dir)
}

proseg_annotations <- rbindlist(lapply(proseg_files, function(f) {
    dt <- fread(
        f,
        select = c("cell_id", "annotation", "region_id"),
        na.strings = c("", "NA")
    )
    dt[, sample_id := parse_sample_id(f)]
    setnames(dt, "annotation", "pathology_annotation")
    dt[, pathology_annotation := normalize_pathology_annotation(pathology_annotation)]
    dt[, merge_id := paste0(cell_id, "__", sample_id)]
    dt[, .(cell_id, sample_id, merge_id, pathology_annotation, region_id)]
}))

if (anyDuplicated(proseg_annotations$merge_id) > 0) {
    duplicate_keys <- unique(proseg_annotations$merge_id[duplicated(proseg_annotations$merge_id)])
    stop(
        "Duplicate annotation keys detected. Example duplicated merge_id values: ",
        paste(head(duplicate_keys, 5), collapse = ", ")
    )
}

object_sample_ids <- sort(unique(original_metadata$sample_id_chr))
annotation_sample_ids <- sort(unique(proseg_annotations$sample_id))
missing_annotation_files <- setdiff(object_sample_ids, annotation_sample_ids)
extra_annotation_files <- setdiff(annotation_sample_ids, object_sample_ids)

if (length(missing_annotation_files) > 0) {
    stop(
        "Missing annotation files for sample_id values: ",
        paste(missing_annotation_files, collapse = ", ")
    )
}
if (length(extra_annotation_files) > 0) {
    warning(
        "Annotation files found for sample_id values not present in BANKSY object: ",
        paste(extra_annotation_files, collapse = ", ")
    )
}

annotation_lookup <- unique(proseg_annotations[, .(
    merge_id, pathology_annotation, region_id
)])
new_metadata <- merge(
    original_metadata,
    annotation_lookup,
    by = "merge_id",
    all.x = TRUE,
    sort = FALSE
)
rownames(new_metadata) <- new_metadata$merge_id
new_metadata <- new_metadata[original_order, ]

if (!all(rownames(new_metadata) == original_order)) {
    stop("Merged metadata row order no longer matches the original BANKSY object.")
}
if (any(is.na(new_metadata$pathology_annotation))) {
    missing_rows <- rownames(new_metadata)[is.na(new_metadata$pathology_annotation)]
    stop(
        "Merge left cells without pathology annotations. Example merge_id values: ",
        paste(head(missing_rows, 5), collapse = ", ")
    )
}

join_summary <- as.data.table(new_metadata)[, .(
    n_cells = .N,
    n_with_pathology_annotation = sum(!is.na(pathology_annotation)),
    n_unannotated = sum(pathology_annotation == "unannotated", na.rm = TRUE),
    n_with_region_id = sum(!is.na(region_id))
), by = sample_id][order(sample_id)]

message("Annotation files consumed: ", length(proseg_files))
message("Annotation sample_id values: ", paste(annotation_sample_ids, collapse = ", "))
print(join_summary)

print(head(new_metadata))
new_metadata$sample_id_chr <- NULL
new_metadata$merge_id <- NULL
new_metadata <- S4Vectors::DataFrame(
    new_metadata,
    row.names = rownames(new_metadata)
)
colData(banksy_embeddings) <- new_metadata

saveRDS(banksy_embeddings, file = annotated_rds_path)

subset_indices <- sample(1:ncol(banksy_embeddings), ncol(banksy_embeddings) * 0.05)
banksy_embeddings_subset <- banksy_embeddings[, subset_indices]
saveRDS(banksy_embeddings_subset, file = annotated_subset_rds_path)
