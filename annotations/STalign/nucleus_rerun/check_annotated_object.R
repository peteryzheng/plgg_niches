# Check the annotated BANKSY object written by annotations/STalign/add_annotations.R from the nucleus rerun.
#
# Inputs : the new object (..._pathology_nucleus.rds) and its 5% subset; the per-sample metas in
#          pathology_annotations_101326_nucleus; the cell-metadata TSV exported from the old (030226-based) object
#          (xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv), for the old labels.
# Checks : (1) every cell's pathology_annotation / region_id in the object equals its meta row (annotation recoded
#          with add_annotations.R's normalisation); (2) the subset is 5% of the cells with the same labels;
#          (3) per sample, how many cells changed label from the old object (all cells, and cells the old map
#          had placed in an annotated region), to compare with qc/label_impact.tsv.
# Output : NEW_DIR/qc/annotated_object_check.tsv; stops on any mismatch in (1) or (2).
# Usage  : Rscript annotations/STalign/nucleus_rerun/check_annotated_object.R   (conda env spatial)
suppressMessages({
    library(data.table)
    library(SummarizedExperiment)
})

data_dir <- "/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data"
new_dir <- file.path(data_dir, "Xenium_annotations/pathology_annotations_101326_nucleus")
obj_stem <- file.path(
    data_dir, "Xenium_annotations",
    "banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134_pathology_nucleus"
)
old_tsv <- file.path(data_dir, "xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv")

# the same recoding add_annotations.R applies (copied from there)
normalize_pathology_annotation <- function(annotation_values) {
    fcase(
        annotation_values == "#000000", "Blood",
        annotation_values == "#800000", "Compact. fibrillary component",
        annotation_values == "#ff0000", "Microvascular proliferation",
        annotation_values == "Compact, eosinophilic component", "Compact. fibrillary component",
        default = annotation_values
    )
}

# labels in the object, keyed like add_annotations.R: <cell_id>__<numeric sample id>
cd_of <- function(path) {
    cd <- as.data.table(as.data.frame(colData(readRDS(path))), keep.rownames = "merge_id")
    cd[, .(merge_id, pathology_annotation = as.character(pathology_annotation), region_id = as.character(region_id))]
}
obj <- cd_of(paste0(obj_stem, ".rds"))
sub <- cd_of(paste0(obj_stem, "_subset.rds"))

# labels in the metas
meta <- rbindlist(lapply(list.files(new_dir, pattern = "^proseg_pathology_.*_proseg_meta\\.csv$", full.names = TRUE), function(f) {
    dt <- fread(f, select = c("cell_id", "annotation", "region_id"), na.strings = c("", "NA"), colClasses = list(character = "region_id"))
    dt[, sample := sub("^proseg_pathology_([A-Z]{2}[0-9]+)_proseg_meta\\.csv$", "\\1", basename(f))]
    dt[, merge_id := paste0(cell_id, "__", sub("^[A-Z]{2}", "", sample))]
    dt[, .(merge_id, sample, meta_annotation = normalize_pathology_annotation(annotation), meta_region_id = region_id)]
}))

# (1) object vs metas, every cell
m <- merge(obj, meta, by = "merge_id", all.x = TRUE)
stopifnot(!anyNA(m$sample))
m[, ok := pathology_annotation == meta_annotation & (region_id == meta_region_id | (is.na(region_id) & is.na(meta_region_id)))]
if (!all(m$ok)) stop("object labels differ from the metas in ", sum(!m$ok), " cells")

# (2) the subset: 5% of cells, same labels
s <- merge(sub, obj, by = "merge_id", suffixes = c("_sub", ""))
stopifnot(nrow(s) == nrow(sub), nrow(sub) == floor(nrow(obj) * 0.05),
          all(s$pathology_annotation_sub == s$pathology_annotation))

# (3) change from the old object's labels
old <- fread(old_tsv, select = c("cell_id", "pathology_annotation"))
setnames(old, c("merge_id", "old_annotation"))
m <- merge(m, old, by = "merge_id", all.x = TRUE)
stopifnot(!anyNA(m$old_annotation))
res <- m[, .(
    n_cells = .N,
    n_old_annotated = sum(old_annotation != "unannotated"),
    n_new_annotated = sum(pathology_annotation != "unannotated"),
    share_all_relabelled = mean(pathology_annotation != old_annotation),
    share_old_annotated_relabelled = mean((pathology_annotation != old_annotation)[old_annotation != "unannotated"])
), by = sample][order(sample)]
res[, `:=`(object_matches_metas = TRUE, subset_cells = nrow(sub))]
fwrite(res, file.path(new_dir, "qc/annotated_object_check.tsv"), sep = "\t")
cat("object cells:", nrow(obj), "| subset cells:", nrow(sub), "| all labels match the metas\n")
print(res[, .(sample, n_cells, n_old_annotated, n_new_annotated,
              all = round(share_all_relabelled, 3), old_annotated = round(share_old_annotated_relabelled, 3))])
