#!/usr/bin/env Rscript
# Reticulate-free spatial-imputation pathway analysis stage (R):
# - Reads Python handoff files (coverage/shared universe/pseudobulk/design)
# - Runs limma DE with ~ sample + group
# - Runs preranked fgsea across GO:BP/REAC/KEGG/WP/Hallmark
# - Writes CSV outputs used by 3_PA_myeloids.qmd
# Note: ORA is intentionally not included.

suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(msigdbr)
})

required_pkgs <- c("data.table", "limma", "msigdbr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for spatial-imputation GSEA. Install it before running this script.")
}

set.seed(5555)
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

detect_workdir <- function() {
    home <- Sys.getenv("HOME")
    if (home %in% c("/Users/youyun", "/Users/youyunzheng")) {
        return(path.expand("~/Documents/HMS/PhD/beroukhimlab/dfci_mount"))
    }
    if (home == "/PHShome/yz762") {
        return("/data/beroukhim1")
    }
    if (home == "/home/yz762") {
        return("/mnt/storage/dept/medonc/beroukhim")
    }
    return("/data/beroukhim1")
}

workdir <- detect_workdir()
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
if (length(file_arg) < 1) {
    stop("Unable to determine script path from --file argument.")
}
script_path <- sub("^--file=", "", file_arg[1])
# Resolve output_dir from script location so this stage reads local handoff files.
output_dir <- dirname(normalizePath(script_path))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_set_size <- 10L
max_set_size <- 500L

coverage_path <- file.path(output_dir, "myeloid_spatial_imputation_py_coverage_summary.csv")
overlap_path <- file.path(output_dir, "myeloid_spatial_imputation_py_pairwise_gene_overlap.csv")
shared_path <- file.path(output_dir, "myeloid_spatial_imputation_py_shared_universe.csv")
pb_path <- file.path(output_dir, "myeloid_spatial_imputation_py_pseudobulk_matrix.csv")
design_path <- file.path(output_dir, "myeloid_spatial_imputation_py_sample_design.csv")
diag_path <- file.path(output_dir, "myeloid_spatial_imputation_py_sample_diagnostics.csv")

# These are produced by myeloid_spatial_imputation_prepare.py in the same folder.
required_inputs <- c(coverage_path, overlap_path, shared_path, pb_path, design_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
    stop("Missing Python handoff files: ", paste(missing_inputs, collapse = "; "))
}

coverage_dt <- fread(coverage_path)
pairwise_overlap <- fread(overlap_path)
shared_universe_dt <- fread(shared_path)
pb_dt <- fread(pb_path)
design_dt <- fread(design_path)
py_diag_dt <- if (file.exists(diag_path)) fread(diag_path) else data.table()

if (!"gene" %in% names(pb_dt)) {
    stop("Pseudobulk matrix must include a 'gene' column.")
}

shared_universe <- unique(toupper(shared_universe_dt$gene))
if (length(shared_universe) == 0L) {
    stop("Shared universe is empty in Python handoff.")
}

pb_dt[, gene := toupper(gene)]
pb_dt <- unique(pb_dt, by = "gene")
pb_dt <- pb_dt[gene %in% shared_universe]
if (nrow(pb_dt) != length(shared_universe)) {
    stop("Pseudobulk matrix gene rows do not match shared universe size.")
}

pb_mat <- as.matrix(pb_dt[, setdiff(names(pb_dt), "gene"), with = FALSE])
rownames(pb_mat) <- pb_dt$gene

if (!all(design_dt$column_id %in% colnames(pb_mat))) {
    stop("Design contains column_id values missing from pseudobulk matrix.")
}
if (!all(colnames(pb_mat) %in% design_dt$column_id)) {
    stop("Pseudobulk matrix contains columns missing from design table.")
}
if (anyDuplicated(design_dt$column_id)) {
    stop("Design contains duplicated column_id values.")
}

ord <- match(colnames(pb_mat), design_dt$column_id)
if (anyNA(ord)) {
    stop("Failed to map pseudobulk columns to design table ordering.")
}
design_dt <- design_dt[ord]
pb_mat <- pb_mat[, design_dt$column_id, drop = FALSE]
if (!all(design_dt$column_id == colnames(pb_mat))) {
    stop("Design rows and pseudobulk columns are not aligned after reindexing.")
}

if (nrow(design_dt) != ncol(pb_mat)) {
    stop("Design rows and pseudobulk columns are not one-to-one.")
}

sample_factor <- factor(design_dt$sample)
group_factor <- factor(design_dt$group, levels = c("Myeloid 1", "Myeloid 2"))
if (anyNA(group_factor)) {
    stop("Design contains groups outside expected values: Myeloid 1 / Myeloid 2.")
}

# Model Myeloid 2 vs Myeloid 1 while adjusting for sample effects.
design <- model.matrix(~ sample_factor + group_factor)
fit <- limma::lmFit(pb_mat, design)
fit <- limma::eBayes(fit, trend = TRUE)
coef_name <- "group_factorMyeloid 2"
if (!coef_name %in% colnames(design)) {
    stop("Unable to find expected limma coefficient: ", coef_name)
}

tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
de_dt <- data.table(
    gene = toupper(rownames(tt)),
    logFC_m2_vs_m1 = as.numeric(tt$logFC),
    t_stat = as.numeric(tt$t),
    pvalue = as.numeric(tt$P.Value),
    fdr = as.numeric(tt$adj.P.Val)
)
de_dt <- unique(de_dt, by = "gene")

if (length(unique(de_dt$gene)) != length(shared_universe)) {
    stop("DE row count must equal shared universe size.")
}
if (!all(de_dt$gene %in% shared_universe)) {
    stop("DE genes include members outside shared universe.")
}

msig_all <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens"))
collection_col <- if ("gs_collection" %in% names(msig_all)) "gs_collection" else "gs_cat"
subcollection_col <- if ("gs_subcollection" %in% names(msig_all)) "gs_subcollection" else "gs_subcat"
gene_col <- if ("gene_symbol" %in% names(msig_all)) "gene_symbol" else "human_gene_symbol"
set_id_col <- if ("gs_id" %in% names(msig_all)) "gs_id" else "gs_name"
set_name_col <- "gs_name"

subset_sources <- function(dt, source_name) {
    if (source_name == "H") {
        return(dt[get(collection_col) == "H"])
    }
    if (source_name == "GO:BP") {
        return(dt[get(collection_col) == "C5" & grepl("GO:BP", get(subcollection_col), ignore.case = TRUE)])
    }
    if (source_name == "REAC") {
        return(dt[get(collection_col) == "C2" & grepl("REACTOME", get(subcollection_col), ignore.case = TRUE)])
    }
    if (source_name == "KEGG") {
        return(dt[get(collection_col) == "C2" & grepl("KEGG", get(subcollection_col), ignore.case = TRUE)])
    }
    if (source_name == "WP") {
        return(dt[get(collection_col) == "C2" & grepl("WIKIPATHWAYS|WP", get(subcollection_col), ignore.case = TRUE)])
    }
    data.table()
}

source_levels <- c("GO:BP", "REAC", "KEGG", "WP", "H")
source_tables <- lapply(source_levels, function(src) {
    src_dt <- subset_sources(msig_all, src)
    if (nrow(src_dt) == 0L) {
        return(data.table())
    }
    unique(src_dt[, .(
        source = src,
        set_id = as.character(get(set_id_col)),
        set_name = as.character(get(set_name_col)),
        gene = toupper(as.character(get(gene_col)))
    )])
})
set_gene_dt <- rbindlist(source_tables, fill = TRUE)
set_gene_dt <- set_gene_dt[!is.na(gene) & gene != ""]
set_gene_dt <- unique(set_gene_dt)

set_gene_dt <- set_gene_dt[gene %in% shared_universe]
set_size_dt <- set_gene_dt[, .(set_size_in_shared_universe = uniqueN(gene)), by = .(source, set_id, set_name)]
set_size_dt <- set_size_dt[
    set_size_in_shared_universe >= min_set_size &
        set_size_in_shared_universe <= max_set_size
]

if (nrow(set_size_dt) == 0L) {
    stop("No gene sets remain after shared-universe intersection and size filtering.")
}

set_gene_filtered <- set_gene_dt[set_size_dt, on = .(source, set_id, set_name)]
set_gene_filtered[, pathway_key := paste(source, set_id, set_name, sep = "|||")]
set_size_dt[, pathway_key := paste(source, set_id, set_name, sep = "|||")]

pathways <- split(set_gene_filtered$gene, set_gene_filtered$pathway_key)
pathways <- lapply(pathways, unique)

if (!all(set_size_dt$set_size_in_shared_universe >= min_set_size & set_size_dt$set_size_in_shared_universe <= max_set_size)) {
    stop("Set-size filter validation failed.")
}

run_fgsea <- function(stats_vec, contrast_label) {
    fg_res <- fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, eps = 0)
    fg_dt <- as.data.table(fg_res)
    if (nrow(fg_dt) == 0L) {
        return(data.table(
            contrast = character(),
            source = character(),
            set_id = character(),
            set_name = character(),
            size = integer(),
            NES = numeric(),
            pval = numeric(),
            padj = numeric(),
            leading_edge = character()
        ))
    }
    parts <- tstrsplit(fg_dt$pathway, "|||", fixed = TRUE)
    fg_dt[, `:=`(
        source = parts[[1]],
        set_id = parts[[2]],
        set_name = parts[[3]],
        contrast = contrast_label,
        leading_edge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    )]
    fg_dt[, .(contrast, source, set_id, set_name, size, NES, pval, padj, leading_edge)]
}

stats_m2 <- de_dt$t_stat
names(stats_m2) <- de_dt$gene
stats_m2 <- sort(stats_m2, decreasing = TRUE)

stats_m1 <- -de_dt$t_stat
names(stats_m1) <- de_dt$gene
stats_m1 <- sort(stats_m1, decreasing = TRUE)

gsea_dt <- rbindlist(list(
    run_fgsea(stats_m2, "Myeloid 2 > Myeloid 1"),
    run_fgsea(stats_m1, "Myeloid 1 > Myeloid 2")
), fill = TRUE)

coverage_out <- coverage_dt[, .(sample, n_tested_genes)]
pairwise_out <- pairwise_overlap[, .(sample_a, sample_b, n_intersection, n_union, jaccard)]
shared_out <- data.table(gene = shared_universe)
limma_out <- de_dt[, .(gene, logFC_m2_vs_m1, t_stat, pvalue, fdr)]
geneset_out <- set_size_dt[, .(source, set_id, set_name, set_size_in_shared_universe)]
gsea_out <- gsea_dt[, .(contrast, source, set_id, set_name, size, NES, pval, padj, leading_edge)]

if (nrow(py_diag_dt) > 0L) {
    sample_diag_out <- merge(py_diag_dt, coverage_out, by = c("sample", "n_tested_genes"), all.y = TRUE)
} else {
    sample_diag_out <- coverage_out[, .(sample, n_tested_genes)]
}

fwrite(coverage_out, file.path(output_dir, "myeloid_spatial_imputation_coverage_summary.csv"))
fwrite(pairwise_out, file.path(output_dir, "myeloid_spatial_imputation_pairwise_gene_overlap.csv"))
fwrite(shared_out, file.path(output_dir, "myeloid_spatial_imputation_shared_universe.csv"))
fwrite(limma_out, file.path(output_dir, "myeloid_spatial_imputation_limma_de.csv"))
fwrite(geneset_out, file.path(output_dir, "myeloid_spatial_imputation_gene_sets_filtered.csv"))
fwrite(gsea_out, file.path(output_dir, "myeloid_spatial_imputation_gsea_results.csv"))
fwrite(sample_diag_out, file.path(output_dir, "myeloid_spatial_imputation_sample_diagnostics.csv"))

metadata_path <- file.path(output_dir, "myeloid_spatial_imputation_metadata.txt")
db_version <- if ("db_version" %in% names(msig_all)) unique(msig_all$db_version) else NA_character_
metadata_lines <- c(
    paste0("run_time_utc: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("n_samples_total: ", nrow(coverage_out)),
    paste0("n_samples_used: ", design_dt[, uniqueN(sample)]),
    paste0("shared_universe_n: ", length(shared_universe)),
    paste0("min_set_size: ", min_set_size),
    paste0("max_set_size: ", max_set_size),
    paste0("msigdb_db_version: ", paste(na.omit(db_version), collapse = ";")),
    paste0("python_handoff_files: ", paste(basename(required_inputs), collapse = ", "))
)
writeLines(metadata_lines, metadata_path)

cat("Saved spatial-imputation outputs to:", output_dir, "\n")
cat("Shared universe size:", length(shared_universe), "\n")
cat("DE genes:", nrow(limma_out), "\n")
cat("Filtered gene sets:", nrow(geneset_out), "\n")
cat("GSEA rows:", nrow(gsea_out), "\n")
