#!/usr/bin/env Rscript
# Morphology-state pathway analysis stage (R), reticulate-free.
# - Reads Python handoff files (pseudobulk / design / shared universe).
# - Per target cell type: limma DE with ~ sample + region (compact = reference),
#   so the region coefficient reads as loose(myxoid)-vs-compact.
# - Preranked fgsea across GO:BP / REAC / KEGG / WP / Hallmark on the moderated t.
# - Tags each gene set as glycan/ECM vs SASP/senescence for downstream focus.
# - Writes CSV outputs consumed by 4_PA_morphology_state.qmd.
#
# Inferential note: the pseudobulk columns are sample x region means, so the
# region effect is estimated across usable samples (n = 3 PA samples with both
# regions; 258482 has no compact annotation). eBayes borrows strength across
# genes, but with ~2 residual df the FDRs are fragile -- always read them
# alongside cross-sample directional consistency (handled in the notebook).
# ORA is intentionally not included.

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
    stop("Package 'fgsea' is required for morphology-state GSEA. Install it before running this script.")
}

set.seed(5555)
options(stringsAsFactors = FALSE)

# Resolve output_dir from the script location so we read local handoff files.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
if (length(file_arg) < 1) {
    stop("Unable to determine script path from --file argument.")
}
script_path <- sub("^--file=", "", file_arg[1])
output_dir <- dirname(normalizePath(script_path))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_set_size <- 10L
max_set_size <- 500L

# Themes for downstream focus; these are the analysis question (myxoid = glycan-
# rich matrix; SASP ties to the BRAF-fusion senescence axis).
glycan_pattern <- "ECM|EXTRACELLULAR_MATRIX|GLYCOSYL|GLYCAN|PROTEOGLYCAN|HYALURON|GLYCOSAMINOGLYCAN|CHONDROITIN|HEPARAN|SIALYL|GALACTOSYL"
sasp_pattern   <- "SENESCEN|SASP|SECRETORY_PHENOTYPE"

coverage_path <- file.path(output_dir, "morphology_state_py_coverage_summary.csv")
overlap_path  <- file.path(output_dir, "morphology_state_py_pairwise_gene_overlap.csv")
shared_path   <- file.path(output_dir, "morphology_state_py_shared_universe.csv")
pb_path       <- file.path(output_dir, "morphology_state_py_pseudobulk_matrix.csv")
design_path   <- file.path(output_dir, "morphology_state_py_sample_design.csv")
counts_path   <- file.path(output_dir, "morphology_state_py_cell_counts.csv")

# Produced by morphology_state_prepare.py in the same folder.
required_inputs <- c(coverage_path, overlap_path, shared_path, pb_path, design_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
    stop("Missing Python handoff files: ", paste(missing_inputs, collapse = "; "))
}

coverage_dt        <- fread(coverage_path)
pairwise_overlap   <- fread(overlap_path)
shared_universe_dt <- fread(shared_path)
pb_dt              <- fread(pb_path)
design_dt          <- fread(design_path)
counts_dt          <- if (file.exists(counts_path)) fread(counts_path) else data.table()

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
if (anyDuplicated(design_dt$column_id)) {
    stop("Design contains duplicated column_id values.")
}

# ---- Expression filter + log transform (REQUIRED before limma) --------------
# The pseudobulk matrix holds arithmetic means of ENVI-decoded expression on a
# LINEAR scale (median ~3e-4). limma's eBayes(trend = TRUE) assumes LOG-scale data
# with a modelled mean-variance trend, so fitting the raw matrix deflates every
# statistic: the p-value distribution comes out humped in the middle instead of
# uniform, and the deflated moderated-t propagates straight into the fgsea ranking.
#
# Two corrections, applied per cell type inside the loop below:
#   1. drop genes with essentially no expression. Unexpressed genes are not
#      independent -- they arrive in large correlated families (olfactory and
#      vomeronasal receptors, keratins) whose near-random ranks cluster together and
#      manufacture spurious enrichment. In preranked GSEA the ranked list defines the
#      null background for EVERY set, so noise genes distort all of them, not just
#      their own. The filter is deliberately permissive so genuinely low-but-real
#      expression survives.
#   2. log2 transform, so the trend model sees the scale it expects.
MIN_MEAN_EXPR <- 0.01      # linear-scale floor, applied within each cell type
LOG_EPS       <- 1e-6      # well below the expressed-gene median (~3e-4)

filter_and_log <- function(mat) {
    keep <- rowMeans(mat) >= MIN_MEAN_EXPR
    log2(mat[keep, , drop = FALSE] + LOG_EPS)
}

# Ratio-to-uniform QC on a p-value vector: under a true null every entry is ~1.
# A hump in the middle bins means the test is conservative (the bug we fixed); a
# spike in the first bin beyond what the biology supports means anti-conservative.
pval_qc <- function(p, label) {
    br <- c(0, .05, .1, .2, .4, .6, .8, 1)
    obs <- as.numeric(table(cut(p, br, include.lowest = TRUE))) / length(p)
    data.table(label = label, bin = paste0("(", head(br, -1), ",", br[-1], "]"),
               frac = obs, expected = diff(br), ratio_to_uniform = obs / diff(br))
}

# ---- Build the MSigDB pathway list once (shared across cell types) ----------
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

# Restrict sets to genes actually testable (shared imputed universe) and size-filter.
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
pathways <- split(set_gene_filtered$gene, set_gene_filtered$pathway_key)
pathways <- lapply(pathways, unique)

# theme tag per set_name (used to focus on the glycan/SASP question downstream)
set_theme <- function(set_name) {
    fifelse(grepl(glycan_pattern, set_name, ignore.case = TRUE), "glycan_ecm",
        fifelse(grepl(sasp_pattern, set_name, ignore.case = TRUE), "sasp_senescence", "other"))
}

run_fgsea <- function(stats_vec, contrast_label, ct_label) {
    # minSize/maxSize re-applied HERE, against the actual ranked vector: the global
    # set_size_dt gate was computed on the unfiltered universe, so after the per-cell-type
    # expression filter a set can fall below the floor in this particular ranking.
    fg_res <- fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, eps = 0,
                                     minSize = min_set_size, maxSize = max_set_size)
    fg_dt <- as.data.table(fg_res)
    if (nrow(fg_dt) == 0L) {
        return(data.table())
    }
    parts <- tstrsplit(fg_dt$pathway, "|||", fixed = TRUE)
    fg_dt[, `:=`(
        source = parts[[1]],
        set_id = parts[[2]],
        set_name = parts[[3]],
        cell_type = ct_label,
        contrast = contrast_label,
        leading_edge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    )]
    fg_dt[, theme := set_theme(set_name)]
    fg_dt[, .(cell_type, contrast, source, set_id, set_name, theme, size, NES, pval, padj, leading_edge)]
}

# ---- Per cell type: limma (~ sample + region) + fgsea -----------------------
cell_types <- unique(design_dt$cell_type)
de_list <- list()
gsea_list <- list()
pqc_list <- list()

for (ct in cell_types) {
    ct_design <- design_dt[cell_type == ct]
    # Need >= 2 usable samples so the paired ~ sample + region design is estimable.
    n_samp <- ct_design[, uniqueN(sample)]
    if (n_samp < 2L) {
        message(sprintf("Skipping cell_type=%s: only %d usable sample(s).", ct, n_samp))
        next
    }
    ct_cols <- ct_design$column_id
    # filter + log2 within this cell type: a gene unexpressed in one lineage may be
    # perfectly real in another, so the gate must not be applied cohort-wide.
    ct_mat <- filter_and_log(pb_mat[, ct_cols, drop = FALSE])

    sample_factor <- factor(ct_design$sample)
    region_factor <- factor(ct_design$region, levels = c("compact", "loose"))  # compact = reference
    if (anyNA(region_factor)) {
        stop("Design contains regions outside expected values: compact / loose (cell_type=", ct, ").")
    }

    design_mat <- model.matrix(~ sample_factor + region_factor)
    fit <- limma::lmFit(ct_mat, design_mat)
    fit <- limma::eBayes(fit, trend = TRUE)
    coef_name <- "region_factorloose"
    if (!coef_name %in% colnames(design_mat)) {
        stop("Unable to find expected limma coefficient: ", coef_name, " (cell_type=", ct, ").")
    }

    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    de_dt <- data.table(
        cell_type = ct,
        gene = toupper(rownames(tt)),
        log2FC_loose_vs_compact = as.numeric(tt$logFC),
        t_stat = as.numeric(tt$t),
        pvalue = as.numeric(tt$P.Value),
        fdr = as.numeric(tt$adj.P.Val)
    )
    de_dt <- unique(de_dt, by = "gene")
    de_list[[ct]] <- de_dt

    # p-value distribution QC for this fit (acceptance test for the scale fix)
    pq <- pval_qc(de_dt$pvalue, ct)
    pq[, `:=`(n_genes_tested = nrow(de_dt), frac_p_lt_0.05 = mean(de_dt$pvalue < 0.05),
              min_pvalue = min(de_dt$pvalue), min_fdr = min(de_dt$fdr))]
    pqc_list[[ct]] <- pq


    # Preranked fgsea on the moderated t (both directions).
    stats_loose <- de_dt$t_stat
    names(stats_loose) <- de_dt$gene
    stats_loose <- sort(stats_loose, decreasing = TRUE)

    stats_compact <- -de_dt$t_stat
    names(stats_compact) <- de_dt$gene
    stats_compact <- sort(stats_compact, decreasing = TRUE)

    gsea_list[[paste0(ct, "_loose")]] <- run_fgsea(stats_loose, "loose > compact", ct)
    gsea_list[[paste0(ct, "_compact")]] <- run_fgsea(stats_compact, "compact > loose", ct)

    message(sprintf("cell_type=%s: %d samples, %d pseudobulk columns.", ct, n_samp, ncol(ct_mat)))
}

if (length(de_list) == 0L) {
    stop("No cell type had >= 2 usable samples; nothing to test.")
}

limma_out <- rbindlist(de_list, fill = TRUE)
gsea_out <- rbindlist(gsea_list, fill = TRUE)

coverage_out <- coverage_dt[, .(sample, n_tested_genes)]
pairwise_out <- pairwise_overlap[, .(sample_a, sample_b, n_intersection, n_union, jaccard)]
shared_out <- data.table(gene = shared_universe)
geneset_out <- set_size_dt[, .(source, set_id, set_name, set_size_in_shared_universe)]
geneset_out[, theme := set_theme(set_name)]

fwrite(coverage_out, file.path(output_dir, "morphology_state_coverage_summary.csv"))
fwrite(pairwise_out, file.path(output_dir, "morphology_state_pairwise_gene_overlap.csv"))
fwrite(shared_out, file.path(output_dir, "morphology_state_shared_universe.csv"))
fwrite(limma_out, file.path(output_dir, "morphology_state_limma_de.csv"))
fwrite(geneset_out, file.path(output_dir, "morphology_state_gene_sets_filtered.csv"))
fwrite(gsea_out, file.path(output_dir, "morphology_state_gsea_results.csv"))
fwrite(rbindlist(pqc_list, fill = TRUE),
       file.path(output_dir, "morphology_state_pvalue_qc.csv"))
if (nrow(counts_dt) > 0L) {
    fwrite(counts_dt, file.path(output_dir, "morphology_state_cell_counts.csv"))
}

metadata_path <- file.path(output_dir, "morphology_state_metadata.txt")
db_version <- if ("db_version" %in% names(msig_all)) unique(msig_all$db_version) else NA_character_
metadata_lines <- c(
    paste0("run_time_utc: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("cell_types_tested: ", paste(names(de_list), collapse = ", ")),
    paste0("n_samples_total: ", nrow(coverage_out)),
    paste0("shared_universe_n: ", length(shared_universe)),
    paste0("expression_filter_min_mean: ", MIN_MEAN_EXPR),
    paste0("transform: log2(x + ", LOG_EPS, ") applied per cell type after filtering"),
    paste0("genes_tested_per_cell_type: ",
           paste(vapply(de_list, nrow, integer(1)), collapse = ", ")),
    paste0("min_set_size: ", min_set_size),
    paste0("max_set_size: ", max_set_size),
    paste0("reference_region: compact"),
    paste0("msigdb_db_version: ", paste(na.omit(db_version), collapse = ";")),
    paste0("python_handoff_files: ", paste(basename(required_inputs), collapse = ", "))
)
writeLines(metadata_lines, metadata_path)

cat("Saved morphology-state outputs to:", output_dir, "\n")
cat("Cell types tested:", paste(names(de_list), collapse = ", "), "\n")
cat("Shared universe size:", length(shared_universe), "\n")
cat("DE rows:", nrow(limma_out), "\n")
cat("Filtered gene sets:", nrow(geneset_out), "\n")
cat("GSEA rows:", nrow(gsea_out), "\n")
