#!/usr/bin/env Rscript
# DLGNT_1 vs DLGNT_2 DE + Hallmark fgsea stage.
# - Reads raw pseudobulk handoffs written by dlgnt12_prepare.py
# - Runs edgeR filterByExpr + TMM + voom/limma per retained cell type
# - Runs preranked Hallmark fgsea on the tested-gene universe for each cell type
# - Writes one DE table, tested-gene table, and fgsea table per cell type

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(limma)
    library(msigdbr)
})

required_pkgs <- c("data.table", "edgeR", "limma", "msigdbr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for dlgnt12 fgsea.")
}

options(stringsAsFactors = FALSE)
set.seed(5555)

`%||%` <- function(a, b) if (is.null(a)) b else a

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
if (length(file_arg) < 1) {
    stop("Unable to determine script path from --file argument.")
}
# Resolve paths relative to the script so the wrapper can run from either local
# or cluster working directories without editing hardcoded analysis paths.
script_path <- sub("^--file=", "", file_arg[1])
script_dir <- dirname(normalizePath(script_path))

get_arg_value <- function(flag, default = NULL) {
    prefix <- paste0(flag, "=")
    hit <- args_all[startsWith(args_all, prefix)]
    if (length(hit) < 1) {
        return(default)
    }
    sub(prefix, "", hit[1], fixed = TRUE)
}

analysis_dir <- normalizePath(get_arg_value("--indir", script_dir), mustWork = FALSE)
handoff_dir <- file.path(analysis_dir, "handoff")
results_dir <- normalizePath(get_arg_value("--outdir", file.path(analysis_dir, "results")), mustWork = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

min_set_size <- as.integer(get_arg_value("--min-pathway-size", "10"))
max_set_size <- as.integer(get_arg_value("--max-pathway-size", "500"))

eligible_path <- file.path(analysis_dir, "eligible_cell_types.tsv")
if (!file.exists(eligible_path)) {
    stop("Missing eligible cell-type table: ", eligible_path)
}

eligible_dt <- fread(eligible_path)
if (nrow(eligible_dt) < 1L) {
    stop("eligible_cell_types.tsv is empty.")
}

# msigdbr has changed some column names across releases, so resolve them once
# here and keep the downstream Hallmark construction version-tolerant.
msig_all <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens"))
collection_col <- if ("gs_collection" %in% names(msig_all)) "gs_collection" else "gs_cat"
gene_col <- if ("gene_symbol" %in% names(msig_all)) "gene_symbol" else "human_gene_symbol"
set_id_col <- if ("gs_id" %in% names(msig_all)) "gs_id" else "gs_name"
hallmark_dt <- unique(msig_all[get(collection_col) == "H", .(
    source = "H",
    set_id = as.character(get(set_id_col)),
    set_name = as.character(gs_name),
    gene = toupper(as.character(get(gene_col)))
)])
hallmark_dt <- hallmark_dt[!is.na(gene) & gene != ""]

run_fgsea <- function(pathways, stats_vec, contrast_label) {
    if (length(pathways) < 1L) {
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

    fg_res <- fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, eps = 0)
    fg_dt <- as.data.table(fg_res)
    if (nrow(fg_dt) < 1L) {
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

summary_rows <- vector("list", nrow(eligible_dt))

for (i in seq_len(nrow(eligible_dt))) {
    cell_type <- eligible_dt$cell_type[[i]]
    cell_type_safe <- eligible_dt$cell_type_safe[[i]]
    counts_path <- file.path(handoff_dir, paste0(cell_type_safe, "_counts.tsv"))
    design_path <- file.path(handoff_dir, paste0(cell_type_safe, "_design.tsv"))

    if (!file.exists(counts_path) || !file.exists(design_path)) {
        stop("Missing handoff files for cell type: ", cell_type)
    }

    counts_dt <- fread(counts_path)
    design_dt <- fread(design_path)
    if (!"gene" %in% names(counts_dt)) {
        stop("Counts file must include a 'gene' column: ", counts_path)
    }

    count_cols <- setdiff(names(counts_dt), "gene")
    mat <- as.matrix(counts_dt[, ..count_cols])
    storage.mode(mat) <- "double"
    rownames(mat) <- toupper(counts_dt$gene)

    # Keep the design row order identical to the pseudobulk matrix columns so
    # sample labels and subtype assignments cannot drift during modeling.
    if (!all(design_dt$column_id %in% colnames(mat)) || !all(colnames(mat) %in% design_dt$column_id)) {
        stop("Counts columns and design column_id values do not match for ", cell_type)
    }
    ord <- match(colnames(mat), design_dt$column_id)
    if (anyNA(ord)) {
        stop("Failed to align design to counts matrix for ", cell_type)
    }
    design_dt <- design_dt[ord]
    if (!all(design_dt$column_id == colnames(mat))) {
        stop("Design rows and counts columns are misaligned for ", cell_type)
    }

    design_dt[, IntDx := factor(IntDx, levels = c("DLGNT_1", "DLGNT_2"))]
    if (anyNA(design_dt$IntDx)) {
        stop("Unexpected IntDx level in design for ", cell_type)
    }

    y <- edgeR::DGEList(counts = mat)
    # The tested-gene universe is defined here and must match the fgsea ranking
    # universe for this cell type.
    keep <- edgeR::filterByExpr(y, group = design_dt$IntDx)
    tested_genes <- rownames(mat)[keep]
    tested_dt <- data.table(cell_type = cell_type, gene = tested_genes)
    fwrite(tested_dt, file.path(results_dir, paste0(cell_type_safe, "_tested_genes.tsv")), sep = "\t")

    y <- y[keep, , keep.lib.sizes = FALSE]
    # TMM adjusts for library-size and composition differences across pseudobulk
    # samples before voom estimates the mean-variance relationship.
    y <- edgeR::calcNormFactors(y, method = "TMM")
    model_dt <- data.frame(IntDx = design_dt$IntDx)
    # The sample is the replicate unit here, so the design is just subtype with
    # DLGNT_1 as the baseline reference level.
    design <- model.matrix(~ IntDx, data = model_dt)
    coef_name <- "IntDxDLGNT_2"
    if (!coef_name %in% colnames(design)) {
        stop("Unable to find expected coefficient for ", cell_type)
    }

    v <- limma::voom(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")

    de_dt <- data.table(
        cell_type = cell_type,
        gene = toupper(rownames(tt)),
        logFC_dlgnt2_vs_dlgnt1 = as.numeric(tt$logFC),
        t_stat = as.numeric(tt$t),
        pvalue = as.numeric(tt$P.Value),
        fdr = as.numeric(tt$adj.P.Val)
    )
    fwrite(de_dt, file.path(results_dir, paste0(cell_type_safe, "_de.tsv")), sep = "\t")

    tested_gene_set <- unique(de_dt$gene)
    # GSEA should only see genes that were actually eligible for DE testing in
    # this cell type, otherwise pathway sizes and rankings become inconsistent.
    hallmark_keep <- hallmark_dt[gene %in% tested_gene_set]
    set_size_dt <- hallmark_keep[, .(set_size_in_tested_genes = uniqueN(gene)), by = .(source, set_id, set_name)]
    set_size_dt <- set_size_dt[
        set_size_in_tested_genes >= min_set_size &
            set_size_in_tested_genes <= max_set_size
    ]
    if (nrow(set_size_dt) > 0L) {
        hallmark_keep <- hallmark_keep[set_size_dt, on = .(source, set_id, set_name)]
        hallmark_keep[, pathway_key := paste(source, set_id, set_name, sep = "|||")]
        pathways <- split(hallmark_keep$gene, hallmark_keep$pathway_key)
        pathways <- lapply(pathways, unique)
    } else {
        pathways <- list()
    }

    stats_dlgnt2 <- de_dt$t_stat
    names(stats_dlgnt2) <- de_dt$gene
    stats_dlgnt2 <- sort(stats_dlgnt2, decreasing = TRUE)

    stats_dlgnt1 <- -de_dt$t_stat
    names(stats_dlgnt1) <- de_dt$gene
    stats_dlgnt1 <- sort(stats_dlgnt1, decreasing = TRUE)

    # Running both signed rankings makes the two directional subtype contrasts
    # explicit without refitting a second linear model.
    gsea_dt <- rbindlist(list(
        run_fgsea(pathways, stats_dlgnt2, "DLGNT_2 > DLGNT_1"),
        run_fgsea(pathways, stats_dlgnt1, "DLGNT_1 > DLGNT_2")
    ), fill = TRUE)
    if (nrow(gsea_dt) > 0L) {
        gsea_dt[, cell_type := cell_type]
    }
    fwrite(gsea_dt, file.path(results_dir, paste0(cell_type_safe, "_gsea.tsv")), sep = "\t")

    summary_rows[[i]] <- data.table(
        cell_type = cell_type,
        cell_type_safe = cell_type_safe,
        n_genes_input = nrow(mat),
        n_tested_genes = length(tested_gene_set),
        n_hallmark_pathways_retained = if (nrow(set_size_dt) > 0L) nrow(set_size_dt) else 0L
    )
}

# Join the model-side gene/pathway counts back onto the prepare-stage sample
# support summary so the notebook has one table per analyzed cell type.
summary_dt <- rbindlist(summary_rows, fill = TRUE)
summary_dt <- merge(
    eligible_dt[, .(
        cell_type,
        cell_type_safe,
        n_samples_dlgnt1,
        n_samples_dlgnt2,
        n_samples_total,
        total_cells_dlgnt1,
        total_cells_dlgnt2
    )],
    summary_dt,
    by = c("cell_type", "cell_type_safe"),
    all.y = TRUE,
    sort = TRUE
)
fwrite(summary_dt, file.path(results_dir, "cell_type_analysis_summary.tsv"), sep = "\t")

metadata_dt <- data.table(
    # data.table() reserves the name `key=` for table keys, so build this as
    # `meta_key` first and then rename it to keep the on-disk schema consistent
    # with prepare_metadata.tsv.
    meta_key = c(
        "analysis_dir",
        "handoff_dir",
        "results_dir",
        "min_pathway_size",
        "max_pathway_size",
        "pathway_source",
        "normalization",
        "analysis_note"
    ),
    value = c(
        analysis_dir,
        handoff_dir,
        results_dir,
        as.character(min_set_size),
        as.character(max_set_size),
        "Hallmark",
        "edgeR filterByExpr + TMM + limma voom/eBayes",
        "Exploratory subtype-only analysis; subtype and location are confounded."
    )
)
setnames(metadata_dt, "meta_key", "key")
fwrite(metadata_dt, file.path(results_dir, "results_metadata.tsv"), sep = "\t")

cat("Saved dlgnt12 limma/fgsea outputs to:", results_dir, "\n")
cat("Cell types analyzed:", nrow(summary_dt), "\n")
