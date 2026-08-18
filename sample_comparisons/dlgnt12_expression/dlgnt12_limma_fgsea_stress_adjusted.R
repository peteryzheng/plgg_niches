#!/usr/bin/env Rscript
# DLGNT_1 vs DLGNT_2 stress-adjusted pseudobulk DE + Hallmark fgsea.
# - Mirrors dlgnt12_limma_fgsea.R but adds a per-pseudobulk dissociation/IEG
#   stress score as a continuous covariate in the limma design.
# - Optionally drops user-specified samples (--drop-samples 267134,...).
# - Filters the IEG gene list from the fgsea ranking universe so the same genes
#   that defined the covariate cannot then drive the Hallmark enrichment.
# - Writes one DE table, tested-gene table, fgsea table, and stress-diagnostic
#   file per cell type, plus a global summary and metadata table.

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

has_flag <- function(flag) {
    any(args_all == flag)
}

analysis_dir <- normalizePath(get_arg_value("--indir", script_dir), mustWork = FALSE)
handoff_dir <- file.path(analysis_dir, "handoff")
results_dir <- normalizePath(
    get_arg_value("--outdir", file.path(analysis_dir, "results_stress_adjusted")),
    mustWork = FALSE
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

min_set_size <- as.integer(get_arg_value("--min-pathway-size", "10"))
max_set_size <- as.integer(get_arg_value("--max-pathway-size", "500"))

# --no-stress-covariate disables the covariate but keeps the rest of the script
# (drop-samples, IEG filtering of GSEA universe) so a baseline-equivalent run
# can be produced for direct comparison against the original results dir.
use_stress_covariate <- !has_flag("--no-stress-covariate")
# Drop IEG genes from the fgsea ranking by default; keep an opt-out for users
# who want to verify that filtering is what changed the GSEA result.
filter_iegs_from_gsea <- !has_flag("--keep-iegs-in-gsea")

drop_samples_arg <- get_arg_value("--drop-samples", "")
drop_samples <- character(0)
if (nzchar(drop_samples_arg)) {
    drop_samples <- trimws(strsplit(drop_samples_arg, ",", fixed = TRUE)[[1]])
    drop_samples <- drop_samples[nzchar(drop_samples)]
}

ieg_list_path <- get_arg_value("--ieg-list", file.path(analysis_dir, "dissociation_iegs.tsv"))
if (!file.exists(ieg_list_path)) {
    stop("Missing IEG/dissociation gene list: ", ieg_list_path)
}
ieg_dt <- fread(ieg_list_path)
if (!"gene" %in% names(ieg_dt) || nrow(ieg_dt) < 1L) {
    stop("IEG list must have a 'gene' column with at least one row: ", ieg_list_path)
}
ieg_genes <- unique(toupper(trimws(ieg_dt$gene)))
ieg_genes <- ieg_genes[nzchar(ieg_genes)]
if (length(ieg_genes) < 5L) {
    stop("IEG list has fewer than 5 valid genes: ", ieg_list_path)
}

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
diagnostic_rows <- vector("list", nrow(eligible_dt))

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

    # Apply --drop-samples before any modeling so the design and the count
    # matrix shrink in lockstep and downstream alignment checks stay valid.
    if (length(drop_samples) > 0L) {
        keep_design <- !(design_dt$sample %in% drop_samples)
        design_dt <- design_dt[keep_design]
        mat <- mat[, design_dt$column_id, drop = FALSE]
    }

    # Skip cell types where the drop list left fewer than 2 samples in either
    # subtype, since limma cannot estimate a contrast without replication.
    n_dlgnt1 <- sum(design_dt$IntDx == "DLGNT_1")
    n_dlgnt2 <- sum(design_dt$IntDx == "DLGNT_2")
    if (n_dlgnt1 < 2L || n_dlgnt2 < 2L) {
        message(sprintf(
            "Skipping %s after drop-samples: DLGNT_1=%d, DLGNT_2=%d (need >=2 each).",
            cell_type, n_dlgnt1, n_dlgnt2
        ))
        next
    }

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

    # Per-pseudobulk stress score: log10(IEG count fraction). Computed on the
    # raw count matrix (before filterByExpr) so the score reflects the full
    # pseudobulk transcriptome, not just the post-filter expressed-gene subset.
    iegs_in_mat <- intersect(ieg_genes, rownames(mat))
    if (length(iegs_in_mat) < 10L) {
        warning(sprintf(
            "Only %d/%d IEG genes intersect counts for %s; stress score may be unstable.",
            length(iegs_in_mat), length(ieg_genes), cell_type
        ))
    }
    if (length(iegs_in_mat) < 1L) {
        stop("No IEG genes intersect counts for ", cell_type, "; cannot compute stress score.")
    }
    ieg_sum <- colSums(mat[iegs_in_mat, , drop = FALSE])
    total_sum <- colSums(mat)
    # 1e-9 epsilon avoids log(0) in the rare case a pseudobulk has zero IEG
    # counts (would only happen with an unusually small handoff column).
    stress_score <- log10(ieg_sum / total_sum + 1e-9)
    design_dt[, stress_score := stress_score]

    # Stress-vs-IntDx diagnostic: how separable are the covariate and the
    # contrast in this cell type? Reported per-pseudobulk plus a one-row
    # cell-type summary so the user can spot near-collinear cases.
    intdx_numeric <- as.integer(design_dt$IntDx) - 1L
    stress_intdx_cor <- if (sd(stress_score) > 0 && sd(intdx_numeric) > 0) {
        suppressWarnings(stats::cor(stress_score, intdx_numeric))
    } else {
        NA_real_
    }
    stress_intdx_p <- tryCatch(
        stats::t.test(stress_score ~ design_dt$IntDx)$p.value,
        error = function(e) NA_real_
    )
    diagnostic_rows[[i]] <- data.table(
        cell_type = cell_type,
        cell_type_safe = cell_type_safe,
        column_id = design_dt$column_id,
        sample = design_dt$sample,
        IntDx = as.character(design_dt$IntDx),
        n_cells = design_dt$n_cells,
        ieg_count_sum = ieg_sum,
        total_count_sum = total_sum,
        ieg_fraction = ieg_sum / total_sum,
        stress_score = stress_score,
        n_iegs_intersected = length(iegs_in_mat),
        stress_intdx_pearson = stress_intdx_cor,
        stress_intdx_t_pvalue = stress_intdx_p,
        adjustment_unstable = !is.na(stress_intdx_cor) & abs(stress_intdx_cor) > 0.9
    )

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
    # Stress score enters as a continuous covariate ahead of IntDx so the
    # IntDxDLGNT_2 coefficient is interpreted as the subtype effect at fixed
    # per-pseudobulk stress. With --no-stress-covariate the model collapses to
    # the original ~ IntDx so this script can also reproduce the baseline.
    if (use_stress_covariate) {
        model_dt <- data.frame(
            stress_score = design_dt$stress_score,
            IntDx = design_dt$IntDx
        )
        design <- model.matrix(~ stress_score + IntDx, data = model_dt)
    } else {
        model_dt <- data.frame(IntDx = design_dt$IntDx)
        design <- model.matrix(~ IntDx, data = model_dt)
    }
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

    # Drop the IEG genes from the fgsea ranking universe so they cannot drive
    # Hallmark enrichment (their adjusted t is partially regressed-out by the
    # covariate, so leaving them in would mix circular and non-circular signal).
    gsea_genes <- unique(de_dt$gene)
    if (filter_iegs_from_gsea) {
        gsea_genes <- setdiff(gsea_genes, ieg_genes)
    }
    de_for_gsea <- de_dt[gene %in% gsea_genes]

    # GSEA should only see genes that were actually eligible for DE testing in
    # this cell type, otherwise pathway sizes and rankings become inconsistent.
    hallmark_keep <- hallmark_dt[gene %in% gsea_genes]
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

    stats_dlgnt2 <- de_for_gsea$t_stat
    names(stats_dlgnt2) <- de_for_gsea$gene
    stats_dlgnt2 <- sort(stats_dlgnt2, decreasing = TRUE)

    stats_dlgnt1 <- -de_for_gsea$t_stat
    names(stats_dlgnt1) <- de_for_gsea$gene
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
        n_samples_dlgnt1_used = n_dlgnt1,
        n_samples_dlgnt2_used = n_dlgnt2,
        n_genes_input = nrow(mat),
        n_tested_genes = length(tested_genes),
        n_genes_in_gsea_universe = length(gsea_genes),
        n_iegs_in_handoff = length(iegs_in_mat),
        n_iegs_filtered_from_gsea = if (filter_iegs_from_gsea) length(intersect(tested_genes, ieg_genes)) else 0L,
        stress_intdx_pearson = stress_intdx_cor,
        stress_intdx_t_pvalue = stress_intdx_p,
        adjustment_unstable = !is.na(stress_intdx_cor) & abs(stress_intdx_cor) > 0.9,
        n_hallmark_pathways_retained = if (nrow(set_size_dt) > 0L) nrow(set_size_dt) else 0L
    )
}

# Drop NULL entries from cell types that were skipped due to the drop-samples
# filter leaving fewer than two samples in a subtype.
summary_dt <- rbindlist(summary_rows[!vapply(summary_rows, is.null, logical(1))], fill = TRUE)
diagnostic_dt <- rbindlist(diagnostic_rows[!vapply(diagnostic_rows, is.null, logical(1))], fill = TRUE)

# Carry the prepare-stage sample-support columns onto the summary table so the
# notebook does not need to re-merge eligible_cell_types.tsv.
if (nrow(summary_dt) > 0L) {
    summary_dt <- merge(
        eligible_dt[, .(
            cell_type,
            cell_type_safe,
            n_samples_dlgnt1,
            n_samples_dlgnt2,
            total_cells_dlgnt1,
            total_cells_dlgnt2
        )],
        summary_dt,
        by = c("cell_type", "cell_type_safe"),
        all.y = TRUE,
        sort = TRUE
    )
}

fwrite(summary_dt, file.path(results_dir, "cell_type_analysis_summary.tsv"), sep = "\t")
fwrite(diagnostic_dt, file.path(results_dir, "stress_adjustment_diagnostics.tsv"), sep = "\t")

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
        "model_formula",
        "stress_covariate_used",
        "iegs_filtered_from_gsea",
        "ieg_list_path",
        "n_iegs_in_list",
        "drop_samples",
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
        if (use_stress_covariate) "~ stress_score + IntDx" else "~ IntDx",
        as.character(use_stress_covariate),
        as.character(filter_iegs_from_gsea),
        ieg_list_path,
        as.character(length(ieg_genes)),
        if (length(drop_samples) > 0L) paste(drop_samples, collapse = ",") else "",
        "Stress-adjusted exploratory subtype analysis; subtype and location remain confounded."
    )
)
setnames(metadata_dt, "meta_key", "key")
fwrite(metadata_dt, file.path(results_dir, "results_metadata.tsv"), sep = "\t")

cat("Saved stress-adjusted dlgnt12 outputs to:", results_dir, "\n")
cat("Cell types analyzed:", nrow(summary_dt), "\n")
if (length(drop_samples) > 0L) {
    cat("Samples dropped:", paste(drop_samples, collapse = ", "), "\n")
}
cat("Stress covariate used:", use_stress_covariate, "\n")
cat("IEGs filtered from GSEA universe:", filter_iegs_from_gsea, "\n")
