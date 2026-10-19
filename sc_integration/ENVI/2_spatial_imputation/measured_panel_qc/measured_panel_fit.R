#!/usr/bin/env Rscript
# Measured-panel limma fit -- ground-truth arm of the imputation QC.
#
# Runs the IDENTICAL statistical treatment as morphology_state.R (expression filter ->
# log2 -> lmFit -> eBayes(trend=TRUE) on ~ sample + region) so that any difference between
# this and the imputed result is attributable to imputation, not to method.
#
# SCALE NOTE: the expression floor cannot be a shared constant across the two arms.
# measured X is library-size normalized (per-cell total 65.4 over 266 genes, so a typical
# gene mean is ~0.25), while the imputed matrix sits near 3e-4. A fixed absolute floor
# would therefore mean completely different things in the two arms. We use a QUANTILE
# floor here so the filter removes a comparable share of each arm's own distribution, and
# the downstream comparison is additionally restricted to genes passing in BOTH arms.
suppressPackageStartupMessages({library(data.table); library(limma)})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- dirname(normalizePath(sub("^--file=", "",
            commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))

EXPR_QUANTILE <- 0.10   # drop the bottom decile of expression within each cell type
LOG_EPS       <- 1e-6

COH <- list(CL = c("nonlipid", "lipid"), PA = c("compact", "loose"))
if (length(args) > 0) COH <- COH[args]   # optional cohort names, as in measured_panel_prepare.py

for (coh in names(COH)) {
    ref <- COH[[coh]][1]; alt <- COH[[coh]][2]
    pb  <- fread(file.path(out_dir, sprintf("measured_%s_pseudobulk_matrix.csv", coh)))
    des <- fread(file.path(out_dir, sprintf("measured_%s_sample_design.csv", coh)))
    genes <- toupper(pb$gene)
    mat <- as.matrix(pb[, setdiff(names(pb), "gene"), with = FALSE]); rownames(mat) <- genes

    de_list <- list()
    for (ct in unique(des$cell_type)) {
        d <- des[cell_type == ct]
        if (uniqueN(d$sample) < 2L) next
        M <- mat[, d$column_id, drop = FALSE]
        keep <- rowMeans(M) >= quantile(rowMeans(M), EXPR_QUANTILE, na.rm = TRUE)
        M <- log2(M[keep, , drop = FALSE] + LOG_EPS)

        sf <- factor(d$sample)
        rf <- factor(d$region, levels = c(ref, alt))
        fit <- eBayes(lmFit(M, model.matrix(~ sf + rf)), trend = TRUE)
        tt  <- topTable(fit, coef = paste0("rf", alt), number = Inf, sort.by = "none")
        de_list[[ct]] <- data.table(cell_type = ct, gene = toupper(rownames(tt)),
                                    log2FC = as.numeric(tt$logFC), t_stat = as.numeric(tt$t),
                                    pvalue = as.numeric(tt$P.Value), fdr = as.numeric(tt$adj.P.Val))
        cat(sprintf("  %s %-14s %d genes kept, %d columns\n", coh, ct, nrow(M), ncol(M)))
    }
    de <- rbindlist(de_list)
    fwrite(de, file.path(out_dir, sprintf("measured_%s_limma_de.csv", coh)))
    cat(sprintf("%s: %d rows | genes FDR<0.25: %d | min p %.2e | min FDR %.4f\n\n",
                coh, nrow(de), sum(de$fdr < 0.25), min(de$pvalue), min(de$fdr)))
}
