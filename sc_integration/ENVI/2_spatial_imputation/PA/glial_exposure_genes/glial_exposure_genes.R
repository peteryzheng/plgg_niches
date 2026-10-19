#!/usr/bin/env Rscript
# Gene-level myeloid glial-adjacency contrast (R aggregation stage) -- pilocytic astrocytoma.
#
# Reads the per-sample, per-gene effect sizes from glial_exposure_genes_prepare.py and
# collapses them to one cohort row per gene, then writes the collaborator handoff.
#
# INFERENTIAL FRAMING (this drives every choice below). The per-cell p-values from the
# Python stage are pseudoreplicated -- thousands of cells per group collapse the SEs, so
# nearly everything would clear q<0.05 within a sample. The real inferential unit is the
# SAMPLE (n = 4 PA). We therefore rank by cross-sample directional consistency plus mean
# standardized effect, exactly as the parent notebook does, and report the sample-level
# t-test q as a descriptive column only. We deliberately do NOT select on FDR: a
# one-sample t-test on 4 betas has 3 df, so clearing BH across ~18,000 genes needs
# p < 3e-6, which no gene can reach -- an FDR-selected table would simply be empty.
# The deliverable is an explicitly RANKED, HYPOTHESIS-GENERATING list.
# (see feedback_percell_pseudoreplication, project_cohort_granularity)
#
# Inputs  (from the Python stage, same directory)
#   - <tag>_per_sample.tsv, <tag>_sample_diagnostics.tsv
#   - ../myeloid_glial_exposure_nested_results.tsv  (to re-derive the 20 signatures the
#     notebook found all-4-consistent; not hardcoded)
#   - ../../myeloid_signatures.csv                  (56 signature memberships)
#
# Outputs
#   - <tag>_cohort_results.tsv          every gene, cohort-level stats + annotation
#   - <tag>_gsea_crosscheck.tsv         preranked fgsea of the 56 signatures on the
#                                       gene-level statistic (does the notebook's
#                                       signature story reproduce from the genes?)
#   - <tag>_spillover_diagnostics.tsv   glial-share vs effect-size relationship
#   - handoff/myeloid_glial_adjacency_gene_list.tsv   the curated ranked list
#   - handoff/glial_exposure_gene_results_full.tsv    full appendix table

suppressPackageStartupMessages({
    library(data.table)
})
if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for the signature cross-check.")
}

set.seed(5555)

# Resolve paths from the script location so handoff files land next to the pipeline.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
if (length(file_arg) < 1) stop("Unable to determine script path from --file argument.")
out_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
pa_dir <- dirname(out_dir)
imp_dir <- dirname(pa_dir)
handoff_dir <- file.path(out_dir, "handoff")
dir.create(handoff_dir, recursive = TRUE, showWarnings = FALSE)

# Output prefix; must match --out-tag passed to the Python stage.
tag <- Sys.getenv("OUT_TAG", unset = "glial_exposure_genes")

per_sample_path <- file.path(out_dir, paste0(tag, "_per_sample.tsv"))
if (!file.exists(per_sample_path)) stop("Missing Python handoff: ", per_sample_path)
per_sample <- fread(per_sample_path, colClasses = list(character = "sample_id"))

n_samples_total <- per_sample[source == "imputed", uniqueN(sample_id)]
cat("Per-sample rows:", nrow(per_sample), "| samples:", n_samples_total, "\n")

# ---- Cohort collapse: directional consistency + mean standardized effect ------------
# `consistent_all` requires the gene to be fittable in EVERY sample and all per-sample
# betas to share a sign. Calibration for the reader: with 4 samples a 4/4 sign run is
# p = 0.125 under a sign test, so ~1/8 of genes reach it by chance -- consistency is a
# ranking device here, not a significance claim.
collapse_cohort <- function(dt) {
    dt[is.finite(beta_m1), .(
        n_samples      = .N,
        n_pos          = sum(beta_m1 > 0),
        n_neg          = sum(beta_m1 < 0),
        mean_beta_m1   = mean(beta_m1),
        min_beta_m1    = min(beta_m1),
        max_beta_m1    = max(beta_m1),
        mean_beta_m2   = mean(beta_m2[is.finite(beta_m2)]),
        # sample-level test on the per-sample betas (descriptive; see header note)
        t_p            = if (.N >= 2) tryCatch(t.test(beta_m1)$p.value,
                                               error = function(e) NA_real_) else NA_real_,
        mean_expr_myeloid = mean(mean_expr_myeloid),
        min_sd_myeloid    = min(sd_myeloid),
        mean_glial_share  = mean(glial_share[is.finite(glial_share)]),
        min_n_other       = min(n_other),
        min_n_pa          = min(n_pa)
    ), by = .(gene, source)]
}
cohort <- collapse_cohort(per_sample)
cohort[, consistent_all := n_samples == n_samples_total & (n_pos == n_samples_total |
                                                           n_neg == n_samples_total)]
cohort[, direction := fifelse(mean_beta_m1 > 0, "up_near_pa_niche_glia",
                              "down_near_pa_niche_glia")]
# % shrinkage of the effect once the neighbouring glia's own MAPK dose is adjusted for;
# small attenuation = the effect is not merely a MAPK proxy.
cohort[, attenuation_pct_m2 := round(100 * (1 - mean_beta_m2 / mean_beta_m1), 1)]
# BH within each source (imputed / panel are separate testing families).
cohort[, t_q := p.adjust(t_p, method = "BH"), by = source]

imputed <- cohort[source == "imputed"]
panel <- cohort[source == "panel"]
cat(sprintf("Imputed genes: %d (%d all-%d-consistent) | panel genes: %d\n",
            nrow(imputed), imputed[consistent_all == TRUE, .N], n_samples_total,
            nrow(panel)))

# ---- On-panel validation arm ---------------------------------------------------------
# The 266 measured Xenium genes are immune to imputation artifacts (though not to
# segmentation spillover). Agreement here is the strongest credibility evidence we have.
panel_ann <- panel[, .(gene,
                       panel_beta_m1 = mean_beta_m1,
                       panel_consistent = consistent_all,
                       panel_direction = direction)]
res <- merge(imputed, panel_ann, by = "gene", all.x = TRUE)
res[, on_xenium_panel := !is.na(panel_beta_m1)]
res[, panel_agrees := on_xenium_panel & (sign(panel_beta_m1) == sign(mean_beta_m1))]

# ---- Signature membership annotation -------------------------------------------------
# Re-derive the notebook's all-4-consistent signature set programmatically instead of
# hardcoding it, so this stays correct if the parent analysis is re-run.
nested <- fread(file.path(pa_dir, "myeloid_glial_exposure_nested_results.tsv"),
                colClasses = list(character = "sample_id"))
sig_cons <- nested[!is.na(Estimate), .(n = .N, n_pos = sum(Estimate > 0),
                                       n_neg = sum(Estimate < 0)), by = signature]
n_sig_samples <- nested[, uniqueN(sample_id)]
consistent_sigs <- sub("_score$", "",
                       sig_cons[n == n_sig_samples & (n_pos == n | n_neg == n), signature])
cat("Notebook all-consistent signatures:", length(consistent_sigs), "\n")

# The 6 literature myeloid/immune sets (everything that is not a Hallmark set).
sig_members <- fread(file.path(imp_dir, "myeloid_signatures.csv"))
sig_members[, Gene := toupper(Gene)]
literature_sigs <- setdiff(unique(sig_members$Signature),
                           grep("^HALLMARK_", unique(sig_members$Signature), value = TRUE))

annotate_sets <- function(genes, sets) {
    sub_dt <- sig_members[Signature %in% sets & Gene %in% genes]
    if (nrow(sub_dt) == 0L) return(rep(NA_character_, length(genes)))
    lookup <- sub_dt[, .(sigs = paste(sort(unique(Signature)), collapse = ";")), by = Gene]
    lookup[match(genes, Gene), sigs]
}
res[, consistent_signatures := annotate_sets(gene, consistent_sigs)]
res[, literature_myeloid_sets := annotate_sets(gene, literature_sigs)]

# ---- Spillover control ---------------------------------------------------------------
# ENVI imputation conditions on spatial context, so a myeloid cell ringed by glia is
# predicted to look more glial. A signature dilutes such genes across 200 members; a
# RANKED GENE LIST actively selects for them. glial_share > 0.5 means the gene is
# expressed more in glia than in myeloid cells, i.e. a spillover-prone gene.
spill_rows <- list()
add_spill <- function(label, dt) {
    if (nrow(dt) < 3L) return(invisible(NULL))
    ct <- suppressWarnings(cor.test(dt$mean_beta_m1, dt$mean_glial_share,
                                    method = "spearman", exact = FALSE))
    spill_rows[[label]] <<- data.table(
        subset = label, n_genes = nrow(dt),
        spearman_rho = unname(ct$estimate), p_value = ct$p.value,
        median_glial_share_up = median(dt[mean_beta_m1 > 0, mean_glial_share], na.rm = TRUE),
        median_glial_share_down = median(dt[mean_beta_m1 < 0, mean_glial_share], na.rm = TRUE),
        frac_up_glial_leaning = dt[mean_beta_m1 > 0, mean(mean_glial_share > 0.5, na.rm = TRUE)]
    )
}
add_spill("all_imputed_genes", res)
add_spill("consistent_only", res[consistent_all == TRUE])
add_spill("top100_abs_effect", res[order(-abs(mean_beta_m1))][1:min(100, .N)])
spillover <- rbindlist(spill_rows, fill = TRUE)
cat("\nSpillover diagnostics:\n"); print(spillover)

# ---- Curated collaborator list -------------------------------------------------------
# Three filters beyond consistency, all for the same reason: beta is STANDARDIZED, so a
# gene with near-zero imputed expression or near-zero within-sample variance gets an
# inflated |beta| and would otherwise dominate the ranking on noise alone.
expr_floor <- quantile(res$mean_expr_myeloid, 0.10, na.rm = TRUE)
sd_floor <- quantile(res$min_sd_myeloid, 0.10, na.rm = TRUE)
MIN_ABS_BETA <- 0.1
cat(sprintf(paste0("\nFilters: mean_expr_myeloid >= %.3g (10th pct), ",
                   "min_sd_myeloid >= %.3g (10th pct), |mean_beta| >= %.2f\n"),
            expr_floor, sd_floor, MIN_ABS_BETA))

curated <- res[consistent_all == TRUE &
               mean_expr_myeloid >= expr_floor &
               min_sd_myeloid >= sd_floor &
               abs(mean_beta_m1) >= MIN_ABS_BETA]
setorder(curated, -consistent_all, direction)
curated <- curated[order(-abs(mean_beta_m1))]
curated[, rank := .I]
cat("Curated list:", nrow(curated), "genes (",
    curated[direction == "up_near_pa_niche_glia", .N], "up /",
    curated[direction == "down_near_pa_niche_glia", .N], "down )\n")

# sanity read-out: are the top-ranked genes near-constant (the failure mode the sd floor
# guards against)? reported so it can be judged rather than assumed away.
cat("Top-50 min_sd_myeloid quantiles:",
    paste(signif(quantile(curated[1:min(50, .N), min_sd_myeloid], c(0, .5, 1)), 3),
          collapse = " / "), "\n")

# ---- Signature cross-check: preranked GSEA on the gene-level statistic ---------------
# If the gene-level contrast is sound, the signatures the notebook found consistent
# should re-emerge in the SAME direction from a ranking that never saw them. Uses the
# same 56 sets as the notebook, restricted to the tested gene universe.
universe <- res$gene
sets_dt <- sig_members[Gene %in% universe]
set_sizes <- sets_dt[, .N, by = Signature][N >= 10L]
sets_dt <- sets_dt[Signature %in% set_sizes$Signature]
pathways <- split(sets_dt$Gene, sets_dt$Signature)
pathways <- lapply(pathways, unique)

stats_vec <- res$mean_beta_m1
names(stats_vec) <- res$gene
stats_vec <- sort(stats_vec[is.finite(stats_vec)], decreasing = TRUE)

fg <- as.data.table(fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, eps = 0))
fg[, leading_edge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
fg[, notebook_consistent := pathway %in% consistent_sigs]
# direction the notebook reported for each signature, to check agreement
nb_dir <- nested[!is.na(Estimate), .(nb_mean_beta = mean(Estimate)), by = signature]
nb_dir[, pathway := sub("_score$", "", signature)]
fg <- merge(fg, nb_dir[, .(pathway, nb_mean_beta)], by = "pathway", all.x = TRUE)
fg[, direction_agrees := sign(NES) == sign(nb_mean_beta)]
fg[, abs_nes := abs(NES)]
setorder(fg, padj, -abs_nes)

agree <- fg[notebook_consistent == TRUE]
cat(sprintf(paste0("\nGSEA cross-check: %d/%d notebook-consistent signatures recovered ",
                   "with the SAME direction (%d at padj<0.05).\n"),
            agree[direction_agrees == TRUE, .N], nrow(agree),
            agree[direction_agrees == TRUE & padj < 0.05, .N]))
print(agree[, .(pathway, size, NES, padj, nb_mean_beta, direction_agrees)])

# ---- Write outputs -------------------------------------------------------------------
full_cols <- c("gene", "direction", "mean_beta_m1", "mean_beta_m2", "attenuation_pct_m2",
               "n_samples", "n_pos", "n_neg", "consistent_all", "min_beta_m1", "max_beta_m1",
               "t_p", "t_q", "mean_expr_myeloid", "min_sd_myeloid", "mean_glial_share",
               "on_xenium_panel", "panel_beta_m1", "panel_agrees",
               "consistent_signatures", "literature_myeloid_sets",
               "min_n_pa", "min_n_other")
full_out <- res[order(-consistent_all, -abs(mean_beta_m1)), ..full_cols]

fwrite(full_out, file.path(out_dir, paste0(tag, "_cohort_results.tsv")), sep = "\t")
fwrite(fg[, .(pathway, size, NES, pval, padj, notebook_consistent, nb_mean_beta,
              direction_agrees, leading_edge)],
       file.path(out_dir, paste0(tag, "_gsea_crosscheck.tsv")), sep = "\t")
fwrite(spillover, file.path(out_dir, paste0(tag, "_spillover_diagnostics.tsv")), sep = "\t")

fwrite(full_out, file.path(handoff_dir, "glial_exposure_gene_results_full.tsv"), sep = "\t")
fwrite(curated[, c("rank", full_cols), with = FALSE],
       file.path(handoff_dir, "myeloid_glial_adjacency_gene_list.tsv"), sep = "\t")

cat("\nWrote:\n  ", file.path(out_dir, paste0(tag, "_cohort_results.tsv")),
    "\n  ", file.path(handoff_dir, "myeloid_glial_adjacency_gene_list.tsv"),
    "\n  ", file.path(handoff_dir, "glial_exposure_gene_results_full.tsv"), "\n")
