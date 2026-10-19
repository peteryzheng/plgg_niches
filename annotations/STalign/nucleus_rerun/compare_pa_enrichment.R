# Compare the pathology-enrichment results built on the old 030226 labels with the reruns on the nucleus labels.
#
# Old side: NEW_DIR/qc/pa_enrichment_old_labels/ (the outputs copied on 2026-10-19 before the rerun; MANIFEST.md5).
# New side: the same files at their canonical places in the main checkout (the notebooks hardcode those paths).
# Analyses (each skipped with a message when its new output is not there yet):
#   A  annotations/linear_models.qmd          beta-binomial, cell types and niches vs region (all histologies)
#   B  annotations/linear_models_niches.qmd   beta-binomial, niches vs region (all histologies; deck PA panel)
#   C  niche_interpretation_proseg.qmd        top_pathology column of niche_dictionary.tsv (reads B)
#   D  PA morphology state                    limma + fgsea, loose vs compact (the old side is the 2026-09-20
#                                             run: same 091126 export and log2 fix, so only the labels differ)
#   E  measured_panel_qc                      measured vs imputed attenuation (PA rows)
# A call is the analysis' own: A/B sign_call (q < 0.25); D tiers q < 0.05 / q < 0.25 with the sign of NES.
# Output: NEW_DIR/qc/pa_enrichment_old_vs_new.tsv (one row per test; D only the themed glycan/SASP sets) and
#         NEW_DIR/qc/pa_enrichment_old_vs_new_summary.tsv (one row per analysis and group).
# Usage : Rscript annotations/STalign/nucleus_rerun/compare_pa_enrichment.R   (conda env spatial)
suppressMessages(library(data.table))

data_dir <- "/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data"
code_dir <- "/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code"
qc_dir <- file.path(data_dir, "Xenium_annotations/pathology_annotations_101326_nucleus/qc")
old_dir <- file.path(qc_dir, "pa_enrichment_old_labels")
ms_dir <- file.path(code_dir, "sc_integration/ENVI/2_spatial_imputation/PA/morphology_state_gsea")

# new outputs count only if written after the snapshot (otherwise they are still the old files)
snapshot_time <- file.mtime(file.path(old_dir, "MANIFEST.md5"))
new_file <- function(path) {
    if (!file.exists(path) || file.mtime(path) <= snapshot_time) {
        message("skip (no new output yet): ", path)
        return(NULL)
    }
    path
}

tests <- list()
summ <- list()

# ---- A and B: one row per (model file, analysis, feature, region) ----
mm_old <- file.path(old_dir, "A_B_mixed_models")
mm_new <- file.path(code_dir, "annotations/annotations/model_outputs/mixed_models")
for (stem in c("linear_models_cell_type", "linear_models_niche_type", "linear_models_niches")) {
    f_new <- new_file(file.path(mm_new, paste0(stem, "_primary_results.tsv")))
    if (is.null(f_new)) next
    cols <- c("analysis_id", "feature", "annotation_region", "estimate_log_odds", "q_value", "sign_call")
    o <- fread(file.path(mm_old, paste0(stem, "_primary_results.tsv")), select = cols, colClasses = list(character = "feature"))
    n <- fread(f_new, select = cols, colClasses = list(character = "feature"))
    m <- merge(o, n, by = c("analysis_id", "feature", "annotation_region"), all = TRUE, suffixes = c("_old", "_new"))
    tests[[stem]] <- m[, .(analysis = stem, group = analysis_id, feature, region = annotation_region,
                           old_estimate = estimate_log_odds_old, new_estimate = estimate_log_odds_new,
                           old_q = q_value_old, new_q = q_value_new, old_call = sign_call_old, new_call = sign_call_new)]
}

# ---- C: the niche dictionary's top pathology region ----
f_new <- new_file(file.path(code_dir, "niche/proseg_output_analysis/niche_dictionary.tsv"))
if (!is.null(f_new)) {
    o <- fread(file.path(old_dir, "C_niche_interpretation/niche_dictionary.tsv"), select = c("niche_type", "top_pathology"))
    n <- fread(f_new, select = c("niche_type", "top_pathology"))
    m <- merge(o, n, by = "niche_type", all = TRUE, suffixes = c("_old", "_new"))
    tests[["niche_dictionary"]] <- m[, .(analysis = "niche_dictionary", group = "all", feature = as.character(niche_type),
                                         region = "top_pathology", old_call = top_pathology_old, new_call = top_pathology_new)]
}

# ---- D: fgsea (themed sets as rows, all sets in the summary) and limma ----
tier <- function(nes, padj) fifelse(is.na(padj), NA_character_,
                                    fifelse(padj < 0.05, paste0(fifelse(nes > 0, "up", "down"), " q<0.05"),
                                            fifelse(padj < 0.25, paste0(fifelse(nes > 0, "up", "down"), " q<0.25"), "ns")))
f_new <- new_file(file.path(ms_dir, "morphology_state_gsea_results.csv"))
if (!is.null(f_new)) {
    cols <- c("cell_type", "contrast", "set_name", "theme", "NES", "padj")
    o <- fread(file.path(old_dir, "D_morphology_state/morphology_state_gsea_results.csv"), select = cols)
    n <- fread(f_new, select = cols)
    # fgsea runs both directions; keep one contrast so each set is counted once
    o <- o[contrast == "loose > compact"]
    n <- n[contrast == "loose > compact"]
    m <- merge(o, n, by = c("cell_type", "set_name", "theme"), all = TRUE, suffixes = c("_old", "_new"))
    m[, `:=`(old_call = tier(NES_old, padj_old), new_call = tier(NES_new, padj_new))]
    tests[["morphology_gsea"]] <- m[theme != "other", .(analysis = "morphology_gsea", group = cell_type,
                                                        feature = set_name, region = paste("loose vs compact;", theme),
                                                        old_estimate = NES_old, new_estimate = NES_new,
                                                        old_q = padj_old, new_q = padj_new, old_call, new_call)]
    summ[["morphology_gsea"]] <- m[, .(analysis = "morphology_gsea_all_sets", n_tests = .N,
                                       n_call_changed = sum(old_call != new_call, na.rm = TRUE),
                                       n_sig_old = sum(padj_old < 0.25, na.rm = TRUE), n_sig_new = sum(padj_new < 0.25, na.rm = TRUE),
                                       estimate_r = cor(NES_old, NES_new, use = "complete.obs")), by = .(group = cell_type)]
}

# the 34-gene curated glycan-ECM panel (copied from 4_PA_morphology_state.qmd, GLYCAN_ECM_GENES)
glycan_genes <- c("HAS1", "HAS2", "HAS3", "HAPLN1", "HAPLN2", "VCAN", "BCAN", "ACAN", "NCAN", "LUM",
                  "SPOCK1", "SPOCK2", "SPOCK3", "BGN", "DCN", "FMOD", "PRELP", "SDC1", "SDC2", "GPC3", "TNC", "FN1", "SPP1",
                  "CHSY1", "CHPF", "CSGALNACT1", "XYLT1", "B3GNT2", "B4GALT1", "MGAT5", "ST6GAL1", "ST3GAL1", "FUT8", "GALNT1")
f_new <- new_file(file.path(ms_dir, "morphology_state_limma_de.csv"))
if (!is.null(f_new)) {
    o <- fread(file.path(old_dir, "D_morphology_state/morphology_state_limma_de.csv"))
    n <- fread(f_new)
    m <- merge(o, n, by = c("cell_type", "gene"), suffixes = c("_old", "_new"))
    summ[["morphology_limma"]] <- m[, .(analysis = "morphology_limma", n_tests = .N,
                                        n_call_changed = sum((fdr_old < 0.05) != (fdr_new < 0.05)),
                                        n_sig_old = sum(fdr_old < 0.05), n_sig_new = sum(fdr_new < 0.05),
                                        estimate_r = cor(log2FC_loose_vs_compact_old, log2FC_loose_vs_compact_new),
                                        glycan_mean_log2FC_old = mean(log2FC_loose_vs_compact_old[gene %in% glycan_genes]),
                                        glycan_mean_log2FC_new = mean(log2FC_loose_vs_compact_new[gene %in% glycan_genes]),
                                        glycan_n_up_old = sum(log2FC_loose_vs_compact_old[gene %in% glycan_genes] > 0),
                                        glycan_n_up_new = sum(log2FC_loose_vs_compact_new[gene %in% glycan_genes] > 0),
                                        glycan_n_genes = sum(gene %in% glycan_genes)), by = .(group = cell_type)]
}

# ---- E: measured vs imputed attenuation, PA rows ----
f_new <- new_file(file.path(code_dir, "sc_integration/ENVI/2_spatial_imputation/measured_panel_qc/measured_vs_imputed_attenuation.tsv"))
if (!is.null(f_new)) {
    # integer count columns become doubles so each table melts into one numeric value column
    read_att <- function(f) { d <- fread(f)[cohort == "PA"]; num <- setdiff(names(d), c("cohort", "lineage"))
                              d[, (num) := lapply(.SD, as.numeric), .SDcols = num] }
    o <- read_att(file.path(old_dir, "E_measured_panel_qc/measured_vs_imputed_attenuation.tsv"))
    n <- read_att(f_new)
    m <- merge(melt(o, id.vars = c("cohort", "lineage"), variable.factor = FALSE),
               melt(n, id.vars = c("cohort", "lineage"), variable.factor = FALSE),
               by = c("cohort", "lineage", "variable"), suffixes = c("_old", "_new"))
    tests[["measured_attenuation"]] <- m[, .(analysis = "measured_attenuation", group = lineage, feature = variable,
                                             region = "loose vs compact", old_estimate = value_old, new_estimate = value_new)]
}

# ---- write ----
res <- rbindlist(tests, fill = TRUE)
if (nrow(res) == 0) stop("no new outputs yet")
res[, call_changed := !is.na(old_call) & !is.na(new_call) & old_call != new_call]
fwrite(res, file.path(qc_dir, "pa_enrichment_old_vs_new.tsv"), sep = "\t")

# A/B summarised per model file and histology, alongside the D summaries
ab <- res[analysis %like% "^linear_models"]
ab <- if (nrow(ab) == 0) NULL else ab[, .(n_tests = .N, n_call_changed = sum(call_changed),
                                          n_sig_old = sum(old_call != "not_significant", na.rm = TRUE),
                                          n_sig_new = sum(new_call != "not_significant", na.rm = TRUE),
                                          estimate_r = cor(old_estimate, new_estimate, use = "complete.obs")),
                                      by = .(analysis, group)]
summ_dt <- rbindlist(c(list(ab), summ), fill = TRUE)
fwrite(summ_dt, file.path(qc_dir, "pa_enrichment_old_vs_new_summary.tsv"), sep = "\t")
print(summ_dt)
cat("\ncalls that changed:\n")
print(res[call_changed == TRUE, .(analysis, group, feature, region, old_call, new_call,
                                  old_estimate = signif(old_estimate, 3), new_estimate = signif(new_estimate, 3),
                                  old_q = signif(old_q, 2), new_q = signif(new_q, 2))])
