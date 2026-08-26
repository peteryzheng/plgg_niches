# Builds a full per-cluster summary table (all lam0.2 cell-type clusters,
# across every resolution tested so far -- 0.2, 0.3, 0.4, 0.5, 1) so
# different cell-type-calling heuristics can be tried quickly against saved
# numbers instead of re-deriving them by hand per cluster. For each cluster:
# size, sample-mixing entropy, the full per-lineage marker_panel score vector
# (not just the winning call) with the margin over the runner-up, and the
# top unbiased (full-transcriptome) DE genes already computed by
# cell_type_marker_ident(). No new clustering -- everything here is derived
# from already-saved objects.
#
# Covers three separate annotated objects (the original res 0.5/1 run, plus
# the two later low-resolution reclustering jobs at res 0.3 and res 0.2/0.4)
# -- add a new list entry below whenever another resolution is tested.

# Path resolution per AGENTS.md so this script works across mac/cluster mounts.
home <- Sys.getenv("HOME")
if (home %in% c("/Users/youyun", "/Users/youyunzheng")) {
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
} else if (home == "/PHShome/yz762") {
    workdir <- "/data/beroukhim1/"
} else if (home == "/home/yz762") {
    workdir <- "/mnt/storage/dept/medonc/beroukhim/"
} else {
    workdir <- "/xchip/beroukhimlab/"
}

source(paste0(workdir, "youyun/plgg/code/helpers/spatial_helper.R"))

DATA_ROOT = paste0(workdir, "youyun/plgg/data")
outputdir = paste0(DATA_ROOT, "/banksy_param_search/k1_15_k2_30_lambda1_0.2_lambda2_0.8_npcs_20_kct_50_resct_0.5,1_kni_50_resni_0.5,1")

# One entry per annotated object -- `ts` must match the timestamp embedded in
# both the object's filename and its cell_type_marker_ident() marker files.
runs = list(
    list(ts = "20260818_005830",
         rds = paste0(outputdir, "/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.5_1_kni_50_resni_0.5_1_20260818_005830.rds")),
    list(ts = "20260824_225621",
         rds = paste0(outputdir, "/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.3_kni_50_resni_0.5_1_20260824_225621.rds")),
    list(ts = "20260824_042548",
         rds = paste0(outputdir, "/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_0.4_kni_50_resni_0.5_1_20260824_042548.rds"))
)

lineage_names = names(marker_panel)
all_rows = list()
score_mats = list()

for (run in runs) {
    print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loading ', run$rds, ' ...'))
    banksy_spe = readRDS(run$rds)
    print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loaded. ncol = ', ncol(banksy_spe)))

    n_samples = length(unique(colData(banksy_spe)$sample_id))
    expr = assay(banksy_spe, "normcounts")
    cell_type_clusters = grep('lam0.2', clusterNames(banksy_spe), value = TRUE)

    for (cc in cell_type_clusters) {
        res_label = sub(".*_res", "res", cc)
        print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Scoring ', cc, ' ...'))

        cl_vec = colData(banksy_spe)[[cc]]
        cl = as.character(cl_vec)

        # Full per-lineage score matrix (mean(auc-0.5) per cluster x lineage) --
        # the argmax across this row is annotate_cell_types()'s actual call.
        sig_wide = score_lineages(expr, cl_vec, return_scores = TRUE)
        score_mats[[res_label]] = sig_wide

        m = as.matrix(sig_wide[, ..lineage_names])
        ord = t(apply(m, 1, function(x) sort(x, decreasing = TRUE)))
        top_lineage = lineage_names[max.col(replace(m, is.na(m), -Inf), ties.method = 'first')]
        margin = ord[, 1] - ord[, 2]

        # Cluster size + sample-mixing entropy (reuse entropy_norm(), same calc
        # used by the review panels).
        by_clust = split(colData(banksy_spe)$sample_id, cl)
        entropy = vapply(by_clust[as.character(sig_wide$cluster)], entropy_norm, numeric(1), n_samples = n_samples)
        n_cells = vapply(by_clust[as.character(sig_wide$cluster)], length, integer(1))

        # Top 5 unbiased DE genes per cluster from the matching full-transcriptome
        # findMarkers() output, ranked by summary.AUC (one-vs-rest specificity)
        # rather than the DataFrame's default row order (scran's `Top` column --
        # best rank in ANY single pairwise comparison). `Top` turned out to be
        # dominated by a handful of genes (e.g. B4GALNT1, COL25A1) recurring
        # across nearly every cluster, a known artifact of that statistic rather
        # than genuine broad markers; summary.AUC gives a cleaner, more
        # comparable cross-check against the marker_panel score above.
        marker_file = paste0(outputdir, "/banksy_clusters_markers_k_geom_15_30_pc_20_", gsub("clust_Harmony_BANKSY_", "", cc), "_", run$ts, ".rds")
        mk = readRDS(marker_file)
        top_de = vapply(as.character(sig_wide$cluster), function(id) {
            if (!id %in% names(mk)) return(NA_character_)
            d = as.data.frame(mk[[id]])
            paste(rownames(d)[order(-d$summary.AUC)][1:5], collapse = ", ")
        }, character(1))

        row_dt = data.table(
            resolution = res_label, cluster = sig_wide$cluster,
            n_cells = n_cells, entropy = entropy,
            top_lineage = top_lineage, margin = margin,
            top_unbiased_genes = top_de
        )
        row_dt = cbind(row_dt, sig_wide[, ..lineage_names])
        all_rows[[paste0(run$ts, "_", cc)]] = row_dt
    }
    rm(banksy_spe); gc()
}

summary_dt = rbindlist(all_rows)
# Sort by resolution (numeric) then margin, so heuristic comparisons across
# the whole sweep read top-to-bottom by resolution first.
summary_dt[, res_numeric := as.numeric(sub("res", "", resolution))]
setorder(summary_dt, res_numeric, margin)
summary_dt[, res_numeric := NULL]

CODE_ROOT = paste0(workdir, "youyun/plgg/code")
out_tsv = paste0(CODE_ROOT, "/niche/banksy_workflow/cluster_summary_stats_full_sweep.tsv")
out_rds = paste0(CODE_ROOT, "/niche/banksy_workflow/cluster_summary_stats_scores_full_sweep.rds")
fwrite(summary_dt, out_tsv, sep = "\t")
saveRDS(score_mats, out_rds)

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Wrote ', out_tsv, ' and ', out_rds))
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Cluster counts per resolution:'))
print(summary_dt[, .N, by = resolution])
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Most borderline calls overall (lowest margin):'))
print(summary_dt[order(margin)][1:15, .(resolution, cluster, n_cells, entropy, top_lineage, margin, top_unbiased_genes)])
