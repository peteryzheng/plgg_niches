# Assembles the final production annotated object at the chosen resolutions
# (lam0.2 cell-type res=0.2, lam0.8 niche res=0.5 -- see param_search.tsv,
# the source of truth for these values) via a validated fast-path rather
# than a fresh `qsub banksy_cohort_qsub.sh` submission.
#
# Provenance chain this script depends on:
#   1. banksy_cohort_qsub.sh -> banksy_cohort.R (2026-08-18, job on argos6)
#      produced the base object with BANKSY/Harmony/UMAP embeddings and the
#      original lam0.2{0.5,1}/lam0.8{0.5,1} Leiden clusters, connectClusters()
#      already applied to all of them:
#        banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.5_1_kni_50_resni_0.5_1_20260818_005830.rds
#   2. celltype_lowres_recluster.R --resolution 0.2,0.4 (SGE job 3625399,
#      2026-08-24) reused that object's Harmony_BANKSY_lam0.2 embedding to
#      cheaply re-cluster lam0.2 at new candidate resolutions (Leiden is
#      deterministic given seed=55555, so this is the actual, only, compute
#      -- BANKSY/Harmony/UMAP are not recomputed). It calls
#      Banksy::clusterBanksy() directly rather than going through
#      banksy_clustering(), so it does NOT call connectClusters() on the
#      new lam0.2 columns (see its own header comment / XENIUM_RUNBOOK.md
#      for why this script is cell-type-axis-only). Produced:
#        banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_0.4_kni_50_resni_0.5_1_20260824_042548.rds
#      (this is the object THIS script loads -- it already carries both the
#      chosen lam0.2 res=0.2 clustering and the untouched, already-connected
#      lam0.8 res=0.5 clustering inherited from the base object).
#
# Why this is equivalent to a fresh full run (not just "close enough"):
#   BANKSY matrix computation, Harmony batch correction, and Leiden
#   clustering are all deterministic given the same input data + fixed seed
#   (55555, unchanged throughout). The one step a full run would additionally
#   apply here -- Banksy::connectClusters() -- is a pure post-hoc integer
#   cluster-ID relabeling (Hungarian-matches each clust_* column's labels
#   against a seed column to maximize overlap); read its source directly to
#   confirm it does not touch cluster membership, counts, or boundaries, and
#   completes in ~seconds (confirmed via the original run's own log). Since
#   annotate_cell_types() maps clusters to cell-type labels via marker
#   expression per cluster (not via the raw integer ID), skipping it upstream
#   has zero effect on any cell-type or niche call. This script closes that
#   one cosmetic gap directly.
#
# Steps: load -> drop the resolutions not chosen as final (lam0.2 res0.4,
# lam0.8 res1) -> connectClusters() -> annotate_cell_types() (uses
# spatial_helper.R's default heuristic, score_lineages_top_gene() with
# min_gene_auc=0.55) -> regenerate QC/review outputs.

workdir <- Sys.getenv("HOME")
if (workdir %in% c("/Users/youyun", "/Users/youyunzheng")) {
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
} else if (workdir == "/PHShome/yz762") {
    workdir <- "/data/beroukhim1/"
} else if (workdir == "/home/yz762") {
    workdir <- "/mnt/storage/dept/medonc/beroukhim/"
} else {
    workdir <- "/xchip/beroukhimlab/"
}
source(paste0(workdir, "youyun/plgg/code/helpers/spatial_helper.R"))

outputdir = paste0(workdir, "youyun/plgg/data/banksy_param_search/k1_15_k2_30_lambda1_0.2_lambda2_0.8_npcs_20_kct_50_resct_0.5,1_kni_50_resni_0.5,1")
rds = paste0(outputdir, "/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_0.4_kni_50_resni_0.5_1_20260824_042548.rds")

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loading ', rds, ' ...'))
banksy_spe = readRDS(rds)
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loaded. ncol = ', ncol(banksy_spe)))

# Drop the resolutions not chosen as final so the final object only carries
# the two production resolutions (param_search.tsv: res_ct=0.2, res_ni=0.5).
drop_cols = c(
    "clust_Harmony_BANKSY_lam0.2_k50_res0.4",
    "cell_type_clust_Harmony_BANKSY_lam0.2_k50_res0.4",
    "clust_Harmony_BANKSY_lam0.8_k50_res1"
)
drop_cols = intersect(drop_cols, colnames(colData(banksy_spe)))
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Dropping non-final-resolution columns: ', paste(drop_cols, collapse = ', ')))
colData(banksy_spe)[, drop_cols] = NULL

# Closes the one real (cosmetic) gap vs. a full banksy_clustering() run --
# see header comment above. Cheap: expect seconds, not minutes.
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Connecting clusters...'))
t0 = Sys.time()
banksy_spe = Banksy::connectClusters(banksy_spe, map_to = clusterNames(banksy_spe)[1])
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | connectClusters() took ', round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), 's'))

current_time = format(Sys.time(), "%Y%m%d_%H%M%S")
k_geom = c(15, 30); lambda = c(0.2, 0.8); npc = 20
k_ct = 50; res_ct = 0.2; k_ni = 50; res_ni = 0.5

# Re-annotate: drops the existing (already-correct) cell_type_ column for
# res0.2 and regenerates it fresh so it's computed post-connectClusters()
# under the same integer cluster IDs saved below (score_lineages_top_gene()
# maps by marker expression per cluster, not raw ID, so the calls themselves
# are unaffected -- this just keeps the saved object internally consistent).
banksy_spe = annotate_cell_types(
    banksy_spe, aname = 'normcounts',
    k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
    k_leiden_celltype = k_ct, resolution_celltype = res_ct,
    k_leiden_niche = k_ni, resolution_niche = res_ni,
    output_dir = outputdir, current_time = current_time
)

for (cc in grep('lam0.2', clusterNames(banksy_spe), value = TRUE)) {
    ct_col = paste0('cell_type_', cc)
    tab = table(colData(banksy_spe)[[ct_col]])
    print(paste0(cc, ': ', paste(names(tab), tab, sep = '=', collapse = ', ')))
}

generate_qc_plots(
    banksy_spe,
    k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
    k_leiden_celltype = k_ct, resolution_celltype = res_ct,
    k_leiden_niche = k_ni, resolution_niche = res_ni,
    output_dir = outputdir, current_time = current_time, seed_val = 55555
)

generate_celltype_review_panels(
    banksy_spe, aname = 'normcounts',
    k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
    k_leiden_celltype = k_ct, resolution_celltype = res_ct,
    k_leiden_niche = k_ni, resolution_niche = res_ni,
    output_dir = outputdir, current_time = current_time, seed_val = 55555
)

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Done. Final timestamp: ', current_time))
