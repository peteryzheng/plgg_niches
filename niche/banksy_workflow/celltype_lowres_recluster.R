library(optparse)

# Re-runs Leiden clustering at new candidate resolution(s) on the FULL
# annotated object -- NOT the 90k subset. A resolution picked as
# "satisfactory" on the subset (res=0.5) turned out to over-segment on the
# full object, so any new resolution has to be tested directly on the full
# object to mean anything. Only Leiden is redone here (BANKSY matrix /
# Harmony / UMAP are already saved on the object and are not recomputed) --
# see XENIUM_RUNBOOK.md / AGENTS.md for the full pipeline context.
#
# CELL-TYPE (lam0.2) AXIS ONLY -- not just by convention, but mechanically:
# the Banksy::clusterBanksy() call below is generic and would recluster
# whichever --lam1 dimred you point it at, but the two steps after it are
# not. annotate_cell_types() and generate_celltype_review_panels() (in
# helpers/spatial_helper.R) both hardcode `min(lambda_vec)` internally to
# find "the" cell-type cluster column, regardless of which lambda was just
# reclustered. Pointing this script at lam0.8 (--lam1 0.8) would silently
# recluster niches correctly but then annotate/plot the OLD, untouched
# lam0.2 clusters instead -- no error, just wrong output. Niches also don't
# get marker-based annotation anywhere in this pipeline, so those two steps
# genuinely don't apply to them. A niche-resolution equivalent would need
# its own script: reuse the clusterBanksy() pattern below, but replace
# annotate_cell_types()/generate_celltype_review_panels() with
# connectClusters() + niche-appropriate (spatial-contiguity) QC instead.

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

option_list = list(
    make_option(c('--annotated_rds'), type = 'character',
        help = 'path to the FULL banksy_clusters_connected_annotated_*.rds file (not the subset)',
        metavar = 'path'),
    # New cell-type (lam0.2) resolution candidate(s) to test, comma-separated.
    make_option(c('--resolution'), type = 'character',
        help = 'comma-separated new lam0.2 resolution value(s) to test, e.g. "0.3"', metavar = 'res_list'),
    make_option(c('--k_ct'), type = 'integer', default = 50,
        help = 'k_neighbors for the new lam0.2 leiden clustering', metavar = 'k'),
    make_option(c('--seed'), type = 'integer', default = 55555, metavar = 'seed'),
    # Original run's params, needed only so cluster_param_string()/output
    # paths/niche columns match the existing run's directory layout.
    make_option(c('--k1'), type = 'integer', default = 15, metavar = 'k1'),
    make_option(c('--k2'), type = 'integer', default = 30, metavar = 'k2'),
    make_option(c('--lam1'), type = 'numeric', default = 0.2,
        help = 'cell-type (smaller) lambda -- do not set to 0.8 to target niches, see header comment: annotate_cell_types()/generate_celltype_review_panels() hardcode min(lambda_vec), so this script is cell-type-axis-only regardless of this flag',
        metavar = 'lambda1'),
    make_option(c('--lam2'), type = 'numeric', default = 0.8,
        help = 'niche (larger) lambda -- kept only for cluster_param_string()/output naming, not reclustered',
        metavar = 'lambda2'),
    make_option(c('--npc'), type = 'integer', default = 20, metavar = 'npc'),
    make_option(c('--k_ni'), type = 'integer', default = 50, metavar = 'k'),
    make_option(c('--res_ni'), type = 'character', default = '0.5,1', metavar = 'res_list'),
    make_option(c('-o', '--outputdir'), type = 'character',
        help = 'banksy_param_search/<param combo> directory containing the run', metavar = 'outputdir')
)
opt = parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$annotated_rds), !is.null(opt$resolution), !is.null(opt$outputdir))

source(paste0(workdir, "youyun/plgg/code/helpers/spatial_helper.R"))

k_geom = c(opt$k1, opt$k2)
lambda = c(opt$lam1, opt$lam2)
res_new = as.numeric(strsplit(trimws(opt$resolution), ',')[[1]])
res_ni  = as.numeric(strsplit(trimws(opt$res_ni), ',')[[1]])
stopifnot(all(!is.na(res_new)))
current_time = format(Sys.time(), "%Y%m%d_%H%M%S")

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loading FULL object: ', opt$annotated_rds, ' ...'))
banksy_spe = readRDS(opt$annotated_rds)
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loaded. ncol = ', ncol(banksy_spe)))

# Drop only the existing lam0.2 clust_*/cell_type_* columns -- leave lam0.8
# (niche) columns untouched so the object still carries the production
# niche calls even though this script only touches cell-type clustering.
drop_cols = grep(
    paste0("^(clust_Harmony_BANKSY_lam", opt$lam1, "|cell_type_clust_Harmony_BANKSY_lam", opt$lam1, ")"),
    colnames(colData(banksy_spe)), value = TRUE
)
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Dropping existing lam', opt$lam1, ' columns: ', paste(drop_cols, collapse = ', ')))
colData(banksy_spe)[, drop_cols] = NULL

# Leiden re-clustering on the already-saved Harmony-corrected BANKSY PCs --
# this is the expensive step (see XENIUM_RUNBOOK.md / PR notes: full-object
# Leiden clustering dominates pipeline runtime, not BANKSY/Harmony/UMAP).
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Clustering lam', opt$lam1, ' at res=', opt$resolution, ' (k=', opt$k_ct, ') on FULL object ...'))
banksy_spe = Banksy::clusterBanksy(
    banksy_spe, dimred = paste0("Harmony_BANKSY_lam", opt$lam1),
    k_neighbors = opt$k_ct, resolution = res_new,
    algo = "leiden", seed = opt$seed
)
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Finished clustering. Cluster counts: ',
    paste(sapply(res_new, function(r) {
        cc = paste0('clust_Harmony_BANKSY_lam', opt$lam1, '_k', opt$k_ct, '_res', r)
        paste0(r, '=', length(unique(colData(banksy_spe)[[cc]])))
    }), collapse = ', ')))

# Annotate the new resolution(s) -- reuses annotate_cell_types() as-is; it
# greps for lam0.2 cluster columns internally, so with the old ones dropped
# it only sees and annotates the new one(s). This also saveRDS()s a fresh
# full-object copy under a new timestamp.
print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Annotating cell types for new resolution(s)...'))
banksy_spe = annotate_cell_types(
    banksy_spe, aname = 'normcounts',
    k_geom_vec = k_geom, lambda_vec = lambda, pc_val = opt$npc,
    k_leiden_celltype = opt$k_ct, resolution_celltype = res_new,
    k_leiden_niche    = opt$k_ni, resolution_niche    = res_ni,
    output_dir = opt$outputdir, current_time = current_time
)

# Cell-type review panels for the new resolution(s), same layout as the
# production run's panels (UMAP pair + entropy-annotated marker AUC heatmap).
generate_celltype_review_panels(
    banksy_spe, aname = 'normcounts',
    k_geom_vec = k_geom, lambda_vec = lambda, pc_val = opt$npc,
    k_leiden_celltype = opt$k_ct, resolution_celltype = res_new,
    k_leiden_niche    = opt$k_ni, resolution_niche    = res_ni,
    output_dir = opt$outputdir, current_time = current_time, seed_val = opt$seed
)

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Done. New timestamp: ', current_time))
