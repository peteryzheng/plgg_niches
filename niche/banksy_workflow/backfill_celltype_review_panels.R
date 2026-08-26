library(optparse)

# One-off/rerunnable backfill: adds generate_celltype_review_panels() output
# (combined cell-type-review panel per cell-type resolution) to a
# banksy_cohort.R run's existing qc_plots_* directory, without re-running the
# rest of the (multi-hour) BANKSY/Harmony/clustering pipeline. Reads the
# already-saved *_clusters_connected_annotated_*.rds object for that run and
# calls the same helper the main pipeline now calls going forward.

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
        help = 'path to a banksy_clusters_connected_annotated_*.rds file from banksy_cohort.R',
        metavar = 'path'),
    make_option(c('--k1'), type = 'integer', default = 15, metavar = 'k1'),
    make_option(c('--k2'), type = 'integer', default = 30, metavar = 'k2'),
    make_option(c('--lam1'), type = 'numeric', default = 0.2, metavar = 'lambda1'),
    make_option(c('--lam2'), type = 'numeric', default = 0.8, metavar = 'lambda2'),
    make_option(c('--npc'), type = 'integer', default = 20, metavar = 'npc'),
    make_option(c('--k_ct'), type = 'integer', default = 50, metavar = 'k'),
    make_option(c('--res_ct'), type = 'character', default = '0.5,1', metavar = 'res_list'),
    make_option(c('--k_ni'), type = 'integer', default = 50, metavar = 'k'),
    make_option(c('--res_ni'), type = 'character', default = '0.5,1', metavar = 'res_list'),
    make_option(c('--seed'), type = 'integer', default = 55555, metavar = 'seed'),
    # current_time must match the timestamp embedded in annotated_rds's filename
    # so cluster_param_string() resolves to the same qc_plots_* directory.
    make_option(c('--current_time'), type = 'character',
        help = 'run timestamp (YYYYMMDD_HHMMSS) matching annotated_rds filename', metavar = 'timestamp'),
    make_option(c('-o', '--outputdir'), type = 'character',
        help = 'banksy_param_search/<param combo> directory containing the run', metavar = 'outputdir')
)
opt = parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$annotated_rds), !is.null(opt$current_time), !is.null(opt$outputdir))

source(paste0(workdir, "youyun/plgg/code/helpers/spatial_helper.R"))

k_geom = c(opt$k1, opt$k2)
lambda = c(opt$lam1, opt$lam2)
res_ct = as.numeric(strsplit(trimws(opt$res_ct), ',')[[1]])
res_ni = as.numeric(strsplit(trimws(opt$res_ni), ',')[[1]])

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Loading ', opt$annotated_rds, ' ...'))
banksy_spe = readRDS(opt$annotated_rds)

generate_celltype_review_panels(
    banksy_spe, aname = 'normcounts',
    k_geom_vec = k_geom, lambda_vec = lambda, pc_val = opt$npc,
    k_leiden_celltype = opt$k_ct, resolution_celltype = res_ct,
    k_leiden_niche    = opt$k_ni, resolution_niche    = res_ni,
    output_dir = opt$outputdir, current_time = opt$current_time, seed_val = opt$seed
)

print(paste0('[', format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), '] | Backfill complete.'))
