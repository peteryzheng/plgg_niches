library(optparse)

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


if(!interactive()) {
    option_list = list(
        make_option(
            c('--k1'), type = 'integer', default = 15,
            help = 'k neighbors for BANKSY to calculate mean neighborhood expression', 
            metavar = 'k1'
        ),
        make_option(
            c('--k2'), type = 'integer', default = 30,
            help = 'k neighbors for BANKSY to calculate neighborhood expression gradient', 
            metavar = 'k2'
        ),
        make_option(
            c('--lam1'), type = 'numeric', default = 0.2,
            help = 'lower lambda for BANKSY to do cell typing', 
            metavar = 'lambda1'
        ),
        make_option(
            c('--lam2'), type = 'numeric', default = 0.8,
            help = 'higher lambda for BANKSY to do domain segmentation', 
            metavar = 'lambda2'
        ),
        make_option(
            c('--npc'), type = 'integer', default = 20,
            help = 'number of PCs to use for BANKSY', 
            metavar = 'npc'
        ),
        # Cell-type (small lambda) leiden params. One k + a comma-separated
        # list of resolutions; clusterBanksy expands the (k x res) grid.
        # k=50 is BANKSY-recommended for ~90k+ cells; res 0.5 sat on the widest
        # plateau in the prior sweep, res 1 is kept for continuity with earlier
        # downstream analyses.
        make_option(
            c('--k_ct'), type = 'integer', default = 50,
            help = 'k_neighbors for cell-type leiden clustering (lam1)',
            metavar = 'k'
        ),
        make_option(
            c('--res_ct'), type = 'character', default = '0.5,1',
            help = 'comma-separated resolutions for cell-type leiden clustering (lam1)',
            metavar = 'res_list'
        ),
        # Niche (large lambda) leiden params. Same k/res as cell type so the
        # two sides are compared on equal-footing grid; downstream analysis
        # picks the resolution that yields contiguous niches per sample.
        make_option(
            c('--k_ni'), type = 'integer', default = 50,
            help = 'k_neighbors for niche leiden clustering (lam2)',
            metavar = 'k'
        ),
        make_option(
            c('--res_ni'), type = 'character', default = '0.5,1',
            help = 'comma-separated resolutions for niche leiden clustering (lam2)',
            metavar = 'res_list'
        ),
        make_option(
            c('--seed'), type = 'integer', default = 55555,
            help = 'seed value for reproducibility', 
            metavar = 'seed'
        ),
        make_option(c("-o", "--outputdir"),
            type = "character", default = paste0(workdir, "coja/Spatial_PLGG/data/Xenium/Xenium_Objects/"),
            help = "Output directory to use.", metavar = "outputdir"
        )
    )

    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)

    k_geom = c(opt$k1, opt$k2)
    lambda = c(opt$lam1, opt$lam2)
    npc = opt$npc
    # Parse comma-separated resolutions into numeric vectors; trimws guards
    # against accidental whitespace when called from a shell wrapper.
    k_ct  = opt$k_ct
    res_ct = as.numeric(strsplit(trimws(opt$res_ct), ',')[[1]])
    k_ni  = opt$k_ni
    res_ni = as.numeric(strsplit(trimws(opt$res_ni), ',')[[1]])
    stopifnot(all(!is.na(res_ct)), all(!is.na(res_ni)))
    seed = opt$seed
    output_dir = opt$outputdir

    # source helper functions
    source(paste0(workdir, "youyun/plgg/code/helpers/spatial_helper.R"))
    current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
    
    # LOAD DATA ===========================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Loading data...'))
    total_se = loading_data(segmentation_method = 'proseg')

    # STAGGER SPATIAL COORDINATES ==========================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Staggering spatial coordinates...'))
    total_se_staggered = stagger_spatial_coords(total_se)

    # SUBSET to only gene expression and not control probes ===============
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Subsetting to gene expression...'))
    total_se_staggered = total_se_staggered[rowData(total_se_staggered)$feature_type == 'Gene Expression',]

    # QC and NORMALIZATION =================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','QC and normalization...'))
    total_se_staggered = QC_and_normalize(total_se_staggered)

    # BANKSY ==============================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Running BANKSY workflow...'))
    total_se_staggered = banksy_workflow(
        total_se_staggered,
        aname = 'normcounts', seed_val = seed, 
        k_geom_vec = k_geom,
        lambda_vec = lambda, 
        pc_val = npc,
        output_dir, current_timestamp
    )

    # Clustering ================================================================
    # banksy_clustering expects separate (k, res) per lambda so cell-type and
    # niche grids can be tuned independently.
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Running clustering on BANKSY output...'))
    total_se_staggered = banksy_clustering(
        total_se_staggered,
        aname = 'normcounts', seed_val = seed,
        k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
        k_leiden_celltype = k_ct, resolution_celltype = res_ct,
        k_leiden_niche    = k_ni, resolution_niche    = res_ni,
        output_dir, current_timestamp
    )

    # Find Cell Type Markers ===================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finding cell type markers...'))
    cell_type_marker_ident(
        total_se_staggered,
        aname = 'normcounts', seed_val = seed,
        k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
        output_dir, current_timestamp
    )

    # Automated Cell Type Annotation ==========================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Annotating cell types...'))
    total_se_staggered = annotate_cell_types(
        total_se_staggered,
        aname = 'normcounts',
        k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
        k_leiden_celltype = k_ct, resolution_celltype = res_ct,
        k_leiden_niche    = k_ni, resolution_niche    = res_ni,
        output_dir, current_timestamp
    )

    # QC Plots =================================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Generating QC plots...'))
    generate_qc_plots(
        total_se_staggered,
        k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
        k_leiden_celltype = k_ct, resolution_celltype = res_ct,
        k_leiden_niche    = k_ni, resolution_niche    = res_ni,
        output_dir, current_timestamp, seed_val = seed
    )

    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Done!'))
}