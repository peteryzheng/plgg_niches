library(optparse)

# local vs UGER
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
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
        make_option(
            c('--k_leiden_lam1'), type = 'character', default = '30,50',
            help = 'comma-separated k neighbors for leiden clustering on lam1 (cell typing)',
        ),
        make_option(
            c('--k_leiden_lam2'), type = 'character', default = '30,50',
            help = 'comma-separated k neighbors for leiden clustering on lam2 (niche calling)',
        ),
        make_option(
            c('--res_lam1'), type = 'character', default = '0.75,1',
            help = 'comma-separated Leiden resolutions to try for lam1 (cell typing)',
            metavar = 'resolution'
        ),
        make_option(
            c('--res_lam2'), type = 'character', default = '0.75,1',
            help = 'comma-separated Leiden resolutions to try for lam2 (niche calling)',
            metavar = 'resolution'
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
    # One Leiden k_neighbors / resolution vector per lambda, in the same
    # order as `lambda` (lam1 = cell typing, lam2 = niche calling) -- these
    # no longer have to share one grid across both lambdas.
    k_leiden_list = list(
        as.numeric(strsplit(opt$k_leiden_lam1, ',')[[1]]),
        as.numeric(strsplit(opt$k_leiden_lam2, ',')[[1]])
    )
    resolution_list = list(
        as.numeric(strsplit(opt$res_lam1, ',')[[1]]),
        as.numeric(strsplit(opt$res_lam2, ',')[[1]])
    )
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
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Running clustering on BANKSY output...'))
    total_se_staggered = banksy_clustering(
        total_se_staggered,
        aname = 'normcounts', seed_val = seed,
        k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
        k_leiden_list = k_leiden_list, resolution_list = resolution_list,
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
        k_leiden_list = k_leiden_list, resolution_list = resolution_list,
        output_dir, current_timestamp
    )

    # QC Plots =================================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Generating QC plots...'))
    generate_qc_plots(
        total_se_staggered,
        k_geom_vec = k_geom, lambda_vec = lambda, pc_val = npc,
        k_leiden_list = k_leiden_list, resolution_list = resolution_list,
        output_dir, current_timestamp, seed_val = seed
    )

    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Done!'))
}