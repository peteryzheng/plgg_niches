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
            c('--kc1'), type = 'integer', default = 30,
            help = 'lower k neighbors for leiden clustering',
        ),
        make_option(
            c('--kc2'), type = 'integer', default = 50,
            help = 'higher k neighbors for leiden clustering',
        ),
        make_option(
            c('--res1'), type = 'numeric', default = 0.75,
            help = 'lower resolution for BANKSY to do clustering', 
            metavar = 'resolution'
        ),
        make_option(
            c('--res2'), type = 'numeric', default = 1,
            help = 'higher resolution for BANKSY to do clustering', 
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
    k_leiden = c(opt$kc1, opt$kc2)
    res = c(opt$res1, opt$res2)
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
        k_leiden_vec = k_leiden, resolution_vec = res,
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

    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Done!'))
}