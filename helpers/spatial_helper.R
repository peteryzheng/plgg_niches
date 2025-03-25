suppressPackageStartupMessages({
    library(rhdf5)
    library(Matrix)
    library(SpatialExperiment)
    library(data.table)
    library(ggplot2)
    library(Banksy)
    library(SummarizedExperiment)
    library(scuttle)
    library(scater)
    library(cowplot)
    library(ggplot2)
    library(harmony)
    library(scran)
    library(Seurat)
})


xenium2SPE = function(data_dir){
    # load gene expression info ========================
    # https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs#feature-matrix
    # the cell feature matrix only include transcripts that pass the default quality value (Q-Score) threshold of Q20.
    exp_file = paste0(data_dir, "/cell_feature_matrix.h5")
    exp_data <- h5read(exp_file, "matrix")
    ## from https://gist.github.com/slowkow/d3c4b77c9bf2a75f6dad4843d7d3aefc
    counts <- sparseMatrix(
        dims = exp_data$shape,
        i = as.numeric(exp_data$indices),
        p = as.numeric(exp_data$indptr),
        x = as.numeric(exp_data$data),
        index1 = FALSE
    )
    colnames(counts) <- exp_data$barcodes
    rownames(counts) <- exp_data$features$name
    tmp = data.frame(exp_data$features)

    # load the spatial info ========================
    pos_info = fread(cmd = paste0(
        'zcat < ',data_dir, "/cells.csv.gz"
    ))

    # all together now ========================
    return(SpatialExperiment(
        assays = list(counts = counts),
        # feature/gene metadata
        rowData = data.frame(exp_data$features),
        # observation/cell metadata
        colData = pos_info,
        spatialCoordsNames = c("x_centroid", "y_centroid")
    ))
}

loading_data = function(segmentation_method = c('proseg','default')){
    if(segmentation_method == 'default'){
        # default segmentation
        xenium_dirs = list.files(paste0(
            workdir,'coja/Spatial_PLGG/data/Xenium/Xenium_Analyzer/'
        ), full.names = TRUE)#[c(1,7)]
    }else if(segmentation_method == 'proseg'){
        # proseg segmentation
        xenium_dirs = system(paste0(
            'find ', workdir,
            'youyun/plgg/data/segmentation/proseg_run_121024/*/*_proseg/outs',
            ' -name outs'
        ), intern = TRUE)
    }else{
        stop('segmentation_method must be one of "proseg" or "default"')
    }

    metadata = data.table(readxl::read_excel(paste0(
        workdir,'coja/Spatial_PLGG/data/metadata/Xenium_PS.xlsx'
    )))[
        ,batch := as.factor(gsub('__.*','',gsub('output-XETG[0-9]+__','',file)))
    ]
    

    total_se = do.call(
        'cbind',
        lapply(xenium_dirs, function(x){
            tmp_se = xenium2SPE(x)
            file_name = ifelse(
                segmentation_method == 'default',
                basename(x), 
                basename(gsub('_proseg/outs','',x))
            )
            sample_id = metadata[file == file_name]$id
            sample_idat = metadata[file == file_name]$idat
            histology = metadata[file == file_name]$mc
            alteration = metadata[file == file_name]$alt
            batch = metadata[file == file_name]$batch
            tmp_se$sample_id = sample_id
            tmp_se$sample_idat = sample_idat
            tmp_se$histology = histology
            tmp_se$alteration = alteration
            tmp_se$batch = batch
            colnames(tmp_se) = paste0(colnames(tmp_se), '__', sample_id)
            return(tmp_se)
        }
    ))
    return(total_se)
}

stagger_spatial_coords = function(spe){
    # Stagger spatial coordinates
    locs <- spatialCoords(spe)
    locs <- cbind(locs, sample_id = factor(spe$sample_id))
    locs_dt <- data.table(locs)
    colnames(locs_dt) <- c("sdimx", "sdimy", "group")
    locs_dt[, sdimx := sdimx - min(sdimx), by = group]
    global_max <- max(locs_dt$sdimx) * 1.5
    locs_dt[, sdimx := sdimx + group * global_max]
    locs <- as.matrix(locs_dt[, 1:2])
    rownames(locs) <- colnames(spe)
    spatialCoords(spe) <- locs
    return(spe)
}

is_outlier = function(metric, nmad = 3, threshold = 'both'){
    if(!threshold %in% c('both','upper','lower')){
        stop('threshold must be one of "both", "upper", or "lower"')
    }
    # inspired by the MAD idea from single cell best practice 
    # https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-cells
    metric = as.numeric(metric)
    mad = median(abs(metric - median(metric)))
    if(threshold == 'both'){
        return(abs(metric - median(metric)) > nmad * mad)
    }else if(threshold == 'upper'){
        return(metric - median(metric) > nmad * mad)
    }else if(threshold == 'lower'){
        return(metric - median(metric) < -nmad * mad)
    }
}

QC_and_normalize = function(spe){
    qcstats <- perCellQCMetrics(spe)
    
    # outliers for log1p total counts -- upper lower
    outliers = is_outlier(log(qcstats$total + 1), nmad = 3, threshold = 'both') | 
        # outliers for log1p genes by count -- lower
        is_outlier(log(qcstats$detected + 1), nmad = 3, threshold = 'lower')

    print(paste0(
        'Removing ', sum(outliers)/length(outliers)*100,
        '% of cells due to low quality'
    ))
    keep <- !outliers
    spe <- spe[, keep]

    # Normalization to mean library size
    spe <- computeLibraryFactors(spe)
    aname <- "normcounts"
    # log-normalize the counts
    assay(spe, aname) <- normalizeCounts(spe, log = FALSE)
    
    return(spe)
}

banksy_workflow = function(
    banksy_spe, aname, seed_val, 
    k_geom_vec, lambda_vec, pc_val,
    output_dir, current_time
){
    set.seed(seed_val)
    subset_indices = sample(1:ncol(banksy_spe), ncol(banksy_spe) * 0.05)
    # current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    param_string = paste0(
        'k_geom_', paste0(k_geom_vec,collapse = '_'), '_',
        'pc_', pc_val, '_',
        'lam_', paste0(lambda_vec,collapse = '_')
    )

    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Computing BANKSY Matrices...'))
    # calculate mean neighborhood and AGF matricies
    banksy_spe <- Banksy::computeBanksy(banksy_spe, assay_name = aname, compute_agf = TRUE, k_geom = k_geom_vec)

    # 0 for non spatial clustering, 0.2 for cell typing, and 0.8 for domain segmentation
    banksy_spe <- Banksy::runBanksyPCA(
        banksy_spe, assay_name = aname,
        use_agf = TRUE, lambda = lambda_vec, 
        npcs = pc_val, seed = seed_val
    )

    # Batch correction =====================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Batch correction...'))
    lapply(lambda_vec, function(x){
        set.seed(seed_val)
        harmony_embedding_PCA_M1 <- RunHarmony(
            data_mat = reducedDim(banksy_spe, paste0("PCA_M1_lam", x)),
            meta_data = colData(banksy_spe),
            vars_use = c('sample_id','batch'),
            do_pca = FALSE,
            max_iter = 50,
            verbose = TRUE
        )
        reducedDim(banksy_spe, paste0("Harmony_BANKSY_lam", x)) <<- harmony_embedding_PCA_M1
    })

    # UMAP ================================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Running UMAP...'))
    lapply(lambda_vec, function(x){
        # run UMAP on the Harmony corrected embeddings
        banksy_spe <<- runBanksyUMAP(banksy_spe, dimred = paste0("Harmony_BANKSY_lam", x))
    })

    # saving the objects after UMAP ========================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Saving UMAP objects...'))
    print(paste0('Saving to: ', output_dir))
    saveRDS(banksy_spe, paste0(
        output_dir,'/banksy_',
        param_string,'_',current_time,'.rds'
    ))
    saveRDS(banksy_spe[, subset_indices], paste0(
        output_dir,'/banksy_subset_',
        param_string,'_',current_time,'.rds'
    ))

    return(banksy_spe)
}


banksy_clustering = function(
    banksy_spe, aname, seed_val, 
    # not a clustering hyperparameter, but we need it for the output file name
    k_geom_vec, lambda_vec, pc_val,
    k_leiden_vec, resolution_vec,
    output_dir, current_time
){
    set.seed(seed_val)
    subset_indices = sample(1:ncol(banksy_spe), ncol(banksy_spe) * 0.05)
    # current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    param_string = paste0(
        'k_geom_', paste0(k_geom_vec,collapse = '_'), '_',
        'pc_', pc_val, '_',
        'lam_', paste0(lambda_vec,collapse = '_'), '_',
        'k_leiden_', paste0(k_leiden_vec,collapse = '_'), '_',
        'res_', paste0(resolution_vec,collapse = '_')
    )
    # Leiden clustering ===================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Clustering...'))
    lapply(lambda_vec, function(x){
        print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Clustering lam', x, ' ...'))
        # running Leiden clustering on the Harmony corrected PCA loadings
        banksy_spe <<- Banksy::clusterBanksy(
            banksy_spe, dimred = paste0("Harmony_BANKSY_lam", x), 
            k_neighbors = k_leiden_vec,
            resolution = resolution_vec, 
            algo = 'leiden',
            seed = seed_val
        )
        print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finished clustering lam', x, ' ...'))
    })
    saveRDS(banksy_spe, paste0(
        output_dir,'/banksy_clusters_',
        param_string,'_',current_time,'.rds'
    ))

    # Connect clusters ====================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Connecting clusters...'))
    banksy_spe <- Banksy::connectClusters(banksy_spe, map_to = clusterNames(banksy_spe)[1])

    # Saving final objects ================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Saving Final Objects...'))
    saveRDS(banksy_spe, paste0(
        output_dir,'/banksy_clusters_connected_',
        param_string,'_',current_time,'.rds'
    ))
    saveRDS(banksy_spe[, subset_indices], paste0(
        output_dir,'/banksy_clusters_connected_subset_',
        param_string,'_',current_time,'.rds'
    ))
    file.remove(paste0(
        output_dir,'/banksy_',
        param_string,'_',current_time,'.rds'
    ))
    file.remove(paste0(
        output_dir,'/banksy_subset_',
        param_string,'_',current_time,'.rds'
    ))
    file.remove(paste0(
        output_dir,'/banksy_clusters_',
        param_string,'_',current_time,'.rds'
    ))
    return(banksy_spe)
}

cell_type_marker_ident = function(
    banksy_spe, aname, seed_val, 
    # not a clustering hyperparameter, but we need it for the output file name
    k_geom_vec, lambda_vec, pc_val,
    output_dir, current_time
){
    # current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    param_string = paste0(
        'k_geom_', paste0(k_geom_vec,collapse = '_'), '_',
        'pc_', pc_val
    )
    # Find markers =======================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Finding markers...'))
    cell_type_clusters = grep(
        paste0('lam',min(as.numeric(lambda_vec))),
        clusterNames(banksy_spe), value = TRUE
    )
    lapply(cell_type_clusters, function(x){
        # lam lowest
        cell_type_markers = findMarkers(
            assay(banksy_spe, "counts"),
            groups = banksy_spe[[x]],
            test.type="wilcox"
        )
        saveRDS(cell_type_markers, paste0(
            output_dir,'/banksy_clusters_markers_',
            param_string,'_', gsub('.*BANKSY_','',x),'_',
            current_time,'.rds'
        ))
    })
}