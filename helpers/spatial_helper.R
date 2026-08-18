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
    library(presto)  # wilcoxauc(): single-cell Wilcoxon-AUC, used in annotate_cell_types()/generate_qc_plots()
    library(pals)    # polychrome()/alphabet() categorical palettes, used in generate_qc_plots()
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
        workdir,'youyun/plgg/data/metadata/Xenium_PS.xlsx'
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

# PLGG marker panel for automated cell-type calling, ported from
# evaluate_leiden_sweep.qmd's validated setup chunk (2026-08-17) -- see that
# qmd for the full curation rationale (panel-coverage checks, why radial_glia
# is trimmed to PAX6+SOX2, etc). Kept as a separate copy rather than shared
# code so this production path doesn't depend on the exploratory qmd; if one
# changes, check whether the other should too.
marker_panel = list(
    astro_glial   = c("AQP4", "APOE", "GJA1"),
    radial_glia   = c("PAX6", "SOX2"),
    opc           = c("OLIG1", "OLIG2", "PDGFRA", "CSPG4", "SOX10", "BCAN", "PTPRZ1"),
    oligo_mature  = c("MOG", "MOBP", "CLDN11", "MAG", "ST18"),
    neuronal      = c("SLC17A7", "SLC17A6", "RORB", "GAD1", "GAD2", "SOX11"),
    myeloid       = c("P2RY12", "CX3CR1", "AIF1", "CD68", "ITGAM", "CD163"),
    t_cell        = c("CD4", "CD2", "TRAC", "IL7R"),
    endothelial   = c("PECAM1", "FLT1"),
    proliferation = c("MKI67", "TOP2A", "PCNA", "CENPF")
)

# Per-cluster top lineage for one cluster column: single-cell Wilcoxon AUC
# (one-vs-rest) -> mean(auc - 0.5) within each lineage's panel-present
# markers -> argmax, ties broken by marker_panel list order. Mirrors
# evaluate_leiden_sweep.qmd's validated 2c.5 logic (plain mean-AUC, no
# single-marker gating).
score_lineages = function(expr_matrix, cluster_vec, panel = marker_panel) {
    genes_all = unlist(panel, use.names = FALSE)
    genes_present = intersect(genes_all, rownames(expr_matrix))
    lineage_of = setNames(rep(names(panel), lengths(panel)), genes_all)

    auc_res = as.data.table(presto::wilcoxauc(
        expr_matrix[genes_present, , drop = FALSE], factor(cluster_vec)
    ))
    setnames(auc_res, c('feature', 'group'), c('gene', 'cluster'))
    auc_res[, lineage := lineage_of[gene]]

    sig_long = auc_res[, .(score = mean(auc - 0.5)), by = .(cluster, lineage)]
    sig_wide = dcast(sig_long, cluster ~ lineage, value.var = 'score')
    lineage_names = names(panel)
    m = as.matrix(sig_wide[, ..lineage_names])
    top = lineage_names[max.col(replace(m, is.na(m), -Inf), ties.method = 'first')]
    data.table(cluster = sig_wide$cluster, top_lineage = top)
}

# Shared filename-building logic between annotate_cell_types() and
# generate_qc_plots(), matching banksy_clustering()'s own param_string
# format exactly (kct_/resct_/kni_/resni_) so all three agree on a run's
# filenames without banksy_clustering() needing to return param_string.
cluster_param_string = function(k_geom_vec, lambda_vec, pc_val, k_leiden_celltype, resolution_celltype, k_leiden_niche, resolution_niche) {
    paste0(
        'k_geom_', paste0(k_geom_vec, collapse = '_'), '_',
        'pc_', pc_val, '_',
        'lam_', paste0(lambda_vec, collapse = '_'), '_',
        'kct_', paste0(k_leiden_celltype, collapse = '_'), '_',
        'resct_', paste0(resolution_celltype, collapse = '_'), '_',
        'kni_', paste0(k_leiden_niche, collapse = '_'), '_',
        'resni_', paste0(resolution_niche, collapse = '_')
    )
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


# Leiden clustering on the Harmony-corrected BANKSY PCs, run once per lambda.
# Cell-type (small lambda) and niche (large lambda) embeddings are clustered at
# DIFFERENT (k_leiden, resolution) grids because the natural granularity is
# different: cell types are finer than spatial niches. The smaller of the two
# lambdas is treated as cell-type, the larger as niche; this assumes lambda_vec
# has length exactly 2, matching the upstream BANKSY workflow.
banksy_clustering = function(
    banksy_spe, aname, seed_val,
    # not a clustering hyperparameter, but needed for the output file name
    k_geom_vec, lambda_vec, pc_val,
    k_leiden_celltype, resolution_celltype,
    k_leiden_niche, resolution_niche,
    output_dir, current_time
){
    stopifnot(length(lambda_vec) == 2)
    set.seed(seed_val)
    subset_indices = sample(1:ncol(banksy_spe), ncol(banksy_spe) * 0.05)
    # Map lambda -> (k_leiden, resolution) grid. Smaller lambda = cell-type,
    # larger lambda = niche. Done once so the lapply below can index by value.
    lam_celltype = min(as.numeric(lambda_vec))
    lam_niche    = max(as.numeric(lambda_vec))
    leiden_grid_by_lambda = list()
    leiden_grid_by_lambda[[as.character(lam_celltype)]] = list(
        k_leiden = k_leiden_celltype, resolution = resolution_celltype
    )
    leiden_grid_by_lambda[[as.character(lam_niche)]] = list(
        k_leiden = k_leiden_niche, resolution = resolution_niche
    )
    # Encode both grids in the output filename so different sweeps don't collide.
    param_string = paste0(
        'k_geom_', paste0(k_geom_vec,collapse = '_'), '_',
        'pc_', pc_val, '_',
        'lam_', paste0(lambda_vec,collapse = '_'), '_',
        'kct_', paste0(k_leiden_celltype,collapse = '_'), '_',
        'resct_', paste0(resolution_celltype,collapse = '_'), '_',
        'kni_', paste0(k_leiden_niche,collapse = '_'), '_',
        'resni_', paste0(resolution_niche,collapse = '_')
    )
    # Leiden clustering ===================================================
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Clustering...'))
    lapply(lambda_vec, function(x){
        grid = leiden_grid_by_lambda[[as.character(x)]]
        print(paste0(
            '[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ',
            'Clustering lam', x, ' with k=', paste0(grid$k_leiden, collapse=','),
            ' res=', paste0(grid$resolution, collapse=','), ' ...'
        ))
        # running Leiden clustering on the Harmony corrected PCA loadings
        banksy_spe <<- Banksy::clusterBanksy(
            banksy_spe, dimred = paste0("Harmony_BANKSY_lam", x),
            k_neighbors = grid$k_leiden,
            resolution = grid$resolution,
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

# Writes an automated cell-type call directly into colData, replacing the
# manual hand-curated lookup-table annotation previously done per-run in
# banksy_clusters_proseg.qmd. Only annotates cell-typing (smaller-lambda)
# cluster columns -- there can be more than one now that it gets multiple
# candidate resolutions -- not niche columns, matching
# cell_type_marker_ident()'s lowest-lambda convention.
annotate_cell_types = function(
    banksy_spe, aname, k_geom_vec, lambda_vec, pc_val,
    k_leiden_celltype, resolution_celltype,
    k_leiden_niche, resolution_niche,
    output_dir, current_time
){
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Annotating cell types...'))
    param_string = cluster_param_string(k_geom_vec, lambda_vec, pc_val,
        k_leiden_celltype, resolution_celltype, k_leiden_niche, resolution_niche)
    cell_type_clusters = grep(
        paste0('lam', min(as.numeric(lambda_vec))),
        clusterNames(banksy_spe), value = TRUE
    )
    expr = assay(banksy_spe, aname)
    for (cc in cell_type_clusters) {
        top_lineage_dt = score_lineages(expr, colData(banksy_spe)[[cc]])
        lin_map = setNames(top_lineage_dt$top_lineage, top_lineage_dt$cluster)
        colData(banksy_spe)[[paste0('cell_type_', cc)]] = lin_map[as.character(colData(banksy_spe)[[cc]])]
    }
    saveRDS(banksy_spe, paste0(
        output_dir,'/banksy_clusters_connected_annotated_',
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
            assay(banksy_spe, aname),
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

# Verification plots for a clustering run, saved directly as PNGs (this
# runs inside the same non-interactive qsub job that produced the clustered
# object, so there's no benefit to a separate rendered notebook -- avoids a
# second load of the multi-million-cell object). Replaces
# banksy_clusters_proseg.qmd's manual verification workflow with the
# equivalent automated views from evaluate_leiden_sweep.qmd: per cluster
# column (both lambdas), a cluster/sample UMAP pair over a sample-mixing
# entropy barchart; for lam1 (cell-typing) columns only, a UMAP of the
# annotate_cell_types() call and the marker AUC heatmap that drives it.
# UMAP scatters are downsampled via plot_frac for render time/file size on
# the full multi-million-cell object; entropy/AUC computations use all cells.
generate_qc_plots = function(
    banksy_spe, k_geom_vec, lambda_vec, pc_val,
    k_leiden_celltype, resolution_celltype,
    k_leiden_niche, resolution_niche,
    output_dir, current_time, plot_frac = 0.1, seed_val = 55555
){
    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','Generating QC plots...'))
    param_string = cluster_param_string(k_geom_vec, lambda_vec, pc_val,
        k_leiden_celltype, resolution_celltype, k_leiden_niche, resolution_niche)
    qc_dir = paste0(output_dir, '/qc_plots_', param_string, '_', current_time)
    dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

    set.seed(seed_val)
    plot_idx = sort(sample(ncol(banksy_spe), round(ncol(banksy_spe) * plot_frac)))

    n_samples = length(unique(colData(banksy_spe)$sample_id))
    entropy_norm = function(x) {
        p = prop.table(table(x))
        p = p[p > 0]
        -sum(p * log(p)) / log(n_samples)
    }

    # pals::polychrome() errors above 36 colors -- fine for samples, but at
    # full-cohort scale some (k, res) combos produce far more clusters than
    # that, so interpolate beyond its max instead of erroring.
    categorical_pal = function(n) {
        if (n <= 36) pals::polychrome(n) else grDevices::colorRampPalette(pals::polychrome(36))(n)
    }

    samp_levels = sort(unique(as.character(colData(banksy_spe)$sample_id)))
    sample_pal  = structure(categorical_pal(length(samp_levels)), names = samp_levels)

    umap_scatter = function(coords_dt, vals, title, palette, show_labels = TRUE) {
        d = copy(coords_dt)[, grp := factor(vals[plot_idx])]
        p = ggplot(d, aes(umap_1, umap_2, colour = grp)) +
            geom_point(size = 0.1, alpha = 0.5) +
            scale_colour_manual(values = palette, name = NULL, na.value = "grey80") +
            guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
            labs(title = title, x = "UMAP 1", y = "UMAP 2") +
            theme_bw()
        if (show_labels) {
            cent = d[, .(x = median(umap_1), y = median(umap_2)), by = grp]
            p = p + geom_text(data = cent, aes(x = x, y = y, label = grp),
                              inherit.aes = FALSE, colour = "black", size = 3.2)
        }
        p
    }

    # Cluster/sample UMAP pair + entropy barchart, per cluster column, both lambdas.
    for (lam in lambda_vec) {
        umap_name = paste0("UMAP_Harmony_BANKSY_lam", lam)
        if (!umap_name %in% reducedDimNames(banksy_spe)) next
        coords = as.data.table(reducedDim(banksy_spe, umap_name))[plot_idx, 1:2]
        setnames(coords, c("umap_1", "umap_2"))

        lam_cluster_cols = grep(paste0('lam', lam), clusterNames(banksy_spe), value = TRUE)
        for (cc in lam_cluster_cols) {
            cl = as.character(colData(banksy_spe)[[cc]])
            cl_levels = sort(unique(cl))
            clpal = structure(categorical_pal(length(cl_levels)), names = cl_levels)
            p_cluster = umap_scatter(coords, cl, paste0(cc, "\nclusters"), clpal)
            p_sample = umap_scatter(coords, as.character(colData(banksy_spe)$sample_id),
                                    "sample_id", sample_pal, show_labels = FALSE)

            by_clust = split(colData(banksy_spe)$sample_id, colData(banksy_spe)[[cc]])
            ent = vapply(by_clust, entropy_norm, numeric(1))
            eb = data.table(cluster = names(ent), entropy = ent)[order(entropy)]
            eb[, cluster := factor(cluster, levels = cluster)]
            p_entropy = ggplot(eb, aes(cluster, entropy, fill = entropy < 0.30)) +
                geom_col() +
                geom_hline(yintercept = 0.30, linetype = "dashed") +
                scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                                  name = "entropy < 0.30") +
                labs(x = "cluster", y = "normalized entropy",
                     title = "per-cluster sample-mixing entropy") +
                theme_bw()

            ggsave(
                paste0(qc_dir, '/umap_sample_entropy_', gsub('clust_Harmony_BANKSY_', '', cc), '.png'),
                cowplot::plot_grid(
                    cowplot::plot_grid(p_cluster, p_sample, nrow = 1),
                    p_entropy, ncol = 1, rel_heights = c(2, 1)
                ),
                width = 13, height = 9, dpi = 150
            )
        }
    }

    # Cell-type UMAP + marker AUC heatmap: lam1 (cell-typing) columns only.
    lam1 = min(as.numeric(lambda_vec))
    cell_type_clusters = grep(paste0('lam', lam1), clusterNames(banksy_spe), value = TRUE)
    umap_name = paste0("UMAP_Harmony_BANKSY_lam", lam1)
    coords = as.data.table(reducedDim(banksy_spe, umap_name))[plot_idx, 1:2]
    setnames(coords, c("umap_1", "umap_2"))

    lineage_pal = structure(pals::alphabet(length(marker_panel)), names = names(marker_panel))
    lineage_umaps = Filter(Negate(is.null), lapply(cell_type_clusters, function(cc) {
        ct_col = paste0('cell_type_', cc)
        if (!ct_col %in% colnames(colData(banksy_spe))) return(NULL)
        umap_scatter(coords, colData(banksy_spe)[[ct_col]],
                    gsub('clust_Harmony_BANKSY_', '', cc), lineage_pal, show_labels = FALSE)
    }))
    if (length(lineage_umaps) > 0) {
        ggsave(
            paste0(qc_dir, '/celltype_umap.png'),
            cowplot::plot_grid(plotlist = lineage_umaps, nrow = 1),
            width = 6 * length(lineage_umaps), height = 5.5, dpi = 150
        )
    }

    genes_all = unlist(marker_panel, use.names = FALSE)
    genes_present = intersect(genes_all, rownames(banksy_spe))
    lineage_of = setNames(rep(names(marker_panel), lengths(marker_panel)), genes_all)
    expr_all = assay(banksy_spe, "normcounts")[genes_present, , drop = FALSE]
    auc_long = rbindlist(lapply(cell_type_clusters, function(cc) {
        res = as.data.table(presto::wilcoxauc(expr_all, factor(colData(banksy_spe)[[cc]])))
        res[, cluster_col := cc]
        res
    }))
    setnames(auc_long, c("feature", "group"), c("gene", "cluster"))
    auc_long[, lineage := lineage_of[gene]]
    auc_long[, gene := factor(gene, levels = genes_all)]
    auc_long[, cluster_col_short := sub("clust_Harmony_BANKSY_", "", cluster_col)]

    p_auc = ggplot(auc_long, aes(x = cluster, y = gene, fill = auc)) +
        geom_tile() +
        facet_grid(lineage ~ cluster_col_short, scales = "free", space = "free") +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                             midpoint = 0.5, limits = c(0, 1), name = "AUC\n(one-vs-rest)") +
        labs(title = "Marker specificity per cluster (single-cell Wilcoxon AUC)",
             x = "cluster", y = NULL) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 6), strip.text.y = element_text(angle = 0))
    ggsave(paste0(qc_dir, '/marker_auc_heatmap.png'), p_auc, width = 10, height = 5, dpi = 150)

    print(paste0('[',format(Sys.time(), "%Y/%m/%d-%H:%M:%S"),'] | ','QC plots saved to: ', qc_dir))
}