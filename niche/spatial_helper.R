library(rhdf5)
library(Matrix)
library(SpatialExperiment)
library(data.table)

xenium2SPE = function(data_dir){
    # load gene expression info ========================
    exp_file = paste0(data_dir, "/cell_feature_matrix.h5")
    exp_data <- h5read(exp_file, "matrix")
    barcodes <- as.character(h5read(exp_file, "matrix/barcodes"))
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

stagger_spatial_coords = function(spe){
    #' Stagger spatial coordinates
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