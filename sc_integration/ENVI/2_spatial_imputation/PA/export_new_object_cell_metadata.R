#!/usr/bin/env Rscript
# Phase 0 -- Export per-cell metadata from the finalized (res_ct=0.2 / res_ni=0.5)
# pathology-annotated BANKSY object as a single source of truth so that the ENVI
# re-scoring (1_PA_MAPK_scores.qmd) and the R niche analysis (2_MAPK_adj.qmd) label
# every cell identically. Reads the FULL object (heavy -- run on a compute node).
#
# Cell type and niche are read DIRECTLY from the object's own finalized colData
# columns. Do NOT join the external banksy_proseg_cell_cluster_annotation.tsv: it
# is stale and its cluster numbers do not map onto this object (joining it
# scrambles every label). See memory reference_new_banksy_object_columns.
#
# Output (large; gitignored): data/xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv

suppressMessages({library(SummarizedExperiment); library(SpatialExperiment); library(data.table)})

# workdir resolves per-machine from HOME (argos = /mnt/storage/dept/medonc/beroukhim/)
home = path.expand('~')
workdir = if (home %in% c('/Users/youyun','/Users/youyunzheng')) file.path(home,'Documents/HMS/PhD/beroukhimlab/dfci_mount') else
          if (home == '/home/yz762') '/mnt/storage/dept/medonc/beroukhim' else
          if (home == '/PHShome/yz762') '/data/beroukhim1' else '/xchip/beroukhimlab'
workdir = paste0(workdir, '/')

obj_path = paste0(workdir, 'youyun/plgg/data/Xenium_annotations/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134_pathology.rds')
out_path = paste0(workdir, 'youyun/plgg/data/xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv')

cat('Reading object:', obj_path, '\n'); flush.console()
o = readRDS(obj_path)
cd = colData(o)
cat('Loaded colData:', nrow(cd), 'cells x', ncol(cd), 'cols\n'); flush.console()

# Finalized labels live in the object's own colData -- read directly (no join).
ct_clust_col = 'clust_Harmony_BANKSY_lam0.2_k50_res0.2'            # raw cell-type cluster id (13 clusters)
ct_type_col  = 'cell_type_clust_Harmony_BANKSY_lam0.2_k50_res0.2' # finalized cell type (10 labels, incl 'myeloid')
ni_col       = 'clust_Harmony_BANKSY_lam0.8_k50_res0.5'           # niche cluster id (17 niches)
stopifnot(all(c(ct_clust_col, ct_type_col, ni_col, 'sample_id', 'histology') %in% colnames(cd)))
alt_col  = if ('alteration' %in% colnames(cd)) 'alteration' else NA_character_
has_path = 'pathology_annotation' %in% colnames(cd)   # STalign pathology-region category (optional)
has_reg  = 'region_id' %in% colnames(cd)              # STalign region UUID (optional)

# Spatial (x, y) centroids for the Phase-2 kNN neighborhood graph. Exporting them
# here makes 2_MAPK_adj.qmd object-free: it can rebuild the spatial graph from this
# TSV instead of re-deserializing the 5 GB object each render. Prefer the
# SpatialExperiment spatialCoords slot; fall back to colData centroid columns.
xy = tryCatch(SpatialExperiment::spatialCoords(o), error = function(e) NULL)
if (!is.null(xy) && ncol(xy) >= 2) {
  x_coord = as.numeric(xy[, 1]); y_coord = as.numeric(xy[, 2])
  cat('Coordinates from spatialCoords():', paste(colnames(xy)[1:2], collapse = ','), '\n')
} else {
  # fall back to conventionally-named centroid columns in colData
  xcol = intersect(c('x_centroid','x','sdimx','CellCenterX','x_global_px','x_location'), colnames(cd))[1]
  ycol = intersect(c('y_centroid','y','sdimy','CellCenterY','y_global_px','y_location'), colnames(cd))[1]
  stopifnot(!is.na(xcol), !is.na(ycol))
  x_coord = as.numeric(cd[[xcol]]); y_coord = as.numeric(cd[[ycol]])
  cat('Coordinates from colData columns:', xcol, ',', ycol, '\n')
}

meta = data.table(
  cell_id    = rownames(cd),                          # 'barcode-1__sampleID' -- matches the ENVI obs_names
  sample_id  = as.character(cd$sample_id),
  histology  = as.character(cd$histology),
  alteration = if (!is.na(alt_col)) as.character(cd[[alt_col]]) else NA_character_,
  ct_cluster = as.character(cd[[ct_clust_col]]),      # raw cell-type cluster id (kept for reference)
  cell_type  = as.character(cd[[ct_type_col]]),       # finalized cell type, read straight from the object
  niche      = as.character(cd[[ni_col]]),            # niche id (1-17)
  x          = x_coord,                               # spatial centroid x (per-sample local px)
  y          = y_coord,                               # spatial centroid y (per-sample local px)
  pathology_annotation = if (has_path) as.character(cd$pathology_annotation) else NA_character_,
  region_id  = if (has_reg)  as.character(cd$region_id) else NA_character_
)
setcolorder(meta, c('cell_id','sample_id','histology','alteration','ct_cluster',
                    'cell_type','niche','x','y','pathology_annotation','region_id'))

fwrite(meta, out_path, sep = '\t')
cat('Wrote', nrow(meta), 'rows ->', out_path, '\n')
cat('Cell-type counts (whole cohort):\n'); print(meta[, .N, cell_type][order(-N)])
cat('NA cell_type cells:', meta[is.na(cell_type), .N], '\n')

# Sanity check the thing that was previously broken: PA-only cell-type fractions.
pa_ids = c('188570','258482','320784','328368')
pa = meta[sample_id %in% pa_ids]
cat('PA-only cell-type counts (n =', nrow(pa), 'cells):\n'); print(pa[, .N, cell_type][order(-N)])
cat('PA myeloid fraction:', round(mean(pa$cell_type == 'myeloid', na.rm=TRUE), 4), '\n')
cat('Top niches by count:\n'); print(head(meta[, .N, niche][order(-N)], 20))
