library(sceasy)
library(Seurat)

use_condaenv(
    condaenv = 'rameen', 
    conda = '/xchip/beroukhimlab/youyun/miniconda3/condabin/conda',
    required = TRUE
)
loompy <- reticulate::import('loompy')

# local vs UGER
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}

# load the data 
gnt_filepath = paste0(workdir, "coja/Spatial_PLGG/data/Seurat_Object/harmony_lambda10_t2.rds")
gnt_sc = readRDS(gnt_filepath)
gnt_sc[["RNA"]] <- as(gnt_sc[["RNA"]], "Assay")

print('converting to anndata')
sceasy::convertFormat(
    gnt_sc, from="seurat", to="anndata",
    outFile=gsub('.rds$', '_sceasy.h5ad', gnt_filepath)
)
print(paste0('Wrote file to ', gsub('.rds$', '_sceasy.h5ad', gnt_filepath)))
