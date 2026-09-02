# Display labels for the automated cell-type (lineage) calls.
#
# The pipeline's cell-type keys are the names of `marker_panel` in
# helpers/spatial_helper.R (astro_glial, radial_glia, ...), assigned by
# score_lineages_top_gene(); "undetermined" is added when the top marker AUC
# is below threshold. Those snake_case keys are machine identifiers reused
# verbatim across many notebooks/TSVs, so we do NOT rename them at the source --
# this is a display-only lookup, applied at plot time.
#
# Kept in its own lightweight file (no library() calls) so lean notebooks can
# source it without pulling in spatial_helper.R's heavy stack (Seurat/Banksy/
# harmony/ComplexHeatmap). If marker_panel gains a lineage, add it here too.
lineage_display = c(
    astro_glial   = "Astrocytes",
    radial_glia   = "Radial glia",
    opc           = "OPC",
    oligo_mature  = "Oligodendrocytes",
    neuronal      = "Neuronal",
    myeloid       = "Myeloid",
    t_cell        = "T cells",
    endothelial   = "Endothelial",
    proliferation = "Proliferating",
    undetermined  = "Undetermined"
)
