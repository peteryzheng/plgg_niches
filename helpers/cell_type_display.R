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

# Single-nucleus label -> the SAME finalized display vocabulary as `lineage_display`,
# so spatial and single-nucleus facets of a shared figure line up in the same columns.
# Keys are the annot_v0 values AFTER the plural normalization done in
# 1_PA_MAPK_scores.qmd (Astrocytes/Neurons/Lymphocytes -> Astrocyte/Neuron/Lymphoid).
# Two sn types have no exact spatial analog and are force-aligned to the closest spatial
# lineage by explicit decision (see 2_MAPK_adj.qmd): Lymphoid -> "T cells" (spatial panel
# only resolves t_cell) and Stromal -> "Endothelial" (spatial panel only resolves
# endothelial). sn has no radial-glia call, so "Radial glia" stays spatial-only.
sn_lineage_display = c(
    Astrocyte        = "Astrocytes",
    OPC              = "OPC",
    Myeloid          = "Myeloid",
    Neuron           = "Neuronal",
    Oligodendrocytes = "Oligodendrocytes",
    Lymphoid         = "T cells",
    Stromal          = "Endothelial",
    Cycling_Cells    = "Proliferating",
    undetermined     = "Undetermined"
)
