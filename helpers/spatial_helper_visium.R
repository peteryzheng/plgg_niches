suppressPackageStartupMessages({
    library(Banksy)
    library(SpatialExperiment)
    library(SummarizedExperiment)
    library(SingleCellExperiment)
    library(scuttle)
    library(scran)
    library(harmony)
    library(Matrix)
    library(data.table)
})

timestamp_message = function(text) {
    cat(sprintf("[%s] | %s\n", format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), text))
}

format_lambda_label = function(lambda) {
    as.character(as.numeric(lambda))
}

parse_numeric_csv = function(x) {
    values = trimws(unlist(strsplit(x, ",")))
    values = values[nzchar(values)]
    as.numeric(values)
}

na_vector_like = function(x, n) {
    if (is.factor(x)) {
        return(factor(rep(NA_character_, n), levels = levels(x)))
    }
    if (is.logical(x)) {
        return(rep(NA, n))
    }
    if (is.integer(x)) {
        return(rep(NA_integer_, n))
    }
    if (is.numeric(x)) {
        return(rep(NA_real_, n))
    }
    if (inherits(x, "Date")) {
        return(as.Date(rep(NA_real_, n), origin = "1970-01-01"))
    }

    rep(NA_character_, n)
}

flag_truthy = function(x) {
    if (is.logical(x)) {
        return(!is.na(x) & x)
    }

    x_chr = tolower(trimws(as.character(x)))
    !is.na(x_chr) & x_chr %in% c("true", "t", "1", "yes", "y")
}

mad_outlier_by_group = function(metric, group, threshold = "both", nmad = 3) {
    if (!threshold %in% c("both", "upper", "lower")) {
        stop("threshold must be one of 'both', 'upper', or 'lower'")
    }

    flags = rep(FALSE, length(metric))
    group = as.character(group)

    for (current_group in unique(group)) {
        idx = which(group == current_group & is.finite(metric))
        if (length(idx) < 3) {
            next
        }

        current_metric = as.numeric(metric[idx])
        current_median = median(current_metric)
        current_mad = median(abs(current_metric - current_median))
        if (is.na(current_mad) || current_mad == 0) {
            next
        }

        if (threshold == "both") {
            flags[idx] = abs(current_metric - current_median) > nmad * current_mad
        } else if (threshold == "upper") {
            flags[idx] = current_metric - current_median > nmad * current_mad
        } else {
            flags[idx] = current_metric - current_median < -nmad * current_mad
        }
    }

    flags
}

derive_sample_id = function(spe, file_path) {
    candidate_fields = c("sample.x", "sample.y")

    for (field in candidate_fields) {
        if (!field %in% colnames(colData(spe))) {
            next
        }

        values = unique(as.character(colData(spe)[[field]]))
        values = values[!is.na(values) & nzchar(values)]
        if (length(values) == 1) {
            return(values[[1]])
        }
    }

    sub("_sce$", "", tools::file_path_sans_ext(basename(file_path)))
}

coerce_to_visium_spe = function(se, file_path) {
    if (!"counts" %in% assayNames(se)) {
        stop(sprintf("counts assay missing in %s", file_path))
    }

    if (!all(c("x", "y") %in% colnames(colData(se)))) {
        stop(sprintf("x and y columns missing in colData for %s", file_path))
    }

    SpatialExperiment(
        assays = list(counts = assay(se, "counts")),
        rowData = rowData(se),
        colData = colData(se),
        spatialCoordsNames = c("x", "y")
    )
}

read_visium_spe = function(file_path) {
    timestamp_message(sprintf("Reading %s", basename(file_path)))
    se = readRDS(file_path)
    spe = coerce_to_visium_spe(se, file_path)

    sample_id = derive_sample_id(spe, file_path)
    tissue_section = if ("tissue_section" %in% colnames(colData(spe))) {
        as.character(spe$tissue_section)
    } else {
        rep("tissue_section_0", ncol(spe))
    }
    tissue_section[is.na(tissue_section) | !nzchar(tissue_section)] = "tissue_section_0"

    spe$visium_sample_id = sample_id
    spe$tissue_section = tissue_section
    spe$visium_section_id = paste(sample_id, tissue_section, sep = "__")
    spe$source_file = basename(file_path)
    spe$spot_id = colnames(spe)

    spe
}

harmonize_coldata = function(objects) {
    all_cols = unique(unlist(lapply(objects, function(spe) {
        colnames(colData(spe))
    })))

    prototypes = vector("list", length(all_cols))
    names(prototypes) = all_cols
    for (colname in all_cols) {
        for (spe in objects) {
            if (colname %in% colnames(colData(spe))) {
                prototypes[[colname]] = colData(spe)[[colname]]
                break
            }
        }
    }

    lapply(objects, function(spe) {
        current_coldata = colData(spe)
        missing_cols = setdiff(all_cols, colnames(current_coldata))
        for (colname in missing_cols) {
            current_coldata[[colname]] = na_vector_like(prototypes[[colname]], ncol(spe))
        }
        colData(spe) = current_coldata[, all_cols, drop = FALSE]
        spe
    })
}

load_visium_cohort = function(input_dir, pattern = "*_sce.rds") {
    file_paths = sort(Sys.glob(file.path(input_dir, pattern)))
    if (length(file_paths) == 0) {
        stop(sprintf("No files matched %s in %s", pattern, input_dir))
    }

    objects = lapply(file_paths, read_visium_spe)
    shared_genes = Reduce(intersect, lapply(objects, rownames))
    if (length(shared_genes) == 0) {
        stop("No shared genes across Visium input files")
    }

    shared_rowdata_cols = Reduce(intersect, lapply(objects, function(spe) {
        colnames(rowData(spe))
    }))
    objects = lapply(objects, function(spe) {
        spe = spe[shared_genes, ]
        rowData(spe) = rowData(spe)[, shared_rowdata_cols, drop = FALSE]
        spe
    })
    objects = harmonize_coldata(objects)
    load_summary = rbindlist(lapply(objects, function(spe) {
        data.table(
            source_file = unique(as.character(spe$source_file)),
            sample_id = unique(as.character(spe$visium_sample_id)),
            n_spots_input = ncol(spe),
            n_sections_input = length(unique(as.character(spe$visium_section_id)))
        )
    }))

    combined = do.call(cbind, objects)
    metadata(combined)$visium_input_files = basename(file_paths)
    metadata(combined)$visium_shared_genes = shared_genes

    list(
        spe = combined,
        load_summary = load_summary
    )
}

filter_and_normalize_visium = function(spe, aname = "normcounts", nmad = 3) {
    counts = assay(spe, "counts")
    total_counts = Matrix::colSums(counts)
    detected_genes = Matrix::colSums(counts > 0)
    sample_id = as.character(spe$visium_sample_id)

    hard_keep = rep(TRUE, ncol(spe))
    if ("in_tissue" %in% colnames(colData(spe))) {
        hard_keep = hard_keep & flag_truthy(spe$in_tissue)
    }
    if ("exclude" %in% colnames(colData(spe))) {
        hard_keep = hard_keep & !flag_truthy(spe$exclude)
    }

    nonzero_keep = total_counts > 0
    candidate_keep = hard_keep & nonzero_keep

    outlier_total = rep(FALSE, ncol(spe))
    outlier_detected = rep(FALSE, ncol(spe))
    if (any(candidate_keep)) {
        outlier_total[candidate_keep] = mad_outlier_by_group(
            log1p(total_counts[candidate_keep]),
            sample_id[candidate_keep],
            threshold = "both",
            nmad = nmad
        )
        outlier_detected[candidate_keep] = mad_outlier_by_group(
            log1p(detected_genes[candidate_keep]),
            sample_id[candidate_keep],
            threshold = "lower",
            nmad = nmad
        )
    }

    qc_keep = candidate_keep & !(outlier_total | outlier_detected)
    if (!any(qc_keep)) {
        stop("No Visium spots remain after QC")
    }

    summary_dt = data.table(
        source_file = as.character(spe$source_file),
        sample_id = sample_id,
        hard_keep = hard_keep,
        nonzero_keep = candidate_keep,
        qc_keep = qc_keep
    )[
        ,
        .(
            n_spots_input = .N,
            n_spots_hard_filter = sum(hard_keep),
            n_spots_nonzero = sum(nonzero_keep),
            n_spots_qc_keep = sum(qc_keep)
        ),
        by = .(source_file, sample_id)
    ]

    filtered = spe[, qc_keep]
    filtered = computeLibraryFactors(filtered)
    assay(filtered, aname) = normalizeCounts(filtered, log = FALSE)

    list(
        spe = filtered,
        qc_summary = summary_dt
    )
}

compute_sectionwise_banksy = function(
    spe,
    aname = "normcounts",
    target_k_geom = 18L,
    min_section_spots = 2L
) {
    split_levels = unique(as.character(spe$visium_section_id))
    split_indices = split(
        seq_len(ncol(spe)),
        factor(as.character(spe$visium_section_id), levels = split_levels)
    )

    section_objects = list()
    section_summary = vector("list", length(split_indices))

    for (i in seq_along(split_indices)) {
        section_id = names(split_indices)[i]
        idx = split_indices[[i]]
        section = spe[, idx]
        n_spots = ncol(section)
        sample_id = unique(as.character(section$visium_sample_id))
        tissue_section = unique(as.character(section$tissue_section))
        source_file = unique(as.character(section$source_file))

        if (n_spots < min_section_spots) {
            section_summary[[i]] = data.table(
                source_file = source_file,
                sample_id = sample_id,
                tissue_section = tissue_section,
                sample_section_id = section_id,
                n_spots = n_spots,
                k_geom_used = NA_integer_,
                status = "dropped_too_small"
            )
            next
        }

        k_geom_used = min(as.integer(target_k_geom), n_spots - 1L)
        timestamp_message(sprintf(
            "Computing BANKSY for %s (%i spots, k_geom=%i)",
            section_id, n_spots, k_geom_used
        ))

        section_objects[[length(section_objects) + 1]] = computeBanksy(
            section,
            assay_name = aname,
            compute_agf = TRUE,
            k_geom = k_geom_used,
            verbose = FALSE
        )

        section_summary[[i]] = data.table(
            source_file = source_file,
            sample_id = sample_id,
            tissue_section = tissue_section,
            sample_section_id = section_id,
            n_spots = n_spots,
            k_geom_used = k_geom_used,
            status = "kept"
        )
    }

    if (length(section_objects) == 0) {
        stop("All Visium sections were dropped before BANKSY")
    }

    combined = suppressWarnings(do.call(cbind, section_objects))
    combined$visium_sample_id = unlist(lapply(section_objects, function(x) {
        as.character(x$visium_sample_id)
    }), use.names = FALSE)
    combined$tissue_section = unlist(lapply(section_objects, function(x) {
        as.character(x$tissue_section)
    }), use.names = FALSE)
    combined$visium_section_id = unlist(lapply(section_objects, function(x) {
        as.character(x$visium_section_id)
    }), use.names = FALSE)
    combined$source_file = unlist(lapply(section_objects, function(x) {
        as.character(x$source_file)
    }), use.names = FALSE)
    metadata(combined)$BANKSY_params = list(
        assay_name = aname,
        M = c(0L, 1L),
        k_geom = as.numeric(target_k_geom),
        spatial_mode = "kNN_median"
    )

    list(
        spe = combined,
        section_summary = rbindlist(section_summary, fill = TRUE)
    )
}

run_visium_banksy = function(
    spe,
    lambda_vec = c(0, 0.2),
    npcs = 20L,
    k_neighbors = 30L,
    resolution = 0.8,
    seed = 55555L,
    harmony_vars = c("visium_sample_id")
) {
    timestamp_message("Running BANKSY PCA")
    spe = runBanksyPCA(
        spe,
        assay_name = "normcounts",
        use_agf = TRUE,
        lambda = lambda_vec,
        npcs = npcs,
        seed = seed
    )

    meta_df = as.data.frame(colData(spe))

    for (lambda in lambda_vec) {
        lambda_label = format_lambda_label(lambda)
        pca_name = paste0("PCA_M1_lam", lambda_label)
        harmony_name = paste0("Harmony_BANKSY_lam", lambda_label)

        timestamp_message(sprintf("Running Harmony for lambda=%s", lambda_label))
        if (all(harmony_vars %in% colnames(meta_df)) &&
            length(unique(meta_df[[harmony_vars[[1]]]])) > 1) {
            reducedDim(spe, harmony_name) = RunHarmony(
                data_mat = reducedDim(spe, pca_name),
                meta_data = meta_df,
                vars_use = harmony_vars,
                do_pca = FALSE,
                max_iter = 50,
                verbose = TRUE
            )
        } else {
            reducedDim(spe, harmony_name) = reducedDim(spe, pca_name)
        }

        timestamp_message(sprintf("Running UMAP for lambda=%s", lambda_label))
        spe = runBanksyUMAP(
            spe,
            dimred = harmony_name,
            seed = seed
        )

        timestamp_message(sprintf("Running Leiden clustering for lambda=%s", lambda_label))
        spe = clusterBanksy(
            spe,
            dimred = harmony_name,
            algo = "leiden",
            k_neighbors = k_neighbors,
            resolution = resolution,
            seed = seed
        )
    }

    spe
}

find_lambda_cluster_column = function(spe, lambda) {
    lambda_label = gsub("\\.", "\\\\.", format_lambda_label(lambda))
    candidates = grep(
        paste0("^clust_Harmony_BANKSY_lam", lambda_label),
        clusterNames(spe),
        value = TRUE
    )
    if (length(candidates) == 0) {
        stop(sprintf("No cluster column found for lambda=%s", format_lambda_label(lambda)))
    }
    candidates[[1]]
}

resolve_visium_cluster_column = function(spe, lambda = NULL, preferred_col = NULL) {
    if (!is.null(preferred_col) &&
        nzchar(preferred_col) &&
        preferred_col %in% colnames(colData(spe))) {
        return(preferred_col)
    }

    if (is.null(lambda)) {
        stop("lambda is required when preferred_col is unavailable")
    }

    find_lambda_cluster_column(spe, lambda)
}

sort_visium_values = function(values) {
    values = unique(as.character(values))
    values = values[!is.na(values) & nzchar(values)]
    if (length(values) == 0) {
        return(character())
    }

    numeric_values = suppressWarnings(as.numeric(values))
    if (all(!is.na(numeric_values))) {
        return(values[order(numeric_values)])
    }

    sort(values)
}

build_visium_palette = function(values) {
    values = sort_visium_values(values)
    if (length(values) == 0) {
        return(setNames(character(), character()))
    }

    hues = seq(15, 375, length.out = length(values) + 1)
    colors = grDevices::hcl(
        h = hues[seq_len(length(values))],
        l = 65,
        c = 100
    )
    setNames(colors, values)
}

visium_sample_order = function(spe) {
    data.table(
        sample_id = as.character(spe$visium_sample_id)
    )[
        ,
        .N,
        by = sample_id
    ][
        order(-N, sample_id)
    ]$sample_id
}

visium_spatial_plot_dt = function(spe, cluster_col, sample_order = NULL) {
    coords = as.data.table(spatialCoords(spe))
    if (ncol(coords) < 2) {
        stop("spatialCoords must contain at least two dimensions")
    }

    coords = coords[, seq_len(2), with = FALSE]
    setnames(coords, colnames(coords), c("x", "y"))

    plot_dt = copy(coords)
    plot_dt[, sample_id := as.character(colData(spe)$visium_sample_id)]
    plot_dt[, tissue_section := as.character(colData(spe)$tissue_section)]
    plot_dt[, cluster := as.character(colData(spe)[[cluster_col]])]

    if (is.null(sample_order)) {
        sample_order = visium_sample_order(spe)
    }

    plot_dt[, sample_id := factor(sample_id, levels = sample_order)]
    plot_dt[, cluster := factor(cluster, levels = sort_visium_values(cluster))]
    plot_dt
}

visium_cluster_composition_dt = function(spe, cluster_col, sample_order = NULL) {
    if (is.null(sample_order)) {
        sample_order = visium_sample_order(spe)
    }

    composition = data.table(
        sample_id = as.character(spe$visium_sample_id),
        cluster = as.character(colData(spe)[[cluster_col]])
    )[
        ,
        .N,
        by = .(sample_id, cluster)
    ][
        ,
        fraction := N / sum(N),
        by = sample_id
    ]

    composition[, sample_id := factor(sample_id, levels = sample_order)]
    composition[, cluster := factor(cluster, levels = sort_visium_values(cluster))]
    composition[]
}

plot_visium_umap = function(
    spe,
    dimred_name,
    color_var,
    title,
    point_size = 0.3,
    alpha = 0.8,
    palette = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plotting")
    }

    embedding = as.data.table(reducedDim(spe, dimred_name))
    if (ncol(embedding) < 2) {
        stop(sprintf("%s must contain at least two dimensions", dimred_name))
    }

    embedding = embedding[, seq_len(2), with = FALSE]
    setnames(embedding, colnames(embedding), c("umap_1", "umap_2"))
    embedding[, value := factor(
        as.character(colData(spe)[[color_var]]),
        levels = sort_visium_values(colData(spe)[[color_var]])
    )]

    if (is.null(palette)) {
        palette = build_visium_palette(levels(embedding$value))
    }

    ggplot2::ggplot(embedding, ggplot2::aes(umap_1, umap_2, color = value)) +
        ggplot2::geom_point(size = point_size, alpha = alpha) +
        ggplot2::scale_color_manual(
            values = palette,
            drop = FALSE,
            na.value = "grey80"
        ) +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title = title,
            color = color_var,
            x = "UMAP 1",
            y = "UMAP 2"
        )
}

plot_visium_spatial_sample_grid = function(
    spe,
    cluster_col,
    title = NULL,
    sample_order = NULL,
    point_size = 0.35,
    alpha = 0.85,
    ncol = 3L,
    palette = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plotting")
    }
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("cowplot is required for plotting")
    }

    plot_dt = visium_spatial_plot_dt(
        spe,
        cluster_col = cluster_col,
        sample_order = sample_order
    )
    sample_order = levels(plot_dt$sample_id)

    if (is.null(palette)) {
        palette = build_visium_palette(levels(plot_dt$cluster))
    }

    sample_plots = lapply(sample_order, function(current_sample) {
        current_dt = plot_dt[sample_id == current_sample]

        ggplot2::ggplot(current_dt, ggplot2::aes(x, y, color = cluster)) +
            ggplot2::geom_point(size = point_size, alpha = alpha) +
            ggplot2::coord_equal() +
            ggplot2::scale_color_manual(
                values = palette,
                drop = FALSE,
                na.value = "grey80"
            ) +
            ggplot2::theme_bw() +
            ggplot2::labs(
                title = sprintf(
                    "%s (n=%s)",
                    current_sample,
                    format(nrow(current_dt), big.mark = ",")
                ),
                x = NULL,
                y = NULL
            ) +
            ggplot2::theme(
                legend.position = "none",
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(hjust = 0.5, size = 10)
            )
    })

    legend = cowplot::get_legend(
        ggplot2::ggplot(plot_dt, ggplot2::aes(x, y, color = cluster)) +
            ggplot2::geom_point(size = 1.5, alpha = 1) +
            ggplot2::scale_color_manual(
                values = palette,
                drop = FALSE,
                na.value = "grey80"
            ) +
            ggplot2::theme_void() +
            ggplot2::theme(legend.position = "bottom") +
            ggplot2::guides(
                color = ggplot2::guide_legend(
                    override.aes = list(size = 3, alpha = 1),
                    nrow = max(1L, ceiling(length(palette) / 8))
                )
            )
    )

    combined_plot = cowplot::plot_grid(
        plotlist = sample_plots,
        ncol = as.integer(ncol)
    )
    combined_plot = cowplot::plot_grid(
        combined_plot,
        legend,
        ncol = 1,
        rel_heights = c(1, 0.14)
    )

    if (!is.null(title) && nzchar(title)) {
        combined_plot = cowplot::plot_grid(
            cowplot::ggdraw() +
                cowplot::draw_label(title, fontface = "bold"),
            combined_plot,
            ncol = 1,
            rel_heights = c(0.08, 1)
        )
    }

    combined_plot
}

find_visium_markers = function(spe, lambda = 0.2) {
    cluster_col = find_lambda_cluster_column(spe, lambda)
    groups = as.factor(colData(spe)[[cluster_col]])
    if (length(unique(groups)) < 2) {
        warning(sprintf("Skipping markers for %s because only one cluster is present", cluster_col))
        return(NULL)
    }

    list(
        cluster_column = cluster_col,
        markers = findMarkers(
            assay(spe, "counts"),
            groups = groups,
            test.type = "wilcox"
        )
    )
}

write_tsv = function(x, output_path) {
    fwrite(as.data.table(x), output_path, sep = "\t")
}

render_visium_report = function(qmd_path, run_dir) {
    quarto_candidates = c("/usr/local/bin/quarto", Sys.which("quarto"))
    quarto_candidates = unique(quarto_candidates[nzchar(quarto_candidates)])
    existing_candidates = quarto_candidates[file.exists(quarto_candidates)]
    quarto_bin = if (length(existing_candidates) > 0) {
        existing_candidates[[1]]
    } else {
        ""
    }
    if (!nzchar(quarto_bin)) {
        warning("quarto not found on PATH; skipping report render")
        return(FALSE)
    }

    render_args = c(
        "render",
        qmd_path,
        "--output-dir", run_dir,
        "-P", paste0("run_dir:", normalizePath(run_dir))
    )
    status = system2(quarto_bin, render_args)
    identical(status, 0L)
}
