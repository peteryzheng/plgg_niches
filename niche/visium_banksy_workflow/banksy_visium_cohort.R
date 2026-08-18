library(optparse)

get_script_path = function() {
    file_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) > 0) {
        return(normalizePath(sub("^--file=", "", file_arg[[1]])))
    }

    normalizePath(getwd())
}

script_path = get_script_path()
script_dir = dirname(script_path)
code_root = normalizePath(file.path(script_dir, "..", ".."))
source(file.path(code_root, "helpers", "spatial_helper_visium.R"))

option_list = list(
    make_option(
        c("-i", "--inputdir"),
        type = "character",
        help = "Directory containing Visium *_sce.rds files"
    ),
    make_option(
        c("-o", "--outputdir"),
        type = "character",
        help = "Directory where BANKSY outputs should be written"
    ),
    make_option(
        c("-p", "--pattern"),
        type = "character",
        default = "*_sce.rds",
        help = "Glob pattern used to match Visium input files [default %default]"
    ),
    make_option(
        c("--k_geom"),
        type = "integer",
        default = 18,
        help = "Target BANKSY k_geom for Visium sections [default %default]"
    ),
    make_option(
        c("--npc"),
        type = "integer",
        default = 20,
        help = "Number of BANKSY PCs [default %default]"
    ),
    make_option(
        c("--lambdas"),
        type = "character",
        default = "0,0.2",
        help = "Comma-separated lambda values [default %default]"
    ),
    make_option(
        c("--k_neighbors"),
        type = "integer",
        default = 30,
        help = "Leiden neighbor count [default %default]"
    ),
    make_option(
        c("--resolution"),
        type = "double",
        default = 0.8,
        help = "Leiden resolution [default %default]"
    ),
    make_option(
        c("--seed"),
        type = "integer",
        default = 55555,
        help = "Random seed [default %default]"
    ),
    make_option(
        c("--skip_report"),
        action = "store_true",
        default = FALSE,
        help = "Skip Quarto QC report rendering"
    )
)

opt = parse_args(OptionParser(option_list = option_list))
if (is.null(opt$inputdir) || is.null(opt$outputdir)) {
    stop("--inputdir and --outputdir are required")
}

lambda_vec = parse_numeric_csv(opt$lambdas)
if (any(is.na(lambda_vec))) {
    stop("Unable to parse --lambdas")
}

current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir = file.path(opt$outputdir, paste0("visium_banksy_", current_timestamp))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

timestamp_message("Loading Visium cohort")
load_result = load_visium_cohort(opt$inputdir, pattern = opt$pattern)

timestamp_message("Filtering and normalizing Visium spots")
filter_result = filter_and_normalize_visium(load_result$spe, aname = "normcounts")

timestamp_message("Computing section-wise BANKSY matrices")
banksy_input = compute_sectionwise_banksy(
    filter_result$spe,
    aname = "normcounts",
    target_k_geom = opt$k_geom
)

timestamp_message("Running Visium BANKSY embedding and clustering")
banksy_output = run_visium_banksy(
    banksy_input$spe,
    lambda_vec = lambda_vec,
    npcs = opt$npc,
    k_neighbors = opt$k_neighbors,
    resolution = opt$resolution,
    seed = opt$seed
)

timestamp_message("Finding lambda 0.2 markers")
marker_result = find_visium_markers(banksy_output, lambda = 0.2)

output_rds = file.path(run_dir, "visium_banksy_output.rds")
markers_rds = file.path(run_dir, "visium_banksy_markers_lam0.2.rds")
qc_summary_tsv = file.path(run_dir, "visium_banksy_qc_summary.tsv")
section_summary_tsv = file.path(run_dir, "visium_banksy_section_summary.tsv")
run_summary_tsv = file.path(run_dir, "visium_banksy_run_summary.tsv")

saveRDS(banksy_output, output_rds)
saveRDS(marker_result, markers_rds)
write_tsv(filter_result$qc_summary, qc_summary_tsv)
write_tsv(banksy_input$section_summary, section_summary_tsv)

run_summary = data.table(
    run_dir = normalizePath(run_dir),
    input_dir = normalizePath(opt$inputdir),
    pattern = opt$pattern,
    n_files = nrow(load_result$load_summary),
    n_samples = length(unique(as.character(banksy_output$visium_sample_id))),
    n_spots_final = ncol(banksy_output),
    n_genes_shared = length(metadata(load_result$spe)$visium_shared_genes),
    n_sections_kept = sum(banksy_input$section_summary$status == "kept"),
    n_sections_dropped = sum(banksy_input$section_summary$status != "kept"),
    k_geom_requested = opt$k_geom,
    lambdas = paste(format_lambda_label(lambda_vec), collapse = ","),
    npc = opt$npc,
    k_neighbors = opt$k_neighbors,
    resolution = opt$resolution,
    seed = opt$seed,
    marker_cluster_column = if (is.null(marker_result)) NA_character_ else marker_result$cluster_column
)
write_tsv(run_summary, run_summary_tsv)

report_rendered = FALSE
if (!opt$skip_report) {
    timestamp_message("Rendering Visium BANKSY QC report")
    report_rendered = render_visium_report(
        file.path(script_dir, "visium_banksy_qc.qmd"),
        run_dir
    )
}

timestamp_message(sprintf("Done. Run directory: %s", normalizePath(run_dir)))
timestamp_message(sprintf("QC report rendered: %s", report_rendered))
