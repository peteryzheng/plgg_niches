# Shared mixed-model utilities for pathology-region enrichment analyses.

build_feature_count_dt = function(
    coldata_dt,
    feature_col,
    sample_col = "sample_id",
    region_col = "pathology_annotation",
    observation_col = region_col,
    extra_cols = c("histology", "alteration")
) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("The data.table package is required.")
    }

    required_cols = unique(c(sample_col, region_col, observation_col, feature_col, extra_cols))
    missing_cols = setdiff(required_cols, colnames(coldata_dt))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }

    dt = data.table::as.data.table(coldata_dt)[
        , ..required_cols
    ]

    dt[
        , `:=`(
            sample_id = as.character(get(sample_col)),
            observation_id = as.character(get(observation_col)),
            pathology_annotation = as.character(get(region_col)),
            feature = data.table::fifelse(
                is.na(get(feature_col)) | trimws(as.character(get(feature_col))) == "",
                "Other",
                as.character(get(feature_col))
            )
        )
    ]

    observation_meta = dt[
        , c(list(n = .N), lapply(.SD, function(x) x[1])),
        by = .(sample_id, observation_id, pathology_annotation),
        .SDcols = extra_cols
    ]

    feature_counts = dt[
        , .(k = .N),
        by = .(sample_id, observation_id, pathology_annotation, feature)
    ]

    all_features = sort(unique(dt$feature))

    expanded = observation_meta[
        , .(feature = all_features),
        by = c("sample_id", "observation_id", "pathology_annotation", extra_cols, "n")
    ]

    out_dt = feature_counts[
        expanded,
        on = .(sample_id, observation_id, pathology_annotation, feature)
    ]

    out_dt[
        is.na(k),
        k := 0L
    ]

    out_dt[
        , `:=`(
            k = as.integer(k),
            n = as.integer(n)
        )
    ]

    data.table::setcolorder(
        out_dt,
        c("sample_id", "observation_id", "pathology_annotation", extra_cols, "feature", "k", "n")
    )

    return(out_dt[])
}

classify_model_status = function(fit_obj, warning_log, coef_dt, expect_random_effect = TRUE) {
    reason_codes = c()
    fatal_reasons = c()
    random_effect_sd = NA_real_

    if (is.null(fit_obj)) {
        return(list(
            model_status = "fail",
            status_reason = "fit_error",
            random_effect_sd = random_effect_sd
        ))
    }

    if (is.null(coef_dt) || nrow(coef_dt) == 0) {
        return(list(
            model_status = "fail",
            status_reason = "no_coefficients",
            random_effect_sd = random_effect_sd
        ))
    }

    if (any(!is.finite(coef_dt$estimate_log_odds)) || any(!is.finite(coef_dt$se))) {
        fatal_reasons = c(fatal_reasons, "non_finite_coefficients")
    }

    conv_code = tryCatch(fit_obj$fit$convergence, error = function(e) NA_integer_)
    if (!is.na(conv_code) && conv_code != 0) {
        fatal_reasons = c(fatal_reasons, "nonzero_convergence_code")
    }

    pd_hess = tryCatch(isTRUE(fit_obj$sdr$pdHess), error = function(e) NA)
    if (!isTRUE(pd_hess)) {
        fatal_reasons = c(fatal_reasons, "non_pd_hessian")
    }

    # Only assess the random-effect SD when a random intercept was actually fit.
    # For between-sample groupings (histology) we deliberately omit the RI, so a
    # missing/near-zero RE SD is expected, not a defect -- don't flag it.
    if (expect_random_effect) {
        random_effect_sd = tryCatch({
            vc = glmmTMB::VarCorr(fit_obj)$cond
            if (length(vc) == 0) {
                NA_real_
            } else {
                first_group = vc[[1]]
                sd_vec = attr(first_group, "stddev")
                if (is.null(sd_vec) || length(sd_vec) == 0) NA_real_ else as.numeric(sd_vec[1])
            }
        }, error = function(e) NA_real_)

        if (is.na(random_effect_sd)) {
            reason_codes = c(reason_codes, "random_effect_sd_na")
        } else if (random_effect_sd < 1e-6) {
            reason_codes = c(reason_codes, "near_zero_random_effect_sd")
        }
    }

    if (length(warning_log) > 0) {
        reason_codes = c(reason_codes, "warnings_present")

        fatal_warning_patterns = c(
            "failed to converge",
            "false convergence",
            "non-positive-definite",
            "non positive definite",
            "hessian",
            "nan",
            "na/nan"
        )
        if (any(vapply(
            fatal_warning_patterns,
            function(pat) any(grepl(pat, warning_log, ignore.case = TRUE)),
            logical(1)
        ))) {
            fatal_reasons = c(fatal_reasons, "fatal_warning_pattern")
        }
    }

    all_reasons = unique(c(fatal_reasons, reason_codes))

    model_status = if (length(fatal_reasons) > 0) {
        "fail"
    } else if (length(all_reasons) > 0) {
        "warn"
    } else {
        "pass"
    }

    status_reason = if (length(all_reasons) == 0) "none" else paste(all_reasons, collapse = ";")

    return(list(
        model_status = model_status,
        status_reason = status_reason,
        random_effect_sd = random_effect_sd
    ))
}

fit_feature_betabinom = function(
    count_dt,
    feature_name,
    ref_region,
    sample_col = "sample_id",
    region_col = "pathology_annotation",
    fixed_effect_col = region_col,      # column supplying the fixed effect (region / histology / indicator)
    random_effect_col = sample_col      # grouping for (1 | .); NULL omits the random intercept entirely
) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop("The glmmTMB package is required for beta-binomial mixed effects models.")
    }
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("The data.table package is required.")
    }

    dt = data.table::as.data.table(count_dt)[feature == feature_name]
    warning_log = c()
    fit_error = NA_character_

    if (nrow(dt) == 0) {
        return(list(
            coef_dt = data.table::data.table(),
            model_status = "fail",
            status_reason = "empty_feature_subset",
            warning_log = "",
            fit_error = "empty_feature_subset",
            random_effect_sd = NA_real_
        ))
    }

    # Fixed effect is always carried internally as `pathology_annotation` so the
    # formula and term-stripping below stay identical regardless of which source
    # column supplies it (a real region, histology, or a binary indicator).
    model_dt = dt[
        , .(
            k = as.integer(k),
            n = as.integer(n),
            pathology_annotation = as.factor(get(fixed_effect_col))
        )
    ]
    # Random-intercept grouping is optional: NULL drops the (1 | .) term (used for
    # between-sample groupings like histology, one row per sample, where a sample
    # RI is unidentifiable -- the beta-binomial dispersion carries that variance).
    if (!is.null(random_effect_col)) {
        model_dt[, sample_id := as.factor(dt[[random_effect_col]])]
    }

    if (ref_region %in% levels(model_dt$pathology_annotation)) {
        model_dt[
            , pathology_annotation := stats::relevel(pathology_annotation, ref = ref_region)
        ]
    }

    # Build the RHS dynamically; defaults reproduce the original literal
    # `pathology_annotation + (1 | sample_id)` byte-for-byte (pathology unchanged).
    formula_rhs = if (!is.null(random_effect_col)) {
        "pathology_annotation + (1 | sample_id)"
    } else {
        "pathology_annotation"
    }
    model_formula = stats::as.formula(paste0("cbind(k, n - k) ~ ", formula_rhs))

    fit_obj = withCallingHandlers(
        tryCatch(
            glmmTMB::glmmTMB(
                model_formula,
                data = model_dt,
                family = glmmTMB::betabinomial(link = "logit")
            ),
            error = function(e) {
                fit_error <<- conditionMessage(e)
                NULL
            }
        ),
        warning = function(w) {
            warning_log <<- c(warning_log, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    if (is.null(fit_obj)) {
        return(list(
            coef_dt = data.table::data.table(),
            model_status = "fail",
            status_reason = "fit_error",
            warning_log = paste(unique(warning_log), collapse = " | "),
            fit_error = fit_error,
            random_effect_sd = NA_real_
        ))
    }

    coef_mat = tryCatch(
        summary(fit_obj)$coefficients$cond,
        error = function(e) NULL
    )

    if (is.null(coef_mat) || nrow(coef_mat) == 0) {
        return(list(
            coef_dt = data.table::data.table(),
            model_status = "fail",
            status_reason = "no_coefficients",
            warning_log = paste(unique(warning_log), collapse = " | "),
            fit_error = fit_error,
            random_effect_sd = NA_real_
        ))
    }

    coef_dt = data.table::as.data.table(coef_mat, keep.rownames = "term")

    estimate_col = "Estimate"
    se_col = names(coef_dt)[grepl("^Std", names(coef_dt))][1]
    stat_col = names(coef_dt)[grepl(" value$", names(coef_dt))][1]
    p_col = names(coef_dt)[grepl("^Pr\\(", names(coef_dt))][1]

    coef_dt[
        , `:=`(
            estimate_log_odds = get(estimate_col),
            se = get(se_col),
            z_or_t = get(stat_col),
            p_value = get(p_col)
        )
    ]

    coef_dt[
        , annotation_region := data.table::fifelse(
            term == "(Intercept)",
            "(Intercept)",
            sub("^pathology_annotation", "", term)
        )
    ]

    coef_dt[
        , `:=`(
            odds_ratio = exp(estimate_log_odds),
            ci_low = exp(estimate_log_odds - 1.96 * se),
            ci_high = exp(estimate_log_odds + 1.96 * se)
        )
    ]

    coef_dt = coef_dt[
        , .(
            term,
            annotation_region,
            estimate_log_odds,
            se,
            z_or_t,
            p_value,
            odds_ratio,
            ci_low,
            ci_high
        )
    ]

    status_obj = classify_model_status(
        fit_obj, warning_log, coef_dt,
        expect_random_effect = !is.null(random_effect_col)
    )

    return(list(
        coef_dt = coef_dt,
        model_status = status_obj$model_status,
        status_reason = status_obj$status_reason,
        warning_log = paste(unique(warning_log), collapse = " | "),
        fit_error = fit_error,
        random_effect_sd = status_obj$random_effect_sd
    ))
}

fit_feature_set = function(
    count_dt,
    feature_col,
    histology_value,
    ref_region,
    filters = list(),
    fdr_scope = c("histology_modality")
) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("The data.table package is required.")
    }

    fdr_scope = match.arg(fdr_scope)
    dt = data.table::as.data.table(data.table::copy(count_dt))

    if (!"feature" %in% names(dt)) {
        if (feature_col %in% names(dt)) {
            data.table::setnames(dt, feature_col, "feature")
        } else {
            stop(sprintf("Feature column '%s' not found.", feature_col))
        }
    }

    required_cols = c("sample_id", "pathology_annotation", "histology", "feature", "k", "n")
    missing_cols = setdiff(required_cols, names(dt))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns in count_dt: %s", paste(missing_cols, collapse = ", ")))
    }

    dt = dt[histology == histology_value]

    if (!is.null(filters$sample_ids_include)) {
        dt = dt[sample_id %in% as.character(filters$sample_ids_include)]
    }
    if (!is.null(filters$sample_ids_exclude)) {
        dt = dt[!sample_id %in% as.character(filters$sample_ids_exclude)]
    }
    if (!is.null(filters$pathology_exclude_regex) && !is.na(filters$pathology_exclude_regex)) {
        dt = dt[!grepl(filters$pathology_exclude_regex, pathology_annotation)]
    }
    if (!is.null(filters$pathology_include_regex) && !is.na(filters$pathology_include_regex)) {
        dt = dt[grepl(filters$pathology_include_regex, pathology_annotation)]
    }

    modality = if (!is.null(filters$modality_label)) {
        as.character(filters$modality_label)
    } else {
        as.character(feature_col)
    }

    filter_log = data.table::data.table(
        histology = histology_value,
        modality = modality,
        stage = character(),
        detail = character(),
        n_removed = integer()
    )

    row_total_before = nrow(dt)
    dt = dt[n >= 30]
    filter_log = rbind(
        filter_log,
        data.table::data.table(
            histology = histology_value,
            modality = modality,
            stage = "row_filter_n",
            detail = "Removed rows with n < 30",
            n_removed = row_total_before - nrow(dt)
        )
    )

    unique_regions_before = unique(dt$pathology_annotation)
    region_counts = unique(dt[, .(sample_id, pathology_annotation)])[,
        .(sample_region_rows = .N),
        by = pathology_annotation
    ]
    kept_regions = region_counts[sample_region_rows >= 2, pathology_annotation]
    dropped_regions = setdiff(unique_regions_before, kept_regions)
    dt = dt[pathology_annotation %in% kept_regions]

    filter_log = rbind(
        filter_log,
        data.table::data.table(
            histology = histology_value,
            modality = modality,
            stage = "pathology_level_filter",
            detail = ifelse(
                length(dropped_regions) == 0,
                "No pathology levels removed",
                paste("Dropped:", paste(dropped_regions, collapse = ", "))
            ),
            n_removed = length(dropped_regions)
        )
    )

    feature_qc = dt[
        , .(
            total_k = sum(k),
            nonzero_samples = uniqueN(sample_id[k > 0])
        ),
        by = feature
    ]

    kept_features = feature_qc[
        total_k >= 20 & nonzero_samples >= 2,
        feature
    ]

    dropped_feature_count = nrow(feature_qc) - length(kept_features)
    dt = dt[feature %in% kept_features]

    filter_log = rbind(
        filter_log,
        data.table::data.table(
            histology = histology_value,
            modality = modality,
            stage = "feature_filter",
            detail = "Applied total_k >= 20 and nonzero_samples >= 2",
            n_removed = dropped_feature_count
        )
    )

    data_checks = data.table::data.table(
        histology = histology_value,
        modality = modality,
        check_name = c(
            "all_n_integer",
            "all_k_integer",
            "k_plus_n_minus_k_equals_n"
        ),
        check_pass = c(
            all(dt$n == as.integer(dt$n)),
            all(dt$k == as.integer(dt$k)),
            all((dt$k + (dt$n - dt$k)) == dt$n)
        )
    )

    if (nrow(dt) == 0 || length(kept_features) == 0) {
        empty_results = data.table::data.table(
            feature = character(),
            histology = character(),
            modality = character(),
            annotation_region = character(),
            term = character(),
            estimate_log_odds = numeric(),
            se = numeric(),
            z_or_t = numeric(),
            p_value = numeric(),
            q_value = numeric(),
            odds_ratio = numeric(),
            ci_low = numeric(),
            ci_high = numeric(),
            model_status = character(),
            status_reason = character(),
            warning_log = character(),
            fit_error = character(),
            random_effect_sd = numeric(),
            sign_call = character()
        )

        model_qc_details = data.table::data.table(
            feature = NA_character_,
            histology = histology_value,
            modality = modality,
            model_status = "fail",
            status_reason = "no_features_after_filters",
            warning_log = "",
            fit_error = "no_features_after_filters",
            random_effect_sd = NA_real_
        )

        model_qc_summary = model_qc_details[, .N, by = .(histology, modality, model_status, status_reason)]

        return(list(
            all_results = empty_results,
            primary_results = empty_results,
            model_qc_summary = model_qc_summary,
            model_qc_details = model_qc_details,
            filter_log = filter_log,
            feature_qc = feature_qc,
            data_checks = data_checks
        ))
    }

    if (data.table::uniqueN(dt$sample_id) < 2) {
        descriptive_base = dt[
            , .(k = sum(k), n = sum(n)),
            by = .(feature, pathology_annotation)
        ]

        if (!(ref_region %in% descriptive_base$pathology_annotation)) {
            ref_region = sort(unique(descriptive_base$pathology_annotation))[1]
        }

        ref_dt = descriptive_base[pathology_annotation == ref_region,
            .(feature, ref_k = k, ref_n = n)
        ]

        descriptive_results = merge(
            descriptive_base,
            ref_dt,
            by = "feature",
            all.x = TRUE
        )[
            pathology_annotation != ref_region
        ][
            , `:=`(
                estimate_log_odds = log((k + 0.5) / (n - k + 0.5)) -
                    log((ref_k + 0.5) / (ref_n - ref_k + 0.5)),
                se = NA_real_,
                z_or_t = NA_real_,
                p_value = NA_real_,
                q_value = NA_real_,
                odds_ratio = exp(estimate_log_odds),
                ci_low = NA_real_,
                ci_high = NA_real_,
                model_status = "warn",
                status_reason = "single_sample_descriptive_only",
                warning_log = "single_sample_descriptive_only",
                fit_error = NA_character_,
                random_effect_sd = NA_real_,
                term = paste0("pathology_annotation", pathology_annotation),
                histology = histology_value,
                modality = modality,
                sign_call = "not_significant"
            )
        ][
            , .(
                feature,
                histology,
                modality,
                annotation_region = pathology_annotation,
                term,
                estimate_log_odds,
                se,
                z_or_t,
                p_value,
                q_value,
                odds_ratio,
                ci_low,
                ci_high,
                model_status,
                status_reason,
                warning_log,
                fit_error,
                random_effect_sd,
                sign_call
            )
        ]

        model_qc_details = descriptive_results[
            , .(
                feature,
                histology,
                modality,
                model_status,
                status_reason,
                warning_log,
                fit_error,
                random_effect_sd
            )
        ]

        model_qc_summary = model_qc_details[, .N, by = .(histology, modality, model_status, status_reason)]

        return(list(
            all_results = descriptive_results,
            primary_results = descriptive_results[0],
            model_qc_summary = model_qc_summary,
            model_qc_details = model_qc_details,
            filter_log = filter_log,
            feature_qc = feature_qc,
            data_checks = data_checks
        ))
    }

    fit_results = lapply(kept_features, function(feature_name) {
        fit_obj = fit_feature_betabinom(
            count_dt = dt,
            feature_name = feature_name,
            ref_region = ref_region,
            sample_col = "sample_id",
            region_col = "pathology_annotation"
        )

        coef_dt = fit_obj$coef_dt
        if (nrow(coef_dt) == 0) {
            coef_dt = data.table::data.table(
                term = NA_character_,
                annotation_region = NA_character_,
                estimate_log_odds = NA_real_,
                se = NA_real_,
                z_or_t = NA_real_,
                p_value = NA_real_,
                odds_ratio = NA_real_,
                ci_low = NA_real_,
                ci_high = NA_real_
            )
        }

        coef_dt[
            , `:=`(
                feature = feature_name,
                histology = histology_value,
                modality = modality,
                model_status = fit_obj$model_status,
                status_reason = fit_obj$status_reason,
                warning_log = fit_obj$warning_log,
                fit_error = fit_obj$fit_error,
                random_effect_sd = fit_obj$random_effect_sd
            )
        ]

        qc_row = data.table::data.table(
            feature = feature_name,
            histology = histology_value,
            modality = modality,
            model_status = fit_obj$model_status,
            status_reason = fit_obj$status_reason,
            warning_log = fit_obj$warning_log,
            fit_error = fit_obj$fit_error,
            random_effect_sd = fit_obj$random_effect_sd
        )

        return(list(coef_dt = coef_dt, qc_row = qc_row))
    })

    all_results = data.table::rbindlist(
        lapply(fit_results, function(x) x$coef_dt),
        use.names = TRUE,
        fill = TRUE
    )

    model_qc_details = data.table::rbindlist(
        lapply(fit_results, function(x) x$qc_row),
        use.names = TRUE,
        fill = TRUE
    )

    all_results[
        , q_value := NA_real_
    ]

    test_rows = all_results[
        annotation_region != "(Intercept)" & is.finite(p_value),
        which = TRUE
    ]

    if (length(test_rows) > 0) {
        if (fdr_scope == "histology_modality") {
            all_results[
                test_rows,
                q_value := stats::p.adjust(p_value, method = "fdr"),
                by = .(histology, modality)
            ]
        }
    }

    all_results[
        , sign_call := data.table::fifelse(
            model_status == "pass" & !is.na(q_value) & q_value < 0.25 & estimate_log_odds > 0,
            "enriched",
            data.table::fifelse(
                model_status == "pass" & !is.na(q_value) & q_value < 0.25 & estimate_log_odds < 0,
                "depleted",
                "not_significant"
            )
        )
    ]

    all_results = all_results[
        , .(
            feature,
            histology,
            modality,
            annotation_region,
            term,
            estimate_log_odds,
            se,
            z_or_t,
            p_value,
            q_value,
            odds_ratio,
            ci_low,
            ci_high,
            model_status,
            status_reason,
            warning_log,
            fit_error,
            random_effect_sd,
            sign_call
        )
    ]

    primary_results = all_results[
        model_status == "pass" & annotation_region != "(Intercept)"
    ]

    model_qc_summary = model_qc_details[
        , .N,
        by = .(histology, modality, model_status, status_reason)
    ][order(histology, modality, model_status)]

    return(list(
        all_results = all_results,
        primary_results = primary_results,
        model_qc_summary = model_qc_summary,
        model_qc_details = model_qc_details,
        filter_log = filter_log,
        feature_qc = feature_qc,
        data_checks = data_checks
    ))
}

summarize_results_for_plot = function(results_dt, primary_only = TRUE) {
    dt = data.table::as.data.table(data.table::copy(results_dt))
    dt = dt[annotation_region != "(Intercept)" & !is.na(annotation_region)]

    if (primary_only) {
        dt = dt[model_status == "pass"]
    }

    if (!"sign_call" %in% names(dt)) {
        dt[
            , sign_call := data.table::fifelse(
                !is.na(q_value) & q_value < 0.05 & estimate_log_odds > 0,
                "enriched",
                data.table::fifelse(
                    !is.na(q_value) & q_value < 0.05 & estimate_log_odds < 0,
                    "depleted",
                    "not_significant"
                )
            )
        ]
    }

    return(dt[])
}


# One-vs-rest histology enrichment via beta-binomial, WITHOUT a random intercept.
# Histology is a between-sample grouping (each sample is exactly one histology,
# one aggregate row per sample), so a sample random intercept is unidentifiable;
# the beta-binomial overdispersion carries the sample-level variance instead.
# For each (feature, histology) it fits cbind(k, n-k) ~ is_histology on the
# per-sample counts (is_histology = this histology vs the rest). Returns the same
# result schema the pathology path produces, so one forest-plot function serves
# both. `count_dt` is the output of build_feature_count_dt() with
# observation_col = "sample_id", region_col = "histology" (one row per sample x
# feature, k = feature cells in the sample, n = total cells in the sample).
fit_feature_histology_set = function(count_dt, sig_threshold = 0.25) {
    dt = data.table::as.data.table(count_dt)
    features = sort(unique(dt$feature))
    histologies = sort(unique(dt$histology))
    rows = list()
    for (f in features) {
        sub_f = dt[feature == f]
        for (h in histologies) {
            sub = data.table::copy(sub_f)
            # binary indicator: this histology vs. the rest (reference = FALSE)
            sub[, is_histology := factor(histology == h, levels = c(FALSE, TRUE))]
            n_samples = sub[histology == h, data.table::uniqueN(sample_id)]
            fit = fit_feature_betabinom(
                count_dt = sub,
                feature_name = f,
                ref_region = "FALSE",
                fixed_effect_col = "is_histology",
                random_effect_col = NULL          # <- no RI (between-sample factor)
            )
            row = data.table::data.table(
                feature = as.character(f), histology = as.character(h),
                n_samples = as.integer(n_samples),
                estimate_log_odds = NA_real_, se = NA_real_, p_value = NA_real_,
                model_status = fit$model_status, status_reason = fit$status_reason,
                warning_log = fit$warning_log, fit_error = fit$fit_error
            )
            if (nrow(fit$coef_dt) > 0) {
                hit = fit$coef_dt[annotation_region == "TRUE"]
                if (nrow(hit) == 1) {
                    row[, `:=`(estimate_log_odds = hit$estimate_log_odds,
                               se = hit$se, p_value = hit$p_value)]
                }
            }
            rows[[length(rows) + 1L]] = row
        }
    }
    res = data.table::rbindlist(rows, fill = TRUE)
    # BH across all (feature, histology) tests in this modality, then call at q<0.25
    res[, q_value := stats::p.adjust(p_value, method = "BH")]
    res[, `:=`(ci_low = estimate_log_odds - 1.96 * se,
               ci_high = estimate_log_odds + 1.96 * se)]
    res[, sign_call := data.table::fifelse(
        model_status == "pass" & !is.na(q_value) & q_value < sig_threshold & estimate_log_odds > 0,
        "enriched",
        data.table::fifelse(
            model_status == "pass" & !is.na(q_value) & q_value < sig_threshold & estimate_log_odds < 0,
            "depleted", "not_significant"
        )
    )]
    return(res[])
}
