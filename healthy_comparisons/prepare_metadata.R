# prepare_metadata.R
#
# Clean and harmonize the integrated healthy+tumor obs metadata for downstream
# cell-type composition modeling (healthy vs tumor, controlling for age and
# brain location).
#
# Input:  combined_healthyTumor_obs_metadata.parquet  (681,384 cells, 27 cols)
# Outputs (two parquets, same harmonization, differ only in neoplastic handling):
#   - combined_healthyTumor_obs_metadata_clean.parquet
#       Neoplastic cells removed. Denominator = non-neoplastic (TME) cells.
#       Used by the primary composition analysis.
#   - combined_healthyTumor_obs_metadata_withneoplastic.parquet
#       Neoplastic retained as a cell-type category. Denominator = all cells.
#       Used by the sensitivity analysis that keeps neoplastic in the denominator
#       so TME fractions stay comparable across entities with different neoplastic
#       call rates (see celltype_composition_model.qmd).
#
# Key decisions (document here for collaborator review):
#   1. Location crosswalk: healthy atlas `region` and tumor `Location_standard`
#      use different ontologies; mapped to a shared `location_broad` (see below).
#   2. Age: `age_combined` parsed to numeric years; fetal weeks retained as
#      fractional years (e.g. 23w -> 0.44 yr); Unclassified -> NA.
#   3. Sex: `Sex` (tumor) and `donor_gender` (healthy) harmonized to `sex_combined`.
#   4. Neoplastic cells: removed from the primary (clean) parquet's denominator,
#      but retained in the withneoplastic parquet as their own category. Both
#      parquets carry per-sample fraction_neoplastic_detected; it is a covariate
#      for the clean analysis only (it becomes collinear with the denominator
#      once neoplastic cells are counted, so the sensitivity analysis omits it).
#   5. Spinal tumors + Pending location: no healthy BCA reference available;
#      excluded from the location-controlled analysis.

suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(stringr)
})

# ---------------------------------------------------------------------------
# Paths — resolve data root from HOME (repo convention, AGENTS.md)
# ---------------------------------------------------------------------------
home = Sys.getenv("HOME")
if (startsWith(home, "/Users/youyun") || startsWith(home, "/Users/youyunzheng")) {
    data_root = file.path(home, "Documents/HMS/PhD/beroukhimlab/dfci_mount")
} else if (home == "/PHShome/yz762") {
    data_root = "/data/beroukhim1"
} else if (home == "/home/yz762") {
    data_root = "/mnt/storage/dept/medonc/beroukhim"
} else {
    data_root = "/xchip/beroukhimlab"
}

sc_dir        = file.path(data_root, "youyun/plgg/data/single_cell")
in_path       = file.path(sc_dir, "combined_healthyTumor_obs_metadata.parquet")
out_path      = file.path(sc_dir, "combined_healthyTumor_obs_metadata_clean.parquet")
out_path_neo  = file.path(sc_dir, "combined_healthyTumor_obs_metadata_withneoplastic.parquet")

cat("Reading:", in_path, "\n")
df = read_parquet(in_path)
cat("Input:", nrow(df), "cells x", ncol(df), "columns\n\n")

# ---------------------------------------------------------------------------
# 1. Age: parse age_combined to numeric years
# ---------------------------------------------------------------------------
# age_combined merges `age` (tumor, e.g. "16yr") with `donor_age` (healthy,
# e.g. "54.43yr", "23w"). Fetal samples use weeks (e.g. "23w"); convert to
# fractional years so the age covariate is continuous. "Unclassified" -> NA.
df = df |>
    mutate(
        age_yr = case_when(
            str_detect(age_combined, "^[0-9.]+yr$") ~
                as.numeric(str_extract(age_combined, "[0-9.]+")),
            # Prenatal samples: weeks of gestation -> fractional year
            str_detect(age_combined, "^[0-9]+w$") ~
                as.numeric(str_extract(age_combined, "[0-9]+")) / 52,
            TRUE ~ NA_real_   # Unclassified, NA, unexpected formats
        )
    )

cat("age_yr summary:\n")
print(summary(df$age_yr))
cat("age_yr NA count:", sum(is.na(df$age_yr)), "\n\n")

# ---------------------------------------------------------------------------
# 2. Location crosswalk: harmonize to location_broad
# ---------------------------------------------------------------------------
# Healthy cells have `region` (BCA atlas anatomy); tumor cells have
# `Location_standard` (clinical gross location). We map both to a shared
# coarser variable so the model can jointly control for brain region.
#
# Crosswalk rationale:
#   Supratentorial: all cortical + subcortical structures above the tentorium
#     cerebelli (cerebral cortex, basal ganglia, amygdala, hippocampus,
#     thalamus/diencephalon, ganglionic eminences).
#   Posterior_fossa: midbrain/brainstem and posterior fossa structures.
#     (BCA has only Midbrain; cerebellum is absent from the healthy atlas.)
#   Intraventricular: choroid plexus (healthy) maps to intraventricular tumors
#     since choroid plexus lines the ventricles.
#   Spinal / Pending: excluded — no healthy BCA reference for spinal cord;
#     Pending location is unknown.

healthy_region_map = c(
    "Cerebral cortex"      = "Supratentorial",
    "Cortex"               = "Supratentorial",   # alternate label, same anatomy
    "Basal ganglia"        = "Supratentorial",
    "Amygdala"             = "Supratentorial",
    "Hippocampus"          = "Supratentorial",
    "Thalamus"             = "Supratentorial",    # diencephalon, above tentorium
    "Ganglionic eminences" = "Supratentorial",    # fetal basal ganglia precursor
    "Midbrain"             = "Posterior_fossa",   # brainstem = posterior fossa
    "Choroid"              = "Intraventricular"   # choroid plexus lines ventricles
)

tumor_location_map = c(
    "Supratentorial"  = "Supratentorial",
    "Posterior fossa" = "Posterior_fossa",
    "Intraventricular" = "Intraventricular"
    # Spinal + Pending intentionally omitted -> will map to NA -> excluded below
)

df = df |>
    mutate(
        location_broad = case_when(
            condition == "healthy" ~ healthy_region_map[as.character(region)],
            condition == "tumor"   ~ tumor_location_map[as.character(Location_standard)],
            TRUE ~ NA_character_
        )
    )

cat("location_broad distribution:\n")
print(table(df$location_broad, useNA = "always"))
cat("\n")

# ---------------------------------------------------------------------------
# 3. Sex: harmonize Sex (tumor) and donor_gender (healthy) -> sex_combined
# ---------------------------------------------------------------------------
# Tumor cells use `Sex`; healthy BCA cells use `donor_gender`. Standardize
# both to Male/Female; NaN / unrecognised values -> NA.
df = df |>
    mutate(
        sex_combined = case_when(
            condition == "tumor"   ~ as.character(Sex),
            condition == "healthy" ~ as.character(donor_gender),
            TRUE ~ NA_character_
        ),
        sex_combined = case_when(
            sex_combined %in% c("Male",   "M") ~ "Male",
            sex_combined %in% c("Female", "F") ~ "Female",
            TRUE ~ NA_character_
        )
    )

cat("sex_combined NA rate by condition:\n")
print(df |> group_by(condition) |> summarise(pct_na_sex = mean(is.na(sex_combined))))
cat("\n")

# ---------------------------------------------------------------------------
# 3.5 Donor proxy: collapse healthy cells to biological-donor units
# ---------------------------------------------------------------------------
# The pseudobulk unit downstream is donor_proxy × location_broad. For healthy
# cells, sample_ID encodes per-library or per-cell GEO accessions, not
# biological donors, so we coarsen here before saving.
#
# tumor:           donor_proxy = sample_ID  (one patient per ID, correct as-is)
# healthy non-GSM: paste(study_prefix, donor_age, donor_gender)
#                  Groups multi-batch / multi-region IDs from the same donor.
#                  E.g. two Tran 2021 runs of the same 54yr Male donor collapse.
# healthy GSM:     paste(region, donor_age, donor_gender)
#                  Smart-seq2 datasets use per-cell GSM IDs with no study
#                  prefix; region is added to prevent merging donors from
#                  different anatomical origins that share age/sex.
df = df |>
    mutate(
        .study_prefix = case_when(
            condition == "tumor"         ~ NA_character_,
            str_detect(sample_ID, "GSM") ~ NA_character_,
            TRUE ~ str_extract(sample_ID, "^BCA_[^_]+")  # e.g. "BCA_Tran"
        ),
        donor_proxy = case_when(
            condition == "tumor" ~ as.character(sample_ID),
            str_detect(sample_ID, "GSM") ~
                paste(as.character(region), as.character(donor_age),
                      as.character(donor_gender), sep = "|"),
            TRUE ~
                paste(.study_prefix, as.character(donor_age),
                      as.character(donor_gender), sep = "|")
        )
    ) |>
    select(-.study_prefix)

cat("donor_proxy unique counts by condition:\n")
print(df |> group_by(condition) |>
    summarise(n_donor_proxies = n_distinct(donor_proxy), .groups = "drop"))
cat("\n")

# ---------------------------------------------------------------------------
# 4. Per-sample fraction_neoplastic_detected (before removing neoplastic)
# ---------------------------------------------------------------------------
# Compute the fraction of sequenced cells classified as Neoplastic per sample.
# This captures variable neoplastic cell detectability across samples and will
# serve as a sensitivity covariate in downstream models.
neoplastic_frac = df |>
    group_by(sample_ID) |>
    summarise(
        n_cells_total       = n(),
        n_neoplastic        = sum(cell_type_harmonised == "Neoplastic"),
        fraction_neoplastic_detected = n_neoplastic / n_cells_total,
        .groups = "drop"
    )

df = df |> left_join(neoplastic_frac |> select(sample_ID, fraction_neoplastic_detected),
                     by = "sample_ID")

cat("fraction_neoplastic_detected (tumor samples only):\n")
print(summary(df$fraction_neoplastic_detected[df$condition == "tumor"]))
cat("\n")

# ---------------------------------------------------------------------------
# 5. Filter location, then split into the two parquets
# ---------------------------------------------------------------------------
# Drop cells with unresolvable location first (Spinal, Pending, unmapped
# regions) — this applies to BOTH parquets. Then:
#   df_neo   = location-resolved cells WITH neoplastic (all-cell denominator)
#   df       = df_neo minus neoplastic (non-neoplastic / TME denominator)
n_before = nrow(df)

df_neo = df |>
    filter(!is.na(location_broad))                     # shared: drops Spinal, Pending, unmapped

df = df_neo |>
    filter(cell_type_harmonised != "Neoplastic")       # clean parquet: TME denominator only

cat(sprintf(
    "Location-resolved: %d -> %d cells (dropped %d unresolvable location)\n",
    n_before, nrow(df_neo), n_before - nrow(df_neo)
))
cat(sprintf(
    "Neoplastic removed for clean parquet: %d -> %d cells (removed %d neoplastic)\n\n",
    nrow(df_neo), nrow(df), nrow(df_neo) - nrow(df)
))

# ---------------------------------------------------------------------------
# 6. Summary checks (on the clean / non-neoplastic set)
# ---------------------------------------------------------------------------
cat("Cell type x condition (clean, non-neoplastic):\n")
print(as.data.frame(table(df$cell_type_harmonised, df$condition)))

cat("\nlocation_broad x condition (clean):\n")
print(table(df$location_broad, df$condition))

cat("\nage_yr NA rate by condition (clean):\n")
print(df |> group_by(condition) |> summarise(pct_na_age = mean(is.na(age_yr))))

cat("\nsex_combined distribution (clean):\n")
print(table(df$sex_combined, df$condition, useNA = "always"))

cat("\nNeoplastic cells retained in withneoplastic parquet (tumor only):\n")
print(table(df_neo$condition[df_neo$cell_type_harmonised == "Neoplastic"]))

# ---------------------------------------------------------------------------
# 7. Save both parquets
# ---------------------------------------------------------------------------
write_parquet(df, out_path)
cat("\nSaved", nrow(df), "rows (clean, no neoplastic) ->", out_path, "\n")

write_parquet(df_neo, out_path_neo)
cat("Saved", nrow(df_neo), "rows (with neoplastic) ->", out_path_neo, "\n")
