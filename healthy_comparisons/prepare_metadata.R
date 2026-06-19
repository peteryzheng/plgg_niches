# prepare_metadata.R
#
# Clean and harmonize the integrated healthy+tumor obs metadata for downstream
# cell-type composition modeling (healthy vs tumor, controlling for age and
# brain location).
#
# Input:  combined_healthyTumor_obs_metadata.parquet  (681,384 cells, 26 cols)
# Output: combined_healthyTumor_obs_metadata_clean.parquet
#
# Key decisions (document here for collaborator review):
#   1. Location crosswalk: healthy atlas `region` and tumor `Location_standard`
#      use different ontologies; mapped to a shared `location_broad` (see below).
#   2. Age: `age_combined` parsed to numeric years; fetal weeks retained as
#      fractional years (e.g. 23w -> 0.44 yr); Unclassified -> NA.
#   3. Sex: `Sex` (tumor) and `donor_gender` (healthy) harmonized to `sex_combined`.
#   4. Neoplastic cells: excluded from composition denominators; per-sample
#      fraction_neoplastic_detected added as a sensitivity covariate.
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

sc_dir    = file.path(data_root, "youyun/plgg/data/single_cell")
in_path   = file.path(sc_dir, "combined_healthyTumor_obs_metadata.parquet")
out_path  = file.path(sc_dir, "combined_healthyTumor_obs_metadata_clean.parquet")

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
# 5. Filter: remove Neoplastic cells and cells with unresolvable location
# ---------------------------------------------------------------------------
n_before = nrow(df)

df = df |>
    filter(cell_type_harmonised != "Neoplastic") |>   # exclude tumor cells from TME denominator
    filter(!is.na(location_broad))                     # drops Spinal, Pending, unmapped regions

cat(sprintf(
    "Filtered: %d -> %d cells (removed %d neoplastic + unresolvable location)\n\n",
    n_before, nrow(df), n_before - nrow(df)
))

# ---------------------------------------------------------------------------
# 6. Summary checks
# ---------------------------------------------------------------------------
cat("Cell type x condition:\n")
print(as.data.frame(table(df$cell_type_harmonised, df$condition)))

cat("\nlocation_broad x condition:\n")
print(table(df$location_broad, df$condition))

cat("\nage_yr NA rate by condition:\n")
print(df |> group_by(condition) |> summarise(pct_na_age = mean(is.na(age_yr))))

cat("\nsex_combined distribution:\n")
print(table(df$sex_combined, df$condition, useNA = "always"))

# ---------------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------------
write_parquet(df, out_path)
cat("\nSaved", nrow(df), "rows ->", out_path, "\n")
