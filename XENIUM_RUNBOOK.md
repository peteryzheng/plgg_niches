# Xenium Runbook — Essential Run Order

This document lists the minimum set of code needed to take the PLGG Xenium
cohort from raw transcripts to the niche / cell-type / pathology /
SC-integration analyses we report on. Each stage describes **how to run
it**, **inputs**, **outputs**, and the **parameters/arguments** that
matter.

> Folder layout in this repo is described in `AGENTS.md`. Read that first
> if any path below is ambiguous.

## Contents
- [0. Environment + path conventions](#0-environment--path-conventions)
- [1. ProSeg segmentation](#1-proseg-segmentation)
- [2. BANKSY cell-type + niche clustering (cohort)](#2-banksy-cell-type--niche-clustering-cohort)
- [3. Pathology annotation integration (STalign)](#3-pathology-annotation-integration-stalign)
- [4. Single-cell integration (ENVI)](#4-single-cell-integration-envi)
- [5. Cohort-level descriptors](#5-cohort-level-descriptors)
- [6. Optional QC / diagnostic notebooks](#6-optional-qc--diagnostic-notebooks)
- [End-to-end run order (critical path)](#end-to-end-run-order-critical-path)
- [Conventions enforced everywhere](#conventions-enforced-everywhere)

---

## 0. Environment + path conventions

### 0a. Conda environment
The single env is `spatial` (R + Python + Quarto). Activate before any
render:
```bash
conda activate spatial
```
Do **not** use `conda run -n spatial quarto render ...` for QMDs — it
triggers a Deno/Quarto version conflict. The QC qsub script is the one
exception (it runs in a clean shell).

### 0b. Cluster filesystem
Active root is `/mnt/storage/dept/medonc/beroukhim/` (new server,
`$HOME=/home/yz762`). `/xchip/beroukhimlab/` is legacy and read-only for
back-reference. Qsub scripts source
`/mnt/storage/dept/medonc/beroukhim/youyun/util/miniforge3/etc/profile.d/conda.sh`
before `conda activate spatial`.

### 0c. Path resolution in code
Every R/QMD/PY entry point derives `workdir` from `$HOME` (see
`AGENTS.md` mapping). Never hard-code one machine's path. Data root
(under `workdir`) is `youyun/plgg/data/` for all inputs/outputs.

### 0d. Shared helpers (`helpers/`)
- `spatial_helper.R` — loading, QC, BANKSY workflow.
- `test_enrichment.R` — legacy categorical/continuous enrichment (Wilcoxon/Fisher);
  still used for the within-sample niche↔cell-type composition test.
- `enrichment_models.R` — **unified beta-binomial GLMM engine** for all
  label-vs-grouping enrichment (cell type / niche × histology **or** pathology
  region). `build_feature_count_dt` → `fit_feature_set` (pathology) /
  `fit_feature_histology_set` (histology). Random intercept `(1|sample_id)` is
  included **only when the grouping varies within a sample** (pathology regions);
  it is omitted for histology (one aggregate row per sample, so a sample RI is
  unidentifiable — the beta-binomial overdispersion carries the sample-level
  variance). Formerly `pathology_mixed_models.R`. Note: histology q-values are
  anti-conservative at small n (2–4 samples/group) — report effect size (logOR +
  CI), treat q as ranking only.
- `cell_type_display.R` — snake_case lineage key → human-readable label map,
  applied at plot time (keeps canonical `marker_panel` keys stable).
- `spatial_helper_visium.R` — Visium-only (not on Xenium critical path).

---

## 1. ProSeg segmentation

Re-segments each Xenium sample from raw transcripts and imports the
result back into a Xeniumranger bundle.

### 1a. Scripts
- `segmentation/proseg.sh` — single-sample driver.
- `segmentation/run_proseg.sh` — UGE/SGE array launcher (one task per
  sample).
- `segmentation/transcripts.txt` — sample list, one transcript
  `transcripts.csv.gz` path per line; `SGE_TASK_ID` indexes into it.

### 1b. How to run
```bash
qsub segmentation/run_proseg.sh
# or interactively for one sample:
bash segmentation/proseg.sh <transcripts.csv.gz> <output_dir> <threads>
```
Resources (in `run_proseg.sh`): `h_rt=48h`, 8 threads, `h_vmem=64G`.

### 1c. Inputs / outputs
- **Inputs**: `transcripts.csv.gz` from a Xenium Analyzer reupload bundle.
- **Outputs** (per sample, under
  `youyun/plgg/data/segmentation/proseg_run_121024/<project_id>/`):
  - `transcript-metadata.csv.gz`, `cell-polygons.geojson.gz` (proseg native).
  - `proseg-transcript-metadata-baysor-import.csv` and
    `proseg-cell-polygons-baysor-import.geojson` (Baysor-compatible).
  - A new Xeniumranger bundle `<project_id>_proseg/` produced by
    `xeniumranger import-segmentation --units microns`.

### 1d. Notes
- `xeniumranger` lives at
  `/xchip/beroukhimlab/youyun/util/xeniumranger/xeniumranger-xenium3.0/`.
- The import step works on compute nodes but not on `dipg`.

### 1e. Optional: QC threshold sweep
- `QC/determine_qc_thresholds.qmd` — MAD-based sweep of total count /
  detected genes thresholds to pick QC cutoffs.
- `QC/determine_qc_threshold.sh` — qsub wrapper (96G, 48h); uses
  `conda run -n spatial quarto render` because no interactive shell
  exists.
- Render locally:
  ```bash
  conda activate spatial
  quarto render QC/determine_qc_thresholds.qmd
  ```

---

## 2. BANKSY cell-type + niche clustering (cohort)

Runs the full BANKSY workflow on the ProSeg-segmented cohort: load → QC
→ normalize → BANKSY embeddings (2 lambdas) → Harmony → Leiden
clustering for both cell type (`lam=0.2`) and niche (`lam=0.8`) →
cell-type markers. Then interpret/annotate the resulting clusters.

### 2a. Run the cohort
- **Driver R**: `niche/banksy_workflow/banksy_cohort.R`
- **Array launcher**: `niche/banksy_workflow/banksy_cohort_qsub.sh`
- **Param grid**: `niche/banksy_workflow/param_search.tsv` — one config
  per row. Columns: `k1 k2 lambda1 lambda2 npcs k_ct res_ct k_ni res_ni`.
  In practice this file holds a single row, always **the current active
  config**; history of prior configs lives in git log/commit messages, not
  accumulated rows (see the header comment in `banksy_cohort_qsub.sh`).
- **Active config (`param_search.tsv`, single row)**:
  `k_geom = (15,30)`, `lambda = (0.2,0.8)`, `npc=20`,
  `k_ct=50, res_ct=0.2`, `k_ni=50, res_ni=0.5`, `seed=55555`. The two
  resolutions here are the **final chosen** values (single, not
  comma-separated). The file is deliberately overwritten to the current
  finals rather than accumulating sweep rows — prior configs live in git,
  not in extra rows (see the `banksy_cohort_qsub.sh` header).
- **Two routes to the production object** — these are *alternatives*, not a
  sequence, and both go through the same automatic `annotate_cell_types()` call
  (so neither needs a separate cell-typing notebook):
  1. **From-scratch (canonical one-step reproduction)** — `qsub
     banksy_cohort_qsub.sh` with the current single-value row. `banksy_cohort.R`
     runs the entire workflow *including* `connectClusters()` (inside
     `banksy_clustering()`) and `annotate_cell_types()`, so this one submission
     yields the final annotated object directly at `res_ct=0.2 / res_ni=0.5`.
     BANKSY, Harmony, and Leiden are all deterministic given the fixed seed
     (55555), so this reproduces a scientifically identical object (~35–55h).
  2. **Fast-path (how the on-disk object was actually built)** — the base run
     was launched at the **sweep** resolutions `res_ct=0.5,1 / res_ni=0.5,1`
     (two values each — hence the base RDS/dir carry `resct_0.5,1 … resni_0.5,1`),
     then §2c re-tuned the cell-type axis down to `res_ct=0.2` on the full
     object and finalized. Only the resolution-independent embeddings
     (BANKSY/Harmony/UMAP) and the niche `res=0.5` clustering carry directly
     from the base run; the cell-type `res=0.2` clustering was produced by the
     §2c recluster, **not** the base run. This is the fast-path the old
     "provenance note" referred to — it avoided a ~35–55h from-scratch rerun.
- **Reproducing the fast-path**: because `param_search.tsv` now holds only the
  final single-value row, this qsub script reproduces route 1 (from-scratch) but
  no longer route 2 — the `0.5,1` sweep row route 2 needed is gone (recover it
  from git if you must re-run the exact fast-path). See §2c and
  `finalize_final_resolution.R`'s header for the full provenance chain.
- **Resources**: 1 slot, `h_vmem=96G`, `h_rt=120h`.
- **Run**:
  ```bash
  qsub niche/banksy_workflow/banksy_cohort_qsub.sh
  ```
- **Helpers used** (all in `helpers/spatial_helper.R`):
  `loading_data(segmentation_method='proseg')`, `stagger_spatial_coords()`,
  `QC_and_normalize()`, `banksy_workflow()`, `banksy_clustering()`,
  `cell_type_marker_ident()`.
- **Inputs**: ProSeg outputs from stage 1 +
  `youyun/plgg/data/metadata/Xenium_PS.xlsx` (sample-to-file mapping
  used by `loading_data`).
- **Outputs** (under `youyun/plgg/data/banksy_param_search/<config_dir>/`):
  - Full BANKSY `SpatialExperiment` RDS:
    `banksy_clusters_connected_k_geom_..._<timestamp>.rds`
  - 5% downsample for fast exploration:
    `banksy_clusters_connected_subset_k_geom_..._<timestamp>.rds`
  - Cell-type marker tables.

### 2b. Optional: Leiden resolution sweep (5% subset)
Reuses the saved BANKSY subset RDS (no re-running BANKSY/Harmony):
- `niche/banksy_workflow/leiden_resolution_sweep.R` — runs Leiden over a
  `k × res` grid for both lambdas; writes
  `leiden_resolution_sweep_counts.tsv` next to itself.
- `niche/banksy_workflow/evaluate_leiden_sweep.qmd` — 3-layer
  evaluation: cluster-count plateaus → biology check on finalists →
  decision table.
- **Run**:
  ```bash
  conda activate spatial
  Rscript niche/banksy_workflow/leiden_resolution_sweep.R
  quarto render niche/banksy_workflow/evaluate_leiden_sweep.qmd
  ```
- **Caveat**: a resolution that looks satisfactory on the subset can
  over-segment on the full object (cluster granularity at fixed resolution
  is not scale-invariant) — any candidate needs to be re-tested on the full
  object (§2c below) before being trusted.

### 2c. Cell-type resolution re-tuning + finalizing the production object (full object)
This is the **fast-path (§2a route 2) that actually built the current
production object** — and the tool set for re-tuning the cell-type (lam0.2)
resolution on the **full** annotated object (not the 5% subset) without
re-running BANKSY/Harmony/UMAP. Skip it only if you rebuilt from scratch via
§2a route 1 (which already finalizes + annotates in one pass).
- `niche/banksy_workflow/celltype_lowres_recluster.R` — cheap Leiden-only
  re-cluster at new resolution candidate(s), reusing the existing
  `Harmony_BANKSY_lam0.2` embedding. **Cell-type axis only** —
  `annotate_cell_types()`/`generate_celltype_review_panels()` hardcode
  `min(lambda_vec)`, so pointing this at lam0.8 silently annotates/plots the
  wrong (untouched) clusters; see the script's own header comment for the
  full mechanism and why a niche equivalent would need its own script.
- `niche/banksy_workflow/cluster_summary_stats.R` — per-cluster summary
  table (size, entropy, all lineage scores, margin, top unbiased genes)
  across one or more annotated objects/resolutions, for comparing
  cell-type-calling heuristics quantitatively rather than eyeballing panels.
- `niche/banksy_workflow/finalize_final_resolution.R` — once a resolution
  is chosen, drops the non-final resolution columns, applies
  `Banksy::connectClusters()` (the one step the cheap recluster above
  skips — see its own header comment), re-annotates with the current
  default cell-type-calling heuristic
  (`score_lineages_top_gene()`/`min_gene_auc=0.55` in
  `helpers/spatial_helper.R`), and regenerates QC/review outputs. This is
  the script that produces the object treated as final production output.
- **Cell typing is automatic here** — `annotate_cell_types()` writes the
  marker-AUC lineage call to `cell_type_clust_Harmony_BANKSY_lam0.2_k50_res0.2`
  on the object. This is *the* cell-type source for all downstream work; the old
  marker-annotation notebook (`banksy_clusters_proseg.qmd`, §2d) and its
  `banksy_proseg_cell_cluster_annotation.tsv` are superseded and no longer read.
  Note `celltype_lowres_recluster.R` only tunes the resolution *candidate*; the
  final production label is the re-annotate in `finalize_final_resolution.R`
  (run after `connectClusters()`).
- **Final production object**:
  `banksy_param_search/k1_15_..._resct_0.5,1_kni_50_resni_0.5,1/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134.rds`
  (the dir name still carries the base run's `0.5,1` sweep suffix; the file
  itself is at the final `resct_0.2/resni_0.5`). Its pathology-merged copy
  (§3f) is `..._20260826_115134_pathology.rds`.

### 2d. Cluster annotation + spatial visualization
All under `niche/proseg_output_analysis/`.

**Cell-type labels are no longer produced by a notebook.** They are assigned
automatically by `annotate_cell_types()` in the R pipeline (§2a route 1 / §2c
finalize) and stored on the object as
`cell_type_clust_Harmony_BANKSY_lam0.2_k50_res0.2`; human-readable lineage
labels come from `helpers/cell_type_display.R` at plot time. **Niche labels /
interpretation** are the "niche dictionary" in `niche_interpretation_proseg.qmd`
(§2e), not a marker-annotation notebook.

**Spatial visualization per sample** (the still-relevant notebooks here):
- `banksy_cell_type_in_space_proseg.qmd` — cell types per sample.
- `banksy_niche_in_space_proseg.qmd` — niches per sample.
- ⚠️ **Not yet rewired to the active object.** Both still `readRDS()` the OLD
  2024-12-30 run (`…lam_0.2_0.8_k_leiden_30_50_res_0.75_1_20241230_232215…`), so
  their maps do **not** reflect the current `resct_0.2 / resni_0.5` object.
  Repoint them at the §2c final object (or the §3f pathology copy) before
  trusting their output. Once rewired, render with:
  ```bash
  conda activate spatial
  quarto render niche/proseg_output_analysis/banksy_cell_type_in_space_proseg.qmd
  ```

**Superseded** (kept for reference only — still load the OLD 2024-12-30 object;
do **not** treat their output as current):
- `banksy_clusters_proseg.qmd` — old marker-based cell-type annotation; produced
  `banksy_proseg_cell_cluster_annotation.tsv`, never regenerated at the final
  resolution and no longer read. Replaced by the automatic `annotate_cell_types()`
  call (§2c).
- `banksy_niches_proseg.qmd` — old niche annotation + cross-sample tests.
  Replaced by `niche_interpretation_proseg.qmd` (§2e).

### 2e. Cell-type / niche enrichment + interpretation (histology axis)
These notebooks need only the BANKSY object (no pathology required) and live in
`niche/proseg_output_analysis/`. `celltype_niche_histology_enrichment.qmd` and
`niche_interpretation_proseg.qmd` are the **current** annotation/interpretation
notebooks on the active object (they read `cell_type_clust_...lam0.2_k50_res0.2`
directly, replacing the superseded §2d marker notebooks).
- `celltype_niche_histology_enrichment.qmd` — **current** cell-type × histology and
  niche × histology enrichment via the unified beta-binomial engine
  (`fit_feature_histology_set`, one-vs-rest, no random intercept), plus per-sample
  composition bars. Writes `celltype_histology_enrichment.tsv`,
  `niche_histology_enrichment.tsv`. Report logOR + CI (q is ranking-only,
  anti-conservative at n=2–4). This replaces the earlier Wilcoxon approach.
- `niche_interpretation_proseg.qmd` — the **niche "dictionary"**: per-niche
  neighborhood cell-type composition (self + kNN neighbors pooled), enrichment vs
  cohort, pathology top-region, enriched histology, and a sample-private
  (candidate-tumor-niche) flag → `niche_dictionary.tsv`.
- `banksy_niches_glm_proseg.qmd`, `banksy_niches_by_histology_proseg.qmd` — older
  GLM / histology-stratified niche analyses (superseded by the two above for the
  final object; kept for reference — note these still load the OLD 2024-12-30
  object, not the active one).
- **Helpers used**: `helpers/enrichment_models.R`, `helpers/cell_type_display.R`,
  `helpers/test_enrichment.R`.
- **Run**:
  ```bash
  conda activate spatial
  quarto render niche/proseg_output_analysis/celltype_niche_histology_enrichment.qmd
  ```

---

## 3. Pathology annotation integration (STalign)

Aligns pathologist GeoJSON polygons (drawn on H&E in QuPath) onto Xenium
space, writes per-cell pathology labels, merges them back into the
BANKSY object, then runs pathology × cell-type / niche enrichment.

### 3a. Per-sample STalign QMDs
`annotations/STalign/pathology_<SAMPLE>_proseg.qmd` — 11 files, one per
cohort sample. These are the source of truth (edits go here).

### 3b. Convert QMD → IPYNB (needed for the cluster step)
The cluster runs `.ipynb` via `jupyter nbconvert`, so convert first:
```bash
conda activate spatial
bash annotations/STalign/convert_qmd_to_ipynb_local.sh           # all
bash annotations/STalign/convert_qmd_to_ipynb_local.sh PA258482  # one
```
Generated `.ipynb` files go to `annotations/STalign/ipynb_generated/`.

### 3c. Array launcher
`annotations/STalign/stalign_qmd_array_qsub.sh` executes the ipynbs and
renders HTML.
```bash
qsub annotations/STalign/stalign_qmd_array_qsub.sh
```
- `-t 1-11`, 64G, 24h.
- Falls back to `legacy_ipynb/` if `ipynb_generated/` is missing.
- Output HTML lands in `annotations/STalign/ipynb_rendered/`.

### 3d. Inputs / outputs
- **Inputs**:
  - GeoJSONs at `youyun/plgg/data/Xenium_annotations/geojsons/`
    (e.g. `230918_Xenium_CytAssist_LGG1_merged.geojson`).
  - JPGs at `youyun/plgg/data/Xenium_annotations/images/`.
  - ProSeg `cells.csv.gz` at
    `youyun/plgg/data/segmentation/proseg_run_121024/<file>/<file>_proseg/outs/`.
  - Sample-id → file mapping: `youyun/plgg/data/metadata/Xenium_PS.xlsx`.
- **Outputs (per sample)**:
  `youyun/plgg/data/Xenium_annotations/pathology_annotations_030226/proseg_pathology_<SAMPLE>_proseg_meta.csv`
  with columns `cell_id, annotation, region_id`.

### 3e. Optional: rebuild the LGG1 merged GeoJSON
`annotations/STalign/build_lgg1_merged_geojson.py` merges multiple
GeoJSONs into the canonical `LGG1_merged.geojson` used by the four LGG1
samples.

### 3f. Merge pathology back into the BANKSY object
- **Script**: `annotations/STalign/add_annotations.R`
- **Run**:
  ```bash
  conda activate spatial
  Rscript annotations/STalign/add_annotations.R
  ```
- Reads every `proseg_pathology_<SAMPLE>_proseg_meta.csv` from 3d, joins
  by `(cell_id, sample_id)`, normalizes pathology labels (e.g. hex
  colors → human names), writes:
  - Annotated full RDS:
    `youyun/plgg/data/Xenium_annotations/banksy_clusters_connected_..._annotated.rds`
  - Annotated 5% subset RDS (re-sampled fresh each run):
    `..._connected_subset_..._annotated.rds`
  - Adds `pathology_annotation` + `region_id` columns to `colData`.

### 3g. Pathology × cell-type / niche enrichment (`annotations/`)
These QMDs consume the **annotated** RDS from 3f. Beta-binomial GLMM
with sample random intercept (`glmmTMB`), modeled per-feature
independently.
- `annotations/linear_models.qmd` — **cell-type** proportions across
  pathology `region_id`s.
- `annotations/linear_models_niches.qmd` — same model on **niche**
  proportions.
- `annotations/linear_models_sample_level.qmd` — sample-level variant
  (sample as the observation).
- **Helpers used**: `helpers/enrichment_models.R`
  (`build_feature_count_dt`, `fit_feature_betabinom`, `fit_feature_set`,
  `summarize_results_for_plot`); `helpers/test_enrichment.R` for the
  pairwise enrichment fallbacks. **CN reference region = `unannotated`** (the
  former `Neuropil zone` ref was tiny/noisy; `unannotated` matches PA/GG and is
  large in both CN samples). Cell-type is taken from the object's own
  `cell_type_<clust>` marker call, not the old res1 annotation TSV.
- **Outputs**: `annotations/cln_associations.csv`,
  `annotations/pa_associations.csv`, and rendered HTML reports.
- **Run**:
  ```bash
  conda activate spatial
  quarto render annotations/linear_models.qmd
  ```

### 3h. Auxiliary annotation notebooks
- `annotations/cancer_cells.qmd` — flag enriched cells in annotated
  regions as cancer cells (earlier-style Fisher → linear-model pipeline).
- `annotations/cell_positions.qmd` — cell-position sanity checks against
  raw Xenium outputs.
- `annotations/CCA_in_annotated_regions.qmd` — CCA-style analysis inside
  pathology regions.

---

## 4. Single-cell integration (ENVI)

Two phases: prepare training AnnDatas, then run ENVI per histology and
use the imputation / co-embedding for downstream analyses.

### 4a. Training data prep — `sc_integration/ENVI/1_training_data/`
- `sce_2_anndata.qmd` — convert the BANKSY `SpatialExperiment` RDS into
  Python `AnnData` (`sceasy`).
- `seurat_2_anndata.qmd` — convert the SN Seurat (`sn_LGG.rds`) to
  AnnData (raw counts; ENVI does its own log1p).
- `h5ad_2_sample_h5ad.qmd` — split the cohort SN AnnData
  (`Extended_sndata_filtered.h5ad`) into per-`<sample_id>_<histology>`
  h5ads using the BANKSY metadata as the sample→histology map.
  - **Output dir**: `youyun/plgg/data/sc_integration/`
  - **Log**: `youyun/plgg/data/sc_integration/h5ad_2_sample_h5ad_log.tsv`

Render each:
```bash
conda activate spatial
quarto render sc_integration/ENVI/1_training_data/h5ad_2_sample_h5ad.qmd
```

### 4b. ENVI training (external)
The ENVI training run itself (writing `*_sn_envi.h5ad` and
`*_spatial_envi.h5ad` to
`youyun/plgg/data/sc_integration/ENVI_results_022626/`) lives outside
this repo. Downstream notebooks consume those h5ads.

### 4c. PA-specific downstream — `sc_integration/ENVI/2_spatial_imputation/PA/`
- `1_PA_MAPK_scores.qmd` — compare MAPK pathway scores between SN and
  spatial-imputed PA samples. Uses
  `sc_integration/ENVI/2_spatial_imputation/MAPK_signatures.csv`.
- `2_MAPK_adj.qmd` — MAPK adjacency analysis.
- `3_PA_myeloids.qmd` — Myeloid 1 vs Myeloid 2 pathway analysis. Triggers
  two parallel GSEA tracks (see `README_myeloid_gsea.md`).

#### 4c-i. Spatial-imputation GSEA track
`PA/spatial_imputation_gsea/` — pathway analysis on the ENVI-imputed
matrix.
- Driver: `myeloid_spatial_imputation_prepare.py` (build pseudobulks) →
  `myeloid_spatial_imputation.R` (limma + GSEA).
- Outputs: `myeloid_spatial_imputation_*` CSVs in the same folder.

#### 4c-ii. SN label-transfer GSEA track
`PA/sn_label_transfer_gsea/` — pathway analysis on SN observed counts
after spatial→SN label transfer in `latent_umap`.
- Same `_prepare.py` + `.R` pattern.
- Co-embedding PNGs in `PA/sn_label_transfer_gsea/coembedding/`.

#### 4c-iii. Gene-universe audit
`PA/PA_imputation_gene_audit_022626.md` documents the canonical
gene-universe extraction logic for ENVI h5ads — read it if any
imputation gene set looks off.

### 4d. Other histologies — single-gene imputation checks
- `sc_integration/ENVI/2_spatial_imputation/GG/GG_CD34.qmd`
- `sc_integration/ENVI/2_spatial_imputation/CN/CN_NKX2-1.qmd`

Render each:
```bash
conda activate spatial
quarto render sc_integration/ENVI/2_spatial_imputation/PA/3_PA_myeloids.qmd
```

---

## 5. Cohort-level descriptors

Run once per cohort update; rerender as needed.

### 5a. Cohort metadata
- `cohort_metadata/cohort_description.qmd` — cohort summary tables and
  figures.

### 5b. Single-nucleus composition
- `sample_comparisons/single_nuc_samples.qmd` — beta-binomial proportion
  modeling of SN cell-type composition against covariates.
- **Input**: `youyun/plgg/data/single_cell/harmony_lambda2_t1.rds`
  (originally pulled from the Terra bucket —
  `gs://fc-secure-351a0c20-.../Extended_sn/`).
- **Outputs** in `sample_comparisons/`:
  - `betabinom_*_results.tsv`
  - `glm_quasibinomial_annot_*_results.tsv`

---

## 6. Optional QC / diagnostic notebooks

Not on the critical path but useful for sanity checks.

### 6a. Unassigned transcripts
`segmentation/1_unassigned_transcripts/unassigned_transcripts.qmd` —
per-sample and per-cell unassigned-transcript fractions after ProSeg.
Primary focus is `dx == "LIPN"`.

### 6b. CRAWDAD multiscale niches
`niche/crawdad/crawdad_niches.qmd` — alternative multiscale niche method
for sanity-checking BANKSY niches on two samples (CLN 87352, PA 258482).

### 6c. Visium-only workflow
`niche/visium_banksy_workflow/` — Visium BANKSY workflow; **not** part
of the Xenium critical path.

---

## End-to-end run order (critical path)

```
1. segmentation/run_proseg.sh                                       (qsub array)
2a. niche/banksy_workflow/banksy_cohort_qsub.sh                     (qsub; single-row param_search.tsv = final res_ct=0.2/res_ni=0.5. Route 1 finalizes + annotates in one pass -> skip 2c)
2b. (opt) niche/banksy_workflow/leiden_resolution_sweep.R + evaluate_leiden_sweep.qmd   (5% subset)
2c. niche/banksy_workflow/{celltype_lowres_recluster,cluster_summary_stats,finalize_final_resolution}.R   (full-object re-tune + finalize -> production object; the fast-path that built the on-disk object)
2d. niche/proseg_output_analysis/banksy_{cell_type,niche}_in_space_proseg.qmd   (spatial viz; NOT yet rewired to active object)
2e. niche/proseg_output_analysis/celltype_niche_histology_enrichment.qmd + niche_interpretation_proseg.qmd   (enrichment + niche dictionary, active object)
3a-c. annotations/STalign/convert_qmd_to_ipynb_local.sh             (local)
      annotations/STalign/stalign_qmd_array_qsub.sh                 (qsub array)
3f.   annotations/STalign/add_annotations.R                         (Rscript)
3g.   annotations/linear_models{,_niches,_sample_level}.qmd         (quarto render)
4a. sc_integration/ENVI/1_training_data/{sce_2_anndata,seurat_2_anndata,h5ad_2_sample_h5ad}.qmd
4b. (external) ENVI training -> ENVI_results_022626/*_envi.h5ad
4c. sc_integration/ENVI/2_spatial_imputation/PA/{1,2,3}*.qmd        (+ GSEA scripts)
5.  cohort_metadata/cohort_description.qmd, sample_comparisons/single_nuc_samples.qmd
```

## Conventions enforced everywhere
- `workdir` derived from `$HOME` (mac, eristwo `/PHShome/yz762`, new
  server `/home/yz762`, legacy UGER). Never hard-code a single absolute
  prefix.
- BANKSY object identity travels through the pipeline by the suffix
  `k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134`.
  Files with this suffix in `banksy_param_search/` (`…_annotated.rds`) and
  `Xenium_annotations/` (`…_pathology.rds`, §3f) are the same **active** cohort
  run; do not mix with other timestamps. The older
  `…res_0.75_1_20241230_232215` suffix is the **superseded** 2024-12-30 run —
  the not-yet-rewired `*_in_space_proseg.qmd` notebooks (§2d) still point at it.
- Generated outputs (`*.html`, `figure-html/`, `libs/`, `._*`) are not
  committed — see `.gitignore`.
