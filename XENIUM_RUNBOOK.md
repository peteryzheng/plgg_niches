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
- **Active config (row 1)**:
  `k_geom = (15,30)`, `lambda = (0.2,0.8)`, `npc=20`,
  `k_ct=50, res_ct=0.2`, `k_ni=50, res_ni=0.5`, `seed=55555`.
- **Provenance note**: the current production annotated object was *not*
  produced by a fresh run of this qsub script at the active config above.
  It was assembled via a validated fast-path
  (`niche/banksy_workflow/finalize_final_resolution.R`) that reuses the
  deterministic BANKSY/Harmony/UMAP embeddings + Leiden clustering already
  computed for this exact config (re-running this script from scratch would
  cost ~35-55h to reproduce a scientifically identical result — BANKSY,
  Harmony, and Leiden are all deterministic given the same seed/data). See
  §2b below for the resolution-retuning + finalization tools and
  `finalize_final_resolution.R`'s header comment for the full provenance
  chain. This qsub script remains the correct from-scratch fallback if that
  equivalence ever needs re-validating or the upstream data changes.
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

### 2b. Optional: Leiden resolution sweep
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
  object (below) before being trusted.

### 2b-2. Cell-type resolution re-tuning + finalizing on the full object
For re-testing lam0.2 (cell-type) resolution candidates on the **full**
annotated object (not the subset) without re-running BANKSY/Harmony/UMAP:
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

### 2c. Cluster annotation + spatial visualization
All under `niche/proseg_output_analysis/`. Render order:
1. `banksy_clusters_proseg.qmd` — annotate cell-type clusters (lam0.2)
   with markers. Produces
   `youyun/plgg/data/xenium_cell_types/banksy_proseg_cell_cluster_annotation.tsv`.
2. `banksy_cell_type_in_space_proseg.qmd` — visualize cell types per
   sample.
3. `banksy_niches_proseg.qmd` — annotate niche clusters (lam0.8) by
   cell-type composition + cross-sample tests.
4. `banksy_niche_in_space_proseg.qmd` — visualize niches per sample.

Render any of these with:
```bash
conda activate spatial
quarto render niche/proseg_output_analysis/banksy_clusters_proseg.qmd
```

### 2d. Cell-type / niche enrichment modeling (histology axis)
These notebooks need only the BANKSY object (no pathology required) and live in
`niche/proseg_output_analysis/`.
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
  final object; kept for reference).
- **Helpers used**: `helpers/enrichment_models.R`, `helpers/cell_type_display.R`,
  `helpers/test_enrichment.R`.
- **Run**:
  ```bash
  conda activate spatial
  quarto render niche/proseg_output_analysis/banksy_niches_glm_proseg.qmd
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
2a. niche/banksy_workflow/banksy_cohort_qsub.sh                     (qsub, row 1 of param_search.tsv)
2b. (opt) niche/banksy_workflow/leiden_resolution_sweep.R + evaluate_leiden_sweep.qmd
2c. niche/proseg_output_analysis/banksy_{clusters,niches,*_in_space}_proseg.qmd   (quarto render)
2d. niche/proseg_output_analysis/banksy_niches_{glm,by_histology}_proseg.qmd
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
  `k_geom_15_30_pc_20_lam_0.2_0.8_k_leiden_30_50_res_0.75_1_20241230_232215`.
  Files with this suffix in `banksy_param_search/` and
  `Xenium_annotations/` are the same cohort run; do not mix with other
  timestamps.
- Generated outputs (`*.html`, `figure-html/`, `libs/`, `._*`) are not
  committed — see `.gitignore`.
