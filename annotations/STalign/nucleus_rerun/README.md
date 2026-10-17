# Nucleus-guided STAlign rerun (H&E ↔ Xenium)

**Branch:** `stalign-nucleus-rerun`, in the worktree `/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code_stalign`.

## Status (2026-10-14): done for all 11 samples

- **New registration:** `DATA/Xenium_annotations/pathology_annotations_101326_nucleus/proseg_pathology_<sample>_proseg_meta.csv`.
  It has the same columns as 030226. `STAlign_x`/`STAlign_y`, `adjust_xy`, `annotation` and
  `region_id` are recomputed; everything else is copied unchanged.
- **Accuracy:** on held-out windows (never used to fit), DAPI nuclei land a median **0.9–1.4 µm** from
  their H&E nuclei in every sample, with p90 ≤ 2.3 µm. The old map was 3.8–85 µm off.
- **Old folder:** `pathology_annotations_030226` is unchanged (md5 of all 11 files checked before and
  after).
- **Annotated BANKSY object (2026-10-16):** `add_annotations.R` now reads the new folder and writes
  `..._20260826_115134_pathology_nucleus.rds` (+ `_subset.rds`); the 030226-based `..._pathology.rds` is kept.
  See "Annotated BANKSY object" below for what was repointed and what still uses the old labels.
- **Write-up with figures:** Notion page "stalign rerun".
- **Original task brief (2026-10-13):** `handoff_2026-10-13.md`.
- **Code:** the pipeline as run is commit `7de415f` on this branch. Later commits add checks and
  documentation only. `NEW_DIR/provenance.json` records the HEAD it was written from and the md5 of
  every script.

### Old vs new pipeline

| | Old (030226; the sample notebooks) | New (101326_nucleus; this folder) |
|---|---|---|
| H&E image (STalign source) | whole-slide JPG, RGB, 1.76 µm/px | SVS hematoxylin channel (colour deconvolution), about 4 then 2 µm/px, one section |
| Xenium image (target) | cell-density map of the ProSeg cell centroids: 10 µm grid, each cell a Gaussian with sd 20 µm | DAPI `morphology_focus` image at the same resolution as the H&E |
| Landmarks | 4–6 placed by hand; 80–373 µm off (per-sample medians) | 20–49 measured from the nucleus pattern of 150 µm windows; affine 0.5–0.9 µm on held-out pairs (GG241670 3.6) |
| Landmark weight σP | 20 µm (STalign default) | 0.5 µm |
| Affine during LDDMM | free (epL, epT at STalign defaults) | frozen (epL = epT = 0) |
| Deformation step epV | 10 | 0.3 |
| Smoothness a | 500 JPG px, about 880 µm | 200 then 150 µm |
| Iterations | 2000 | 150 + 100 |
| Pixel re-weighting | on (STalign always does it, from iteration 50) | on, the same |
| Held-out nucleus error, median | 3.8–85 µm | 0.9–1.4 µm |

- **Why the old map was off:** its step sizes were too large and the pixel re-weighting kicked or trapped
  the map. Its blurred cell-density image would have limited it to about 5 µm in any case. See "Why the
  old map was off".
- **Why the new map is right:** its start is already accurate (nucleus landmarks), stiff landmarks hold
  it there, steps are small with the affine frozen, and nucleus-level images pin it sharply. The
  re-weighting still runs but has little room to move the map. See "Why the new map ends right".

### Reproduce from scratch (for a new user)

Starting from the raw H&E and Xenium data, a new user needs the inputs and software below and about
3 h on 4 cores. `DATA = /mnt/storage/dept/medonc/beroukhim/youyun/plgg/data`; `NEW_DIR` is the
folder above.

**Inputs.** The md5 of each is in `NEW_DIR/provenance_files.tsv` (roles `input_*`).

| Input | Path | Used for |
|---|---|---|
| H&E whole-slide images | `DATA/Xenium_annotations/images/230918_Xenium_CytAssist_{LGG1,LGG2}.svs` | H&E nuclei in each window; the hematoxylin image for LDDMM |
| Xenium DAPI | `DATA/segmentation/proseg_run_121024/<folder>/<folder>_proseg/outs/morphology_focus/morphology_focus_0000.ome.tif` | DAPI nuclei in each window; the DAPI image for LDDMM |
| Xenium folder lookup | `DATA/metadata/Xenium_PS.xlsx` | sample → `<folder>` |
| Old registration (030226) | `DATA/Xenium_annotations/pathology_annotations_030226/proseg_pathology_<s>_proseg_meta.csv` | the cell table (cell ids, Xenium x/y, every column not recomputed) and the starting map (below) |
| Pathology polygons | `DATA/Xenium_annotations/geojsons/<file named in the notebook>` | region labels |
| Sample notebooks | `annotations/STalign/pathology_<s>_proseg.qmd` (this repo) | which GeoJSON each sample uses, and its colour → label `recode_map` |

**The 030226 metas are a fixed input; don't regenerate them.**
- The old notebooks don't reproduce at this scale: the same run ends 26 µm apart on 1 vs 2 CPU threads.
- The pipeline uses the old map only as a lookup: where to read each window's DAPI, and how to turn a
  measured offset into a Xenium point. Each landmark's offset is measured from nuclei, so the old map's
  error does not carry over.
- The result is checked against nuclei on held-out windows, which doesn't depend on the old map.

**Software.** conda env `spatial`; the full export is `NEW_DIR/environment_spatial.yml`.
- Versions that are actually imported (`provenance.json` → `software`): Python 3.10.19, torch 2.10.0
  (conda-forge CPU build, MKL, AVX512 path), numpy 1.23.4, scipy 1.13.1, scikit-image 0.24.0, tifffile
  2025.5.10, imagecodecs 2024.12.30, pandas 2.0.0, geopandas 1.1.2, shapely 2.0.7, Pillow 9.5.0.
- STalign: `pip install git+https://github.com/JEFworks-Lab/STalign@b2068edc98974efa54537eca194736e177bbe11d`.
  Every installed STalign file matches that commit by md5 (`STalign.py` ef918cea…).
- `pip list` shows torch 2.0.0 because of a leftover pip metadata folder; the torch that loads is 2.10.0.
  `provenance.json` files written before 2026-10-15 recorded 2.0.0 by mistake.

**Paths.** The scripts use absolute paths on this cluster: `DATA` and `ROOT` in `nucleus_core.py`, and
`PY` in `run_all.sh`. On another system, edit those three.

**Run.** Use a 4-slot allocation, e.g. `qrsh -pe smp 4 -l h_rt=12:00:00,mem_free=32G` from argos-qsub1.
```bash
git checkout <commit>; conda activate spatial
# registration from the raw inputs into an empty folder; nothing cached is reused (~3 h)
OUT_DIR=<empty dir> bash annotations/STalign/nucleus_rerun/run_all.sh          # PARALLEL=2 THREADS=2
# compare with the production outputs -> <empty dir>/compare_to_production.tsv
python annotations/STalign/nucleus_rerun/provenance.py --compare <empty dir>
```
- Expect identical meta files with THREADS=2 on this hardware. With a different thread count, cell
  positions differ by ≤ 0.002 µm (PA320784 test below). Other CPUs are untested.
- Without `OUT_DIR`, `run_all.sh` writes to `NEW_DIR` itself and reuses the cached old-map window
  tables there.
- Each sample log prints one `ERROR 1: PROJ: proj_create_from_database` line, as the production logs do.
  It comes from reading the GeoJSON; the polygons are in pixel units, so no projection is used.
- **Verification (2026-10-15, done):** all 11 samples were rebuilt this way into
  `NEW_DIR/qc/repro_from_scratch/` (about 3 h; two groups of PARALLEL=1 THREADS=2). All 11 metas are
  byte-identical to production. So are the recomputed old-map window tables, the landmark/held-out
  roles, the after-registration window tables and the transforms (`compare_to_production.tsv`: every
  difference 0).

### Annotated BANKSY object (2026-10-16)

The new metas reach the analyses through one object. `annotations/STalign/add_annotations.R` copies each cell's
`annotation` (recoded) and `region_id` from the per-sample metas into the colData of the final BANKSY
SpatialExperiment (`banksy_param_search/.../..._20260826_115134.rds`). Only labels are copied, not coordinates.

- **Now:** it reads `pathology_annotations_101326_nucleus` and writes, in `DATA/Xenium_annotations/`:
  - `banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134_pathology_nucleus.rds`
    (md5 `9f0cfb6f…`, 1,802,005 cells);
  - `..._pathology_nucleus_subset.rds` (md5 `05044d4f…`, 90,100 cells, `set.seed(1)`; the 030226 subset was
    not seeded).
  - The 030226-based `..._pathology.rds` and `_subset.rds` are kept, unchanged.
- **Check** (`check_annotated_object.R` → `qc/annotated_object_check.tsv`): every cell's label in the object
  equals its meta row, and the subset carries the same labels. Relabelled against the old object, per
  sample: 0.1–16% of all cells, and 0.2–75% of the cells the old map put in an annotated region. That matches
  `qc/label_impact.tsv`; the counts are a little lower because the object holds only QC-passing cells.
- **Repointed to the new object** (not re-rendered): `annotations/linear_models.qmd`,
  `annotations/linear_models_niches.qmd`, `niche/proseg_output_analysis/niche_interpretation_proseg.qmd`.
  `linear_models_niches.qmd` reads the subset from a local Mac path when run on the Mac, so copy the new
  subset there first.
- **Still on the old labels or coordinates** (not in this branch; untracked or uncommitted in the main checkout):
  - `sc_integration/ENVI/2_spatial_imputation/PA/export_new_object_cell_metadata.R`, which writes
    `xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826.tsv` from the old object. The PA MAPK and
    morphology notebooks read that TSV.
  - The CL lipid / vacuole scripts under `sc_integration/ENVI/2_spatial_imputation/CL/`, which read
    `STAlign_x/_y` from the 030226 metas directly.
  - `XENIUM_RUNBOOK.md`, which still names the 030226 folder and `..._pathology.rds`.
- **Not affected:** the ENVI training h5ads (`sce_2_anndata.qmd` reads the unlabelled BANKSY object).

### All commands (production outputs, figures, diagnostics)

```bash
# 1. registration, all 11 samples (~2.3 h on 4 slots). Per sample: old-map window check (cached in
#    qc/windows_before_<s>.tsv) -> nucleus landmarks -> affine -> 2-stage LDDMM -> new meta + validation
bash annotations/STalign/nucleus_rerun/run_all.sh                 # PARALLEL=2 THREADS=2 by default
#    affine-only reference variant (written to qc/trial_affine_all/)
python annotations/STalign/nucleus_rerun/stalign_nucleus.py <sample> --niter1 0 --niter2 0 --tag affine_all
# 2. label changes, figures
python annotations/STalign/nucleus_rerun/label_impact.py   $NEW_DIR $NEW_DIR/qc $NEW_DIR/qc/trial_affine_all
python annotations/STalign/nucleus_rerun/cohort_figures.py $NEW_DIR $NEW_DIR/qc $NEW_DIR/figures $NEW_DIR/qc/trial_affine_all
# 3. diagnostics behind the explanations (outputs in qc/detection, qc/landmarks, qc/replay_old, qc/frame, qc/image_term)
python annotations/STalign/nucleus_rerun/diag_detection.py         # one-to-one H&E/DAPI peak accounting (~15 min)
python annotations/STalign/nucleus_rerun/diag_detection.py --threshold-area   # DAPI threshold from the window only
python annotations/STalign/nucleus_rerun/diag_jpg_frame.py         # JPG vs SVS, polygon outlines vs tissue edges
python annotations/STalign/nucleus_rerun/stalign_nucleus.py PA320784 --threads 1 --tag repro_1thread   # determinism
python annotations/STalign/nucleus_rerun/diag_old_landmarks.py     # the notebooks' hand landmarks vs truth
for s in PA320784 CL42112 GG240468; do for l in original true true_grid; do   # old pipeline replays, ~6-11 min each
    python annotations/STalign/nucleus_rerun/replay_old_stalign.py $s $l --threads 2; done; done
python annotations/STalign/nucleus_rerun/replay_old_stalign.py PA320784 original --threads 1 --tag 1thread
# step-size tests (1 thread each, OMP_NUM_THREADS=1): notebook steps with the affine recorded, affine frozen, 10x smaller
for s in PA320784 CL42112; do
    python annotations/STalign/nucleus_rerun/replay_old_stalign.py $s true --threads 1 --tag traj_1thread
    for l in true original; do python annotations/STalign/nucleus_rerun/replay_old_stalign.py $s $l --threads 1 \
        --epL 2e-9 --epT 2e-2 --epV 1 --tag smallstep_1thread
        # the same with STalign's pixel re-weighting switched off
        python annotations/STalign/nucleus_rerun/replay_old_stalign.py $s $l --threads 1 --epL 2e-9 --epT 2e-2 --epV 1 \
            --reweight-start 100000 --tag smallstep_noreweight_1thread; done; done
for s in PA320784 CL42112 GG240468; do
    python annotations/STalign/nucleus_rerun/replay_old_stalign.py $s true --threads 1 --epL 0 --epT 0 --tag frozen_affine_1thread; done
python annotations/STalign/nucleus_rerun/replay_old_stalign.py --figure
python annotations/STalign/nucleus_rerun/replay_old_stalign.py --trajectory
python annotations/STalign/nucleus_rerun/replay_old_stalign.py --conditions   # end error by start and setting
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python annotations/STalign/nucleus_rerun/diag_image_term.py <sample> [...]   # ~8 min per sample
python annotations/STalign/nucleus_rerun/diag_image_term.py --figure
# 4. from-scratch rerun and comparison (above)
OUT_DIR=$NEW_DIR/qc/repro_from_scratch bash annotations/STalign/nucleus_rerun/run_all.sh
python annotations/STalign/nucleus_rerun/provenance.py --compare $NEW_DIR/qc/repro_from_scratch
# 5. the annotated BANKSY object from the new metas (~16 min, 14 GB RAM), its check, then the provenance manifest
Rscript annotations/STalign/add_annotations.R
Rscript annotations/STalign/nucleus_rerun/check_annotated_object.R
python annotations/STalign/nucleus_rerun/provenance.py            # records HEAD, md5s, versions, conda env export
```
(GG240468 was replayed with `original` and `true` only.)

- `NEW_DIR/qc/transform_<s>.npz` stores the fitted transform: affine `A`, velocity field `v` on grid
  `xv0`/`xv1`, and the run arguments.
- `NEW_DIR/provenance.json` records the commit, script md5s, software versions, constants, seeds and
  per-sample arguments. `provenance_files.tsv` lists the md5 of every input and output.
- Windows are drawn with seed 0 and split into landmark and held-out halves with seed 1.

| Script | Role |
|---|---|
| `nucleus_core.py` | Nucleus detection, FFT shift search and image readers. Vendored verbatim from the CL lipid session (`he_image.py`, `registration_nuclei_check.py`, `lipid_detector_v2.py`; source md5s in its header) |
| `nucleus_windows.py` | Window drawing and the per-window nucleus check, for any meta file |
| `stalign_nucleus.py` | The pipeline for one sample |
| `run_all.sh` | Runs all samples |
| `label_impact.py` | Per-cell displacement and label changes, old vs new |
| `cohort_figures.py` | Offsets dumbbell, error field, per-slide overlays (H&E: SendUserFile / Notion only) |
| `overlay_before_after.py` | Three held-out windows for one sample |
| `diag_detection.py` | One-to-one accounting of H&E vs DAPI peaks, peak strength, cross-stain signal, example windows |
| `diag_old_landmarks.py` | The original notebooks' hand landmarks vs their true positions; figures |
| `replay_old_stalign.py` | Reruns the original notebook pipeline with hand, exact or grid-rounded landmarks, optionally with other step sizes or the pixel re-weighting off; records the affine per iteration |
| `diag_jpg_frame.py` | Checks that the JPG is the SVS at level 0 / 7, and the polygon outlines against tissue edges |
| `diag_image_term.py` | How sharply each H&E/Xenium image pair's match pins the map (R² of STAlign's contrast fit vs offset) |
| `check_annotated_object.R` | Checks the annotated BANKSY object against the new metas; label changes vs the old object |
| `provenance.py` | Provenance manifest; `--compare` checks a rerun against the production outputs |

### Method (final)

1. **Old-map check.** For each sample, draw 100 windows of 150 µm at random cells: at least 100 cells
   each, at least 250 µm apart. Measure the old map's nucleus offset in each window.
2. **Landmarks and held-out set.** Confident windows (peak z ≥ 5) are split at random into halves.
   - In a landmark window, the H&E centre is paired with the Xenium point the old map puts there,
     corrected by the measured offset.
   - A robust affine fit (`L_T_from_points` plus outlier rejection) gives the start.
   - See "Landmarks used by this run" below.
3. **LDDMM**, with the affine held fixed (epL = epT = 0) so only the velocity field adapts.
   - Images: H&E hematoxylin from the SVS against Xenium DAPI (log), both scaled 0–1.
   - Stage 1: about 4 µm/px, a = 200 µm, 150 iterations.
   - Stage 2: about 2 µm/px, a = 150 µm, 100 iterations.
   - Other settings: epV 0.3, sigmaM 0.5, sigmaP 0.5 (landmark term kept), sigmaR 5e5, nt 3, float32
     on CPU.
4. **Cells.** Map each cell to H&E and relabel it with the notebook's own rule: the first polygon in
   GeoJSON order that contains the cell, its colour recoded with the notebook's `recode_map`.
5. **Validation.** Re-measure the same windows through the new map. Report held-out windows
   separately from landmark windows.

**Why LDDMM rather than affine only (user decision 2026-10-14).**
- LDDMM keeps the original method and allows local deformation.
- Affine alone is slightly tighter in 10/11 samples: held-out median 0.56–0.90 µm against 0.90–1.41
  µm for LDDMM.
- LDDMM is clearly better in GG241670: 1.00 µm (p90 1.53) against 2.25 µm (p90 4.28). There, 23
  windows the old map could not measure sit 2.5–5.6 µm off under affine and about 1 µm under LDDMM.
- The two variants assign a different label to ≤ 0.6% of cells in every sample.

### Landmarks used by this run

No landmark is placed by hand. `stalign_nucleus.py` finds them in the tissue:
1. Up to 100 window centres are drawn at random cells (seed 0). Each window is 150 µm, holds at least
   100 cells, and is at least 250 µm from the others, so small sections get fewer.
2. In each window, nuclei are called in the H&E and in the DAPI (next section). The DAPI nuclei are
   placed on the H&E through the **old** map. An FFT search over ±150 µm finds the offset d that best
   overlays the two point patterns.
3. A window is confident when the correlation peak has z ≥ 5. A random half of the confident windows
   (seed 1) become landmarks; the other half are held out and used only for scoring.
4. Each landmark pairs the H&E window centre c with the Xenium point that truly lies under it: the
   old map's preimage of c + d. The old map only serves as a lookup between the two images, and its
   error there is exactly what d measures.
   - **One window gives one landmark, not one per nucleus.** The landmark is a point pair at the
     window centre. Its offset d comes from the whole nucleus pattern of the window: a median 78–203
     H&E nuclei per landmark window (by sample; minimum 28), overlaid on the DAPI nuclei.
   - **Why not individual nuclei.** Single-nucleus pairs are ambiguous: the stains call peaks
     1.17–1.57× apart in number, and a mapped nucleus can sit closer to a neighbour than to itself. A
     pattern of about 100 nuclei has one clear best overlay. That overlay is found to 0.25 µm, and its
     z-score says whether it is real.
   - **Assumption:** within a 150 µm window, the old map's error is a pure shift. After that shift, matched
     nuclei sit a median 1.1–1.7 µm apart (per-sample medians of `residual_um` in landmark windows). That
     bounds what the assumption plus nucleus-centre noise can cost.
5. `robust_affine` fits L, T and drops pairs whose residual exceeds max(15 µm, median + 3 robust SD).
   Only GG241670 lost any (5).
6. LDDMM keeps every remaining pair in its landmark term, with σP = 0.5 µm.

| Sample | Windows | Confident | Landmarks (rejected) | Held out | Affine residual at landmarks, median µm | Affine on held-out pairs, median µm |
|---|---|---|---|---|---|---|
| CL42112 | 100 | 90 | 45 | 45 | 0.9 | 0.7 |
| CL67068 | 100 | 66 | 33 | 33 | 1.3 | 0.8 |
| CL87352 | 100 | 83 | 41 | 42 | 0.6 | 0.6 |
| CN51420 | 100 | 73 | 36 | 37 | 0.6 | 0.8 |
| CN53398 | 100 | 55 | 27 | 28 | 0.6 | 0.8 |
| GG240468 | 54 | 53 | 26 | 27 | 0.7 | 0.6 |
| GG241670 | 100 | 77 | 33 (5) | 39 | 3.3 | 3.6 |
| PA188570 | 77 | 75 | 37 | 38 | 0.6 | 0.7 |
| PA258482 | 100 | 98 | 49 | 49 | 0.7 | 0.6 |
| PA320784 | 40 | 40 | 20 | 20 | 0.7 | 0.8 |
| PA328368 | 90 | 88 | 44 | 44 | 0.6 | 0.5 |

(From `qc/run_<sample>.log`; window roles are in `qc/windows_before_roles_<sample>.tsv`.)
- Where the old map was more than 150 µm off, no confident shift exists, so those windows never become
  landmarks (for example the 200–600 µm corner of CL67068). The new map is re-measured there as an
  independent check (Results).
- "Affine on held-out pairs" scores the held-out pairs measured on the old map. GG241670 includes some
  false correlation peaks there. The definitive score is the re-measurement on the new map (Results).

### Frames: SVS, JPG and polygons

- **The polygons are in SVS level-0 px.** Their coordinates span the SVS level 0: up to 85,866 ×
  43,370 px on LGG1, whose level 0 is 89,640 × 46,004. The JPG is 12,806 px wide.
- **The only polygon step is a fixed 1/7 scaling.** The notebooks scale the polygons by 1/7 into the
  JPG frame (`s.geometry.scale(xfact=1/7, yfact=1/7, origin=(0,0))`) and label each cell by the
  polygon containing its `STAlign_x/_y`. The polygons are never moved; STAlign is applied to cells.
- **1/7 is exact.** The JPG is the SVS at level 0 / 7. At 66 patches (6 per section), phase correlation
  of the JPG against the SVS gives offsets ≤ 0.73 µm (`diag_jpg_frame.py`, `qc/frame/jpg_vs_svs.tsv`).
- **The new map never uses the JPG.** It is fitted in SVS level-0 µm, from the SVS hematoxylin. Cell
  positions are written as level 0 / 7 to keep the `STAlign_x/_y` columns, and labels use polygons / 7,
  so the containment test is effectively done in SVS level 0.
- **Polygon outlines vs tissue edges** (`qc/frame/polygon_edges.tsv`). The test shifts the outlines
  and finds where they best follow the JPG's tissue edges.
  - In PA188570, PA320784, PA328368, CL42112 and CL87352, 66–92% of outline points lie within 26 µm
    of an edge. There the best shift is 4–12 µm and gains ≤ 0.11 px of mean distance, so it is flat
    around zero.
  - Samples that share a GeoJSON file don't share an offset, which a frame error would cause.
  - In GG240468, PA258482 and CN51420, few outline points follow an edge (18–36%), so the test is not
    informative there.
  - GG241670's six polygons (110 annotated cells) fit best 40 µm away. PA188570 and PA328368, from
    the same file, don't. So this reflects how those outlines were drawn, not a frame shift.
  - How precisely hand-drawn outlines follow the tissue is a property of the annotation, separate
    from the registration.

### How nuclei are called in the window check (exact)

This is `nucleus_windows.measure()`, with `peaks()` and `best_shift()` from `nucleus_core.py`.

| Step | H&E | Xenium DAPI |
|---|---|---|
| Image | SVS full resolution (0.2522 µm/px), averaged in 2×2 blocks → RGB at 0.504 µm/px | `morphology_focus` full resolution (0.2125 µm/px), raw intensity |
| Area read | The 150 µm window | The Xenium area the map under test puts under the window, padded by 150 µm per side (about 450 µm square), because the shift search needs it |
| Channel | Hematoxylin optical density from colour deconvolution (`skimage.color.rgb2hed`, channel 0) | The DAPI image itself |
| Smoothing | Gaussian, σ = 1 µm (1.98 px) | Gaussian, σ = 1 µm (4.7 px) |
| Threshold | Otsu threshold of the smoothed window, all pixels including glass | Otsu threshold of the smoothed, padded read area |
| Peaks | `peak_local_max(min_distance=7 px)`: pixels above the threshold that are the maximum of the (2·7+1)² px square around them. Two peaks both survive only if they are more than 7 px apart along row or column, i.e. ≥ 8 px = **4.03 µm** | Same rule, `min_distance=18 px`: ≥ 19 px = **4.04 µm** |

Why the two columns differ:
- **Peak spacing: they don't, in effect.** Both use `NUC_MIN_DIST_UM = 4.0`, converted to whole
  pixels with `int()`: 7.94 → 7 H&E px and 18.8 → 18 DAPI px. Because a peak must beat every pixel
  within that many px, survivors are ≥ 8 px (4.03 µm) and ≥ 19 px (4.04 µm) apart. Only the pixel
  grid differs.
- **Area read: the DAPI side is the search area.** The H&E window is the pattern being looked for.
  The DAPI read must hold every place it could be, so it adds the ±150 µm search range on each side.
  - The side effect is a DAPI threshold taken over about 9× the area.
  - `diag_detection.py --threshold-area` reran the accounting with the DAPI threshold taken from the
    window alone (`qc/detection/threshold_area.tsv`). That threshold is 0.5–12% higher (sample
    medians).
  - DAPI peak counts change by −2.0% to +0.3%, the H&E/DAPI ratio goes from 1.17–1.57 to 1.17–1.60,
    and the H&E paired share moves by ≤ 1 point. So the area difference does not explain the H&E
    surplus.

Then:
- DAPI peaks are mapped into the H&E window through the map under test.
- "Matched" = the share of H&E peaks whose **nearest** DAPI peak is ≤ 2.5 µm away. This is not one-to-one: two H&E peaks can share one DAPI peak.
- `diag_detection.py` does the one-to-one accounting. Each peak is either paired (mutual nearest, ≤ 2.5 µm), a duplicate, a near miss (2.5–6 µm), or has no partner within 6 µm. Its outputs are in `qc/detection/`.

### Results (held-out windows unless stated)

| Sample | Old median (p90) µm | New median (p90) µm | Affine only | Nuclei matched old → new [ceiling] | Old-unmeasurable windows now confident (median µm) | Cell shift old → new, median / p90 µm | Annotated cells relabelled (n annotated) | All cells relabelled |
|---|---|---|---|---|---|---|---|---|
| CL42112 | 78.0 (96.5) | 1.06 (1.65) | 0.71 | 0.08 → 0.45 [0.49] | 8/10 (0.75) | 77 / 99 | 43% (21,943) | 11.7% |
| CL67068 | 30.2 (56.1) | 1.41 (2.00) | 0.56 | 0.11 → 0.59 [0.68] | 33/34 (1.06) | 44 / 308 | 50% (38,264) | 9.8% |
| CL87352 | 29.8 (52.8) | 0.97 (1.35) | 0.56 | 0.10 → 0.50 [0.52] | 12/17 (0.67) | 29 / 62 | 25% (29,161) | 1.8% |
| CN51420 | 52.7 (58.3) | 1.35 (2.26) | 0.75 | 0.11 → 0.61 [0.68] | 23/27 (1.46) | 56 / 114 | 75% (8,077) | 5.2% |
| CN53398 | 34.2 (44.3) | 0.97 (1.48) | 0.71 | 0.12 → 0.56 [0.57] | 32/45 (1.03) | 34 / 55 | 59% (1,002) | 0.7% |
| GG240468 (control) | 3.8 (7.1) | 1.06 (1.55) | 0.56 | 0.07 → 0.72 [0.74] | 1/1 (0.75) | 4 / 9 | 0.2% (49,204) | 0.4% |
| GG241670 | 85.1 (143.1) | 1.00 (1.53) | 2.25 | 0.05 → 0.70 [0.74] | 23/23 (0.90) | 91 / 188 | 67% (110) | 0.1% |
| PA188570 | 34.0 (40.3) | 1.18 (1.59) | 0.77 | 0.08 → 0.70 [0.76] | 0/2 (too few nuclei) | 34 / 42 | 18% (19,900) | 16.1% |
| PA258482 | 16.6 (25.0) | 0.90 (1.70) | 0.75 | 0.07 → 0.68 [0.72] | 1/2 (6.3, torn tissue edge) | 15 / 24 | 2.5% (23,029) | 1.3% |
| PA320784 | 28.1 (31.8) | 1.06 (1.58) | 0.90 | 0.05 → 0.68 [0.73] | 0/0 | 25 / 31 | 8.5% (11,262) | 6.8% |
| PA328368 | 33.2 (43.7) | 1.03 (1.56) | 0.56 | 0.06 → 0.59 [0.63] | 2/2 (1.35) | 32 / 46 | 5.5% (48,174) | 7.9% |

How to read the columns:
- **Nuclei matched** is the share of H&E nuclei that have a DAPI nucleus within 2.5 µm.
  - The **ceiling** is the same share at each window's best shift, so it is the most any map could
    reach.
  - The ceiling sits below 1 because the hematoxylin channel calls 1.17–1.57× more peaks than DAPI.
    `diag_detection.py` puts every peak into exactly one class, so each modality sums to 100% and both
    sides count the same pairs:

    | | H&E peaks | DAPI peaks |
    |---|---|---|
    | Paired one-to-one (mutual nearest, ≤ 2.5 µm) | 47–74% | 74–88% |
    | Near miss (nearest other-stain peak 2.5–6 µm) | 10–22% | 9–22% |
    | No partner within 6 µm | 13–33% | 1.6–7% |
    | Duplicate (two peaks on one partner) | ≤ 2 peaks per sample | 0–1 |

  - **Unpartnered peaks are weaker than paired ones, in both stains.** Pooled medians are 1.62 vs
    2.19× the H&E threshold, and 1.71 vs 2.23× the DAPI threshold. They are the borderline calls.
  - **Near misses are the same nucleus with a different centre.** The DAPI under a near-miss H&E peak
    is 0.67–0.91 of a matched nucleus's.
  - **H&E peaks with no partner are a mix, and the mix depends on the tissue.** Measured as DAPI
    relative to a matched nucleus in the same window:
    - In PA/GG sections it is 0.04–0.61 (sample medians). The DAPI image is dark there, so these are
      mostly hematoxylin staining with little DNA.
    - In the dense CL/CN sections it is 0.58–0.77, and many look like nuclei in the H&E. The DAPI there
      is dim or diffuse, so no separate peak is called. Plausible causes: nuclei dim or out of focus in
      the single 2D DAPI image, or neighbours merged by the 1 µm smoothing and 4 µm peak spacing.
      This is not resolved.
  - These extra peaks lower the match share but not the shift, because the shift is the best overlay
    of the whole point pattern.
  - *Correction:* an earlier version said "most extra H&E peaks are not nuclei". That holds for PA/GG,
    not for CL/CN.
- **Old-unmeasurable windows** gave no confident shift on the old map (z < 5).
  - Usually this is because the old map was more than the ±150 µm search range off there. For
    example, a 200–600 µm region in CL67068 holds both cerebellum layers and most of the infiltration
    zone.
  - These windows were never landmarks, so their new-map offsets are an independent check.
  - It also means the "old median" column **understates** the old error in CL67068, GG241670, CN51420,
    CN53398 and CL87352. The cell-shift column covers every cell.
- **Annotated cells relabelled** counts cells the old map placed inside a polygon whose label differs
  under the new map. GG241670 (110 cells) and CN53398 (1,002 cells) have very small denominators.
- **Per-label retention** is in `qc/label_retention.tsv`. Examples of cells keeping their label, with
  the share of the new population that are new entrants:
  - Lipidized cells: CL42112 57% (33%), CL67068 71% (29%), CL87352 75% (20%).
  - CN51420 perivascular pseudorosettes: 24% kept.
  - CL67068 cerebellum molecular / granule layer: 7% / 13% kept.

### Why the new map ends right, and how reproducible it is

- **The accuracy comes from the landmarks.** There are 20–49 per sample, each from the nucleus
  pattern of a 150 µm window. The affine fitted to them alone is already 0.5–0.9 µm on held-out
  pairs (3.6 µm in GG241670; table above). The old notebooks started from 4–6 hand landmarks that
  were 80–373 µm off.
- **The LDDMM is kept from drifting**, which the old one did even from exact landmarks (see "Why the old
  map was off"):
  - the affine is held fixed (epL = epT = 0);
  - the step is small (epV 0.3, against 10);
  - the landmark term is 40× tighter (σP 0.5 µm, against 20 µm), over 20–49 landmarks instead of 4–6;
  - deformation is local (a = 150–200 µm, against about 880 µm);
  - it runs 250 iterations, against 2000.
- **The pixel re-weighting is not switched off here.** `stalign_nucleus.py` calls STalign's standard
  LDDMM, which re-estimates the weights from iteration 50 of each stage. The stiff landmarks and the
  accurate start leave it little room.
- **Its image term still pulls a little.** In 10 of 11 samples, held-out medians are 0.90–1.41 µm
  with LDDMM against 0.56–0.90 µm for the affine alone. The landmarks hold this under 1 µm. Whether the
  re-weighting causes this pull is untested. GG241670 is the exception: LDDMM takes it from 2.25 to
  1.00 µm.
- **Reproducibility (PA320784 rerun from the committed code):**
  - Same thread count (2): byte-identical meta (`qc/trial_repro/`).
  - 1 thread instead of 2 (`qc/trial_repro_1thread/`): cell positions differ by at most 0.002 µm
    (median 0.0004 µm), no label changes, identical validation. The md5 differs only through those
    last digits.
  - The old pipeline under the same test moved 26 µm.
- **Code vs outputs.** `nucleus_core.py` and `nucleus_windows.py` were made self-contained (the
  measurement code copied in verbatim) while the full run was in progress. The rerun from that
  version gave a byte-identical meta.

### Known quirks

- **CL42112 window 72** (old-unmeasurable) reads 204 µm at z 6.2 under LDDMM. Its zero-shift match
  share is 0.30, the same as under affine (0.56 µm there). It is a false correlation peak in a sparse
  window (44 H&E nuclei), not a misregistration.
- **PA258482 window 4** reads 6.3 µm, the same under both variants. It sits on a section edge where
  DAPI nuclei extend past the torn or retracted H&E tissue.
- **GG241670 colour `#800000`** has no entry in the notebook's `recode_map`, so 79 cells kept the raw
  hex label in 030226. The rerun reproduces this as it is.
- **LDDMM energy** can rise between iterations because STalign re-estimates its artefact/background
  weights. Accuracy is judged by the held-out windows, not by the energy.
- **GG240468** started at 3.8 µm and ends at 1.06 µm. As the no-change control it confirms the new
  pipeline keeps a near-correct map correct: only 0.2% of annotated cells change label.

### Why the old map was off (tested 2026-10-14 to 2026-10-16)

**In short.** Two settings of the old LDDMM caused the tens of µm.
- Its step sizes were so large that the map never settled.
- STalign's pixel re-weighting, from iteration 50 on, kicks the map and can trap it in a wrong place.

With both changed, the old images and even the old hand landmarks end 4.5–4.9 µm from the truth (2
samples tested). The hand landmarks were far off but hardly mattered. The old images point to the right
place, but they are too blurred to pin the map better than about 5 µm, so they hold it only weakly.

**The old pipeline** (`annotations/STalign/pathology_<s>_proseg.qmd`; see "Old vs new pipeline" above):
- H&E side (STalign's source): the whole-slide JPG, RGB, 1.76 µm/px.
- Xenium side (target): a cell-density image of the ProSeg cell centroids (`cells.csv.gz`). Each centroid
  is spread as a Gaussian with sd 20 µm on a 10 µm grid (`STalign.rasterize(dx=10)`). No DAPI is used.
- Start: an affine from 4–6 hand-placed landmarks (`L_T_from_points`).
- LDDMM settings:
  - 2000 iterations; a = 500 JPG px (about 880 µm); sigmaM 0.15, sigmaB 0.10, sigmaA 0.11.
  - Step sizes: epV 10, with epL / epT at STalign's defaults (2e-8 / 0.2, which LDDMM divides by 10).
  - Landmark term on, with σP = 20 µm (default). Pixel re-weighting on (STalign always does it).

**How it was tested.** `replay_old_stalign.py` reruns that pipeline exactly and changes one thing at a
time, scoring each run by every cell's distance from the new map.
- Landmarks: the notebook's hand landmarks, or "exact" ones (the same Xenium points, with their true H&E
  positions from the new map).
- Step sizes: notebook, affine frozen (epL = epT = 0), or all 10× smaller (epL 2e-9, epT 2e-2, epV 1).
- Pixel re-weighting: on, or off (`--reweight-start 100000`, so every pixel keeps weight 0.5).

Median cell error at the end, µm (`qc/replay_old/replay_summary.tsv`, `fig_replay_conditions.png`):

| Sample | Start: affine from landmarks | Notebook step sizes | Affine frozen | All steps 10× smaller | 10× smaller, re-weighting off | Stored 030226 map |
|---|---|---|---|---|---|---|
| CL42112 | hand, 70 µm | 65.2 | | 4.2 | **4.9** | 76.8 |
| CL42112 | exact, 0.8 µm | 78.1 (2 threads), 34.0 (1 thread) | 55.9 | 5.9 | **4.5** | |
| PA320784 | hand, 243 µm | 10.4 (2 threads), 22.8 (1 thread) | | 124 | **4.5** | 25.3 |
| PA320784 | exact, 0.2 µm | 32.1 (2 threads), 51.5 (1 thread) | 6.3 | 17.5 | **4.7** | |
| GG240468 | hand, 115 µm | 5.4 | | | | 4.4 |
| GG240468 | exact, 0.2 µm | 3.5 | 7.0 | | | |

(Runs tagged `_1thread` used 1 CPU thread. Exact landmarks rounded to the notebooks' 50 grid: PA320784
43.5 → 19.0 µm, CL42112 47.0 → 78.7 µm.)

**1. The step sizes were too large, so the map never settles** (`fig_replay_trajectory.png`).
- Each iteration, STalign moves each part of the map by step size × the slope of its energy: epT for
  the shift, epL for rotation and stretch, epV for the deformation. The step sizes are unitless
  multipliers. How far the map actually moves depends on how steep the energy is for these images, so
  defaults tuned on other data need not fit.
- From exact landmarks, the map jumps from the first iteration. The landmarks move a median 18–27 µm
  per iteration in PA320784 and CL42112. The affine alone moves the PA320784 section a median 34 µm
  per iteration, up to 134 µm.
- Meanwhile the image term gets worse, rising 37% (PA320784) and 3× (CL42112) within 10 iterations. The
  optimizer is overshooting, not finding a better fit.
- It bounces like this for all 2000 iterations, so the final map is wherever iteration 2000 lands.
  That is why the same run ends 26 µm apart on 1 vs 2 threads, and why exact landmarks end 32–78 µm off.
- GG240468 bounces least (median 5.8 µm per iteration). That, not its landmarks, is why it came out
  right.
- Both parts of the map bounce. Freezing the affine takes PA320784 from 51.5 to 6.3 µm. But in CL42112
  the deformation alone still moves the landmarks a median 12.6 µm per iteration and ends 56 µm off.

**2. STalign's pixel re-weighting kicks the map and can trap it.**
- *What it is.* STalign assumes some pixels of the target image should not be matched: background, and
  artefacts with no counterpart in the other image. For each pixel it estimates the chance of being
  matched tissue, artefact or background (WM / WA / WB). Only WM enters the image term.
  - Until iteration 50 every pixel has WM = 0.5.
  - From iteration 50 on, every 5 iterations, the weights are re-estimated from the current overlay.
    Pixels the overlay explains badly are down-weighted.
- *Tested by switching it off*, with 10× smaller steps. Off, the four runs (2 samples × hand / exact
  landmarks) all end at 4.5–4.9 µm. On, they end at 4.2, 5.9, 17.5 and 124 µm.
- *The kick.* With it on, the landmarks jump at the second re-estimate (iteration 55): PA320784 from 2.3
  to 23 µm, CL42112 from 7.3 to 34 µm. With it off there is no jump (3.6 and 7.6 µm).
- *The trap.* PA320784 starts from its hand landmarks 230 µm off. In the first 50 iterations the affine
  brings the section to 120 µm off, either way.
  - With the re-weighting on, it stays 126–129 µm off for the remaining 1950 iterations. Presumably the
    pixels that don't fit are written off as artefact or background, so the wrong position stops being
    penalised (the weights themselves were not inspected).
  - With it off, the cells end 4.5 µm off.

**3. Nothing held the map at the right place.**
- *The landmarks are slack.* With σP = 20 µm and 4–5 landmarks, the landmark term is at most 1.3% of
  the energy at any iteration of any replay.
  - That was harmless for the hand landmarks, which were wrong anyway (below).
  - Several of the notebooks' Xenium landmarks sit at the section edge (CL42112: 67 and 75 µm from the
    nearest cell; PA320784: 81 µm). There nothing constrains the map. With the re-weighting off they
    drift 13–36 µm while the cells end 4.5 µm off. So the landmark distance in
    `fig_replay_trajectory.png` overstates the cell error of the small-step runs; the table is the cell
    error.
- *The old images point to the right place, but only weakly* (`diag_image_term.py`,
  `fig_image_term.png`). The scan shifts the map rigidly off the true position and computes the R² of
  STAlign's contrast fit over tissue pixels.
  - The old pair's best match is 0–4.5 µm from the truth in all 11 samples, while the stored map is
    3.5–89 µm off.
  - But 5 µm off costs only 0.7–3.4% of the match (10 µm: 1.7–6.8%). For the new hematoxylin–DAPI
    pair it costs 30–78%.
  - A shallow bowl pulls back weakly, so it can't stop the bouncing or the re-weighting.
- *The shallowness comes from the 20 µm blur, not from cell counts vs DAPI.* DAPI blurred the same way
  is just as shallow (5 µm costs 0.6–2.2%). The same DAPI at 10 µm pixels without the blur loses 7–41%.
- *So about 5 µm is the limit of the old images.* With the optimizer fixed they reach 4.5–4.9 µm
  (p90 10–19 µm). About 1 µm needs nucleus-level images and many tight landmarks, as in the rerun.

**4. The hand landmarks were far off, but not what set the final error** (`diag_old_landmarks.py`,
`qc/landmarks/`).
- There are 48 landmarks in total. Placed on the new map, they are off by a median 189 µm (29–630 µm).
  Per-sample medians are 80–373 µm.
- It isn't the 50 grid.
  - In 9 of 11 notebooks every coordinate is a multiple of 50: JPG px on the H&E side (88 µm) and µm
    on the Xenium side.
  - Rounding a perfectly placed pair to that grid moves it at most 97 µm, typically 37 µm (simulated).
    Yet 77% of landmarks are off by more than 97 µm.
  - The grid point nearest each true position was only 4–56 µm away.
  - The two notebooks not on the grid (CL87352, GG240468) were as far off (187 and 178 µm).
  - So the landmarks were misidentified features, not rounded ones.
- They set only the start. From exact landmarks the old pipeline still ends 32–78 µm off, and with the
  optimizer fixed the hand landmarks end 4.5–4.9 µm off.

**5. Ruled out.**
- JPG resolution (1.76 µm/px).
- A JPG/SVS mismatch. The JPG is the SVS at level 0 / 7 to within 0.73 µm at 66 patches across all 11
  sections (`diag_jpg_frame.py`, `qc/frame/jpg_vs_svs.tsv`), so "exact" landmarks in the replays are
  exact in the JPG the old pipeline registered.

**Corrections to earlier versions of this explanation.**
- The 2026-10-13 brief (`handoff_2026-10-13.md`) said the old setup "can't resolve anything under
  about 100 µm" and has "no nucleus-level features to lock onto". The second part is right. But the
  old images do point to the truth; the optimizer settings caused the error.
- Earlier versions said the old image match "is not best at the true alignment", and that each section
  has its own old-pipeline error whatever the start. Both are wrong. The match peaks at the truth, and
  the error changes from run to run (CL42112 from exact landmarks: 34 or 78 µm).
- The 2026-10-15 version said small steps "can't cover a far start" and that the re-weighting "adds
  6–18 µm". The re-weighting-off test (2026-10-16) shows that the far start (PA320784) was stuck because
  of the re-weighting, not the step size.
- Before 2026-10-16, landmark distances derived from the landmark energy were multiplied by 1.76, as if
  σP were 20 JPG px (35 µm). LDDMM measures landmark distances in the Xenium frame, so σP is 20 µm.
  The numbers were 1.76× too large (for example "31–48 µm per iteration", now 18–27).
  `replay_old_stalign.py --trajectory` is fixed, and the figure and numbers above are corrected.

### Why the 5 µm target

- Neighbouring DAPI nuclei sit a median 8.7–10.3 µm apart. Beyond about half that spacing, a mapped
  point is closer to the next nucleus than its own, so per-cell H&E features attach to the wrong cell.
- Region labels tolerate more, as long as errors stay well under region width. The annotated regions
  here are often only 40–140 µm wide (CN51420 pseudorosettes 43 µm, CL42112 lipidized islands 135 µm,
  effective width).
- The CL touching rule (3 µm) needs about 1–2 µm. The final maps reach that.
