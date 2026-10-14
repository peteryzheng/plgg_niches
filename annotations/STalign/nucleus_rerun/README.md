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
- **Not done (user's call):** rerunning `annotations/STalign/add_annotations.R` on the new folder.
- **Write-up with figures:** Notion page "stalign rerun".
- **Code:** the pipeline as run is commit `7de415f` on this branch. Later commits add checks and
  documentation only. `NEW_DIR/provenance.json` records the HEAD it was written from and the md5 of
  every script.

### Reproduce

Use the conda env `spatial`. `run_all.sh` calls its python directly. `NEW_DIR` is the folder above.

```bash
# 1. registration, all 11 samples (~2.3 h on 4 slots). Per sample: old-map window check (cached in
#    qc/windows_before_<s>.tsv) -> nucleus landmarks -> affine -> 2-stage LDDMM -> new meta + validation
bash annotations/STalign/nucleus_rerun/run_all.sh                 # PARALLEL=2 THREADS=2 by default
#    affine-only reference variant (written to qc/trial_affine_all/)
python annotations/STalign/nucleus_rerun/stalign_nucleus.py <sample> --niter1 0 --niter2 0 --tag affine_all
# 2. label changes, figures, provenance manifest
python annotations/STalign/nucleus_rerun/label_impact.py   $NEW_DIR $NEW_DIR/qc $NEW_DIR/qc/trial_affine_all
python annotations/STalign/nucleus_rerun/cohort_figures.py $NEW_DIR $NEW_DIR/qc $NEW_DIR/figures $NEW_DIR/qc/trial_affine_all
# 3. diagnostics behind the explanations (outputs in qc/detection, qc/landmarks, qc/replay_old, qc/frame)
python annotations/STalign/nucleus_rerun/diag_detection.py         # one-to-one H&E/DAPI peak accounting (~15 min)
python annotations/STalign/nucleus_rerun/diag_detection.py --threshold-area   # DAPI threshold from the window only
python annotations/STalign/nucleus_rerun/diag_jpg_frame.py         # JPG vs SVS, polygon outlines vs tissue edges
python annotations/STalign/nucleus_rerun/stalign_nucleus.py PA320784 --threads 1 --tag repro_1thread   # determinism
python annotations/STalign/nucleus_rerun/diag_old_landmarks.py     # the notebooks' hand landmarks vs truth
for s in PA320784 CL42112 GG240468; do for l in original true true_grid; do   # old pipeline replays, ~6-11 min each
    python annotations/STalign/nucleus_rerun/replay_old_stalign.py $s $l --threads 2; done; done
python annotations/STalign/nucleus_rerun/replay_old_stalign.py PA320784 original --threads 1 --tag 1thread
python annotations/STalign/nucleus_rerun/replay_old_stalign.py --figure
python annotations/STalign/nucleus_rerun/provenance.py            # after committing; records HEAD + md5s
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
| `replay_old_stalign.py` | Reruns the original notebook pipeline with hand, exact or grid-rounded landmarks |
| `diag_jpg_frame.py` | Checks that the JPG is the SVS at level 0 / 7, and the polygon outlines against tissue edges |
| `provenance.py` | Provenance manifest |

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
- **The LDDMM is kept from drifting**, which the old one did even from exact landmarks:
  - the affine is held fixed (epL = epT = 0);
  - the step is small (epV 0.3, against 10);
  - the landmark term is 40× tighter (σP 0.5 µm, against 20 µm), over 20–49 landmarks instead of 4–6;
  - deformation is local (a = 150–200 µm, against about 880 µm);
  - it runs 250 iterations, against 2000.
- **Its image term still pulls a little.** In 10 of 11 samples, held-out medians are 0.90–1.41 µm
  with LDDMM against 0.56–0.90 µm for the affine alone. This is the same tendency as the old pipeline,
  held under 1 µm by the landmarks. GG241670 is the exception: LDDMM takes it from 2.25 to 1.00 µm.
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

### Why the old map was off (verified 2026-10-14; corrects the handoff text below)

**The hand landmarks were far off, and not because of the 50 grid** (`diag_old_landmarks.py`,
`qc/landmarks/`).
- There are 48 landmarks in total. Placed on the new map, they are off by a median 189 µm (29–630 µm).
  Per-sample medians are 80–373 µm.
- In 9 of 11 notebooks every coordinate is a multiple of 50: JPG px on the H&E side (88 µm) and µm on
  the Xenium side.
- Rounding a perfectly placed pair to that grid moves it at most 97 µm, typically 37 µm (simulated).
  Yet 77% of landmarks are off by more than 97 µm.
- The grid point nearest each true position was only 4–56 µm away.
- The two notebooks not on the grid (CL87352, GG240468) were as far off (187 and 178 µm).
- So the landmarks were misidentified features, not rounded ones.

**The landmarks didn't decide the final map either. The old LDDMM did, and it doesn't converge to the
truth** (`replay_old_stalign.py`, `qc/replay_old/`). The original notebook pipeline was replayed:
whole-slide JPG vs 10 µm cell-density raster, a = 500 JPG px, 2000 iterations, sigmaM 0.15, epV 10,
landmark term on.

| Sample | Landmarks | Start (affine from landmarks) | End (after old LDDMM) | 030226 |
|---|---|---|---|---|
| PA320784 | hand | 243 µm | 10.4 µm (2 threads), 22.8 µm (1 thread) | 25.3 µm |
| PA320784 | exact | 0.2 µm | 32.1 µm | |
| PA320784 | exact, rounded to grid | 43.5 µm | 19.0 µm | |
| CL42112 | hand | 70.3 µm | 65.2 µm | 76.8 µm |
| CL42112 | exact | 0.8 µm | 78.1 µm | |
| CL42112 | exact, rounded to grid | 47.0 µm | 78.7 µm | |
| GG240468 | hand | 115 µm | 5.4 µm | 4.4 µm |
| GG240468 | exact | 0.2 µm | 3.5 µm | |

- **Exact landmarks don't help.** They give a near-perfect start (0.2–0.8 µm), and the old LDDMM
  then pulls the map 32–78 µm away in PA320784 and CL42112. Its image term matches H&E colour to cell
  counts per 10 µm pixel, and that match is not best at the true alignment.
- **Each section has its own old-pipeline error, whatever the start:** about 4 µm in GG240468, 10–32
  µm in PA320784, 65–79 µm in CL42112. That is why GG240468 came out right. It was the section, not
  its landmarks.
- **The old pipeline isn't reproducible at this scale.** The same PA320784 run on 1 vs 2 CPU threads
  ends 26 µm apart, from identical inputs. Replays with the notebooks' own landmarks land 6–89 µm from
  030226, so the stored map is one draw from a wide spread.
- **JPG resolution (1.76 µm/px) was not the cause.**
- **Nor was a JPG/SVS mismatch.** The JPG is the SVS at level 0 / 7 to within 0.73 µm at 66 patches
  across all 11 sections (`diag_jpg_frame.py`, `qc/frame/jpg_vs_svs.tsv`), so "exact" landmarks in
  the replays are exact in the JPG the old pipeline registered.
- **Not established:** why the image term's optimum is off. One hypothesis is that H&E darkness tracks
  cell density inconsistently across tissue types, plus edge and background handling, but it hasn't
  been tested.

### Why the 5 µm target

- Neighbouring DAPI nuclei sit a median 8.7–10.3 µm apart. Beyond about half that spacing, a mapped
  point is closer to the next nucleus than its own, so per-cell H&E features attach to the wrong cell.
- Region labels tolerate more, as long as errors stay well under region width. The annotated regions
  here are often only 40–140 µm wide (CN51420 pseudorosettes 43 µm, CL42112 lipidized islands 135 µm,
  effective width).
- The CL touching rule (3 µm) needs about 1–2 µm. The final maps reach that.

---

## Original handoff brief (2026-10-13)

This brief was written for the Claude session that ran this work. That work started in a different
session (the CL lipid-map work, run from the main checkout `/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code`).
A session in this directory has a different project-memory folder, so everything needed is written
here.

## Problem (found 2026-10-13)

The stored H&E↔Xenium registration is in `data/Xenium_annotations/pathology_annotations_030226/proseg_pathology_<sample>_proseg_meta.csv`.
Its columns `STAlign_x` / `STAlign_y` are the per-cell H&E position in JPG pixels: row and column,
JPG = SVS level 0 ÷ 7. This registration is off by tens of µm in 10 of 11 samples.

**How it was measured.** The H&E is a post-Xenium stain of the same slide, so every H&E nucleus should
land on a Xenium DAPI nucleus.
- Each sample gets 30 windows of 150 µm.
- H&E nuclei: peaks of the hematoxylin channel. DAPI nuclei: peaks of `morphology_focus`.
- The DAPI points are mapped into the H&E through the stored coordinates.
- The best shift per window comes from FFT cross-correlation (±150 µm), with a peak z-score.

| Sample | Confident windows | Median offset (µm) | p90 (µm) | Spread (µm) | Annotated cells whose label changes |
|---|---|---|---|---|---|
| CL42112 | 26/29 | 75 | 94 | 21 | 41% |
| CL67068 | 21/28 | 25 | 50 | 33 | 7% |
| CL87352 | 24/30 | 26 | 44 | 31 | 13% |
| CN51420 | 20/30 | 54 | 58 | 11 | 70% |
| CN53398 | 13/30 | 38 | 44 | 17 | 32% |
| GG240468 | 30/30 | **3** | 7 | 4 | 0% (positive control) |
| GG241670 | 13/30 | 62 | 98 | 75 | 14% |
| PA188570 | 29/30 | 34 | 42 | 10 | 18% |
| PA258482 | 30/30 | 16 | 22 | 13 | 1% |
| PA320784 | 30/30 | 26 | 32 | 6 | 9% |
| PA328368 | 29/30 | 35 | 45 | 11 | 6% |

- **Reading the table.** A small spread means a near-rigid offset (CL42112, CN51420, the PA samples).
  A large spread means a local, non-rigid error (CL67068, CL87352, GG241670).
- **"Label changes"** = the share of pathology-annotated cells whose region label changes when the
  map moves by the median offset.
- **Source tables:** `<main>/sc_integration/ENVI/2_spatial_imputation/CL/segmentation_qc/registration/`,
  in `offsets_summary.tsv`, `annotation_impact.tsv` and `offsets_<sample>.tsv` (per window).
- **Picture:** `R0_nuclei_before_after.jpg` in the same folder.

**Why it happened.** The original notebooks (`annotations/STalign/pathology_<sample>_proseg.qmd`):
- align the whole-slide **JPG at 1.76 µm/px**. LGG1 and LGG2 each hold several sections on one slide.
- align it to a **10 µm cell-density raster** of ProSeg centroids (`STalign.rasterize(..., dx=10)`).
- start from **4–6 hand-placed landmarks** via `L_T_from_points`.
- run LDDMM with **default smoothness `a = 500`**, plus niter 2000, sigmaM 0.15, sigmaB 0.10,
  sigmaA 0.11 and epV 10.

That setup is fine at mm scale. It can't resolve anything under about 100 µm, and the cell-density
raster has no nucleus-level features to lock onto. *(2026-10-14: corrected in "Why the old map was off" above. Replays show the old LDDMM ends 3.5–78 µm
from the truth even from exact landmarks.)*

**Why it matters.**
- Cell-level pathology labels near region edges are affected.
- So is everything in the CL lipid work below about 100 µm: the touching rule (nucleus within 3 µm
  of a lipid pocket) and the rim profiles.
- Treat the current labels as reliable at mm scale only.

## Decision (user)

Rerun STAlign itself rather than build a custom correction field. User: "I feel like this is what
STAlign was meant to do. so we shouldn't reinvent the wheel here."

## Pilot design (CL42112 first)

1. **Images.** Both sides show nuclei.
   - H&E: the hematoxylin channel (`skimage.color.rgb2hed(rgb)[..., 0]`) from the **SVS**, not the
     JPG, at about 4 µm/px, then 2 µm/px. Crop to this sample's section.
   - Xenium: the DAPI `morphology_focus` image at the same resolution.
2. **Initial affine from nucleus-matched landmarks.**
   - Use windows where nuclei match confidently (FFT shift, peak z-score ≥ about 5).
   - Each gives a pair: H&E window centre ↔ that point mapped into Xenium µm with its measured
     offset applied.
   - Use **half** the windows for the fit; hold the other half out for scoring.
   - Feed the pairs to `STalign.L_T_from_points`.
3. **LDDMM, coarse to fine.**
   - Smoothness `a` about 150–200 µm, in the image's units, instead of 500.
   - Keep the landmark term (`pointsI` / `pointsJ`).
   - First at 4 µm/px, then refine at 2 µm/px.
4. **Output.**
   - Write a new `_meta.csv` per sample, with the same columns as the 030226 files (`STAlign_x`,
     `STAlign_y` in JPG px and the region labels), to a new folder such as
     `data/Xenium_annotations/pathology_annotations_<date>_nucleus/`.
   - **Never modify `pathology_annotations_030226`.** Record its md5s before you start and check them
     at the end.
5. **Validation.**
   - Rerun the nucleus check against the new coordinates, on the **held-out** windows only.
   - Target: median ≤ 5 µm. Report p90, the match fraction, and before vs after.
   - Make an overlay picture: H&E with DAPI nuclei on top, before and after.
6. **⛔ Gate.** Show the user the CL42112 result before running any other sample.
   - After approval, run all 11 samples.
   - GG240468 is the no-change control: it is already at 3 µm and must stay there.
   - Then rerun `annotations/STalign/add_annotations.R` against the new folder. The user decides when.

## Code to reuse

**The nucleus check** is `<main>/sc_integration/ENVI/2_spatial_imputation/CL/segmentation_qc/registration_nuclei_check.py`.
- It is **untracked** in the main checkout and absent from this branch. Import it by absolute path:
  ```python
  sys.path.insert(0, "/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/sc_integration/ENVI/2_spatial_imputation/CL/segmentation_qc")
  from registration_nuclei_check import DAPI, peaks, raster, best_shift, folder, GROUP, DATA, DAPI_UM, MATCH_UM
  from he_image import SVS, MPP0, JPG_TO_L0          # SVS(path).read_scaled(r0, c0, h, w, factor)
  ```
- Don't edit files in the main checkout; another session is working there.
- `run()` / `run_cells()` hard-code the 030226 folder and, for CL samples, the lipid candidate masks.
  Write your own window loop (copy `run_cells`, which needs no lipid data) with the metadata path as
  a parameter.
- Convention: `best_shift` returns d such that DAPI_mapped ≈ H&E + d, in µm, as (row, col).

**STalign** (conda env `spatial`):
- `LDDMM(xI, I, xJ, J, pointsI, pointsJ, L, T, a=500.0, niter=5000, ...)`
- `L_T_from_points`
- `transform_points_target_to_source`
- `build_transform`, `rasterize`, `normalize`

## Data paths (base `DATA = /mnt/storage/dept/medonc/beroukhim/youyun/plgg/data`)

| What | Path |
|---|---|
| SVS (level 0 = 0.2522 µm/px) | `DATA/Xenium_annotations/images/230918_Xenium_CytAssist_{LGG1,LGG2}.svs`; the slide group per sample is `GROUP` in the check |
| JPG (level 0 ÷ 7) | the same name with `.jpg` |
| Pathology polygons | `DATA/Xenium_annotations/geojsons/230918_Xenium_CytAssist_{group}.geojson`, in level-0 px; the notebooks scale them by 1/7 to JPG |
| Xenium DAPI (0.2125 µm/px, tiled JPEG2000) | `DATA/segmentation/proseg_run_121024/<folder>/<folder>_proseg/outs/morphology_focus/morphology_focus_0000.ome.tif`; `<folder>` comes from `folder(sample)` via `DATA/metadata/Xenium_PS.xlsx` |
| ProSeg cells | `.../outs/cells.csv.gz` in the same folder |
| Old registration (read-only) | `DATA/Xenium_annotations/pathology_annotations_030226/proseg_pathology_<sample>_proseg_meta.csv` |

**Frames.**
- Xenium `x_centroid` (µm) runs along H&E rows (`STAlign_x`).
- `y_centroid` runs against H&E columns (`STAlign_y`): correlations +1.0 and −1.0 respectively.
  Check this per sample.

## Compute and standing rules

- **Cluster.** This node (argos4) **cannot submit jobs.**
  - Write an SGE script (`#$ -pe smp 4`; 4-slot jobs start sooner on the fragmented queue) for the
    user to `qsub` from **argos-qsub1**.
  - The pilot at 4 µm/px may fit in an interactive run.
- **Python.** `python` is Python 2; activate `conda activate spatial` first.
- **Rules:**
  - Don't run `quarto render`; the user renders.
  - Commit only when asked. When asked, commit just this folder's new files on this branch, e.g.
    `annotations: nucleus-guided STAlign rerun`.
  - Follow AGENTS.md style: snake_case, a header comment on every script, and comments on every
    non-obvious step.
  - Don't commit generated outputs.
- **Patient H&E** images go to the user only via SendUserFile or Notion. Never publish them as an
  Artifact or upload them to consumer cloud services.
- **Explanations:** plain language, and reason through confounders, e.g. whether a window "improved"
  only because it was a landmark window.
