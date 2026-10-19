#!/usr/bin/env python3
# Provenance manifest for the nucleus-guided STAlign rerun (pathology_annotations_101326_nucleus).
#
# Writes, next to the new meta files:
#   provenance.json      : code commit (+ whether the pipeline folder was clean), md5 of every pipeline
#                          script, software versions, fixed constants, and each sample's run arguments
#                          (read back from its transform_<sample>.npz, i.e. what actually ran)
#   provenance_files.tsv : role, sample, path, bytes, md5 for every input and output file
#   environment_spatial.yml : `conda env export` of the env that ran it
# It also records, when present, the downstream annotated BANKSY object (annotations/STalign/add_annotations.R)
# and its check (check_annotated_object.R), the per-cell metadata TSV exported from it, the PA pathology-enrichment
# reruns on the new labels and their comparison with the old-label outputs (compare_pa_enrichment.R).
# Usage: python provenance.py          (conda env spatial; after run_all.sh, from a committed checkout)
#        python provenance.py --compare <dir>   compare a rerun (run_all.sh with OUT_DIR=<dir>) with the
#                                               production outputs; writes <dir>/compare_to_production.tsv
# Assumes: inputs are those stalign_nucleus.py reads (old 030226 meta, SVS, GeoJSON + recode map from the
#          sample notebook, Xenium DAPI, Xenium_PS.xlsx for the folder lookup); outputs in NEW_DIR(/qc).
import hashlib, importlib, json, platform, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nucleus_core, nucleus_windows
from stalign_nucleus import NEW_DIR, NOTEBOOKS, N_WIN, WIN_SPACING_UM, SVS_PAGE, DAPI_LEVEL, notebook_spec
from nucleus_windows import DATA, GROUP, Z_CONF, MIN_CELLS, meta_path, dapi_path

SCRIPTS = ["nucleus_core.py", "nucleus_windows.py", "stalign_nucleus.py", "run_all.sh", "label_impact.py", "cohort_figures.py",
           "overlay_before_after.py", "diag_detection.py", "diag_old_landmarks.py", "replay_old_stalign.py", "diag_jpg_frame.py",
           "diag_image_term.py", "check_annotated_object.R", "compare_pa_enrichment.R", "provenance.py"]
# downstream scripts outside this folder that consume the new labels (repo-relative)
DOWNSTREAM = ["annotations/STalign/add_annotations.R",
              "sc_integration/ENVI/2_spatial_imputation/PA/export_new_object_cell_metadata.R",
              "sc_integration/ENVI/2_spatial_imputation/PA/morphology_state_gsea/morphology_state_prepare.py",
              "sc_integration/ENVI/2_spatial_imputation/PA/morphology_state_gsea/morphology_state.R",
              "sc_integration/ENVI/2_spatial_imputation/measured_panel_qc/measured_panel_prepare.py",
              "sc_integration/ENVI/2_spatial_imputation/measured_panel_qc/measured_panel_fit.R",
              "sc_integration/ENVI/2_spatial_imputation/measured_panel_qc/compare_measured_imputed.py",
              "annotations/STalign/plot_endothelial_vessel_registration.py"]
REPO = HERE.parents[2]
# distribution name -> module name. Versions are read from the imported module: the env holds a stale pip
# torch-2.0.0.dist-info next to the conda torch 2.10.0 that is actually imported, so package metadata misleads.
PACKAGES = {"torch": "torch", "numpy": "numpy", "scipy": "scipy", "scikit-image": "skimage", "tifffile": "tifffile",
            "imagecodecs": "imagecodecs", "pandas": "pandas", "geopandas": "geopandas", "shapely": "shapely", "Pillow": "PIL"}
# STalign has no version of its own ("1.0"); it was installed from the GitHub main zip. Every installed package file
# is identical (md5) to upstream main b2068ed (2023-12-11, still main on 2026-10-15); STalign.py last changed at 5837b03.
STALIGN_REPO, STALIGN_COMMIT = "https://github.com/JEFworks-Lab/STalign", "b2068edc98974efa54537eca194736e177bbe11d"


def md5(p, block=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        while chunk := f.read(block):
            h.update(chunk)
    return h.hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True).stdout.strip()


def software():
    """Versions of what is actually imported, the CPU path torch uses, and STalign's source identity."""
    import torch, STalign.STalign as st
    out = dict(python=sys.version.split()[0], platform=platform.platform(),
               **{d: importlib.import_module(m).__version__ for d, m in PACKAGES.items()})
    out.update(torch_cpu_capability=torch.backends.cpu.get_cpu_capability(), torch_threads_default=torch.get_num_threads(),
               STalign=dict(installed="1.0, from the GitHub main zip", repo=STALIGN_REPO, commit=STALIGN_COMMIT,
                            STalign_py_md5=md5(st.__file__)))
    return out


def export_env():
    """`conda env export` of the running env (envs/<name>/bin/python -> <base>/bin/conda) to NEW_DIR."""
    prefix = Path(sys.prefix)
    r = subprocess.run([str(prefix.parents[1] / "bin/conda"), "env", "export", "-p", str(prefix)], capture_output=True, text=True)
    (NEW_DIR / "environment_spatial.yml").write_text(r.stdout if r.returncode == 0 else f"# conda env export failed:\n# {r.stderr}")


def compare(other):
    """A rerun in `other` (same layout as NEW_DIR) against the production outputs, per sample: identical meta file?
    largest cell displacement, label differences, and the largest differences in the window tables and the transform."""
    rows = []
    for s in GROUP:
        a, b = other / f"proseg_pathology_{s}_proseg_meta.csv", NEW_DIR / f"proseg_pathology_{s}_proseg_meta.csv"
        if not a.exists():
            continue
        r = dict(sample=s, meta_identical=md5(a) == md5(b))
        ma, mb = pd.read_csv(a), pd.read_csv(b)
        assert (ma.cell_id.values == mb.cell_id.values).all(), f"{s}: cell order differs"
        d = np.hypot(ma.STAlign_x - mb.STAlign_x, ma.STAlign_y - mb.STAlign_y) * nucleus_core.JPG_TO_L0 * nucleus_core.MPP0
        r.update(cell_shift_max_um=d.max(), cell_shift_median_um=d.median(),
                 annotation_differs=int((ma.annotation.fillna("") != mb.annotation.fillna("")).sum()),
                 region_id_differs=int((ma.region_id.astype(str) != mb.region_id.astype(str)).sum()))
        for tab in ("windows_before", "windows_after"):
            ta, tb = (pd.read_csv(x / "qc" / f"{tab}_{s}.tsv", sep="\t") for x in (other, NEW_DIR))
            c = ["row_l0", "col_l0", "shift_row_um", "shift_col_um", "peak_z", "match_at_zero"]
            r[f"{tab}_max_abs_diff"] = np.nanmax(np.abs(ta[c].values - tb[c].values)) if ta.shape == tb.shape else np.nan
        ra = pd.read_csv(other / "qc" / f"windows_before_roles_{s}.tsv", sep="\t").role
        r["roles_identical"] = (ra.values == pd.read_csv(NEW_DIR / "qc" / f"windows_before_roles_{s}.tsv", sep="\t").role.values).all()
        za, zb = (np.load(x / "qc" / f"transform_{s}.npz") for x in (other, NEW_DIR))
        r.update(transform_A_max_abs_diff=np.abs(za["A"] - zb["A"]).max(), transform_v_max_abs_diff_um=np.abs(za["v"] - zb["v"]).max())
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(other / "compare_to_production.tsv", sep="\t", index=False)
    print(t.to_string(index=False))


def main():
    if sys.argv[1:2] == ["--compare"]:
        return compare(Path(sys.argv[2]))
    qc = NEW_DIR / "qc"
    export_env()
    # every file the pipeline read or wrote, per sample (shared inputs listed once per sample for clarity)
    files = []
    for s in GROUP:
        gj, _ = notebook_spec(s)
        files += [("input_old_meta", s, meta_path(s)),
                  ("input_he_svs", s, DATA / f"Xenium_annotations/images/230918_Xenium_CytAssist_{GROUP[s]}.svs"),
                  ("input_geojson", s, gj),
                  ("input_notebook_recode_map", s, NOTEBOOKS / f"pathology_{s}_proseg.qmd"),
                  ("input_xenium_dapi", s, dapi_path(s)),
                  ("input_xenium_folder_lookup", s, DATA / "metadata/Xenium_PS.xlsx"),
                  ("output_meta", s, NEW_DIR / f"proseg_pathology_{s}_proseg_meta.csv"),
                  ("output_transform", s, qc / f"transform_{s}.npz"),
                  ("qc_windows_old_map", s, qc / f"windows_before_{s}.tsv"),
                  ("qc_windows_roles", s, qc / f"windows_before_roles_{s}.tsv"),
                  ("qc_windows_new_map", s, qc / f"windows_after_{s}.tsv"),
                  ("qc_summary", s, qc / f"summary_{s}.tsv"),
                  ("qc_log", s, qc / f"run_{s}.log")]
    # cohort-level QC tables and figures (label_impact.py, cohort_figures.py)
    files += [("qc_label_impact", "all", qc / f) for f in ["label_impact.tsv", "label_retention.tsv"]]
    files += [("figure", "all", f) for f in sorted((NEW_DIR / "figures").glob("*.jpg")) + sorted((NEW_DIR / "figures").glob("*.png"))]
    # diagnostics of the old map and of the nucleus check (diag_detection.py, diag_old_landmarks.py, replay_old_stalign.py,
    # diag_jpg_frame.py) and the 1-thread rerun of PA320784 (determinism across thread counts)
    for sub in ("detection", "landmarks", "replay_old", "frame", "trial_repro_1thread", "image_term"):
        files += [(f"diag_{sub}", "all", f) for f in sorted((qc / sub).glob("*")) if f.suffix in (".tsv", ".gz", ".png", ".jpg", ".npz", ".csv")]
    # the from-scratch rerun's comparison with these outputs, and the environment snapshot
    files += [("repro_from_scratch", "all", qc / "repro_from_scratch/compare_to_production.tsv"),
              ("environment", "all", NEW_DIR / "environment_spatial.yml")]
    # downstream: the BANKSY object labelled from these metas by add_annotations.R, its 5% subset, and the check
    obj = DATA / "Xenium_annotations/banksy_clusters_connected_annotated_k_geom_15_30_pc_20_lam_0.2_0.8_kct_50_resct_0.2_kni_50_resni_0.5_20260826_115134_pathology_nucleus"
    if Path(str(obj) + ".rds").exists():                        # (the name has dots, so no with_suffix)
        files += [("annotated_object", "all", Path(str(obj) + ".rds")),
                  ("annotated_object_subset", "all", Path(str(obj) + "_subset.rds")),
                  ("annotated_object_check", "all", qc / "annotated_object_check.tsv"),
                  ("annotated_object_log", "all", qc / "add_annotations_20261016.log")]
    # downstream: the per-cell metadata TSV exported from that object, the old-label snapshot, the PA reruns, the comparison
    tsv = DATA / "xenium_cell_types/cell_metadata_resct0.2_resni0.5_20260826_pathology_nucleus.tsv"
    if tsv.exists():
        files += [("cell_metadata_tsv", "all", tsv), ("cell_metadata_tsv_log", "all", qc / "export_cell_metadata_20261019.log")]
    files += [(role, "all", qc / name) for role, name in [
        ("pa_old_label_snapshot", "pa_enrichment_old_labels/MANIFEST.md5"),
        ("pa_morphology_rerun_log", "morphology_state_rerun_20261019.log"),
        ("pa_measured_panel_rerun_log", "measured_panel_rerun_20261019.log"),
        ("pa_enrichment_comparison", "pa_enrichment_old_vs_new.tsv"),
        ("pa_enrichment_comparison_summary", "pa_enrichment_old_vs_new_summary.tsv")] if (qc / name).exists()]
    paths = sorted({p for _, _, p in files})
    missing = [str(p) for p in paths if not Path(p).exists()]
    with ThreadPoolExecutor(4) as ex:                       # ~9 GB of images: hash in parallel
        sums = dict(zip(paths, ex.map(lambda p: md5(p) if Path(p).exists() else "MISSING", paths)))
    tab = pd.DataFrame([dict(role=r, sample=s, path=str(p), bytes=Path(p).stat().st_size if Path(p).exists() else np.nan,
                             md5=sums[p]) for r, s, p in files])
    tab.to_csv(NEW_DIR / "provenance_files.tsv", sep="\t", index=False)

    # the arguments each sample actually ran with, as stored in its transform file
    run_args = {}
    for s in GROUP:
        t = qc / f"transform_{s}.npz"
        if t.exists():
            run_args[s] = json.loads(str(np.load(t)["args"]))
    man = dict(
        created=time.strftime("%F %T"),
        code=dict(repo_head=git("rev-parse", "HEAD"), branch=git("rev-parse", "--abbrev-ref", "HEAD"),
                  pipeline_folder_clean=git("status", "--porcelain", "--", ".") == "",
                  folder="annotations/STalign/nucleus_rerun",
                  script_md5={f: md5(HERE / f) for f in SCRIPTS if (HERE / f).exists()},
                  downstream_script_md5={f: md5(REPO / f) for f in DOWNSTREAM if (REPO / f).exists()}),
        commands=["bash annotations/STalign/nucleus_rerun/run_all.sh   (PARALLEL=2 THREADS=2 by default)",
                  "python annotations/STalign/nucleus_rerun/label_impact.py <NEW_DIR> <NEW_DIR>/qc <NEW_DIR>/qc/trial_affine_all",
                  "python annotations/STalign/nucleus_rerun/cohort_figures.py <NEW_DIR> <NEW_DIR>/qc <NEW_DIR>/figures <NEW_DIR>/qc/trial_affine_all",
                  "python annotations/STalign/nucleus_rerun/diag_detection.py",
                  "python annotations/STalign/nucleus_rerun/diag_detection.py --threshold-area",
                  "python annotations/STalign/nucleus_rerun/diag_jpg_frame.py",
                  "python annotations/STalign/nucleus_rerun/stalign_nucleus.py PA320784 --threads 1 --tag repro_1thread",
                  "python annotations/STalign/nucleus_rerun/diag_old_landmarks.py",
                  "python annotations/STalign/nucleus_rerun/replay_old_stalign.py <sample> {original,true,true_grid} [--threads 1 --tag 1thread]",
                  "python annotations/STalign/nucleus_rerun/replay_old_stalign.py --figure",
                  "python annotations/STalign/nucleus_rerun/diag_image_term.py <sample> [...]; python annotations/STalign/nucleus_rerun/diag_image_term.py --figure",
                  "OUT_DIR=<NEW_DIR>/qc/repro_from_scratch bash annotations/STalign/nucleus_rerun/run_all.sh",
                  "python annotations/STalign/nucleus_rerun/provenance.py --compare <NEW_DIR>/qc/repro_from_scratch",
                  "Rscript annotations/STalign/add_annotations.R   (annotated BANKSY object from these metas)",
                  "Rscript annotations/STalign/nucleus_rerun/check_annotated_object.R",
                  "Rscript sc_integration/ENVI/2_spatial_imputation/PA/export_new_object_cell_metadata.R   (per-cell TSV with the new labels)",
                  "(in PA/morphology_state_gsea/) python morphology_state_prepare.py; Rscript morphology_state.R",
                  "(in measured_panel_qc/) python measured_panel_prepare.py PA; Rscript measured_panel_fit.R PA; python compare_measured_imputed.py",
                  "python annotations/STalign/plot_endothelial_vessel_registration.py",
                  "quarto render annotations/linear_models.qmd, annotations/linear_models_niches.qmd, niche/proseg_output_analysis/niche_interpretation_proseg.qmd, sc_integration/ENVI/2_spatial_imputation/PA/4_PA_morphology_state.qmd   (by the user)",
                  "Rscript annotations/STalign/nucleus_rerun/compare_pa_enrichment.R",
                  "python annotations/STalign/nucleus_rerun/provenance.py"],
        new_dir=str(NEW_DIR),
        software=software(),
        constants=dict(N_WIN=N_WIN, WIN_SPACING_UM=WIN_SPACING_UM, MIN_CELLS=MIN_CELLS, Z_CONF=Z_CONF,
                       WIN_UM=nucleus_core.WIN_UM, SHIFT_UM=nucleus_core.SHIFT_UM, MATCH_UM=nucleus_core.MATCH_UM,
                       NUC_MIN_DIST_UM=nucleus_core.NUC_MIN_DIST_UM, MPP0=nucleus_core.MPP0, JPG_TO_L0=nucleus_core.JPG_TO_L0,
                       DAPI_UM=nucleus_core.DAPI_UM, SVS_PAGE=SVS_PAGE, DAPI_LEVEL=DAPI_LEVEL,
                       window_seed=0, landmark_split_seed=1),
        run_args=run_args,
        missing_files=missing,
    )
    (NEW_DIR / "provenance.json").write_text(json.dumps(man, indent=2))
    print(f"wrote {NEW_DIR / 'provenance.json'} and provenance_files.tsv ({len(tab)} rows, {len(missing)} missing)")


if __name__ == "__main__":
    main()
