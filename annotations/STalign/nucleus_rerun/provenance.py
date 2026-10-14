#!/usr/bin/env python3
# Provenance manifest for the nucleus-guided STAlign rerun (pathology_annotations_101326_nucleus).
#
# Writes, next to the new meta files:
#   provenance.json      : code commit (+ whether the pipeline folder was clean), md5 of every pipeline
#                          script, software versions, fixed constants, and each sample's run arguments
#                          (read back from its transform_<sample>.npz, i.e. what actually ran)
#   provenance_files.tsv : role, sample, path, bytes, md5 for every input and output file
# Usage: python provenance.py          (conda env spatial; after run_all.sh, from a committed checkout)
# Assumes: inputs are those stalign_nucleus.py reads (old 030226 meta, SVS, GeoJSON + recode map from the
#          sample notebook, Xenium DAPI, Xenium_PS.xlsx for the folder lookup); outputs in NEW_DIR(/qc).
import hashlib, json, platform, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nucleus_core, nucleus_windows
from stalign_nucleus import NEW_DIR, NOTEBOOKS, N_WIN, WIN_SPACING_UM, SVS_PAGE, DAPI_LEVEL, notebook_spec
from nucleus_windows import DATA, GROUP, Z_CONF, MIN_CELLS, meta_path, dapi_path

SCRIPTS = ["nucleus_core.py", "nucleus_windows.py", "stalign_nucleus.py", "run_all.sh", "label_impact.py", "cohort_figures.py",
           "overlay_before_after.py", "diag_detection.py", "diag_old_landmarks.py", "replay_old_stalign.py", "provenance.py"]
PACKAGES = ["STalign", "torch", "numpy", "scipy", "scikit-image", "tifffile", "imagecodecs", "pandas", "geopandas", "shapely"]


def md5(p, block=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        while chunk := f.read(block):
            h.update(chunk)
    return h.hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True).stdout.strip()


def main():
    qc = NEW_DIR / "qc"
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
    # diagnostics of the old map and of the nucleus check (diag_detection.py, diag_old_landmarks.py, replay_old_stalign.py)
    for sub in ("detection", "landmarks", "replay_old"):
        files += [(f"diag_{sub}", "all", f) for f in sorted((qc / sub).glob("*")) if f.suffix in (".tsv", ".gz", ".png", ".jpg", ".npz")]
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
                  script_md5={f: md5(HERE / f) for f in SCRIPTS if (HERE / f).exists()}),
        commands=["bash annotations/STalign/nucleus_rerun/run_all.sh   (PARALLEL=2 THREADS=2 by default)",
                  "python annotations/STalign/nucleus_rerun/label_impact.py <NEW_DIR> <NEW_DIR>/qc <NEW_DIR>/qc/trial_affine_all",
                  "python annotations/STalign/nucleus_rerun/cohort_figures.py <NEW_DIR> <NEW_DIR>/qc <NEW_DIR>/figures <NEW_DIR>/qc/trial_affine_all",
                  "python annotations/STalign/nucleus_rerun/diag_detection.py",
                  "python annotations/STalign/nucleus_rerun/diag_old_landmarks.py",
                  "python annotations/STalign/nucleus_rerun/replay_old_stalign.py <sample> {original,true,true_grid} [--threads 1 --tag 1thread]",
                  "python annotations/STalign/nucleus_rerun/replay_old_stalign.py --figure",
                  "python annotations/STalign/nucleus_rerun/provenance.py"],
        new_dir=str(NEW_DIR),
        software=dict(python=sys.version.split()[0], platform=platform.platform(),
                      **{p: version(p) for p in PACKAGES}),
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
