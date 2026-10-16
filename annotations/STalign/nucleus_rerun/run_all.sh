#!/usr/bin/env bash
# Run the nucleus-guided STAlign registration (landmark affine + LDDMM) for all 11 STAlign samples.
# Each sample: old-map window check (cached) -> landmark affine -> LDDMM -> new meta + held-out
# validation (stalign_nucleus.py). PARALLEL samples run at once with THREADS torch threads each; keep
# PARALLEL x THREADS <= the allocation's slots (LDDMM is CPU-bound). Per-sample stdout goes to
# <OUT_DIR>/qc/stdout_<sample>.log. OUT_DIR defaults to the production folder (pathology_annotations_101326_nucleus);
# set it to an empty folder to rerun everything from the raw inputs and compare (provenance.py --compare).
# Usage: [OUT_DIR=<dir>] bash run_all.sh [sample ...]     (default: all 11; conda env spatial's python is used directly)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PY=/mnt/storage/dept/medonc/beroukhim/youyun/util/miniforge3/envs/spatial/bin/python
OUT_DIR=${OUT_DIR:-/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data/Xenium_annotations/pathology_annotations_101326_nucleus}
QC=$OUT_DIR/qc
PARALLEL=${PARALLEL:-2}
THREADS=${THREADS:-2}
SAMPLES=("$@")
[ ${#SAMPLES[@]} -eq 0 ] && SAMPLES=(CL42112 CL67068 CL87352 CN51420 CN53398 GG240468 GG241670 PA188570 PA258482 PA320784 PA328368)
mkdir -p "$QC"
export PY HERE OUT_DIR QC THREADS PYTHONWARNINGS=ignore
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$PARALLEL" -I{} bash -c \
    '"$PY" "$HERE/stalign_nucleus.py" {} --threads "$THREADS" --out-dir "$OUT_DIR" > "$QC/stdout_{}.log" 2>&1 && echo "done {}" || echo "FAILED {} (see $QC/stdout_{}.log)"'
