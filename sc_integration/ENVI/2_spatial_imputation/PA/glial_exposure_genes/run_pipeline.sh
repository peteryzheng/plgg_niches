#!/usr/bin/env bash
# Run the gene-level myeloid glial-adjacency pipeline (Python streaming -> R aggregation).
#
# Like the CL morphology-state pipeline this does NOT need a big-RAM node: the Python stage
# streams obsm['imputation'] one gene column at a time, so peak memory is ~one column plus
# the dense 266-gene panel (a few hundred MB). It is I/O bound -- ~26 GB of sequential
# reads across the 4 PA ENVI objects -- so run it where the data is LOCAL (argos), never
# over the sshfs mount.
#
# Run from the spatial2 tmux on argos. Logs are teed to the shared FS next to this script
# (not /tmp). Verify by artifacts, not exit code.
#
# Usage:
#   bash run_pipeline.sh                 # ENVI_results_091126, 18,077 genes
#
# The old ~2.2k per-sample-HVG sensitivity arm is retired: its ENVI_results/ objects were
# deleted in the 2026-09-11 export migration, and the current export carries one uniform
# 18,077-gene universe across every sample.
set -o pipefail
cd "$(dirname "$0")" || exit 1

ENVI_DIR="ENVI_results_091126"
OUT_TAG="glial_exposure_genes"
EXTRA=""
export OUT_TAG

echo "=== PREPARE start $(date) host=$(hostname) envi=${ENVI_DIR} tag=${OUT_TAG} ==="
python glial_exposure_genes_prepare.py --envi-dir "${ENVI_DIR}" --out-tag "${OUT_TAG}" \
    ${EXTRA} 2>&1 | tee "prepare_${OUT_TAG}.log"
pex=${PIPESTATUS[0]}
echo "PREPARE_EXIT=${pex}"
if [ "${pex}" -ne 0 ]; then
    echo "prepare failed; aborting before R stage."
    exit "${pex}"
fi

# TEST A: closed-form vs a real least-squares fit (arithmetic check).
echo "=== VERIFY TEST A start $(date) ==="
python verify_glial_exposure_genes.py test-a \
    --per-sample "${OUT_TAG}_per_sample.tsv" --envi-dir "${ENVI_DIR}" \
    2>&1 | tee "verify_a_${OUT_TAG}.log"

echo "=== RSTAGE start $(date) ==="
Rscript glial_exposure_genes.R 2>&1 | tee "rstage_${OUT_TAG}.log"
rex=${PIPESTATUS[0]}
echo "RSTAGE_EXIT=${rex}"
echo "=== DONE $(date) ==="
