#!/usr/bin/env bash
# Run the stress-adjusted DLGNT_1 vs DLGNT_2 pseudobulk DE/fgsea stage.
# Reuses the existing handoff/ TSVs produced by dlgnt12_prepare.py - this
# wrapper only re-runs the modeling stage with the IEG covariate enabled and
# optional --drop-samples / --no-stress-covariate / --keep-iegs-in-gsea flags.
# Pass any of those flags through positional args, e.g.:
#   bash run_dlgnt12_expression_stress_adjusted.sh --drop-samples=267134
#   bash run_dlgnt12_expression_stress_adjusted.sh --no-stress-covariate

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Run from repo root so relative paths inside any downstream tooling resolve
# the same way on local machines and the cluster.
cd "${ROOT_DIR}"

run_cmd() {
    # Echo the fully quoted command before execution so cluster logs show the
    # exact R invocation that was used for this run.
    printf '+'
    for arg in "$@"; do
        printf ' %q' "${arg}"
    done
    printf '\n'
    "$@"
}

R_CMD=(conda run -n spatial Rscript "${SCRIPT_DIR}/dlgnt12_limma_fgsea_stress_adjusted.R" "$@")

echo "Running stress-adjusted DLGNT_1 vs DLGNT_2 modeling from ${ROOT_DIR}"
run_cmd "${R_CMD[@]}"
