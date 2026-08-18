#!/usr/bin/env bash
# Run the full DLGNT_1 vs DLGNT_2 pseudobulk workflow in one command.
# This wrapper keeps the notebook read-only while avoiding manual step-by-step execution.
# The prepare stage resolves the input H5AD from the repo-standard runtime-dependent workdir.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Run from repo root so relative paths inside Quarto and any downstream helper
# code resolve the same way on local machines and the cluster.
cd "${ROOT_DIR}"

run_cmd() {
    # Echo the fully quoted command before execution so cluster logs show the
    # exact Python, R, and Quarto invocations that were used.
    printf '+'
    for arg in "$@"; do
        printf ' %q' "${arg}"
    done
    printf '\n'
    "$@"
}

# Wrapper arguments are forwarded to the Python prepare stage, which is the only
# stage in this workflow that exposes user-facing CLI options.
PYTHON_CMD=(conda run -n spatial python "${SCRIPT_DIR}/dlgnt12_prepare.py" "$@")
R_CMD=(conda run -n spatial Rscript "${SCRIPT_DIR}/dlgnt12_limma_fgsea.R")
QUARTO_CMD=(conda run -n spatial quarto render "${SCRIPT_DIR}/dlgnt12_expression.qmd")

echo "Running DLGNT_1 vs DLGNT_2 pseudobulk workflow from ${ROOT_DIR}"
run_cmd "${PYTHON_CMD[@]}"
run_cmd "${R_CMD[@]}"
run_cmd "${QUARTO_CMD[@]}"
