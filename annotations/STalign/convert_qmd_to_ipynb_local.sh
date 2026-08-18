#!/usr/bin/env bash
set -euo pipefail

# Convert STalign QMD notebooks to IPYNB for interactive/local use.
# Default behavior: convert all pathology_*_proseg.qmd files in this directory.
# You can also pass one or more qmd paths explicitly.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/ipynb_generated"
mkdir -p "${OUT_DIR}"

if ! command -v jupytext >/dev/null 2>&1; then
  echo "ERROR: jupytext not found in PATH." >&2
  echo "Activate your env first, e.g. conda activate spatial" >&2
  exit 1
fi

declare -a QMD_FILES

if [ "$#" -eq 0 ] || [ "${1:-}" = "--all" ]; then
  mapfile -t QMD_FILES < <(find "${SCRIPT_DIR}" -maxdepth 1 -type f -name 'pathology_*_proseg.qmd' | sort)
else
  for arg in "$@"; do
    if [ -f "${arg}" ]; then
      QMD_FILES+=("${arg}")
    elif [ -f "${SCRIPT_DIR}/${arg}" ]; then
      QMD_FILES+=("${SCRIPT_DIR}/${arg}")
    else
      echo "ERROR: QMD file not found: ${arg}" >&2
      exit 1
    fi
  done
fi

if [ "${#QMD_FILES[@]}" -eq 0 ]; then
  echo "ERROR: No QMD files found to convert." >&2
  exit 1
fi

echo "Converting ${#QMD_FILES[@]} QMD files to ${OUT_DIR}"
for qmd in "${QMD_FILES[@]}"; do
  base="$(basename "${qmd}" .qmd)"
  out="${OUT_DIR}/${base}.ipynb"
  echo "  - ${qmd} -> ${out}"
  jupytext --from quarto --to ipynb "${qmd}" -o "${out}"
done

echo "Done."
