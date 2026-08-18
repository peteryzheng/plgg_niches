#!/bin/bash -l
#$ -N stalign_ipynb
#$ -cwd
#$ -t 1-11
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=24:00:00
#$ -l h_vmem=64G
#$ -o /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/annotations/STalign/logs
#$ -e /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/annotations/STalign/logs

set -euo pipefail

# Initialize conda from the explicit Miniforge installation you provided
MINIFORGE_PATH="/mnt/storage/dept/medonc/beroukhim/youyun/util/miniforge3"
if [ -f "${MINIFORGE_PATH}/etc/profile.d/conda.sh" ]; then
  . "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
else
  echo "Warning: ${MINIFORGE_PATH}/etc/profile.d/conda.sh not found; will try Miniforge env python directly." >&2
fi
conda activate spatial

CODE_ROOT="/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code"
IPYNB_PRIMARY_DIR="${CODE_ROOT}/annotations/STalign/ipynb_generated"
IPYNB_FALLBACK_DIR="${CODE_ROOT}/annotations/STalign/legacy_ipynb"
HTML_DIR="${CODE_ROOT}/annotations/STalign/ipynb_rendered"
mkdir -p "${HTML_DIR}"

NOTEBOOK_BASENAMES=(
  "pathology_CL42112_proseg"
  "pathology_CL67068_proseg"
  "pathology_CL87352_proseg"
  "pathology_CN51420_proseg"
  "pathology_CN53398_proseg"
  "pathology_GG240468_proseg"
  "pathology_GG241670_proseg"
  "pathology_PA188570_proseg"
  "pathology_PA258482_proseg"
  "pathology_PA320784_proseg"
  "pathology_PA328368_proseg"
)

TASK_INDEX=$((SGE_TASK_ID - 1))
if [ "${TASK_INDEX}" -lt 0 ] || [ "${TASK_INDEX}" -ge "${#NOTEBOOK_BASENAMES[@]}" ]; then
  echo "ERROR: SGE_TASK_ID=${SGE_TASK_ID} is out of range for ${#NOTEBOOK_BASENAMES[@]} notebooks." >&2
  exit 1
fi

BASENAME="${NOTEBOOK_BASENAMES[${TASK_INDEX}]}"
IPYNB_FILE="${IPYNB_PRIMARY_DIR}/${BASENAME}.ipynb"
if [ ! -f "${IPYNB_FILE}" ]; then
  IPYNB_FILE="${IPYNB_FALLBACK_DIR}/${BASENAME}.ipynb"
fi
if [ ! -f "${IPYNB_FILE}" ]; then
  echo "ERROR: No notebook found for ${BASENAME}. Looked in:" >&2
  echo "  ${IPYNB_PRIMARY_DIR}" >&2
  echo "  ${IPYNB_FALLBACK_DIR}" >&2
  exit 1
fi

echo "[$(date)] Host: $(hostname)"
echo "[$(date)] Task ${SGE_TASK_ID} executing: ${IPYNB_FILE}"
echo "[$(date)] CONDA_PREFIX: ${CONDA_PREFIX:-unset}"
echo "[$(date)] JUPYTER_BIN: $(which jupyter)"

echo "[$(date)] Executing/rendering via nbconvert -> ${HTML_DIR}"
jupyter nbconvert \
  --to html \
  --execute \
  --ExecutePreprocessor.timeout=-1 \
  --output "${BASENAME}" \
  --output-dir "${HTML_DIR}" \
  "${IPYNB_FILE}"

echo "[$(date)] Completed task ${SGE_TASK_ID}: ${IPYNB_FILE}"
