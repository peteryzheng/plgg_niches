#!/bin/bash -l
#$ -N extract_metadata
#$ -cwd
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=4:00:00
#$ -l h_vmem=64G
#$ -o /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/healthy_comparisons/logs
#$ -e /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/healthy_comparisons/logs

# Extract .obs metadata from the 35G integrated h5ad into a parquet file.
# Only .obs is read into RAM (backed mode); 64G is well above the required headroom.

set -euo pipefail

MINIFORGE_PATH="/mnt/storage/dept/medonc/beroukhim/youyun/util/miniforge3"
if [ -f "${MINIFORGE_PATH}/etc/profile.d/conda.sh" ]; then
    . "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
else
    echo "Warning: ${MINIFORGE_PATH}/etc/profile.d/conda.sh not found; relying on PATH conda." >&2
fi
conda activate spatial

CODE_ROOT="/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code"

echo "[$(date)] Host: $(hostname)"
echo "[$(date)] CONDA_PREFIX: ${CONDA_PREFIX:-unset}"
echo "[$(date)] python: $(which python)"

mkdir -p "${CODE_ROOT}/healthy_comparisons/logs"

python "${CODE_ROOT}/healthy_comparisons/extract_metadata.py"

echo "[$(date)] Done."
