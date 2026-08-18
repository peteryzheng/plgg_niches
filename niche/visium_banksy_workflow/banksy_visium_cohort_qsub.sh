#!/bin/bash -l
#$ -l h_rt=120:00:00
#$ -cwd
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=128G
#$ -o /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/niche/visium_banksy_workflow/logs
#$ -e /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/niche/visium_banksy_workflow/logs
#$ -N visium_banksy

set -euo pipefail

MINIFORGE_PATH="/mnt/storage/dept/medonc/beroukhim/youyun/util/miniforge3"
if [ -f "${MINIFORGE_PATH}/etc/profile.d/conda.sh" ]; then
  . "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
else
  echo "Warning: ${MINIFORGE_PATH}/etc/profile.d/conda.sh not found; will try env binaries directly." >&2
fi
conda activate spatial

CODE_ROOT="/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code"
DATA_ROOT="/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data"
SCRIPT_ROOT="${CODE_ROOT}/niche/visium_banksy_workflow"
LOG_ROOT="${SCRIPT_ROOT}/logs"

inputdir=${1:-${DATA_ROOT}/visium/rds_spe}
outputdir=${2:-${DATA_ROOT}/visium/banksy}

mkdir -p "${LOG_ROOT}"
mkdir -p "$outputdir"

echo "[$(date)] Host: $(hostname)"
echo "[$(date)] CONDA_PREFIX: ${CONDA_PREFIX:-unset}"
echo "[$(date)] RSCRIPT_BIN: $(which Rscript)"
echo "[$(date)] Input dir: ${inputdir}"
echo "[$(date)] Output dir: ${outputdir}"

echo "[$(date)] Running: Rscript ${SCRIPT_ROOT}/banksy_visium_cohort.R \
    --inputdir ${inputdir} \
    --outputdir ${outputdir} \
    --k_geom 18 \
    --npc 20 \
    --lambdas 0,0.2 \
    --k_neighbors 30 \
    --resolution 0.8 \
    --seed 55555"

Rscript "${SCRIPT_ROOT}/banksy_visium_cohort.R" \
    --inputdir "$inputdir" \
    --outputdir "$outputdir" \
    --k_geom 18 \
    --npc 20 \
    --lambdas 0,0.2 \
    --k_neighbors 30 \
    --resolution 0.8 \
    --seed 55555

echo "[$(date)] Completed Visium BANKSY cohort run"
