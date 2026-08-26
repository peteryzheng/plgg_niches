#!/bin/bash -l
#$ -N banksy_cohort
#$ -cwd
#$ -t 1
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=120:00:00
#$ -l h_vmem=96G
#$ -o /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/niche/banksy_workflow/logs
#$ -e /mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/niche/banksy_workflow/logs

# task 1 = current chosen config (k_geom 15/30, lam 0.2/0.8, npc 20,
# k_ct=50, res_ct=0.2, k_ni=50, res_ni=0.5) -- the only row in
# param_search.tsv. Older exploratory sweep rows used a different,
# incompatible column layout (kc1/kc2/res1/res2 shared across both lambdas)
# and were removed rather than migrated; see git history if needed.
#
# res_ct=0.2/res_ni=0.5 were chosen after an extensive resolution sweep
# (see cluster_summary_stats.R and celltype_lowres_recluster.R) -- the
# production annotated object was actually assembled via a validated
# fast-path (finalize_final_resolution.R) that reuses the deterministic
# BANKSY/Harmony/UMAP embeddings + Leiden clustering already computed for
# this config rather than re-running this qsub script from scratch (which
# would take ~35-55h to reproduce a scientifically identical result). This
# script remains the correct from-scratch fallback if that equivalence
# ever needs re-validating or the upstream data changes.

set -euo pipefail

# Activate the spatial conda env via the explicit Miniforge install on the
# new server (matches annotations/STalign/stalign_qmd_array_qsub.sh).
MINIFORGE_PATH="/mnt/storage/dept/medonc/beroukhim/youyun/util/miniforge3"
if [ -f "${MINIFORGE_PATH}/etc/profile.d/conda.sh" ]; then
    . "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
else
    echo "Warning: ${MINIFORGE_PATH}/etc/profile.d/conda.sh not found; relying on PATH conda." >&2
fi
conda activate spatial

CODE_ROOT="/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code"
DATA_ROOT="/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data"

# Pull one row from param_search.tsv per array task. Columns are
# whitespace-delimited; res_ct / res_ni are comma-separated lists (no
# spaces) so awk treats each as a single token, e.g. "0.5,1".
parameter_file_path="${CODE_ROOT}/niche/banksy_workflow/param_search.tsv"
line=$(sed -n -e "${SGE_TASK_ID}p" "${parameter_file_path}")
echo "[$(date)] Host: $(hostname)"
echo "[$(date)] Task ${SGE_TASK_ID} param line: ${line}"
# columns: k1 k2 lambda1 lambda2 npcs k_ct res_ct k_ni res_ni
# example: 15 30 0.2 0.8 20 50 0.5,1 50 0.5,1
k1=$(echo "${line}" | awk '{print $1}')
k2=$(echo "${line}" | awk '{print $2}')
lambda1=$(echo "${line}" | awk '{print $3}')
lambda2=$(echo "${line}" | awk '{print $4}')
npcs=$(echo "${line}" | awk '{print $5}')
k_ct=$(echo "${line}" | awk '{print $6}')
res_ct=$(echo "${line}" | awk '{print $7}')
k_ni=$(echo "${line}" | awk '{print $8}')
res_ni=$(echo "${line}" | awk '{print $9}')

outputdir="${DATA_ROOT}/banksy_param_search/k1_${k1}_k2_${k2}_lambda1_${lambda1}_lambda2_${lambda2}_npcs_${npcs}_kct_${k_ct}_resct_${res_ct}_kni_${k_ni}_resni_${res_ni}"
mkdir -p "${outputdir}"

echo "[$(date)] CONDA_PREFIX: ${CONDA_PREFIX:-unset}"
echo "[$(date)] Rscript: $(which Rscript)"
echo "[$(date)] Output dir: ${outputdir}"

Rscript "${CODE_ROOT}/niche/banksy_workflow/banksy_cohort.R" \
    --k1 "${k1}" --k2 "${k2}" --lam1 "${lambda1}" --lam2 "${lambda2}" --npc "${npcs}" \
    --k_ct "${k_ct}" --res_ct "${res_ct}" --k_ni "${k_ni}" --res_ni "${res_ni}" \
    --seed 55555 \
    -o "${outputdir}"

echo "[$(date)] Completed task ${SGE_TASK_ID}"
