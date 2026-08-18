#!/bin/bash
#$ -l h_rt=120:00:00
#$ -t 109
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=96G
#$ -o '/xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/logs'
#$ -e '/xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/logs'
#$ -N banksy_cohort

# task 109 is default param (k1=15,k2=30,lam1=0.2,lam2=0.8,npc=20 -- same as
# the old task 14 default -- but with per-lambda Leiden k/res picked from the
# evaluate_leiden_sweep.qmd analysis instead of one shared kc/res pair)
set -e

export PATH="/xchip/beroukhimlab/youyun/miniconda3/bin:$PATH"

# go into /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/param_search.tsv to get each individual parameter
parameter_file_path=/xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/param_search.tsv
line=$(sed -n -e "${SGE_TASK_ID}p" $parameter_file_path)
echo $line
# k1  k2  lambda1  lambda2  npcs  k_leiden_lam1  k_leiden_lam2  res_lam1  res_lam2
# k_leiden_lam1/k_leiden_lam2/res_lam1/res_lam2 are comma-separated (e.g.
# "0.1,0.2") since cell typing (lam1) and niche calling (lam2) can each test
# a different set of Leiden k-neighbors / resolutions rather than sharing one
# grid across both lambdas.
# 5	10	0.1	0.8	20	30,50	30,50	0.75,1	0.75,1
k1=$(echo $line | awk '{print $1}')
k2=$(echo $line | awk '{print $2}')
lambda1=$(echo $line | awk '{print $3}')
lambda2=$(echo $line | awk '{print $4}')
npcs=$(echo $line | awk '{print $5}')
k_leiden_lam1=$(echo $line | awk '{print $6}')
k_leiden_lam2=$(echo $line | awk '{print $7}')
res_lam1=$(echo $line | awk '{print $8}')
res_lam2=$(echo $line | awk '{print $9}')


outputdir=/xchip/beroukhimlab/youyun/plgg/data/banksy_param_search/k1_${k1}_k2_${k2}_lambda1_${lambda1}_lambda2_${lambda2}_npcs_${npcs}_kleidenlam1_${k_leiden_lam1}_kleidenlam2_${k_leiden_lam2}_reslam1_${res_lam1}_reslam2_${res_lam2}
mkdir -p $outputdir

# echo the run command

echo "/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n spatial --live-stream \
	Rscript /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/banksy_cohort.R \
	--k1 $k1 --k2 $k2 --lam1 $lambda1 --lam2 $lambda2 --npc $npcs \
	--k_leiden_lam1 $k_leiden_lam1 --k_leiden_lam2 $k_leiden_lam2 \
	--res_lam1 $res_lam1 --res_lam2 $res_lam2 \
	--seed 55555 \
	-o $outputdir"

/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n spatial --live-stream \
    Rscript /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/banksy_cohort.R \
	--k1 $k1 --k2 $k2 --lam1 $lambda1 --lam2 $lambda2 --npc $npcs \
	--k_leiden_lam1 "$k_leiden_lam1" --k_leiden_lam2 "$k_leiden_lam2" \
	--res_lam1 "$res_lam1" --res_lam2 "$res_lam2" \
	--seed 55555 \
	-o $outputdir
    



