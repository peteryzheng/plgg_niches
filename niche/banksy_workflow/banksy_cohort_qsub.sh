#!/bin/bash
#$ -l h_rt=120:00:00
#$ -t 14
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_vmem=96G
#$ -o '/xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/logs'
#$ -e '/xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/logs'
#$ -N banksy_cohort

#task 14 is default param
set -e

export PATH="/xchip/beroukhimlab/youyun/miniconda3/bin:$PATH"

# go into /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/param_search.tsv to get each individual parameter
parameter_file_path=/xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/param_search.tsv
line=$(sed -n -e "${SGE_TASK_ID}p" $parameter_file_path)
echo $line
# k1    k2 lambda1 lambda2  npcs   kc1   kc2      res1  res2
# 5	10	0.1	0.8	20	20	30	0.75	1
# OR
# 5 10 0.1 0.8 20 20 30 0.75 1
k1=$(echo $line | awk '{print $1}')
k2=$(echo $line | awk '{print $2}')
lambda1=$(echo $line | awk '{print $3}')
lambda2=$(echo $line | awk '{print $4}')
npcs=$(echo $line | awk '{print $5}')
kc1=$(echo $line | awk '{print $6}')
kc2=$(echo $line | awk '{print $7}')
res1=$(echo $line | awk '{print $8}')
res2=$(echo $line | awk '{print $9}')


outputdir=/xchip/beroukhimlab/youyun/plgg/data/banksy_param_search/k1_${k1}_k2_${k2}_lambda1_${lambda1}_lambda2_${lambda2}_npcs_${npcs}_kc1_${kc1}_kc2_${kc2}_res1_${res1}_res2_${res2}
mkdir -p $outputdir

# echo the run command

echo "/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n spatial --live-stream \
	Rscript /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/banksy_cohort.R \
	--k1 $k1 --k2 $k2 --lam1 $lambda1 --lam2 $lambda2 --npc $npcs \
	--kc1 $kc1 --kc2 $kc2 --res1 $res1 --res2 $res2 \
	--seed 55555 \
	-o $outputdir"

/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n spatial --live-stream \
    Rscript /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_workflow/banksy_cohort.R \
	--k1 $k1 --k2 $k2 --lam1 $lambda1 --lam2 $lambda2 --npc $npcs \
	--kc1 $kc1 --kc2 $kc2 --res1 $res1 --res2 $res2 \
	--seed 55555 \
	-o $outputdir
    



