#!/bin/bash
#$ -l h_rt=120:00:00
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=64G
#$ -o '/xchip/beroukhimlab/youyun/plgg/code/niche/logs'
#$ -e '/xchip/beroukhimlab/youyun/plgg/code/niche/logs'
#$ -N banksy_cohort

set -e

export PATH="/xchip/beroukhimlab/youyun/miniconda3/bin:$PATH"

/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n spatial --live-stream \
    Rscript /xchip/beroukhimlab/youyun/plgg/code/niche/banksy_cohort.R 
