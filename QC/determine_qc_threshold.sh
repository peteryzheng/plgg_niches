#!/bin/bash
#$ -l h_vmem=96G
#$ -l h_rt=48:00:00
#$ -pe smp 1
#$ -binding linear:1

export PATH="/xchip/beroukhimlab/youyun/miniconda3/bin:$PATH"

/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n spatial --live-stream \
    quarto render /xchip/beroukhimlab/youyun/plgg/code/QC/determine_qc_thresholds.qmd