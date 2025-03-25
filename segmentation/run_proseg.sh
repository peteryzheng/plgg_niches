#!/bin/bash
#$ -l h_rt=48:00:00
#$ -t 2-11
#$ -pe smp 8
#$ -binding linear:8
#$ -l h_vmem=64G
#$ -o '/xchip/beroukhimlab/youyun/plgg/data/segmentation/logs'
#$ -e '/xchip/beroukhimlab/youyun/plgg/data/segmentation/logs'
#$ -N proseg

export PATH="/home/unix/youyun/.cargo/bin:$PATH"

# the #SGE_TASK_ID# is the task id
# that line in the file /xchip/beroukhimlab/youyun/plgg/code/segmentation/transcripts.txt is the input transcript gz path
transcript_gz_path=$(sed -n "$SGE_TASK_ID p" /xchip/beroukhimlab/youyun/plgg/code/segmentation/transcripts.txt)
project_id=$(basename $(dirname $transcript_gz_path))
output_dir='/xchip/beroukhimlab/youyun/plgg/data/segmentation/proseg_run_121024'
mkdir -p $output_dir/$project_id

bash /xchip/beroukhimlab/youyun/plgg/code/segmentation/proseg.sh \
	$transcript_gz_path $output_dir/$project_id 8
