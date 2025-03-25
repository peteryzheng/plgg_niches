#! /usr/bin/bash

transcript_gz_path=$1
output_dir=$2
thread_num=$3

# parent directory of the transcript gz file
original_xenium_directory=$(dirname $transcript_gz_path)
project_id=$(basename $original_xenium_directory)

# Create output directory
mkdir -p $output_dir
cd $output_dir

# Running proseg in output directory
proseg --xenium --nthreads $thread_num $transcript_gz_path

proseg-to-baysor \
    transcript-metadata.csv.gz \
    cell-polygons.geojson.gz \
    --output-transcript-metadata proseg-transcript-metadata-baysor-import.csv \
    --output-cell-polygons proseg-cell-polygons-baysor-import.geojson

# interestingly this works on a compute node interactive session but not on dipg
/xchip/beroukhimlab/youyun/util/xeniumranger/xeniumranger-xenium3.0/xeniumranger \
    import-segmentation --id "${project_id}_proseg" \
    --xenium-bundle $original_xenium_directory \
    --viz-polygons proseg-cell-polygons-baysor-import.geojson \
    --transcript-assignment proseg-transcript-metadata-baysor-import.csv \
    --units microns