#!/bin/sh

map=$1
i=$SLURM_ARRAY_TASK_ID

if [[ -z $SLURM_ARRAY_TASK_ID ]]; then
  i=$2
fi

if [[ -z $map || -z $i ]]; then
  echo "Usage: _kmer_based_dj_counting_job.sh <map_file> <array_task_id>"
  echo "  map_file: Tab-delimited file with sample names and input files"
  echo "  array_task_id: SLURM array task ID to select the line from the map file"
  exit 1
fi

sample=$(sed -n "${i}p" $map | awk '{print $1}')
input=$(sed -n "${i}p" $map | awk '{print $2}')

$tools/DJCounter/scripts/kmer_based_dj_counting.sh $sample $input hg38
