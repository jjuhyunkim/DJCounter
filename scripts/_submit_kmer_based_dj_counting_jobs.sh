#! /bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: sh _submit_kmer_based_dj_counting_jobs.sh sample_input.map"
  echo "  sample_input.map: Tab-delimited file with sample name and input file(s) in one line"
  exit -1
fi

NUM_LINES=`wc -l $1 | awk '{print $1}'`

cpus=16
mem=120g
name=kmer_dj_count
script=$tools/DJCounter/scripts/_kmer_based_dj_counting_job.sh
args=$1
partition=quick
walltime=4:00:00
path=`pwd`
local="--gres=lscratch:120"
array="--array=1-$NUM_LINES"

mkdir -p logs
log=logs/$name.%A_%a.log

set -x
sbatch -J $name \
  --cpus-per-task=$cpus --mem=$mem \
  --partition=$partition \
  $array $local \
  -D $path --time=$walltime --error=$log \
  --output=$log $script $args

