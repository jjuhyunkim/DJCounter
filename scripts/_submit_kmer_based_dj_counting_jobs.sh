#! /bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: sh _submit_kmer_based_dj_counting_jobs.sh sample_input.map [-ref hg38] [-array array_idx]"
  echo "  sample_input.map: Tab-delimited file with sample name and input file(s) in one line"
  echo "                    Format: <sample_name> <input.bam|input.cram|input.fq.gz>"
  echo "  -ref: (optional) Reference genome used in the cram input file. Use hg38 for GRCh38_full_analysis_set_plus_decoy_hla.fa"
  echo "        If not provided, assumes cram has all the sequences."
  echo "  -array: (optional) Comma separated list of line numbers to submit as array job."
  exit -1
fi


# Default values
map=$1
shift 1

ref=""
array_idx=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -ref)
            ref="$2"
            shift 2
            ;;
        -array)
            array_idx="$2"
            shift 2
            ;;
        *)
            # Handle non-option arguments or break the loop
            echo "Unknown option: $1"
            exit 1
    esac
done


args="$map"
if [ -n "$ref" ]; then
    echo "Reference genome: $ref"
    args="$map $ref"
fi


NUM_LINES=`wc -l $map | awk '{print $1}'`
array="--array=1-$NUM_LINES"
if [ -n "$array_idx" ]; then
    echo "Array job indices: $array_idx"
    array="--array=$array_idx"
fi

cpus=16
mem=120g
name=kmer_dj_count
script=$tools/DJCounter/scripts/_kmer_based_dj_counting_job.sh
args="$args"
partition=quick
walltime=4:00:00
path=`pwd`
local="--gres=lscratch:300"

mkdir -p logs
log=logs/$name.%A_%a.log

set -x
sbatch -J $name \
  --cpus-per-task=$cpus --mem=$mem \
  --partition=$partition \
  $array $local \
  -D $path --time=$walltime --error=$log \
  --output=$log $script $args

