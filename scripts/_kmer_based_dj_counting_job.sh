#!/bin/sh

map=$1
if [[ -z $map ]]; then
  echo "Usage: _kmer_based_dj_counting_job.sh <map_file> [-ref hg38] [-array array_task_id]"
  echo "  map_file: Tab-delimited file with sample names and input files"
  echo "  -ref: GRCh38_full_analysis_set_plus_decoy_hla.fa will be used if ref is hg38."
  echo "  -array: SLURM array task ID to select the line from the map file"
  exit 1
fi

shift 1
ref=""
i=$SLURM_ARRAY_TASK_ID

while [[ $# -gt 0 ]]; do
    case "$1" in
        -ref)
            ref="$2"
            shift 2
            ;;
        -array)
            i="$2"
            shift 2
            ;;
        *)
            # Handle non-option arguments or break the loop
            echo "Unknown option: $1"
            exit 1
    esac
done


if [[ -z $i ]]; then
  echo "No SLURM_ARRAY_TASK_ID provided or -array provided. Exit."
fi

sample=$(sed -n "${i}p" $map | awk '{print $1}')
input=$(sed -n "${i}p" $map | awk '{print $2}')

$tools/DJCounter/scripts/kmer_based_dj_counting.sh $sample $input $ref
