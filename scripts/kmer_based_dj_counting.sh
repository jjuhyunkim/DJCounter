#!/bin/sh

sample=$1
input=$2
ref=$3
DJ_TARGET=$tools/DJCounter/resources/DJtarget.meryl

if [[ -z $sample || -z $input ]]; then
  echo "Usage: kmer_based_dj_counting.sh <sample_name> <input.bam|input.cram|input.fq.gz> [ref]"
  echo "  sample_name: Sample identifier"
  echo "  input.bam|input.cram|input.fq.gz: Input sequencing reads in BAM or FASTQ format (gz or not)."
  echo "  For paired-end reads, provide files as a comma separated list e.g. \"input1.fq.gz,input2.fq.gz\""
  echo "  ref: (Optional) Reference genome used in the bam/cram input file."
  exit 1
fi

cpus=$SLURM_CPUS_PER_TASK
if [ -z "$cpus" ]; then
  cpus=24
fi

if [[ -z $SLURM_MEM_PER_NODE ]]; then
  mem=48
else
  # Convert MB to GB
  mem=$(((SLURM_MEM_PER_NODE/1024)-10))
fi

if [[ -z $SLURM_JOB_ID ]]; then
  tmp="."
else
  tmp="/lscratch/$SLURM_JOB_ID"
fi

if [[ -s DJcounts/${sample}_DJ_count.txt ]]; then
  echo "DJ count file already exists for ${sample}. Nothing to do."
  exit 0
fi

set -e
set -o pipefail
set -x

mkdir -p hist DJcounts

if [[ -d $tmp/${sample}.k31.meryl ]]; then
  echo "Kmer database already exists for ${sample}, skipping counting."
else
  input=$(echo $input | tr ',' ' ')
  echo "Counting kmers for ${sample} from ${input}"
  if [[ "${input}" == *".bam" || "${input}" == *".cram" ]]; then
    echo "Input file is a BAM/CRAM file."
	  if ! [[ -s $input ]]; then
	    echo "$input file is empty."
	    exit 1
	  fi
    module load samtools/1.21
    ref=""
    if [[ -z $ref ]] || [[ $ref == "" ]]; then
      echo "No reference provided for BAM/CRAM input. Assuming bam/cram has all the sequences."
    elif [[ $ref == "hg38" ]]; then
      ref="--reference $tools/DJCounter/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    fi
    # mkfifo ${sample}_fq_pipe
    # $samtools fastq -@ ${cpus} $ref $input > ${sample}_fq_pipe &
    samtools fastq -@ ${cpus} $ref $input | pigz -c - > $tmp/${sample}.fq.gz
	  # meryl count k=31 threads=${cpus} memory=${mem} output ${sample}.k31.meryl ${sample}_fq_pipe || exit -1
    meryl count k=31 threads=${cpus} memory=${mem} output $tmp/${sample}.k31.meryl $tmp/${sample}.fq.gz || exit -1
    #rm ${sample}_fq_pipe
  else
    meryl count k=31 threads=${cpus} memory=${mem} output $tmp/${sample}.k31.meryl ${input}
  fi
  meryl histogram $tmp/${sample}.k31.meryl > hist/${sample}.k31.hist
fi

if [[ -s hist/${sample}.DJ.hist ]]; then
  echo "DJ meryl already exists for ${sample}, skipping intersect."
else
  meryl intersect threads=${cpus} memory=${mem} \
    $tmp/${sample}.k31.meryl $DJ_TARGET output hist/${sample}.DJ.meryl
fi
meryl histogram hist/${sample}.DJ.meryl  > hist/${sample}.DJ.hist

count_mid=`cat hist/${sample}.DJ.hist | awk '{ count+=$NF; } END {print count/2}'`
med=`cat hist/${sample}.DJ.hist | \
  awk -v count_mid=${count_mid} '{cnt_sum+=$NF; if (cnt_sum > count_mid) {print $(NF-1); exit;} }'`

peak2=`java -jar -Xmx256m $MERQURY/eval/kmerHistToPloidyDepth.jar hist/${sample}.k31.hist |\
  tail -n1 | awk '{print $2}'`

# print Sample, DJmedCov, PeakCP2, Peak_Est
echo -e "${sample}\t${med}\t${peak2}" |
  awk -F "\t" '{print $1"\t"$2"\t"$3"\t"(($2*2)/$3)}' \
  > DJcounts/${sample}_DJ_count.txt

cat DJcounts/${sample}_DJ_count.txt
