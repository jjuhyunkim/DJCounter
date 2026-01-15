#!/bin/sh

sample=$1
input=$2
DJ_TARGET=$tools/DJcounter/resources/DJtarget.meryl

if [[ -z $sample || -z $input ]]; then
  echo "Usage: kmer_based_dj_counting.sh <sample_name> <input.bam|input.fq.gz>"
  echo "  sample_name: Sample identifier"
  echo "  input.bam|input.fq.gz: Input sequencing reads in BAM or FASTQ format (gz or not)."
  echo "  For paired-end reads, provide files as a comma separated list e.g. \"input1.fq.gz,input2.fq.gz\""
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
  mem=$((SLURM_MEM_PER_NODE/1024))
fi

if [[ -s ${sample}_DJ_count.txt ]]; then
  echo "DJ count file already exists for ${sample}. Nothing to do."
  echo
  exit 0
fi

set -e
set -o pipefail

if [[ -d $sample.k31.meryl ]]; then
  echo "Kmer database already exists for ${sample}, skipping counting."
else
  input=$(echo $input | tr ',' ' ')
  echo "Counting kmers for ${sample} from ${input}"
  if [[ "${input}" == *".bam" ]]; then
    echo "Input file is a BAM file."
	  if ! [[ -s $input ]]; then
	    echo "BAM file is empty."
	    exit 1
	  fi
    mkfifo ${sample}_fq_pipe
    module load samtools/1.21
    samtools fastq -@ ${cpus} $input > ${sample}_fq_pipe &
	  meryl count k=31 threads=${cpus} memory=${mem} output ${sample}.k31.meryl ${sample}_fq_pipe || exit -1
    rm ${sample}_fq_pipe
  else
    meryl count k=31 threads=${cpus} memory=${mem} output ${sample}.k31.meryl ${input}
  fi
  meryl histogram ${sample}.k31.meryl > ${sample}.k31.hist
fi

if [[ -s ${sample}.DJ.meryl ]]; then
  echo "DJ meryl already exists for ${sample}, skipping intersect."
else
  meryl intersect threads=${cpus} memory=${mem} ${sample}.k31.meryl $DJ_TARGET output ${sample}.DJ.meryl
fi
meryl histogram ${sample}.DJ.meryl  > ${sample}.DJ.hist

count_mid=`cat ${sample}.DJ.hist | awk '{ count+=$NF; } END {print count/2}'`
med=`cat ${sample}.DJ.hist | \
  awk -v count_mid=${count_mid} '{cnt_sum+=$NF; if (cnt_sum > count_mid) {print $(NF-1); exit;} }'`

peak2=`java -jar -Xmx256m $MERQURY/eval/kmerHistToPloidyDepth.jar ${sample}.k31.hist | tail -n1 | awk '{print $2}'`

# print Sample, DJmedCov, PeakCP2, Peak_Est
echo -e "${sample}\t${med}\t${peak2}" |
  awk -F "\t" '{print $1"\t"$2"\t"$3"\t"(($2*2)/$3)}' \
  > ${sample}_DJ_count.txt
