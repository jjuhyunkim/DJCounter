#! /bin/bash
sample=$1
prefix=$2
bam=$3

threads=10
covLen=136405
bed=ukb-dj.bed
outFile=out.txt

ml samtools

fragmentSize=$(samtools view -@ $threads $bam | head -10000 | awk '{print length($10)}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 

dj=$(samtools view -@ $threads -c -L $bed $bam)

bgCov=$(samtools idxstats -@ $threads $bam | grep -E '^chr(1?[0-9]|2[0-2])\s'| awk -v frag=$fragmentSize '{print $3/$2*frag}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 

djCount=$(echo "scale=5; 2 * $dj * $fragmentSize / $covLen / $bgCov" | bc)

echo -e "${sample}_${prefix}\t$fragmentSize\t$bgCov\t$dj\t$djCount" >> $outFile
