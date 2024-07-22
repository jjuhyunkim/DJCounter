# DJ Counter

## Purpose
This repository provides tools for calculating the number of DJ counts in the genome aligned to the GRCh38 (broad reference).

## Principle
The DJ counts are determined by assessing the coverage of specific regions across the GRCh38 reference genome.

## Main workflow
1. Calculating the background coverage
2. Calculating the coverage of DJ regions scattered across GRCh38 reference
3. Calculating the number of DJ counts

### Input files and settings
BED files are provided in this repository (link here).
Aligned BAM files should be based on the Broad GRCh38 reference. You can download the reference here.
```
# This is an example
sample=GM04890
threads=10
covLen=136405

bam="/data/01.broad_hg38/$sample/$sample.dedup.bam" # BAM or CRAM file
outdir="/data/01.broad_hg38/$sample" # The output directory
prefix="$sample" # Prefix for output files
bed="/data/01.broad_hg38/uk.dj.bed" # BED file for DJ counting
idxstats="$outdir/$prefix.idxstat" # If you already have one, specify its name here. If not, you can create it.
```

### Calculating the bacgkround coverage 
We calculate background coverage using autosomal chromosomes from the GRCh38 reference.
This is achieved by analyzing aligned BAM files with tools like samtools depth or samtools idxstats.
Ensure that your BAM files contain the necessary contigs for accurate DJ coverage calculation.
Background coverage is computed as the median depth across autosomal chromosomes.
```
## Calculate the read length
fragmentSize=$(samtools view $bam | head -10000 | awk '{print length($10)}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
echo $fragmentSize
# Check the idxstats 
# fragmentSize=$(samtools view $bam | head -10000 | awk '{print length($10)}' | datamash median 1)
if [ ! -f $idxstats ]; then
	samtools idxstats -@ $threads $bam > $idxstats
fi

## Calculate the total background
bgCov=$(cat $idxstats | head -22 | awk -v frag=$fragmentSize '{print $3/$2*frag}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
echo $bgCov
```

### Calculating the coverage of DJ regions scattered across GRCh38 reference
The provided BED file contains regions highly similar to DJ regions on CHM13 chromosome 13, excluding high variable or repeat regions.
Coverage calculation utilizes samtools depth, with computational time typically under 5 minutes using 10 threads, depending on BAM file size.
```
# Calculate the DJ regions' coverage
samtools depth -@ $threads -b $bed $bam > $outdir/$prefix.dj_coverage_results.depth
sum=$(awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }' "$outdir/$prefix.dj_coverage_results.depth")
echo $sum
```

### Calculating the number of DJ counts
DJ counts are derived by dividing the total depth on DJ regions in GRCh38 by the background coverage.
The total length of DJ regions used in this analysis is fixed at 136,405 bp.
Results are presented in diploid genome bases by multiplying by 2.
```
djCount=$(echo "scale=5; 2 * $sum / $covLen / $bgCov" | bc)
echo -e "$prefix\t$djCount" > $outdir/$prefix.dj_hg38.txt
```

###  Expectable output
Normal human samples typically yield around 10 copies of DJ counts.
Robertsonian samples usually show approximately ~8 copies, with occasional deviations to ~11 or ~9.

<img src="https://github.com/user-attachments/assets/9212dabb-593f-4040-bebc-494a74301fa0" width="200">

### Logs




