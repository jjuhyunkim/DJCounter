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
BED files are provided in this repository (link [here](https://github.com/jjuhyunkim/DJCounter/raw/main/ukb-dj.bed)).
Aligned BAM files should be based on the Broad GRCh38 reference. You can download the reference from [Broad Github](https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz).
```bash
# This is an example
sample=Sample01
threads=10

bam="/data/01.broad_hg38/$sample/$sample.dedup.bam" # BAM or CRAM file
outdir="/data/01.broad_hg38/$sample" # The output directory
prefix="$sample" # Prefix for output files
bed="/data/01.broad_hg38/uk.dj.bed" # BED file for DJ counting
```

### Calculating the bacgkround coverage 
We calculate background coverage using autosomal chromosomes from the GRCh38 reference.
This is achieved by analyzing aligned BAM files with tools like samtools depth or samtools idxstats.
Ensure that your BAM files contain the necessary contigs for accurate DJ coverage calculation.
Background coverage is computed as the median depth across autosomal chromosomes.
```bash
## Calculate the read length
fragmentSize=$(samtools view $bam | head -10000 | awk '{print length($10)}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
echo $fragmentSize

## Calculate the total background
bgCov=$(samtools coverage $bam | grep -E '^chr(1?[0-9]|2[0-2])\s' | awk '{print $7}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
echo $bgCov
```

### Calculating the coverage of DJ regions scattered across GRCh38 reference
The provided BED file contains regions highly similar to DJ regions on CHM13 chromosome 13, excluding high variable or repeat regions.
Coverage calculation utilizes `samtools depth`, with computational time typically under 5~8 minutes using 10 threads, depending on BAM file size.
```bash
# Calculate the DJ regions' coverage
sum=$(samtools depth -@ $threads -b $bed $bam | awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }')
echo $sum
```

### Calculating the number of DJ counts
DJ counts are derived by dividing the total depth on DJ regions in GRCh38 by the background coverage.
The total length of DJ regions used in this analysis is fixed at 136,405 bp.
Results are presented in diploid genome bases by multiplying by 2.
```bash
covLen=136405
djCount=$(echo "scale=5; 2 * $sum / $covLen / $bgCov" | bc)
echo -e "$prefix\t$djCount" > $outdir/$prefix.dj_hg38.txt
```

###  Expectable output
The output file `$outdir/$prefix.dj_hg38.txt` contains two columns separated by tabs:
1. **$prefix:** This column contains identifiers or names associated with each estimation of DJ count.
2. **Estimation of DJ count based on diploid genome:** This column provides the calculated DJ count values adjusted for diploid genome context.
```bash
Sample01	8.30800
```

Normal human samples typically yield around 10 copies of DJ counts, with occasional deviations to ~11 or ~9.
Robertsonian samples usually show approximately ~8 copies.

<img src="https://github.com/user-attachments/assets/9212dabb-593f-4040-bebc-494a74301fa0" width="200">


## Changing logs
V0.1(2024-07-17)
* first commit

V0.2(2024-07-25)
* Changing the background coverage estimation methods from samtools idxstats to samtools coverage.
* Removing the step of saving temporary files; instead, we assign everything to variables.
