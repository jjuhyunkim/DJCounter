# Mapping based approach
1. Calculae background coverage
2. Calculae coverage of DJ target regions scattered across a given reference
3. Calculae the number of DJs

## Input files and settings
Target DJ region used in the assessment are provided in this repository for
 [GRCh38/hg38](https://github.com/jjuhyunkim/DJCounter/blob/main/roi/GRCh38/ukb-dj.bed) or
 [GRCh37/hg19](https://github.com/jjuhyunkim/DJCounter/blob/main/roi/hg19/DJ.bed).

For GRCh38, the DJ sequences are on `chr21`, `chr17_GL000205v2_random` and `chrUn_GL000195v1`. Our pipeline is optimized for processing aligned BAM files to the Broad reference version of the GRCh38. You can download the reference from [Broad Github](https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz).

```bash
# This is an example
sample=Sample01
threads=10

bam="/data/01.broad_hg38/$sample/$sample.dedup.bam" # BAM or CRAM file
outdir="/data/01.broad_hg38/$sample" # The output directory
prefix="$sample" # Prefix for output files
bed="/data/01.broad_hg38/uk.dj.bed" # BED file for DJ counting
```

## Bacgkround coverage
We calculate background coverage using autosomal chromosomes from the reference.
This is achieved by analyzing aligned BAM files with `samtools coverage` and `samtools depth`.
Ensure that your BAM files contain the necessary contigs for accurate DJ coverage calculation.
Background coverage is computed as the median depth across autosomal chromosomes.

```bash
## Get the read length from the first 10k reads
fragmentSize=$(samtools view $bam | head -10000 | awk '{print length($10)}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}')
echo $fragmentSize

## Calculate the total background
bgCov=$(for i in {1..22}; do
  samtools coverage -r chr${i} "$bam" | sed -n 2p | awk -v frag="$fragmentSize" '{print $4/$3*frag}'
done | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]};')
echo $bgCov
```

## Target DJ coverage
The provided BED file contains regions highly similar to DJ regions on CHM13 chromosome 13, excluding high variable or repeat regions.
Coverage calculation utilizes `samtools depth`, with computational time typically under 5~8 minutes using 10 threads, depending on BAM file size.
```bash
# Calculate the DJ regions' coverage
sum=$(samtools depth -@ $threads -b $bed $bam | awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }')
echo $sum
```

### Copy number estimate of the DJ
DJ counts are derived by dividing the total depth on DJ regions in GRCh38 by the background coverage.
The total length of DJ regions used in this analysis is fixed at 136,405 bp.
Results are presented in diploid genome bases by multiplying by 2.
```bash
covLen=136405
djCount=$(echo "scale=5; 2 * $sum / $covLen / $bgCov" | bc)
echo -e "$prefix\t$fragmentSize\t$bgCov\t$sum\t$djCount" > $outdir/$prefix.dj_hg38.txt
```

###  Expected output
The output file `$outdir/$prefix.dj_hg38.txt` contains two columns separated by tabs:
1. **$prefix:** This column contains identifiers or names associated with each estimation of DJ count.
2. **Estimation of DJ count based on diploid genome:** This column provides the calculated DJ count values adjusted for diploid genome context.
```bash
Sample01	8.30800
```

Normal human samples typically yield around 10 copies of DJ counts, with occasional deviations to ~11 or ~9.
Robertsonian samples usually show approximately ~8 copies.

<img src="https://github.com/user-attachments/assets/9212dabb-593f-4040-bebc-494a74301fa0" width="200">

