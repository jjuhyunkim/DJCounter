# Mapping based approach
1. Calculae background coverage
2. Calculae coverage of DJ target regions scattered across a given reference
3. Calculae the number of DJs

## Input files and settings
Target DJ region used in the assessment are provided in this repository for
 [GRCh38/hg38](https://github.com/jjuhyunkim/DJCounter/blob/main/roi/GRCh38/ukb-dj.bed) or
 [GRCh37/hg19](https://github.com/jjuhyunkim/DJCounter/blob/main/roi/hg19/DJ.bed).

For GRCh38, the DJ sequences are on `chr21`, `chr17_GL000205v2_random` and `chrUn_GL000195v1`. Our pipeline is optimized for processing aligned BAM files to the Broad reference version of the GRCh38. You can download the reference from [Broad Github](https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz).


## Wrapper script 
```bash
#############################################################################
.DDD..   .JJJ.    .CCCC.   .OOOO.   .U...U.  .N...N.  .TTTTT. .EEEEE.  .RRRR.
.D...D     .J.   .C....   .O....O.  .U...U.  .NN..N.   ..T..  .E.....  .R...R
.D...D     .J.   .C....   .O....O.  .U...U.  .N.N.N.   ..T..  .EEE...  .RRRR.
.D...D  .J..J.   .C....   .O....O.  .U...U.  .N..NN.   ..T..  .E.....  .R..R.
.DDD..   .JJ..    .CCCC.   .OOOO.   ..UUU..  .N...N.   ..T..  .EEEEE.  .R...R
#############################################################################

Usage: calCounts.sh --sample <sample> --bam <bam_file> --ref <reference> 
  --sample: sample name
  --bam: path to BAM file
  --ref: reference name (GRCh38 by default)
  --noGap: exclude gap regions in calculations (False by default)
  --filter: filter flag option for calculations (3332 by default, excluding UNMAP,SECONDARY,DUP,SUPPLEMENTARY)
  --threads: number of threads to use (10 by default)
  --accurate: use accurate background calculation (False by default, take more time but more accurate, especially for samples with uneven coverage across chromosomes)
  --fast: use fast mode (False by default, using all reads without filtering, which might be faster but less accurate)
```

###  Expected output
The output file `$outdir/$sample.$ref.tg.$filter_condition.$gap_condition.txt` contains five columns separated by tabs:
1. **sample:** This column contains identifiers or names associated with each estimation of DJ count.
2. **ref:** reference name
3. **roi:** region of interest.
4. **Estimation of DJ count based on diploid genome:** This column provides the calculated DJ count values adjusted for diploid genome context.

```bash
Sample01    GRCh38  DJ_filt 11.01608
```

## Bacgkround coverage
We calculate background coverage using autosomal chromosomes from the reference genome. This is achieved by analyzing aligned BAM files with `samtools idxstat` or `samtools view`, applying any user-specified filtering conditions, multiplying by the fragment size, and dividing by the corresponding chromosome length.

Ensure that your BAM files contain all contigs listed in `DJ_filt.bed` to enable accurate DJ coverage calculation. This checking step is also included in the script. Background coverage is computed as the median depth across autosomal chromosomes to minimize the impact of potential aneuploidy in individual chromosomes.

```
Background Coverage = median over autosomes of:

(N_c × L_fragment) / L_c

where :
- **Nc** : Number of reads mapped to the chromosome
- **Lc** : length of the chromosome
```

## Calculating the diploid DJ counts
We calculated the DJ counts for a diploid genome using the equation below.

```
norm_count = (2 × tgCount ) / (covLen × bgCov)

Where:  
- tgCount: the number of reads aligned to the target region  
- covLen: the original length of DJ on CHM13 used to normalize tgCount  
- bgCov: background autosomal coverage used for normalization  
```

## Expected results

Normal human samples typically yield around 10 copies of DJ counts, with occasional deviations to ~11 or ~9.
Robertsonian samples usually show approximately ~8 copies.

<img src="https://github.com/user-attachments/assets/9212dabb-593f-4040-bebc-494a74301fa0" width="200">

