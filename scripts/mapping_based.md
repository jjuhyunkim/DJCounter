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
Usage: DJCounter/scripts/calCounts.sh --sample <sample> --ref <reference> --bam <bam_file>
  --sample: sample name
  --ref: reference name (hg19 or GRCh38)
  --bam: path to BAM file
  --noGap: exclude gap regions in calculations (False by default)
  --filter: filter option for calculations (3332 by default)
  --threads: number of threads to use (10 by default)
  --targetlist: list of target regions (DJ_filt by default)
```

###  Expected output
The output file `$outdir/$sample.$ref.tg.$filter_condition.$gap_condition.txt` contains five columns separated by tabs:
1. **sample:** This column contains identifiers or names associated with each estimation of DJ count.
2. **ref:** reference name
3. **roi:** region of interest.
4. **background coverage:** background coverage calculated from autosome
5. **Estimation of DJ count based on diploid genome:** This column provides the calculated DJ count values adjusted for diploid genome context.
```bash
Sample01    GRCh38  DJ_filt 33.9158 11.01608
```

## Details
```bash
# This is an input example
sample=Sample01
threads=10

bam="/data/01.broad_hg38/$sample/$sample.dedup.bam" # BAM or CRAM file
outdir="/data/01.broad_hg38/$sample" # The output directory
prefix="$sample" # Prefix for output files
bed="/roi/DJ_filt.bed" # BED file for DJ counting
```

## Bacgkround coverage
We calculate background coverage using autosomal chromosomes from the reference genome. This is achieved by analyzing aligned BAM files with `samtools idxstat` or `samtools view`, applying any user-specified filtering conditions, multiplying by the fragment size, and dividing by the corresponding chromosome length.

Ensure that your BAM files contain all contigs listed in `DJ_filt.bed` to enable accurate DJ coverage calculation. This checking step is also included in the script. Background coverage is computed as the median depth across autosomal chromosomes to minimize the impact of potential aneuploidy in individual chromosomes.

$$
\text{Background Coverage}
=
\operatorname{median}_{c \in \text{autosomes}}
\left(
\frac{N_c \times L_{\text{fragment}}}{L_c}
\right)
$$

where :
- **Nc** : Number of reads mapped to the chromosome
- **Lfragment** : length of fragments
- **Lc** : length of the chromosome


## Calculating the diploid DJ counts
We calculated the DJ counts for a diploid genome using the equation below.

$$
\text{norm\_count}
=
\frac{2 \times \text{tgCount} \times \text{fragmentSize}}
{\text{covLen} \times \text{bgCov}}
$$

Where:  
- **tgCount**: the number of reads aligned to the target region  
- **covLen**: the original length of DJ on CHM13 used to normalize tgCount  
- **bgCov**: background autosomal coverage used for normalization  
- **fragmentSize**: average fragment size of the sequencing library

## Expected results

Normal human samples typically yield around 10 copies of DJ counts, with occasional deviations to ~11 or ~9.
Robertsonian samples usually show approximately ~8 copies.

<img src="https://github.com/user-attachments/assets/9212dabb-593f-4040-bebc-494a74301fa0" width="200">

