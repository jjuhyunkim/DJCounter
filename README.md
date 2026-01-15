# DJ Counter

Count the number of distal junctions (DJ) of the ribosomal DNA (rDNA)

## Purpose
This repository provides tools for estimating the copy number of the DJs in a genome from sequencing reads.

## Principle
Copy number of DJs are determined by 1) the sequencing coverage in mapped reads or 2) k-mer multiplicity in raw reads.

## Main workflow
Copy number of the DJ can be estimated with the following approaches:
1. [Mapping based approach](/scripts/mapping_based.md)
2. [K-mer based, reference-free approach](/scripts/kmer_based.md)

### Mapping based estimates

The mapping based approach is recommended when the reads are already aligned to one of the following references.

1. GRCh38/hg38 [Homo_sapiens_assembly38.fasta.gz](https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/) (1000 Genomes Project Broad ver. Suitable for UKBioBank)

   Requires `chr21`, `chr17_GL000205v2_random` and `chrUn_GL000195v1`.

2. GRCh37/hg19 [human_g1k_v37.fasta.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/) (1000 Genomes Project ver)

   Requires `chr7_gl000195_random` and `chr17_gl000205_random`.

3. T2T-CHM13/hs1 (Will be updated soon)

For hg38 or hg19, check with `samtools` to confirm the sequence exists in the header; such as using `samtools view -H in.bam | grep chr17_GL000205v2_random`.

Read coverage is assessed on the mapped BAM file for the target DJ region and compared against the background coverage collected from autosomes.

### K-mer based approach

The k-mer based approach is reference-free.

This approach is recommended when reads are aligned to hg38 or hg19 _without any decoy sequences_ or are in its raw FASTQ form. A collected set of target k-mers are pre-built to query the k-mer multiplicity of the DJ and is compared against the single / 2-copy copy number estimates inferred from the k-mer multiplicity histogram.

## Change logs
<details>
<summary>v0.1(2024-07-17)</summary>
* first commit
</details>

<details>
<summary>v0.2(2024-07-25)</summary>
* Changing the background coverage estimation methods from samtools idxstats to samtools coverage.<br />
* Removing the step of saving temporary files; instead, we assign everything to variables.
</details>

<details>
<summary>v0.2.1(2024-07-29)</summary>
* Add background and fragment size to the output file. <br />
* Fix the command line used for calculating the background to ensure it works correctly.
</details>

<details>
<summary>v0.2.2(2025-11-26)</summary>
* Add BED file for roi on hg19 <br />
</details>




