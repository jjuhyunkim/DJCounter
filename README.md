<h1 align="center">🧬 DJCounter</h1>

<p align="center">
  <em>Estimate the copy number of ribosomal-DNA distal junctions (DJs) from sequencing data.</em>
</p>

---

## Overview

`DJCounter` estimates how many copies of the ribosomal-DNA **distal junction (DJ)** are present in a human genome from short-read sequencing data. It supports two complementary modes:

| Mode | Input | When to use |
| ---- | ----- | ----------- |
| **Mapping-based** | aligned BAM/CRAM | reads already aligned to GRCh38 / GRCh37 / CHM13|
| **K-mer based** *(reference-free)* | raw FASTQ (or BAM/CRAM) | raw reads or BAM file |

- Mapping-based: DJ copy number is derived from sequencing **coverage** in the target region, normalized to autosomal background.
- K-mer based: DJ copy number is derived from the **k-mer multiplicity** of a curated DJ-specific 31-mer set, normalized to the 2-copy peak in the read k-mer histogram.

Typical human samples yield ~10 DJ copies and Robertsonian samples typically show ~8.

## Install

DJCounter is a shell/data pipeline that wraps external bioinformatics tools. Any of the following work on Linux and macOS (arm64/x86_64):

### Conda / mamba (recommended, lightweight)

1. Use the environment file shipped with this repo

```bash
git clone https://github.com/marbl/DJCounter && cd DJCounter
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate djcounter
# unpack the bundled DJ k-mer database once
tar -xzf resources/DJtarget.meryl.tar.gz -C resources
# put the wrappers on $PATH (matches the bioconda install)
export PATH="$PWD/bin:$PATH"
```

2. When a bioconda version becomes available

```bash
mamba install -c bioconda -c conda-forge djcounter
```

### Docker / Singularity (via BioContainers)

```bash
# Docker
docker run --rm -v "$PWD":/data quay.io/biocontainers/djcounter:1.1--hdfd78af_0 \
    djcounter-kmer -sample Sample01 -input /data/reads.fq.gz

# Singularity / Apptainer
singularity exec \
    https://depot.galaxyproject.org/singularity/djcounter:1.1--hdfd78af_0 \
    djcounter-kmer -sample Sample01 -input reads.fq.gz
```

* Replace `<build>` with the tag shown on [quay.io/biocontainers/djcounter](https://quay.io/repository/biocontainers/djcounter?tab=tags).

### Manual install

Assuming the rest of the [Dependencies](#dependencies) below are available, download and set environment paths as needed. Make sure `samtools`, `java`, and `bc` are on `$PATH`.

```bash
# DJCounter - nothing to install
git clone https://github.com/marbl/DJCounter && cd DJCounter
tar -xzf resources/DJtarget.meryl.tar.gz -C resources # for k-mer counting DJs
export PATH="$PWD/bin:$PATH"
cd ../

# Merqury - nothing to install
git clone https://github.com/marbl/merqury
export MERQURY=$PWD/merqury

# Meryl - start with runnable binaries
wget https://github.com/marbl/meryl/releases/download/v1.4.2/meryl-1.4.2.Linux-amd64.tar.xz # replace to your system-compatible binary
tar -xJf meryl-1.4.2.Linux-amd64.tar.xz
export PATH=$PWD/meryl-1.4.2/bin:$PATH
```


## Quick start

Here we are testing with a sub-sampled 5x HG002 bam file.

Note the results are slightly less accurate due to its low coverage.

### 0. Download a test bam

Download the following files:
```
https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/DJCounter_2026/test/hg002_5x.bam
https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/DJCounter_2026/test/hg002_5x.bam.bai
```

### 1. Mapping-based

Suitable when your BAM/CRAM is aligned to one of the supported references (see [References](#supported-references)).

```bash
djcounter-map \
    --sample  hg002_5x_map \
    --bam     hg002_5x.bam \
    --ref     GRCh38 \
    --threads 10
```

Output: `hg002_5x_map.GRCh38.filter_3332.incGap.fast-accurate.byWGS.readCount.tg.txt`

```
sample           ref      roi  tgCount   bgCov  DJ_count
hg002_5x_map  GRCh38  DJ_filt    61321  .03237  10.94933
```
* ref: reference
* roi: region of interest, bed file used for collecting coverage under [roi](roi)
* tgCount: target read count
* bgCov: background coverage
* DJ_count: estimated DJ counts

📘 Details: [scripts/mapping_based.md](scripts/mapping_based.md)

### 2. K-mer based (reference-free)

Suitable for BAM/CRAMs mapped to a non-supported reference or raw FASTQ files.

```bash
djcounter-kmer -sample hg002_5x_kmer -input hg002_5x.bam -tmp .

# also works with fq.gz files, i.e.:
djcounter-kmer -sample hg002_5x_kmer -input hg002_5x.R1.fq.gz,hg002_5x.R2.fq.gz -tmp .
```

Output: DJcounts/hg002_5x_kmer_DJ_count.txt

```
sample          tgMult  bgMult  DJ_count
hg002_5x_kmer       22       4        11
```
* tgMult: target k-mer multiplicity, median multiplicity of the target DJ k-mers
* bgMult: background k-mer multiplicity, estimated for 2-copy peak
* DJ_count: estimated DJ counts

Plot the distribution across many samples:

```bash
cat DJcounts/*_DJ_count.txt > DJ_counts.txt
Rscript scripts/plot_dist.R
```

📘 Details: [scripts/kmer_based.md](scripts/kmer_based.md)

## How it works

### Mapping-based

```
DJ_count = (2 × tgCount) / (covLen × bgCov)

  tgCount : reads aligned to the DJ target regions
  covLen  : DJ length on CHM13 used to normalize tgCount
  bgCov   : background autosomal coverage
```

### K-mer based

1. Count all 31-mers in the input (`meryl count k=31`).
2. Intersect with the curated `DJtarget.meryl` set (52,227 distinct k-mers; 26,140,589 occurrences) and read the median frequency from its histogram.
3. Use Merqury's `kmerHistToPloidyDepth.jar` to estimate the 2-copy peak from the read k-mer histogram.

```
DJ_count = (2 × tgMult) / bgMult

  tgMult.  : target k-mer multiplicity, median multiplicity of the target DJ k-mers
  bgMult.  : background k-mer multiplicity, estimated for 2-copy peak
```

## Supported references for mapping-based counting

| Build       | Required sequences | Notes |
| ----------- | ------------------ | ----- |
| **GRCh38 / hg38** | `chr21`, `chrUn_GL000220v1`, `chr17_GL000205v2_random`, `chr22_KI270733v1_random`, `chrUn_GL000195v1` | [Broad ver.](https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large/) (UK Biobank) or [1KGP NYGC ver.](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa) |
| **T2T-CHM13 / hs1** | `chr13`, `chr14`, `chr15`, `chr21`, `chr22` | [Masked T2T-CHM13](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/), using all-but-one target masked, sample-matched [XX (chrY masked)](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/chm13v2.0_masked_DJ_5S_rDNA_PHR_noY_wi_rCRS.fa) or [XY (PAR masked)](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa) sex chromosome compliment |

Verify your BAM contains the required contigs:

```bash
samtools view -H sample.bam | grep chr17_GL000205v2_random
```

## Repository layout

```
DJCounter/
├── bin/             # Thin wrappers (djcounter-kmer, djcounter-map)
├── scripts/         # Pipeline scripts and per-mode docs
│   ├── calCounts.sh                 # mapping based
│   ├── kmer_based_dj_counting.sh    # kmer based
│   ├── mapping_based.md
│   └── kmer_based.md
├── resources/       # Pre-built DJ k-mer database & references
│   └── DJtarget.meryl.tar.gz
├── roi/             # Target BED files
│   ├── GRCh38/
│   ├── hg19/
│   └── CHM13/
├── recipes/         # Bioconda recipe (also produces Docker/Singularity via BioContainers)
│   └── djcounter/
├── environment.yml  # Conda/mamba env for one-command install
└── paper/           # jupyter notebook for generating plots
```

## Dependencies

- [`samtools`](https://www.htslib.org/) ≥ 1.21 — mapping-based mode
- [`meryl`](https://github.com/marbl/meryl) ≥ 1.4.2 — k-mer mode
- [`merqury`](https://github.com/marbl/merqury) — only `eval/kmerHistToPloidyDepth.jar`; set `$MERQURY` to the clone path
- Java runtime (for the Merqury jar)
- `R` (optional, for plotting)

## Changelog

| Version | Date | Changes |
| ------- | ---- | ------- |
| **v1.1** | 2026-07-22 | Update meryl ver to accept BAM/CRAM for k-mer mode |
| v1.0 | 2026-03-08 | Finalized hg38 and k-mer modes |
| v0.2.2  | 2025-11-26 | Added BED file for ROI on hg19 |
| v0.2.1  | 2024-07-29 | Output background and fragment size; fixed background command |
| v0.2    | 2024-07-25 | `samtools idxstats` → `samtools coverage` for background; removed temp files |
| v0.1    | 2024-07-17 | First commit |

## Citation

Please use this [paper](https://doi.org/10.64898/2026.03.08.710242) to cite DJCounter:

Rhie, A., Kim, J., Rodriguez-Algarra, F. et al. Biobank-scale genotyping of Robertsonian translocations reveals hidden structural variation on the human acrocentric chromosomes. bioRxiv (2026). https://doi.org/10.64898/2026.03.08.710242
