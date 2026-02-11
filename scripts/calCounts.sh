#!/bin/bash
ml samtools
ml bedtools

set -e 

# Parse named arguments
sample=""
refname=""
bam=""
fasta=true
noGap=false
accurate=false
targetlist="DJ_filt"
filter_flag="3332"
threads=10

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            sample="$2"
            shift 2
            ;;
        --ref)
            refname="$2"
            shift 2
            ;;
        --bam)
            bam="$2"
            shift 2
            ;;
        --noGap)
            noGap=true
            shift 1
            ;;
        --accurate)
            accurate=true
            shift 1
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --targetlist)
            targetlist="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --sample <sample> --ref <reference> --bam <bam_file> [--noGap --filter <filter> --threads <threads> --targetlist <targetlist>]"
            exit 1
            ;;
    esac
done
 
if [ -z "$sample" ] || [ -z "$refname" ] || [ -z "$bam" ]; then
    echo "Usage: $0 --sample <sample> --ref <reference> --bam <bam_file>"
    echo "  --sample: sample name"
    echo "  --ref: reference name (hg19 or GRCh38)"
    echo "  --bam: path to BAM file"
    # echo "  --noGap: exclude gap regions in calculations (False by default)"
    # echo "  --filter: filter option for calculations (3332 by default)"
    echo "  --threads: number of threads to use (10 by default)"
    # echo "  --targetlist: list of target regions (DJ_filt by default)"
    exit 1
fi
prefix=${sample}.${refname}

echo -e "calculating median background coverage"
if [ "$noGap" = true ]; then
    gap_info=noGap
else 
    gap_info=incGap
fi

if [ "$filter_flag" != "" ]; then
    filter_info=filter_$filter_flag
else
    filter_info="noFilter"
fi

## check the parameters
echo -e "Sample: $sample"
echo -e "Reference: $refname"
echo -e "BAM file: $bam"
echo -e "No Gap: $noGap"
echo -e "Filter: $filter_flag"
echo -e "Threads: $threads"
echo -e "Prefix: $prefix"
echo -e "Target List: $targetlist"

roi_dir=$(dirname $(realpath $0))/../roi
echo "ROI directory: $roi_dir"

if [ "$refname" == "hg19" ]; then
    # ref=$roi_dir/hg19/hg19.p13.plusMT.full_analysis_set.fa.gz
    targetDir=$roi_dir/hg19/
    if [ "$noGap" == true ]; then
        autosomeBed=$roi_dir/hg19/autosome.nogap.bed
        backgroundLen=$roi_dir/hg19/autosome.nogap.len
    else
        autosomeBed=$roi_dir/hg19/autosome.bed
        backgroundLen=$roi_dir/hg19/autosome.len
    fi
elif [ "$refname" == "GRCh38" ]; then 
    # ref=$roi_dir/GRCh38/Homo_sapiens_assembly38.fasta
    targetDir=$roi_dir/GRCh38/
    if [ "$noGap" == true ]; then
        autosomeBed=$roi_dir/GRCh38/autosome.nogap.bed
        backgroundLen=$roi_dir/GRCh38/autosome.nogap.len
    else
        autosomeBed=$roi_dir/GRCh38/autosome.bed
        backgroundLen=$roi_dir/GRCh38/autosome.len
    fi
else
    echo "Unknown reference name: $refname"
    exit 1
fi
echo "Using reference: $refname"

# Check whether all required contigs are present in the BAM header.
for target in $targetlist; do
    if [ ! -f "$roi_dir/$refname/$target.bed" ]; then
        echo -e "$target.bed is not detected in the $roi_dir/$refname"
        exit 1
    else
        cut -f 1 "$roi_dir/$refname/$target.bed" | sort | uniq > tg.contigs.tmp
        contigs_bed=$(wc -w < tg.contigs.tmp)
        echo -e "Number of contigs in the target bed ($target.bed): $contigs_bed"
        contigs_bam=$(samtools view -H "$bam" | grep -wf tg.contigs.tmp | wc -l)
        echo -e "Number of contigs in the BAM header: ${contigs_bam}"
        rm tg.contigs.tmp
        if [ "$contigs_bam" -ne "$contigs_bed" ]; then
            echo "Some contigs in $target.bed are not present in the BAM header."
            exit 1
        else 
            echo "All contigs in $target.bed are present in the BAM header."
        fi
    fi
done

# FRAGEMENT SIZE
# fragmentSize=$(samtools view "$bam" | head -10000 | awk '{print length($10)}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
# echo $fragmentSize > fragmentsize.txt
# echo "fragment size : $fragmentSize"

if [ ! -f ${bam}.bai ]; then
	samtools index -@ $threads $bam
fi
if [ -f "${prefix}.$filter_info.$gap_info.readCount.txt" ] && [ $(wc -l < "${prefix}.$filter_info.$gap_info.readCount.txt") -lt 22 ]; then
    echo "Removing incorrect filter.readCount.txt file"
    rm "${prefix}.$filter_info.$gap_info.readCount.txt"
fi

# BACKGROUND PROCESSING
if [ "$filter_flag" == "" ] && [ ! -f "${prefix}.noFilter.${gap_info}.readCount.txt" ]; then
    samtools idxstats -@ $threads $bam > ${prefix}.noFilter.idxstats
    cut -f 1,3 ${prefix}.noFilter.idxstats > ${prefix}.noFilter.${gap_info}.readCount.txt

elif [ "$filter_flag" != "" ] && [ ! -f "${prefix}.$filter_info.$gap_info.readCount.txt" ]; then
    if [ $accurate == false ]; then
        echo -e "Calculating background read count without per-chromosome processing"
        autosomeReadCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L "$autosomeBed" $bam)
        bgCov=$(echo "scale=5; $autosomeReadCount / 2875001522" | bc)
    else
        echo -e "Calculating background read count with per-chromosome processing"
        echo -e "This might take longer time but will be more accurate, especially for samples with uneven coverage across chromosomes."
        for chr in {1..22}; do
            echo "Processing chromosome $chr"
            sed -n "${chr}p" "$autosomeBed" > tmp.bed
            readCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L tmp.bed "$bam")
            echo -e "chr${chr}\t$readCount" >> "${prefix}.$filter_info.$gap_info.readCount.txt"
        done
        rm tmp.bed
        bgCov=$(join -t $'\t' -1 1 -2 1  "${prefix}.$filter_info.$gap_info.readCount.txt" "$backgroundLen" | awk '{print $2/$3}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
    fi
fi


#readCount=$(samtools view -@ $threads -F0x1504 -c -L "$autosomeBed"  "$bam")
#echo -e "$readCount" > autosome.readCount.txt
#echo "scale=6; $fragmentSize * $readCount / 2745187818" | bc >> background.allAutosome_noDup.txt        


# bgCov=$(join -t $'\t' -1 1 -2 1  "${prefix}.$filter_info.$gap_info.readCount.txt" "$backgroundLen" | awk '{print $2/$3}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
# echo $bgCov > background.$filter_info.$gap_info.txt
echo "Background coverage (bgCov): $bgCov"
# CALCULATE DEPTH AND COUNTS FOR TARGET REGIONS
rm -rf "${prefix}.tg.$filter_info.$gap_info.txt"

echo -e "Reading background coverage: background.$filter_info.$gap_info.txt"
# bgCov=$(cat background.$filter_info.$gap_info.txt)

# targetlist=$(ls "${targetDir}"/*.bed 2>/dev/null | xargs -n 1 basename | sed -e 's/\.bed$//' | tr '\n' ' ')
for target in $targetlist; do
    bed="${targetDir}/$target.bed"
    len="${targetDir}/$target.len"
    if [ ! -f "$bed" ]; then
        echo "Warning: BED file $bed does not exist, skipping target $target."
        continue
    fi
   
    if [ ! -s "$bed" ]; then
        echo "Warning: BED file $bed is missing or empty, skipping."
        continue
    fi
    covLen=$(awk '{sum += $3 - $2} END {print sum}' "$bed")
    if [ -z "$covLen" ] || [ "$covLen" -eq 0 ]; then
        echo "Warning: covLen is zero for $bed, skipping."
        continue
    fi
    echo "Processing target: $target"
    # if [ ! -f "${target}.depth" ]; then
    #    echo "Calculating depth for $target"
    #    samtools depth -@ "$threads" -b "$bed" "$bam" > "${target}.depth"
    #else
    #    echo "Depth file already exists for $target, skipping calculation."
    #fi
    #tgCount_depth=$(awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }' "${target}.depth")
    # echo "$tgCount_depth"
    # if [ ! -f "${target}.$filter_info.$gap_info.counts" ]; then
    if [ "$filter_flag" != "" ]; then
        echo "Calculating counts for $target with filter: ${filter_flag})"
        tgCount_count=$(samtools view -@ "$threads" -F ${filter_flag} -c -L "$bed" "$bam")
    else
        echo "Calculating counts for $target including all alignments - no filter" 
        tgCount_count=$(samtools view -@ "$threads" -c -L "$bed" "$bam") 
    fi

    # echo -e "Reading counts for $target : ${target}.$filter_info.$gap_info.counts"
    # tgCount_count=$(cat "${target}.$filter_info.$gap_info.counts")
    # echo "$tgCount_count"
    # Estimate the DJ count
    if [ -z "$bgCov" ] || [ "$bgCov" = "0" ]; then
        #echo "Warning: bgCov is zero, skipping normalization for $target."
        norm_depth="NA"
        norm_count="NA"
    elif [ ! -f "$len" ]; then
        #echo "Warning: Length file $len does not exist, skipping target $target."
        norm_count="NA"
        norm_depth="NA"
    else
        covLen=$(cat "$len")
        # norm_depth=$(echo "scale=5; 2 * $tgCount_depth / $covLen / $bgCov" | bc)
        # norm_count=$(echo "scale=5; 2 * $tgCount_count * $fragmentSize / $covLen / $bgCov" | bc)
        norm_count=$(echo "scale=5; 2 * $tgCount_count / $covLen / $bgCov" | bc)
    fi
    # echo -e "$prefix\t${target}\t${norm_depth}\t${norm_count}" >> "$prefix.tg.$filter_info.$gap_info.txt"
    echo -e "$sample\t$refname\t${target}\t${norm_count}" >> "$prefix.tg.$filter_info.$gap_info.txt"
done