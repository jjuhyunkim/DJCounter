#!/bin/bash
# if module could be used, load the required modules
if command -v module &> /dev/null; then
    ml samtools
    ml bedtools
fi

set -e 

echo "#############################################################################"
echo ".DDD..   .JJJ.    .CCCC.   .OOOO.   .U...U.  .N...N.  .TTTTT. .EEEEE.  .RRRR."
echo ".D...D     .J.   .C....   .O....O.  .U...U.  .NN..N.   ..T..  .E.....  .R...R"
echo ".D...D     .J.   .C....   .O....O.  .U...U.  .N.N.N.   ..T..  .EEE...  .RRRR."
echo ".D...D  .J..J.   .C....   .O....O.  .U...U.  .N..NN.   ..T..  .E.....  .R..R."
echo ".DDD..   .JJ..    .CCCC.   .OOOO.   ..UUU..  .N...N.   ..T..  .EEEEE.  .R...R"
echo "#############################################################################"

# Parse named arguments
sample=""
bam=""
refname="GRCh38"
noGap=false
fast=false
accurate=false
targetlist="DJ_filt"
filter_flag="3332"
threads=10
clean=false

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
        --fast)
            fast=true
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
        -clean)
            clean=true
            shift 1
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --sample <sample> --bam <bam_file> --ref <reference> [--noGap --filter <filter> --threads <threads> --targetlist <targetlist>]"
            exit 1
            ;;
    esac
done
 
if [ -z "$sample" ] || [ -z "$refname" ] || [ -z "$bam" ]; then
    echo "Usage: $0 --sample <sample> --bam <bam_file> --ref <reference> "
    echo "  --sample: sample name"
    echo "  --bam: path to BAM file"
    echo "  --ref: reference name (GRCh38 by default)"
    echo "  --noGap: exclude gap regions in calculations (False by default)"
    echo "  --filter: filter option for calculations (3332 by default)"
    echo "  --threads: number of threads to use (10 by default)"
    echo "  --accurate: use accurate background calculation (False by default, take more time but more accurate, especially for samples with uneven coverage across chromosomes)"
    echo "  --fast: use fast mode (False by default, using all reads without filtering, which might be faster but less accurate)"
    # echo "  --targetlist: list of target regions (DJ_filt by default)"
    exit 1
fi
prefix=${sample}.${refname}

if [ "$fast" = true ] && [ "$accurate" = true ]; then
    echo "Fast mode and accurate mode cannot be enabled at the same time. Please choose one of them."
    exit 1
fi

if [ "$fast" = true ]; then
    filter_flag=""
    echo "Fast mode enabled: using all reads without filtering."
fi

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
echo ""
echo "Parameters:"
echo -e "Sample: $sample"
echo -e "Reference: $refname"
echo -e "BAM file: $bam"
echo -e "Filter: $filter_flag"
echo -e "Threads: $threads"
echo -e "Prefix: $prefix"
echo -e "Target List: $targetlist"
echo -e "Accurate mode: $accurate"
echo -e "Fast mode: $fast"
echo -e "Gap info: $gap_info"
echo -e "Clean mode: $clean"
echo ""

roi_dir=$(dirname $(realpath $0))/../roi
# echo "ROI directory: $roi_dir"

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
autosomeLen=$(awk '{sum += $3 - $2} END {print sum}' $autosomeBed)

# Check whether all required contigs are present in the BAM header.
for target in $targetlist; do
    if [ ! -f "$roi_dir/$refname/$target.bed" ]; then
        echo -e "$target.bed is not detected in the $roi_dir/$refname"
        exit 1
    else
        cut -f 1 "$roi_dir/$refname/$target.bed" | sort | uniq > tg.contigs.tmp
        contigs_bed=$(wc -w < tg.contigs.tmp)
        # echo -e "Number of contigs in the target bed ($target.bed): $contigs_bed"
        contigs_bam=$(samtools view -H "$bam" | grep -wf tg.contigs.tmp | wc -l)
        # echo -e "Number of contigs in the BAM header: ${contigs_bam}"
        rm tg.contigs.tmp
        if [ "$contigs_bam" -ne "$contigs_bed" ]; then
            echo "Some contigs in $target.bed are not present in the BAM header."
            exit 1
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

if [ "$accurate" = true ]; then
    accurate_info="byChr"
else 
    accurate_info="wholeGenome"
fi

# BACKGROUND PROCESSING
if [ -f "${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt" ]; then
    if [ ! -s "${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt" ]; then
        echo "Removing incorrect bg read count file"
        rm "${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt"
    fi
fi

if [ ! -f "${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt" ] ; then
    if [ "$fast" = true ]; then
        samtools idxstats -@ $threads $bam > ${prefix}.noFilter.idxstats
        cut -f 1,3 ${prefix}.noFilter.idxstats > ${prefix}.noFilter.${gap_info}.readCount.txt
        autosomeReadCount=$(awk '{sum += $2} END {print sum}' "${prefix}.noFilter.${gap_info}.readCount.txt")
        rm ${prefix}.noFilter.idxstats ${prefix}.noFilter.${gap_info}.readCount.txt
        bgCov=$(echo "scale=5; $autosomeReadCount / $autosomeLen" | bc)

    elif [ "$filter_flag" != "" ] && [ "$accurate" = false ];then 
        if [ ! -f "${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt" ]; then
            echo -e "Calculating background read count without per-chromosome processing"
            autosomeReadCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L "$autosomeBed" $bam)
            bgCov=$(echo "scale=5; $autosomeReadCount / $autosomeLen" | bc)
        fi
        
    elif [ "$filter_flag" != "" ] && [ "$accurate" = true ]; then
        if [ -f "${prefix}.$filter_info.$gap_info.readCount.txt" ] && [ $(wc -l < "${prefix}.$filter_info.$gap_info.readCount.txt") -lt 22 ]; then
            echo "Removing incorrect filter.readCount.txt file"
            rm "${prefix}.$filter_info.$gap_info.readCount.txt"
        fi
        
        if [ ! -f "${prefix}.$filter_info.$gap_info.readCount.txt" ]; then
        echo -e "Calculating background read count with per-chromosome processing"
        echo -e "This might take longer time but will be more accurate, especially for samples with uneven coverage across chromosomes."
            
            for chr in {1..22}; do
                echo "Processing chromosome $chr"
                sed -n "${chr}p" "$autosomeBed" > tmp.bed
                readCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L tmp.bed "$bam")
                echo -e "chr${chr}\t$readCount" >> "${prefix}.$filter_info.$gap_info.readCount.txt"
            done
            rm tmp.bed
        fi
        bgCov=$(join -t $'\t' -1 1 -2 1  "${prefix}.$filter_info.$gap_info.readCount.txt" "$backgroundLen" | awk '{print $2/$3}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
    fi
    echo $bgCov > "${prefix}.$filter_info.${gap_info}.${accurate_info}.readCount.bg.txt"
else
    # echo -e "Reading background coverage from file: ${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt"
    bgCov=$(cat "${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt")
fi

# echo "Background coverage (bgCov): $bgCov"

rm -rf "${prefix}.tg.$filter_info.$gap_info.txt"

# echo -e "Reading background coverage: background.$filter_info.$gap_info.txt"
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
    # echo "Processing target: $target"
    # if [ ! -f "${target}.depth" ]; then
    #    echo "Calculating depth for $target"
    #    samtools depth -@ "$threads" -b "$bed" "$bam" > "${target}.depth"
    #else
    #    echo "Depth file already exists for $target, skipping calculation."
    #fi
    #tgCount_depth=$(awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }' "${target}.depth")
    # echo "$tgCount_depth"
    # if [ ! -f "${target}.$filter_info.$gap_info.counts" ]; then

    if [ -f "${target}.$filter_info.$gap_info.counts" ]; then
        # remove if the filesize is 1 or less, which indicates an error in the previous run
        if [ ! -s "${target}.$filter_info.$gap_info.counts" ]; then
            echo "Warning: Count file ${target}.$filter_info.$gap_info.counts is empty, removing it and recalculating."
            rm "${target}.$filter_info.$gap_info.counts"
        else
            echo "Count file already exists for $target, skipping calculation."
            tgCount_count=$(cat "${target}.$filter_info.$gap_info.counts")
        fi
    else
        if [ "$filter_flag" != "" ]; then
            echo "Calculating counts for $target with filter: ${filter_flag}"
            tgCount_count=$(samtools view -@ "$threads" -F ${filter_flag} -c -L "$bed" "$bam")
        else
            echo "Calculating counts for $target including all alignments - no filter" 
            tgCount_count=$(samtools view -@ "$threads" -c -L "$bed" "$bam") 
        fi
        echo "$tgCount_count" > "${target}.$filter_info.$gap_info.counts"
    fi

    if [ "$clean" = true ]; then
        echo "Clean mode enabled: removing intermediate files for $target"
        rm -f "${target}.$filter_info.$gap_info.counts"
    fi

    # echo -e "Reading counts for $target : ${target}.$filter_info.$gap_info.counts"
    # tgCount_count=$(cat "${target}.$filter_info.$gap_info.counts")
    # echo "$tgCount_count"
    # Estimate the DJ count
    if [ -z "$bgCov" ] || [ "$bgCov" = "0" ]; then
        #echo "Warning: bgCov is zero, skipping normalization for $target."
        # norm_depth="NA"
        norm_count="NA"
    elif [ ! -f "$len" ]; then
        #echo "Warning: Length file $len does not exist, skipping target $target."
        # norm_count="NA"
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