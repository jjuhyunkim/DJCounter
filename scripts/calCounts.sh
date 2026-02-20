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
mode="fast-accurate"
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
        --mode)
            mode="$2"
            shift 2
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
            echo "Usage: $0 --sample <sample> --bam <bam_file> --ref <reference> --mode <mode> [--noGap --filter <filter> --threads <threads> --targetlist <targetlist>]"
            exit 1
            ;;
    esac
done
 
if [ -z "$sample" ] || [ -z "$refname" ] || [ -z "$bam" ]; then
    echo "Usage: $0 --sample <sample> --bam <bam_file> --ref <reference> "
    echo "  --sample: sample name"
    echo "  --bam: path to BAM file"
    echo "  --ref: reference name (GRCh38 by default, CHM13 and hg19 are also supported)"
    echo "  --noGap: exclude gap regions in calculations (False by default)"        
    echo "  --threads: number of threads to use (10 by default)"
    echo "  --mode: calculation mode (fast-accurate by default, options: fast, fast-accurate, high-accurate, high-accurate-coverage)"
    # echo "  --targetlist: list of target regions (DJ_filt by default)"
    exit 1
fi
prefix=${sample}.${refname}

# check if the mode is valid
if [ "$mode" != "fast" ] && [ "$mode" != "fast-accurate" ]  && [ "$mode" != "high-accurate" ] && [ "$mode" != "high-accurate-coverage" ] ; then
    echo "Invalid mode: $mode. Valid options are: fast, fast-accurate, high-accurate, high-accurate-coverage"
    exit 1
fi

if [ "$mode" = "fast" ]; then
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
echo -e "Mode: $mode"
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
elif [ "$refname" == "CHM13" ]; then
    targetDir=$roi_dir/CHM13/
    if [ "$noGap" == true ]; then
        autosomeBed=$roi_dir/CHM13/autosome.nogap.bed
        backgroundLen=$roi_dir/CHM13/autosome.nogap.len
    else
        autosomeBed=$roi_dir/CHM13/autosome.bed
        backgroundLen=$roi_dir/CHM13/autosome.len
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

if [ ! -f ${bam}.bai ]; then
	samtools index -@ $threads $bam
fi

if [ "$mode" = "high-accurate" ]; then
    accurate_info="byChr"
else 
    accurate_info="wholeGenome"
fi

# BACKGROUND PROCESSING
if [ -f "${prefix}.$filter_info.$gap_info.$mode.readCount.bg.txt" ]; then
    if [ ! -s "${prefix}.$filter_info.$gap_info.$mode.readCount.bg.txt" ]; then
        echo "Removing incorrect bg read count file"
        rm "${prefix}.$filter_info.$gap_info.$mode.readCount.bg.txt"
    fi
fi

if [ ! -f "${prefix}.$filter_info.$gap_info.${mode}.readCount.bg.txt" ] ; then
    if [ "$mode" = "fast" ]; then
        samtools idxstats -@ $threads $bam > ${prefix}.noFilter.idxstats
        awk '{print $1,$3}' ${prefix}.noFilter.idxstats > ${prefix}.noFilter.${gap_info}.readCount.txt
        autosomeReadCount=$(awk '{sum += $2} END {print sum}' "${prefix}.noFilter.${gap_info}.readCount.txt")
        rm ${prefix}.noFilter.idxstats ${prefix}.noFilter.${gap_info}.readCount.txt
        bgCov=$(echo "scale=5; $autosomeReadCount / $autosomeLen" | bc)
        echo $bgCov > "${prefix}.$filter_info.${gap_info}.${mode}.readCount.bg.txt"

    elif [ "$filter_flag" != "" ] && [ "$mode" = "fast-accurate" ];then 
        if [ ! -f "${prefix}.$filter_info.$gap_info.${mode}.readCount.bg.txt" ]; then
            echo -e "Calculating background read count without per-chromosome processing"
            autosomeReadCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L "$autosomeBed" $bam)
            bgCov=$(echo "scale=5; $autosomeReadCount / $autosomeLen" | bc)
            echo $bgCov > "${prefix}.$filter_info.${gap_info}.${mode}.readCount.bg.txt"
        fi
        
    elif [ "$filter_flag" != "" ] && [ "$mode" = "high-accurate-read" ]; then
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
        echo $bgCov > "${prefix}.$filter_info.${gap_info}.${mode}.readCount.bg.txt"
    
    elif [ "$filter_flag" != "" ] && [ "$mode" = "high-accurate-coverage" ]; then
        echo "MODE : $mode"
        if [ ! -f $prefix.cov ]; then
            samtools coverage $bam > $prefix.cov
        fi
        bgCov=$(grep -E '^chr([1-9]|1[0-9]|2[0-2])([[:space:]]|$)' $prefix.cov | cut -f 7 | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}')
    fi
    
else
    # echo -e "Reading background coverage from file: ${prefix}.$filter_info.$gap_info.$accurate_info.readCount.bg.txt"
    bgCov=$(cat "${prefix}.$filter_info.$gap_info.$mode.readCount.bg.txt")
fi

rm -rf "${prefix}.tg.$filter_info.$mode.$gap_info.txt"

for target in $targetlist; do
    bed="${targetDir}/$target.bed"
    len="${targetDir}/$target.len"
    
    if [ ! -f "$bed" ] || [ ! -s "$bed" ]; then
        echo "Warning: BED file $bed does not exist or is empty, skipping target $target."
        continue
    fi

    covLen=$(awk '{sum += $3 - $2} END {print sum}' "$bed")
    if [ -z "$covLen" ] || [ "$covLen" -eq 0 ]; then
        echo "Warning: covLen is zero for $bed, skipping."
        continue
    fi

    # echo "Processing target: $target"
    if [ "$mode" = "high-accurate-coverage" ]; then
        if [ ! -f "${target}.depth" ]; then
            echo "Calculating depth for $target"
            samtools depth -@ "$threads" -b "$bed" "$bam" > "${target}.depth"
        elif [ -f "${target}.depth" ]; then
            echo "Depth file already exists for $target, skipping calculation."
        fi
        tgCount_count=$(awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }' "${target}.depth")
    else 
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
    fi

    if [ "$clean" = true ]; then
        echo "Clean mode enabled: removing intermediate files for $target"
        rm -f "${target}.$filter_info.$gap_info.counts"
    fi

    if [ -z "$bgCov" ] || [ "$bgCov" = "0" ]; then
        norm_count="NA"
    elif [ ! -f "$len" ]; then
        norm_count="NA"
    else
        covLen=$(cat "$len")
        norm_count=$(echo "scale=5; 2 * $tgCount_count / $covLen / $bgCov" | bc)
    fi
    echo -e "$sample\t$refname\t${target}\t${norm_count}" >> "$prefix.tg.$filter_info.$mode.$gap_info.txt"
done