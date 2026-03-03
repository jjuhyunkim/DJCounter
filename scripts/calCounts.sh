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
log=""

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
        --clean)
            clean=true
            shift 1
            ;;
        --log)
            log="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --sample <sample> --bam <bam_file> --ref <reference> --mode <mode> [--noGap --filter <filter> --threads <threads> --targetlist <targetlist> --log <log_file>]"
            exit 1
            ;;
    esac
done
 
if [ -z "$sample" ] || [ -z "$refname" ] || [ -z "$bam" ]; then
    echo "Usage: $0 --sample <sample> --bam <bam_file> --ref <reference> "
    echo "  --sample: sample name"
    echo "  --bam: path to BAM file"
    echo "  --ref: reference name (GRCh38 by default, CHM13 and hg19 are also supported)"
    echo "  --noGap: exclude gap regions in calculations (False by default) --ref CHM13 will automatically set this to false"
    echo "  --threads: number of threads to use (10 by default)"
    echo "  --mode: calculation mode (fast-accurate by default, options: fast, fast-accurate, high-accurate, high-accurate-coverage)"
    echo "  --log: path to log file (optional)"
    exit 1
fi
prefix=${sample}.${refname}

# check if the mode is valid
if [ "$mode" != "fast" ] && [ "$mode" != "fast-accurate" ] && [ "$mode" != "high-accurate" ]  &&  [ "$mode" != "high-accurate-coverage" ] ; then
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
    if [ "$targetlist" == "DJ_filt" ]; then
        targetlist="DJ_on_chr13"
    fi
    noGap=false
    autosomeBed=$roi_dir/CHM13/autosome.bed
    backgroundLen=$roi_dir/CHM13/autosome.len

else
    echo "Unknown reference name: $refname"
    exit 1
fi

autosomeLen=$(awk '{sum += $3 - $2} END {print sum}' $autosomeBed)

# Check whether all required contigs are present in the BAM header.
for target in $targetlist
do
    if [ ! -s "$roi_dir/$refname/$target.bed" ]; then
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

# Index BAM if the index file does not exist
if [ ! -f ${bam}.bai ] && [ ! -f ${bam}.crai ] && [ ! -f ${bam}.csi ]; then
	samtools index -@ $threads $bam
fi

# Preset parameters based on the mode
if [ "$mode" = "high-accurate" ] || [ "$mode" == "high-accurate-coverage" ]; then
    accurate_info="byChr"
else 
    accurate_info="byWGS"
fi

if [ "$mode" = "fast" ]; then
    filter_info="noFilter"
fi

# Shorten output file
OUT=${prefix}.${filter_info}.${gap_info}.${mode}.${accurate_info}

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
echo -e "Output prefix: $OUT"
echo ""



# Collect background read count for normalization per modes

if [ -s "$OUT.readCount.bg.txt" ]; then
    echo "Removing prior bg read count file"
    rm "$OUT.readCount.bg.txt"
fi

if [ "$mode" = "fast" ]; then
    if [ ! -f "$OUT.idxstats" ]; then
        samtools idxstats -@ "$threads" "$bam" > "$OUT.idxstats"
    fi
    bgCov=$(grep -E '^chr([1-9]|1[0-9]|2[0-2])([[:space:]]|$)' "$OUT.idxstats" | awk '{print $3/$2}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}')
    echo $bgCov > "$OUT.readCount.bg.txt"
    
elif [ "$mode" = "fast-accurate" ];then
    echo -e "Calculating background read count without per-chromosome processing"
    autosomeReadCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L "$autosomeBed" $bam)
    bgCov=$(echo "scale=5; $autosomeReadCount / $autosomeLen" | bc)
    echo $bgCov > "$OUT.readCount.bg.txt"
    
elif [ "$mode" = "high-accurate" ]; then
    if [ ! -f "$OUT.readCount.txt" ]; then
        echo -e "Calculating background read count with per-chromosome processing"
        echo -e "This might take longer time but will be more accurate, especially for samples with uneven coverage across chromosomes."
        
        for chr in {1..22}; do
            echo "Processing chromosome $chr"
            sed -n ${chr}p $autosomeBed > $OUT.tmp.bed
            readCount=$(samtools view -@ ${threads} -F ${filter_flag} -c -L $OUT.tmp.bed "$bam")
            echo -e "chr${chr}\t$readCount" >> $OUT.readCount.txt
        done
        rm $OUT.tmp.bed
    fi
    bgCov=$(join -t $'\t' -1 1 -2 1  "$OUT.readCount.txt" "$backgroundLen" | awk '{print $2/$3}'| sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}') 
    echo $bgCov > "$OUT.readCount.bg.txt"

elif [ "$mode" = "high-accurate-coverage" ]; then
    if [ ! -s $OUT.cov ]; then
        samtools coverage $bam > $OUT.cov
    fi
    bgCov=$(grep -E '^chr([1-9]|1[0-9]|2[0-2])([[:space:]]|$)' $OUT.cov | cut -f 7 | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)];}')
    echo $bgCov > "$OUT.readCount.bg.txt"
fi

echo "Background coverage (bgCov): $bgCov"

if [ -s $OUT.readCount.tg.txt ]; then
    echo "Removing prior target read count file"
    rm $OUT.readCount.tg.txt
fi

for target in $targetlist
do
    echo "Target: $target"
    bed="${targetDir}/$target.bed"
    len="${targetDir}/$target.len"
    
    if [ ! -s "$bed" ]; then
        echo "Warning: BED file $bed does not exist or is empty, skipping target $target."
        continue
    fi

    # covLen=$(awk '{sum += $3 - $2} END {print sum}' "$bed")
    covLen=$(cat "$len")
    if [ -z "$covLen" ] || [ "$covLen" -eq 0 ]; then
        echo "Warning: covLen is zero for $bed, skipping."
        continue
    fi

    # echo "Processing target: $target"
    if [ "$mode" = "high-accurate-coverage" ]; then
        # Using depth instead of count here
        if [ ! -s "$OUT.$target.counts" ]; then
            echo "Calculating depth for $target"
            samtools depth -@ "$threads" -b "$bed" "$bam" > "$OUT.$target.counts"
        elif [ -s "$OUT.$target.counts" ]; then
            echo "Depth file already exists for $target, skipping calculation."
        fi
        tgCount_count=$(awk 'BEGIN { SUM=0 } { SUM+=$3 } END { print SUM }' "$OUT.$target.counts")
        echo $tgCount_count > $OUT.$target.counts
    else
        # Using readCount for other modes
        if [ -s "$OUT.$target.counts" ]; then
            echo "Count file already exists for $target, skipping calculation."
            tgCount_count=$(cat "$OUT.$target.counts")
        else
            if [ "$filter_info" != "noFilter" ]; then
                # for fast-accurate and high-accurate modes (filter enabled by default)
                echo "Calculating counts for $target with filter: ${filter_flag}"
                tgCount_count=$(samtools view -@ "$threads" -F ${filter_flag} -c -L "$bed" "$bam")
            else
                # for fast mode
                echo "Calculating counts for $target including all alignments (no filter)"
                tgCount_count=$(samtools view -@ "$threads" -c -L "$bed" "$bam")
            fi
            echo $tgCount_count > $OUT.$target.counts
        fi
    fi

    if [ "$clean" = true ]; then
        echo "Clean mode enabled: removing intermediate files for $target"
        rm -f "$OUT.$target.counts"
    fi

    if [ -z "$bgCov" ] || [ "$bgCov" = "0" ]; then
        norm_count="NA"
    elif [ ! -f "$len" ]; then
        norm_count="NA"
    else
        norm_count=$(echo "scale=5; 2 * $tgCount_count / $covLen / $bgCov" | bc)
    fi
    # Output format: sample, refname, target, target read count, background coverage, normalized count
    echo -e "$sample\t$refname\t${target}\t${tgCount_count}\t${bgCov}\t${norm_count}" >> $OUT.readCount.tg.txt
done
