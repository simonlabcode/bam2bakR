#!/bin/bash
#
# To make filtered bam files and tracks with different numbers of mutations

# Get parameters
# This allows to be script run as: sbatch tracks.sh configFile sampleName
    if [[ -n $1 ]]; then config=$1; fi
    if [[ -n $2 ]]; then sample=$2; fi

    source $config

# Set normalization value
normVal='1'

## Dependencies
# Python
# Pysam
# Igvtools
# parallel
# Samtools
# Star

cpus=$1
sample=$2
input1=$3 # counts.csv files
input2=$4 # _tl.bam files
output=$5
mut_tracks=$6
TLcount2tracks=$7


# Test if count files exist

    echo '* Creating files for each level of counting.'

    python $TLcount2tracks \
        -i "$input1" \
        -s $sample \
        --mutType $mut_tracks


# Sort the starting .bam file by coordinate
    samtools sort -@ "$cpus" -o "$sample"_sort.bam "$input2"
    samtools index -@ "$cpus" "$sample"_sort.bam


# Make .chrom.sizes file from .bam file header (alternative to .genome file for toTDF)
    samtools view -H "$sample"_sort.bam \
        | awk -v OFS="\t" ' $1 ~ /^@SQ/ {split($2, chr, ":")
                                         split($3, size, ":")
                                         print chr[2], size[2]}' > "$chrom_sizes"



    muts=$(echo $mut_tracks | tr ',' ' ')

    for b in $muts; do
        if [ $b == "GA" ]; then
            colVal[0]='200,200,200'
            colVal[1]='120,188,230'
            colVal[2]='65,125,195'
            colVal[3]='36,110,182'
            colVal[4]='27,78,165'
            colVal[5]='18,50,120'
        else
            colVal[0]='200,200,200'
            colVal[1]='250,150,150'
            colVal[2]='250,0,0'
            colVal[3]='150,0,0'
            colVal[4]='100,0,0'
            colVal[5]='50,0,0'
        fi

        # Make track headers
        for count in $(seq 0 5); do

            echo "track type=bedGraph name=\" $sample $b $count "plus \" description=\" $sample $b $count "positive strand\" visibility=full autoScale=off windowingFunction=mean viewLimits=0:10 color="${colVal[$count]} " altColor=10,10,10 priority=20" > "$sample"."$b"."$count".pos.bedGraph
            echo "track type=bedGraph name=\" $sample $b $count "minus \" description=\" $sample $b $count "minus strand\" visibility=full negateValues=on autoScale=off windowingFunction=mean viewLimits=-10:0 color=10,10,10 altColor="${colVal[$count]} " priority=20" > "$sample"."$b"."$count".min.bedGraph


        done
    done


    # Filter the reads
    parllel -j 1 samtools view -@ "$cpus" \
                                     -b \
                                     -N "$sample"_{1}_{2}_reads.txt \
                                     -o "$sample"_{1}_{2}.bam \
                                     "$sample"_sort.bam ::: $muts \
                                                        ::: $(seq 0 5)

        # Make tracks
        parallel -j "$cpus" STAR \
                                --runMode inputAlignmentsFromBAM \
                                --inputBAMfile "$sample"_{1}_{2}.bam \
                                --outWigType bedGraph \
                                --outWigNorm None \
                                --outWigStrand Stranded \
                                --outFileNamePrefix {1}_{2}_ ::: $muts \
                                                               ::: $(seq 0 5)



        # Take only unique component of track
        parallel -j "$cpus" "awk -v norm=${normVal} \
                                        '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
                                        {1}_{2}_Signal.Unique.{4}.out.bg \
                                        >> ${sample}.{1}.{2}.{5}.bedGraph" ::: $muts \
                                                                                         ::: $(seq 0 5) \
                                                                                         ::: + - \
                                                                                         :::+ str1 str2 \
                                                                                         :::+ pos min

        rm *.bg



    # Make tdf files from the tracks
    parallel -j "$cpus" igvtools toTDF \
                                -f mean,max \
                                "$sample".{1}.{2}.{3}.bedGraph \
                                "$output" \
                                "$chrom_sizes" ::: $muts \
                                             ::: $(seq 0 5) \
                                             ::: pos min






    echo "* TDF track files created for sample $sample."

    rm -f *.bam
    rm -f *_reads.txt
    rm -f *.bedGraph
