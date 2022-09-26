#!/bin/bash
#
# To make filtered bam files and tracks with different numbers of mutations

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
input1=$3
input2=$4
mut_tracks=$5
output=$6
genome_fasta=$7
WSL_b=$8


# Test if count files exist

    echo '* Creating files for each level of counting.'

    python ./workflow/scripts/count_to_tracks.py \
        -i $input1 \
        -s $sample


# Sort the starting .bam file by coordinate
    samtools sort -@ "$cpus" -o "$sample"_sort.bam "$input2"
    samtools index -@ "$cpus" "$sample"_sort.bam


# Make .chrom.sizes file from .bam file header (alternative to .genome file for toTDF)
    chrom_sizes=$(echo ${genome_fasta##*/} | cut -f 1 -d '.')".chrom.sizes"
    if [ ! -f $chrom_sizes ]; then
            samtools view -H "$sample"_sort.bam \
                | awk -v OFS="\t" ' $1 ~ /^@SQ/ {split($2, chr, ":")
                                                 split($3, size, ":")
                                                 print chr[2], size[2]}' > "$chrom_sizes"
    fi


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
    parallel -j 1 samtools view -@ "$cpus" \
                                     -b \
                                     -N "$sample"_{1}_{2}_reads.txt \
                                     -o "$sample"_{1}_{2}.bam \
                                     "$sample"_sort.bam ::: $muts \
                                                        ::: $(seq 0 5)


    # if [ $WSL_b ]; then
    #
    #     for iter in $(seq 0 5);
    #     do
    #         STAR \
    #             --runMode inputAlignmentsFromBAM \
    #             --runThreadN $cpus \
    #             --inputBAMfile "$sample"_"$muts"_"$iter".bam \
    #             --outWigType bedGraph \
    #             --outWigNorm None \
    #             --outWigStrand Stranded \
    #             --outFileNamePrefix ./"$sample"_"$muts"_"$iter"_
    #     done
    #
    #
    #
    #     # Take only unique component of track
    #     parallel -j "$cpus" "awk -v norm=${normVal} \
    #                                     '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
    #                                     ${sample}_{1}_{2}_Signal.Unique.{4}.out.bg \
    #                                     >> ${sample}.{1}.{2}.{5}.bedGraph" ::: $muts \
    #                                                                        ::: $(seq 0 5) \
    #                                                                        ::: + - \
    #                                                                        :::+ str1 str2 \
    #                                                                        :::+ pos min
    #
    #         #rm *.bg
    #
    # else
    #     # # Make tracks
    #     # parallel -j "$cpus" STAR \
    #     #                         --runMode inputAlignmentsFromBAM \
    #     #                         --inputBAMfile "$sample"_{1}_{2}.bam \
    #     #                         --outWigType bedGraph \
    #     #                         --outWigNorm None \
    #     #                         --outWigStrand Stranded \
    #     #                         --outFileNamePrefix ./"$sample"_{1}_{2}_ ::: $muts \
    #     #                                                                  ::: $(seq 0 5)
    #
    #      parallel -j "$cpus" "bedtools genomecov \
    #                                      -ibam ${sample}_{1}_{2}.bam \
    #                                      -bg \
    #                                      -pc \
    #                                      -strand {3} \
    #                                  | awk -v norm=${normVal} \
    #                                      '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
    #                                      >> ${sample}.{1}.{2}.{4}.bedGraph" ::: $muts \
    #                                                                                       ::: $(seq 0 5) \
    #                                                                                       ::: + - \
    #                                                                                       :::+ pos min
    #
    # fi

    parallel -j "$cpus" "bedtools genomecov \
                                    -ibam ${sample}_{1}_{2}.bam \
                                    -bg \
                                    -pc \
                                    -strand {3} \
                                | awk -v norm=${normVal} \
                                    '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
                                    >> ${sample}.{1}.{2}.{4}.bedGraph" ::: $muts \
                                                                                     ::: $(seq 0 5) \
                                                                                     ::: + - \
                                                                                     :::+ pos min




    # Make tdf files from the tracks
    parallel -j "$cpus" igvtools toTDF \
                                -f mean,max \
                                "$sample".{1}.{2}.{3}.bedGraph \
                                "$sample".{1}.{2}.{3}.tdf \
                                "$chrom_sizes" ::: $muts \
                                             ::: $(seq 0 5) \
                                             ::: pos min

    touch "$output"

    mv *.tdf ./results/tracks/

    rm -f *.bam
    rm -f *_reads.txt
    rm -f *.bedGraph
    rm -f *.bai
    rm -f *.out
    rm -f *.chrom.sizes
    rm -f igv*

    echo "* TDF track files created for sample $sample."
