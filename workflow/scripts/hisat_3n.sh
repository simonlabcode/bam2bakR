#!/bin/bash

# Source the paths and variables:
cpus=$1 # number of cpus
sample=$2 # sample names
format=$3 # Paired or single end
reads=$4 # Forward or reverse
chr_tag=$5
HISAT_3N_INDEX=$6
HISAT_3N_PATH=$7
mut_tracks=$8
yale=$9

if [ "$format" = "PE" ]; then
    input1=${10}
    input2=${11}
    output=${12}

else
    input1=${10}
    output=${11}
fi

if [ "$yale" = "True" ]; then
    module purge
    module load StdEnv
    module load HISAT-3N
fi



    echo "* Aligning reads with HISAT-3n for sample $sample"

    if [ "$chr_tag" = 'True' ]; then
        echo "* chr tag will be included in alignment for " $sample
    fi


    if [[ "$format" = "PE" ]]; then
        if [[ "$reads" = "F" ]]; then
        
            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ "$chr_tag" = 'True' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness FR \
                | samtools view -Sbh -o "$output"
        else
            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ "$chr_tag" = 'True' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness RF \
                | samtools view -Sbh -o "$output"
        fi

    elif [[ "$format" = "SE" ]]; then
        if [[ "$reads" = "F" ]]; then

            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -U "$input1" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ "$chr_tag" = 'True' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness F \
                | samtools view -Sbh -o "$output"
        else

            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -U "$input1" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ "$chr_tag" = 'True' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness R \
                | samtools view -Sbh -o "$output"

        fi
    else
        echo "! No PE/SE FORMAT parameter recognized for " $sample
        exit 1
    fi &&


    echo "* Alignment script finished for " $sample
