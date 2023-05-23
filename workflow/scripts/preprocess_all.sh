#!/bin/bash

## Maybe I can have input be folders containing fastq files (of either SE or PE variety).
## Throw out step where I check if gz or not (force them to be fastq files)
## But other problem I have is that for whatever reasons, files get renamed
    ## I guess this is needed to ensure consistency of naming

## Realized I can get file names from fastq folder with some bash
    ## names=($(ls $fastq_dir)); names[0] is one of the read pair, names[1] is other if it exists
    cpus=$1 # number of cpus
    sample=$2 # sample names
    fastq_dir=$3 # input fastqs
    format=$4 # Paired or single end



if [ "$format" = "PE" ]; then

    output1=$5
    output2=$6
	adapter1=$7
	adapter2=$8

    # Create array of fastq file names
    declare -a fastqs
    for file in ${fastq_dir}/*.fastq
    do
        fastqs=(${fastqs[*]} "$file")
    done

    # remove duplicate fastq reads:
            #fastuniq \
             #   -i <(echo "${fastqs[0]}"; echo "${fastqs[1]}") \
             #   -o "$sample"_1u.fastq \
             #   -p "$sample"_2u.fastq &&

             #echo "* fastquniq finished for sample " $sample



    # trim reads
            echo "* Running cutadapt in parallel mode for sample $sample"


                cutadapt \
                    -a "$adapter1" \
                    -A "$adapter2" \
                    --minimum-length=20 \
                    --cores="$cpus" \
                    -o "$output1" \
                    -p "$output2"\
                    <(echo "${fastqs[0]}"; echo "${fastqs[1]}") &&
                echo "* cutadapt finished for " $sample





elif [ "$format" = "SE" ]; then

    output1=$5 # output1
	adapter=$6

    # Create array of fastq file names
    declare -a fastqs
    for file in ${fastq_dir}/*.fastq
    do
        fastqs=(${fastqs[*]} "$file")
    done

    # trim reads
            echo "* Running cutadapt in parallel mode for sample $sample"


                cutadapt \
                    -a "$adapter" \
                    --minimum-length=20 \
                    --cores="$cpus" \
                    -o "$output1" \
                    <(echo "${fastqs[0]}") &&
                echo "* cutadapt finished for " $sample


else
    echo "!!! format must be PE or SE !!!"
fi
