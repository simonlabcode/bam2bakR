#!/bin/bash

# Source the paths and variables:
    cd "$MASTER_DIR"
    source $1
    sample=$2

    cd "$MASTER_DIR"/"$sample".dir


# align trimmed reads to the  genome
    if [[ -n $BOWTIE2_MOD ]]; then module load ${BOWTIE2_MOD}; fi
    if [[ -n $BISMARK_MOD ]]; then module load ${BISMARK_MOD}; fi
    if [[ -n $SAMTOOLS_MOD ]]; then module load ${SAMTOOLS_MOD}; fi

    echo "* Aligning reads with Bowtie2 and Bismark for sample $sample"

    if [[ "$READS" == "F" ]]; then
        $BISMARK \
            $BISMARK_INDEX \
            --multicore `echo $(( $cpus/4 + 1 ))` \
            --se "$sample".t.fastq \
            --local \
            -N 0 \
            -L 20 \
            --un
    elif [[ "$READS" = "FR" ]]; then
        $BISMARK \
            $BISMARK_INDEX \
            --multicore `echo $(( $cpus/4 + 1 ))` \
            -1 "$sample".t.r1.fastq \
            -2 "$sample".t.r2.fastq \
            --local \
            -N 0 \
            -L 20 \
            --un
    elif [[ "$READS" = "RF" ]]; then
        $BISMARK \
            $BISMARK_INDEX \
            --multicore `echo $(( $cpus/4 + 1 ))` \
            -1 "$sample".t.r2.fastq \
            -2 "$sample".t.r1.fastq \
            --local \
            -N 0 \
            -L 20 \
            --un
    else
        echo "! No FR/RF/F method recognized for " $sample
        exit 1
    fi


    if [[ "$FORMAT" = "SE" ]]; then
        rm "$sample".t.fastq
        mv "$sample".t_bismark_bt2.bam "$sample"Aligned.out.bam

    else
        rm "$sample".t.r[12].fastq
        mv "$sample".t.r2_bismark_bt2_pe.bam "$sample"Aligned.out.bam
    fi


    echo "* Alignment script finished for " $sample
