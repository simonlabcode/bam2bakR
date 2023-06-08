#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4
output2=$5
annotation=$6
strand=$7
count_script=$8
flattened=$9


if [ "$strand" = "F" ]; then

    strand="yes"

else

    strand="reverse"

fi

    # Will create the ./results/htseq
    touch "$output2"


    if [ "$flattened" = "True" ]; then
    

        samtools view -h -@ "$cpus" "$input" \
        | parallel \
            -L 2 \
            -j "$cpus" \
            --linebuffer \
            --joblog "$sample"_htseq_parallel.log \
            --roundrobin \
            --header '(@.*\n)*' \
            --pipe python $count_script \
                        -f sam \
                        --samout ./results/htseq/"$sample"_htseq.{#}_temp.sam \
                        -t aggregate_gene,exonic_part,exonic_part \
                        -i gene_id,transcripts,exon_id \
                        -m union,union,intersection-strict \
                        -s "$strand" \
                        -c ./results/htseq/"$sample"_GF_htseq.{#}_temp.txt,./results/htseq/"$sample"_EF_htseq.{#}_temp.txt,./results/htseq/"$sample"_XF_htseq.{#}_temp.txt \
                        - \
                        "$annotation"

    else


        samtools view -h -@ "$cpus" "$input" \
        | parallel \
            -L 2 \
            -j "$cpus" \
            --linebuffer \
            --joblog "$sample"_htseq_parallel.log \
            --roundrobin \
            --header '(@.*\n)*' \
            --pipe python $count_script \
                        -f sam \
                        --samout ./results/htseq/"$sample"_htseq.{#}_temp.sam \
                        -t transcript,exon,exon \
                        -i gene_id,gene_id,gene_id \
                        -m union,union,intersection-strict \
                        -s "$strand" \
                        -c ./results/htseq/"$sample"_GF_htseq.{#}_temp.txt,./results/htseq/"$sample"_EF_htseq.{#}_temp.txt,./results/htseq/"$sample"_XF_htseq.{#}_temp.txt \
                        - \
                        "$annotation"

    fi


    # Error catching
    if [ $(cut -f 7 "$sample"_htseq_parallel.log | grep '1' | wc -l ) -eq 0 ]; then
        echo '* HTSeq counting finished for sample' $sample
    else
        echo '!!! HTSeq counting failed!'
        exit 1
    fi



# Combine outputs from individual jobs
    # .sam file
    samtools view -H "$input" \
        | cat - ./results/htseq/"$sample"_htseq.*_temp.sam \
        | samtools sort \
            -@ "$cpus" \
            -n \
            -o "$output" -

    echo "* HTSeq .sam files merged for sample $sample"

### Need to make this a separate rule!
    # gene counts
    parallel -j "$cpus" "awk -v 'OFS=\t' '
                                          {
                                                gene[\$1] += \$2
                                          }
                                          END {
                                                for (name in gene) {
                                                    print name, gene[name]
                                                }
                                          }' ./results/htseq/${sample}_{1}_htseq.*_temp.txt  \
                            | LC_COLLATE=C sort -k1,1 > ./results/htseq/{2}.${sample}_htseq.txt" ::: GF EF XF \
                                                                                 :::+ gene exon_w_overlaps mature


    echo "* HTSeq counts files merged for sample $sample"

    mv ${sample}_htseq_parallel.log ./results/htseq/
    rm ./results/htseq/${sample}*temp*