#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4
annotation=$5
mutcnt=$6


    samtools view -h -@ "$cpus" "$input" \
        | parallel \
            -L 2 \
            -j "$cpus" \
            --linebuffer \
            --joblog htseq_parallel.log \
            --roundrobin \
            --header '(@.*\n)*' \
            --pipe python $mutcnt \
                        -f sam \
                        --samout "$sample"_htseq.{#}_temp.sam \
                        -t gene,exon,exon \
                        -i gene_id,gene_id,gene_id \
                        -m union,union,intersection-strict \
                        -c "$sample"_GF_htseq.{#}_temp.txt,"$sample"_EF_htseq.{#}_temp.txt,"$sample"_XF_htseq.{#}_temp.txt \
                        - \
                        "$annotation"

    # Error catching
    if [ $(cut -f 7 htseq_parallel.log | grep '1' | wc -l ) -eq 0 ]; then
        echo '* HTSeq counting finished for sample' $sample
    else
        echo '!!! HTSeq counting failed!'
        exit 1
    fi



# Combine outputs from individual jobs
    # .sam file
    samtools view -H "$input" \
        | cat - "$sample"_htseq.*_temp.sam \
        | samtools sort \
            -@ "$cpus" \
            -n \
            -o "$output" -

    rm "${sample}"_htseq.*_temp.sam
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
                                          }' ${sample}_{1}_htseq.*_temp.txt  \
                            | LC_COLLATE=C sort -k1,1 > {2}.${sample}_htseq.txt" ::: GF EF XF \
                                                                                 :::+ gene exon_w_overlaps mature

    rm "${sample}"_GF_htseq.*temp.txt
	rm "${sample}"_EF_htseq.*temp.txt
	rm "${sample}"_XF_htseq.*temp.txt

    echo "* HTSeq counts files merged for sample $sample"

	mv exon_w_overlaps.${sample}_htseq.txt ./results/bams/
	mv gene.${sample}_htseq.txt ./results/bams
	mv mature.${sample}_htseq.txt ./results/bams/
