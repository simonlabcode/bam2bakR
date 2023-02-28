#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
input_snp=$4
output=$5
output2=$6
fragment_size=$7
minqual=$8
mut_tracks=$9
format=${10}
strand=${11}

# Create results/counts/
touch "$output2"

# Infer a decent fragment size
num_reads=$(samtools view -@ "$cpus" -c "$input")
fragment_size=$(echo "scale=0; $num_reads/$cpus | bc)
(( fragment_size++ ))


# Spread the work load so all cpus are working at all times
    declare $(samtools view -@ "$cpus"  "$input" \
    			| wc -l \
    			| awk -v fragment_size="$fragment_size" \
                      -v cpus="$cpus" '{
                                        nFrag = int($1 / fragment_size)                         # Calculate to how many fragments must the .bam be split
                                        if ($1 % fragment_size != 0) {nFrag = nFrag + 1}        # If there are some left over reads then increase number of fractions

                                        nBatchRun = int(nFrag / cpus)                           # In how many batches will the fragments be processed
                                        if (nFrag % cpus != 0) { nBatchRun = nBatchRun + 1}     # If there are some leftover fragments that will not fill up all cpus, then incerease number of batches

                                        newFragmentNumber = nBatchRun * cpus

                                        newFragmentSize = int($1 / newFragmentNumber)
                                        if ($1 % newFragmentNumber != 0) {newFragmentSize = newFragmentSize + 2}    # if there are still some reads leftover, then increase fragment size by 1 read-pair

                                        print "newFragmentNumber="newFragmentNumber
                                        print "newFragmentSize="newFragmentSize
                                    }')

    echo "* The number of fragments for sample $sample is $newFragmentNumber"
    echo "* The fragment size is set to $newFragmentSize alignments per fragment"



    for f in $(seq 1 $newFragmentNumber); do
        samtools view -H "$input" > ./results/counts/"$f"_"$sample".sam
    done &&


    samtools view "$input" \
    	| awk \
    		-v fragment_size="$newFragmentSize" \
    		-v sample="$sample" \
    		-f "./workflow/scripts/fragment_sam.awk" &&

    for f in $(seq 1 $newFragmentNumber); do
        samtools view -@ "$cpus" -o ./results/counts/"$f"_"$sample"_frag.bam ./results/counts/"$f"_"$sample".sam
        rm ./results/counts/"$f"_"$sample".sam
    done &&

    echo "* Aligned .sam file fragmented for sample $sample"


# remove any remaining sam files:

    #rm -f *.sam

#!/bin/bash
# Main sript for mutation calling in samples

## Need to add mutation call script to scripts!

# Call mutations
    parallel -j $cpus "python ./workflow/scripts/mut_call.py -b {1} \
                                              --mutType $mut_tracks \
                                              --minQual $minqual \
											  --SNPs "./results/snps/snp.txt" \
                                              --strandedness $strand \
                                              --reads $format" ::: ./results/counts/*_"$sample"_frag.bam \



    echo "** Mutations called for sample $sample"


# Combine output from fragments into single file
    # 1) _count.csv files
    awk 'FNR > 1 || NR == 1' ./results/counts/*_"$sample"_frag_counts.csv \
        | pigz -p $cpus > "$output"

    rm ./results/counts/*_"$sample"_frag_counts.csv



    echo "** Results fragments merged into final files"



	rm -f ./results/counts/*_"$sample"_frag.bam

	echo '* Cleaning up fragmented .bam files'
