#!/bin/bash

# number of cpus
# number of control samples
# list of control samples

control_samples=("$@")

unset control_samples[0]
unset control_samples[1]
unset control_samples[2]
unset control_samples[3]



cpus=$1
nsamps=$2
output_txt=$3
genome_fasta=$4




if [ $nsamps > 0 ]; then
    # Loop through control samples:
        for cs in ${control_samples[@]}; do
            name=$(echo "$cs" | cut -d '/' -f 3 | rev | cut -c8- | rev)
            samtools sort -@ "$cpus" -o "$name"_sort.bam "$cs"
            samtools index -@ "$cpus" "$name"_sort.bam
        done

        for cs in ${control_samples[@]}; do
            NAMES+=($(echo "$cs" | cut -d '/' -f 3 | rev | cut -c8- | rev))
        done

    # Parallelize SNPs calling. Each chromosome in each .bam file is processed as separate job
        # Note: This approach does not give the same snp.txt result. In 2199712 SNPs there were 16 different.
        # {1} : path to fasta file
        # {2} : samtools view -H ${control_samples[0]}_sort.bam | awk ' $1 == "@SQ" {split($2,a,":"); print a[2]}' : Extracts chromosome names from .bam file header
        # {3} : ${control_samples[@]/%/_sort.bam}                                                    : Appends "_sort.bam" to the end of control names and prints
        parallel -j "$cpus" "samtools mpileup -uf {1} \
                                                   -r {2} {3} \
                                | bcftools call -mv" ::: $genome_fasta \
                                                      ::: $(samtools view -H ${NAMES[0]}_sort.bam \
                                                                    | awk ' $1 == "@SQ" {split($2,a,":"); print a[2]}') \
                                                      ::: ${NAMES[@]/%/_sort.bam} > snp.vcf


        # Note: Easier and also fast option would be:  bcftools mpileup --threads $cpus -f $genome_fasta "$cs"_sort.bam | bcftools call --threads $cpus-mv > snp-"$cs".vcf


    # Clean this up for all possible mutations:
        # Note: Filtering now inludes cases where both alleles are mutated
        awk '$1 !~ /^#/ && length($4) == 1 {if (length($5) == 1) {print $4":"$5":"$1":"$2}
                                            else if (length($5) == 3 && $5 ~ /,/) {split($5, mut, ",")
                                                                                   print $4":"mut[1]":"$1":"$2
                                                                                   print $4":"mut[2]":"$1":"$2}
                                            }' snp.vcf \
            | sort \
            | uniq > $output_txt

        echo '* SNPs called and snp.txt generated'
fi
