#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4
format=$5


samtools fixmate -@ "$cpus" "$input" "$sample"_fixed_mate.bam

	if [ "$format" = "NU" ]; then
	samtools view -@ "$cpus" -h "$sample"_fixed_mate.bam \
		| awk awk '$1 ~ /^@/ {print}
                   (($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163)) || (($2 == 355 || $2 == 403) || ($2 == 339 || $2 == 419))  {print}' > "$sample".f.sam
    elif [ "$format" = "PE" ]; then
	    samtools view -@ "$cpus" -q 2 -h "$sample"_fixed_mate.bam \
            | awk '$1 ~ /^@/ {print}
                  ($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163) {print}' > "$sample".f.sam
    elif [ "$format" = "SE" ]; then
        samtools view -@ "$cpus" -q 2 -h "$sample"_fixed_mate.bam \
            | awk '$1 ~ /^@/ {print}
                  ($2 == 0 || $2 == 16) {print}' > "$sample".f.sam
	else
		samtools view -@ "$cpus" -h "$sample"_fixed_mate.bam > "$sample".f.sam
	fi &&
    echo "* Reads filtered for sample $sample"
	



samtools sort -@ "$cpus" -n -o "$output" "$sample".f.sam

rm "$sample".f.sam
rm "$sample"_fixed_mate.bam