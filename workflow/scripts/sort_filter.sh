#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4
output2=$5
output3=$6
format=$7

# Make sure bam file is query sorted rather than coordinate sorted
samtools sort -@ "$cpus" -n "$input" | samtools fixmate -@ "$cpus" - - | samtools view -@ "$cpus" -b - > "$output2"

	if [ "$format" = "NU" ]; then
	samtools view -@ "$cpus" -h "$output2" \
		| awk awk '$1 ~ /^@/ {print}
                   (($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163)) || (($2 == 355 || $2 == 403) || ($2 == 339 || $2 == 419))  {print}' > "$output3"
    elif [ "$format" = "PE" ]; then
	    samtools view -@ "$cpus" -q 2 -h "$output2" \
            | awk '$1 ~ /^@/ {print}
                  ($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163) {print}' > "$output3"
    elif [ "$format" = "SE" ]; then
        samtools view -@ "$cpus" -q 2 -h "$output2" \
            | awk '$1 ~ /^@/ {print}
                  ($2 == 0 || $2 == 16) {print}' > "$output3"
	else
		samtools view -@ "$cpus" -h "$sample"_fixed_mate.bam > "$output3"
	fi &&
    echo "* Reads filtered for sample $sample"
	



samtools sort -@ "$cpus" -n -O bam -o "$output" "$output3"
