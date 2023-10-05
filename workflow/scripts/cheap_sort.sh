#!/bin/bash

# Variables
input_file=$1
output_file=$2
sort_column=$3 # Column number to sort by
lines=$4
gzipped=$5
dir=$6

if [ "$gzipped" = "TRUE" ]; then
    
    pigz -d -c $input_file > $dir/input_decompressed.csv

    # Extract header
    header=$(head -n 1 $dir/input_decompressed.csv)
    tail -n +2 $dir/input_decompressed.csv > $dir/no_header.csv

    #rm -f ./results/read_to_transcripts/input_decompressed.csv


else
    
    # Extract header
    header=$(head -n 1 $input_file)
    tail -n +2 $input_file > $dir/no_header.csv


fi


# Stage 1: Split large CSV (without header) into smaller files
split -l $lines $dir/no_header.csv $dir/chunk_ # Adjust the number of lines (10000) as needed

# Stage 2: Sort each smaller file
for file in $dir/chunk_*
do
    sort -t, -k$sort_column -V -o $file $file # Adjust delimiter (-t) and sort column (-k) as needed
done

# Stage 3: Merge sorted files into a single sorted file (without header)
sort -m -t, -k$sort_column -V -o $dir/sorted_no_header.csv $dir/chunk_*

# Add header back to the sorted file
echo $header > $output_file
cat $dir/sorted_no_header.csv >> $output_file

# Clean up temporary files
#rm -f $dir/sorted_no_header.csv
#rm -f $dir/chunk_*
#rm -f $dir/no_header.csv