__author__ = "Isaac Vock"
__copyright__ = "Copyright 2023, Isaac Vock"
__email__ = "isaac.vock@yale.edu"
__license__ = "MIT"


import pysam
import sys

input_bam = snakemake.input.input_bam
output_bam = snakemake.output.output_bam

def remove_tags(input_bam, output_bam):
    # Open the input BAM file
    input_file = pysam.AlignmentFile(input_bam, "rb")

    # Create a new BAM file for writing
    output_file = pysam.AlignmentFile(output_bam, "wb", header=input_file.header)

    # Iterate over each read in the input BAM file
    for read in input_file:
        # Remove the 'jI' and 'jM' tags, if present
        read.tags = [(tag, value) for tag, value in read.tags if tag not in ['jI', 'jM', 'GF', 'XF', 'EF']]

        # Write the modified read to the output BAM file
        output_file.write(read)

    # Close the BAM files
    input_file.close()
    output_file.close()

# Execute the tag removal function
remove_tags(input_bam, output_bam)