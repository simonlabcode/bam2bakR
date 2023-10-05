import csv
import argparse
import datetime


genomic_counts = snakemake.input.get("muts")
transcript_counts = snakemake.input.get("transcripts")

output = snakemake.output.get("merged")




transcripts_dict = {}

# Open sorted files and output file
with open(genomic_counts, 'r') as g_obj, open(transcript_counts, 'r') as t_obj, open(output, 'w', newline='') as output_file:

    reader_t = csv.reader(t_obj)
    reader_g = csv.reader(g_obj)
    writer = csv.writer(output_file)
    
    # Get the header and write to the output
    header_t = next(reader_t)
    header_g = next(reader_g)
    writer.writerow(header_g + [header_t[1]])
    
    # Initialize variables
    row_t = next(reader_t, None)
    for row_g in reader_g:
        # Loop until qnames match or file1 ends
        while row_t and row_t[0] < row_g[0]:
            row_t = next(reader_t, None)
        # Write matched row or row with 'NA'
        if row_t and row_t[0] == row_g[0]:
            writer.writerow(row_g + [row_t[1]])
            row_t = next(reader_t, None)
        else:
            writer.writerow(row_g + ['NA'])


