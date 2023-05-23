__author__ = "Isaac Vock"
__copyright__ = "Copyright 2023, Isaac Vock"
__email__ = "isaac.vock@yale.edu"
__license__ = "MIT"


import pysam
import csv
import argparse
import datetime

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is python implementation of TimeLapse mutation calling')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-b', '--bam', type=str, required=True, metavar = 'in_file.bam',
                    help='Bam file to process')

args = parser.parse_args()

inputName = args.bam.split('.bam')[0]           # name without .bam suffix

samfile = pysam.AlignmentFile(args.bam, "rb")

myfile = open(inputName + '_rsem.csv', 'w', newline = '')
wr = csv.writer(myfile)

header = ['qname', 'TF', 'pt']
wr.writerow(header)

for r in samfile:
    r_info = [r.query_name, r.reference_name, r.get_tag('ZW')]
    wr.writerow(r_info)

myfile.close()
