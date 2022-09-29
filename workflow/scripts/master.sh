#!/bin/bash

# Source the paths and variables:
cpus=$1
cBout=$2
keepcols=$3
mut_tracks=$4

#day=$(date +"%y%m%d")

# create base count name: "TC" => "nT"
base=$(echo $mut_tracks | awk -v OFS="," '
                               {
                                    split($1,a,",")
                               }
                               END {
                                   for (base in a) {
                                       $base = "n"substr(a[base],1,1)
                                   }
                                   print $0
                               }')


keepcols=${keepcols}","${base}","${mut_tracks}

#cd ./results/counts

# Read all _counts.csv.gz files and save them as master-DATE.csv.gz and cB-DATE.csv.gz
parallel -j 1 --compress --plus "cat <(echo Filename:{1%_counts.csv.gz}) <(pigz -d -k -c -p $cpus {1})" ::: ./results/counts/*_counts* \
    | awk -v OFS="," '
            $1 ~ /Filename/ {
                split($1, sample, ":")
                next
            }
            NR == 2 {
                header = $0
                print "sample", $0
                next
            }
            $0 == header {
                next
            }
            {
                print sample[2], $0
            }' \
    | awk -v colNames="$keepcols" '
            BEGIN {
                FS=","
                ncol=split(colNames,tmp)
                for (i in tmp) {
                    names[tmp[i]]
                }
            }
            NR == 1 {
                for (i=1; i<=NF; i++) {
                    if ($i in names) {
                        f[++nf] = i
                    }
                }
            }
            {
                out = ""
                for (i=1; i<=ncol; i++) {
                    out = (out=="" ? "" : out ",") $(f[i])
                }
                print out
            } ' \
    | awk -v OFS="," -v FS="," '
            NR == 1 {
                print $0, "n"
                next
            }
            NR == 2 {
                sample = $1
            }
            sample != $1 {

                for (row in count) {
                    print row, count[row]
                }
                sample = $1
                delete count
            }
            {
                ++count[$0]
            }
            END {
                for (row in count) {
                    print row, count[row]
                }
            } '\
    | pigz -p $cpus > "$cBout"


# Explanation of the data flow:
#
#   cat (echo [filename as identiefier of data origin] + gzip [decompress file]) [use Parallels to create a single stream of alternating filename line + csv lines]
#    |
#    --> awk [use filename as a value for new column named "sample"; keep only the very first .csv header]
#         |
#         --> tee  --> gzip [compress to master file]
#              |
#              --> awk [look which column names in header match keepcols and print only those]
#                   |
#                   --> awk [calculate occurences of reads with the same values]
#                        |
#                        --> gzip [compress output to cB file]

echo "** cB file created: cB.csv.gz"

rm -f *temp*
